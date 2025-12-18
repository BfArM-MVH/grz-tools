import json
import logging
import queue
import subprocess
import threading

daemon_logger = logging.getLogger("daemon_monitor")

finish_sentinel = object()
submission_queue = queue.Queue()

daemon_keepalive_marker = "<results>/daemon/daemon_keepalive_marker.marker"
if "daemon" in sys.argv:
    submission_queue.put(daemon_keepalive_marker)


def _run_grzctl_db_command(*args):
    """Helper to run grzctl db commands from monitoring thread"""
    try:
        db_config_path = config["config_paths"]["db"]
        full_cmd = ["grzctl", "db", "--config-file", db_config_path, *args]
        result = subprocess.run(
            full_cmd, check=True, text=True, capture_output=True, timeout=120
        )
        return result
    except subprocess.CalledProcessError as e:
        daemon_logger.error(
            "A grzctl command failed: %s\n%s", " ".join(full_cmd), e.stderr
        )
        return None
    except Exception as e:
        daemon_logger.error(
            "An unexpected error occurred while running a grzctl command: %s", e
        )
        return None


def monitor_and_queue_submissions(shutdown_event):
    """
    Runs in background thread, replicates discovery logic (scan → sync → select)
    to find new submissions and puts the respective target paths onto snakemake queue.
    """
    interval = int(config["monitor"]["interval"])
    inbox_configs = config["config_paths"]["inbox"]
    already_queued = set()

    daemon_logger.info(
        "Starting monitoring thread with a scan interval of %s seconds.",
        interval,
    )
    print(f"Starting monitoring thread with a scan interval of {interval} seconds.")

    try:
        while not shutdown_event.is_set():
            if (
                not submission_queue.empty()
                and submission_queue.queue[0] is finish_sentinel
            ):
                daemon_logger.info(
                    "Finish sentinel detected in queue. Shutting down monitoring thread."
                )
                break
            daemon_logger.info(
                "Checking for new submissions in all configured inboxes."
            )

            # scan all inboxes
            all_s3_submissions = []
            for submitter, inboxes in inbox_configs.items():
                for inbox, config_path in inboxes.items():
                    try:
                        result = subprocess.run(
                            [
                                "grzctl",
                                "list",
                                "--config-file",
                                config_path,
                                "--json",
                                "--show-cleaned",
                            ],
                            check=True,
                            text=True,
                            capture_output=True,
                            timeout=120,
                        )
                        submissions = json.loads(result.stdout)
                        for submission in submissions:
                            submission["origin"] = {
                                "submitter_id": submitter,
                                "inbox": inbox,
                            }
                        all_s3_submissions.extend(submissions)
                    except Exception as e:
                        daemon_logger.error(
                            "Failed to scan the inbox for submitter '%s', inbox '%s': %s",
                            submitter,
                            inbox,
                            e,
                        )

            # sync with db
            initial_db_submission_list = _run_grzctl_db_command("list", "--json")
            if initial_db_submission_list is None:
                daemon_logger.critical(
                    "The database is unavailable. Shutting down monitoring."
                )
                submission_queue.put(finish_sentinel)
                return

            db_states = {
                e["id"]: e.get("latest_state", {}).get("state", "").casefold()
                for e in json.loads(initial_db_submission_list.stdout or "[]")
            }

            for submission in all_s3_submissions:
                submission_id, s3_state = (
                    submission["submission_id"],
                    submission["state"],
                )

                if s3_state == "complete":
                    target_db_state = "uploaded"
                elif s3_state == "incomplete":
                    target_db_state = "uploading"
                else:
                    daemon_logger.debug(
                        "Skipping submission '%s' because its S3-state is '%s'.",
                        submission_id,
                        s3_state,
                    )
                    continue

                if submission_id not in db_states:
                    daemon_logger.info(
                        "Found a new submission '%s' and registered it in the database with state '%s'.",
                        submission_id,
                        target_db_state,
                    )
                    _run_grzctl_db_command("submission", "add", submission_id)
                    _run_grzctl_db_command(
                        "submission", "update", submission_id, target_db_state
                    )
                elif (
                    db_states.get(submission_id) == "uploading"
                    and target_db_state == "uploaded"
                ):
                    daemon_logger.info(
                        "Submission '%s' has completed its upload. The database state is being updated to 'uploaded'.",
                        submission_id,
                    )
                    _run_grzctl_db_command(
                        "submission", "update", submission_id, "uploaded"
                    )

            # select pending submissions
            db_submissions_list = _run_grzctl_db_command("list", "--json")
            if db_submissions_list is None:
                daemon_logger.critical(
                    "The database became unavailable after the sync operation. Shutting down monitoring."
                )
                submission_queue.put(finish_sentinel)
                return

            final_db_states = {
                e["id"]: e.get("latest_state", {}).get("state", "").casefold()
                for e in json.loads(db_submissions_list.stdout or "[]")
            }
            sub_id_to_origin = {
                s["submission_id"]: s["origin"] for s in all_s3_submissions
            }

            new_targets_found = []
            for submission_id, state in final_db_states.items():
                if state == "uploaded":
                    if origin := sub_id_to_origin.get(submission_id):
                        target_path = f"<results>/{origin['submitter_id']}/{origin['inbox']}/{submission_id}/processed"
                        if target_path not in already_queued:
                            new_targets_found.append(target_path)

            if new_targets_found:
                daemon_logger.info(
                    "Found %d new submission(s) ready for processing; adding them to the queue.",
                    len(new_targets_found),
                )
                for target in new_targets_found:
                    submission_queue.put(target)
                    already_queued.add(target)
            else:
                daemon_logger.info(
                    "No new submissions were found that are ready for processing in this cycle."
                )

            shutdown_event.wait(timeout=interval)

    except (KeyboardInterrupt, SystemExit):
        daemon_logger.info(
            "Monitoring thread received an interrupt signal. Signalling workflow to finish processing."
        )
    finally:
        daemon_logger.info(
            "Placing the finish sentinel on the queue to allow a graceful shutdown."
        )


monitor_thread = None
if "daemon" in sys.argv:
    monitor_thread = threading.Thread(
        target=monitor_and_queue_submissions, args=(shutdown_event,), daemon=True
    )


rule daemon_keepalive:
    localrule: True
    output:
        marker=touch(temp(daemon_keepalive_marker)),
    input:
        rules.init_db.output.marker,
    params:
        delay=lambda wildcards: int(config["monitor"]["interval"]) * 2 + 3,
    shell:
        "sleep {params.delay};"


rule daemon:
    """
    Consumes submissions from monitoring queue and sends them off for processing.
    """
    localrule: True
    input:
        rules.daemon_keepalive.output.marker,
        from_queue(submission_queue, finish_sentinel=finish_sentinel),
    default_target: True
