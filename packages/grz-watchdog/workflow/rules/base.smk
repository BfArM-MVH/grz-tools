import os
from snakemake.io import temp, directory, touch


rule init_db:
    """
    Initialize the submission database if it doesn't yet exist.
    """
    input:
        db_config_path=cfg_path("config_paths/db"),
    output:
        marker=touch("<results>/db/init.marker"),
    log:
        stdout="<logs>/db/init.stdout.log",
        stderr="<logs>/db/init.stderr.log",
    priority: 2
    resources:
        db_handles=1,
    script:
        "../scripts/init_db.sh"


rule scan_inbox:
    """
    Scan an inbox for new submissions.
    """
    input:
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
    output:
        submissions=temp("<results>/scan_inbox/{submitter_id}/{inbox}/submissions.json"),
    log:
        stdout="<logs>/scan_inbox/{submitter_id}/{inbox}.stdout.log",
        stderr="<logs>/scan_inbox/{submitter_id}/{inbox}.stderr.log",
    benchmark:
        "<benchmarks>/scan_inbox/{submitter_id}/{inbox}/benchmark.tsv"
    priority: 2
    params:
        s3_access_key=os.environ.get("GRZ_S3__ACCESS_KEY"),
        s3_secret=os.environ.get("GRZ_S3__SECRET"),
    script:
        "../scripts/scan_inbox.py"


rule filter_single_submission:
    """
    Takes a JSON file from scan_inbox and filters it to contain only the single submission specified by the wildcards.
    """
    input:
        submissions=rules.scan_inbox.output.submissions,
    output:
        filtered_submission=temp(
            "<results>/scan_inbox/filtered/{submitter_id}/{inbox}/{submission_id}.json"
        ),
    log:
        stderr="<logs>/scan_inbox/filtered/{submitter_id}/{inbox}/{submission_id}.stderr.log",
    priority: 1
    params:
        submission_id=lambda wildcards: wildcards.submission_id,
    shell:
        """(
        jq --arg sid "{params.submission_id}" '[.[] | select(.submission_id == $sid)]' {input.submissions} > {output.filtered_submission}
        ) >{log.stderr} 2>&1
        """


rule sync_db_from_inbox:
    """
    Scan an inbox and register all new/updated submissions in the database.
    Sets state to 'uploading' or 'uploaded' mirroring S3 status, and transitions
    'uploading' → 'uploaded' when an upload completes.
    Runs once per inbox per workflow cycle (temp output is cleaned up on completion).
    """
    input:
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        db_config_path=cfg_path("config_paths/db"),
        db_initialized=ancient(rules.init_db.output.marker),
    output:
        marker=touch(
            temp("<results>/sync_db_from_inbox/{submitter_id}/{inbox}/synced.marker")
        ),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/sync_db_from_inbox.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/sync_db_from_inbox.stderr.log",
    priority: 2
    resources:
        db_handles=1,
    script:
        "../scripts/sync_db_from_inbox.sh"


checkpoint select_submissions:
    """
    Determine which submissions are to be processed.
    """
    input:
        scans=expand(
            "<results>/scan_inbox/{submitter}/{inbox}/submissions.json",
            zip,
            submitter=map(itemgetter(0), ALL_INBOX_PAIRS),
            inbox=map(itemgetter(1), ALL_INBOX_PAIRS),
        ),
        db_config_path=cfg_path("config_paths/db"),
    output:
        submissions_batch=temp("<results>/select_submissions/pending_batch.txt"),
    log:
        stdout="<logs>/select_submissions.stdout.log",
        stderr="<logs>/select_submissions.stderr.log",
    priority: 2
    params:
        batch_limit=config["batch"].get("limit", None),
    script:
        "../scripts/select_submissions.py"


rule pending:
    default_target: True
    input:
        all_pending_submissions,
    priority: 3


rule metadata:
    """
    Fetch a submission's metadata from the target S3 bucket.
    Useful for estimating resource usage.
    """
    input:
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        metadata=temp(
            "<results>/{submitter_id}/{inbox}/{submission_id}/metadata/metadata.json"
        ),
        timestamp=temp("<results>/{submitter_id}/{inbox}/{submission_id}/timestamp.txt"),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/download.metadata.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/download.metadata.stderr.log",
    priority: 1
    params:
        s3_access_key=register_s3_access_key,
        s3_secret=register_s3_secret,
        s3_bucket=get_s3_bucket,
        s3_metadata_key=get_s3_metadata_key,
        s3_endpoint_url=get_endpoint_url,
    shell:
        """(
        s5cmd --endpoint-url {params.s3_endpoint_url} cp s3://{params.s3_bucket}/{params.s3_metadata_key} {output.metadata}
        s5cmd --endpoint-url {params.s3_endpoint_url} head s3://{params.s3_bucket}/{params.s3_metadata_key} | jq -r '.last_modified | fromdate | strftime("%Y-%m-%d")' > {output.timestamp}
        ) >{log.stdout} 2> {log.stderr}"""


rule download:
    """
    Download a submission from S3 to the local filesystem.
    """
    input:
        db_synced_marker=rules.sync_db_from_inbox.output.marker,
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        metadata=rules.metadata.output.metadata,
        db_config_path=cfg_path("config_paths/db"),
    output:
        base_dir=temp(
            directory("<results>/{submitter_id}/{inbox}/{submission_id}/downloaded")
        ),
        metadata_dir=temp(
            directory(
                "<results>/{submitter_id}/{inbox}/{submission_id}/downloaded/metadata"
            )
        ),
        encrypted_files_dir=temp(
            directory(
                "<results>/{submitter_id}/{inbox}/{submission_id}/downloaded/encrypted_files"
            )
        ),
        progress_log="<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_download.cjson",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/download.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/download.stderr.log",
    benchmark:
        "<benchmarks>/download/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    priority: 1
    resources:
        disk=estimate_download_size,
        runtime=estimate_download_runtime,
        db_handles=1,
    params:
        s3_access_key=os.environ.get("GRZ_S3__ACCESS_KEY"),
        s3_secret=os.environ.get("GRZ_S3__SECRET"),
    script:
        "../scripts/download.sh"


rule decrypt:
    """
    Decrypt a downloaded submission.
    """
    input:
        base_dir=rules.download.output.base_dir,
        encrypted_files_dir=rules.download.output.encrypted_files_dir,
        metadata_dir=rules.download.output.metadata_dir,
        download_progress_log=rules.download.output.progress_log,
        metadata=rules.metadata.output.metadata,
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        base_dir=temp(
            directory("<results>/{submitter_id}/{inbox}/{submission_id}/decrypted")
        ),
        files_dir=temp(
            directory(
                "<results>/{submitter_id}/{inbox}/{submission_id}/decrypted/files"
            )
        ),
        progress_log="<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_decrypt.cjson",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/decrypt.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/decrypt.stderr.log",
    benchmark:
        "<benchmarks>/decrypt/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    priority: 1
    resources:
        disk=estimate_decrypt_size,
        runtime=estimate_decrypt_runtime,
        db_handles=1,
    params:
        grz_private_key_passphrase=os.environ.get(
            "GRZ_KEYS__GRZ_PRIVATE_KEY_PASSPHRASE"
        ),
    script:
        "../scripts/decrypt.sh"


checkpoint validate:
    """
    Validate a decrypted submission.
    """
    input:
        metadata=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
            rules.metadata.output.metadata,
        ),
        base_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
            rules.decrypt.output.base_dir,
        ),
        files_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
            rules.decrypt.output.files_dir,
        ),
        decrypt_progress_log=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
            rules.decrypt.output.progress_log,
        ),
        inbox_config_path=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
            cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        ),
        db_config_path=cfg_path("config_paths/db"),
    output:
        checksum_log=temp(
            "<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_checksum_validation.cjson"
        ),
        seq_data_log="<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_sequencing_data_validation.cjson",
        validation_flag="<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
        validation_errors="<results>/{submitter_id}/{inbox}/{submission_id}/validation_errors.txt",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/validate.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/validate.stderr.log",
    benchmark:
        "<benchmarks>/validate/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    priority: 2
    threads: 4
    resources:
        runtime=estimate_validate_runtime,
        db_handles=1,
    script:
        "../scripts/validate.sh"


rule consent:
    """
    Check if a submission is research consented.
    """
    input:
        metadata=rules.metadata.output.metadata,
    output:
        consent_flag="<results>/{submitter_id}/{inbox}/{submission_id}/consent_flag",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/consent.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/consent.stderr.log",
    priority: 1
    script:
        "../scripts/consent.sh"


rule re_encrypt:
    """
    Re-encrypt a submission using the target public key (depending on whether research consent was given).
    """
    input:
        metadata=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.metadata.output.metadata,
        ),
        files_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.decrypt.output.files_dir,
        ),
        base_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.decrypt.output.base_dir,
        ),
        decrypt_progress_log=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.decrypt.output.progress_log,
        ),
        validation_flag=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.validate.output.validation_flag,
        ),
        validation_checksum_log=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.validate.output.checksum_log,
        ),
        validation_seq_data_log=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files",
            rules.validate.output.seq_data_log,
        ),
        consent_flag=rules.consent.output.consent_flag,
        consented_config_path=cfg_path("config_paths/archive/consented"),
        nonconsented_config_path=cfg_path("config_paths/archive/nonconsented"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        encrypted_files_dir=temp(
            directory(
                "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files"
            )
        ),
        encryption_log="<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_encrypt.cjson",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/re_encrypt.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/re_encrypt.stderr.log",
    benchmark:
        "<benchmarks>/re_encrypt/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    priority: 1
    resources:
        disk=estimate_re_encrypt_size,
        runtime=estimate_encrypt_runtime,
        db_handles=1,
    script:
        "../scripts/re_encrypt.sh"


rule archive:
    """
    Archive a submission to the target s3 bucket (depending on whether research consent was given).
    """
    input:
        metadata=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/archived",
            rules.metadata.output.metadata,
        ),
        re_encrypted_files_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/archived",
            rules.re_encrypt.output.encrypted_files_dir,
        ),
        consent_flag=rules.consent.output.consent_flag,
        progress_logs_to_archive=[
            rules.download.output.progress_log,
            rules.decrypt.output.progress_log,
            rules.validate.output.checksum_log,
            rules.validate.output.seq_data_log,
            rules.re_encrypt.output.encryption_log,
        ],
        consented_config_path=cfg_path("config_paths/archive/consented"),
        nonconsented_config_path=cfg_path("config_paths/archive/nonconsented"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        marker=touch("<results>/{submitter_id}/{inbox}/{submission_id}/archived"),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/archive.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/archive.stderr.log",
    benchmark:
        "<benchmarks>/archive/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    priority: 2
    resources:
        runtime=estimate_archive_runtime,
        db_handles=1,
    script:
        "../scripts/archive.sh"


rule generate_pruefbericht:
    """
    Generate a Prüfbericht for a submission.
    """
    input:
        metadata=rules.metadata.output.metadata,
        timestamp=rules.metadata.output.timestamp,
        validation_flag=rules.validate.output.validation_flag,
        archived_marker=rules.archive.output.marker,
        db_config_path=cfg_path("config_paths/db"),
    output:
        pruefbericht="<results>/{submitter_id}/{inbox}/{submission_id}/pruefbericht.json",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/generate_pruefbericht.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/generate_pruefbericht.stderr.log",
    priority: 2
    script:
        "../scripts/generate_pruefbericht.sh"


rule submit_pruefbericht:
    """
    Report a Prüfbericht to BfArM.
    """
    input:
        pruefbericht=rules.generate_pruefbericht.output.pruefbericht,
        pruefbericht_config_path=cfg_path("config_paths/pruefbericht"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        answer="<results>/{submitter_id}/{inbox}/{submission_id}/pruefbericht_answer",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/submit_pruefbericht.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/submit_pruefbericht.stderr.log",
    priority: 2
    resources:
        db_handles=1,
    params:
        custom_ca_cert=lambda _: (
            "/workdir/config/cert.pem"
            if os.environ.get("GRZ_PRUEFBERICHT_MOCK", False)
            else ""
        ),
    script:
        "../scripts/submit_pruefbericht.sh"


rule setup_qc_workflow:
    output:
        workflow_dir=directory("<resources>/GRZ_QC_Workflow/"),
        pipeline="<resources>/GRZ_QC_Workflow/main.nf",
    log:
        stdout="<logs>/qc/setup_qc_workflow.stdout.log",
        stderr="<logs>/qc/setup_qc_workflow.stderr.log",
    priority: 0
    params:
        revision=get_qc_workflow_revision,
    shell:
        """
        (
          git clone https://github.com/BfArM-MVH/GRZ_QC_Workflow.git {output.workflow_dir}
          pushd {output.workflow_dir}
          git reset --hard {params.revision}
          popd
        ) >{log.stdout} 2> {log.stderr}
        """


qc_prepare_references_mode = (
    config.get("qc", {}).get("prepare-qc", {}).get("mode", "download")
)
if qc_prepare_references_mode == "download":

    rule download_qc_workflow_references:
        output:
            references_dir=get_qc_workflow_references_directory(),
            launch_dir=directory("<resources>/shared_qc_launchdir"),
        priority: 0
        resources:
            disk="50G",
        shell:
            """
            mkdir -p {output.references_dir}
            pushd {output.references_dir}
            wget --recursive --no-parent -nH --cut-dirs=3 -R "index.html*" https://www.cmm.in.tum.de/public/grz/reference/
            popd
            mkdir -p {output.launch_dir}
            """

else:

    rule prepare_qc_workflow_references:
        input:
            # marked as ancient so as not to re-generate references if not needed
            # if you need to re-generate references, you can delete the references directory manually
            workflow_dir=ancient(rules.setup_qc_workflow.output.workflow_dir),
            pipeline=ancient(rules.setup_qc_workflow.output.pipeline),
            custom_configs=lambda wc: config.get("qc", {})
            .get("prepare-qc", {})
            .get("configs", []),
        output:
            references_dir=get_qc_workflow_references_directory(),
            launch_dir=directory("<resources>/shared_qc_launchdir"),
            work_dir=temp(directory("<resources>/prepare_qc_workflow/work")),
        log:
            stdout="<logs>/qc/prepare_qc_workflow_references.stdout.log",
            stderr="<logs>/qc/prepare_qc_workflow_references.stderr.log",
        benchmark:
            "<benchmarks>/prepare_qc_workflow_references/benchmark.tsv"
        priority: 0
        handover: True
        threads: 1
        resources:
            mem="90G",
            disk="60G",
            runtime="3h",
            tmpdir=get_nextflow_tmpdir,
        params:
            profiles=get_prepare_qc_nextflow_profiles,
            configs=get_prepare_qc_nextflow_configs,
            extra=get_prepare_qc_nextflow_extra_params,
        script:
            "../scripts/prepare_qc.sh"


rule qc:
    """
    Perform QC on a submission using the QC nextflow pipeline.
    """
    input:
        metadata=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/qc/success.marker",
            rules.metadata.output.metadata,
        ),
        base_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/qc/success.marker",
            rules.decrypt.output.base_dir,
        ),
        files_dir=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/qc/success.marker",
            rules.decrypt.output.files_dir,
        ),
        validation_flag=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/qc/success.marker",
            rules.validate.output.validation_flag,
        ),
        db_config_path=cfg_path("config_paths/db"),
        workflow_dir=ancient(rules.setup_qc_workflow.output.workflow_dir),
        pipeline=ancient(rules.setup_qc_workflow.output.pipeline),
        reference_path=get_qc_workflow_references_directory(),
        custom_configs=lambda wc: config.get("qc", {})
        .get("run-qc", {})
        .get("configs", []),
    output:
        out_dir=temp(
            directory("<results>/{submitter_id}/{inbox}/{submission_id}/qc/out")
        ),
        work_dir=temp(
            directory("<results>/{submitter_id}/{inbox}/{submission_id}/qc/work")
        ),
        marker=touch(
            "<results>/{submitter_id}/{inbox}/{submission_id}/qc/success.marker"
        ),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/qc.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/qc.stderr.log",
    benchmark:
        "<benchmarks>/qc/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    priority: 0
    handover: True
    resources:
        runtime=estimate_qc_runtime,
        mem=estimate_qc_memory,
        disk=estimate_qc_disk,
        tmpdir=get_nextflow_tmpdir,
        db_handles=1,
    params:
        profiles=get_run_qc_nextflow_profiles,
        configs=get_run_qc_nextflow_configs,
        extra=get_run_qc_nextflow_extra_params,
    script:
        "../scripts/run_qc.sh"


rule process_qc_results:
    input:
        qc_results=rules.qc.output.out_dir,
        qc_done=rules.qc.output.marker,
        db_config_path=cfg_path("config_paths/db"),
    output:
        marker=touch(
            "<results>/{submitter_id}/{inbox}/{submission_id}/qc/processed.marker"
        ),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/process_qc_results/stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/process_qc_results/stderr.log",
    priority: 1
    params:
        report_csv=lambda wildcards, input: Path(input.qc_results) / "report.csv",
        qc_workflow_version=lambda wildcards: config["qc"]["revision"],
    script:
        "../scripts/process_qc_results.sh"


rule clean:
    """
    Clean up a submission after reporting the Prüfbericht / performing QC, or after a failed validation.
    This entails removing the data from the inbox S3 bucket.
    It will, however, leave the respective metadata.json in place to prevent an LE from uploading the same
    submission_id multiple times.
    """
    input:
        db_config_path=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
            cfg_path("config_paths/db"),
        ),
        ready_marker=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
            get_cleanup_prerequisite,
        ),
        validation_flag=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
            rules.validate.output.validation_flag,
        ),
        validation_errors=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
            rules.validate.output.validation_errors,
        ),
        consent_flag=payload(
            "<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
            rules.consent.output.consent_flag,
        ),
        inbox_config_path=anchor(
            "<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
            cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        ),
    output:
        clean_results="<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}.stderr.log",
    benchmark:
        "<benchmarks>/clean/{submitter_id}/{inbox}/{submission_id}/{qc_status}/benchmark.tsv"
    priority: 2
    resources:
        db_handles=1,
    params:
        mode=config.get("auto-cleanup", "none"),
    script:
        "../scripts/clean.sh"


rule finalize_fail:
    """
    Target rule for each _failed_ submission.
    """
    input:
        unpack(get_failed_finalize_inputs),
    output:
        target=touch(
            "<results>/{submitter_id}/{inbox}/{submission_id}/target/failed_validation.{qc_status}"
        ),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/finalize_fail/{qc_status}.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/finalize_fail/{qc_status}.stderr.log",
    priority: 1
    resources:
        db_handles=1,
    shell:
        """
        (
        echo "Submission {wildcards.submission_id} failed validation."
        grzctl db --config-file {input.db_config_path} submission update --ignore-error-state {wildcards.submission_id} error --data '{{"reason": "validation failed"}}'
        echo "Submission {wildcards.submission_id} processing finished due to validation failure." > {output.target}
        ) > {log.stdout} 2> {log.stderr}
        """


rule finalize_success:
    """
    Target rule for each _valid_ submission.
    This explicitly states what needs to happen (which inputs are required)
    in order for a submission to be considered "finished".
    Will update the submission status in the database accordingly.
    """
    input:
        unpack(get_successful_finalize_inputs),
    output:
        target="<results>/{submitter_id}/{inbox}/{submission_id}/target/{qc_status}",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/finalize_success/{qc_status}.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/finalize_success/{qc_status}.stderr.log",
    priority: 1
    resources:
        db_handles=1,
    shell:
        """
        (
            CLEANED=$(cat {input.clean_results})
            if [[ "$CLEANED" == "true" ]]; then
                echo "Submission {wildcards.submission_id} successfully finalized ({wildcards.qc_status})." > {output.target}
            else
                echo "Submission {wildcards.submission_id} failed to finalize ({wildcards.qc_status})." > {output.target}
                exit 1
            fi
            grzctl db --config-file {input.db_config_path} submission update {wildcards.submission_id} finished
        ) > {log.stdout} 2> {log.stderr}
        """


rule process_submission:
    input:
        get_final_submission_target,
        db_config_path=cfg_path("config_paths/db"),
    output:
        touch("<results>/{submitter_id}/{inbox}/{submission_id}/processed"),
    priority: 3
