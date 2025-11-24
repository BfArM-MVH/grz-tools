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
    benchmark:
        "<benchmarks>/scan_inbox/{submitter_id}/{inbox}/benchmark.tsv"
    params:
        s3_access_key=os.environ.get("GRZ_S3__ACCESS_KEY"),
        s3_secret=os.environ.get("GRZ_S3__SECRET"),
    log:
        stdout="<logs>/scan_inbox/{submitter_id}/{inbox}.stdout.log",
        stderr="<logs>/scan_inbox/{submitter_id}/{inbox}.stderr.log",
    priority: 2
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
    params:
        submission_id=lambda wildcards: wildcards.submission_id,
    log:
        stderr="<logs>/scan_inbox/filtered/{submitter_id}/{inbox}/{submission_id}.stderr.log",
    priority: 1
    shell:
        """(
        jq --arg sid "{params.submission_id}" '[.[] | select(.submission_id == $sid)]' {input.submissions} > {output.filtered_submission}
        ) >{log.stderr} 2>&1
        """


rule check_and_register:
    """
    Checks S3 for a 'complete' submission status and ensures the
    submission is registered in the DB with the state 'uploaded'.
    """
    input:
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        db_config_path=cfg_path("config_paths/db"),
        db_initialized=rules.init_db.output.marker,
    output:
        marker=touch(
            "<results>/{submitter_id}/{inbox}/{submission_id}/registered.marker"
        ),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/check_and_register.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/check_and_register.stderr.log",
    priority: 2
    script:
        "../scripts/check_and_register.py"


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
    params:
        batch_limit=config["batch"].get("limit", None),
    log:
        stdout="<logs>/select_submissions.stdout.log",
        stderr="<logs>/select_submissions.stderr.log",
    priority: 2
    script:
        "../scripts/select_submissions.py"


rule pending:
    input:
        all_pending_submissions,
    default_target: True
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
    params:
        s3_access_key=register_s3_access_key,
        s3_secret=register_s3_secret,
        s3_bucket=get_s3_bucket,
        s3_metadata_key=get_s3_metadata_key,
        s3_endpoint_url=get_endpoint_url,
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/download.metadata.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/download.metadata.stderr.log",
    priority: 1
    shell:
        """(
        s5cmd --endpoint-url {params.s3_endpoint_url} cp s3://{params.s3_bucket}/{params.s3_metadata_key} {output.metadata}
        ) >{log.stdout} 2> {log.stderr}"""


rule download:
    """
    Download a submission from S3 to the local filesystem.
    """
    input:
        registered_marker=rules.check_and_register.output.marker,
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
    benchmark:
        "<benchmarks>/download/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    params:
        s3_access_key=os.environ.get("GRZ_S3__ACCESS_KEY"),
        s3_secret=os.environ.get("GRZ_S3__SECRET"),
    resources:
        disk=estimate_download_size,
        runtime=estimate_download_runtime,
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/download.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/download.stderr.log",
    priority: 1
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
    benchmark:
        "<benchmarks>/decrypt/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    params:
        grz_private_key_passphrase=os.environ.get(
            "GRZ_KEYS__GRZ_PRIVATE_KEY_PASSPHRASE"
        ),
    resources:
        disk=estimate_decrypt_size,
        runtime=estimate_decrypt_runtime,
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/decrypt.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/decrypt.stderr.log",
    priority: 1
    script:
        "../scripts/decrypt.sh"


checkpoint validate:
    """
    Validate a decrypted submission.

    This is a checkpoint because the downstream workflow (success path vs. failure path) depends on its output.
    """
    input:
        base_dir=rules.decrypt.output.base_dir,
        files_dir=rules.decrypt.output.files_dir,
        decrypt_progress_log=rules.decrypt.output.progress_log,
        metadata=rules.metadata.output.metadata,
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        checksum_log=temp(
            "<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_checksum_validation.cjson"
        ),
        seq_data_log="<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_sequencing_data_validation.cjson",
        validation_flag="<results>/{submitter_id}/{inbox}/{submission_id}/validation_flag",
        validation_errors="<results>/{submitter_id}/{inbox}/{submission_id}/validation_errors.txt",
    benchmark:
        "<benchmarks>/validate/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/validate.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/validate.stderr.log",
    threads: 4
    resources:
        runtime=estimate_validate_runtime,
    priority: 2
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
        metadata=rules.metadata.output.metadata,
        base_dir=rules.decrypt.output.base_dir,
        files_dir=rules.decrypt.output.files_dir,
        decrypt_progress_log=rules.decrypt.output.progress_log,
        consent_flag=rules.consent.output.consent_flag,
        validation_flag=rules.validate.output.validation_flag,
        consented_config_path=cfg_path("config_paths/archive/consented"),
        nonconsented_config_path=cfg_path("config_paths/archive/nonconsented"),
        db_config_path=cfg_path("config_paths/db"),
        validation_checksum_log=rules.validate.output.checksum_log,
        validation_seq_data_log=rules.validate.output.seq_data_log,
    output:
        encrypted_files_dir=temp(
            directory(
                "<results>/{submitter_id}/{inbox}/{submission_id}/re-encrypted/encrypted_files"
            )
        ),
        encryption_log="<results>/{submitter_id}/{inbox}/{submission_id}/progress_logs/progress_encrypt.cjson",
    benchmark:
        "<benchmarks>/re_encrypt/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    resources:
        disk=estimate_re_encrypt_size,
        runtime=estimate_encrypt_runtime,
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/re_encrypt.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/re_encrypt.stderr.log",
    priority: 1
    script:
        "../scripts/re_encrypt.sh"


rule archive:
    """
    Archive a submission to the target s3 bucket (depending on whether research consent was given).
    """
    input:
        metadata=rules.metadata.output.metadata,
        re_encrypted_files_dir=rules.re_encrypt.output.encrypted_files_dir,
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
    benchmark:
        "<benchmarks>/archive/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    resources:
        runtime=estimate_archive_runtime,
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/archive.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/archive.stderr.log",
    priority: 2
    script:
        "../scripts/archive.sh"


rule generate_pruefbericht:
    """
    Generate a Prüfbericht for a submission.
    """
    input:
        metadata=rules.metadata.output.metadata,
        validation_flag=get_validation_state,
        archived_marker=ancient(rules.archive.output.marker),
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
    params:
        custom_ca_cert=lambda _: (
            "/workdir/config/cert.pem"
            if os.environ.get("GRZ_PRUEFBERICHT_MOCK", False)
            else ""
        ),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/submit_pruefbericht.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/submit_pruefbericht.stderr.log",
    priority: 2
    script:
        "../scripts/submit_pruefbericht.sh"


rule setup_qc_workflow:
    output:
        workflow_dir=directory("<resources>/GRZ_QC_Workflow/"),
        pipeline="<resources>/GRZ_QC_Workflow/main.nf",
    params:
        revision=get_qc_workflow_revision,
    log:
        stdout="<logs>/qc/setup_qc_workflow.stdout.log",
        stderr="<logs>/qc/setup_qc_workflow.stderr.log",
    priority: 0
    shell:
        """
        (
          git clone --branch {params.revision} --depth 1 https://github.com/BfArM-MVH/GRZ_QC_Workflow.git {output.workflow_dir}
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
        resources:
            disk="50G",
        priority: 0
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
        benchmark:
            "<benchmarks>/prepare_qc_workflow_references/benchmark.tsv"
        params:
            profiles=get_prepare_qc_nextflow_profiles,
            configs=get_prepare_qc_nextflow_configs,
            extra=get_prepare_qc_nextflow_extra_params,
        threads: 1
        resources:
            mem="90G",
            disk="60G",
            runtime="3h",
            tmpdir=get_nextflow_tmpdir,
        log:
            stdout="<logs>/qc/prepare_qc_workflow_references.stdout.log",
            stderr="<logs>/qc/prepare_qc_workflow_references.stderr.log",
        priority: 0
        script:
            "../scripts/prepare_qc.sh"


rule qc:
    """
    Perform QC on a submission using the QC nextflow pipeline.
    """
    input:
        base_dir=rules.decrypt.output.base_dir,
        files_dir=rules.decrypt.output.files_dir,
        metadata=rules.metadata.output.metadata,
        validation_flag=get_validation_state,
        db_config_path=cfg_path("config_paths/db"),
        workflow_dir=ancient(rules.setup_qc_workflow.output.workflow_dir),
        pipeline=ancient(rules.setup_qc_workflow.output.pipeline),
        launch_dir="<resources>/shared_qc_launchdir",
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
    benchmark:
        "<benchmarks>/qc/{submitter_id}/{inbox}/{submission_id}/benchmark.tsv"
    resources:
        runtime=estimate_qc_runtime,
        mem=estimate_qc_memory,
        disk=estimate_qc_disk,
        tmpdir=get_nextflow_tmpdir,
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/qc.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/qc.stderr.log",
    params:
        profiles=get_run_qc_nextflow_profiles,
        configs=get_run_qc_nextflow_configs,
        extra=get_run_qc_nextflow_extra_params,
    priority: 0
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
    params:
        report_csv=lambda wildcards, input: Path(input.qc_results) / "report.csv",
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/process_qc_results/stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/process_qc_results/stderr.log",
    priority: 1
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
        ready_marker=get_cleanup_prerequisite,
        validation_flag=get_validation_state,
        validation_errors=rules.validate.output.validation_errors,
        consent_flag=rules.consent.output.consent_flag,
        inbox_config_path=cfg_path("config_paths/inbox/{submitter_id}/{inbox}"),
        db_config_path=cfg_path("config_paths/db"),
    output:
        clean_results="<results>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}",
    benchmark:
        "<benchmarks>/clean/{submitter_id}/{inbox}/{submission_id}/{qc_status}/benchmark.tsv"
    params:
        mode=config.get("auto-cleanup", "none"),
    log:
        stdout="<logs>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}.stdout.log",
        stderr="<logs>/{submitter_id}/{inbox}/{submission_id}/clean/{qc_status}.stderr.log",
    priority: 2
    script:
        "../scripts/clean.sh"


# rule taetigkeitsbericht:
#     """
#     Update the Taetigkeitsbericht with information from the submission DB.
#     """
#     output:
#         taetigkeitsbericht="<results>/{submitter_id}/{inbox}/{submission_id}/taetigkeitsbericht.pdf",
#     params:
#         db_config_path=cfg_path("config_paths/db"),
#     log:
#         "<logs>/{submitter_id}/{inbox}/{submission_id}/taetigkeitsbericht.log",
#     shell:
#         """
#         # (generate taetigkeitsbericht using submission DB via {input.db_config_path}) 2> {log}
#         echo "TODO" > {output.taetigkeitsbericht} 2> {log}
#         """


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
    shell:
        """
        (
        echo "Submission {wildcards.submission_id} failed validation."
        grzctl db --config-file {input.db_config_path} submission update {wildcards.submission_id} error --data '{{"reason": "validation failed"}}'
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
