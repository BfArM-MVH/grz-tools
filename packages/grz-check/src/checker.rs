use crate::checks::bam::BamCheckJob;
use crate::checks::common;
use crate::checks::fastq::{PairedFastqJob, SingleFastqJob};
use crate::checks::raw::RawJob;
use crate::checks::{bam, fastq, raw};
use anyhow::Context;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use serde::Serialize;
use sha2::{Digest, Sha256};
use std::error::Error as StdError;
use std::fmt;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::sync::{Arc, Mutex};

#[derive(Debug, Copy, Clone, Eq, PartialEq, Serialize)]
pub struct Stats {
    pub num_records: u64,
    pub read_length: Option<usize>,
}

#[derive(Debug, Serialize, Clone)]
pub struct FileReport {
    pub path: PathBuf,
    pub stats: Option<Stats>,
    pub sha256: Option<String>,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

impl FileReport {
    pub fn new(
        path: &Path,
        stats: Option<Stats>,
        errors: Vec<String>,
        warnings: Vec<String>,
    ) -> Self {
        Self {
            path: path.to_path_buf(),
            stats,
            sha256: None,
            errors,
            warnings,
        }
    }

    pub fn new_with_error(path: &Path, error: String) -> Self {
        Self {
            path: path.to_path_buf(),
            stats: None,
            sha256: None,
            errors: vec![error],
            warnings: vec![],
        }
    }

    pub fn with_sha256(mut self, sha256: Option<String>) -> Self {
        self.sha256 = sha256;
        self
    }

    pub fn is_ok(&self) -> bool {
        self.errors.is_empty()
    }
}

#[derive(Debug, Serialize, Clone)]
pub struct PairReport {
    pub fq1_report: FileReport,
    pub fq2_report: FileReport,
    pub pair_errors: Vec<String>,
}

impl PairReport {
    fn is_ok(&self) -> bool {
        self.fq1_report.is_ok() && self.fq2_report.is_ok() && self.pair_errors.is_empty()
    }
}

#[derive(Debug)]
pub enum Job {
    SingleFastq(SingleFastqJob),
    PairedFastq(PairedFastqJob),
    Bam(BamCheckJob),
    Raw(RawJob),
}

#[derive(Debug)]
enum CheckResult {
    PairedFastq(PairReport),
    SingleFastq(FileReport),
    Bam(FileReport),
    Raw(FileReport),
}

impl CheckResult {
    fn is_error(&self) -> bool {
        match self {
            CheckResult::PairedFastq(r) => !r.is_ok(),
            CheckResult::SingleFastq(r) => !r.is_ok(),
            CheckResult::Bam(r) => !r.is_ok(),
            CheckResult::Raw(r) => !r.is_ok(),
        }
    }
    fn primary_path(&self) -> &Path {
        match self {
            CheckResult::PairedFastq(r) => &r.fq1_report.path,
            CheckResult::SingleFastq(r) => &r.path,
            CheckResult::Bam(r) => &r.path,
            CheckResult::Raw(r) => &r.path,
        }
    }
}

#[allow(clippy::large_enum_variant)]
#[derive(Debug)]
enum StopReason {
    Error(CheckResult),
    Interrupted,
}

#[derive(Debug)]
pub struct EarlyExitError(StopReason);

impl fmt::Display for EarlyExitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.0 {
            StopReason::Error(_) => write!(f, "A validation error occurred, exiting."),
            StopReason::Interrupted => write!(f, "Operation was interrupted by the user."),
        }
    }
}

impl StdError for EarlyExitError {}

fn filename(path: impl AsRef<Path>) -> String {
    path.as_ref()
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_string()
}

fn process_job(
    (m, main_pb, style): &mut (MultiProgress, ProgressBar, ProgressStyle),
    job: Job,
) -> CheckResult {
    match job {
        Job::SingleFastq(job) => {
            let pb = m.add(ProgressBar::new(job.size));
            pb.set_style(style.clone());
            pb.set_prefix("FASTQ");
            let report = fastq::check_single_fastq(&job.path, job.length_check, &pb, main_pb);
            if report.is_ok() {
                pb.finish_with_message(format!("✓ OK    {}", filename(&job.path)));
            } else {
                pb.abandon_with_message(format!("✗ ERROR {}", filename(&job.path)));
            }
            CheckResult::SingleFastq(report)
        }
        Job::PairedFastq(job) => {
            let fq1_pb = m.add(ProgressBar::new(job.fq1_size));
            fq1_pb.set_style(style.clone());
            fq1_pb.set_prefix("FASTQ R1");

            let fq2_pb = m.add(ProgressBar::new(job.fq2_size));
            fq2_pb.set_style(style.clone());
            fq2_pb.set_prefix("FASTQ R2");

            let fq1_setup = common::setup_file_reader(&job.fq1_path, &fq1_pb, main_pb, true);
            let fq2_setup = common::setup_file_reader(&job.fq2_path, &fq2_pb, main_pb, true);

            let report = match (fq1_setup, fq2_setup) {
                (Ok((reader1, hasher1)), Ok((reader2, hasher2))) => {
                    let (fq1_outcome, fq2_outcome, pair_errors) =
                        match fastq::process_paired_readers(
                            reader1,
                            reader2,
                            job.fq1_length_check,
                            job.fq2_length_check,
                        ) {
                            Ok(result) => result,
                            Err(e) => {
                                let outcome1 = common::CheckOutcome {
                                    errors: vec![e.clone()],
                                    ..Default::default()
                                };
                                let outcome2 = common::CheckOutcome {
                                    errors: vec![e],
                                    ..Default::default()
                                };
                                return CheckResult::PairedFastq(PairReport {
                                    fq1_report: FileReport::new(
                                        &job.fq1_path,
                                        None,
                                        outcome1.errors,
                                        outcome1.warnings,
                                    ),
                                    fq2_report: FileReport::new(
                                        &job.fq2_path,
                                        None,
                                        outcome2.errors,
                                        outcome2.warnings,
                                    ),
                                    pair_errors: vec![
                                        "Parsing error during paired fastq check.".to_string(),
                                    ],
                                });
                            }
                        };

                    let finalize = |hasher: Arc<Mutex<Sha256>>| match Arc::try_unwrap(hasher) {
                        Ok(mutex) => Some(format!("{:x}", mutex.into_inner().unwrap().finalize())),
                        Err(_) => None,
                    };
                    let cs1 = finalize(hasher1);
                    let cs2 = finalize(hasher2);

                    let fq1_report = FileReport::new(
                        &job.fq1_path,
                        fq1_outcome.stats,
                        fq1_outcome.errors,
                        fq1_outcome.warnings,
                    )
                    .with_sha256(cs1);
                    let fq2_report = FileReport::new(
                        &job.fq2_path,
                        fq2_outcome.stats,
                        fq2_outcome.errors,
                        fq2_outcome.warnings,
                    )
                    .with_sha256(cs2);

                    PairReport {
                        fq1_report,
                        fq2_report,
                        pair_errors,
                    }
                }
                (Err(e1), Ok((_r2, _h2))) => {
                    let fq1_report = FileReport::new_with_error(&job.fq1_path, e1.to_string());
                    let fq2_report = FileReport::new(
                        &job.fq2_path,
                        None,
                        vec![format!(
                            "R1 ({:?}) failed to parse; check aborted.",
                            &job.fq1_path
                        )],
                        vec![],
                    );
                    PairReport {
                        fq1_report,
                        fq2_report,
                        pair_errors: vec![],
                    }
                }
                (Ok((_r1, _h1)), Err(e2)) => {
                    let fq1_report = FileReport::new(
                        &job.fq1_path,
                        None,
                        vec![format!(
                            "R2 ({:?}) failed to parse; check aborted.",
                            &job.fq2_path
                        )],
                        vec![],
                    );
                    let fq2_report = FileReport::new_with_error(&job.fq2_path, e2.to_string());
                    PairReport {
                        fq1_report,
                        fq2_report,
                        pair_errors: vec![],
                    }
                }
                (Err(e1), Err(e2)) => {
                    let fq1_report = FileReport::new_with_error(&job.fq1_path, e1.to_string());
                    let fq2_report = FileReport::new_with_error(&job.fq2_path, e2.to_string());
                    PairReport {
                        fq1_report,
                        fq2_report,
                        pair_errors: vec![],
                    }
                }
            };

            let fq1_filename = filename(&job.fq1_path);
            let fq2_filename = filename(&job.fq2_path);
            finish_pb(fq1_pb, fq1_filename, &report.fq1_report);
            finish_pb(fq2_pb, fq2_filename, &report.fq1_report);

            CheckResult::PairedFastq(report)
        }
        Job::Bam(job) => {
            let pb = m.add(ProgressBar::new(job.size));
            pb.set_style(style.clone());
            pb.set_prefix("BAM");
            let filename = filename(&job.path);
            let report = bam::check_bam(&job.path, &pb, main_pb);
            finish_pb(pb, filename, &report);
            CheckResult::Bam(report)
        }
        Job::Raw(job) => {
            let pb = m.add(ProgressBar::new(job.size));
            pb.set_style(style.clone());
            pb.set_prefix("OTHER");
            let report = raw::check_raw(&job.path, &pb, main_pb);
            let filename = filename(&job.path);
            finish_pb(pb, filename, &report);
            CheckResult::Raw(report)
        }
    }
}

fn finish_pb(pb: ProgressBar, filename: String, report: &FileReport) {
    if report.is_ok() {
        pb.finish_with_message(format!("✓ OK    {filename}"));
    } else {
        pb.abandon_with_message(format!("✗ ERROR {filename}"));
    }
}

#[allow(clippy::result_large_err)]
fn process_jobs(
    jobs: Vec<Job>,
    continue_on_error: bool,
    shutdown_flag: Arc<AtomicBool>,
    mpb: MultiProgress,
    main_pb: ProgressBar,
    file_style: ProgressStyle,
    writer: Arc<Mutex<BufWriter<fs::File>>>,
) -> Result<(), EarlyExitError> {
    if continue_on_error {
        let num_failed_jobs = Arc::new(AtomicUsize::new(0));

        jobs.into_par_iter().for_each_with(
            (
                mpb,
                main_pb.clone(),
                file_style,
                writer,
                num_failed_jobs.clone(),
            ),
            |(mpb, main_pb, style, writer, num_failed), job| {
                if shutdown_flag.load(Ordering::Relaxed) {
                    return;
                }

                let report = process_job(&mut (mpb.clone(), main_pb.clone(), style.clone()), job);

                if report.is_error() {
                    num_failed.fetch_add(1, Ordering::SeqCst);
                }

                let mut writer_guard = writer.lock().unwrap();
                if let Err(e) = write_jsonl_report_entry(&report, &mut *writer_guard) {
                    eprintln!(
                        "Failed to write report line for {:?}: {}",
                        report.primary_path(),
                        e
                    );
                }
            },
        );

        let final_fail_count = num_failed_jobs.load(Ordering::SeqCst);
        if shutdown_flag.load(Ordering::SeqCst) {
            main_pb.abandon_with_message("✗ Operation cancelled by user.");
        } else if final_fail_count > 0 {
            main_pb.abandon_with_message(format!(
                "✗ Processing complete. {final_fail_count} pairs/files failed."
            ));
        } else {
            main_pb.finish_with_message("✓ All checks passed!");
        }

        Ok(())
    } else {
        jobs.into_par_iter().try_for_each_with(
            (mpb, main_pb, file_style, writer),
            |(mpb, main_pb, style, writer), job| {
                if shutdown_flag.load(Ordering::Relaxed) {
                    return Err(EarlyExitError(StopReason::Interrupted));
                }
                let report = process_job(&mut (mpb.clone(), main_pb.clone(), style.clone()), job);

                let mut writer_guard = writer.lock().unwrap();
                if let Err(e) = write_jsonl_report_entry(&report, &mut *writer_guard) {
                    eprintln!(
                        "Failed to write report line for {:?}: {}",
                        report.primary_path(),
                        e
                    );
                }
                writer_guard.flush().ok();
                drop(writer_guard);

                if report.is_error() {
                    Err(EarlyExitError(StopReason::Error(report)))
                } else {
                    Ok(())
                }
            },
        )
    }
}

use std::sync::{LazyLock, Once};
static SHUTDOWN_FLAG: LazyLock<Arc<AtomicBool>> =
    LazyLock::new(|| Arc::new(AtomicBool::new(false)));

static SET_HANDLER_ONCE: Once = Once::new();

fn setup_signal_handler() -> anyhow::Result<()> {
    let mut result = Ok(());

    SET_HANDLER_ONCE.call_once(|| {
        let handler_flag = SHUTDOWN_FLAG.clone();
        let set_handler_result = ctrlc::set_handler(move || {
            if handler_flag.swap(true, Ordering::SeqCst) {
                eprintln!("\nSecond interrupt received, exiting immediately.");
                std::process::exit(130);
            }
            eprintln!("\nCtrl+C received, shutting down gracefully…");
        });

        if let Err(e) = set_handler_result {
            result = Err(e).context("Error setting Ctrl-C handler");
        }
    });

    result
}

pub fn run_check(
    jobs: Vec<Job>,
    total_bytes: u64,
    output: &Path,
    continue_on_error: bool,
    show_progress: Option<bool>,
) -> anyhow::Result<()> {
    setup_signal_handler()?;
    let shutdown_flag = SHUTDOWN_FLAG.clone();

    let mpb = MultiProgress::new();
    match show_progress {
        Some(true) => {
            mpb.set_draw_target(ProgressDrawTarget::stderr());
        }
        Some(false) => {
            mpb.set_draw_target(ProgressDrawTarget::hidden());
        }
        _ => {}
    }

    let file_style = ProgressStyle::with_template(
        "{prefix:8.bold} ▕{bar:50.cyan/blue}▏ {bytes:>10}/{total_bytes:<10} ({bytes_per_sec:>12}, ETA: {eta:>6}) {wide_msg}"
    )?.progress_chars("█▒░");
    let main_pb = mpb.add(ProgressBar::new(total_bytes));
    main_pb.set_style(file_style.clone());
    main_pb.set_prefix("Overall");

    let writer = Arc::new(Mutex::new(BufWriter::new(
        fs::File::create(output)
            .with_context(|| format!("Failed to create report file at {}", output.display()))?,
    )));

    let processing_result = process_jobs(
        jobs,
        continue_on_error,
        shutdown_flag.clone(),
        mpb.clone(),
        main_pb.clone(),
        file_style,
        writer.clone(),
    );

    if let Ok(mutex) = Arc::try_unwrap(writer) {
        if let Ok(mut writer_guard) = mutex.into_inner() {
            writer_guard
                .flush()
                .context("Failed to perform final flush of report file")?;
        }
    }
    mpb.clear()?;

    match processing_result {
        Ok(()) => {
            if shutdown_flag.load(Ordering::Relaxed) {
                main_pb.abandon_with_message("✗ Operation cancelled by user.");
                std::process::exit(130);
            } else if !continue_on_error {
                main_pb.finish_with_message("✓ All checks passed!");
            }
        }
        Err(EarlyExitError(reason)) => match reason {
            StopReason::Error(failed_report) => {
                main_pb.abandon_with_message(format!(
                    "✗ Error in {}. See report: {}",
                    failed_report.primary_path().display(),
                    output.display()
                ));
                anyhow::bail!(
                    "A validation error occurred in {}. Aborting.",
                    failed_report.primary_path().display()
                );
            }
            StopReason::Interrupted => {
                main_pb.abandon_with_message("✗ Operation cancelled by user.");
                std::process::exit(130);
            }
        },
    }

    Ok(())
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "snake_case")]
struct FastqReport<'a> {
    path: &'a Path,
    status: &'a str,
    num_records: Option<u64>,
    read_length: Option<usize>,
    checksum: Option<&'a String>,
    errors: Vec<String>,
    warnings: &'a [String],
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "snake_case")]
struct BamReport<'a> {
    path: &'a Path,
    status: &'a str,
    num_records: Option<u64>,
    checksum: Option<&'a String>,
    errors: &'a [String],
    warnings: &'a [String],
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "snake_case")]
struct RawReport<'a> {
    path: &'a Path,
    status: &'a str,
    checksum: Option<&'a String>,
    errors: &'a [String],
    warnings: &'a [String],
}

#[derive(Debug, Serialize)]
#[serde(tag = "check_type", content = "data", rename_all = "snake_case")]
enum JsonReport<'a> {
    Fastq(FastqReport<'a>),
    Bam(BamReport<'a>),
    Raw(RawReport<'a>),
}

fn write_jsonl_report_entry<W: Write>(result: &CheckResult, writer: &mut W) -> anyhow::Result<()> {
    match result {
        CheckResult::PairedFastq(pair_report) => {
            let is_pair_error = !pair_report.pair_errors.is_empty();

            let r1 = &pair_report.fq1_report;
            let r2 = &pair_report.fq2_report;
            let file_reports = [r1, r2];

            for file_report in file_reports {
                let mut errors = file_report.errors.clone();
                if is_pair_error {
                    errors.extend(pair_report.pair_errors.clone());
                }
                let status = if file_report.is_ok() && !is_pair_error {
                    "OK"
                } else {
                    "ERROR"
                };

                let report = JsonReport::Fastq(FastqReport {
                    path: &file_report.path,
                    status,
                    num_records: file_report.stats.map(|s| s.num_records),
                    read_length: file_report.stats.and_then(|s| s.read_length),
                    checksum: file_report.sha256.as_ref(),
                    errors,
                    warnings: &file_report.warnings,
                });
                serde_json::to_writer(&mut *writer, &report)?;
                writer.write_all(b"\n")?;
            }
        }
        CheckResult::SingleFastq(report) => {
            let json_report = JsonReport::Fastq(FastqReport {
                path: &report.path,
                status: if report.is_ok() { "OK" } else { "ERROR" },
                num_records: report.stats.map(|s| s.num_records),
                read_length: report.stats.and_then(|s| s.read_length),
                checksum: report.sha256.as_ref(),
                errors: report.errors.clone(),
                warnings: &report.warnings,
            });
            serde_json::to_writer(&mut *writer, &json_report)?;
            writer.write_all(b"\n")?;
        }
        CheckResult::Bam(report) => {
            let json_report = JsonReport::Bam(BamReport {
                path: &report.path,
                status: if report.is_ok() { "OK" } else { "ERROR" },
                num_records: report.stats.map(|s| s.num_records),
                checksum: report.sha256.as_ref(),
                errors: &report.errors,
                warnings: &report.warnings,
            });
            serde_json::to_writer(&mut *writer, &json_report)?;
            writer.write_all(b"\n")?;
        }
        CheckResult::Raw(report) => {
            let json_report = JsonReport::Raw(RawReport {
                path: &report.path,
                status: if report.is_ok() { "OK" } else { "ERROR" },
                checksum: report.sha256.as_ref(),
                errors: &report.errors,
                warnings: &report.warnings,
            });
            serde_json::to_writer(&mut *writer, &json_report)?;
            writer.write_all(b"\n")?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Result, anyhow};
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use noodles::bam;

    use crate::checks::fastq::ReadLengthCheck;
    use noodles::sam::alignment::io::Write as SamWrite;
    use noodles::sam::alignment::record::Flags;
    use noodles::sam::alignment::record::cigar::op::{Kind, Op};
    use noodles::sam::alignment::record_buf;
    use noodles::sam::alignment::record_buf::QualityScores;
    use noodles::sam::header::record::value::map::ReadGroup;
    use noodles::sam::{Header, header::record::value::Map};
    use serde::Deserialize;
    use std::io::{BufRead, BufReader, Write};
    use tempfile::tempdir;

    fn create_gzipped_fastq(path: &Path, content: &str) -> Result<()> {
        let file = fs::File::create(path)?;
        let mut writer = GzEncoder::new(file, Compression::default());
        writer.write_all(content.as_bytes())?;
        writer.finish()?;
        Ok(())
    }

    struct TestFiles {
        _tempdir: tempfile::TempDir,
        pub dir: PathBuf,
    }

    impl TestFiles {
        fn new() -> Result<Self> {
            let tempdir = tempdir()?;
            let dir = tempdir.path().to_path_buf();

            // Case 1: Valid pair
            create_gzipped_fastq(
                &dir.join("ok_r1.fastq.gz"),
                "@SEQ1\nACGT\n+\nFFFF\n@SEQ2\nTGCA\n+\nFFFF\n",
            )?;
            create_gzipped_fastq(
                &dir.join("ok_r2.fastq.gz"),
                "@SEQ1\nAAAA\n+\nFFFF\n@SEQ2\nTTTT\n+\nFFFF\n",
            )?;

            // Case 2: Inconsistent length in a single file
            create_gzipped_fastq(
                &dir.join("badlen.fastq.gz"),
                "@SEQ1\nACGT\n+\nFFFF\n@SEQ2\nTCG\n+\nFFF\n",
            )?;

            // Case 3: Mismatched read counts in a pair
            create_gzipped_fastq(
                &dir.join("counts1.fastq.gz"),
                "@SEQ1\nACGT\n+\nFFFF\n@SEQ2\nTGCA\n+\nFFFF\n",
            )?;
            create_gzipped_fastq(&dir.join("counts2.fastq.gz"), "@SEQ1\nACGT\n+\nFFFF\n")?;

            // Case 4: Malformed file (seq len != qual len)
            create_gzipped_fastq(&dir.join("malformed.fastq.gz"), "@SEQ1\nACGT\n+\nFF\n")?;

            // Case 5: Empty file
            create_gzipped_fastq(&dir.join("empty.fastq.gz"), "")?;

            // Case 6: Read length 5 for R2 (for checking different read lengths in the same pair)
            create_gzipped_fastq(
                &dir.join("ok_r2_len5.fastq.gz"),
                "@SEQ1\nAAAAA\n+\nFFFFF\n@SEQ2\nTTTTA\n+\nFFFFF\n",
            )?;

            Ok(Self {
                _tempdir: tempdir,
                dir,
            })
        }
    }

    #[allow(dead_code)]
    #[derive(Deserialize, Debug, Clone)]
    #[serde(rename_all = "snake_case")]
    struct TestFastqReportData {
        path: PathBuf,
        status: String,
        num_records: Option<u64>,
        read_length: Option<usize>,
        checksum: Option<String>,
        errors: Vec<String>,
        warnings: Vec<String>,
    }

    #[allow(dead_code)]
    #[derive(Deserialize, Debug, Clone)]
    #[serde(rename_all = "snake_case")]
    struct TestBamReportData {
        path: PathBuf,
        status: String,
        num_records: Option<u64>,
        checksum: Option<String>,
        errors: Vec<String>,
        warnings: Vec<String>,
    }

    #[allow(dead_code)]
    #[derive(Deserialize, Debug, Clone)]
    #[serde(rename_all = "snake_case")]
    struct TestRawReportData {
        path: PathBuf,
        status: String,
        checksum: Option<String>,
        errors: Vec<String>,
        warnings: Vec<String>,
    }

    #[derive(Deserialize, Debug, Clone)]
    #[serde(tag = "check_type", content = "data", rename_all = "snake_case")]
    enum TestReport {
        Fastq(TestFastqReportData),
        Bam(TestBamReportData),
        Raw(TestRawReportData),
    }

    fn read_jsonl_report(report_path: &Path) -> Result<Vec<TestReport>> {
        let file = fs::File::open(report_path)?;
        let reader = BufReader::new(file);
        reader
            .lines()
            .map(|line| {
                let line = line?;
                serde_json::from_str::<TestReport>(&line).map_err(|e| anyhow!(e))
            })
            .collect()
    }

    #[test]
    fn test_valid_pair_with_different_lengths() -> Result<()> {
        let fixture = TestFiles::new()?;
        let output = fixture.dir.join("report.jsonl");

        let fq1_path = fixture.dir.join("ok_r1.fastq.gz");
        let fq2_path = fixture.dir.join("ok_r2_len5.fastq.gz");
        let fq1_size = fs::metadata(&fq1_path)?.len();
        let fq2_size = fs::metadata(&fq2_path)?.len();
        let total_bytes = fq1_size + fq2_size;

        let jobs = vec![Job::PairedFastq(PairedFastqJob {
            fq1_path,
            fq2_path,
            fq1_length_check: ReadLengthCheck::Fixed(4),
            fq2_length_check: ReadLengthCheck::Fixed(5),
            fq1_size,
            fq2_size,
        })];

        run_check(jobs, total_bytes, &output, false, Some(false))?;

        let mut records = read_jsonl_report(&output)?;
        records.sort_by(|a, b| match (a, b) {
            (TestReport::Fastq(d1), TestReport::Fastq(d2)) => d1.path.cmp(&d2.path),
            _ => panic!("Unexpected report types"),
        });

        assert_eq!(records.len(), 2);
        if let TestReport::Fastq(data) = &records[0] {
            assert!(data.path.ends_with("ok_r1.fastq.gz"));
            assert_eq!(data.status, "OK");
            assert_eq!(data.read_length, Some(4));
        } else {
            panic!("Expected a Fastq report for R1");
        }

        if let TestReport::Fastq(data) = &records[1] {
            assert!(data.path.ends_with("ok_r2_len5.fastq.gz"));
            assert_eq!(data.status, "OK");
            assert_eq!(data.read_length, Some(5));
        } else {
            panic!("Expected a Fastq report for R2");
        }
        Ok(())
    }

    #[test]
    fn test_multiple_inputs_with_continue_on_error() -> Result<()> {
        let fixture = TestFiles::new()?;
        let output = fixture.dir.join("report.jsonl");

        let mut jobs = Vec::new();
        let mut total_bytes = 0;

        let p1f1_path = fixture.dir.join("counts1.fastq.gz");
        let p1f2_path = fixture.dir.join("counts2.fastq.gz");
        let p1f1_size = fs::metadata(&p1f1_path)?.len();
        let p1f2_size = fs::metadata(&p1f2_path)?.len();
        total_bytes += p1f1_size + p1f2_size;
        jobs.push(Job::PairedFastq(PairedFastqJob {
            fq1_path: p1f1_path,
            fq2_path: p1f2_path,
            fq1_length_check: ReadLengthCheck::Auto,
            fq2_length_check: ReadLengthCheck::Auto,
            fq1_size: p1f1_size,
            fq2_size: p1f2_size,
        }));

        let p2f1_path = fixture.dir.join("ok_r1.fastq.gz");
        let p2f2_path = fixture.dir.join("ok_r2.fastq.gz");
        let p2f1_size = fs::metadata(&p2f1_path)?.len();
        let p2f2_size = fs::metadata(&p2f2_path)?.len();
        total_bytes += p2f1_size + p2f2_size;
        jobs.push(Job::PairedFastq(PairedFastqJob {
            fq1_path: p2f1_path,
            fq2_path: p2f2_path,
            fq1_length_check: ReadLengthCheck::Fixed(4),
            fq2_length_check: ReadLengthCheck::Fixed(4),
            fq1_size: p2f1_size,
            fq2_size: p2f2_size,
        }));

        let s1_path = fixture.dir.join("badlen.fastq.gz");
        let s1_size = fs::metadata(&s1_path)?.len();
        total_bytes += s1_size;
        jobs.push(Job::SingleFastq(SingleFastqJob {
            path: s1_path,
            length_check: ReadLengthCheck::Auto,
            size: s1_size,
        }));

        run_check(jobs, total_bytes, &output, true, Some(false))?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 5);

        let find_report = |recs: &[TestReport], suffix: &str| -> TestReport {
            recs.iter()
                .find(|r| match r {
                    TestReport::Fastq(d) => d.path.ends_with(suffix),
                    _ => false,
                })
                .unwrap_or_else(|| panic!("Report for file ending in '{suffix}' not found"))
                .clone()
        };

        if let TestReport::Fastq(data) = find_report(&records, "counts1.fastq.gz") {
            assert_eq!(data.status, "ERROR");
            assert!(
                data.errors
                    .iter()
                    .any(|e| e.contains("Mismatched read counts"))
            );
        }

        if let TestReport::Fastq(data) = find_report(&records, "counts2.fastq.gz") {
            assert_eq!(data.status, "ERROR");
            assert!(
                data.errors
                    .iter()
                    .any(|e| e.contains("Mismatched read counts"))
            );
        }

        if let TestReport::Fastq(data) = find_report(&records, "badlen.fastq.gz") {
            assert_eq!(data.status, "ERROR");
            assert!(
                data.errors
                    .iter()
                    .any(|e| e.contains("Found inconsistent read length"))
            );
        }

        if let TestReport::Fastq(data) = find_report(&records, "ok_r1.fastq.gz") {
            assert_eq!(data.status, "OK");
            assert!(data.errors.is_empty());
        }
        if let TestReport::Fastq(data) = find_report(&records, "ok_r2.fastq.gz") {
            assert_eq!(data.status, "OK");
            assert!(data.errors.is_empty());
        }

        Ok(())
    }

    #[test]
    fn test_valid_bam_check() -> Result<()> {
        let dir = tempdir()?;
        let bam_path = dir.path().join("test.bam");

        let header = Header::builder()
            .add_read_group("rg0", Map::<ReadGroup>::default())
            .build();
        let mut writer = bam::io::Writer::new(fs::File::create(&bam_path)?);
        writer.write_header(&header)?;
        let record = record_buf::Builder::default()
            .set_name("r0")
            .set_flags(Flags::UNMAPPED)
            .set_sequence(b"ACGT".into())
            .set_quality_scores(QualityScores::from(vec![1, 1, 1, 1]))
            .build();

        writer.write_alignment_record(&header, &record)?;
        drop(writer);

        let output = dir.path().join("report.jsonl");
        let bam_size = fs::metadata(&bam_path)?.len();
        let jobs = vec![Job::Bam(BamCheckJob {
            path: bam_path,
            size: bam_size,
        })];

        run_check(jobs, bam_size, &output, true, Some(false))?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 1);

        if let TestReport::Bam(data) = &records[0] {
            assert_eq!(data.status, "OK");
            assert_eq!(data.num_records, Some(1));
            assert!(data.errors.is_empty());
            assert!(
                data.warnings.iter().any(|w| w.contains(
                    "Detected a header in BAM file, ensure it contains no private information!"
                )),
                "Expected to find BAM header warning"
            );
        } else {
            panic!("Expected a Bam report");
        }
        Ok(())
    }

    #[test]
    fn test_checksum_only() -> Result<()> {
        let dir = tempdir()?;
        let file_path = dir.path().join("raw.txt");
        let content = "some file contents";
        fs::write(&file_path, content)?;

        let expected_checksum = "cf57fcf9d6d7fb8fd7d8c30527c8f51026aa1d99ad77cc769dd0c757d4fe8667";

        let output = dir.path().join("report.jsonl");

        let file_size = fs::metadata(&file_path)?.len();
        let jobs = vec![Job::Raw(RawJob {
            path: file_path,
            size: file_size,
        })];

        run_check(jobs, file_size, &output, true, Some(false))?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 1);
        if let TestReport::Raw(data) = &records[0] {
            assert_eq!(data.status, "OK");
            assert_eq!(data.checksum.as_deref(), Some(expected_checksum));
            assert!(data.errors.is_empty());
        } else {
            panic!("Expected a Checksum report");
        }
        Ok(())
    }

    #[test]
    fn test_bam_with_multiple_secondary_alignments() -> Result<()> {
        let dir = tempdir()?;
        let bam_path = dir.path().join("secondary.bam");
        let header = Header::default();
        let mut writer = bam::io::Writer::new(fs::File::create(&bam_path)?);
        writer.write_header(&header)?;

        let rec1 = record_buf::Builder::default()
            .set_name("rec1")
            .set_flags(Flags::empty())
            .build();
        let rec2 = record_buf::Builder::default()
            .set_name("rec2_secondary")
            .set_flags(Flags::SECONDARY)
            .build();
        let rec3 = record_buf::Builder::default()
            .set_name("rec3_secondary")
            .set_flags(Flags::SECONDARY)
            .build();

        writer.write_alignment_record(&header, &rec1)?;
        writer.write_alignment_record(&header, &rec2)?;
        writer.write_alignment_record(&header, &rec3)?;
        drop(writer);

        let output = dir.path().join("report.jsonl");
        let bam_size = fs::metadata(&bam_path)?.len();
        let jobs = vec![Job::Bam(BamCheckJob {
            path: bam_path,
            size: bam_size,
        })];
        run_check(jobs, bam_size, &output, true, Some(false))?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 1);
        if let TestReport::Bam(data) = &records[0] {
            assert_eq!(data.status, "OK");
            assert_eq!(data.num_records, Some(3));
            assert_eq!(data.warnings.len(), 1);
            assert_eq!(
                data.warnings[0],
                "File contains 2 secondary alignment(s). First detected at record #2 ('rec2_secondary')."
            );
        } else {
            panic!("Expected a BAM report");
        }
        Ok(())
    }

    #[test]
    fn test_bam_with_multiple_hard_clipped_alignments() -> Result<()> {
        let dir = tempdir()?;
        let bam_path = dir.path().join("hardclip.bam");
        let header = Header::default();
        let mut writer = bam::io::Writer::new(fs::File::create(&bam_path)?);
        writer.write_header(&header)?;

        let cigar_hard_clip: record_buf::Cigar =
            [Op::new(Kind::HardClip, 5), Op::new(Kind::Match, 4)]
                .into_iter()
                .collect();
        let rec1 = record_buf::Builder::default()
            .set_name("rec1_hardclip")
            .set_flags(Flags::empty())
            .set_cigar(cigar_hard_clip.clone())
            .set_sequence(b"ACGT".into())
            .build();
        let rec2 = record_buf::Builder::default()
            .set_name("rec2_noclip")
            .set_flags(Flags::empty())
            .build();
        let rec3 = record_buf::Builder::default()
            .set_name("rec3_hardclip")
            .set_flags(Flags::empty())
            .set_cigar(cigar_hard_clip)
            .set_sequence(b"TGCA".into())
            .build();

        writer.write_alignment_record(&header, &rec1)?;
        writer.write_alignment_record(&header, &rec2)?;
        writer.write_alignment_record(&header, &rec3)?;
        drop(writer);

        let output = dir.path().join("report.jsonl");
        let bam_size = fs::metadata(&bam_path)?.len();
        let jobs = vec![Job::Bam(BamCheckJob {
            path: bam_path,
            size: bam_size,
        })];
        run_check(jobs, bam_size, &output, true, Some(false))?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 1);
        if let TestReport::Bam(data) = &records[0] {
            assert_eq!(data.status, "OK");
            assert_eq!(data.num_records, Some(3));
            assert_eq!(data.warnings.len(), 1);
            assert_eq!(
                data.warnings[0],
                "File contains 2 primary alignment(s) with hard-clipped bases. First detected at record #1 ('rec1_hardclip')."
            );
        } else {
            panic!("Expected a BAM report");
        }
        Ok(())
    }

    #[test]
    fn test_bam_with_mixed_warnings() -> Result<()> {
        let dir = tempdir()?;
        let bam_path = dir.path().join("mixed.bam");
        let header = Header::default();
        let mut writer = bam::io::Writer::new(fs::File::create(&bam_path)?);
        writer.write_header(&header)?;

        let cigar_hard_clip: record_buf::Cigar = [Op::new(Kind::HardClip, 5)].into_iter().collect();
        let rec1 = record_buf::Builder::default()
            .set_name("rec1_hardclip")
            .set_flags(Flags::empty())
            .set_cigar(cigar_hard_clip.clone())
            .build();
        let rec2 = record_buf::Builder::default()
            .set_name("rec2_secondary")
            .set_flags(Flags::SECONDARY)
            .build();
        let rec3 = record_buf::Builder::default()
            .set_name("rec3_hardclip")
            .set_flags(Flags::empty())
            .set_cigar(cigar_hard_clip)
            .build();
        let rec4 = record_buf::Builder::default()
            .set_name("rec4_secondary")
            .set_flags(Flags::SECONDARY)
            .build();

        writer.write_alignment_record(&header, &rec1)?;
        writer.write_alignment_record(&header, &rec2)?;
        writer.write_alignment_record(&header, &rec3)?;
        writer.write_alignment_record(&header, &rec4)?;
        drop(writer);

        let output = dir.path().join("report.jsonl");
        let bam_size = fs::metadata(&bam_path)?.len();
        let jobs = vec![Job::Bam(BamCheckJob {
            path: bam_path,
            size: bam_size,
        })];
        run_check(jobs, bam_size, &output, true, Some(false))?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 1);
        if let TestReport::Bam(data) = &records[0] {
            assert_eq!(data.status, "OK");
            assert_eq!(data.num_records, Some(4));
            assert_eq!(data.warnings.len(), 2);
            assert!(data.warnings.contains(&"File contains 2 secondary alignment(s). First detected at record #2 ('rec2_secondary').".to_string()));
            assert!(data.warnings.contains(&"File contains 2 primary alignment(s) with hard-clipped bases. First detected at record #1 ('rec1_hardclip').".to_string()));
        } else {
            panic!("Expected a BAM report");
        }
        Ok(())
    }
}
