use crate::checks::bam::BamCheckJob;
use crate::checks::fastq::{FastqCheckJob, ReadLengthCheck};
use crate::checks::raw::RawJob;
use crate::checks::{bam, fastq, raw};
use anyhow::Context;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use serde::Serialize;
use std::error::Error as StdError;
use std::fmt;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

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
    pub fq2_report: Option<FileReport>,
    pub pair_errors: Vec<String>,
}

impl PairReport {
    fn is_error(&self) -> bool {
        !self.fq1_report.is_ok()
            || !self.pair_errors.is_empty()
            || self.fq2_report.as_ref().is_some_and(|r2| !r2.is_ok())
    }
}

#[derive(Debug)]
enum Job {
    Fastq(FastqCheckJob),
    Bam(BamCheckJob),
    Raw(RawJob),
}

#[derive(Debug)]
enum CheckResult {
    PairedFastq(PairReport),
    SingleFastq(FileReport),
    Bam(FileReport),
    Checksum(FileReport),
}

impl CheckResult {
    fn is_error(&self) -> bool {
        match self {
            CheckResult::PairedFastq(r) => r.is_error(),
            CheckResult::SingleFastq(r) => !r.is_ok(),
            CheckResult::Bam(r) => !r.is_ok(),
            CheckResult::Checksum(r) => !r.is_ok(),
        }
    }
    fn primary_path(&self) -> &Path {
        match self {
            CheckResult::PairedFastq(r) => &r.fq1_report.path,
            CheckResult::SingleFastq(r) => &r.path,
            CheckResult::Bam(r) => &r.path,
            CheckResult::Checksum(r) => &r.path,
        }
    }
}

#[derive(Debug)]
pub struct EarlyExitError(CheckResult);

impl fmt::Display for EarlyExitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "A validation error occurred, exiting.")
    }
}

impl StdError for EarlyExitError {}

fn create_jobs(
    paired_raw: &[String],
    single_raw: &[String],
    bam_raw: &[PathBuf],
    raw: &[PathBuf],
) -> anyhow::Result<(Vec<Job>, u64)> {
    let mut jobs = Vec::new();
    let mut total_bytes: u64 = 0;

    let parse_len = |len_str: &str| -> anyhow::Result<ReadLengthCheck> {
        let len_val: i64 = len_str
            .parse()
            .context("Invalid read length. Must be an integer.")?;
        Ok(match len_val {
            v if v < 0 => ReadLengthCheck::Skip,
            0 => ReadLengthCheck::Auto,
            v => ReadLengthCheck::Fixed(v as usize),
        })
    };

    for chunk in paired_raw.chunks_exact(4) {
        let fq1_path = PathBuf::from(&chunk[0]);
        let fq2_path = PathBuf::from(&chunk[1]);
        let fq1_length_check = parse_len(&chunk[2]).with_context(|| {
            format!(
                "Invalid read length '{}' for file '{}'",
                &chunk[2], &chunk[0]
            )
        })?;
        let fq2_length_check = parse_len(&chunk[3]).with_context(|| {
            format!(
                "Invalid read length '{}' for file '{}'",
                &chunk[3], &chunk[1]
            )
        })?;
        let fq1_size = fs::metadata(&fq1_path)?.len();
        let fq2_size = fs::metadata(&fq2_path)?.len();
        total_bytes += fq1_size + fq2_size;
        jobs.push(Job::Fastq(FastqCheckJob {
            fq1: fq1_path,
            fq2: Some(fq2_path),
            fq1_length_check,
            fq2_length_check: Some(fq2_length_check),
            fq1_size,
            fq2_size: Some(fq2_size),
        }));
    }

    for chunk in single_raw.chunks_exact(2) {
        let path = PathBuf::from(&chunk[0]);
        let fq1_length_check = parse_len(&chunk[1]).with_context(|| {
            format!(
                "Invalid read length '{}' for file '{}'",
                &chunk[1], &chunk[0]
            )
        })?;
        let fq1_size = fs::metadata(&path)?.len();
        total_bytes += fq1_size;
        jobs.push(Job::Fastq(FastqCheckJob {
            fq1: path,
            fq2: None,
            fq1_length_check,
            fq2_length_check: None,
            fq1_size,
            fq2_size: None,
        }));
    }

    for path_str in bam_raw {
        let path = PathBuf::from(path_str);
        let size = fs::metadata(&path)?.len();
        total_bytes += size;
        jobs.push(Job::Bam(BamCheckJob { path, size }));
    }

    for path_str in raw {
        let path = PathBuf::from(path_str);
        let size = fs::metadata(&path)
            .with_context(|| format!("Could not get metadata for {}", path.display()))?
            .len();
        total_bytes += size;
        jobs.push(Job::Raw(RawJob { path, size }));
    }

    Ok((jobs, total_bytes))
}

pub fn run_check(
    paired_raw: Vec<String>,
    single_raw: Vec<String>,
    bam_raw: Vec<PathBuf>,
    raw: Vec<PathBuf>,
    output: &Path,
    continue_on_error: bool,
    show_progress: Option<bool>,
) -> anyhow::Result<()> {
    if paired_raw.is_empty() && single_raw.is_empty() && bam_raw.is_empty() && raw.is_empty() {
        anyhow::bail!(
            "No input files provided. Use --fastq-single, --fastq-paired, --bam, or --raw."
        );
    }

    let (jobs, total_bytes) = create_jobs(&paired_raw, &single_raw, &bam_raw, &raw)?;

    let mpb = MultiProgress::new();
    match show_progress {
        Some(true) => {
            mpb.set_draw_target(ProgressDrawTarget::stderr());
        }
        Some(false) => {
            mpb.set_draw_target(ProgressDrawTarget::hidden());
        }
        _ => { // keep default
        }
    }

    let file_style = ProgressStyle::with_template(
        "{prefix:15.bold} [{bar:40.cyan/blue}] {bytes}/{total_bytes} ({bytes_per_sec}, {eta}) {wide_msg}",
    )?.progress_chars("#>-");
    let main_pb = mpb.add(ProgressBar::new(total_bytes));
    main_pb.set_style(file_style.clone());
    main_pb.set_prefix("Overall");

    let process_job = |(m, main_pb, style): &mut (MultiProgress, ProgressBar, ProgressStyle),
                       job: Job|
     -> CheckResult {
        match job {
            Job::Fastq(job) => {
                let fq1_pb = m.add(ProgressBar::new(job.fq1_size));
                fq1_pb.set_style(style.clone());
                fq1_pb.set_prefix("Checking R1");
                let fq1_report =
                    fastq::check_single_fastq(&job.fq1, job.fq1_length_check, &fq1_pb, main_pb);
                if fq1_report.is_ok() {
                    fq1_pb.finish_with_message("✓ OK");
                } else {
                    fq1_pb.abandon_with_message("✗ ERROR");
                }

                if let (Some(fq2_path), Some(fq2_size), Some(fq2_len_check)) =
                    (&job.fq2, job.fq2_size, job.fq2_length_check)
                {
                    let fq2_pb = m.add(ProgressBar::new(fq2_size));
                    fq2_pb.set_style(style.clone());
                    fq2_pb.set_prefix("Checking R2");
                    let fq2_report =
                        fastq::check_single_fastq(fq2_path, fq2_len_check, &fq2_pb, main_pb);
                    if fq2_report.is_ok() {
                        fq2_pb.finish_with_message("✓ OK");
                    } else {
                        fq2_pb.abandon_with_message("✗ ERROR");
                    }

                    let mut pair_errors = Vec::new();
                    if let (Some(fq1_stats), Some(fq2_stats)) = (fq1_report.stats, fq2_report.stats)
                    {
                        if fq1_stats.num_records != fq2_stats.num_records {
                            pair_errors.push(format!(
                                "Mismatched read counts: {} vs {}",
                                fq1_stats.num_records, fq2_stats.num_records
                            ));
                        }
                        if matches!(job.fq1_length_check, ReadLengthCheck::Auto)
                            && job
                                .fq2_length_check
                                .is_some_and(|c| matches!(c, ReadLengthCheck::Auto))
                            && fq1_stats.read_length != fq2_stats.read_length
                        {
                            pair_errors.push(format!(
                                "Mismatched read lengths in auto-detect mode: {:?} vs {:?}",
                                fq1_stats.read_length, fq2_stats.read_length
                            ));
                        }
                    }
                    CheckResult::PairedFastq(PairReport {
                        fq1_report,
                        fq2_report: Some(fq2_report),
                        pair_errors,
                    })
                } else {
                    CheckResult::SingleFastq(fq1_report)
                }
            }
            Job::Bam(job) => {
                let pb = m.add(ProgressBar::new(job.size));
                pb.set_style(style.clone());
                pb.set_prefix("Checking BAM");
                let report = bam::check_bam(&job.path, &pb, main_pb);
                if report.is_ok() {
                    pb.finish_with_message("✓ OK");
                } else {
                    pb.abandon_with_message("✗ ERROR");
                }
                CheckResult::Bam(report)
            }
            Job::Raw(job) => {
                let pb = m.add(ProgressBar::new(job.size));
                pb.set_style(style.clone());
                pb.set_prefix("Checksum");
                let report = raw::check_raw(&job.path, &pb, main_pb);
                if report.is_ok() {
                    pb.finish_with_message("✓ OK");
                } else {
                    pb.abandon_with_message("✗ ERROR");
                }
                CheckResult::Checksum(report)
            }
        }
    };

    let mut output_writer = fs::File::create(output)
        .with_context(|| format!("Failed to create report file at {}", output.display()))?;

    if continue_on_error {
        let reports: Vec<CheckResult> = jobs
            .into_par_iter()
            .map_with((mpb.clone(), main_pb.clone(), file_style), process_job)
            .collect();

        let num_failed_jobs = reports.iter().filter(|r| r.is_error()).count();
        write_jsonl_report(&reports, &mut output_writer)?;

        if num_failed_jobs > 0 {
            main_pb.abandon_with_message(format!(
                "Processing complete. {} pairs/files failed. See report: {}",
                num_failed_jobs,
                output.display()
            ));
        } else {
            main_pb.finish_with_message(format!(
                "All checks passed! Report written to {}",
                output.display()
            ));
        }
    } else {
        let result: Result<Vec<CheckResult>, EarlyExitError> = jobs
            .into_par_iter()
            .map_with((mpb.clone(), main_pb.clone(), file_style), |state, job| {
                let report = process_job(state, job);
                if report.is_error() {
                    Err(EarlyExitError(report))
                } else {
                    Ok(report)
                }
            })
            .collect();

        match result {
            Ok(reports) => {
                write_jsonl_report(&reports, &mut output_writer)?;
                main_pb.finish_with_message(format!(
                    "All checks passed! Report written to {}",
                    output.display()
                ));
            }
            Err(EarlyExitError(failed_report)) => {
                let error_msg = format!(
                    "An error occurred in {}. See report for details. Aborting due to fail-fast mode.",
                    failed_report.primary_path().display()
                );
                write_jsonl_report(&[failed_report], &mut output_writer)?;
                main_pb
                    .abandon_with_message(format!("Error found. See report: {}", output.display()));
                mpb.clear()?;
                anyhow::bail!(error_msg);
            }
        }
    }

    mpb.clear()?;
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

fn write_jsonl_report<W: Write>(results: &[CheckResult], writer: &mut W) -> anyhow::Result<()> {
    for result in results {
        match result {
            CheckResult::PairedFastq(pair_report) => {
                let is_pair_error = !pair_report.pair_errors.is_empty();

                // Process R1
                let r1 = &pair_report.fq1_report;
                let mut fq1_errors = r1.errors.clone();
                if is_pair_error {
                    fq1_errors.extend(pair_report.pair_errors.clone());
                }
                let status = if r1.is_ok() && !is_pair_error {
                    "OK"
                } else {
                    "ERROR"
                };

                let report1 = JsonReport::Fastq(FastqReport {
                    path: &r1.path,
                    status,
                    num_records: r1.stats.map(|s| s.num_records),
                    read_length: r1.stats.and_then(|s| s.read_length),
                    checksum: r1.sha256.as_ref(),
                    errors: fq1_errors,
                    warnings: &r1.warnings,
                });
                serde_json::to_writer(&mut *writer, &report1)?;
                writer.write_all(b"\n")?;

                // Process R2, if it exists
                if let Some(r2) = &pair_report.fq2_report {
                    let mut fq2_errors = r2.errors.clone();
                    if is_pair_error {
                        fq2_errors.extend(pair_report.pair_errors.clone());
                    }
                    let status = if r2.is_ok() && !is_pair_error {
                        "OK"
                    } else {
                        "ERROR"
                    };

                    let report2 = JsonReport::Fastq(FastqReport {
                        path: &r2.path,
                        status,
                        num_records: r2.stats.map(|s| s.num_records),
                        read_length: r2.stats.and_then(|s| s.read_length),
                        checksum: r2.sha256.as_ref(),
                        errors: fq2_errors,
                        warnings: &r2.warnings,
                    });
                    serde_json::to_writer(&mut *writer, &report2)?;
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
            CheckResult::Checksum(report) => {
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
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::{Result, anyhow};
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use noodles::bam;

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

        fn path_str(&self, filename: &str) -> String {
            self.dir.join(filename).to_string_lossy().to_string()
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
        let paired = vec![
            fixture.path_str("ok_r1.fastq.gz"),
            fixture.path_str("ok_r2_len5.fastq.gz"),
            "4".to_string(),
            "5".to_string(),
        ];

        run_check(paired, vec![], vec![], vec![], &output, false, Some(false))?;

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

        let all_paired = vec![
            fixture.path_str("counts1.fastq.gz"),
            fixture.path_str("counts2.fastq.gz"),
            "0".to_string(),
            "0".to_string(),
            fixture.path_str("ok_r1.fastq.gz"),
            fixture.path_str("ok_r2.fastq.gz"),
            "4".to_string(),
            "4".to_string(),
        ];

        let all_single = vec![fixture.path_str("badlen.fastq.gz"), "0".to_string()];

        run_check(
            all_paired,
            all_single,
            vec![],
            vec![],
            &output,
            true,
            Some(false),
        )?;

        let records = read_jsonl_report(&output)?;
        assert_eq!(records.len(), 5);

        let find_report = |recs: &[TestReport], suffix: &str| -> TestReport {
            recs.iter()
                .find(|r| match r {
                    TestReport::Fastq(d) => d.path.ends_with(suffix),
                    _ => false,
                })
                .unwrap_or_else(|| panic!("Report for file ending in '{}' not found", suffix))
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

        run_check(
            vec![],
            vec![],
            vec![bam_path],
            vec![],
            &output,
            true,
            Some(false),
        )?;

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

        run_check(
            vec![],
            vec![],
            vec![],
            vec![file_path],
            &output,
            true,
            Some(false),
        )?;

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
        run_check(
            vec![],
            vec![],
            vec![bam_path],
            vec![],
            &output,
            true,
            Some(false),
        )?;

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
        run_check(
            vec![],
            vec![],
            vec![bam_path],
            vec![],
            &output,
            true,
            Some(false),
        )?;

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
        run_check(
            vec![],
            vec![],
            vec![bam_path],
            vec![],
            &output,
            true,
            Some(false),
        )?;

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
