use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{MultiProgress, ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use serde::Serialize;
use std::error::Error as StdError;
use std::fmt;
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};

#[derive(Debug, clap::Parser)]
#[command(author, version, about)]
struct Args {
    #[command(subcommand)]
    cmd: Commands,

    /// Flag to show progress bars during processing.
    #[arg(long)]
    show_progress: Option<bool>,
}

#[derive(Debug, Parser)]
enum Commands {
    /// Checks integrity and consistency of FASTQ files.
    ///
    /// Use --paired or --single for each sample. These flags can be used multiple times.
    ///
    /// By default, the tool will exit immediately after the first error is found.
    /// Use --continue-on-error to check all files regardless of errors.
    Fastq {
        /// A paired-end sample. Provide FQ1, FQ2, FQ1 read length, and FQ2 read length.
        ///
        /// Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check.
        /// Example: --paired fq1.fastq.gz fq2.fastq.gz 150 151
        #[arg(
            long,
            action = clap::ArgAction::Append,
            num_args = 4,
            value_names = ["FQ1_PATH", "FQ2_PATH", "FQ1_READ_LEN", "FQ2_READ_LEN"]
        )]
        paired: Vec<String>,

        /// A single-end sample. Provide the file path and read length.
        ///
        /// Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check.
        /// Example: --single sample.fastq.gz 0
        #[arg(
            long,
            action = clap::ArgAction::Append,
            num_args = 2,
            value_names = ["FQ_PATH", "READ_LEN"]
        )]
        single: Vec<String>,

        /// Path to write the output TSV report.
        #[arg(long, required = true)]
        output: PathBuf,

        /// Continue processing all files even if an error is found.
        #[arg(long, action = clap::ArgAction::SetTrue)]
        continue_on_error: bool,

        /// Number of threads to use for processing.
        #[arg(long)]
        threads: Option<usize>,
    },
}

#[derive(Debug, Copy, Clone, Eq, PartialEq, Serialize)]
struct Stats {
    num_reads: u64,
    read_length: usize,
}

#[derive(Debug, Serialize, Clone)]
struct FileReport {
    path: PathBuf,
    stats: Option<Stats>,
    errors: Vec<String>,
}

impl FileReport {
    fn is_ok(&self) -> bool {
        self.errors.is_empty()
    }
}

#[derive(Debug, Serialize, Clone)]
struct PairReport {
    fq1_report: FileReport,
    fq2_report: Option<FileReport>,
    pair_errors: Vec<String>,
}

impl PairReport {
    fn is_error(&self) -> bool {
        !self.fq1_report.is_ok()
            || !self.pair_errors.is_empty()
            || self.fq2_report.as_ref().is_some_and(|r2| !r2.is_ok())
    }
}

#[derive(Debug, Copy, Clone)]
enum ReadLengthCheck {
    Fixed(usize),
    Auto,
    Skip,
}

#[derive(Debug)]
struct EarlyExitError(PairReport);

impl fmt::Display for EarlyExitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "A validation error occurred, exiting.")
    }
}

impl StdError for EarlyExitError {}

fn check_single_fastq(path: &Path, length_check: ReadLengthCheck, pb: &ProgressBar) -> FileReport {
    pb.set_message(
        path.file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string(),
    );

    let file = match fs::File::open(path) {
        Ok(file) => file,
        Err(e) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                errors: vec![format!("Failed to open file for reading: {}", e)],
            };
        }
    };

    let progress_reader = pb.wrap_read(file);
    let decompressed_reader = match niffler::get_reader(Box::new(BufReader::new(progress_reader))) {
        Ok((reader, _)) => reader,
        Err(e) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                errors: vec![format!("Failed to decompress file: {}", e)],
            };
        }
    };

    let mut reader = noodles::fastq::io::Reader::new(BufReader::new(decompressed_reader));
    let mut records = reader.records();

    let first_record = match records.next() {
        Some(Ok(rec)) => rec,
        Some(Err(e)) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                errors: vec![format!("Failed to parse first record: {}", e)],
            };
        }
        None => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                errors: vec!["File is empty. Expected at least one record.".to_string()],
            };
        }
    };

    let file_read_length = first_record.sequence().len();
    let mut num_reads: u64 = 1;

    if let ReadLengthCheck::Fixed(expected) = length_check {
        if expected != file_read_length {
            return FileReport {
                path: path.to_path_buf(),
                stats: Some(Stats {
                    num_reads,
                    read_length: file_read_length,
                }),
                errors: vec![format!(
                    "Provided read length ({}) does not match first read's length ({}).",
                    expected, file_read_length
                )],
            };
        }
    }

    if !matches!(length_check, ReadLengthCheck::Skip) {
        for (i, record) in records.enumerate() {
            num_reads += 1;
            let record = match record {
                Ok(rec) => rec,
                Err(e) => {
                    return FileReport {
                        path: path.to_path_buf(),
                        stats: Some(Stats {
                            num_reads,
                            read_length: file_read_length,
                        }),
                        errors: vec![format!("Failed to parse record #{}: {}", i + 2, e)],
                    };
                }
            };

            if record.sequence().len() != file_read_length {
                return FileReport {
                    path: path.to_path_buf(),
                    stats: Some(Stats {
                        num_reads,
                        read_length: file_read_length,
                    }),
                    errors: vec![format!(
                        "Found inconsistent read length at record #{}. Expected {}, but got {}.",
                        i + 2,
                        file_read_length,
                        record.sequence().len()
                    )],
                };
            }
        }
    } else {
        num_reads += records.count() as u64;
    }

    FileReport {
        path: path.to_path_buf(),
        stats: Some(Stats {
            num_reads,
            read_length: file_read_length,
        }),
        errors: vec![],
    }
}

#[derive(Debug)]
struct FastqCheckJob {
    fq1: PathBuf,
    fq2: Option<PathBuf>,
    fq1_length_check: ReadLengthCheck,
    fq2_length_check: Option<ReadLengthCheck>,
    fq1_size: u64,
    fq2_size: Option<u64>,
}

fn create_jobs_from_raw_input(
    paired_raw: &[String],
    single_raw: &[String],
) -> Result<(Vec<FastqCheckJob>, u64)> {
    let mut jobs = Vec::new();
    let mut total_bytes: u64 = 0;

    let parse_len = |len_str: &str| -> Result<ReadLengthCheck> {
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
        let fq1_str = &chunk[0];
        let fq2_str = &chunk[1];
        let len1_str = &chunk[2];
        let len2_str = &chunk[3];

        let r1 = PathBuf::from(fq1_str);
        let r2 = PathBuf::from(fq2_str);

        let fq1_length_check = parse_len(len1_str).with_context(|| {
            format!("Invalid read length '{}' for file '{}'", len1_str, fq1_str)
        })?;
        let fq2_length_check = parse_len(len2_str).with_context(|| {
            format!("Invalid read length '{}' for file '{}'", len2_str, fq2_str)
        })?;

        let fq1_size = fs::metadata(&r1)
            .with_context(|| format!("Could not get metadata for {}", r1.display()))?
            .len();
        let fq2_size = fs::metadata(&r2)
            .with_context(|| format!("Could not get metadata for {}", r2.display()))?
            .len();
        total_bytes += fq1_size + fq2_size;

        jobs.push(FastqCheckJob {
            fq1: r1,
            fq2: Some(r2),
            fq1_length_check,
            fq2_length_check: Some(fq2_length_check),
            fq1_size,
            fq2_size: Some(fq2_size),
        });
    }

    for chunk in single_raw.chunks_exact(2) {
        let path_str = &chunk[0];
        let len_str = &chunk[1];

        let path = PathBuf::from(path_str);
        let fq1_length_check = parse_len(len_str).with_context(|| {
            format!("Invalid read length '{}' for file '{}'", len_str, path_str)
        })?;

        let fq1_size = fs::metadata(&path)
            .with_context(|| format!("Could not get metadata for {}", path.display()))?
            .len();
        total_bytes += fq1_size;

        jobs.push(FastqCheckJob {
            fq1: path,
            fq2: None,
            fq1_length_check,
            fq2_length_check: None,
            fq1_size,
            fq2_size: None,
        });
    }

    Ok((jobs, total_bytes))
}

fn run_check(
    paired_raw: Vec<String>,
    single_raw: Vec<String>,
    output: &Path,
    continue_on_error: bool,
    show_progress: Option<bool>,
) -> Result<()> {
    if paired_raw.is_empty() && single_raw.is_empty() {
        anyhow::bail!(
            "No input files provided. Use --paired or --single to specify files to check."
        );
    }

    let (jobs, total_bytes) = create_jobs_from_raw_input(&paired_raw, &single_raw)?;

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
                       job: FastqCheckJob|
     -> PairReport {
        let fq1_pb = m.add(ProgressBar::new(job.fq1_size));
        fq1_pb.set_style(style.clone());
        fq1_pb.set_prefix("Checking R1");
        let fq1_report = check_single_fastq(&job.fq1, job.fq1_length_check, &fq1_pb);
        if fq1_report.is_ok() {
            fq1_pb.finish_with_message("✓ OK");
        } else {
            fq1_pb.abandon_with_message("✗ ERROR");
        }
        main_pb.inc(job.fq1_size);

        let mut fq2_report = None;
        if let (Some(fq2_path), Some(fq2_size), Some(fq2_len_check)) =
            (&job.fq2, job.fq2_size, job.fq2_length_check)
        {
            let fq2_pb = m.add(ProgressBar::new(fq2_size));
            fq2_pb.set_style(style.clone());
            fq2_pb.set_prefix("Checking R2");
            let report = check_single_fastq(fq2_path, fq2_len_check, &fq2_pb);
            if report.is_ok() {
                fq2_pb.finish_with_message("✓ OK");
            } else {
                fq2_pb.abandon_with_message("✗ ERROR");
            }
            main_pb.inc(fq2_size);
            fq2_report = Some(report);
        }

        let mut pair_errors = Vec::new();
        if let (Some(fq1_stats), Some(fq2_report_val)) = (fq1_report.stats, fq2_report.as_ref()) {
            if let Some(fq2_stats) = fq2_report_val.stats {
                if fq1_stats.num_reads != fq2_stats.num_reads {
                    pair_errors.push(format!(
                        "Mismatched read counts: {} vs {}",
                        fq1_stats.num_reads, fq2_stats.num_reads
                    ));
                }
                if matches!(job.fq1_length_check, ReadLengthCheck::Auto)
                    && job
                        .fq2_length_check
                        .is_some_and(|c| matches!(c, ReadLengthCheck::Auto))
                    && fq1_stats.read_length != fq2_stats.read_length
                {
                    pair_errors.push(format!(
                        "Mismatched read lengths in auto-detect mode: {} vs {}",
                        fq1_stats.read_length, fq2_stats.read_length
                    ));
                }
            }
        }
        PairReport {
            fq1_report,
            fq2_report,
            pair_errors,
        }
    };

    if continue_on_error {
        let reports: Vec<PairReport> = jobs
            .into_par_iter()
            .map_with((mpb.clone(), main_pb.clone(), file_style), process_job)
            .collect();

        let num_failed_jobs = reports.iter().filter(|r| r.is_error()).count();
        write_report(&reports, output)?;

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
        // fail-fast
        let result: Result<Vec<PairReport>, EarlyExitError> = jobs
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
                write_report(&reports, output)?;
                main_pb.finish_with_message(format!(
                    "All checks passed! Report written to {}",
                    output.display()
                ));
            }
            Err(EarlyExitError(failed_report)) => {
                let error_msg = format!(
                    "An error occurred in {}. See report for details. Aborting due to fail-fast mode.",
                    failed_report.fq1_report.path.display()
                );
                write_report(&[failed_report], output)?;
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

fn write_report(reports: &[PairReport], output_path: &Path) -> Result<()> {
    let mut writer = csv::WriterBuilder::default()
        .delimiter(b'\t')
        .from_path(output_path)
        .with_context(|| {
            format!(
                "Failed to create report file at '{}'",
                output_path.display()
            )
        })?;

    writer.write_record(["path", "status", "num_reads", "read_length", "errors"])?;

    for report in reports {
        let fq1_status = if report.fq1_report.is_ok() && report.pair_errors.is_empty() {
            "OK"
        } else {
            "ERROR"
        };
        let mut all_errors = report.fq1_report.errors.clone();
        all_errors.extend(report.pair_errors.clone());

        let fq1_stats = report.fq1_report.stats.as_ref();
        writer.write_record(&[
            report.fq1_report.path.to_string_lossy().to_string(),
            fq1_status.to_string(),
            fq1_stats
                .map(|s| s.num_reads.to_string())
                .unwrap_or_else(|| "N/A".to_string()),
            fq1_stats
                .map(|s| s.read_length.to_string())
                .unwrap_or_else(|| "N/A".to_string()),
            all_errors.join("; "),
        ])?;

        if let Some(fq2_report) = &report.fq2_report {
            let fq2_status = if fq2_report.is_ok() && report.pair_errors.is_empty() {
                "OK"
            } else {
                "ERROR"
            };
            let mut all_errors = fq2_report.errors.clone();
            all_errors.extend(report.pair_errors.clone());

            let fq2_stats = fq2_report.stats.as_ref();
            writer.write_record(&[
                fq2_report.path.to_string_lossy().to_string(),
                fq2_status.to_string(),
                fq2_stats
                    .map(|s| s.num_reads.to_string())
                    .unwrap_or_else(|| "N/A".to_string()),
                fq2_stats
                    .map(|s| s.read_length.to_string())
                    .unwrap_or_else(|| "N/A".to_string()),
                all_errors.join("; "),
            ])?;
        }
    }
    writer.flush()?;
    Ok(())
}

fn main() -> Result<()> {
    let args = Args::parse();

    match args.cmd {
        Commands::Fastq {
            paired,
            single,
            output,
            threads,
            continue_on_error,
        } => {
            if let Some(num_threads) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build_global()
                    .context("Failed to set up Rayon thread pool")?;
            }
            run_check(
                paired,
                single,
                &output,
                continue_on_error,
                args.show_progress,
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;
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

    fn read_report_status(report_path: &Path) -> Result<Vec<(String, String, String)>> {
        let mut reader = csv::ReaderBuilder::default()
            .delimiter(b'\t')
            .from_path(&report_path)?;
        reader
            .records()
            .map(|r| {
                let record = r?;
                Ok((
                    record[0].to_string(),
                    record[1].to_string(),
                    record[4].to_string(),
                ))
            })
            .collect()
    }

    #[test]
    fn test_valid_pair_with_different_lengths() -> Result<()> {
        let fixture = TestFiles::new()?;
        let output = fixture.dir.join("report.tsv");
        let paired = vec![
            fixture.path_str("ok_r1.fastq.gz"),
            fixture.path_str("ok_r2_len5.fastq.gz"),
            "4".to_string(),
            "5".to_string(),
        ];

        run_check(paired, vec![], &output, false, Some(false))?;

        let statuses = read_report_status(&output)?;
        assert_eq!(statuses.len(), 2);
        assert_eq!(statuses[0].1, "OK");
        assert_eq!(statuses[1].1, "OK");
        Ok(())
    }

    #[test]
    fn test_multiple_inputs_with_continue_on_error() -> Result<()> {
        let fixture = TestFiles::new()?;
        let output = fixture.dir.join("report.tsv");

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

        run_check(all_paired, all_single, &output, true, Some(false))?;

        let statuses = read_report_status(&output)?;
        assert_eq!(statuses.len(), 5);

        assert_eq!(statuses[0].0, fixture.path_str("counts1.fastq.gz"));
        assert_eq!(statuses[0].1, "ERROR");
        assert!(statuses[0].2.contains("Mismatched read counts"));

        assert_eq!(statuses[2].0, fixture.path_str("ok_r1.fastq.gz"));
        assert_eq!(statuses[2].1, "OK");

        assert_eq!(statuses[4].0, fixture.path_str("badlen.fastq.gz"));
        assert_eq!(statuses[4].1, "ERROR");
        Ok(())
    }
}
