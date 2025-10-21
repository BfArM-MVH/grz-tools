use anyhow::{Context, Result};
use clap::{ArgGroup, Parser};
use std::fs;
use std::path::PathBuf;

use crate::checker::Job;
use crate::checks::bam::BamCheckJob;
use crate::checks::fastq::{PairedFastqJob, ReadLengthCheck, SingleFastqJob};
use crate::checks::raw::RawJob;

mod checker;
mod checks;
mod progress;
mod sha256;

/// Checks integrity of sequencing files (FASTQ, BAM).
///
/// Use --fastq-paired for paired-end FASTQ, --fastq-single for single-end FASTQ,
/// --bam for BAM files, or --raw for only calculating checksums of any file.
/// These flags can be used multiple times.
///
/// By default, the tool will exit immediately after the first error is found.
/// Use --continue-on-error to check all files regardless of errors.
#[derive(Debug, clap::Parser)]
#[command(author, version, about)]
#[command(group(
    ArgGroup::new("input_files")
        .required(true)
        .multiple(true)
))]
struct Args {
    /// Flag to show progress bars during processing.
    #[arg(long, global = true)]
    show_progress: Option<bool>,

    /// A paired-end FASTQ sample. Provide FQ1, FQ2, and minimum mean read length.
    /// Read Length: >0 for fixed, <0 to skip length check.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        allow_hyphen_values = true,
        num_args = 3,
        value_names = ["FQ1_PATH", "FQ2_PATH", "MIN_MEAN_READ_LEN"],
        group = "input_files"
    )]
    fastq_paired: Vec<String>,

    /// A single-end FASTQ sample. Provide the file path and minimum mean read length.
    /// Read Length: >0 for fixed, <0 to skip length check.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        allow_hyphen_values = true,
        num_args = 2,
        value_names = ["FQ_PATH", "MIN_MEAN_READ_LEN"],
        group = "input_files"
    )]
    fastq_single: Vec<String>,

    /// A single BAM file to validate.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        num_args = 1,
        value_names = ["BAM_PATH"],
        group = "input_files"
    )]
    bam: Vec<PathBuf>,

    /// A file for which to only calculate the SHA256 checksum, skipping all other validation.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        num_args = 1,
        value_names = ["FILE_PATH"],
        group = "input_files"
    )]
    raw: Vec<PathBuf>,

    /// Path to write the output JSONL report.
    #[arg(long, required = true)]
    output: PathBuf,

    /// Continue processing all files even if an error is found.
    #[arg(long, action = clap::ArgAction::SetTrue)]
    continue_on_error: bool,

    /// Number of threads to use for processing.
    #[arg(long)]
    threads: Option<usize>,
}

fn create_jobs(
    paired_raw: &[String],
    single_raw: &[String],
    bam_raw: &[PathBuf],
    raw: &[PathBuf],
) -> Result<(Vec<Job>, u64)> {
    let mut jobs = Vec::new();
    let mut total_bytes: u64 = 0;

    let parse_len = |len_str: &str| -> Result<ReadLengthCheck> {
        let len_val: i64 = len_str
            .parse()
            .context("Invalid read length. Must be an integer.")?;
        Ok(match len_val {
            v if v < 0 => ReadLengthCheck::Skip,
            v => ReadLengthCheck::Fixed(v as usize),
        })
    };

    for chunk in paired_raw.chunks_exact(3) {
        let fq1_path = PathBuf::from(&chunk[0]);
        let fq2_path = PathBuf::from(&chunk[1]);
        let length_check =
            parse_len(&chunk[2]).with_context(|| format!("Invalid read length '{}'", &chunk[2]))?;
        let fq1_size = fs::metadata(&fq1_path)?.len();
        let fq2_size = fs::metadata(&fq2_path)?.len();
        total_bytes += fq1_size + fq2_size;
        jobs.push(Job::PairedFastq(PairedFastqJob {
            fq1_path,
            fq2_path,
            length_check,
            fq1_size,
            fq2_size,
        }));
    }

    for chunk in single_raw.chunks_exact(2) {
        let path = PathBuf::from(&chunk[0]);
        let length_check = parse_len(&chunk[1]).with_context(|| {
            format!(
                "Invalid read length '{}' for file '{}'",
                &chunk[1], &chunk[0]
            )
        })?;
        let size = fs::metadata(&path)?.len();
        total_bytes += size;
        jobs.push(Job::SingleFastq(SingleFastqJob {
            path,
            length_check,
            size,
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

fn main() -> Result<()> {
    let args = Args::parse();

    let Args {
        fastq_paired,
        fastq_single,
        bam,
        raw,
        output,
        threads,
        continue_on_error,
        show_progress,
    } = args;

    if let Some(num_threads) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .context("Failed to set up Rayon thread pool")?;
    }

    let (jobs, total_bytes) = create_jobs(&fastq_paired, &fastq_single, &bam, &raw)?;

    checker::run_check(jobs, total_bytes, &output, continue_on_error, show_progress)?;

    Ok(())
}
