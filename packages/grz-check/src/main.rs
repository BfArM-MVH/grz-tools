use anyhow::{Context, Result};
use clap::Parser;
use std::path::PathBuf;

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
struct Args {
    /// Flag to show progress bars during processing.
    #[arg(long, global = true)]
    show_progress: Option<bool>,

    /// A paired-end FASTQ sample. Provide FQ1, FQ2, FQ1 read length, and FQ2 read length.
    /// Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        allow_hyphen_values = true,
        num_args = 4,
        value_names = ["FQ1_PATH", "FQ2_PATH", "FQ1_READ_LEN", "FQ2_READ_LEN"]
    )]
    fastq_paired: Vec<String>,

    /// A single-end FASTQ sample. Provide the file path and read length.
    /// Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        allow_hyphen_values = true,
        num_args = 2,
        value_names = ["FQ_PATH", "READ_LEN"]
    )]
    fastq_single: Vec<String>,

    /// A single BAM file to validate.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        num_args = 1,
        value_names = ["BAM_PATH"]
    )]
    bam: Vec<PathBuf>,

    /// A file for which to only calculate the SHA256 checksum, skipping all other validation.
    #[arg(
        long,
        action = clap::ArgAction::Append,
        num_args = 1,
        value_names = ["FILE_PATH"]
    )]
    raw: Vec<PathBuf>,

    /// Path to write the output TSV report.
    #[arg(long, required = true)]
    output: PathBuf,

    /// Continue processing all files even if an error is found.
    #[arg(long, action = clap::ArgAction::SetTrue)]
    continue_on_error: bool,

    /// Number of threads to use for processing.
    #[arg(long)]
    threads: Option<usize>,
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
    {
        if let Some(num_threads) = threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(num_threads)
                .build_global()
                .context("Failed to set up Rayon thread pool")?;
        }
        checker::run_check(
            fastq_paired,
            fastq_single,
            bam,
            raw,
            &output,
            continue_on_error,
            show_progress,
        )?;
    }

    Ok(())
}
