//! grz-check library for validating sequencing files (FASTQ, BAM)
//!
//! This library provides validation functionality for sequencing files,
//! exposed both as a Rust library and Python bindings via PyO3.

pub mod checker;
pub mod checks;
pub mod progress;
pub mod sha256;

#[cfg(feature = "python")]
pub mod python;

// Re-export core types for library users
pub use checker::{FileReport, PairReport, Stats};
pub use checks::bam::{BamCheckJob, validate_bam_data};
pub use checks::fastq::{PairedFastqJob, ReadLengthCheck, SingleFastqJob, validate_fastq_data};
pub use checks::raw::RawJob;
