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
pub use checker::{FileReport, Stats, PairReport};
pub use checks::fastq::{ReadLengthCheck, SingleFastqJob, PairedFastqJob};
pub use checks::bam::BamCheckJob;
pub use checks::raw::RawJob;
