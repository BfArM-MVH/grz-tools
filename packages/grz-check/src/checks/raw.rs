use crate::checker::FileReport;
use crate::checks::common::{CheckOutcome, check_file};
use indicatif::ProgressBar;
use std::io;
use std::path::{Path, PathBuf};

pub fn check_raw(path: &Path, file_pb: &ProgressBar, global_pb: &ProgressBar) -> FileReport {
    check_file(path, file_pb, global_pb, false, |reader| {
        match io::copy(reader, &mut io::sink()) {
            Ok(_) => Ok(CheckOutcome::default()),
            Err(e) => Err(format!("Failed to read file: {e}")),
        }
    })
}

#[derive(Debug)]
pub struct RawJob {
    pub path: PathBuf,
    pub size: u64,
}
