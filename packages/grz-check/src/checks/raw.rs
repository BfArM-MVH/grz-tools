use crate::checker::{DataSource, FileReport};
use crate::checks::common::{CheckOutcome, check_data};
use indicatif::ProgressBar;
use std::io;

pub fn check_raw(input: DataSource, file_pb: &ProgressBar, global_pb: &ProgressBar) -> FileReport {
    check_data(input, file_pb, global_pb, false, |reader| {
        match io::copy(reader, &mut io::sink()) {
            Ok(_) => Ok(CheckOutcome::default()),
            Err(e) => Err(format!("Failed to read file: {e}")),
        }
    })
}

#[derive(Debug)]
pub struct RawJob {
    pub input: DataSource,
}
