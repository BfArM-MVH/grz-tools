use crate::checker::{FileReport, Stats};
use crate::progress::DualProgressReader;
use crate::sha256::SharedHashingReader;
use anyhow::Context;
use indicatif::ProgressBar;
use sha2::{Digest, Sha256};
use std::fs;
use std::io::{BufReader, Read};
use std::path::Path;
use std::sync::{Arc, Mutex};
#[derive(Debug, Default)]
pub struct CheckOutcome {
    pub stats: Option<Stats>,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

type ReaderAndHasher = (Box<dyn Read>, Arc<Mutex<Sha256>>);

pub fn setup_file_reader(
    path: &Path,
    file_pb: &ProgressBar,
    global_pb: &ProgressBar,
    decompress: bool,
) -> anyhow::Result<ReaderAndHasher> {
    file_pb.set_message(format!(
        "~ CHECK {}",
        path.file_name().unwrap_or_default().to_string_lossy()
    ));

    let file = fs::File::open(path)
        .with_context(|| format!("Failed to open file for reading: {}", path.display()))?;

    let hasher = Arc::new(Mutex::new(Sha256::new()));
    let hashing_reader = SharedHashingReader::new(BufReader::new(file), hasher.clone());
    let progress_reader =
        DualProgressReader::new(hashing_reader, file_pb.clone(), global_pb.clone());

    let reader: Box<dyn Read> = if decompress {
        let (decompressed_reader, _) = niffler::get_reader(Box::new(progress_reader))
            .with_context(|| format!("Failed to decompress file: {}", path.display()))?;
        decompressed_reader
    } else {
        Box::new(progress_reader)
    };

    Ok((reader, hasher))
}

pub fn check_file<F>(
    path: &Path,
    file_pb: &ProgressBar,
    global_pb: &ProgressBar,
    decompress: bool,
    logic: F,
) -> FileReport
where
    F: FnOnce(&mut dyn Read) -> Result<CheckOutcome, String>,
{
    let (mut reader, hasher) = match setup_file_reader(path, file_pb, global_pb, decompress) {
        Ok(setup) => setup,
        Err(e) => return FileReport::new_with_error(path, e.to_string()),
    };

    let outcome = match logic(&mut reader) {
        Ok(outcome) => outcome,
        Err(error_msg) => {
            return FileReport::new_with_error(path, error_msg);
        }
    };

    // Ensure the reader is fully consumed, such that the hasher can finalize
    drop(reader);

    let checksum = match Arc::try_unwrap(hasher) {
        Ok(mutex) => {
            let final_hasher = mutex.into_inner().unwrap();
            Some(format!("{:x}", final_hasher.finalize()))
        }
        Err(_) => {
            let mut final_report = FileReport::new(path, outcome.stats, vec![], outcome.warnings);
            final_report
                .errors
                .push("Failed to finalize checksum: hasher is still in use.".to_string());
            return final_report;
        }
    };

    FileReport::new(path, outcome.stats, outcome.errors, outcome.warnings).with_sha256(checksum)
}
