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

// Gzip magic bytes
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];

#[derive(Debug, Default)]
pub struct CheckOutcome {
    pub stats: Option<Stats>,
    pub errors: Vec<String>,
    pub warnings: Vec<String>,
}

type ReaderAndHasher = (Box<dyn Read>, Arc<Mutex<Sha256>>);

/// Check if a file has gzip magic bytes.
/// Returns an error message if the file has a .gz extension but doesn't have gzip magic bytes.
fn check_gzip_magic(path: &Path) -> Option<String> {
    // Only check files with .gz extension
    if !path.extension().map_or(false, |ext| ext == "gz") {
        return None;
    }

    // Try to read the first 2 bytes
    let mut file = match fs::File::open(path) {
        Ok(f) => f,
        Err(e) => return Some(format!("Failed to open file for magic bytes check: {}", e)),
    };

    let mut magic_bytes = [0u8; 2];
    match file.read_exact(&mut magic_bytes) {
        Ok(_) => {
            if magic_bytes != GZIP_MAGIC {
                return Some(format!(
                    "File has .gz extension but does not have gzip magic bytes (expected 0x1f 0x8b, found 0x{:02x} 0x{:02x})",
                    magic_bytes[0], magic_bytes[1]
                ));
            }
        }
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => {
            return Some("File has .gz extension but is too small to contain gzip header".to_string());
        }
        Err(e) => {
            return Some(format!("Failed to read magic bytes: {}", e));
        }
    }

    None
}

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

    // Check gzip magic bytes if file has .gz extension
    if let Some(error_msg) = check_gzip_magic(path) {
        anyhow::bail!(error_msg);
    }

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
