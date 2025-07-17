use crate::checker::{FileReport, Stats};
use crate::checks::common::{CheckOutcome, check_file};
use indicatif::ProgressBar;
use std::io::BufReader;
use std::path::{Path, PathBuf};

#[derive(Debug, Copy, Clone)]
pub enum ReadLengthCheck {
    Fixed(usize),
    Auto,
    Skip,
}

pub fn check_single_fastq(
    path: &Path,
    length_check: ReadLengthCheck,
    file_pb: &ProgressBar,
    global_pb: &ProgressBar,
) -> FileReport {
    check_file(path, file_pb, global_pb, true, |reader| {
        let mut fastq_reader = noodles::fastq::io::Reader::new(BufReader::new(reader));
        let mut records = fastq_reader.records();

        let first_record = match records.next() {
            Some(Ok(rec)) => rec,
            Some(Err(e)) => return Err(format!("Failed to parse first record: {e}")),
            None => {
                return Ok(CheckOutcome {
                    errors: vec!["File is empty. Expected at least one record.".to_string()],
                    ..Default::default()
                });
            }
        };

        let file_read_length = first_record.sequence().len();
        let mut num_reads: u64 = 1;

        if let ReadLengthCheck::Fixed(expected) = length_check {
            if expected != file_read_length {
                return Ok(CheckOutcome {
                    stats: Some(Stats {
                        num_records: num_reads,
                        read_length: Some(file_read_length),
                    }),
                    errors: vec![format!(
                        "Provided read length ({}) does not match first read's length ({}).",
                        expected, file_read_length
                    )],
                    warnings: vec![],
                });
            }
        }

        if !matches!(length_check, ReadLengthCheck::Skip) {
            for (i, record) in records.enumerate() {
                num_reads += 1;
                let record = match record {
                    Ok(rec) => rec,
                    Err(e) => return Err(format!("Failed to parse record #{i}: {e}")),
                };

                if record.sequence().len() != file_read_length {
                    return Ok(CheckOutcome {
                        stats: Some(Stats {
                            num_records: num_reads,
                            read_length: Some(file_read_length),
                        }),
                        errors: vec![format!(
                            "Found inconsistent read length at record #{}. Expected {}, but got {}.",
                            i,
                            file_read_length,
                            record.sequence().len()
                        )],
                        warnings: vec![],
                    });
                }
            }
        } else {
            num_reads += records.count() as u64;
        }

        Ok(CheckOutcome {
            stats: Some(Stats {
                num_records: num_reads,
                read_length: Some(file_read_length),
            }),
            errors: vec![],
            warnings: vec![],
        })
    })
}

#[derive(Debug)]
pub struct FastqCheckJob {
    pub fq1: PathBuf,
    pub fq2: Option<PathBuf>,
    pub fq1_length_check: ReadLengthCheck,
    pub fq2_length_check: Option<ReadLengthCheck>,
    pub fq1_size: u64,
    pub fq2_size: Option<u64>,
}
