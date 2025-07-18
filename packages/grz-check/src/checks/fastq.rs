use crate::checker::{FileReport, Stats};
use crate::checks::common::{CheckOutcome, check_file};
use indicatif::ProgressBar;
use itertools::EitherOrBoth::{Both, Left, Right};
use itertools::Itertools;
use noodles::fastq;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};

#[derive(Debug, Copy, Clone)]
pub enum ReadLengthCheck {
    Fixed(usize),
    Auto,
    Skip,
}

#[derive(Debug)]
pub struct SingleFastqJob {
    pub path: PathBuf,
    pub length_check: ReadLengthCheck,
    pub size: u64,
}

#[derive(Debug)]
pub struct PairedFastqJob {
    pub fq1_path: PathBuf,
    pub fq2_path: PathBuf,
    pub fq1_length_check: ReadLengthCheck,
    pub fq2_length_check: ReadLengthCheck,
    pub fq1_size: u64,
    pub fq2_size: u64,
}

struct FastqCheckProcessor {
    length_check: ReadLengthCheck,
    num_records: u64,
    detected_read_length: Option<usize>,
    errors: Vec<String>,
}

impl FastqCheckProcessor {
    fn new(length_check: ReadLengthCheck) -> Self {
        Self {
            length_check,
            num_records: 0,
            detected_read_length: None,
            errors: Vec::new(),
        }
    }

    fn is_ok(&self) -> bool {
        self.errors.is_empty() && self.num_records > 0
    }

    fn process_record(
        &mut self,
        record: Result<fastq::Record, std::io::Error>,
        file_id: &str,
    ) -> Result<(), String> {
        self.num_records += 1;

        let record = record.map_err(|e| {
            format!(
                "Failed to parse {} record #{}: {}",
                file_id, self.num_records, e
            )
        })?;

        let current_len = record.sequence().len();

        if let Some(expected_len) = self.detected_read_length {
            if !matches!(self.length_check, ReadLengthCheck::Skip) && current_len != expected_len {
                self.errors.push(format!(
                    "Found inconsistent read length at record #{}. Expected {}, but got {}.",
                    self.num_records, expected_len, current_len
                ));
            }
        } else {
            self.detected_read_length = Some(current_len);
            if let ReadLengthCheck::Fixed(expected) = self.length_check {
                if current_len != expected {
                    self.errors.push(format!(
                        "Provided read length ({expected}) does not match first read's length ({current_len})."
                    ));
                }
            }
        }
        Ok(())
    }

    fn finalize(mut self) -> CheckOutcome {
        if self.num_records == 0 && self.is_ok() {
            self.errors
                .push("File is empty. Expected at least one record.".to_string());
        }

        CheckOutcome {
            stats: if self.num_records > 0 {
                Some(Stats {
                    num_records: self.num_records,
                    read_length: self.detected_read_length,
                })
            } else {
                None
            },
            errors: self.errors,
            warnings: vec![],
        }
    }
}

pub fn check_single_fastq(
    path: &Path,
    length_check: ReadLengthCheck,
    file_pb: &ProgressBar,
    global_pb: &ProgressBar,
) -> FileReport {
    check_file(path, file_pb, global_pb, true, |reader| {
        let mut fastq_reader = fastq::io::Reader::new(BufReader::new(reader));
        let mut processor = FastqCheckProcessor::new(length_check);

        for record_res in fastq_reader.records() {
            processor.process_record(record_res, "record")?;
            if !processor.is_ok() {
                break;
            }
        }

        Ok(processor.finalize())
    })
}

pub fn process_paired_readers<R1, R2>(
    reader1: R1,
    reader2: R2,
    length_check1: ReadLengthCheck,
    length_check2: ReadLengthCheck,
) -> Result<(CheckOutcome, CheckOutcome, Vec<String>), String>
where
    R1: Read,
    R2: Read,
{
    let mut fq1_reader = fastq::io::Reader::new(BufReader::new(reader1));
    let mut fq2_reader = fastq::io::Reader::new(BufReader::new(reader2));

    let mut fq1_processor = FastqCheckProcessor::new(length_check1);
    let mut fq2_processor = FastqCheckProcessor::new(length_check2);
    let mut pair_errors = Vec::new();

    for result in fq1_reader.records().zip_longest(fq2_reader.records()) {
        match result {
            Both(r1_res, r2_res) => {
                fq1_processor.process_record(r1_res, "R1")?;
                fq2_processor.process_record(r2_res, "R2")?;
            }
            Left(r1_res) => {
                fq1_processor.process_record(r1_res, "R1")?;
                pair_errors
                    .push("Mismatched read counts: R1 has more records than R2.".to_string());
            }
            Right(r2_res) => {
                fq2_processor.process_record(r2_res, "R2")?;
                pair_errors
                    .push("Mismatched read counts: R2 has more records than R1.".to_string());
            }
        }
        if !fq1_processor.is_ok() || !fq2_processor.is_ok() || !pair_errors.is_empty() {
            break;
        }
    }

    let outcome1 = fq1_processor.finalize();
    let outcome2 = fq2_processor.finalize();

    Ok((outcome1, outcome2, pair_errors))
}
