use crate::checker::{FileReport, Stats};
use crate::checks::common::{CheckOutcome, check_file};
use indicatif::ProgressBar;
use itertools::EitherOrBoth::{Both, Left, Right};
use itertools::Itertools;
use noodles::fastq;
use std::io::BufRead;
use std::path::{Path, PathBuf};

#[derive(Debug, Copy, Clone)]
pub enum ReadLengthCheck {
    Fixed(usize),
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
    pub length_check: ReadLengthCheck,
    pub fq1_size: u64,
    pub fq2_size: u64,
}

struct FastqCheckProcessor {
    length_check: ReadLengthCheck,
    num_records: u64,
    total_read_length: u64,
    errors: Vec<String>,
}

impl FastqCheckProcessor {
    fn new(length_check: ReadLengthCheck) -> Self {
        Self {
            length_check,
            num_records: 0,
            total_read_length: 0,
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

        if record.sequence().len() != record.quality_scores().len() {
            return Err(format!(
                "Failed to parse {} record #{}: sequence length ({}) does not match quality scores length ({})",
                file_id,
                self.num_records,
                record.sequence().len(),
                record.quality_scores().len()
            ));
        }

        self.total_read_length = self
            .total_read_length
            .checked_add(
                u64::try_from(record.sequence().len())
                    .expect("Single FASTQ record length should fit in u64"),
            )
            .expect("Total length of all reads should fit in u64");

        Ok(())
    }

    fn finalize(mut self) -> CheckOutcome {
        if self.num_records == 0 && self.is_ok() {
            self.errors
                .push("File is empty. Expected at least one record.".to_string());
        }

        let mean_read_length = (self.total_read_length as f64) / (self.num_records as f64);

        match self.length_check {
            ReadLengthCheck::Fixed(min_mean_read_length) => {
                // if mean_read_length is NaN (num_records is zero) then following conditional will
                // be false and the error correctly not reported, since the empty file error was
                // already recorded above.
                if mean_read_length < (min_mean_read_length as f64) {
                    self.errors.push(format!(
                        "Mean read length ({}) is less than minimum required ({})",
                        mean_read_length, min_mean_read_length
                    ))
                }
            }
            ReadLengthCheck::Skip => (),
        };

        CheckOutcome {
            stats: if self.num_records > 0 {
                Some(Stats {
                    num_records: self.num_records,
                    total_read_length: Some(self.total_read_length),
                })
            } else {
                None
            },
            errors: self.errors,
            warnings: vec![],
        }
    }
}

/// Validate FASTQ data from any BufRead source.
/// This is the core validation logic, independent of file I/O.
pub fn validate_fastq_data<R: BufRead>(
    reader: R,
    length_check: ReadLengthCheck,
) -> Result<CheckOutcome, String> {
    let mut fastq_reader = fastq::io::Reader::new(reader);
    let mut processor = FastqCheckProcessor::new(length_check);

    for record_res in fastq_reader.records() {
        processor.process_record(record_res, "record")?;
        if !processor.is_ok() {
            break;
        }
    }

    Ok(processor.finalize())
}

pub fn check_single_fastq(
    path: &Path,
    length_check: ReadLengthCheck,
    file_pb: &ProgressBar,
    global_pb: &ProgressBar,
) -> FileReport {
    check_file(path, file_pb, global_pb, true, |reader| {
        validate_fastq_data(reader, length_check)
    })
}

pub fn process_paired_readers<R1, R2>(
    reader1: R1,
    reader2: R2,
    length_check: ReadLengthCheck,
) -> Result<(CheckOutcome, CheckOutcome, Vec<String>), String>
where
    R1: BufRead,
    R2: BufRead,
{
    let mut fq1_reader = fastq::io::Reader::new(reader1);
    let mut fq2_reader = fastq::io::Reader::new(reader2);

    let mut fq1_processor = FastqCheckProcessor::new(length_check);
    let mut fq2_processor = FastqCheckProcessor::new(length_check);
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

#[cfg(test)]
mod tests {
    use super::*;

    const VALID_FASTQ: &[u8] = b"@read1\nACGT\n+\n!!!!\n";
    const INVALID_FASTQ_LEN: &[u8] = b"@read1\nACGT\n+\n!!!\n";
    const VALID_FASTQ_2_RECORDS: &[u8] = b"@read1\nACGT\n+\n!!!!\n@read2\nACGT\n+\n!!!!\n";

    #[test]
    fn test_valid_single_fastq() {
        let outcome = validate_fastq_data(VALID_FASTQ, ReadLengthCheck::Skip).unwrap();

        assert!(outcome.errors.is_empty());
        let stats = outcome
            .stats
            .expect("Stats should be present for non-empty file");
        assert_eq!(stats.num_records, 1);
        assert_eq!(stats.total_read_length, Some(4));
    }

    #[test]
    fn test_mismatched_sequence_and_quality() {
        let err = validate_fastq_data(INVALID_FASTQ_LEN, ReadLengthCheck::Skip).unwrap_err();

        assert!(err.contains("sequence length (4) does not match quality scores length (3)"));
    }

    #[test]
    fn test_mean_read_length_pass() {
        let outcome = validate_fastq_data(VALID_FASTQ, ReadLengthCheck::Fixed(4)).unwrap();
        assert!(outcome.errors.is_empty());
    }

    #[test]
    fn test_mean_read_length_fail() {
        let outcome = validate_fastq_data(VALID_FASTQ, ReadLengthCheck::Fixed(5)).unwrap();

        assert_eq!(
            outcome.errors,
            vec!["Mean read length (4) is less than minimum required (5)"]
        );
    }

    #[test]
    fn test_invalid_fastq_syntax() {
        let invalid_syntax: &[u8] = b"NOT_A_FASTQ_RECORD\nACGT\n+\n!!!!\n";
        let res = validate_fastq_data(invalid_syntax, ReadLengthCheck::Skip);

        assert!(res.is_err());
        assert!(
            res.unwrap_err()
                .contains("Failed to parse record record #1")
        );
    }

    #[test]
    fn test_paired_readers_valid() {
        let (out1, out2, pair_errors) =
            process_paired_readers(VALID_FASTQ, VALID_FASTQ, ReadLengthCheck::Skip).unwrap();

        assert!(out1.errors.is_empty());
        assert!(out2.errors.is_empty());
        assert!(pair_errors.is_empty());

        assert_eq!(out1.stats.unwrap().num_records, 1);
        assert_eq!(out2.stats.unwrap().num_records, 1);
    }

    #[test]
    fn test_paired_readers_r1_longer() {
        let (out1, out2, pair_errors) =
            process_paired_readers(VALID_FASTQ_2_RECORDS, VALID_FASTQ, ReadLengthCheck::Skip)
                .unwrap();

        assert_eq!(
            pair_errors,
            vec!["Mismatched read counts: R1 has more records than R2."]
        );

        assert_eq!(out1.stats.unwrap().num_records, 2);
        assert_eq!(out2.stats.unwrap().num_records, 1);
    }

    #[test]
    fn test_paired_readers_r2_longer() {
        let (_out1, _out2, pair_errors) =
            process_paired_readers(VALID_FASTQ, VALID_FASTQ_2_RECORDS, ReadLengthCheck::Skip)
                .unwrap();

        assert_eq!(
            pair_errors,
            vec!["Mismatched read counts: R2 has more records than R1."]
        );
    }
}
