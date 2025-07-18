use crate::checker::{FileReport, Stats};
use crate::checks::common::{CheckOutcome, check_file};
use indicatif::ProgressBar;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use std::io::BufReader;
use std::path::{Path, PathBuf};

pub fn check_bam(path: &Path, file_pb: &ProgressBar, global_pb: &ProgressBar) -> FileReport {
    check_file(path, file_pb, global_pb, false, |reader| {
        let mut bam_reader = bam::io::Reader::new(BufReader::new(reader));
        let header = match bam_reader.read_header() {
            Ok(h) => h,
            Err(e) => return Err(format!("Failed to read BAM header: {e}")),
        };

        let mut warnings = Vec::new();
        if !header.reference_sequences().is_empty()
            || !header.read_groups().is_empty()
            || (header.programs().roots().count() != 0)
            || !header.comments().is_empty()
        {
            warnings.push(
                "Detected a header in BAM file, ensure it contains no private information!"
                    .to_string(),
            );
        }

        let mut num_records = 0;
        let mut secondary_alignment_count: u64 = 0;
        let mut first_secondary_warning_details: Option<(u64, String)> = None;
        let mut hard_clip_count: u64 = 0;
        let mut first_hard_clip_warning_details: Option<(u64, String)> = None;

        for (i, result) in bam_reader.records().enumerate() {
            let record = match result {
                Ok(rec) => rec,
                Err(e) => return Err(format!("Failed to parse record #{}: {}", i + 1, e)),
            };
            num_records += 1;

            if record.flags().is_secondary() {
                secondary_alignment_count += 1;
                if first_secondary_warning_details.is_none() {
                    first_secondary_warning_details = Some((
                        num_records,
                        record.name().map(|n| n.to_string()).unwrap_or_default(),
                    ));
                }
            }

            if !record.flags().is_secondary()
                && record
                    .cigar()
                    .iter()
                    .any(|op| op.is_ok_and(|op| op.kind() == Kind::HardClip))
            {
                hard_clip_count += 1;
                if first_hard_clip_warning_details.is_none() {
                    first_hard_clip_warning_details = Some((
                        num_records,
                        record.name().map(|n| n.to_string()).unwrap_or_default(),
                    ));
                }
            }
        }

        if num_records == 0 {
            return Ok(CheckOutcome {
                errors: vec!["File is empty. Expected at least one record.".to_string()],
                ..Default::default()
            });
        }

        if let Some((rec_num, read_name)) = first_secondary_warning_details {
            warnings.push(format!(
                "File contains {secondary_alignment_count} secondary alignment(s). First detected at record #{rec_num} ('{read_name}')."
            ));
        }

        if let Some((rec_num, read_name)) = first_hard_clip_warning_details {
            warnings.push(format!(
                "File contains {hard_clip_count} primary alignment(s) with hard-clipped bases. First detected at record #{rec_num} ('{read_name}')."
            ));
        }

        Ok(CheckOutcome {
            stats: Some(Stats {
                num_records,
                read_length: None,
            }),
            errors: vec![],
            warnings,
        })
    })
}

#[derive(Debug)]
pub struct BamCheckJob {
    pub path: PathBuf,
    pub size: u64,
}
