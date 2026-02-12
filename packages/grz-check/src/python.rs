use pyo3::prelude::*;
use std::io::{BufReader, Read, Result as IoResult};
use std::path::PathBuf;

use crate::checker::FileReport;
use crate::checks::fastq::{check_single_fastq, process_paired_readers, ReadLengthCheck};
use crate::checks::bam::check_bam;

/// Adapter to use Python file-like objects as Rust Read trait
struct PyFileRead {
    obj: PyObject,
}

impl PyFileRead {
    fn new(obj: PyObject) -> Self {
        Self { obj }
    }
}

impl Read for PyFileRead {
    fn read(&mut self, buf: &mut [u8]) -> IoResult<usize> {
        Python::with_gil(|py| {
            let bytes = self.obj.call_method1(py, "read", (buf.len(),))?;
            let bytes: &[u8] = bytes.extract(py)?;
            buf[..bytes.len()].copy_from_slice(bytes);
            Ok(bytes.len())
        }).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e.to_string()))
    }
}

/// Validation report exposed to Python
#[pyclass(name = "ValidationReport")]
#[derive(Clone)]
pub struct PyValidationReport {
    #[pyo3(get)]
    pub path: String,
    #[pyo3(get)]
    pub num_records: Option<u64>,
    #[pyo3(get)]
    pub mean_read_length: Option<f64>,
    #[pyo3(get)]
    pub sha256: Option<String>,
    #[pyo3(get)]
    pub errors: Vec<String>,
    #[pyo3(get)]
    pub warnings: Vec<String>,
    #[pyo3(get)]
    pub is_valid: bool,
}

impl From<FileReport> for PyValidationReport {
    fn from(report: FileReport) -> Self {
        let is_valid = report.is_ok();
        Self {
            path: report.path.to_string_lossy().to_string(),
            num_records: report.stats.map(|s| s.num_records),
            mean_read_length: report.stats.and_then(|s| s.mean_read_length()),
            sha256: report.sha256,
            errors: report.errors,
            warnings: report.warnings,
            is_valid,
        }
    }
}

/// Helper to parse read length check from Python input
fn parse_read_length_check(min_mean_read_length: Option<i64>) -> ReadLengthCheck {
    match min_mean_read_length {
        None => ReadLengthCheck::Skip,
        Some(n) if n < 0 => ReadLengthCheck::Skip,
        Some(n) => ReadLengthCheck::Fixed(n as usize),
    }
}

/// Helper to create a no-op progress bar (not used in Python bindings)
fn noop_progress_bar() -> indicatif::ProgressBar {
    indicatif::ProgressBar::hidden()
}

/// Validate a single-end FASTQ file
///
/// Args:
///     path: Path to the FASTQ file
///     min_mean_read_length: Minimum mean read length required (>0), or -1/None to skip check
///
/// Returns:
///     ValidationReport with validation results
#[pyfunction]
#[pyo3(signature = (path, min_mean_read_length=None))]
fn validate_fastq_single(
    py: Python,
    path: &str,
    min_mean_read_length: Option<i64>,
) -> PyResult<PyValidationReport> {
    let path = PathBuf::from(path);
    let length_check = parse_read_length_check(min_mean_read_length);

    // Run validation (may take time, so release GIL)
    let report = py.allow_threads(|| {
        let pb = noop_progress_bar();
        let main_pb = noop_progress_bar();
        check_single_fastq(&path, length_check, &pb, &main_pb)
    });

    Ok(PyValidationReport::from(report))
}

/// Validate paired-end FASTQ from Python file-like objects
///
/// Args:
///     r1: File-like object with .read() method (e.g., io.BytesIO, gzip.open())
///     r2: File-like object with .read() method
///     min_mean_read_length: Minimum mean read length (>0), or -1/None to skip
#[pyfunction]
#[pyo3(signature = (r1, r2, min_mean_read_length=None))]
fn validate_fastq_paired(
    py: Python,
    r1: PyObject,
    r2: PyObject,
    min_mean_read_length: Option<i64>,
) -> PyResult<(PyValidationReport, PyValidationReport)> {
    let length_check = parse_read_length_check(min_mean_read_length);

    let (outcome1, outcome2, _pair_errors) = py.allow_threads(|| {
        let py_reader1 = PyFileRead::new(r1);
        let py_reader2 = PyFileRead::new(r2);
        // Use 32KB buffer to reduce GIL crossings
        let reader1 = BufReader::with_capacity(32 * 1024, py_reader1);
        let reader2 = BufReader::with_capacity(32 * 1024, py_reader2);
        process_paired_readers(reader1, reader2, length_check)
    }).map_err(|e: String| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    // Build reports (no path/sha256 for streams)
    let r1_report = PyValidationReport {
        path: "".to_string(),
        num_records: outcome1.stats.map(|s| s.num_records),
        mean_read_length: outcome1.stats.and_then(|s| s.mean_read_length()),
        sha256: None,
        errors: outcome1.errors.clone(),
        warnings: outcome1.warnings,
        is_valid: outcome1.errors.is_empty(),
    };

    let r2_report = PyValidationReport {
        path: "".to_string(),
        num_records: outcome2.stats.map(|s| s.num_records),
        mean_read_length: outcome2.stats.and_then(|s| s.mean_read_length()),
        sha256: None,
        errors: outcome2.errors.clone(),
        warnings: outcome2.warnings,
        is_valid: outcome2.errors.is_empty(),
    };

    Ok((r1_report, r2_report))
}

/// Validate a BAM file
///
/// Args:
///     path: Path to the BAM file
///
/// Returns:
///     ValidationReport with validation results
#[pyfunction]
fn validate_bam(py: Python, path: &str) -> PyResult<PyValidationReport> {
    let path = PathBuf::from(path);

    let report = py.allow_threads(|| {
        let pb = noop_progress_bar();
        let main_pb = noop_progress_bar();
        check_bam(&path, &pb, &main_pb)
    });

    Ok(PyValidationReport::from(report))
}

/// Calculate SHA256 checksum of a file
///
/// Args:
///     path: Path to the file
///
/// Returns:
///     Hex-encoded SHA256 checksum string
#[pyfunction]
fn calculate_checksum(py: Python, path: &str) -> PyResult<String> {
    use crate::checks::raw::check_raw;

    let path = PathBuf::from(path);

    let report = py.allow_threads(|| {
        let pb = noop_progress_bar();
        let main_pb = noop_progress_bar();
        check_raw(&path, &pb, &main_pb)
    });

    report.sha256.ok_or_else(|| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("Failed to compute checksum")
    })
}

/// The grz_check Python module
#[pymodule]
fn grz_check(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyValidationReport>()?;
    m.add_function(wrap_pyfunction!(validate_fastq_single, m)?)?;
    m.add_function(wrap_pyfunction!(validate_fastq_paired, m)?)?;
    m.add_function(wrap_pyfunction!(validate_bam, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_checksum, m)?)?;
    Ok(())
}
