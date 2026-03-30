use pyo3::prelude::*;
use pyo3::types::PyBytes;
use sha2::{Digest, Sha256};
use std::io::{self, BufReader, Read};
use std::path::PathBuf;

use crate::checker::FileReport;
use crate::checks::bam::validate_bam_data;
use crate::checks::common::CheckOutcome;
use crate::checks::fastq::{
    check_single_fastq, process_paired_readers, validate_fastq_data, ReadLengthCheck,
};

/// Wrapper for Python file-like objects that implements Rust's Read trait.
/// This is a custom implementation compatible with PyO3 0.23's Bound API.
pub struct PyFileLikeObject {
    obj: PyObject,
}

impl PyFileLikeObject {
    /// Create a new PyFileLikeObject from a PyObject.
    /// The object must have a .read() method.
    pub fn new(obj: &Bound<'_ , PyAny>) -> PyResult<Self> {
        // Verify the object has a read method
        if !obj.hasattr("read")? {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected object with .read() method",
            ));
        }
        Ok(Self {
            obj: obj.clone().into_pyobject(obj.py())?.into(),
        })
    }
}

impl Read for PyFileLikeObject {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        Python::with_gil(|py| -> PyResult<usize> {
            let result = self.obj.call_method1(py, "read", (buf.len(),))?;
            let bytes = result.bind(py).downcast::<PyBytes>()?;
            let data = bytes.as_bytes();
            let len = data.len();
            buf[..len].copy_from_slice(data);
            Ok(len)
        })
        .map_err(|e: PyErr| io::Error::new(io::ErrorKind::Other, e.to_string()))
    }
}

impl<'py> FromPyObject<'py> for PyFileLikeObject {
    fn extract_bound(obj: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::new(obj)
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

impl From<CheckOutcome> for PyValidationReport {
    fn from(outcome: CheckOutcome) -> Self {
        let is_valid = outcome.errors.is_empty();
        Self {
            path: "".to_string(),
            num_records: outcome.stats.map(|s| s.num_records),
            mean_read_length: outcome.stats.and_then(|s| s.mean_read_length()),
            sha256: None,
            errors: outcome.errors,
            warnings: outcome.warnings,
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

/// Validate a single-end FASTQ file from a path
///
/// This is a thin wrapper around the Rust file-based validation.
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

    let report = py.allow_threads(|| {
        let pb = noop_progress_bar();
        let main_pb = noop_progress_bar();
        check_single_fastq(&path, length_check, &pb, &main_pb)
    });

    Ok(PyValidationReport::from(report))
}

/// Validate a single-end FASTQ from a Python file-like object
///
/// This is a thin wrapper around `validate_fastq_data` from the Rust library.
///
/// Args:
///     source: File-like object with .read() method (e.g., io.BytesIO, gzip.open(), open())
///     min_mean_read_length: Minimum mean read length required (>0), or -1/None to skip check
///
/// Returns:
///     ValidationReport with validation results
#[pyfunction]
#[pyo3(signature = (source, min_mean_read_length=None))]
fn validate_fastq_single_stream(
    py: Python,
    source: PyFileLikeObject,
    min_mean_read_length: Option<i64>,
) -> PyResult<PyValidationReport> {
    let length_check = parse_read_length_check(min_mean_read_length);

    let outcome = py
        .allow_threads(|| validate_fastq_data(source, length_check))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    Ok(PyValidationReport::from(outcome))
}

/// Validate paired-end FASTQ from Python file-like objects (stream-based)
///
/// This is a thin wrapper around `process_paired_readers` from the Rust library.
///
/// Args:
///     r1: File-like object with .read() method (e.g., io.BytesIO, gzip.open())
///     r2: File-like object with .read() method
///     min_mean_read_length: Minimum mean read length (>0), or -1/None to skip
#[pyfunction]
#[pyo3(signature = (r1, r2, min_mean_read_length=None))]
fn validate_fastq_paired_stream(
    py: Python,
    r1: PyFileLikeObject,
    r2: PyFileLikeObject,
    min_mean_read_length: Option<i64>,
) -> PyResult<(PyValidationReport, PyValidationReport)> {
    let length_check = parse_read_length_check(min_mean_read_length);

    let (outcome1, outcome2, _pair_errors) = py
        .allow_threads(|| {
            // Use 128KB buffer for better performance
            let reader1 = BufReader::with_capacity(128 * 1024, r1);
            let reader2 = BufReader::with_capacity(128 * 1024, r2);
            process_paired_readers(reader1, reader2, length_check)
        })
        .map_err(|e: String| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    let r1_report = PyValidationReport::from(outcome1);
    let r2_report = PyValidationReport::from(outcome2);

    Ok((r1_report, r2_report))
}

/// Validate paired-end FASTQ files from paths
///
/// This opens the files and validates them using the same logic as `validate_fastq_paired`.
///
/// Args:
///     r1_path: Path to the first FASTQ file (R1)
///     r2_path: Path to the second FASTQ file (R2)
///     min_mean_read_length: Minimum mean read length (>0), or -1/None to skip
///
/// Returns:
///     Tuple of (r1_report, r2_report) ValidationReports
#[pyfunction]
#[pyo3(signature = (r1_path, r2_path, min_mean_read_length=None))]
fn validate_fastq_paired_paths(
    py: Python,
    r1_path: &str,
    r2_path: &str,
    min_mean_read_length: Option<i64>,
) -> PyResult<(PyValidationReport, PyValidationReport)> {
    use std::fs::File;

    let length_check = parse_read_length_check(min_mean_read_length);

    let (outcome1, outcome2, _pair_errors) = py
        .allow_threads(|| {
            let file1 = File::open(r1_path)
                .map_err(|e| format!("Failed to open R1 file: {e}"))?;
            let file2 = File::open(r2_path)
                .map_err(|e| format!("Failed to open R2 file: {e}"))?;

            // Auto-detect and decompress (e.g. .gz files) using niffler
            let (reader1, _fmt1) = niffler::get_reader(Box::new(BufReader::with_capacity(128 * 1024, file1)))
                .map_err(|e| format!("Failed to decompress R1 file: {e}"))?;
            let (reader2, _fmt2) = niffler::get_reader(Box::new(BufReader::with_capacity(128 * 1024, file2)))
                .map_err(|e| format!("Failed to decompress R2 file: {e}"))?;

            process_paired_readers(reader1, reader2, length_check)
        })
        .map_err(|e: String| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    let r1_report = PyValidationReport::from(outcome1);
    let r2_report = PyValidationReport::from(outcome2);

    Ok((r1_report, r2_report))
}

/// Validate a BAM file from a path
///
/// This is a thin wrapper around the Rust file-based validation.
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
        crate::checks::bam::check_bam(&path, &pb, &main_pb)
    });

    Ok(PyValidationReport::from(report))
}

/// Validate a BAM file from a Python file-like object
///
/// This is a thin wrapper around `validate_bam_data` from the Rust library.
///
/// Args:
///     source: File-like object with .read() method (e.g., io.BytesIO, open())
///
/// Returns:
///     ValidationReport with validation results
#[pyfunction]
fn validate_bam_stream(py: Python, source: PyFileLikeObject) -> PyResult<PyValidationReport> {
    let outcome = py
        .allow_threads(|| validate_bam_data(source))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    Ok(PyValidationReport::from(outcome))
}

/// Internal helper to calculate checksum from a path
fn calculate_file_checksum(path: &std::path::Path) -> Result<String, std::io::Error> {
    use std::fs::File;

    let file = File::open(path)?;
    let mut reader = BufReader::new(file);
    let mut hasher = Sha256::new();

    std::io::copy(&mut reader, &mut hasher)?;

    Ok(format!("{:x}", hasher.finalize()))
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
    let path = PathBuf::from(path);

    let checksum = py
        .allow_threads(|| calculate_file_checksum(&path))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    Ok(checksum)
}

/// Internal helper to calculate checksum from any Read source
fn calculate_stream_checksum<R: Read>(mut reader: R) -> Result<String, std::io::Error> {
    let mut hasher = Sha256::new();
    std::io::copy(&mut reader, &mut hasher)?;
    Ok(format!("{:x}", hasher.finalize()))
}

/// Calculate SHA256 checksum from a Python file-like object
///
/// Args:
///     source: File-like object with .read() method
///
/// Returns:
///     Hex-encoded SHA256 checksum string
#[pyfunction]
fn calculate_checksum_stream(py: Python, source: PyFileLikeObject) -> PyResult<String> {
    let checksum = py
        .allow_threads(|| calculate_stream_checksum(source))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    Ok(checksum)
}

/// The grz_check Python module
#[pymodule]
fn grz_check(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyValidationReport>()?;
    m.add_function(wrap_pyfunction!(validate_fastq_single, m)?)?;
    m.add_function(wrap_pyfunction!(validate_fastq_single_stream, m)?)?;
    m.add_function(wrap_pyfunction!(validate_fastq_paired_stream, m)?)?;
    m.add_function(wrap_pyfunction!(validate_fastq_paired_paths, m)?)?;
    m.add_function(wrap_pyfunction!(validate_bam, m)?)?;
    m.add_function(wrap_pyfunction!(validate_bam_stream, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_checksum, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_checksum_stream, m)?)?;
    Ok(())
}