use pyo3::buffer::PyBuffer;
use pyo3::prelude::*;
use pyo3::types::PyBytes;
use sha2::{Digest, Sha256};
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::PathBuf;

use crate::checker::FileReport;
use crate::checks::bam::validate_bam_data;
use crate::checks::common::CheckOutcome;
use crate::checks::fastq::{
    ReadLengthCheck, check_single_fastq, process_paired_readers, validate_fastq_data,
};

const STREAM_BUF_SIZE: usize = 4 * 1024 * 1024; // 4 MB

/// Wrapper for Python file-like objects that implements Rust's Read trait.
pub struct PyFileLikeObject {
    obj: PyObject,
}

impl PyFileLikeObject {
    /// Create a new PyFileLikeObject from a PyObject.
    /// The object must have a .read() method.
    pub fn new(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        // Verify the object has a read method
        if !obj.hasattr("read")? {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected object with .read() method",
            ));
        }
        // If the object has a string mode attribute (e.g. files opened with open()), verify
        // it is binary. Skip non-string modes: gzip.GzipFile exposes an integer mode, and
        // io.BytesIO has no mode attribute at all — both are valid binary sources.
        if let Ok(mode) = obj.getattr("mode") {
            if let Ok(mode_str) = mode.extract::<String>() {
                if !mode_str.contains('b') {
                    return Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!(
                        "File must be opened in binary mode (e.g. 'rb'), got mode '{mode_str}'"
                    )));
                }
            }
        }
        Ok(Self {
            obj: obj.clone().unbind(),
        })
    }
}

impl Read for PyFileLikeObject {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        Python::with_gil(|py| -> PyResult<usize> {
            let result = self.obj.bind(py).call_method1("read", (buf.len(),))?;
            let bytes = result.downcast::<PyBytes>()?;
            let data = bytes.as_bytes();
            let len = data.len();
            buf[..len].copy_from_slice(data);
            Ok(len)
        })
        .map_err(|e: PyErr| {
            io::Error::new(
                io::ErrorKind::Other,
                format!("Failed to read from Python file-like object: {e}"),
            )
        })
    }
}

/// Represents a Python input ready to be consumed by the Rust validator.
/// `Send`-safe so it can cross `py.allow_threads` boundaries.
enum StreamInput {
    /// Zero-copy view into a Python buffer-protocol object
    /// (bytes, bytearray, memoryview, mmap, numpy u8 arrays, array.array('B', …)).
    /// The underlying memory is never copied into a Rust `Vec`.
    PyBuf(PyBuffer<u8>),
    /// File-like object with .read() method.
    FileLike(PyFileLikeObject),
}

/// Try to extract a Python object as a buffer-protocol object, falling back to file-like.
/// Buffer-protocol inputs are held as a zero-copy `PyBuffer<u8>` view; the Rust side
/// reads directly from Python-owned memory via a `Cursor` over the buffer's slice.
fn extract_stream_input(obj: &Bound<'_, PyAny>) -> PyResult<StreamInput> {
    if let Ok(buffer) = PyBuffer::<u8>::get(obj) {
        if !buffer.is_c_contiguous() {
            return Err(PyErr::new::<pyo3::exceptions::PyBufferError, _>(
                "Buffer must be C-contiguous (flat, no strides); strided or multi-dimensional views are not supported",
            ));
        }
        return Ok(StreamInput::PyBuf(buffer));
    }
    Ok(StreamInput::FileLike(PyFileLikeObject::new(obj)?))
}

/// Borrow a `PyBuffer<u8>`'s memory as a plain byte slice.
///
/// SAFETY: The returned slice is only valid while `buffer` is alive and the
/// underlying Python memory is not mutated. `PyBuffer` holds a `Py_buffer`
/// view that pins the memory (Python cannot resize/realloc the object while
/// the view is held). Inside `py.allow_threads` the GIL is released, so
/// Python-level mutation is blocked; callers passing `bytes`, `memoryview`,
/// or `mmap(ACCESS_READ)` get an effectively immutable buffer. Writable
/// buffers (bytearray, mmap ACCESS_WRITE) must not be concurrently mutated
/// from other OS threads for the duration of the borrow.
unsafe fn buffer_as_slice(buffer: &PyBuffer<u8>) -> &[u8] {
    // SAFETY: see function doc. `PyBuffer<u8>` is constructed via PyBuffer::<u8>::get
    // which guarantees itemsize == 1 and format 'B'; `is_c_contiguous` is checked at
    // extraction time, so the buffer is a flat u8 run of `item_count()` bytes.
    unsafe { std::slice::from_raw_parts(buffer.buf_ptr() as *const u8, buffer.item_count()) }
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

/// Run `f` with a `BufRead` over `input`, applying niffler auto-decompression
/// (gzip / bgzip / bzip2 / xz / zstd) transparently.
///
/// For `PyBuf` inputs this is zero-copy: the Rust side constructs a
/// `Cursor<&[u8]>` directly over Python-owned memory — no bytes are copied
/// into a `Vec` first. A single 4 MB outer `BufReader` is kept for noodles'
/// line-based FASTQ parser.
fn with_fastq_reader<F, T>(input: StreamInput, f: F) -> Result<T, String>
where
    F: FnOnce(&mut dyn BufRead) -> Result<T, String>,
{
    match input {
        StreamInput::PyBuf(buffer) => {
            // SAFETY: `buffer` outlives `reader` below (local binding order);
            // see `buffer_as_slice` doc for the non-mutation contract.
            let slice = unsafe { buffer_as_slice(&buffer) };
            let (decompressed, _format) = niffler::get_reader(Box::new(Cursor::new(slice)))
                .map_err(|e| format!("Failed to detect or decompress stream format: {e}"))?;
            let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, decompressed);
            let result = f(&mut reader);
            drop(reader);
            drop(buffer);
            result
        }
        StreamInput::FileLike(file_like) => {
            let (decompressed, _format) = niffler::get_reader(Box::new(file_like))
                .map_err(|e| format!("Failed to detect or decompress stream format: {e}"))?;
            let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, decompressed);
            f(&mut reader)
        }
    }
}

/// Validate a single-end FASTQ file from a path
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

/// Validate a single-end FASTQ from a Python stream or buffer.
///
/// Accepts bytes / bytearray / memoryview / mmap (zero-copy, no memory copy
/// into Rust) or any file-like object with .read(). Compression (gzip, bgzip)
/// is detected and decompressed transparently on the Rust side — no need to
/// wrap in gzip.open().
#[pyfunction]
#[pyo3(signature = (source, min_mean_read_length=None))]
fn validate_fastq_single_stream(
    py: Python,
    source: &Bound<'_, PyAny>,
    min_mean_read_length: Option<i64>,
) -> PyResult<PyValidationReport> {
    let length_check = parse_read_length_check(min_mean_read_length);
    let input = extract_stream_input(source)?;

    let outcome = py
        .allow_threads(|| {
            with_fastq_reader(input, |reader| validate_fastq_data(reader, length_check))
        })
        .map_err(|e: String| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    Ok(PyValidationReport::from(outcome))
}

/// Validate paired-end FASTQ from Python streams or buffers.
///
/// Accepts bytes / bytearray / memoryview / mmap (zero-copy, no memory copy
/// into Rust) or any file-like object with .read(). Compression (gzip, bgzip)
/// is detected and decompressed transparently on the Rust side — no need to
/// wrap in gzip.open().
#[pyfunction]
#[pyo3(signature = (r1, r2, min_mean_read_length=None))]
fn validate_fastq_paired_stream(
    py: Python,
    r1: &Bound<'_, PyAny>,
    r2: &Bound<'_, PyAny>,
    min_mean_read_length: Option<i64>,
) -> PyResult<(PyValidationReport, PyValidationReport)> {
    let length_check = parse_read_length_check(min_mean_read_length);
    let input1 = extract_stream_input(r1)?;
    let input2 = extract_stream_input(r2)?;

    let (outcome1, outcome2, _pair_errors) = py
        .allow_threads(|| {
            with_fastq_reader(input1, |reader1| {
                with_fastq_reader(input2, |reader2| {
                    process_paired_readers(reader1, reader2, length_check)
                })
            })
        })
        .map_err(|e: String| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    let r1_report = PyValidationReport::from(outcome1);
    let r2_report = PyValidationReport::from(outcome2);

    Ok((r1_report, r2_report))
}

/// Validate paired-end FASTQ files from paths
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
            let file1 = File::open(r1_path).map_err(|e| format!("Failed to open R1 file: {e}"))?;
            let file2 = File::open(r2_path).map_err(|e| format!("Failed to open R2 file: {e}"))?;

            let (decompressed1, _fmt1) = niffler::get_reader(Box::new(file1))
                .map_err(|e| format!("Failed to decompress R1 file: {e}"))?;
            let (decompressed2, _fmt2) = niffler::get_reader(Box::new(file2))
                .map_err(|e| format!("Failed to decompress R2 file: {e}"))?;

            let reader1 = BufReader::with_capacity(STREAM_BUF_SIZE, decompressed1);
            let reader2 = BufReader::with_capacity(STREAM_BUF_SIZE, decompressed2);

            process_paired_readers(reader1, reader2, length_check)
        })
        .map_err(|e: String| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    let r1_report = PyValidationReport::from(outcome1);
    let r2_report = PyValidationReport::from(outcome2);

    Ok((r1_report, r2_report))
}

/// Validate a BAM file from a path
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

/// Validate a BAM file from a Python stream or buffer.
///
/// Accepts bytes / bytearray / memoryview / mmap (zero-copy, no memory copy
/// into Rust) or any file-like object with .read(). Noodles handles bgzf
/// decompression internally.
#[pyfunction]
fn validate_bam_stream(py: Python, source: &Bound<'_, PyAny>) -> PyResult<PyValidationReport> {
    let input = extract_stream_input(source)?;

    let outcome = py
        .allow_threads(|| match input {
            StreamInput::PyBuf(buffer) => {
                // SAFETY: `buffer` lives through the end of this arm, so the
                // slice (and Cursor borrowing it) remain valid for validate_bam_data.
                let slice = unsafe { buffer_as_slice(&buffer) };
                let result = validate_bam_data(&mut Cursor::new(slice));
                drop(buffer);
                result
            }
            StreamInput::FileLike(file_like) => {
                let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, file_like);
                validate_bam_data(&mut reader)
            }
        })
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e))?;

    Ok(PyValidationReport::from(outcome))
}

/// Internal helper to calculate checksum from a path
fn calculate_file_checksum(path: &std::path::Path) -> Result<String, std::io::Error> {
    use std::fs::File;

    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, file);
    let mut hasher = Sha256::new();

    std::io::copy(&mut reader, &mut hasher)?;

    Ok(format!("{:x}", hasher.finalize()))
}

/// Calculate SHA256 checksum of a file
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

/// Calculate SHA256 checksum from a Python stream or buffer.
///
/// Accepts bytes / bytearray / memoryview / mmap (zero-copy, no memory copy
/// into Rust) or any file-like object with .read(). Computes checksum on raw
/// (uncompressed) bytes.
#[pyfunction]
fn calculate_checksum_stream(py: Python, source: &Bound<'_, PyAny>) -> PyResult<String> {
    let input = extract_stream_input(source)?;

    let checksum = py
        .allow_threads(|| match input {
            StreamInput::PyBuf(buffer) => {
                // SAFETY: `buffer` outlives the cursor (local drop order).
                let slice = unsafe { buffer_as_slice(&buffer) };
                let result = calculate_stream_checksum(Cursor::new(slice));
                drop(buffer);
                result
            }
            StreamInput::FileLike(file_like) => {
                let reader = BufReader::with_capacity(STREAM_BUF_SIZE, file_like);
                calculate_stream_checksum(reader)
            }
        })
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
