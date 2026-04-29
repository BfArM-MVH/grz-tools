use pyo3::buffer::PyBuffer;
use pyo3::prelude::*;
use sha2::{Digest, Sha256};
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use crate::checker::FileReport;
use crate::checks::bam::{check_bam, validate_bam_data};
use crate::checks::common::{CheckOutcome, READ_BUF_SIZE};
use crate::checks::fastq::{
    ReadLengthCheck, check_single_fastq, process_paired_readers, validate_fastq_data,
};

const STREAM_BUF_SIZE: usize = READ_BUF_SIZE;

/// sha2 dropped std::io::Write implementation, so we work around it.
struct HashWriter<'a, D>(&'a mut D);

impl<'a, D: Digest> std::io::Write for HashWriter<'a, D> {
    #[inline]
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.0.update(buf);
        Ok(buf.len())
    }

    #[inline]
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

/// Wrapper for Python file-like objects that implements Rust's Read trait.
pub struct PyFileLikeObject {
    obj: Py<PyAny>,
}

impl PyFileLikeObject {
    /// Create a new PyFileLikeObject from a Py<PyAny>.
    /// The object must have a .read() method.
    pub fn new(obj: &Bound<'_, PyAny>) -> PyResult<Self> {
        // Verify the object has a read method
        if !obj.hasattr("read")? {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected object with .read() method",
            ));
        }
        // If the object has a string mode attribute (e.g. files opened with open()), verify
        // it is binary.
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

/// Implements the Read trait for Python file-like objects.
///
/// Each Rust read call translates to one Python .read() call and copies
/// the returned bytes into the Rust buffer.
/// For inputs that expose the buffer protocol, the zero-copy     
/// path in `buffer_as_slice` avoids this cost.
impl Read for PyFileLikeObject {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        Python::attach(|py| -> PyResult<usize> {
            let result = self.obj.bind(py).call_method1("read", (buf.len(),))?;
            let data: &[u8] = result.extract()?;
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

enum StreamInput {
    /// Read in place — no copy.                                              
    PyBuf(PyBuffer<u8>),
    /// Read via Python's .read() — copies each chunk.           
    FileLike(PyFileLikeObject),
}

/// Extracts a stream input from a Python object.
///
/// This function is indirectly called from Python.
/// Prefer the fast path (read directly from a buffer) over the slow one (.read() loop).
fn extract_stream_input(obj: &Bound<'_, PyAny>) -> PyResult<StreamInput> {
    if let Ok(buffer) = PyBuffer::<u8>::get(obj) {
        if !buffer.readonly() {
            return Err(PyErr::new::<pyo3::exceptions::PyBufferError, _>(
                "Writable buffers (bytearray, memoryview of bytearray, mmap ACCESS_WRITE) are not supported; \
                 pass bytes, memoryview(bytes), or mmap(..., access=mmap.ACCESS_READ) instead",
            ));
        }
        if !buffer.is_c_contiguous() {
            return Err(PyErr::new::<pyo3::exceptions::PyBufferError, _>(
                "Input buffer must be a contiguous byte sequence; \
                 strided or multi-dimensional views (e.g. arr[::2], arr[:, 0]) are not supported",
            ));
        }
        return Ok(StreamInput::PyBuf(buffer));
    }
    Ok(StreamInput::FileLike(PyFileLikeObject::new(obj)?))
}

/// Return the Python buffer's memory as a plain bytes.
/// The caller must keep the buffer alive for as long as the byte slice is used.
/// In function `extract_stream_input` it is verified that the buffer is read-only
unsafe fn buffer_as_slice(buffer: &PyBuffer<u8>) -> &[u8] {
    let ptr = buffer.buf_ptr() as *const u8;
    let len = buffer.item_count();
    // `buffer` comes from Python, so buf_ptr() can be null (even for length 0).
    if len == 0 || ptr.is_null() {
        return &[];
    }
    unsafe { std::slice::from_raw_parts(ptr, len) }
}

/// Validation report exposed to Python
#[pyclass(name = "ValidationReport", from_py_object)]
#[derive(Clone)]
pub struct PyValidationReport {
    #[pyo3(get, set)]
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

struct StreamHasher<R> {
    inner: R,
    hasher: Arc<Mutex<Sha256>>,
}

impl<R: Read> Read for StreamHasher<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let n = self.inner.read(buf)?;
        if n > 0 {
            self.hasher.lock().unwrap().update(&buf[..n]);
        }
        Ok(n)
    }
}

fn with_fastq_reader<F, T>(input: StreamInput, f: F) -> Result<(T, String), String>
where
    F: FnOnce(&mut dyn BufRead) -> Result<T, String>,
{
    let hasher = Arc::new(Mutex::new(Sha256::new()));

    let result = match input {
        StreamInput::PyBuf(buffer) => {
            // The buffer is dropped after the reader, so the slice stays valid.
            let slice = unsafe { buffer_as_slice(&buffer) };
            let hash_reader = StreamHasher {
                inner: Cursor::new(slice),
                hasher: hasher.clone(),
            };
            let (decompressed, _) = niffler::get_reader(Box::new(hash_reader))
                .map_err(|e| format!("Failed to detect or decompress stream format: {e}"))?;
            let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, decompressed);
            f(&mut reader)?
        }
        StreamInput::FileLike(file_like) => {
            let hash_reader = StreamHasher {
                inner: file_like,
                hasher: hasher.clone(),
            };
            let (decompressed, _) = niffler::get_reader(Box::new(hash_reader))
                .map_err(|e| format!("Failed to detect or decompress stream format: {e}"))?;
            let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, decompressed);
            f(&mut reader)?
        }
    };

    let final_hash = match Arc::try_unwrap(hasher) {
        Ok(mutex) => hex::encode(mutex.into_inner().unwrap().finalize()),
        Err(_) => return Err("Failed to finalize checksum: hasher is still in use.".to_string()),
    };

    Ok((result, final_hash))
}

/// Validate a single-end FASTQ file from a path
#[pyfunction]
#[pyo3(signature = (path, min_mean_read_length=None))]
fn validate_fastq_single(
    py: Python,
    path: PathBuf,
    min_mean_read_length: Option<i64>,
) -> PyResult<PyValidationReport> {
    let length_check = parse_read_length_check(min_mean_read_length);

    let report = py.detach(|| {
        let pb = noop_progress_bar();
        let main_pb = noop_progress_bar();
        check_single_fastq(&path, length_check, &pb, &main_pb)
    });

    Ok(PyValidationReport::from(report))
}

/// Validate a single-end FASTQ from an mmap or file-like object.
#[pyfunction]
#[pyo3(signature = (source, min_mean_read_length=None))]
fn validate_fastq_single_stream(
    py: Python,
    source: &Bound<'_, PyAny>,
    min_mean_read_length: Option<i64>,
) -> PyResult<PyValidationReport> {
    let length_check = parse_read_length_check(min_mean_read_length);
    let input = extract_stream_input(source)?;

    let (outcome, checksum) = py
        .detach(|| -> Result<_, String> {
            with_fastq_reader(input, |reader| validate_fastq_data(reader, length_check))
        })
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)?;

    let mut report = PyValidationReport::from(outcome);
    report.sha256 = Some(checksum);
    Ok(report)
}

/// Validate paired-end FASTQ from mmaps or file-like objects.
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

    let (inner_res, hash1) = py
        .detach(|| -> Result<_, String> {
            with_fastq_reader(input1, |reader1| {
                with_fastq_reader(input2, |reader2| {
                    process_paired_readers(reader1, reader2, length_check)
                })
            })
        })
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)?;

    let ((outcome1, outcome2, pair_errors), hash2) = inner_res;

    let mut r1_report = PyValidationReport::from(outcome1);
    let mut r2_report = PyValidationReport::from(outcome2);

    r1_report.sha256 = Some(hash1);
    r2_report.sha256 = Some(hash2);

    r1_report.errors.extend(pair_errors.clone());
    r2_report.errors.extend(pair_errors.clone());
    r1_report.is_valid &= r1_report.errors.is_empty();
    r2_report.is_valid &= r2_report.errors.is_empty();

    Ok((r1_report, r2_report))
}

/// Helper to set up a buffered reader and a hasher for a given file path.
fn get_file_reader_and_hasher(
    path: &std::path::Path,
) -> Result<(BufReader<Box<dyn Read>>, Arc<Mutex<Sha256>>), String> {
    use std::fs::File;
    let file = File::open(path).map_err(|e| format!("Failed to open {}: {e}", path.display()))?;
    let hasher = Arc::new(Mutex::new(Sha256::new()));
    let hr = StreamHasher {
        inner: file,
        hasher: hasher.clone(),
    };
    let (decompressed, _) = niffler::get_reader(Box::new(hr)).map_err(|e| e.to_string())?;
    Ok((
        BufReader::with_capacity(STREAM_BUF_SIZE, decompressed),
        hasher,
    ))
}

/// Validate paired-end FASTQ files from paths
#[pyfunction]
#[pyo3(signature = (r1_path, r2_path, min_mean_read_length=None))]
fn validate_fastq_paired_paths(
    py: Python,
    r1_path: PathBuf,
    r2_path: PathBuf,
    min_mean_read_length: Option<i64>,
) -> PyResult<(PyValidationReport, PyValidationReport)> {
    let length_check = parse_read_length_check(min_mean_read_length);

    let (outcome1, outcome2, pair_errors, hash1, hash2) = py
        .detach(|| -> Result<_, String> {
            let (mut reader1, hasher1) = get_file_reader_and_hasher(&r1_path)?;
            let (mut reader2, hasher2) = get_file_reader_and_hasher(&r2_path)?;

            let res = process_paired_readers(&mut reader1, &mut reader2, length_check)?;

            // Explicitly drop readers to release the hasher Arcs
            drop(reader1);
            drop(reader2);

            let h1 = hex::encode(
                Arc::try_unwrap(hasher1)
                    .unwrap()
                    .into_inner()
                    .unwrap()
                    .finalize(),
            );
            let h2 = hex::encode(
                Arc::try_unwrap(hasher2)
                    .unwrap()
                    .into_inner()
                    .unwrap()
                    .finalize(),
            );

            Ok((res.0, res.1, res.2, h1, h2))
        })
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)?;

    let mut r1_report = PyValidationReport::from(outcome1);
    let mut r2_report = PyValidationReport::from(outcome2);
    r1_report.sha256 = Some(hash1);
    r2_report.sha256 = Some(hash2);
    r1_report.errors.extend(pair_errors.clone());
    r2_report.errors.extend(pair_errors.clone());
    r1_report.is_valid &= r1_report.errors.is_empty();
    r2_report.is_valid &= r2_report.errors.is_empty();

    Ok((r1_report, r2_report))
}

/// Validate a BAM file from a path
#[pyfunction]
fn validate_bam(py: Python, path: PathBuf) -> PyResult<PyValidationReport> {
    let report = py.detach(|| {
        let pb = noop_progress_bar();
        let main_pb = noop_progress_bar();
        check_bam(&path, &pb, &main_pb)
    });

    Ok(PyValidationReport::from(report))
}

/// Validate a BAM from an mmap or file-like object.
#[pyfunction]
fn validate_bam_stream(py: Python, source: &Bound<'_, PyAny>) -> PyResult<PyValidationReport> {
    let input = extract_stream_input(source)?;

    let (outcome, checksum) = py
        .detach(|| -> Result<_, String> {
            let hasher = Arc::new(Mutex::new(Sha256::new()));

            let res = match input {
                StreamInput::PyBuf(buffer) => {
                    // The buffer is dropped after the cursor, so the slice stays valid.
                    let slice = unsafe { buffer_as_slice(&buffer) };
                    let hash_reader = StreamHasher {
                        inner: Cursor::new(slice),
                        hasher: hasher.clone(),
                    };
                    let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, hash_reader);
                    validate_bam_data(&mut reader)
                }
                StreamInput::FileLike(file_like) => {
                    let hash_reader = StreamHasher {
                        inner: file_like,
                        hasher: hasher.clone(),
                    };
                    let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, hash_reader);
                    validate_bam_data(&mut reader)
                }
            };

            let final_hash = hex::encode(
                Arc::try_unwrap(hasher)
                    .unwrap()
                    .into_inner()
                    .unwrap()
                    .finalize(),
            );
            Ok((res.map_err(|e| e.to_string())?, final_hash))
        })
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)?;

    let mut report = PyValidationReport::from(outcome);
    report.sha256 = Some(checksum);
    Ok(report)
}

#[pyfunction]
fn validate_raw(py: Python, path: PathBuf) -> PyResult<PyValidationReport> {
    let checksum = py
        .detach(|| -> Result<_, String> {
            calculate_file_checksum(&path).map_err(|e| e.to_string())
        })
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)?;

    Ok(PyValidationReport {
        path: path.to_string_lossy().to_string(),
        num_records: None,
        mean_read_length: None,
        sha256: Some(checksum),
        errors: vec![],
        warnings: vec![],
        is_valid: true,
    })
}

#[pyfunction]
fn validate_raw_stream(py: Python, source: &Bound<'_, PyAny>) -> PyResult<PyValidationReport> {
    let input = extract_stream_input(source)?;

    let checksum = py
        .detach(|| -> Result<_, String> {
            match input {
                StreamInput::PyBuf(buffer) => {
                    let slice = unsafe { buffer_as_slice(&buffer) };
                    calculate_stream_checksum(Cursor::new(slice)).map_err(|e| e.to_string())
                }
                StreamInput::FileLike(file_like) => {
                    let reader = BufReader::with_capacity(STREAM_BUF_SIZE, file_like);
                    calculate_stream_checksum(reader).map_err(|e| e.to_string())
                }
            }
        })
        .map_err(PyErr::new::<pyo3::exceptions::PyIOError, _>)?;

    Ok(PyValidationReport {
        path: "".to_string(),
        num_records: None,
        mean_read_length: None,
        sha256: Some(checksum),
        errors: vec![],
        warnings: vec![],
        is_valid: true,
    })
}

fn calculate_file_checksum(path: &std::path::Path) -> Result<String, std::io::Error> {
    use std::fs::File;

    let file = File::open(path)?;
    let mut reader = BufReader::with_capacity(STREAM_BUF_SIZE, file);
    let mut hasher = Sha256::new();

    std::io::copy(&mut reader, &mut HashWriter(&mut hasher))?;
    Ok(hex::encode(hasher.finalize()))
}

/// Calculate SHA256 checksum of a file
#[pyfunction]
fn calculate_checksum(py: Python, path: PathBuf) -> PyResult<String> {
    let checksum = py
        .detach(|| calculate_file_checksum(&path))
        .map_err(|e| PyErr::new::<pyo3::exceptions::PyIOError, _>(e.to_string()))?;

    Ok(checksum)
}

/// Internal helper to calculate checksum from any Read source
fn calculate_stream_checksum<R: Read>(mut reader: R) -> Result<String, std::io::Error> {
    let mut hasher = Sha256::new();
    std::io::copy(&mut reader, &mut HashWriter(&mut hasher))?;
    Ok(hex::encode(hasher.finalize()))
}

/// SHA256 checksum of an mmap or file-like object (no decompression).
#[pyfunction]
fn calculate_checksum_stream(py: Python, source: &Bound<'_, PyAny>) -> PyResult<String> {
    let input = extract_stream_input(source)?;

    let checksum = py
        .detach(|| match input {
            StreamInput::PyBuf(buffer) => {
                // The buffer is dropped after the cursor, so the slice stays valid.
                let slice = unsafe { buffer_as_slice(&buffer) };
                calculate_stream_checksum(Cursor::new(slice))
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
    m.add_function(wrap_pyfunction!(validate_raw, m)?)?;
    m.add_function(wrap_pyfunction!(validate_raw_stream, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_checksum, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_checksum_stream, m)?)?;
    Ok(())
}
