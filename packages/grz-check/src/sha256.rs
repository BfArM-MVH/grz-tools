use sha2::{Digest, Sha256};
use std::io::{self, Read};
use std::sync::{Arc, Mutex};

pub struct SharedHashingReader<R: Read> {
    inner: R,
    hasher: Arc<Mutex<Sha256>>,
}

impl<R: Read> SharedHashingReader<R> {
    pub fn new(inner: R, hasher: Arc<Mutex<Sha256>>) -> Self {
        Self { inner, hasher }
    }
}

impl<R: Read> Read for SharedHashingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let bytes_read = self.inner.read(buf)?;
        if bytes_read > 0 {
            self.hasher.lock().unwrap().update(&buf[..bytes_read]);
        }
        Ok(bytes_read)
    }
}
