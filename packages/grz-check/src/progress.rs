use indicatif::ProgressBar;

pub(crate) struct DualProgressReader<R: std::io::Read> {
    inner: R,
    specific_pb: ProgressBar,
    global_pb: ProgressBar,
}

impl<R: std::io::Read> DualProgressReader<R> {
    pub fn new(inner: R, pb1: ProgressBar, pb2: ProgressBar) -> Self {
        Self {
            inner,
            specific_pb: pb1,
            global_pb: pb2,
        }
    }
}

impl<R: std::io::Read> std::io::Read for DualProgressReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let bytes_read = self.inner.read(buf)?;
        if bytes_read > 0 {
            let n = bytes_read as u64;
            self.specific_pb.inc(n);
            self.global_pb.inc(n);
        }
        Ok(bytes_read)
    }
}
