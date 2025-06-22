"""I/O utilities."""

import io
import json
import logging
from typing import TextIO

log = logging.getLogger(__name__)


def read_multiple_json(input_file: TextIO):
    """
    Read multiple JSON objects from a text stream.
    :param input_file:
    """
    for line in input_file:
        if line.strip():
            yield json.loads(line)
        else:
            raise ValueError("Encountered blank line while reading jsonl, invalid.")


class TqdmIOWrapper(io.RawIOBase):
    """
    Wrapper to record reads and writes in a tqdm progress bar.

    Example:
        very_long_input = io.StringIO("0123456789abcdef" * 100000)
        total_size = len(very_long_input.getvalue())

        batch_size = 10 ** 4  # 10kb
        with open("/dev/null", "w") as fd:
            with TqdmIOWrapper(fd, tqdm(
                    total=total_size,
                    desc="Printing: ",
                    unit="B",
                    unit_scale=True,
                    # unit_divisor=1024,  # make use of standard units e.g. KB, MB, etc.
                    miniters=1,
            )) as pbar_fd:
                while chunk := very_long_input.read(batch_size):
                    pbar_fd.write(chunk)

                    time.sleep(0.1)
    """

    def __init__(self, io_buf: io.RawIOBase, progress_bar):
        """

        :param io_buf: the buffer to wrap
        :param progress_bar: tqdm progress bar
        """
        self.io_buf = io_buf
        self.callback = progress_bar.update

    def write(self, data):
        """Write data to the buffer and update the progress bar"""
        nbytes_written = self.io_buf.write(data)
        if nbytes_written:
            self.callback(nbytes_written)
        return nbytes_written

    def read(self, size=-1) -> bytes | None:
        """Read data from the buffer and update the progress bar"""
        data = self.io_buf.read(size)
        if data:
            self.callback(len(data))

        return data

    def readinto(self, buffer, /):
        """Read data into a buffer and update the progress bar"""
        nbytes_written = self.io_buf.readinto(buffer)
        if nbytes_written:
            self.callback(nbytes_written)

        return nbytes_written

    def flush(self):
        """Flush the buffer"""
        # Ensure all data is flushed to the underlying binary IO object
        self.io_buf.flush()

    def close(self):
        """Close the buffer"""
        # Close the underlying binary IO object
        self.io_buf.close()
