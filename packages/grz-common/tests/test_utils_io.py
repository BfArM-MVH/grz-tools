import hashlib

import pytest
from grz_common.utils.io import HashingIOWrapper, TqdmIOWrapper
from tqdm.auto import tqdm


def wrap_with_tqdm(io_wrapper, total_size):
    """Helper to wrap an IO wrapper with tqdm progress bar"""
    return TqdmIOWrapper(io_wrapper, tqdm(total=total_size, disable=True))


class TestHashingIOWrapper:
    ALGORITHMS = ["md5", "sha1", "sha224", "sha256", "sha384", "sha512"]

    @pytest.mark.parametrize("algorithm", ALGORITHMS)
    @pytest.mark.parametrize("with_tqdm", [False, True], ids=["plain", "with_tqdm"])
    def test_read_with_algorithm(self, tmp_path, algorithm, with_tqdm):
        """Test that reading a file calculates the correct hash for each algorithm"""
        test_data = b"Hello, World! This is a test file for hashing."
        test_file = tmp_path / "test_file.bin"
        test_file.write_bytes(test_data)

        expected_hash = hashlib.new(algorithm, test_data).hexdigest()

        with open(test_file, "rb") as fd:
            hashing_fd = HashingIOWrapper(fd, algorithm)
            wrapped_fd = wrap_with_tqdm(hashing_fd, len(test_data)) if with_tqdm else hashing_fd
            with wrapped_fd:
                data = wrapped_fd.read()
                assert data == test_data
            calculated_hash = hashing_fd.hexdigest()

        assert calculated_hash == expected_hash

    @pytest.mark.parametrize("algorithm", ALGORITHMS)
    @pytest.mark.parametrize("with_tqdm", [False, True], ids=["plain", "with_tqdm"])
    def test_write_with_algorithm(self, tmp_path, algorithm, with_tqdm):
        """Test that writing data calculates the correct hash for each algorithm"""
        test_data = b"Hello, World! This is written data."
        test_file = tmp_path / "output_file.bin"

        expected_hash = hashlib.new(algorithm, test_data).hexdigest()

        with open(test_file, "wb") as fd:
            hashing_fd = HashingIOWrapper(fd, algorithm)
            wrapped_fd = wrap_with_tqdm(hashing_fd, len(test_data)) if with_tqdm else hashing_fd
            with wrapped_fd:
                wrapped_fd.write(test_data)
                calculated_hash = hashing_fd.hexdigest()

        assert calculated_hash == expected_hash
        assert test_file.read_bytes() == test_data

    @pytest.mark.parametrize(
        "chunk_size",
        [100, 1000, 10000],
        ids=["small_chunks", "medium_chunks", "large_chunks"],
    )
    @pytest.mark.parametrize("with_tqdm", [False, True], ids=["plain", "with_tqdm"])
    def test_read_chunked(self, tmp_path, chunk_size, with_tqdm):
        """Test that reading in chunks calculates the correct hash"""
        test_data = b"A" * 10000 + b"B" * 10000 + b"C" * 10000
        test_file = tmp_path / "test_file.bin"
        test_file.write_bytes(test_data)

        expected_hash = hashlib.sha256(test_data).hexdigest()

        with open(test_file, "rb") as fd:
            hashing_fd = HashingIOWrapper(fd, "sha256")
            wrapped_fd = wrap_with_tqdm(hashing_fd, len(test_data)) if with_tqdm else hashing_fd
            with wrapped_fd:
                chunks = []
                while chunk := wrapped_fd.read(chunk_size):
                    chunks.append(chunk)

                assert b"".join(chunks) == test_data
                calculated_hash = hashing_fd.hexdigest()

        assert calculated_hash == expected_hash

    @pytest.mark.parametrize(
        "chunks",
        [
            [b"Single chunk"],
            [b"First chunk. ", b"Second chunk. ", b"Third chunk."],
            [b"A" * 1000, b"B" * 1000, b"C" * 1000, b"D" * 1000],
        ],
        ids=["single_chunk", "three_chunks", "four_large_chunks"],
    )
    @pytest.mark.parametrize("with_tqdm", [False, True], ids=["plain", "with_tqdm"])
    def test_write_chunked(self, tmp_path, chunks, with_tqdm):
        """Test that writing in chunks calculates the correct hash"""
        test_data = b"".join(chunks)
        test_file = tmp_path / "output_file.bin"

        expected_hash = hashlib.sha256(test_data).hexdigest()

        with open(test_file, "wb") as fd:
            hashing_fd = HashingIOWrapper(fd, "sha256")
            wrapped_fd = wrap_with_tqdm(hashing_fd, len(test_data)) if with_tqdm else hashing_fd
            with wrapped_fd:
                for chunk in chunks:
                    wrapped_fd.write(chunk)
                calculated_hash = hashing_fd.hexdigest()

        assert calculated_hash == expected_hash
        assert test_file.read_bytes() == test_data

    @pytest.mark.parametrize(
        "test_data",
        [b"", b"x", b"A" * 1000000],
        ids=["empty", "single_byte", "large_file"],
    )
    @pytest.mark.parametrize("with_tqdm", [False, True], ids=["plain", "with_tqdm"])
    def test_various_file_sizes(self, tmp_path, test_data, with_tqdm):
        """Test hashing files of various sizes"""
        test_file = tmp_path / "test_file.bin"
        test_file.write_bytes(test_data)

        expected_hash = hashlib.sha256(test_data).hexdigest()

        with open(test_file, "rb") as fd:
            hashing_fd = HashingIOWrapper(fd, "sha256")
            wrapped_fd = wrap_with_tqdm(hashing_fd, len(test_data)) if with_tqdm else hashing_fd
            with wrapped_fd:
                data = wrapped_fd.read()
                assert data == test_data
                calculated_hash = hashing_fd.hexdigest()

        assert calculated_hash == expected_hash

    @pytest.mark.parametrize("with_tqdm", [False, True], ids=["plain", "with_tqdm"])
    def test_readinto(self, tmp_path, with_tqdm):
        """Test the readinto method"""
        test_data = b"Test data for readinto"
        test_file = tmp_path / "test_file.bin"
        test_file.write_bytes(test_data)

        expected_hash = hashlib.sha256(test_data).hexdigest()

        with open(test_file, "rb") as fd:
            hashing_fd = HashingIOWrapper(fd, "sha256")
            wrapped_fd = wrap_with_tqdm(hashing_fd, len(test_data)) if with_tqdm else hashing_fd
            with wrapped_fd:
                buffer = bytearray(len(test_data))
                nbytes = wrapped_fd.readinto(buffer)
                assert nbytes == len(test_data)
                assert bytes(buffer) == test_data
                calculated_hash = hashing_fd.hexdigest()

        assert calculated_hash == expected_hash
