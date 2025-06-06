"""Hash calculation utilities."""

import hashlib
import logging
from os import PathLike
from os.path import getsize
from pathlib import Path

from tqdm.auto import tqdm

from ..constants import TQDM_SMOOTHING

log = logging.getLogger(__name__)


def calculate_sha256(file_path: str | PathLike, chunk_size=2**16, progress=True) -> str:
    """
    Calculate the sha256 value of a file in chunks

    :param file_path: path to the file
    :param chunk_size: Chunk size in bytes
    :param progress: Print progress
    :return: calculated sha256 value of file_path
    """
    file_path = Path(file_path)
    total_size = getsize(file_path)
    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        if progress and (total_size > chunk_size):
            with tqdm(
                total=total_size,
                unit="B",
                unit_scale=True,
                desc=f"Calculating SHA256 {file_path.name}",
                smoothing=TQDM_SMOOTHING,
            ) as pbar:
                while chunk := f.read(chunk_size):
                    sha256_hash.update(chunk)
                    pbar.update(len(chunk))
        else:
            while chunk := f.read(chunk_size):
                sha256_hash.update(chunk)
    return sha256_hash.hexdigest()
