import io
from typing import Any

from .base import ObserverStream


class TqdmObserver(ObserverStream):
    """
    Updates a shared TQDM progress bar.
    """

    def __init__(self, source: io.BufferedIOBase, pbar: Any):
        super().__init__(source)
        self.pbar = pbar

    def observe(self, chunk: bytes) -> None:
        self.pbar.update(len(chunk))
