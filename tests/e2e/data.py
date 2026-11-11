"""The submissions of grz-mini-test-data, downloaded once and cached."""

import hashlib
import shutil
import tarfile
import tempfile
import urllib.request
from pathlib import Path

VERSION = "0.2.3"
_URL = "https://github.com/twrightsman/grz-mini-test-data/releases/download/version%2F{version}/{name}.tgz"

DATASETS = {
    "panel": "fd9240569455939fd19334cbc71fedfdcb53fd33f83d1c9ac439707681a87eab",
    "wgs": "addca402c066f8778f6463c9feae4034c418ad8f2c780c3ef680b9284438a44d",
    "wgs_lr": "7b3e121a602876a3e38f127695314b39fd2c991180b959e0aedd6544e2894028",
}
"""SHA-256 of each release tarball, which holds one submission (``metadata/`` and ``files/``)."""


def fetch(name: str, cache_dir: Path) -> Path:
    """Return the unpacked submission *name*, downloading it into *cache_dir* on first use.

    Parallel test workers may race here, so each writes into a temporary path and moves it into place.

    :param name: Key of :data:`DATASETS`.
    :param cache_dir: Directory that keeps the unpacked submissions between test runs.
    :returns: The submission directory, with ``metadata/metadata.json`` and ``files/``.
    :raises ValueError: If the download does not match its pinned SHA-256.
    """
    target = cache_dir / VERSION / name
    if target.is_dir():
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=target.parent) as tmp:
        tarball = Path(tmp) / f"{name}.tgz"
        with urllib.request.urlopen(_URL.format(version=VERSION, name=name), timeout=60) as response:
            with tarball.open("wb") as f:
                shutil.copyfileobj(response, f)
        digest = hashlib.sha256(tarball.read_bytes()).hexdigest()
        if digest != DATASETS[name]:
            raise ValueError(
                f"{tarball.name} of grz-mini-test-data {VERSION} has SHA-256 {digest}, expected {DATASETS[name]}"
            )
        unpacked = Path(tmp) / name
        with tarfile.open(tarball) as tar:
            tar.extractall(unpacked, filter="data")
        try:
            unpacked.rename(target)
        except OSError:
            if not target.is_dir():
                raise
    return target
