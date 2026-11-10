"""Fixtures of the end-to-end tests: SeaweedFS in a container, a BfArM mock, and one GRZ per test."""

import os
from collections.abc import Iterator
from pathlib import Path

import boto3
import pytest
from grz_common.models.version import VERSION_FILE_KEY, VersionFile

from . import data, services
from .harness import Grz

_ENGINE_ENV = dict(os.environ)
"""The environment before the root conftest replaces ``HOME``; rootless podman keeps its storage under the real one."""


@pytest.fixture(scope="session")
def s3_endpoint(tmp_path_factory) -> Iterator[str]:
    """SeaweedFS with the GRZ's buckets, and the version policy that grz-cli reads from the inbox."""
    with services.seaweedfs(tmp_path_factory.mktemp("seaweedfs"), _ENGINE_ENV) as endpoint:
        s3 = boto3.client(
            "s3",
            endpoint_url=endpoint,
            aws_access_key_id=services.ADMIN[0],
            aws_secret_access_key=services.ADMIN[1],
        )
        s3.put_object(Bucket=services.INBOX, Key=VERSION_FILE_KEY, Body=VersionFile.read_bundled_text().encode())
        yield endpoint


@pytest.fixture(scope="session")
def bfarm_certificate(tmp_path_factory) -> tuple[Path, Path]:
    return services.self_signed_certificate(tmp_path_factory.mktemp("bfarm"))


@pytest.fixture(scope="session")
def bfarm_mock(bfarm_certificate) -> Iterator[services.BfarmMock]:
    mock = services.BfarmMock(*bfarm_certificate)
    yield mock
    mock.close()


@pytest.fixture
def grz(tmp_path, s3_endpoint, migrated_db_connection, bfarm_mock, bfarm_certificate, request) -> Grz:
    """A GRZ with a fresh database, on SQLite and on PostgreSQL."""
    bfarm_mock.reset()
    cache = request.config.cache.mkdir("grz-mini-test-data")
    return Grz(
        tmp_path,
        s3_endpoint,
        migrated_db_connection,
        bfarm_mock,
        bfarm_certificate[0],
        lambda name: data.fetch(name, cache),
    )
