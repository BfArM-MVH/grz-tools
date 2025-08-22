from typing import Annotated

from grz_common.models.base import IgnoringBaseSettings
from grz_common.models.keys import KeyConfigModel
from grz_common.models.s3 import S3ConfigModel
from pydantic import Field

from .db import DbModel
from .pruefbericht import PruefberichtModel


class ArchiveConfig(S3ConfigModel):
    pass


class DownloadConfig(S3ConfigModel):
    pass


class DecryptConfig(KeyConfigModel):
    pass


class CleanConfig(S3ConfigModel):
    pass


class PruefberichtConfig(IgnoringBaseSettings):
    pruefbericht: PruefberichtModel


class DbConfig(IgnoringBaseSettings):
    db: DbModel


class ListConfig(S3ConfigModel):
    # invalid DbModel will be silently ignored as dict
    db: Annotated[DbModel | dict | None, Field(union_mode="left_to_right")] = None
