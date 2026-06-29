import logging
from typing import Annotated

from grz_common.models.identifiers import IdentifiersConfigModel
from grz_common.models.keys import KeyConfigModel
from grz_common.models.s3 import S3ConfigModel, S3Options
from pydantic import Field

log = logging.getLogger(__name__)


class UploadConfig(S3ConfigModel):
    pass


class EncryptConfig(KeyConfigModel):
    pass


class ValidateConfig(IdentifiersConfigModel):
    s3: Annotated[S3Options | None, Field(default=None)] = None
