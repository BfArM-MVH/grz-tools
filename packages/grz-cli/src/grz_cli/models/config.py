import logging

from grz_common.models.keys import KeyConfigModel
from grz_common.models.s3 import S3ConfigModel

log = logging.getLogger(__name__)


class UploadConfig(S3ConfigModel):
    pass


class EncryptConfig(KeyConfigModel):
    pass
