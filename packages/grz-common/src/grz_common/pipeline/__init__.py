"""Modular streaming pipeline components."""

from .base import (
    PipelineContext,
    PipelineError,
    StreamObserver,
    StreamSink,
    StreamSource,
    StreamStage,
    StreamTransformer,
)
from .compressors import GzipDecompressor
from .constants import MULTIPART_DEFAULT_PART_SIZE, MULTIPART_MAX_PARTS, MULTIPART_MIN_PART_SIZE
from .crypto import Crypt4GHDecryptor, Crypt4GHEncryptor
from .operations import (
    DecryptOperation,
    DownloadOperation,
    EncryptOperation,
    UploadOperation,
    ValidateOperation,
)
from .processing import FileProcessingResult, GrzctlProcessPipeline
from .s3 import S3Downloader, S3MultipartUploader
from .s3_keys import S3KeyBuilder
from .utils import SignalManager, abort_all_stages, drain_queue, finalize_stages_in_order, safe_join_thread
from .validators import BamValidator, FastqValidator, RawChecksumValidator

__all__ = [  # noqa: RUF022
    "BamValidator",
    "Crypt4GHDecryptor",
    "Crypt4GHEncryptor",
    "DecryptOperation",
    "DownloadOperation",
    "EncryptOperation",
    "FastqValidator",
    "FileProcessingResult",
    "GrzctlProcessPipeline",
    "GzipDecompressor",
    "MULTIPART_DEFAULT_PART_SIZE",
    "MULTIPART_MAX_PARTS",
    "MULTIPART_MIN_PART_SIZE",
    "PipelineContext",
    "PipelineError",
    "RawChecksumValidator",
    "S3Downloader",
    "S3KeyBuilder",
    "S3MultipartUploader",
    "SignalManager",
    "StreamObserver",
    "StreamSink",
    "StreamSource",
    "StreamStage",
    "StreamTransformer",
    "UploadOperation",
    "ValidateOperation",
    "abort_all_stages",
    "drain_queue",
    "finalize_stages_in_order",
    "safe_join_thread",
]
