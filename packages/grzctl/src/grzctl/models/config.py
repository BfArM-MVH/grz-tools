import logging
import sys
from typing import Annotated

from grz_common.models.base import IgnoringBaseModel, IgnoringBaseSettings
from grz_common.models.identifiers import IdentifiersModel
from grz_common.models.keys import KeyModel
from grz_common.models.s3 import S3ConnectionBase, S3Options
from pydantic import Field, model_validator
from pydantic import ValidationError as PydanticValidationError

log = logging.getLogger(__name__)

from .db import DbModel
from .pruefbericht import PruefberichtModel


class InboxConfig(S3ConnectionBase):
    """
    Configuration for a specific inbox.
    Includes connection details and the private key needed to decrypt its contents.
    """

    private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key used to decrypt files from this inbox."""

    private_key_passphrase: Annotated[str | None, Field(default=None)] = None
    """Passphrase to the GRZ private key used to decrypt files from this inbox."""


class InboxTarget(IgnoringBaseModel):
    """
    Fully resolved source configuration.
    Encapsulates everything needed to read and decrypt from a specific inbox.
    """

    s3: S3Options
    """Fully resolved S3 options, including the bucket name."""

    private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key used to decrypt files from this inbox."""

    private_key_passphrase: Annotated[str | None, Field(default=None)] = None
    """Passphrase to the GRZ private key used to decrypt files from this inbox."""


class ProcessS3Options(IgnoringBaseModel):
    """
    Root S3 configuration for the process command.
    Enforces that at least one LE and one inbox are defined.
    """

    inboxes: Annotated[
        dict[str, Annotated[dict[str, InboxConfig], Field(min_length=1)]],
        Field(min_length=1),
    ]
    """
    Mapping: LE-Id -> BucketName -> InboxConfig.
    """


class ProcessKeyConfigModel(IgnoringBaseSettings):
    """Key configuration for the process command."""

    grz_private_key_path: Annotated[str, Field(min_length=1)]
    """Path to the GRZ private key for decryption."""

    consented_archive_public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of consented submissions."""

    non_consented_archive_public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of non-consented submissions."""


class DetailedQcModel(IgnoringBaseSettings):
    local_storage: Annotated[str, Field(min_length=1)]
    """Path to local storage for detailed QC staging."""

    salt: str
    """Salt to use for deterministic determination of submissions selected for detailed QC."""

    target_percentage: Annotated[float, Field(ge=0.0, le=100.0)] = 2.0
    """Target percentage of submissions selected for detailed QC per month."""


class ArchiveTarget(IgnoringBaseModel):
    """Encapsulates everything needed to write to a specific archive."""

    s3: S3Options
    """S3 connection details and bucket for this archive."""

    public_key_path: Annotated[str, Field(min_length=1)]
    """Path to the public key for re-encryption of files destined for this archive."""


class InterrogationConfig(IgnoringBaseModel):
    """Configuration for the staging interrogation bucket."""

    s3: S3Options
    """S3 connection details and bucket for the staging interrogation bucket."""

    keep_failed: bool = False
    """If true, leaves the failed submission files in the interrogation bucket. Otherwise deletes them."""


class ArchivesConfig(IgnoringBaseModel):
    """Configuration for consented and non-consented archives."""

    consented: ArchiveTarget
    """Target definition for consented submissions."""

    non_consented: ArchiveTarget
    """Target definition for non-consented submissions."""

    interrogation: InterrogationConfig
    """Target definition for the intermediate interrogation bucket."""

    @model_validator(mode="after")
    def check_endpoints_match(self) -> "ArchivesConfig":
        consented_endpoint = str(self.consented.s3.endpoint_url) if self.consented.s3.endpoint_url else None
        non_consented_endpoint = str(self.non_consented.s3.endpoint_url) if self.non_consented.s3.endpoint_url else None
        interrogation_endpoint = str(self.interrogation.s3.endpoint_url) if self.interrogation.s3.endpoint_url else None

        if interrogation_endpoint != consented_endpoint:
            log.warning(
                "Interrogation bucket endpoint (%s) differs from consented archive endpoint (%s). "
                "Server-side copying might be slow or fail.",
                interrogation_endpoint,
                consented_endpoint,
            )

        if interrogation_endpoint != non_consented_endpoint:
            log.warning(
                "Interrogation bucket endpoint (%s) differs from non-consented archive endpoint (%s). "
                "Server-side copying might be slow or fail.",
                interrogation_endpoint,
                non_consented_endpoint,
            )

        return self

    @model_validator(mode="after")
    def check_buckets_are_unique(self) -> "ArchivesConfig":
        buckets = {self.consented.s3.bucket, self.non_consented.s3.bucket, self.interrogation.s3.bucket}
        if len(buckets) != 3:
            raise ValueError("consented, non-consented and interrogation buckets must be distinct.")
        return self


def _require_section(config, field_name: str, section_name: str):
    """Resolve an optional config section, exiting with a clear message if missing."""
    value = getattr(config, field_name)
    if value is None:
        log.error(
            "Configuration section '%s' is required for this command. Add a '%s' section to your config.",
            section_name,
            section_name,
        )
        sys.exit(1)
    return value


class GrzctlConfig(IgnoringBaseSettings):
    """Unified configuration for all grzctl commands."""

    s3: ProcessS3Options
    """Configuration for S3 inbox connections."""

    archives: ArchivesConfig
    """Configuration for consented and non-consented archives."""

    db: DbModel | None = None
    """Database configuration for submission tracking."""

    pruefbericht: PruefberichtModel | None = None
    """Configuration for Prüfbericht submission."""

    detailed_qc: DetailedQcModel | None = None
    """Configuration for detailed QC selection and staging."""

    keys: KeyModel | None = None
    """Key configuration for encryption/decryption commands."""

    identifiers: IdentifiersModel | None = None
    """Identifiers for the GRZ and LE."""

    @model_validator(mode="before")
    @classmethod
    def _drop_incomplete_optional_sections(cls, data: dict) -> dict:
        """Drop optional sections that are present but incomplete (e.g. from partial env vars)."""
        optional_sections = {
            "db": DbModel,
            "pruefbericht": PruefberichtModel,
            "detailed_qc": DetailedQcModel,
        }
        for field_name, model_cls in optional_sections.items():
            if field_name in data and data[field_name] is not None:
                try:
                    model_cls.model_validate(data[field_name])  # type: ignore[attr-defined]
                except PydanticValidationError:
                    data[field_name] = None
        return data

    def require_db(self) -> DbModel:
        """Return the db section, or exit if missing."""
        return _require_section(self, "db", "db")

    def require_pruefbericht(self) -> PruefberichtModel:
        """Return the pruefbericht section, or exit if missing."""
        return _require_section(self, "pruefbericht", "pruefbericht")

    def require_detailed_qc(self) -> DetailedQcModel:
        """Return the detailed_qc section, or exit if missing."""
        return _require_section(self, "detailed_qc", "detailed_qc")

    def require_keys(self) -> KeyModel:
        """Return the keys section, or exit if missing."""
        return _require_section(self, "keys", "keys")

    def require_identifiers(self) -> IdentifiersModel:
        """Return the identifiers section, or exit if missing."""
        return _require_section(self, "identifiers", "identifiers")

    def resolve_inbox_by_submission_id(self, submission_id: str, bucket: str | None = None) -> InboxTarget:
        """Resolve an inbox by extracting the LE ID from the submission ID.

        The LE ID is everything before the first ``_`` in the submission ID.
        """
        le_id = submission_id.split("_", maxsplit=1)[0]

        if le_id not in self.s3.inboxes:
            available = ", ".join(self.s3.inboxes.keys())
            log.error("Submitter ID '%s' not found in configuration. Available: %s", le_id, available)
            sys.exit(1)

        submitter_inboxes = self.s3.inboxes[le_id]

        if bucket:
            if bucket not in submitter_inboxes:
                log.error("Inbox bucket '%s' not found for '%s'.", bucket, le_id)
                sys.exit(1)
            bucket_name = bucket
            inbox_config = submitter_inboxes[bucket]
        elif len(submitter_inboxes) == 1:
            bucket_name, inbox_config = next(iter(submitter_inboxes.items()))
        else:
            available_buckets = ", ".join(submitter_inboxes.keys())
            log.error(
                "Multiple inboxes found for '%s' (%s). Please specify --inbox-bucket.",
                le_id,
                available_buckets,
            )
            sys.exit(1)

        s3_options = S3Options(bucket=bucket_name, **inbox_config.model_dump())
        return InboxTarget(
            s3=s3_options, **inbox_config.model_dump(include={"private_key_path", "private_key_passphrase"})
        )

    def resolve_inbox_by_bucket(self, bucket: str | None = None) -> S3Options:
        """Resolve S3Options by iterating all inboxes.

        If *bucket* is specified, find it across all LEs.
        Otherwise require exactly one inbox total.
        """
        if bucket:
            for _le_id, buckets in self.s3.inboxes.items():
                if bucket in buckets:
                    inbox_config = buckets[bucket]
                    return S3Options(bucket=bucket, **inbox_config.model_dump())
            log.error("Inbox bucket '%s' not found in any configured LE.", bucket)
            sys.exit(1)

        all_buckets = []
        for le_id, buckets in self.s3.inboxes.items():
            for bucket_name, inbox_config in buckets.items():
                all_buckets.append((le_id, bucket_name, inbox_config))

        if len(all_buckets) == 1:
            _, bucket_name, inbox_config = all_buckets[0]
            return S3Options(bucket=bucket_name, **inbox_config.model_dump())

        bucket_list = ", ".join(f"{le}/{b}" for le, b, _ in all_buckets)
        log.error("Multiple inboxes configured (%s). Please specify --inbox-bucket.", bucket_list)
        sys.exit(1)
