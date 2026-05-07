import enum
import logging
from importlib.metadata import PackageNotFoundError, version

from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PublicKey
from grz_db.models.base import VerifiableLog

log = logging.getLogger(__name__)


def _get_versions() -> dict[str, str]:
    """Return versions of grzctl and key dependencies."""

    def get_version(package_name: str) -> str:
        try:
            return version(package_name)
        except PackageNotFoundError:
            return "package-not-found"

    return {
        "grzctl": get_version("grzctl"),
        "grz-cli": get_version("grz-cli"),
        "grz-common": get_version("grz-common"),
        "grz-db": get_version("grz-db"),
        "grz-pydantic-models": get_version("grz-pydantic-models"),
    }


class SignatureStatus(enum.StrEnum):
    """Enum for signature status."""

    VERIFIED = "Verified"
    FAILED = "Failed"
    ERROR = "Error"
    UNKNOWN = "Unknown"

    def rich_display(self, comment: str | None) -> str:
        """Displays the signature status in rich format."""
        match self:
            case "Verified":
                return "[green]Verified[/green]" if comment is None else f"[green]Verified ({comment})[/green]"
            case "Failed":
                return "[red]Failed[/red]"
            case "Error":
                return "[red]Error[/red]"
            case "Unknown" | _:
                return "[yellow]Unknown Key[/yellow]"


def _verify_signature(
    public_keys: dict[str, Ed25519PublicKey], expected_key_comment: str, verifiable_log: VerifiableLog
) -> tuple[SignatureStatus, str | None]:
    signature_status = SignatureStatus.UNKNOWN
    verifying_key_comment = None
    if public_key := public_keys.get(expected_key_comment):
        try:
            signature_status = SignatureStatus.VERIFIED if verifiable_log.verify(public_key) else SignatureStatus.FAILED
        except Exception as e:
            signature_status = SignatureStatus.ERROR
            log.error(e)
    else:
        log.debug("Found no key with matching username in comment, trying all keys")
        for comment, public_key in public_keys.items():
            try:
                if verifiable_log.verify(public_key):
                    signature_status = SignatureStatus.VERIFIED
                    verifying_key_comment = comment
                    # stop trying after first verification success
                    break
            except Exception as e:
                signature_status = SignatureStatus.ERROR
                log.error(e)

    return signature_status, verifying_key_comment
