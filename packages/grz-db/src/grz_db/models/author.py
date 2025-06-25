import logging

from cryptography.hazmat.primitives.asymmetric.types import PrivateKeyTypes

log = logging.getLogger(__name__)


class Author:
    def __init__(self, name: str, private_key_bytes: bytes, private_key_passphrase: str | None):
        self.name = name
        self.private_key_bytes = private_key_bytes
        self.private_key_passphrase = private_key_passphrase

    def private_key(self) -> PrivateKeyTypes:
        from functools import partial
        from getpass import getpass

        from cryptography.hazmat.primitives.serialization import load_ssh_private_key

        passphrase = self.private_key_passphrase
        passphrase_callback = (lambda: passphrase) if passphrase else None

        if not passphrase:
            passphrase_callback = partial(getpass, prompt=f"Passphrase for GRZ DB author ({self.name}'s) private key: ")
        log.info(f"Loading private key of {self.name}…")
        private_key = load_ssh_private_key(
            self.private_key_bytes,
            password=passphrase_callback().encode("utf-8"),
        )
        return private_key
