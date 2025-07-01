import functools
import logging

from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

log = logging.getLogger(__name__)


class Author:
    def __init__(self, name: str, private_key_bytes: bytes, private_key_passphrase: str | None):
        self.name = name
        self.private_key_bytes = private_key_bytes
        self.private_key_passphrase = private_key_passphrase

    # cache to avoid asking for passphrase multiple times if needed
    @functools.cache  # noqa: B019
    def private_key(self) -> Ed25519PrivateKey:
        from functools import partial
        from getpass import getpass

        from cryptography.hazmat.primitives.serialization import load_ssh_private_key

        passphrase = self.private_key_passphrase
        if passphrase:
            passphrase_callback = lambda: passphrase
        else:
            passphrase_callback = partial(getpass, prompt=f"Passphrase for GRZ DB author ({self.name}'s) private key: ")

        log.info(f"Loading private key of {self.name}…")
        try:
            private_key = load_ssh_private_key(
                self.private_key_bytes,
                password=passphrase_callback().encode("utf-8"),
            )
        except ValueError as e:
            if "Corrupt data: broken checksum" in str(e):
                raise ValueError("Could not load private key, likely incorrect passphrase supplied.") from e
            else:
                raise e
        if not isinstance(private_key, Ed25519PrivateKey):
            raise TypeError(f"private_key must be an Ed25519PrivateKey. Got {type(private_key)}")
        return private_key
