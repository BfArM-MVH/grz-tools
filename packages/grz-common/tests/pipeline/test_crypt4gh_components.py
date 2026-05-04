"""Tests for the Crypt4GH pipeline components."""

import os
from io import BytesIO

from crypt4gh import SEGMENT_SIZE
from grz_common.pipeline.components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor
from nacl.public import PrivateKey


def generate_keypair() -> tuple[bytes, bytes]:
    """Generate a Crypt4GH keypair (private, public)."""
    private = PrivateKey.generate()
    return bytes(private), bytes(private.public_key)


class TestCrypt4GHEncryptor:
    """Tests for the Crypt4GHEncryptor stage."""

    def test_encrypt_small_data(self):
        """Test encryption of data smaller than one segment."""
        private_key, public_key = generate_keypair()

        plaintext = b"Hello, World!"

        with (
            BytesIO(plaintext) as f,
            Crypt4GHEncryptor(f, sender_privkey=private_key, recipient_pubkey=public_key) as encryptor,
        ):
            encrypted = encryptor.read(-1)

        # Encrypted data should include header + encrypted content
        assert len(encrypted) > len(plaintext)
        # Check for crypt4gh magic
        assert encrypted[:8] == b"crypt4gh"

    def test_encrypt_large_data(self):
        """Test encryption of data spanning multiple segments."""
        private_key, public_key = generate_keypair()

        plaintext = os.urandom(SEGMENT_SIZE * 2 + 1000)

        with (
            BytesIO(plaintext) as f,
            Crypt4GHEncryptor(f, sender_privkey=private_key, recipient_pubkey=public_key) as encryptor,
        ):
            encrypted = encryptor.read(-1)

        assert encrypted[:8] == b"crypt4gh"
        assert len(encrypted) > len(plaintext)


class TestCrypt4GHDecryptor:
    """Tests for the Crypt4GHDecryptor stage."""

    def test_decrypt_small_data(self):
        """Test decryption of data smaller than one segment."""
        private_key, public_key = generate_keypair()

        plaintext = b"Hello, World! This is a test message."

        # First encrypt
        with (
            BytesIO(plaintext) as f,
            Crypt4GHEncryptor(f, sender_privkey=private_key, recipient_pubkey=public_key) as encryptor,
        ):
            encrypted = encryptor.read(-1)

        # Then decrypt
        with (
            BytesIO(encrypted) as f,
            Crypt4GHDecryptor(f, private_key=private_key) as decryptor,
        ):
            decrypted = decryptor.read(-1)

        assert decrypted == plaintext

    def test_decrypt_large_data(self):
        """Test decryption of data spanning multiple segments."""
        private_key, public_key = generate_keypair()

        # Create data larger than one segment
        plaintext = os.urandom(SEGMENT_SIZE * 3 + 500)

        # First encrypt
        with (
            BytesIO(plaintext) as f,
            Crypt4GHEncryptor(f, sender_privkey=private_key, recipient_pubkey=public_key) as encryptor,
        ):
            encrypted = encryptor.read(-1)

        # Then decrypt
        with (
            BytesIO(encrypted) as f,
            Crypt4GHDecryptor(f, private_key=private_key) as decryptor,
        ):
            decrypted = decryptor.read(-1)

        assert decrypted == plaintext
