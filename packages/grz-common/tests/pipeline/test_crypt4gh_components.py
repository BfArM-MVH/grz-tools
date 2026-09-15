"""Tests for the Crypt4GH pipeline components."""

import os
from io import BytesIO

import crypt4gh.lib
import pytest
from crypt4gh import SEGMENT_SIZE
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from grz_common.exceptions import DecryptionError
from grz_common.pipeline.components.crypt4gh import Crypt4GHDecryptor, Crypt4GHEncryptor


def generate_keypair() -> tuple[bytes, bytes]:
    """Generate a Crypt4GH keypair (private, public) as raw X25519 keys."""
    private = X25519PrivateKey.generate()
    private_bytes = private.private_bytes(
        encoding=serialization.Encoding.Raw,
        format=serialization.PrivateFormat.Raw,
        encryption_algorithm=serialization.NoEncryption(),
    )
    public_bytes = private.public_key().public_bytes(
        encoding=serialization.Encoding.Raw,
        format=serialization.PublicFormat.Raw,
    )
    return private_bytes, public_bytes


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

    def test_decrypt_rejects_a_tampered_segment(self):
        """A changed byte fails the segment's authentication instead of producing wrong plaintext."""
        sender_private, _ = generate_keypair()
        recipient_private, recipient_public = generate_keypair()
        encrypted = BytesIO()
        crypt4gh.lib.encrypt([(0, sender_private, recipient_public)], BytesIO(os.urandom(1000)), encrypted)
        tampered = bytearray(encrypted.getvalue())
        tampered[-1] ^= 1  # last byte of the segment's MAC

        with (
            BytesIO(bytes(tampered)) as f,
            Crypt4GHDecryptor(f, private_key=recipient_private) as decryptor,
            pytest.raises(DecryptionError),
        ):
            decryptor.read(-1)


class TestCrypt4GHInterop:
    """The components encrypt segments themselves, so they must stay compatible with the reference crypt4gh."""

    @pytest.mark.parametrize("size", [0, 1, SEGMENT_SIZE, 2 * SEGMENT_SIZE + 12345])
    def test_reference_decrypts_what_the_encryptor_writes(self, size):
        sender_private, _ = generate_keypair()
        recipient_private, recipient_public = generate_keypair()
        plaintext = os.urandom(size)

        with (
            BytesIO(plaintext) as f,
            Crypt4GHEncryptor(f, sender_privkey=sender_private, recipient_pubkey=recipient_public) as encryptor,
        ):
            encrypted = encryptor.read(-1)

        decrypted = BytesIO()
        crypt4gh.lib.decrypt([(0, recipient_private, None)], BytesIO(encrypted), decrypted)
        assert decrypted.getvalue() == plaintext

    @pytest.mark.parametrize("size", [0, 1, SEGMENT_SIZE, 2 * SEGMENT_SIZE + 12345])
    def test_decryptor_reads_what_the_reference_writes(self, size):
        sender_private, _ = generate_keypair()
        recipient_private, recipient_public = generate_keypair()
        plaintext = os.urandom(size)
        encrypted = BytesIO()
        crypt4gh.lib.encrypt([(0, sender_private, recipient_public)], BytesIO(plaintext), encrypted)

        with (
            BytesIO(encrypted.getvalue()) as f,
            Crypt4GHDecryptor(f, private_key=recipient_private) as decryptor,
        ):
            decrypted = decryptor.read(-1)

        assert decrypted == plaintext


@pytest.mark.parametrize("size", [0, 1, SEGMENT_SIZE - 1, SEGMENT_SIZE, SEGMENT_SIZE + 1, 3 * SEGMENT_SIZE])
def test_reencryption_does_not_grow_the_object(size):
    """Re-encrypting must not make the object larger: grzctl process sizes its upload parts by the inbox object."""
    submitter_private, _ = generate_keypair()
    grz_private, grz_public = generate_keypair()
    _, archive_public = generate_keypair()

    # the inbox object, encrypted by the reference crypt4gh implementation
    inbox = BytesIO()
    crypt4gh.lib.encrypt([(0, submitter_private, grz_public)], BytesIO(os.urandom(size)), inbox)

    with (
        BytesIO(inbox.getvalue()) as f,
        Crypt4GHDecryptor(f, private_key=grz_private) as decryptor,
        Crypt4GHEncryptor(decryptor, recipient_pubkey=archive_public) as encryptor,
    ):
        reencrypted = encryptor.read(-1)

    assert len(reencrypted) <= len(inbox.getvalue())
