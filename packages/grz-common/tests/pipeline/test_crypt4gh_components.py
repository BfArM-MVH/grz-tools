"""Tests for the Crypt4GH pipeline components."""

import io
import os

from grz_common.pipeline.base import PipelineContext
from grz_common.pipeline.crypto import (
    SEGMENT_SIZE,
    Crypt4GHDecryptor,
    Crypt4GHEncryptor,
)
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

        encryptor = Crypt4GHEncryptor(
            recipient_public_key=public_key,
            sender_private_key=private_key,
        )
        context = PipelineContext()
        encryptor.initialize(context)

        plaintext = b"Hello, World!"
        encrypted = encryptor.process(plaintext)
        encrypted += encryptor.flush()
        encryptor.finalize()

        # Encrypted data should include header + encrypted content
        assert len(encrypted) > len(plaintext)
        # Check for crypt4gh magic
        assert encrypted[:8] == b"crypt4gh"

    def test_encrypt_large_data(self):
        """Test encryption of data spanning multiple segments."""
        private_key, public_key = generate_keypair()

        encryptor = Crypt4GHEncryptor(
            recipient_public_key=public_key,
            sender_private_key=private_key,
        )
        context = PipelineContext()
        encryptor.initialize(context)

        # Create data larger than one segment
        plaintext = os.urandom(SEGMENT_SIZE * 2 + 1000)

        # Feed in chunks
        chunk_size = 8192
        result = io.BytesIO()
        for i in range(0, len(plaintext), chunk_size):
            chunk = plaintext[i : i + chunk_size]
            encrypted = encryptor.process(chunk)
            result.write(encrypted)

        result.write(encryptor.flush())
        encryptor.finalize()

        encrypted = result.getvalue()
        assert encrypted[:8] == b"crypt4gh"
        assert len(encrypted) > len(plaintext)


class TestCrypt4GHDecryptor:
    """Tests for the Crypt4GHDecryptor stage."""

    def test_decrypt_small_data(self):
        """Test decryption of data smaller than one segment."""
        private_key, public_key = generate_keypair()

        # First encrypt
        encryptor = Crypt4GHEncryptor(
            recipient_public_key=public_key,
            sender_private_key=private_key,
        )
        context = PipelineContext()
        encryptor.initialize(context)

        plaintext = b"Hello, World! This is a test message."
        encrypted = encryptor.process(plaintext)
        encrypted += encryptor.flush()
        encryptor.finalize()

        # Then decrypt
        decryptor = Crypt4GHDecryptor(private_key=private_key)
        context2 = PipelineContext()
        decryptor.initialize(context2)

        decrypted = decryptor.process(encrypted)
        decrypted += decryptor.flush()
        decryptor.finalize()

        assert decrypted == plaintext

    def test_decrypt_large_data(self):
        """Test decryption of data spanning multiple segments."""
        private_key, public_key = generate_keypair()

        # Create data larger than one segment
        plaintext = os.urandom(SEGMENT_SIZE * 3 + 500)

        # Encrypt
        encryptor = Crypt4GHEncryptor(
            recipient_public_key=public_key,
            sender_private_key=private_key,
        )
        context = PipelineContext()
        encryptor.initialize(context)

        encrypted_buffer = io.BytesIO()
        chunk_size = 8192
        for i in range(0, len(plaintext), chunk_size):
            encrypted = encryptor.process(plaintext[i : i + chunk_size])
            encrypted_buffer.write(encrypted)
        encrypted_buffer.write(encryptor.flush())
        encryptor.finalize()

        encrypted = encrypted_buffer.getvalue()

        # Decrypt in chunks
        decryptor = Crypt4GHDecryptor(private_key=private_key)
        context2 = PipelineContext()
        decryptor.initialize(context2)

        decrypted_buffer = io.BytesIO()
        for i in range(0, len(encrypted), chunk_size):
            decrypted = decryptor.process(encrypted[i : i + chunk_size])
            decrypted_buffer.write(decrypted)
        decrypted_buffer.write(decryptor.flush())
        decryptor.finalize()

        assert decrypted_buffer.getvalue() == plaintext

    def test_roundtrip_with_different_chunk_sizes(self):
        """Test that encryption/decryption works with various chunk sizes."""
        private_key, public_key = generate_keypair()

        plaintext = os.urandom(SEGMENT_SIZE * 2 + 12345)

        for encrypt_chunk_size in [1024, 8192, 65536, 100000]:
            for decrypt_chunk_size in [1024, 8192, 65536, 100000]:
                # Encrypt
                encryptor = Crypt4GHEncryptor(
                    recipient_public_key=public_key,
                    sender_private_key=private_key,
                )
                context = PipelineContext()
                encryptor.initialize(context)

                encrypted_buffer = io.BytesIO()
                for i in range(0, len(plaintext), encrypt_chunk_size):
                    encrypted = encryptor.process(plaintext[i : i + encrypt_chunk_size])
                    encrypted_buffer.write(encrypted)
                encrypted_buffer.write(encryptor.flush())
                encrypted = encrypted_buffer.getvalue()

                # Decrypt
                decryptor = Crypt4GHDecryptor(private_key=private_key)
                context2 = PipelineContext()
                decryptor.initialize(context2)

                decrypted_buffer = io.BytesIO()
                for i in range(0, len(encrypted), decrypt_chunk_size):
                    decrypted = decryptor.process(encrypted[i : i + decrypt_chunk_size])
                    decrypted_buffer.write(decrypted)
                decrypted_buffer.write(decryptor.flush())

                assert decrypted_buffer.getvalue() == plaintext, (
                    f"Failed with encrypt_chunk={encrypt_chunk_size}, decrypt_chunk={decrypt_chunk_size}"
                )
