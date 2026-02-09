"""Crypt4GH encryption/decryption pipeline stages."""

from __future__ import annotations

import io
import os

import crypt4gh.header
import crypt4gh.lib
from nacl.bindings import (
    crypto_aead_chacha20poly1305_ietf_decrypt as decrypt_block,
)
from nacl.bindings import (
    crypto_aead_chacha20poly1305_ietf_encrypt as encrypt_block,
)
from nacl.public import PrivateKey

from .base import PipelineError, StreamTransformer

# Crypt4GH constants
NONCE_LENGTH = 12
MAC_LENGTH = 16
SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE  # 65536 bytes of plaintext
CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE  # 65536 + 12 + 16 = 65564 bytes


class Crypt4GHDecryptor(StreamTransformer):
    """Streaming Crypt4GH decryptor."""

    def __init__(self, private_key: bytes, name: str | None = None):
        """
        Initialize the decryptor.

        :param private_key: The recipient's Crypt4GH private key
        :param name: Stage name for logging
        """
        super().__init__(name or "Crypt4GHDecryptor")
        self._private_key = private_key
        self._session_key: bytes | None = None
        self._header_parsed = False
        self._header_buffer = bytearray()
        self._segment_buffer = bytearray()
        self._header_length = 0
        self._bytes_in = 0
        self._bytes_out = 0

    def _parse_header(self) -> bool:
        """
        Attempt to parse the header from the buffer.

        :returns: True if header was successfully parsed
        """
        if len(self._header_buffer) < 16:
            return False  # need at least magic + version + packet count

        try:
            header_stream = io.BytesIO(bytes(self._header_buffer))
            keys = [(0, self._private_key, None)]

            # validate magic
            magic = header_stream.read(8)
            if magic != b"crypt4gh":
                raise ValueError("Not a valid Crypt4GH file (invalid magic)")

            version = int.from_bytes(header_stream.read(4), byteorder="little")
            if version != 1:
                raise ValueError(f"Unsupported Crypt4GH version: {version}")

            packet_count = int.from_bytes(header_stream.read(4), byteorder="little")

            # read all header packets
            for _ in range(packet_count):
                if header_stream.tell() + 4 > len(self._header_buffer):
                    return False  # need more data

                packet_length = int.from_bytes(header_stream.read(4), byteorder="little")
                if header_stream.tell() + packet_length - 4 > len(self._header_buffer):
                    return False  # need more data

                header_stream.read(packet_length - 4)  # skip packet data for now

            self._header_length = header_stream.tell()

            # reset and parse the complete header
            header_stream = io.BytesIO(bytes(self._header_buffer[: self._header_length]))
            session_keys, _ = crypt4gh.header.deconstruct(infile=header_stream, keys=keys, sender_pubkey=None)

            if len(session_keys) != 1:
                raise ValueError(f"Expected 1 session key, got {len(session_keys)}")

            self._session_key = session_keys[0]
            self._header_parsed = True

            # move remaining data to segment buffer
            self._segment_buffer.extend(self._header_buffer[self._header_length :])
            self._header_buffer.clear()

            self._log.debug(f"Parsed Crypt4GH header ({self._header_length} bytes)")
            return True

        except ValueError:
            raise
        except Exception:
            return False  # need more data

    def _decrypt_segment(self, segment: bytes) -> bytes:
        """
        Decrypt a single Crypt4GH segment.

        :param segment: Encrypted segment (CIPHER_SEGMENT_SIZE bytes)
        :returns: Decrypted plaintext
        """
        if not segment:
            return b""

        if self._session_key is None:
            raise RuntimeError("Session key not initialized")

        nonce = segment[:NONCE_LENGTH]
        ciphertext = segment[NONCE_LENGTH:]

        try:
            return decrypt_block(ciphertext, None, nonce, self._session_key)
        except Exception as e:
            raise PipelineError(f"Decryption failed: {e}", self.name, e) from e

    def process(self, data: bytes) -> bytes:
        """
        Process encrypted input data.

        :param data: Encrypted bytes (may include header on first call)
        :returns: Decrypted bytes (may be empty if buffering)
        """
        if not data:
            return b""

        self._bytes_in += len(data)

        # parse header first
        if not self._header_parsed:
            self._header_buffer.extend(data)
            if not self._parse_header():
                return b""  # need more data
            # header parsed, segment_buffer has remaining data
        else:
            self._segment_buffer.extend(data)

        # decrypt complete segments
        result = bytearray()
        while len(self._segment_buffer) >= CIPHER_SEGMENT_SIZE:
            segment = bytes(self._segment_buffer[:CIPHER_SEGMENT_SIZE])
            self._segment_buffer = self._segment_buffer[CIPHER_SEGMENT_SIZE:]
            decrypted = self._decrypt_segment(segment)
            result.extend(decrypted)
            self._bytes_out += len(decrypted)

        return bytes(result)

    def flush(self) -> bytes:
        """
        Process any remaining buffered data.

        :returns: Final decrypted bytes
        """
        if not self._header_parsed:
            raise PipelineError("No data processed - header never parsed", self.name)

        if not self._segment_buffer:
            return b""

        # decrypt final (possibly partial) segment
        decrypted = self._decrypt_segment(bytes(self._segment_buffer))
        self._bytes_out += len(decrypted)
        self._segment_buffer.clear()
        return decrypted

    def finalize(self) -> None:
        """Record decryption stats in context."""
        self.context.bytes_decrypted = self._bytes_out
        self.context.metadata["crypt4gh_header_length"] = self._header_length
        self._log.debug(f"Decrypted {self._bytes_in} -> {self._bytes_out} bytes")

    @property
    def header_parsed(self) -> bool:
        """Return whether the header has been parsed."""
        return self._header_parsed

    @property
    def session_key(self) -> bytes | None:
        """Return the session key (for re-encryption)."""
        return self._session_key

    @property
    def bytes_in(self) -> int:
        """Return bytes received."""
        return self._bytes_in

    @property
    def bytes_out(self) -> int:
        """Return bytes produced."""
        return self._bytes_out


class Crypt4GHEncryptor(StreamTransformer):
    """Streaming Crypt4GH encryptor."""

    def __init__(
        self,
        recipient_public_key: bytes,
        sender_private_key: bytes | None = None,
        session_key: bytes | None = None,
        name: str | None = None,
    ):
        """
        Initialize the encryptor.

        :param recipient_public_key: Recipient's Crypt4GH public key
        :param sender_private_key: Sender's private key for signing (optional)
        :param session_key: Session key to use (optional, generated if not provided)
        :param name: Stage name for logging
        """
        super().__init__(name or "Crypt4GHEncryptor")
        self._recipient_public_key = recipient_public_key
        self._sender_private_key = sender_private_key or bytes(PrivateKey.generate())
        self._session_key = session_key or os.urandom(32)
        self._header: bytes | None = None
        self._header_emitted = False
        self._segment_buffer = bytearray()
        self._bytes_in = 0
        self._bytes_out = 0

    def _build_header(self) -> bytes:
        """Build the Crypt4GH header."""
        keys = [(0, self._sender_private_key, self._recipient_public_key)]
        header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
        header_packets = crypt4gh.header.encrypt(header_content, keys)
        return crypt4gh.header.serialize(header_packets)

    def _encrypt_segment(self, plaintext: bytes) -> bytes:
        """
        Encrypt a single segment of plaintext.

        :param plaintext: Plaintext (up to SEGMENT_SIZE bytes)
        :returns: Encrypted segment (nonce + ciphertext + MAC)
        """
        nonce = os.urandom(NONCE_LENGTH)
        ciphertext = encrypt_block(plaintext, None, nonce, self._session_key)
        return nonce + ciphertext

    def get_header(self) -> bytes:
        """
        Get the Crypt4GH header.

        Can be called before or after processing data.
        """
        if self._header is None:
            self._header = self._build_header()
        return self._header

    def process(self, data: bytes) -> bytes:
        """
        Encrypt plaintext data.

        :param data: Plaintext bytes to encrypt
        :returns: Encrypted bytes (may include header on first call)
        """
        if not data:
            return b""

        self._bytes_in += len(data)
        self._segment_buffer.extend(data)

        result = bytearray()

        # emit header on first output
        if not self._header_emitted:
            result.extend(self.get_header())
            self._header_emitted = True
            self._bytes_out += len(self._header)  # type: ignore[arg-type]

        # encrypt complete segments
        while len(self._segment_buffer) >= SEGMENT_SIZE:
            segment = bytes(self._segment_buffer[:SEGMENT_SIZE])
            self._segment_buffer = self._segment_buffer[SEGMENT_SIZE:]
            encrypted = self._encrypt_segment(segment)
            result.extend(encrypted)
            self._bytes_out += len(encrypted)

        return bytes(result)

    def flush(self) -> bytes:
        """
        Encrypt any remaining buffered data.

        :returns: Final encrypted bytes
        """
        result = bytearray()

        # emit header if not yet emitted (empty file case)
        if not self._header_emitted:
            result.extend(self.get_header())
            self._header_emitted = True
            self._bytes_out += len(self._header)  # type: ignore[arg-type]

        # encrypt remaining data
        if self._segment_buffer:
            encrypted = self._encrypt_segment(bytes(self._segment_buffer))
            result.extend(encrypted)
            self._bytes_out += len(encrypted)
            self._segment_buffer.clear()

        return bytes(result)

    def finalize(self) -> None:
        """Record encryption stats in context."""
        self.context.metadata["crypt4gh_bytes_encrypted"] = self._bytes_in
        self._log.debug(f"Encrypted {self._bytes_in} -> {self._bytes_out} bytes")

    @property
    def session_key(self) -> bytes:
        """Return the session key."""
        return self._session_key

    @property
    def bytes_in(self) -> int:
        """Return bytes received."""
        return self._bytes_in

    @property
    def bytes_out(self) -> int:
        """Return bytes produced."""
        return self._bytes_out
