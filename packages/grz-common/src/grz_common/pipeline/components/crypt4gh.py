import io
import os

import crypt4gh.header
import crypt4gh.lib
from crypt4gh.sodium import chacha20poly1305_decrypt, chacha20poly1305_encrypt
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey
from grz_common.exceptions import DecryptionError

from . import StreamConfigurationError, Transformer


class Crypt4GHDecryptor(Transformer):
    """
    Crypt4GH Decryption Transformer.

    Usage:
        - Standalone: Crypt4GHDecryptor(source_stream, private_key=key)
        - Pipeline:   source >> Crypt4GHDecryptor(private_key=key) >> sink
    """

    CIPHER_DIFF = crypt4gh.lib.CIPHER_DIFF
    CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE

    def __init__(self, source: io.BufferedIOBase | None = None, /, *, private_key: bytes):
        """
        :param source: Optional positional-only source stream.
        :param private_key: Required keyword-only private key for decryption.
        """
        super().__init__(source)
        self._private_key = private_key
        self._session_keys: bytes | None = None
        self._header_parsed = False
        self._buffer = bytearray()
        self._out_buffer = bytearray(self.CIPHER_SEGMENT_SIZE)

    def _fill_buffer(self) -> bytes:
        if not self._header_parsed:
            self._read_header()

        if self.source is None:
            raise StreamConfigurationError("Stream source not set", stage=self.__class__.__name__)

        while len(self._buffer) < self.CIPHER_SEGMENT_SIZE:
            chunk = self.source.read(self.CIPHER_SEGMENT_SIZE)
            if not chunk:
                break
            self._buffer.extend(chunk)

        if not self._buffer:
            return b""

        ciphersegment = bytes(self._buffer[: self.CIPHER_SEGMENT_SIZE])
        del self._buffer[: self.CIPHER_SEGMENT_SIZE]

        if len(ciphersegment) <= self.CIPHER_DIFF:
            raise ValueError("Truncated cipher segment")

        if not self._session_keys:
            raise ValueError("No session keys found in Crypt4GH header")

        segment_len = len(ciphersegment) - self.CIPHER_DIFF
        errors = []
        out_view = memoryview(self._out_buffer)[:segment_len]

        for key in self._session_keys:
            try:
                out_len = chacha20poly1305_decrypt(out_view, ciphersegment, bytes(key))
                return bytes(out_view[:out_len])
            except Exception as e:
                errors.append(repr(e))

        raise DecryptionError(f"Decryption failed: {errors}")

    def _read_header(self) -> None:
        try:
            keys = [(0, self._private_key, None)]
            # Decrypt header to extract session keys using crypt4gh library
            session_keys, edit_list = crypt4gh.header.deconstruct(self.source, keys, sender_pubkey=None)

            if edit_list is not None:
                raise ValueError("Edit lists in Crypt4GH headers are not supported!")

            if not session_keys:
                raise ValueError("No session keys found in Crypt4GH header")

            self._session_keys = session_keys
            self._header_parsed = True
        except Exception as e:
            raise OSError(f"Crypt4GH Header Error: {e}") from e


class Crypt4GHEncryptor(Transformer):
    """
    Crypt4GH Encryption Transformer.
    """

    NONCE_LENGTH = 12
    SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE
    CIPHER_DIFF = crypt4gh.lib.CIPHER_DIFF

    def __init__(
        self,
        source: io.BufferedIOBase | None = None,
        /,
        *,
        recipient_pubkey: bytes,
        sender_privkey: bytes | None = None,
    ):
        """
        :param source: Optional positional-only source stream.
        :param recipient_pubkey: Required keyword-only public key of the recipient.
        :param sender_privkey: Optional keyword-only private key of the sender.
        """
        super().__init__(source)
        self._recipient_pubkey = recipient_pubkey
        self._sender_privkey = sender_privkey or (
            X25519PrivateKey.generate().private_bytes(
                encoding=serialization.Encoding.Raw,
                format=serialization.PrivateFormat.Raw,
                encryption_algorithm=serialization.NoEncryption(),
            )
        )
        self._session_key = os.urandom(32)
        self._header_sent = False
        self._buffer = bytearray()
        self._out_buffer = bytearray(self.SEGMENT_SIZE + self.CIPHER_DIFF)

    def _fill_buffer(self) -> bytes:
        if not self._header_sent:
            self._header_sent = True
            return self._compose_header()

        if self.source is None:
            raise StreamConfigurationError("Stream source not set", stage=self.__class__.__name__)

        # ensure we always get SEGMENT_SIZE long segments to minimize crypt4gh overhead
        while len(self._buffer) < self.SEGMENT_SIZE:
            chunk = self.source.read(self.SEGMENT_SIZE)
            if not chunk:
                break
            self._buffer.extend(chunk)

        if not self._buffer:
            return b""

        segment_len = min(self.SEGMENT_SIZE, len(self._buffer))
        segment = bytes(self._buffer[:segment_len])
        del self._buffer[:segment_len]

        out_length = segment_len + self.CIPHER_DIFF
        out_view = memoryview(self._out_buffer)[:out_length]
        out_len = chacha20poly1305_encrypt(out_view, segment, self._session_key)

        return bytes(out_view[:out_len])

    def _compose_header(self) -> bytes:
        keys = [(0, self._sender_privkey, self._recipient_pubkey)]
        header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
        header_packets = crypt4gh.header.encrypt(header_content, keys)
        return crypt4gh.header.serialize(header_packets)
