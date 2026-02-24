import io
import os

import crypt4gh.header
import crypt4gh.lib
from nacl.bindings import crypto_aead_chacha20poly1305_ietf_encrypt
from nacl.public import PrivateKey

from . import Transformer


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

    def _fill_buffer(self) -> bytes:
        if not self._header_parsed:
            self._read_header()

        if self.source is None:
            raise RuntimeError("Stream source not set")

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

        return crypt4gh.lib.decrypt_block(ciphersegment, self._session_keys)

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
        self._sender_privkey = sender_privkey or bytes(PrivateKey.generate())
        self._session_key = os.urandom(32)
        self._header_sent = False
        self._buffer = bytearray()

    def _fill_buffer(self) -> bytes:
        if not self._header_sent:
            self._header_sent = True
            return self._compose_header()

        if self.source is None:
            raise RuntimeError("Stream source not set")

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

        nonce = os.urandom(self.NONCE_LENGTH)
        ciphertext = crypto_aead_chacha20poly1305_ietf_encrypt(segment, None, nonce, self._session_key)
        return nonce + ciphertext

    def _compose_header(self) -> bytes:
        keys = [(0, self._sender_privkey, self._recipient_pubkey)]
        header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
        header_packets = crypt4gh.header.encrypt(header_content, keys)
        return crypt4gh.header.serialize(header_packets)
