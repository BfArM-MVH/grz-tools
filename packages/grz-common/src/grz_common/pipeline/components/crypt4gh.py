import os

import crypt4gh.header
import crypt4gh.lib
from cryptography.exceptions import InvalidTag
from cryptography.hazmat.primitives import serialization
from cryptography.hazmat.primitives.asymmetric.x25519 import X25519PrivateKey

# Segments go through cryptography's ChaCha20Poly1305, not through crypt4gh's libsodium binding:
# it releases the GIL while it computes, and it is faster.
from cryptography.hazmat.primitives.ciphers.aead import ChaCha20Poly1305
from grz_common.exceptions import DecryptionError, EncryptionError

from . import Readable, StreamConfigurationError, Transformer


def _read_segment(source: Readable, size: int) -> bytes:
    """Read one segment of ``size`` bytes, or fewer only at the end of the stream.

    Crypt4GH requires every segment except the last one to be full.
    A raw source may return fewer bytes than asked before its end, so this reads again for the missing bytes.
    Only an empty read ends the stream.

    :param source: Stream to read the segment from.
    :param size: Length of a full segment.
    :returns: The segment, or ``b""`` at the end of the stream.
    """
    segment = source.read(size)
    if not 0 < len(segment) < size:
        return segment
    # a bytearray grows in place, while adding to bytes would copy the segment again on every short read
    buffer = bytearray(segment)
    while len(buffer) < size:
        rest = source.read(size - len(buffer))
        if not rest:
            break
        buffer += rest
    return bytes(buffer)


class Crypt4GHDecryptor(Transformer):
    """
    Crypt4GH Decryption Transformer.

    Usage:
        - Standalone: Crypt4GHDecryptor(source_stream, private_key=key)
        - Pipeline:   source >> Crypt4GHDecryptor(private_key=key) >> sink
    """

    NONCE_LENGTH = 12
    CIPHER_DIFF = crypt4gh.lib.CIPHER_DIFF
    CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE

    def __init__(self, source: Readable | None = None, /, *, private_key: bytes):
        """
        :param source: Optional positional-only source stream.
        :param private_key: Required keyword-only private key for decryption.
        """
        super().__init__(source)
        self._private_key = private_key
        self._ciphers: list[ChaCha20Poly1305] = []  # one per session key in the header
        self._header_parsed = False

    def _fill_buffer(self) -> bytes:
        if not self._header_parsed:
            self._read_header()

        if self.source is None:
            raise StreamConfigurationError("Stream source not set", stage=self.__class__.__name__)

        ciphersegment = _read_segment(self.source, self.CIPHER_SEGMENT_SIZE)
        if not ciphersegment:
            return b""

        if len(ciphersegment) <= self.CIPHER_DIFF:
            raise DecryptionError("Truncated cipher segment")

        nonce = ciphersegment[: self.NONCE_LENGTH]
        ciphertext = memoryview(ciphersegment)[self.NONCE_LENGTH :]
        errors = []

        for cipher in self._ciphers:
            try:
                return cipher.decrypt(nonce, ciphertext, None)
            except InvalidTag as e:
                errors.append(repr(e))

        raise DecryptionError(f"Decryption failed: {errors}")

    def _read_header(self) -> None:
        # crypt4gh raises ValueError for a malformed header or a key that decrypts no header packet;
        # other errors, such as a failed read from the source, propagate unchanged
        try:
            keys = [(0, self._private_key, None)]
            # Decrypt header to extract session keys using crypt4gh library
            session_keys, edit_list = crypt4gh.header.deconstruct(self.source, keys, sender_pubkey=None)

            if edit_list is not None:
                raise ValueError("Edit lists in Crypt4GH headers are not supported!")

            if not session_keys:
                raise ValueError("No session keys found in Crypt4GH header")

            self._ciphers = [ChaCha20Poly1305(bytes(key)) for key in session_keys]
            self._header_parsed = True
        except ValueError as e:
            if str(e) == "No supported encryption method":
                # crypt4gh raises this if the key opens no packet of the header
                raise DecryptionError("the private key does not open its Crypt4GH header") from e
            raise DecryptionError(f"Crypt4GH header error: {e}") from e


class Crypt4GHEncryptor(Transformer):
    """
    Crypt4GH Encryption Transformer.
    """

    NONCE_LENGTH = 12
    SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE
    CIPHER_DIFF = crypt4gh.lib.CIPHER_DIFF

    def __init__(
        self,
        source: Readable | None = None,
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
        self._cipher = ChaCha20Poly1305(self._session_key)
        self._header_sent = False

    def _fill_buffer(self) -> bytes:
        if not self._header_sent:
            self._header_sent = True
            return self._compose_header()

        if self.source is None:
            raise StreamConfigurationError("Stream source not set", stage=self.__class__.__name__)

        segment = _read_segment(self.source, self.SEGMENT_SIZE)
        if not segment:
            return b""

        nonce = os.urandom(self.NONCE_LENGTH)
        return nonce + self._cipher.encrypt(nonce, segment, None)

    def _compose_header(self) -> bytes:
        keys = [(0, self._sender_privkey, self._recipient_pubkey)]
        try:
            header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
            header_packets = crypt4gh.header.encrypt(header_content, keys)
            return crypt4gh.header.serialize(header_packets)
        except Exception as e:
            raise EncryptionError(f"Crypt4GH header cannot be encrypted: {e}") from e
