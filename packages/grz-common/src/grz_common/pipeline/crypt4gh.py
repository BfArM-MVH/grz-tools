import io
import os

import crypt4gh.header
import crypt4gh.lib
from nacl.bindings import crypto_aead_chacha20poly1305_ietf_decrypt, crypto_aead_chacha20poly1305_ietf_encrypt
from nacl.public import PrivateKey

from .base import TransformStream

# Constants
NONCE_LENGTH = 12
SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE  # 64KB
CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE  # 64KB + overhead
INT_WIDTH = 4  # Byte width for integers in header (version, count, length)
CRYPT4GH_MAGIC = b"crypt4gh"


class Crypt4GHDecryptor(TransformStream):
    def __init__(self, source: io.BufferedIOBase, private_key: bytes):
        super().__init__(source)
        self._private_key = private_key
        self._session_key: bytes | None = None
        self._header_parsed = False
        self._header_buffer = bytearray()
        self._cipher_residue = bytearray()
        self._header_length = 0

    def _fill_buffer(self) -> None:
        # Phase 1: Parse Header
        if not self._header_parsed:
            chunk = self.source.read(io.DEFAULT_BUFFER_SIZE)
            if not chunk:
                # If we hit EOF before parsing header, and buffer has ANY data, it's corrupt.
                # If buffer is empty, it's an empty file (technically invalid for c4gh).
                if len(self._header_buffer) > 0:
                    raise IOError("Unexpected EOF while reading Crypt4GH header")
                self._eof = True
                return

            self._header_buffer.extend(chunk)
            self._try_parse_header()
            return

        # Phase 2: Decrypt Body
        if len(self._cipher_residue) < CIPHER_SEGMENT_SIZE:
            chunk = self.source.read(io.DEFAULT_BUFFER_SIZE)
            if not chunk:
                self._eof = True
                self._flush_residue()
                return
            self._cipher_residue.extend(chunk)

        # Process full blocks
        while len(self._cipher_residue) >= CIPHER_SEGMENT_SIZE:
            segment = self._cipher_residue[:CIPHER_SEGMENT_SIZE]
            self._cipher_residue = self._cipher_residue[CIPHER_SEGMENT_SIZE:]
            self._output_buffer.extend(self._decrypt_segment(segment))

    def _try_parse_header(self) -> None:
        # Need at least Magic(8) + Ver(4) + Count(4) = 16 bytes
        if len(self._header_buffer) < 16:
            return

        try:
            with io.BytesIO(self._header_buffer) as f:
                # Strict Magic Check
                magic = f.read(8)
                if magic != CRYPT4GH_MAGIC:
                    # Fail fast if this isn't a Crypt4GH file
                    raise ValueError(f"Invalid Crypt4GH magic bytes: {magic!r}")

                _ = int.from_bytes(f.read(INT_WIDTH), "little")
                packet_count = int.from_bytes(f.read(INT_WIDTH), "little")

                for _ in range(packet_count):
                    # Check length field available
                    if f.tell() + INT_WIDTH > len(self._header_buffer):
                        return

                    pkt_len = int.from_bytes(f.read(INT_WIDTH), "little")

                    # Check packet body available
                    if f.tell() + pkt_len - INT_WIDTH > len(self._header_buffer):
                        return

                    f.seek(pkt_len - INT_WIDTH, os.SEEK_CUR)

                self._header_length = f.tell()

            # We have the full header. Parse it.
            header_bytes = self._header_buffer[: self._header_length]
            keys = [(0, self._private_key, None)]

            session_keys, _ = crypt4gh.header.deconstruct(io.BytesIO(header_bytes), keys, sender_pubkey=None)

            self._session_key = session_keys[0]
            self._header_parsed = True

            # Move remaining bytes (encrypted body) to residue
            self._cipher_residue.extend(self._header_buffer[self._header_length :])
            self._header_buffer = None

        except Exception as e:
            # Propagate specific validation errors
            raise IOError(f"Crypt4GH Header Error: {e}") from e

    def _decrypt_segment(self, segment: bytes) -> bytes:
        nonce = segment[:NONCE_LENGTH]
        return crypto_aead_chacha20poly1305_ietf_decrypt(segment[NONCE_LENGTH:], None, nonce, self._session_key)

    def _flush_residue(self) -> None:
        if self._cipher_residue:
            self._output_buffer.extend(self._decrypt_segment(bytes(self._cipher_residue)))
            self._cipher_residue.clear()


class Crypt4GHEncryptor(TransformStream):
    def __init__(self, source: io.BufferedIOBase, recipient_pubkey: bytes, sender_privkey: bytes | None = None):
        super().__init__(source)
        self._recipient_pubkey = recipient_pubkey
        self._sender_privkey = sender_privkey or bytes(PrivateKey.generate())
        self._session_key = os.urandom(32)
        self._header_sent = False
        self._plain_residue = bytearray()

    def _fill_buffer(self) -> None:
        if not self._header_sent:
            keys = [(0, self._sender_privkey, self._recipient_pubkey)]
            header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
            header_packets = crypt4gh.header.encrypt(header_content, keys)
            self._output_buffer.extend(crypt4gh.header.serialize(header_packets))
            self._header_sent = True

        chunk = self.source.read(io.DEFAULT_BUFFER_SIZE)
        if not chunk:
            self._eof = True
            if self._plain_residue:
                self._output_buffer.extend(self._encrypt_segment(bytes(self._plain_residue)))
                self._plain_residue.clear()
            return

        self._plain_residue.extend(chunk)
        while len(self._plain_residue) >= SEGMENT_SIZE:
            block = self._plain_residue[:SEGMENT_SIZE]
            self._plain_residue = self._plain_residue[SEGMENT_SIZE:]
            self._output_buffer.extend(self._encrypt_segment(block))

    def _encrypt_segment(self, plaintext: bytes) -> bytes:
        nonce = os.urandom(NONCE_LENGTH)
        return nonce + crypto_aead_chacha20poly1305_ietf_encrypt(plaintext, None, nonce, self._session_key)
