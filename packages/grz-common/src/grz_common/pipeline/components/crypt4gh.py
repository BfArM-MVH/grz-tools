import io
import os

import crypt4gh.header
import crypt4gh.lib
from nacl.bindings import crypto_aead_chacha20poly1305_ietf_encrypt
from nacl.public import PrivateKey

from . import TransformStream


class Crypt4GHDecryptor(TransformStream):
    CIPHER_DIFF = crypt4gh.lib.CIPHER_DIFF
    CIPHER_SEGMENT_SIZE = crypt4gh.lib.CIPHER_SEGMENT_SIZE

    def __init__(self, source: io.BufferedIOBase, private_key: bytes):
        super().__init__(source)
        self._private_key = private_key
        self._session_keys: bytes | None = None
        self._header_parsed = False
        self._cipher_residue = bytearray()

    def _fill_buffer(self) -> bytes:
        if not self._header_parsed:
            self._read_header()

        return self._decrypt_next_segment()

    def _read_header(self) -> None:
        """Parse the Crypt4GH header to extract the session key for decryption."""
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

    def _decrypt_next_segment(self) -> bytes:
        """Decrypt and return the next segment, or None if EOF."""
        ciphersegment = self.source.read(self.CIPHER_SEGMENT_SIZE)
        if not ciphersegment:
            # EOF reached
            return b""
        elif len(ciphersegment) <= self.CIPHER_DIFF:
            # This means we have a truncated segment:
            # Even an empty plaintext segment would have CIPHER_DIFF bytes of overhead
            raise ValueError("Truncated cipher segment")
        else:
            return crypt4gh.lib.decrypt_block(ciphersegment, self._session_keys)


class Crypt4GHEncryptor(TransformStream):
    NONCE_LENGTH = 12
    SEGMENT_SIZE = crypt4gh.lib.SEGMENT_SIZE

    def __init__(self, source: io.BufferedIOBase, recipient_pubkey: bytes, sender_privkey: bytes | None = None):
        super().__init__(source)
        self._recipient_pubkey = recipient_pubkey
        self._sender_privkey = sender_privkey or bytes(PrivateKey.generate())
        self._session_key = os.urandom(32)
        self._header_sent = False
        self._plain_residue = bytearray()

    def _fill_buffer(self) -> bytes:
        if not self._header_sent:
            header = self._compose_header()
            self._header_sent = True

            return header

        return self._encrypt_next_segment()

    def _compose_header(self) -> bytes:
        keys = [(0, self._sender_privkey, self._recipient_pubkey)]

        header_content = crypt4gh.header.make_packet_data_enc(0, self._session_key)
        header_packets = crypt4gh.header.encrypt(header_content, keys)

        return crypt4gh.header.serialize(header_packets)

    def _encrypt_next_segment(self) -> bytes:
        segment = self.source.read(self.SEGMENT_SIZE)

        if len(segment) == 0:
            # EOF reached
            return b""

        nonce = os.urandom(self.NONCE_LENGTH)
        ciphertext = crypto_aead_chacha20poly1305_ietf_encrypt(segment, None, nonce, self._session_key)
        return nonce + ciphertext
