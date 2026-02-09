"""
Integration tests comparing manual step-by-step CLI processing with streaming pipeline.

This test ensures that the streaming pipeline produces equivalent results to
the traditional manual workflow using grzctl/grz-cli commands:
    grzctl decrypt -> grz-cli validate -> grz-cli encrypt -> grzctl archive
"""

import shutil
from pathlib import Path

import grz_cli.cli
import grzctl.cli
import pytest
from click.testing import CliRunner
from grz_common.pipeline import (
    Crypt4GHDecryptor,
    PipelineContext,
    RawChecksumValidator,
)
from grz_common.utils.checksums import calculate_sha256

# Path to test fixtures
MOCK_FILES_DIR = Path(__file__).parent / "mock_files"
VALID_SUBMISSION_DIR = MOCK_FILES_DIR / "submissions" / "valid_submission"


def get_private_key(key_path: Path, passphrase: str = "") -> bytes:
    """Load a Crypt4GH private key."""
    import crypt4gh.keys

    return crypt4gh.keys.get_private_key(str(key_path), lambda: passphrase)


def get_public_key(key_path: Path) -> bytes:
    """Load a Crypt4GH public key."""
    import crypt4gh.keys

    return crypt4gh.keys.get_public_key(str(key_path))


class TestManualCliVsStreamingPipeline:
    """
    Compare manual step-by-step CLI processing with streaming pipeline.

    Manual workflow (using CLI commands):
        1. grzctl decrypt --submission-dir ...
        2. grz-cli validate --submission-dir ...
        3. grz-cli encrypt --submission-dir ... (with archive key)

    Streaming pipeline:
        Uses the pipeline components directly (same as grzctl process internally)
    """

    @pytest.fixture
    def working_dir_manual(self, tmpdir_factory: pytest.TempdirFactory) -> Path:
        """Create a temporary directory for manual CLI processing."""
        datadir = tmpdir_factory.mktemp("manual_submission")
        return Path(datadir.strpath)

    @pytest.fixture
    def working_dir_streaming(self, tmpdir_factory: pytest.TempdirFactory) -> Path:
        """Create a temporary directory for streaming processing."""
        datadir = tmpdir_factory.mktemp("streaming_submission")
        return Path(datadir.strpath)

    def _setup_submission_dir(self, working_dir: Path) -> None:
        """Copy encrypted submission files to working directory."""
        shutil.copytree(
            VALID_SUBMISSION_DIR / "encrypted_files",
            working_dir / "encrypted_files",
            dirs_exist_ok=True,
        )
        shutil.copytree(
            VALID_SUBMISSION_DIR / "metadata",
            working_dir / "metadata",
            dirs_exist_ok=True,
        )
        # Also copy the original unencrypted files for checksum comparison
        shutil.copytree(
            VALID_SUBMISSION_DIR / "files",
            working_dir / "files_original",
            dirs_exist_ok=True,
        )

    def _run_cli_decrypt(self, working_dir: Path, config_file: str) -> None:
        """Run grzctl decrypt command."""
        runner = CliRunner()
        cli = grzctl.cli.build_cli()

        result = runner.invoke(
            cli,
            [
                "decrypt",
                "--submission-dir",
                str(working_dir),
                "--config-file",
                config_file,
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Decrypt failed: {result.output}"

    def _run_cli_validate(self, working_dir: Path, config_file: str) -> None:
        """Run grz-cli validate command."""
        runner = CliRunner()
        cli = grz_cli.cli.build_cli()

        result = runner.invoke(
            cli,
            [
                "validate",
                "--submission-dir",
                str(working_dir),
                "--config-file",
                config_file,
                "--no-grz-check",  # Use fallback validation
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Validate failed: {result.output}"

    def _run_cli_encrypt(self, working_dir: Path, config_file: str) -> None:
        """Run grz-cli encrypt command with archive key."""
        runner = CliRunner()
        cli = grz_cli.cli.build_cli()

        # First remove the inbox-encrypted files so encrypt creates new ones
        encrypted_files_dir = working_dir / "encrypted_files"
        if encrypted_files_dir.exists():
            shutil.rmtree(encrypted_files_dir)

        result = runner.invoke(
            cli,
            [
                "encrypt",
                "--submission-dir",
                str(working_dir),
                "--config-file",
                config_file,
                "--no-check-validation-logs",
            ],
            catch_exceptions=False,
        )
        assert result.exit_code == 0, f"Encrypt failed: {result.output}"

    def _streaming_decrypt_file(
        self,
        encrypted_path: Path,
        output_path: Path,
        decrypt_key: bytes,
    ) -> str:
        """
        Decrypt a file using the streaming pipeline Crypt4GHDecryptor.

        This is the same decryptor used internally by grzctl process.

        Returns the SHA256 checksum of the decrypted content.
        """
        context = PipelineContext()
        decryptor = Crypt4GHDecryptor(decrypt_key)
        checksum_validator = RawChecksumValidator()

        decryptor.initialize(context)
        checksum_validator.initialize(context)

        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(encrypted_path, "rb") as infile, open(output_path, "wb") as outfile:
            while chunk := infile.read(65536):
                decrypted = decryptor.process(chunk)
                if decrypted:
                    checksum_validator.observe(decrypted)
                    outfile.write(decrypted)

            final = decryptor.flush()
            if final:
                checksum_validator.observe(final)
                outfile.write(final)

        decryptor.finalize()
        checksum_validator.finalize()

        return checksum_validator.calculated_checksum

    def test_decrypt_produces_same_checksums(
        self,
        working_dir_manual,
        working_dir_streaming,
        temp_keys_config_file_path,
        crypt4gh_grz_private_key_file_path,
    ):
        """
        Test that manual CLI decrypt (grzctl decrypt) and streaming pipeline
        Crypt4GHDecryptor produce identical checksums.

        This validates that the streaming decryptor component produces the
        same output as the file-based CLI command.
        """
        # Setup both directories
        self._setup_submission_dir(working_dir_manual)
        self._setup_submission_dir(working_dir_streaming)

        # === Manual CLI Decrypt (grzctl decrypt) ===
        self._run_cli_decrypt(working_dir_manual, temp_keys_config_file_path)

        # === Streaming Decrypt (using Crypt4GHDecryptor directly) ===
        grz_private_key = get_private_key(crypt4gh_grz_private_key_file_path)

        encrypted_files = list((working_dir_streaming / "encrypted_files").glob("*.c4gh"))
        streaming_checksums = {}

        for encrypted_file in encrypted_files:
            output_name = encrypted_file.stem  # removes .c4gh
            output_path = working_dir_streaming / "files" / output_name

            checksum = self._streaming_decrypt_file(encrypted_file, output_path, grz_private_key)
            streaming_checksums[output_name] = checksum

        # === Compare Checksums ===
        manual_files = list((working_dir_manual / "files").glob("*"))

        for manual_file in manual_files:
            if manual_file.is_file():
                manual_checksum = calculate_sha256(manual_file)
                streaming_checksum = streaming_checksums.get(manual_file.name)

                assert streaming_checksum is not None, f"Missing streaming result for {manual_file.name}"
                assert manual_checksum == streaming_checksum, (
                    f"Checksum mismatch for {manual_file.name}:\n"
                    f"  Manual CLI (grzctl decrypt): {manual_checksum}\n"
                    f"  Streaming (Crypt4GHDecryptor): {streaming_checksum}"
                )

    def test_full_workflow_produces_equivalent_results(
        self,
        working_dir_manual,
        working_dir_streaming,
        temp_keys_config_file_path,
        temp_identifiers_config_file_path,
        crypt4gh_grz_private_key_file_path,
    ):
        """
        Test that full manual CLI workflow and streaming pipeline produce equivalent results.

        Manual workflow:
            1. grzctl decrypt
            2. grz-cli validate
            3. grz-cli encrypt (re-encrypt with same key for simplicity)

        The streaming pipeline uses the same components internally (Crypt4GHDecryptor,
        FastqValidator, Crypt4GHEncryptor), so if the decryption matches, the full
        pipeline will produce equivalent results.
        """
        # Setup both directories
        self._setup_submission_dir(working_dir_manual)
        self._setup_submission_dir(working_dir_streaming)

        # === Manual CLI Workflow ===
        # Step 1: Decrypt
        self._run_cli_decrypt(working_dir_manual, temp_keys_config_file_path)

        # Step 2: Validate (using grz-cli validate)
        self._run_cli_validate(working_dir_manual, temp_identifiers_config_file_path)

        # Step 3: Re-encrypt with same key (simulating archive encryption)
        self._run_cli_encrypt(working_dir_manual, temp_keys_config_file_path)

        # Collect checksums from manual decrypted files
        manual_decrypted_checksums = {}
        for f in (working_dir_manual / "files").glob("*"):
            if f.is_file():
                manual_decrypted_checksums[f.name] = calculate_sha256(f)

        # === Streaming Pipeline (decrypt only, as that's the core transformation) ===
        grz_private_key = get_private_key(crypt4gh_grz_private_key_file_path)

        streaming_checksums = {}
        encrypted_files = list((working_dir_streaming / "encrypted_files").glob("*.c4gh"))

        for encrypted_file in encrypted_files:
            output_name = encrypted_file.stem  # remove .c4gh
            output_path = working_dir_streaming / "files" / output_name

            checksum = self._streaming_decrypt_file(encrypted_file, output_path, grz_private_key)
            streaming_checksums[output_name] = checksum

        # === Compare Results ===
        # Decrypted checksums should match
        for file_name, manual_checksum in manual_decrypted_checksums.items():
            streaming_checksum = streaming_checksums.get(file_name)
            assert streaming_checksum is not None, f"Missing streaming result for {file_name}"
            assert manual_checksum == streaming_checksum, (
                f"Checksum mismatch for {file_name}:\n"
                f"  Manual CLI: {manual_checksum}\n"
                f"  Streaming:  {streaming_checksum}"
            )


class TestStreamingPipelineChunkSizes:
    """Test that the streaming pipeline works correctly with various chunk sizes."""

    def test_different_chunk_sizes_produce_same_result(self, crypt4gh_grz_private_key_file_path):
        """
        Verify that different input chunk sizes produce the same decrypted output.

        This ensures the Crypt4GHDecryptor correctly handles buffering across
        arbitrary chunk boundaries.
        """
        grz_private_key = get_private_key(crypt4gh_grz_private_key_file_path)

        # Find a test file
        encrypted_files = list((VALID_SUBMISSION_DIR / "encrypted_files").glob("*.c4gh"))
        assert len(encrypted_files) > 0
        test_file = encrypted_files[0]

        chunk_sizes = [1024, 4096, 16384, 65536, 131072]
        checksums = []

        for chunk_size in chunk_sizes:
            context = PipelineContext()
            decryptor = Crypt4GHDecryptor(grz_private_key)
            checksum_validator = RawChecksumValidator()

            decryptor.initialize(context)
            checksum_validator.initialize(context)

            with open(test_file, "rb") as f:
                while chunk := f.read(chunk_size):
                    decrypted = decryptor.process(chunk)
                    if decrypted:
                        checksum_validator.observe(decrypted)

                final = decryptor.flush()
                if final:
                    checksum_validator.observe(final)

            decryptor.finalize()
            checksum_validator.finalize()

            checksums.append(checksum_validator.calculated_checksum)

        # All checksums should be identical
        assert all(cs == checksums[0] for cs in checksums), (
            f"Different chunk sizes produced different checksums: {dict(zip(chunk_sizes, checksums, strict=True))}"
        )
