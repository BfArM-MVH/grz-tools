# packages/grz-check-provider/hatch_hooks.py
import json
import os
import shutil
import subprocess
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface


class CustomBuildHook(BuildHookInterface):
    """
    A Hatch build hook that compiles a Rust binary and includes it in the wheel.
    """

    def _get_cargo_target_dir(self, rust_project_path: Path) -> str:
        """
        Asks Cargo where it will place build artifacts.

        This respects `CARGO_TARGET_DIR` and local `.cargo/config.toml` files.
        """
        try:
            # --no-deps is faster as we only need workspace info.
            # --format-version 1 ensures stable JSON output.
            command = ["cargo", "metadata", "--no-deps", "--format-version", "1"]
            self.app.display_info(f"Querying cargo metadata with: `{' '.join(command)}`")

            result = subprocess.run(
                command,
                cwd=str(rust_project_path),
                check=True,
                capture_output=True,
                text=True,
            )
            metadata = json.loads(result.stdout)
            target_dir = metadata["target_directory"]
            self.app.display_info(f"Cargo target directory detected at: {target_dir}")
            return target_dir

        except (subprocess.CalledProcessError, FileNotFoundError, KeyError) as e:
            self.app.display_warning(
                f"Could not determine cargo target directory via `cargo metadata`: {e}. "
                f"Falling back to default '{rust_project_path / 'target'}'"
            )
            # Fallback for safety, though it shouldn't be needed in normal environments.
            return str(rust_project_path / "target")

    def initialize(self, version, build_data):
        # This hook is only for building wheels, not sdists
        if self.target_name != "wheel":
            return

        self.app.display_info("Custom Rust Build Hook (grz-check)")

        root_dir = Path(self.root)
        rust_project_path = root_dir / "rust_src"
        binary_name = "grz-check"
        target_in_src_dir = root_dir / "src" / "grz_check_provider" / "bin"

        cargo_target_directory = Path(self._get_cargo_target_dir(rust_project_path))

        build_target = os.environ.get("CARGO_BUILD_TARGET")

        if build_target:
            self.app.display_info(f"Cross-compilation detected for target: {build_target}")
            build_command_base = ["cross"]
            build_args = ["--target", build_target]
            source_binary_path = cargo_target_directory / build_target / "release" / binary_name
        else:
            self.app.display_info("Native build detected (CARGO_BUILD_TARGET not set).")
            build_command_base = ["cargo"]
            build_args = []
            source_binary_path = cargo_target_directory / "release" / binary_name

        if not shutil.which(build_command_base[0]):
            raise FileNotFoundError(f"'{build_command_base[0]}' not found. Is it installed and in PATH?")

        build_command = [*build_command_base, "build", "--release", *build_args]

        self.app.display_info(f"Running command: {' '.join(build_command)}")
        try:
            subprocess.run(build_command, cwd=str(rust_project_path), check=True)
        except subprocess.CalledProcessError as e:
            self.app.display_critical("Rust build failed")
            raise e

        self.app.display_info("Rust build successful")

        target_in_src_dir.mkdir(parents=True, exist_ok=True)
        final_binary_path = target_in_src_dir / binary_name

        self.app.display_info(f"Copying binary from '{source_binary_path}' to '{final_binary_path}'")
        if not source_binary_path.exists():
            raise FileNotFoundError(
                f"The compiled rust binary was not found at the expected path: {source_binary_path}. "
                "Please check the build logs."
            )
        shutil.copy(source_binary_path, final_binary_path)

        # tell hatch to include the binary in the wheel.
        relative_artifact_path = final_binary_path.relative_to(root_dir / "src")
        build_data["force_include"][str(final_binary_path)] = str(relative_artifact_path)

        self.app.display_info(f"Included '{final_binary_path}' in wheel at '{relative_artifact_path}'")
