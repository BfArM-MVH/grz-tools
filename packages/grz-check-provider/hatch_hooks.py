import os
import platform
import shutil
import subprocess
from pathlib import Path

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

RUST_GLIBC_TARGET_MAP = {
    "x86_64": "x86_64-unknown-linux-gnu",
    "aarch64": "aarch64-unknown-linux-gnu",
}
RUST_MUSL_TARGET_MAP = {
    "x86_64": "x86_64-unknown-linux-musl",
    "aarch64": "aarch64-unknown-linux-musl",
}


class CustomBuildHook(BuildHookInterface):
    def initialize(self, version, build_data):
        if self.target_name != "wheel":
            return

        print("--- Starting Rust cross-compilation for grz-check ---")

        # Get target architecture from cibuildwheel or fall back to host machine
        target_arch_str = os.environ.get("CIBW_ARCHS", platform.machine())
        print(f"--- Build architecture detected: {target_arch_str} ---")

        if "musl" in target_arch_str:
            print("--- MUSL target detected ---")
            arch_key = target_arch_str.split("_", 1)[1]  # "musllinux_x86_64" -> "x86_64"
            rust_target = RUST_MUSL_TARGET_MAP.get(arch_key)
        else:
            print("--- GLIBC target detected ---")
            rust_target = RUST_GLIBC_TARGET_MAP.get(target_arch_str)

        if not rust_target:
            raise RuntimeError(f"Unsupported architecture for Rust build: {target_arch_str}")

        print(f"--- Rust target triple: {rust_target} ---")

        if not shutil.which("cargo"):
            raise FileNotFoundError("Cargo command not found. Please install Rust.")

        root_dir = Path(self.root)
        rust_project_path = root_dir / "rust_src"
        rust_target_dir = rust_project_path / "target"
        build_command = ["cargo", "build", "--release", f"--target={rust_target}", f"--target-dir={rust_target_dir}"]

        try:
            print(f"--- Running command: {' '.join(build_command)} ---")
            subprocess.run(build_command, cwd=str(rust_project_path), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print("!!! Rust build failed !!!")
            print(f"--- STDERR ---\n{e.stderr}")
            print(f"--- STDOUT ---\n{e.stdout}")
            raise e

        print("--- Rust build successful ---")

        # The binary name from your Cargo.toml
        binary_name = "grz-check"
        source_path = rust_target_dir / rust_target / "release" / binary_name

        # This is the path inside the Python package's source tree
        target_dir = root_dir / "src" / "grz_check_provider" / "bin"
        target_dir.mkdir(exist_ok=True)
        target_path = target_dir / binary_name

        print(f"--- Copying binary from {source_path} to {target_path} ---")
        shutil.copy(source_path, target_path)

        # Explicitly tell hatch to include the compiled binary in the wheel
        build_data["force_include"][str(target_path)] = str(target_path.relative_to(root_dir / "src"))
