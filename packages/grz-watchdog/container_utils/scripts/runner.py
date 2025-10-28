import os
import subprocess
import sys
from pathlib import Path


def main():
    try:
        script_path = Path(__file__).resolve()
        project_root = script_path.parent.parent.parent.parent.parent
    except (IndexError, FileNotFoundError):
        print("Error: Could not determine project root directory.", file=sys.stderr)
        sys.exit(1)

    workdir = project_root / "packages/grz-watchdog/workdir"
    snakefile_path = project_root / "packages/grz-watchdog/workflow/Snakefile"
    snakemake_profile_path = project_root / "packages/grz-watchdog/workflow/profiles/default"
    pixi_manifest_path = project_root / "packages/grz-watchdog/pixi.toml"

    workdir = Path(os.environ.get("GRZ_WATCHDOG_WORKDIR", workdir))
    snakefile_path = Path(os.environ.get("GRZ_WATCHDOG_SNAKEFILE", snakefile_path))
    snakemake_profile_path = Path(os.environ.get("GRZ_WATCHDOG_SNAKEMAKE_PROFILE", snakemake_profile_path))

    workdir.mkdir(parents=True, exist_ok=True)

    snakemake_args = sys.argv[1:]

    command = [
        "pixi",
        "run",
        "--manifest-path",
        str(pixi_manifest_path),
        "--",
        "snakemake",
        "--snakefile",
        str(snakefile_path),
        "--workflow-profile",
        str(snakemake_profile_path),
        "--directory",
        str(workdir),
        *snakemake_args,
    ]

    print("---- grz-watchdog runner ----", flush=True)
    print(f"project root:     {project_root}", flush=True)
    print(f"work directory:   {workdir}", flush=True)
    print(f"pixi manifest:    {pixi_manifest_path}", flush=True)
    print(f"forwarded args:   {' '.join(snakemake_args)}", flush=True)
    print("------------------------------", flush=True)

    try:
        process = subprocess.run(command, check=True)  # noqa: S603
        sys.exit(process.returncode)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)
    except FileNotFoundError:
        print("Error: 'pixi' not found. Is pixi installed and in your PATH?", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
