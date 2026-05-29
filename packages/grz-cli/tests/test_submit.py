from pathlib import Path

import click.testing
import yaml


MINIMAL_CONFIG = {
    "s3": {
        "endpoint_url": "https://example.invalid",
        "bucket": "test-bucket",
        "access_key": "AKIA-test",
        "secret": "test-secret",
    },
}


def _write_config(tmp_path: Path) -> Path:
    config_path = tmp_path / "config.yaml"
    config_path.write_text(yaml.dump(MINIMAL_CONFIG))
    return config_path


def _make_submission_dir(tmp_path: Path) -> Path:
    submission_dir = tmp_path / "submission"
    for sub in ("metadata", "files", "encrypted_files", "logs"):
        (submission_dir / sub).mkdir(parents=True)
    return submission_dir


def test_submit_invokes_subcommands(tmp_path, monkeypatch):
    """submit should parse config + run validate, encrypt, upload in order."""
    from grz_cli.cli import build_cli
    from grz_cli.commands import submit as submit_module

    calls: list[str] = []

    def fake_validate(*args, **kwargs):
        calls.append("validate")

    def fake_encrypt(*args, **kwargs):
        calls.append("encrypt")

    def fake_upload(*args, **kwargs):
        calls.append("upload")

    monkeypatch.setattr(submit_module.validate, "callback", fake_validate)
    monkeypatch.setattr(submit_module.encrypt, "callback", fake_encrypt)
    monkeypatch.setattr(submit_module.upload, "callback", fake_upload)
    monkeypatch.setattr(submit_module, "check_version_and_exit_if_needed", lambda _s3: None)

    config_path = _write_config(tmp_path)
    submission_dir = _make_submission_dir(tmp_path)

    runner = click.testing.CliRunner()
    result = runner.invoke(
        build_cli(),
        [
            "submit",
            "--submission-dir",
            str(submission_dir),
            "--config-file",
            str(config_path),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
    assert calls == ["validate", "encrypt", "upload"]


def test_submit_accepts_multiple_config_files(tmp_path, monkeypatch):
    """Regression: repeated --config-file once produced a tuple that crashed UploadConfig.from_path."""
    from grz_cli.cli import build_cli
    from grz_cli.commands import submit as submit_module

    monkeypatch.setattr(submit_module.validate, "callback", lambda *a, **kw: None)
    monkeypatch.setattr(submit_module.encrypt, "callback", lambda *a, **kw: None)
    monkeypatch.setattr(submit_module.upload, "callback", lambda *a, **kw: None)
    monkeypatch.setattr(submit_module, "check_version_and_exit_if_needed", lambda _s3: None)

    base_config_path = tmp_path / "base.yaml"
    base_config_path.write_text(yaml.dump({"s3": {"endpoint_url": "https://example.invalid", "bucket": "x"}}))

    override_config_path = tmp_path / "override.yaml"
    override_config_path.write_text(yaml.dump({"s3": {"access_key": "k", "secret": "s"}}))

    submission_dir = _make_submission_dir(tmp_path)

    runner = click.testing.CliRunner()
    result = runner.invoke(
        build_cli(),
        [
            "submit",
            "--submission-dir",
            str(submission_dir),
            "--config-file",
            str(base_config_path),
            "--config-file",
            str(override_config_path),
        ],
        catch_exceptions=False,
    )

    assert result.exit_code == 0, result.output
