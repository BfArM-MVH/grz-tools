"""The S3 permissions in the README of grzctl are the ones the end-to-end tests grant.

This test needs no container, so it runs with the default test run.
"""

import re
from pathlib import Path

from . import services

_README = Path(__file__).parents[2] / "packages" / "grzctl" / "README.md"
_CREDENTIALS = {
    "le": "`s3` of the LE's grz-cli config",
    "inbox": "`leistungserbringer.<LE>.inbox_buckets.<inbox>`",
    "interrogation": "`archives.interrogation.s3`",
    "consented": "`archives.consented.s3`",
    "non_consented": "`archives.non_consented.s3`",
}
_BUCKETS = {
    services.INBOX: "inbox",
    services.ARCHIVE_INTERROGATION: "interrogation",
    services.ARCHIVE_CONSENTED: "consented archive",
    services.ARCHIVE_NON_CONSENTED: "non-consented archive",
}


def _documented() -> dict[tuple[str, str], set[str]]:
    section = _README.read_text().split("## S3 permissions\n", 1)[1].split("\n## ", 1)[0]
    rows = [line for line in section.splitlines() if line.startswith("| `")]
    documented = {}
    for row in rows:
        credentials, bucket, actions = (cell.strip() for cell in row.strip("|").split("|"))
        documented[credentials, bucket] = set(re.findall(r"`(s3:\w+)`", actions))
    return documented


def test_the_readme_lists_the_permissions_of_the_end_to_end_tests():
    granted = {
        (_CREDENTIALS[role], _BUCKETS[bucket]): set(actions)
        for role, buckets in services.PERMISSIONS.items()
        for bucket, actions in buckets.items()
    }
    assert _documented() == granted
