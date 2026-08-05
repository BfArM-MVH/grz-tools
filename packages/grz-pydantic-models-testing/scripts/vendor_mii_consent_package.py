#!/usr/bin/env python3
"""
Vendor the MII consent artefacts the consent tests check the model against.

Downloads one published version of ``de.medizininformatikinitiative.kerndatensatz.consent`` and
writes the CodeSystems and the Consent profile into ``example_terminology``, named after the version
stated inside each resource. Run this by hand when the MII publishes a release; the test suite never
reaches the network itself.

Usage::

    python scripts/vendor_mii_consent_package.py 2026.0.1

Then run the consent tests: whatever the model does not yet cover fails by file name.
"""

import argparse
import io
import json
import sys
import tarfile
import urllib.request
from pathlib import Path

REGISTRY = "https://packages2.fhir.org/packages/de.medizininformatikinitiative.kerndatensatz.consent"

#: Artefact inside the FHIR package -> file name prefix used in example_terminology.
ARTEFACTS = {
    "CodeSystem-MiiConsentPolicyCodeSystem.json": "mii-cs-consent-policy",
    "CodeSystem-MiiConsentVersionModuleCodeSystem.json": "mii-cs-consent-version-modules",
    "Profile_MII_Consent_Einwilligung.json": "mii-pr-consent-einwilligung",
}

TERMINOLOGY_DIR = Path(__file__).resolve().parent.parent / "src/grz_pydantic_models_testing/example_terminology"


def download(package_version: str) -> dict[str, dict]:
    """Fetch a package and return its artefacts of interest, keyed by file name.

    :param package_version: published package version, e.g. '2026.0.1'.
    :returns: parsed JSON of every artefact the package ships, missing ones omitted.
    """
    url = f"{REGISTRY}/{package_version}"
    print(f"fetching {url}")
    with urllib.request.urlopen(url) as response:  # noqa: S310 - fixed https registry
        archive = io.BytesIO(response.read())

    found = {}
    with tarfile.open(fileobj=archive, mode="r:gz") as tar:
        shipped = set(tar.getnames())
        for name in ARTEFACTS:
            # not every package ships every artefact: the version and module CodeSystem appears in
            # 2026.0.0 only, and the 2026.0.1 release candidates dropped it again
            if f"package/{name}" not in shipped:
                continue
            with tar.extractfile(f"package/{name}") as member:
                found[name] = json.load(member)
    return found


def vendor(artefacts: dict[str, dict], force: bool) -> int:
    """Write artefacts into the terminology directory, reporting what changed.

    :param artefacts: parsed artefacts keyed by their file name in the package.
    :param force: overwrite a vendored file whose contents differ from upstream.
    :returns: process exit code.
    """
    conflicts = 0
    for name, prefix in ARTEFACTS.items():
        resource = artefacts.get(name)
        if resource is None:
            print(f"  {prefix}: not shipped by this package")
            continue

        target = TERMINOLOGY_DIR / f"{prefix}.{resource['version']}.json"
        serialized = json.dumps(resource, indent=2, ensure_ascii=False) + "\n"

        if target.exists():
            if json.loads(target.read_text()) == resource:
                print(f"  {target.name}: already vendored, unchanged")
                continue
            conflicts += 1
            if not force:
                print(f"  {target.name}: DIFFERS from upstream, kept. Re-run with --force to overwrite")
                continue
            print(f"  {target.name}: differed from upstream, overwritten")
        else:
            print(f"  {target.name}: added")

        target.write_text(serialized)

    return 1 if conflicts and not force else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("package_version", help="published package version, e.g. 2026.0.1")
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite a vendored file whose contents differ, i.e. upstream republished a version",
    )
    args = parser.parse_args()

    exit_code = vendor(download(args.package_version), args.force)
    print(f"\nvendored into {TERMINOLOGY_DIR}")
    print("next: uv run pytest packages/grz-pydantic-models")
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
