import hashlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]


def sha256(path):
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    return hasher.hexdigest()


def test_current_v1_analysis_plan_still_matches_original_g0_freeze():
    manifest = json.loads((ROOT / "data/manifests/confirmatory_v1/freeze_manifest.json").read_text())
    expected = manifest["artifact_bindings"]["analysis_plan"]
    path = ROOT / expected["path"]
    assert path.stat().st_size == expected["bytes"]
    assert sha256(path) == expected["sha256"]


def test_current_v1_design_document_is_byte_identical_to_its_g0_version():
    path = ROOT / "docs/确认性设计_v1.md"
    assert sha256(path) == "d5570ee4c8e3505eec212460165edffa86b5296b9ecc1f1de9f390c2457ca92f"


def test_all_original_g0_bindings_revalidate_when_full_artifacts_are_present():
    manifest = json.loads((ROOT / "data/manifests/confirmatory_v1/freeze_manifest.json").read_text())
    bindings = list(manifest["artifact_bindings"].values()) + manifest["canonical_split_bindings"]
    unavailable = [record["path"] for record in bindings if not (ROOT / record["path"]).exists()]
    if unavailable:
        pytest.skip("full local G0 artifacts unavailable in this checkout: " + ", ".join(unavailable))
    for record in bindings:
        path = ROOT / record["path"]
        assert path.stat().st_size == record["bytes"]
        assert sha256(path) == record["sha256"]
