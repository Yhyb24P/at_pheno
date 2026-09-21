import importlib.util
from pathlib import Path

VALID = {
    "trait_id": "t1", "trait_name": "DTF-s", "study": "1001 Genomes & easyGWAS",
    "source_publication_doi": "10.1105/tpc.16.00551", "arapheno_doi": "10.21958/phenotype:703",
    "measurement_definition": "days until visible flowering buds",
    "environment": "chamber 16C 16:8 RBD x4", "experimental_unit": "plant",
    "replicate_design": "4-replicate randomized block", "database_value_level": "published accession-level value",
    "upstream_aggregation": "not documented", "local_aggregation": "none",
    "transformation_of_selected_target": "none", "finite_panel_n": "9",
    "duplicate_accession_n": "0", "duplicate_resolution": "none",
    "target_status": "PUBLISHED_ACCESSION_VALUE_USABLE", "candidate_role": "reference",
    "values_sha256": "a" * 64, "primary_evidence": "snapshot_703.json", "notes": "",
}

PENDING = {**VALID, "trait_id": "t2", "trait_name": "rosette-s",
           "duplicate_accession_n": "2",
           "duplicate_resolution": "two equal-value duplicate rows; source resolution documented, not re-frozen",
           "target_status": "REQUIRES_DUPLICATE_RESOLUTION"}


def _lib():
    path = Path(__file__).parents[1] / "scripts" / "trait_resolution.py"
    spec = importlib.util.spec_from_file_location("trait_resolution", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_tsv(tmp_path, rows):
    lib = _lib()
    lines = ["\t".join(lib.COLUMNS)]
    for row in rows:
        lines.append("\t".join(row[c] for c in lib.COLUMNS))
    path = tmp_path / "registry.tsv"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path, lib


def test_minimal_fixture_roundtrips_and_reports_no_violations(tmp_path):
    path, lib = _write_tsv(tmp_path, [VALID, PENDING, {**VALID, "trait_id": "t3", "target_status": "UNRESOLVED_SOURCE"}])
    rows = lib.load_tsv(path)
    assert [r["trait_id"] for r in rows] == ["t1", "t2", "t3"]
    assert lib.validate(rows) == []
    assert [r["target_status"] for r in lib.usable_rows(rows)] == ["PUBLISHED_ACCESSION_VALUE_USABLE"]
    assert [r["trait_id"] for r in lib.duplicate_pending_rows(rows)] == ["t2"]


def test_usable_with_duplicate_accessions_is_rejected(tmp_path):
    path, lib = _write_tsv(tmp_path, [{**VALID, "duplicate_accession_n": "2", "duplicate_resolution": "pending"}])
    errors = lib.validate(lib.load_tsv(path))
    assert len(errors) == 1 and "USABLE with duplicate accessions" in errors[0]


def test_duplicate_rows_require_documented_resolution(tmp_path):
    path, lib = _write_tsv(tmp_path, [{**PENDING, "duplicate_resolution": "  "}, {**VALID, "duplicate_accession_n": "0", "duplicate_resolution": "none"}])
    errors = lib.validate(lib.load_tsv(path))
    assert any("t2" in e for e in errors) and "documented duplicate_resolution" in "; ".join(errors)


def test_invalid_status_and_missing_sha_are_rejected(tmp_path):
    path, lib = _write_tsv(tmp_path, [
        {**VALID, "trait_id": "a", "target_status": "MAYBE_USABLE"},
        {**VALID, "trait_id": "b", "values_sha256": "0" * 32},
        {**VALID, "trait_id": "c", "target_status": "REQUIRES_DUPLICATE_RESOLUTION", "duplicate_accession_n": "0"},
    ])
    errors = lib.validate(lib.load_tsv(path))
    assert any("a" in e and "invalid target_status" in e for e in errors)
    assert any("b" in e and "64-hex" in e for e in errors)
    assert len(errors) == 2


def test_real_registry_is_internally_consistent():
    lib = _lib()
    rows = lib.load_tsv(Path(__file__).parents[1] / "data/manifests/trait_resolution_v1.tsv")
    assert len(rows) == 7
    assert lib.validate(rows) == []
    roles = {r["trait_id"]: r["candidate_role"] for r in rows}
    # Mapping assertions per the delivery spec: 704=RL, 705=CL, 390=rosetteDM, 748=stomata density
    assert roles["704"] == "morphology_sensitivity_not_independent_replication"
    assert roles["705"] == "morphology_sensitivity_not_independent_replication"
    assert {r["trait_id"] for r in lib.usable_rows(rows)} == {"703", "704", "705", "748"}
    assert [r["trait_id"] for r in lib.duplicate_pending_rows(rows)] == ["390", "261", "262"]
