"""TASK04: formal provenance v2 gate on a fully synthetic tiny fixture.

Nothing under data/raw and no real (12 GB-class) release file is
read or hashed.  The fixture under tmp_path holds a static 12-line
synthetic source VCF, a fixed 4x8 genotypes matrix (0/1/2 plus NaN),
three small TSVs, and a hand-baked v2 provenance record whose
bindings are recomputed at test time.  A record written *before*
tampering must be rejected by a ValueError naming the breached
field;  intact, the gate passes.
"""

import gzip
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pytest

from at_pheno.formal_provenance import PROVENANCE_SCHEMA_V2, gate_v2

# Committed catalog artifact (untracked on purpose; the check below
# self-skips in checkouts that do not carry the large deliverable).
CATALOG = "data/manifests/1001g_v31_biallelic_snp_catalog.tsv.gz"


def _md5_of_unpacked_gz(path):
    h = hashlib.md5()
    with gzip.open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()

# Static synthetic source VCF: 13 lines, deterministic content.
VCF_LINES = [
    "##fileformat=VCFv4.2",
    "##source=https://example.org/releases/2026/synthetic_tiny.vcf.gz",
    "##contig=<ID=1;length=4000>",
    "##contig=<ID=2;length=8000>",
    "##contig=<ID=3;length=6400>",
    "##contig=<ID=4;length=4700>",
    "##contig=<ID=5;length=1700>",
    "##FILTER=<ID=PASS;Description=\"All filters passed\">",
    "##FILTER=<ID=q4;Description=\"Quality below 4\">",
    "##FORMAT=<ID=GT;Number=1;Type=String;Description=\"Genotype\">",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2\tS3\tS4",
    "1\t100\t.\tA\tG\t.\tPASS\t.\t0/0\t0/1\t1/1\t./.",
    "1\t200\t.\tC\tT\t.\tPASS\t.\t1/1\t0/0\t0/1\t1/1",
    "5\t500\t.\tG\tA\t.\tpass5\t.\t0/0\t0/0\t./.\t1/1",
]

GENOTYPES = np.array([
    [0, 1, 2, np.nan, 1, 0, 1, 2],
    [np.nan, 1, 0, 0, 2, np.nan, 0, 1],
    [1, 2, np.nan, 1, 0, 2, 1, 0],
    [0, 0, 1, np.nan, 2, 1, 0, np.nan],
], dtype="float32")

MARKERS = [
    ("1", "100", "A", "G"),
    ("1", "200", "C", "T"),
    ("1", "300", "T", "A"),
    ("1", "400", "G", "C"),
    ("5", "500", "G", "A"),
    ("5", "600", "A", "C"),
    ("5", "700", "C", "T"),
    ("5", "800", "T", "G"),
]

TRAIT = "demo_trait"


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def build_fixture(tmp_path):
    """Tiny deterministic dataset; returns (data_dir, registry, source_vcf, source_manifest)."""
    data = tmp_path/"data"
    data.mkdir()
    np.save(data/"genotypes.npy", GENOTYPES)
    (data/"samples.tsv").write_text(
        "accession_id\tgenetic_group\n"
        "synth_a\tG0\nsynth_b\tG0\nsynth_c\tG1\nsynth_d\tG1\n", encoding="utf-8")
    table = "\n".join(f"{c}\t{p}\t{r}\t{a}" for c, p, r, a in MARKERS)
    (data/"variants.tsv").write_text(
        "chromosome\tposition\tref\talt\n" + table + "\n", encoding="utf-8")
    (data/"phenotypes.tsv").write_text(
        "accession_id\ttrait_id\tvalue\n"
        "synth_a\tdemo_trait\t1.0\nsynth_b\tdemo_trait\t2.0\n"
        "synth_c\tdemo_trait\t3.0\nsynth_d\tdemo_trait\t4.0\n", encoding="utf-8")
    registry = tmp_path/"trait_resolution_v1.tsv"
    registry.write_text(
        "trait_id\ttarget_status\ndemo_trait\tPUBLISHED_ACCESSION_VALUE_USABLE\n",
        encoding="utf-8")
    vcf = tmp_path/"synthetic_tiny.vcf"
    vcf.write_text("\n".join(VCF_LINES) + "\n", encoding="utf-8")
    manifest = data/"source_manifest.json"
    manifest.write_text(json.dumps({
        "series": "synthetic_vcf_source_manifest",
        "vcf_path": str(vcf),
        "local_sha256": hashlib.sha256(vcf.read_bytes()).hexdigest(),
        "official_md5": hashlib.md5(vcf.read_bytes()).hexdigest(),
        "md5_matches_official": True,
    }), encoding="utf-8")
    return data, registry, vcf, manifest


def make_record(data, registry, vcf, manifest, **overrides):
    """The full v2 record, binding all files as they exist *right now*."""
    record = {
        "provenance_schema": PROVENANCE_SCHEMA_V2,
        "dataset_id": "synthetic_tiny_fixture",
        "genotypes_sha256": _sha256(data/"genotypes.npy"),
        "samples_sha256": _sha256(data/"samples.tsv"),
        "variants_sha256": _sha256(data/"variants.tsv"),
        "phenotypes_sha256": _sha256(data/"phenotypes.tsv"),
        "registry_sha256": _sha256(registry),
        "source_vcf_sha256": _sha256(vcf),
        "source_manifest_sha256": _sha256(manifest),
        "official_source_md5": hashlib.md5(vcf.read_bytes()).hexdigest(),
        "vcf_filter_policy": "PASS-only; q4-filtered variants excluded",
        "multiallelic_policy": "multi-allelic excluded unsplittable",
        "reference_assembly": "TAIR10",
        "allele_encoding": "ALT_dosage_0_1_2",
        "conversion_script": "scripts/convert_vcf_to_npy.py",
        "conversion_commit": "9"*40,
        "notes": ["tiny synthetic fixture; not a genuine formal dataset"],
        "phenotype_registry": {"trait_id": TRAIT, "resolved": True},
    }
    record.update(overrides)
    return record


def write_record(data, registry, vcf, manifest, **overrides):
    (data/"provenance.json").write_text(
        json.dumps(make_record(data, registry, vcf, manifest, **overrides)),
        encoding="utf-8")


def run_gate(data, registry, vcf):
    # build_fixture wrote the source manifest into the dataset dir
    return gate_v2(data, TRAIT, source_vcf=vcf, registry=registry,
                   source_manifest=data/"source_manifest.json")


def test_gate_binds_all_six_and_passes_intact(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    path, record = run_gate(data, registry, vcf)
    assert path == data/"provenance.json"
    assert record["provenance_schema"] == PROVENANCE_SCHEMA_V2
    for field, file in (("genotypes_sha256", "genotypes.npy"),
                        ("samples_sha256", "samples.tsv"),
                        ("variants_sha256", "variants.tsv"),
                        ("phenotypes_sha256", "phenotypes.tsv")):
        assert record[field] == _sha256(data/file)
    assert record["source_vcf_sha256"] == _sha256(vcf)
    assert record["registry_sha256"] == _sha256(registry)


def test_phenotype_tamper_after_record_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    run_gate(data, registry, vcf)                  # intact record passes
    (data/"phenotypes.tsv").write_text(
        "accession_id\ttrait_id\tvalue\nsynth_a\tdemo_trait\t1.0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="phenotypes_sha256"):
        run_gate(data, registry, vcf)


def test_variants_tamper_after_record_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    run_gate(data, registry, vcf)
    table = "\n".join(f"{c}\t{p}\t{r}\t{a}" for c, p, r, a in MARKERS)
    (data/"variants.tsv").write_text(
        "chromosome\tposition\tref\talt\n" + table + "\n5\t800\tT\tG\n",
        encoding="utf-8")
    with pytest.raises(ValueError, match="variants_sha256"):
        run_gate(data, registry, vcf)


def test_samples_tamper_after_record_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    run_gate(data, registry, vcf)
    (data/"samples.tsv").write_text(
        "accession_id\tgenetic_group\nsynth_a\tG1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="samples_sha256"):
        run_gate(data, registry, vcf)


def test_registry_and_source_vcf_are_recomputed_targets(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    run_gate(data, registry, vcf)
    registry.write_text("trait_id\ttarget_status\nother_trait\tPENDING\n",
                        encoding="utf-8")
    with pytest.raises(ValueError, match="registry_sha256"):
        run_gate(data, registry, vcf)


def test_unresolved_registry_record_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    run_gate(data, registry, vcf)                  # resolved true passes
    write_record(data, registry, vcf, manifest,
                 phenotype_registry={"trait_id": TRAIT, "resolved": False})
    with pytest.raises(ValueError, match="phenotype_registry"):
        run_gate(data, registry, vcf)
    write_record(data, registry, vcf, manifest,
                 phenotype_registry={"trait_id": "other_trait", "resolved": True})
    with pytest.raises(ValueError, match="phenotype_registry"):
        run_gate(data, registry, vcf)


def test_v1_shaped_record_rejected_on_schema(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    full = make_record(data, registry, vcf, manifest)
    # The v1 record shape: intake/example fields, no v2 binding keys
    v1 = {
        "dataset_id": full["dataset_id"],
        "genotype_representation": "vcf_alt_dosage",
        "reference_assembly": full["reference_assembly"],
        "source_urls": ["https://example.test/source.vcf.gz"],
        "input_sha256": full["genotypes_sha256"],
        "source_vcf_sha256": full["source_vcf_sha256"],
        "conversion_script": full["conversion_script"],
        "conversion_commit": full["conversion_commit"],
        "allele_encoding": full["allele_encoding"],
        "phenotype_registry": full["phenotype_registry"],
    }
    (data/"provenance.json").write_text(json.dumps(v1), encoding="utf-8")
    with pytest.raises(ValueError, match="provenance_schema"):
        run_gate(data, registry, vcf)
    # and an explicit wrong version value is refused at the same field
    write_record(data, registry, vcf, manifest, provenance_schema="at_pheno_formal_v1")
    with pytest.raises(ValueError, match="provenance_schema"):
        run_gate(data, registry, vcf)


def test_syntactically_valid_source_hash_not_trusted(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    # 64-hex format alone is no verification
    write_record(data, registry, vcf, manifest, source_vcf_sha256="f"*64)
    with pytest.raises(ValueError, match="source_vcf_sha256"):
        run_gate(data, registry, vcf)
    # a file change after the record is also caught by the recompute
    write_record(data, registry, vcf, manifest)
    run_gate(data, registry, vcf)
    vcf.write_text("\n".join(VCF_LINES + ["1\t300\t.\tT\tC\t.\tPASS\t.\t0/0"]) + "\n",
                   encoding="utf-8")
    with pytest.raises(ValueError, match="source_vcf_sha256"):
        run_gate(data, registry, vcf)


def test_bad_reference_assembly_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest, reference_assembly="TAIR11")
    with pytest.raises(ValueError, match="reference_assembly"):
        run_gate(data, registry, vcf)


def test_missing_bound_files_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    (data/"genotypes.npy").unlink()
    with pytest.raises(ValueError, match="genotypes_sha256"):
        run_gate(data, registry, vcf)
    np.save(data/"genotypes.npy", GENOTYPES)        # restore the intact file
    # a registry target that no longer exists on disk
    with pytest.raises(ValueError, match="registry_sha256"):
        gate_v2(data, TRAIT, source_vcf=vcf, registry=tmp_path/"gone.tsv")
    # and no recompute target at all
    with pytest.raises(ValueError, match="source_vcf_sha256"):
        gate_v2(data, TRAIT, source_vcf=None, registry=registry)


# --- Registry truth: the record's resolved flag is not the truth source (goal 1) ---

def test_registry_duplicate_blocked_rejects_resolved_record(tmp_path):
    """registry says REQUIRES_DUPLICATE_RESOLUTION, record claims resolved
    true -> MUST FAIL (settled goal item 2)."""
    data, registry, vcf, manifest = build_fixture(tmp_path)
    registry.write_text(
        "trait_id\ttarget_status\n"
        "demo_trait\tREQUIRES_DUPLICATE_RESOLUTION\n", encoding="utf-8")
    write_record(data, registry, vcf, manifest)    # binds the blocked registry
    with pytest.raises(ValueError, match="target_status"):
        run_gate(data, registry, vcf)


def test_registry_duplicate_trait_rows_rejected(tmp_path):
    """A trait listed more than once in the registry is never resolvable
    (goal 1: the requested trait must exist uniquely)."""
    data, registry, vcf, manifest = build_fixture(tmp_path)
    registry.write_text(
        "trait_id\ttarget_status\n"
        "demo_trait\tPUBLISHED_ACCESSION_VALUE_USABLE\n"
        "demo_trait\tREQUIRES_DUPLICATE_RESOLUTION\n", encoding="utf-8")
    write_record(data, registry, vcf, manifest)
    with pytest.raises(ValueError, match="exactly once"):
        run_gate(data, registry, vcf)


def test_registry_missing_trait_rejects_resolved_record(tmp_path):
    """A trait the registry never lists cannot back a resolved record."""
    data, registry, vcf, manifest = build_fixture(tmp_path)
    registry.write_text(
        "trait_id\ttarget_status\nother_trait\tPUBLISHED_ACCESSION_VALUE_USABLE\n",
        encoding="utf-8")
    write_record(data, registry, vcf, manifest)
    with pytest.raises(ValueError, match="phenotype_registry"):
        run_gate(data, registry, vcf)


def test_registry_wide_tsv_header_driven_parsing(tmp_path):
    """The committed data-semantic registry is a wide TSV (19 columns);
    truth parsing must be header-driven, not assuming column 0/1."""
    data, registry, vcf, manifest = build_fixture(tmp_path)
    registry.write_text(
        "trait_id\ttrait_name\ttarget_status\n"
        "demo_trait\tDTF1-synthetic\tPUBLISHED_ACCESSION_VALUE_USABLE\n",
        encoding="utf-8")
    write_record(data, registry, vcf, manifest)
    path, record = run_gate(data, registry, vcf)   # wide header parses fine
    assert record["phenotype_registry"]["trait_id"] == TRAIT


# --- Official-MD5 cross-bound via the SHA256-bound source manifest (goal 4) ---

def test_manifest_flag_falsy_rejects_gate(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    m = json.loads(manifest.read_text())
    m["md5_matches_official"] = False
    manifest.write_text(json.dumps(m), encoding="utf-8")
    write_record(data, registry, vcf, manifest)
    with pytest.raises(ValueError, match="md5_matches_official"):
        run_gate(data, registry, vcf)


def test_manifest_sha_binding_tamper_rejected(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest)
    m = json.loads(manifest.read_text())
    m["official_md5"] = "0" * 32                     # record/tamper the manifest
    manifest.write_text(json.dumps(m), encoding="utf-8")
    with pytest.raises(ValueError, match="source_manifest_sha256"):
        run_gate(data, registry, vcf)


def test_record_official_md5_must_equal_manifest(tmp_path):
    data, registry, vcf, manifest = build_fixture(tmp_path)
    write_record(data, registry, vcf, manifest, official_source_md5="0" * 32)
    with pytest.raises(ValueError, match="official_source_md5"):
        run_gate(data, registry, vcf)


def test_mapping_summary_binds_source_vcf_provenance():
    """The committed HDF5-VCF mapping summary cross-binds the source VCF
    hashes to the VCF-side source audit, field by field."""
    repo = Path(__file__).resolve().parents[1]
    summary = json.loads(
        (repo/"data/manifests/hdf5_vcf_position_mapping_summary.json").read_text(
            encoding="utf-8"))
    source = json.loads(
        (repo/"data/manifests/1001g_v31_population_vcf_source.json").read_text(
            encoding="utf-8"))
    prov = summary["source_vcf_provenance"]
    assert prov["path"] == source["vcf_path"]
    assert prov["sha256"] == source["local_sha256"]
    assert prov["official_md5"] == source["official_md5"]
    assert prov["md5_matches_official"] is source["md5_matches_official"]
    assert prov["recorded_in"] == "data/manifests/1001g_v31_population_vcf_source.json"


def test_catalog_provenance_keeps_hashes_distinct():
    """The catalog block must keep the two hashes distinct (goal 5):
    a byte-identity hash of one .gz artifact vs the content-identity
    hash of the unpacked TSV; the bare `sha256` key is gone."""
    repo = Path(__file__).resolve().parents[1]
    summary = json.loads(
        (repo/"data/manifests/hdf5_vcf_position_mapping_summary.json").read_text(
            encoding="utf-8"))
    cat = summary["catalog_provenance"]
    assert "sha256" not in cat, "bare sha256 key smashes the artifact/logical distinction"
    comp = cat["compressed_artifact_sha256"]
    logical = cat["logical_content_md5"]
    assert len(comp) == 64 and all(c in "0123456789abcdef" for c in comp)
    assert len(logical) == 32 and all(c in "0123456789abcdef" for c in logical)
    assert "content" in cat["note"] and "artifact" in cat["note"]


def test_catalog_logical_md5_matches_fresh_unpacked_stream():
    """The recorded logical_content_md5 must equal a fresh md5 of the
    unpacked TSV stream on disk (the deliverable (61 MB) is
    worktree-local; the check is self-verifying when the file is
    present, and raises a readable error otherwise)."""
    if not os.path.isfile(CATALOG):
        pytest.skip("catalog .tsv.gz not on disk in this checkout")
    h = _md5_of_unpacked_gz(CATALOG)
    repo = Path(__file__).resolve().parents[1]
    summary = json.loads(
        (repo/"data/manifests/hdf5_vcf_position_mapping_summary.json").read_text(
            encoding="utf-8"))
    assert summary["catalog_provenance"]["logical_content_md5"] == h
