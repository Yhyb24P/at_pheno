"""TASK04: formal provenance v2 gate on a fully synthetic tiny fixture.

Nothing under data/raw and no real (12 GB-class) release file is
read or hashed.  The fixture under tmp_path holds a static 12-line
synthetic source VCF, a fixed 4x8 genotypes matrix (0/1/2 plus NaN),
three small TSVs, and a hand-baked v2 provenance record whose
bindings are recomputed at test time.  A record written *before*
tampering must be rejected by a ValueError naming the breached
field;  intact, the gate passes.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from at_pheno.formal_provenance import PROVENANCE_SCHEMA_V2, gate_v2

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
    """Tiny deterministic dataset; returns (data_dir, registry_tsv, source_vcf)."""
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
    return data, registry, vcf


def make_record(data, registry, vcf, **overrides):
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


def write_record(data, registry, vcf, **overrides):
    (data/"provenance.json").write_text(
        json.dumps(make_record(data, registry, vcf, **overrides)), encoding="utf-8")


def run_gate(data, registry, vcf):
    return gate_v2(data, TRAIT, source_vcf=vcf, registry=registry)


def test_gate_binds_all_six_and_passes_intact(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
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
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
    run_gate(data, registry, vcf)                  # intact record passes
    (data/"phenotypes.tsv").write_text(
        "accession_id\ttrait_id\tvalue\nsynth_a\tdemo_trait\t1.0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="phenotypes_sha256"):
        run_gate(data, registry, vcf)


def test_variants_tamper_after_record_rejected(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
    run_gate(data, registry, vcf)
    table = "\n".join(f"{c}\t{p}\t{r}\t{a}" for c, p, r, a in MARKERS)
    (data/"variants.tsv").write_text(
        "chromosome\tposition\tref\talt\n" + table + "\n5\t800\tT\tG\n",
        encoding="utf-8")
    with pytest.raises(ValueError, match="variants_sha256"):
        run_gate(data, registry, vcf)


def test_samples_tamper_after_record_rejected(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
    run_gate(data, registry, vcf)
    (data/"samples.tsv").write_text(
        "accession_id\tgenetic_group\nsynth_a\tG1\n", encoding="utf-8")
    with pytest.raises(ValueError, match="samples_sha256"):
        run_gate(data, registry, vcf)


def test_registry_and_source_vcf_are_recomputed_targets(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
    run_gate(data, registry, vcf)
    registry.write_text("trait_id\ttarget_status\nother_trait\tPENDING\n",
                        encoding="utf-8")
    with pytest.raises(ValueError, match="registry_sha256"):
        run_gate(data, registry, vcf)


def test_unresolved_registry_record_rejected(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
    run_gate(data, registry, vcf)                  # resolved true passes
    write_record(data, registry, vcf,
                 phenotype_registry={"trait_id": TRAIT, "resolved": False})
    with pytest.raises(ValueError, match="phenotype_registry"):
        run_gate(data, registry, vcf)
    write_record(data, registry, vcf,
                 phenotype_registry={"trait_id": "other_trait", "resolved": True})
    with pytest.raises(ValueError, match="phenotype_registry"):
        run_gate(data, registry, vcf)


def test_v1_shaped_record_rejected_on_schema(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    full = make_record(data, registry, vcf)
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
    write_record(data, registry, vcf, provenance_schema="at_pheno_formal_v1")
    with pytest.raises(ValueError, match="provenance_schema"):
        run_gate(data, registry, vcf)


def test_syntactically_valid_source_hash_not_trusted(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    # 64-hex format alone is no verification
    write_record(data, registry, vcf, source_vcf_sha256="f"*64)
    with pytest.raises(ValueError, match="source_vcf_sha256"):
        run_gate(data, registry, vcf)
    # a file change after the record is also caught by the recompute
    write_record(data, registry, vcf)
    run_gate(data, registry, vcf)
    vcf.write_text("\n".join(VCF_LINES + ["1\t300\t.\tT\tC\t.\tPASS\t.\t0/0"]) + "\n",
                   encoding="utf-8")
    with pytest.raises(ValueError, match="source_vcf_sha256"):
        run_gate(data, registry, vcf)


def test_bad_reference_assembly_rejected(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf, reference_assembly="TAIR11")
    with pytest.raises(ValueError, match="reference_assembly"):
        run_gate(data, registry, vcf)


def test_missing_bound_files_rejected(tmp_path):
    data, registry, vcf = build_fixture(tmp_path)
    write_record(data, registry, vcf)
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
