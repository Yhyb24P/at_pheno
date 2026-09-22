import json
import sys
from argparse import Namespace

import numpy as np
import pytest

from at_pheno.cli import demo, load_dataset, main, read_table, run, sha256, write_json, write_table
from at_pheno.core import random_folds


def _convert_demo_to_formal_int8(data):
    """Convert a synthetic pilot fixture into the formal storage contract."""
    values = np.load(data/"genotypes.npy")
    formal = np.where(np.isnan(values), -1, values).astype(np.int8)
    np.save(data/"genotypes.npy", formal)


def test_outer_test_labels_cannot_change_own_predictions_or_tuning(tmp_path):
    data = tmp_path/"data"
    demo(data)
    config = tmp_path/"test.toml"
    config.write_text('seed=7\nouter_folds=3\ninner_folds=2\ndensities=[16,"all"]\nalphas=[0.1,1.0]\nmin_call_rate=0.9\nmin_maf=0.05\nmin_mac=0\nblock_size=32\nhash_salt="test"\n')
    args = Namespace(data=data, trait="synthetic", config=config, out=tmp_path/"a", protocol="iid", splits=None)
    run(args)
    first = read_table(args.out/"predictions.tsv", ["accession_id", "prediction"])
    ids = [r["accession_id"] for r in read_table(data/"samples.tsv", ["accession_id"])]
    heldout = {s for s, f in zip(ids, random_folds(ids, 3, 7)) if f == 0}
    phenos = read_table(data/"phenotypes.tsv", ["accession_id", "value"])
    for row in phenos:
        if row["accession_id"] in heldout:
            row["value"] = str(float(row["value"])+10000)
    write_table(data/"phenotypes.tsv", phenos, list(phenos[0]))
    args.out = tmp_path/"b"
    run(args)
    second = read_table(args.out/"predictions.tsv", ["accession_id", "prediction"])
    a = [(r["accession_id"], r["density"], r["prediction"], r["alpha"]) for r in first if r["fold"] == "0"]
    b = [(r["accession_id"], r["density"], r["prediction"], r["alpha"]) for r in second if r["fold"] == "0"]
    assert a == b
    for density in {r["density"] for r in first}:
        assert len({r["accession_id"] for r in first if r["density"] == density}) == 60


def test_replicate_rows_are_not_silently_averaged(tmp_path):
    data = tmp_path/"data"
    demo(data)
    rows = read_table(data/"phenotypes.tsv", ["accession_id", "value"])
    write_table(data/"phenotypes.tsv", rows+[rows[0]], list(rows[0]))
    with pytest.raises(ValueError, match="Repeated phenotype"):
        load_dataset(data, "synthetic")


def test_unsorted_physical_marker_manifest_is_rejected(tmp_path):
    data = tmp_path/"data"
    demo(data)
    rows = read_table(data/"variants.tsv", ["chromosome", "position"])
    rows[0], rows[1] = rows[1], rows[0]
    write_table(data/"variants.tsv", rows, list(rows[0]))
    with pytest.raises(ValueError, match="strictly ordered"):
        load_dataset(data, "synthetic")


def test_group_outer_and_inner_splits_never_mix_declared_groups(tmp_path):
    data = tmp_path/"data"
    demo(data)
    config = tmp_path/"test.toml"
    config.write_text('seed=7\nouter_folds=3\ninner_folds=2\ndensities=[16]\nalphas=[0.1,1.0]\nmin_call_rate=0.9\nmin_maf=0.05\nmin_mac=0\nblock_size=32\nhash_salt="test"\n')
    out = tmp_path/"group"
    run(Namespace(data=data, trait="synthetic", config=config, out=out, protocol="group", splits=None, mode="pilot"))
    outer = read_table(out/"splits.tsv", ["fold", "genetic_group"])
    assert all(len({r["fold"] for r in outer if r["genetic_group"] == group}) == 1
               for group in {r["genetic_group"] for r in outer})
    for path in out.glob("inner_splits_*.tsv"):
        rows = read_table(path, ["accession_id", "fold"])
        groups = {r["accession_id"]: r["genetic_group"] for r in outer}
        assert all(len({r["fold"] for r in rows if groups[r["accession_id"]] == group}) == 1
                   for group in {groups[r["accession_id"]] for r in rows})


def test_formal_mode_requires_alt_dosage_provenance(tmp_path):
    data = tmp_path/"data"
    demo(data)
    config = tmp_path/"test.toml"
    config.write_text('seed=7\nouter_folds=3\ninner_folds=2\ndensities=[16]\nalphas=[0.1]\nmin_call_rate=0.9\nmin_maf=0.05\nmin_mac=0\nblock_size=32\nhash_salt="test"\n')
    with pytest.raises(ValueError, match="Formal provenance"):
        run(Namespace(data=data, trait="synthetic", config=config, out=tmp_path/"formal", protocol="iid", splits=None, mode="formal"))


def _v2_record(data, source, registry, manifest, genotypes_sha256="0"*64):
    return {
        "provenance_schema": "at_pheno_formal_v2", "dataset_id": "test",
        "genotypes_sha256": genotypes_sha256,
        "samples_sha256": sha256(data/"samples.tsv"),
        "variants_sha256": sha256(data/"variants.tsv"),
        "phenotypes_sha256": sha256(data/"phenotypes.tsv"),
        "registry_sha256": sha256(registry),
        "source_vcf_sha256": sha256(source),
        "source_manifest_sha256": sha256(manifest),
        "official_source_md5": "0"*32,
        "vcf_filter_policy": "PASS-only nuclear biallelic SNPs",
        "multiallelic_policy": "multi-allelic excluded unsplittable",
        "reference_assembly": "TAIR10", "allele_encoding": "ALT_dosage_0_1_2",
        "conversion_script": "convert.py", "conversion_commit": "a"*40,
        "notes": ["tiny synthetic pilot record"],
        "phenotype_registry": {"trait_id": "synthetic", "resolved": True}}


def test_formal_mode_binds_actual_genotype_hash_and_resolved_trait(tmp_path):
    data = tmp_path/"data"
    demo(data)
    _convert_demo_to_formal_int8(data)
    source = tmp_path/"source.vcf"
    source.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\n1\t100\t.\tA\tG\t.\tPASS\tGT\t0/0\n")
    registry = tmp_path/"trait_resolution_v1.tsv"
    registry.write_text("trait_id\ttarget_status\nsynthetic\tPUBLISHED_ACCESSION_VALUE_USABLE\n")
    manifest = tmp_path/"source_manifest.json"
    write_json(manifest, {"vcf_path": str(source),
                          "local_sha256": sha256(source),
                          "official_md5": "0" * 32,
                          "md5_matches_official": True})
    # v1 record (no provenance_schema key) is no longer accepted by the formal gate
    v1 = {"dataset_id": "test", "genotype_representation": "vcf_alt_dosage",
          "reference_assembly": "TAIR10", "source_urls": ["https://example.test/source.vcf.gz"],
          "input_sha256": "0" * 64, "source_vcf_sha256": "1" * 64,
          "conversion_script": "convert.py", "conversion_commit": "a" * 40,
          "allele_encoding": "ALT_dosage_0_1_2",
          "phenotype_registry": {"trait_id": "synthetic", "resolved": True}}
    write_json(data/"provenance.json", v1)
    config = tmp_path/"test.toml"
    config.write_text('seed=7\nouter_folds=3\ninner_folds=2\ndensities=[16]\nalphas=[0.1]\nmin_call_rate=0.9\nmin_maf=0.05\nmin_mac=0\nblock_size=32\nhash_salt="test"\n')
    args = Namespace(data=data, trait="synthetic", config=config, out=tmp_path/"formal", protocol="iid",
                     splits=None, mode="formal", source_vcf=source, registry=registry,
                     source_manifest=manifest)
    with pytest.raises(ValueError, match="provenance_schema"):
        run(args)
    # v2 record with a wrong genotypes binding -> recompute-not-trust rejection
    write_json(data/"provenance.json", _v2_record(data, source, registry, manifest))
    with pytest.raises(ValueError, match="does not bind"):
        run(args)
    write_json(data/"provenance.json", _v2_record(data, source, registry, manifest,
                                                   genotypes_sha256=sha256(data/"genotypes.npy")))
    run(args)
    assert (args.out/"provenance.json").is_file()


def test_formal_cli_fixture_argparse_run_gate_v2(tmp_path, monkeypatch):
    """Real argparse -> run -> gate_v2 wiring through the CLI itself
    (goal 1): every flag, including the new --source-manifest, arrives
    via the real parser, not a hand-built Namespace."""
    data = tmp_path/"data"
    demo(data)
    _convert_demo_to_formal_int8(data)
    config = tmp_path/"cli.toml"
    config.write_text('seed=7\nouter_folds=3\ninner_folds=2\ndensities=[16]\nalphas=[0.1]\nmin_call_rate=0.9\nmin_maf=0.05\nmin_mac=0\nblock_size=32\nhash_salt="cli"\n')
    source = tmp_path/"source.vcf"
    source.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\n1\t100\t.\tA\tG\t.\tPASS\tGT\t0/0\n")
    registry = tmp_path/"trait_resolution_v1.tsv"
    registry.write_text("trait_id\ttarget_status\nsynthetic\tPUBLISHED_ACCESSION_VALUE_USABLE\n")
    manifest = tmp_path/"source_manifest.json"
    write_json(manifest, {"series": "cli_fixture_source_manifest", "vcf_path": str(source),
                          "local_sha256": sha256(source),
                          "official_md5": "0" * 32,
                          "md5_matches_official": True})
    write_json(data/"provenance.json",
               _v2_record(data, source, registry, manifest,
                          genotypes_sha256=sha256(data/"genotypes.npy")))
    out = tmp_path/"cli_formal_out"
    monkeypatch.setattr(sys, "argv", [
        "at_pheno", "run",
        "--data", str(data), "--trait", "synthetic",
        "--config", str(config), "--out", str(out),
        "--protocol", "iid", "--mode", "formal",
        "--source-vcf", str(source),
        "--registry", str(registry), "--source-manifest", str(manifest)])
    main()
    assert (out/"provenance.json").is_file()
    prov = json.loads((out/"provenance.json").read_text(encoding="utf-8"))
    assert prov["mode"] == "formal"
    assert prov["status"] == "formal_data_contract_passed_not_confirmatory"
    # The run log binds the re-verified registry TSV and the SHA256-bound
    # source manifest; the raw source VCF is deliberately NOT digested
    # by the experiment run (its identity is a build-stage check).
    assert prov["input_sha256"][str(registry)] == sha256(registry)
    assert prov["input_sha256"][str(manifest)] == sha256(manifest)
    assert str(source) not in prov["input_sha256"]
    # The v2 record travels through the run, still cross-bound to the
    # gate-re-verified SHA256-bound source manifest record:
    assert prov["formal_dataset_provenance"]["source_manifest_sha256"] == sha256(manifest)
    assert prov["formal_dataset_provenance"]["provenance_schema"] == "at_pheno_formal_v2"
