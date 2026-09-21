from argparse import Namespace
import pytest

from at_pheno.cli import demo, load_dataset, read_table, run, write_table
from at_pheno.core import random_folds


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
