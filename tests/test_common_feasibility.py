import importlib.util
import csv
import json
from pathlib import Path


def module():
    path = Path(__file__).parents[1] / "scripts" / "audit_common_feasibility_gpu.py"
    spec = importlib.util.spec_from_file_location("common_feasibility", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def write_tsv(path, header, records):
    path.write_text("\t".join(header) + "\n" + "\n".join("\t".join(row) for row in records) + "\n")


def test_frozen_partitions_uses_split_membership_only(tmp_path):
    mod = module()
    split = tmp_path / "splits"
    split.mkdir()
    trait_rows = [{"trait_id": "1", "trait_name": "Trait", "role": "PRIMARY"}]
    header = ["trait_id", "trait_name", "accession_id", "fold"]
    outer = [["1", "Trait", "a", "0"], ["1", "Trait", "b", "0"], ["1", "Trait", "c", "1"], ["1", "Trait", "d", "1"]]
    write_tsv(split / "Trait_outer_blocked.tsv", header, outer)
    write_tsv(split / "Trait_outer_iid.tsv", header, outer)
    for protocol in ("blocked", "iid"):
        for fold in (0, 1):
            write_tsv(split / f"Trait_inner_{protocol}_outer{fold}.tsv", header, outer)
    write_tsv(split / "group_logo_manifest.tsv", ["trait_id", "trait_name", "accession_id", "held_out_group", "logo_fold", "is_test"],
              [["1", "Trait", "a", "g1", "g1", "False"], ["1", "Trait", "b", "g1", "g1", "False"],
               ["1", "Trait", "c", "g1", "g1", "True"], ["1", "Trait", "d", "g1", "g1", "True"]])
    partitions = mod.frozen_partitions(split, trait_rows, {name: i for i, name in enumerate("abcd")})
    assert len(partitions) == 13  # 2 protocols * (2 outer + 2*2 inner) + 1 LOGO
    assert all(set(item["train_indices"]) <= {0, 1, 2, 3} for item in partitions)
    logo = [item for item in partitions if item["protocol"] == "logo"]
    assert logo[0]["train_indices"] == [0, 1]


def test_g0_1_primary_feasibility_and_baseline_are_frozen():
    root = Path(__file__).parents[1]
    spec = json.loads((root / "data/manifests/confirmatory_v1_1/baseline_spec_v1.json").read_text())
    assert spec["kernel"]["formula"] == "K = ZZ^T / sum_j[2 p_j (1-p_j)]"
    assert spec["hyperparameter_selection"]["alpha_hard_log10_cap"] == [-12, 12]
    with (root / "data/manifests/confirmatory_v1_1/common_feasibility.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    primary = [row for row in rows if row["trait_id"] == "703" and row["protocol"] == "blocked"]
    assert len(rows) == 237 and len(primary) == 25
    assert all(row["COMMON_250K_feasible"] == "True" for row in primary)
    amendment = json.loads((root / "data/manifests/confirmatory_v1_1/analysis_plan_v1_1.json").read_text())
    assert amendment["density_availability"]["COMMON_1M"].startswith("PREDECLARED_UNAVAILABLE")
