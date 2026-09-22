"""Genotype-only COMMON feasibility audit across frozen G0.1 train partitions.

Run this file with a CUDA-enabled Python that has numpy and torch.  It reads
only genotype storage, frozen sample IDs and split manifests; phenotype files
and phenotype values are neither opened nor parsed.
"""

import argparse
import csv
import json
from pathlib import Path

import numpy as np


def rows(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def frozen_partitions(split_dir, trait_rows, sample_index):
    """Return only training membership masks derived from frozen split files."""
    result = []
    for trait in trait_rows:
        label = trait["trait_name"].replace(" ", "_")
        for protocol in ("blocked", "iid"):
            outer = rows(split_dir / f"{label}_outer_{protocol}.tsv")
            by_fold = {int(row["fold"]): row for row in outer}
            for held in sorted(by_fold):
                train = [sample_index[row["accession_id"]] for row in outer if int(row["fold"]) != held]
                result.append({"trait_id": trait["trait_id"], "trait_name": trait["trait_name"], "protocol": protocol,
                               "level": "outer", "held_out": str(held), "train_indices": train})
                inner = rows(split_dir / f"{label}_inner_{protocol}_outer{held}.tsv")
                for inner_held in sorted({int(row["fold"]) for row in inner}):
                    train_inner = [sample_index[row["accession_id"]] for row in inner if int(row["fold"]) != inner_held]
                    result.append({"trait_id": trait["trait_id"], "trait_name": trait["trait_name"], "protocol": protocol,
                                   "level": "inner", "held_out": f"outer{held}_inner{inner_held}", "train_indices": train_inner})
    logo = rows(split_dir / "group_logo_manifest.tsv")
    by_logo = {}
    for row in logo:
        by_logo.setdefault((row["trait_id"], row["trait_name"], row["held_out_group"]), []).append(row)
    for (trait_id, trait_name, held), records in sorted(by_logo.items()):
        train = [sample_index[row["accession_id"]] for row in records if row["is_test"] == "False"]
        result.append({"trait_id": trait_id, "trait_name": trait_name, "protocol": "logo", "level": "outer",
                       "held_out": held, "train_indices": train})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genotypes", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--trait-panel", type=Path, required=True)
    parser.add_argument("--splits", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--partition-batch", type=int, default=32)
    parser.add_argument("--marker-block", type=int, default=8192)
    args = parser.parse_args()
    import torch
    if args.device.startswith("cuda") and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested but unavailable")
    samples = rows(args.samples)
    ids = [row["accession_id"] for row in samples]
    sample_index = {value: index for index, value in enumerate(ids)}
    if len(ids) != len(sample_index):
        raise ValueError("sample IDs are not unique")
    trait_rows = [row for row in rows(args.trait_panel) if row["role"] in {"PRIMARY", "INDEPENDENT_REPLICATION", "SECONDARY_SENSITIVITY"}]
    partitions = frozen_partitions(args.splits, trait_rows, sample_index)
    x = np.load(args.genotypes, mmap_mode="r", allow_pickle=False)
    if x.dtype != np.int8 or x.shape[0] != len(ids):
        raise ValueError("formal int8 genotype matrix/sample mismatch")
    weights = np.zeros((len(partitions), len(ids)), dtype=np.float32)
    train_n = np.zeros(len(partitions), dtype=np.float32)
    for index, partition in enumerate(partitions):
        weights[index, partition["train_indices"]] = 1
        train_n[index] = len(partition["train_indices"])
        if train_n[index] < 2:
            raise ValueError(f"training partition too small: {partition}")
    call_count = np.zeros(len(partitions), dtype=np.int64)
    maf_count = np.zeros(len(partitions), dtype=np.int64)
    common_count = np.zeros(len(partitions), dtype=np.int64)
    device = torch.device(args.device)
    for marker_start in range(0, x.shape[1], args.marker_block):
        # np.memmap slices are read-only.  Copy before torch conversion so the
        # tensor contract cannot accidentally expose writable storage.
        source = np.array(x[:, marker_start:marker_start + args.marker_block], copy=True)
        if np.any((source < -1) | (source > 2)):
            raise ValueError("formal genotype contract violated")
        genotype = torch.from_numpy(source).to(device)
        called = (genotype >= 0).to(torch.float32)
        dosage = torch.clamp_min(genotype, 0).to(torch.float32)
        for start in range(0, len(partitions), args.partition_batch):
            stop = min(start + args.partition_batch, len(partitions))
            w = torch.from_numpy(weights[start:stop]).to(device)
            calls = w @ called
            total = w @ dosage
            rate_ok = calls >= torch.from_numpy(.95 * train_n[start:stop, None]).to(device)
            p = total / torch.clamp_min(2 * calls, 1)
            maf_ok = (calls > 0) & (torch.minimum(p, 1 - p) >= .05) & (total > 0) & (total < 2 * calls)
            common = rate_ok & maf_ok
            call_count[start:stop] += rate_ok.sum(dim=1).cpu().numpy().astype(np.int64)
            maf_count[start:stop] += maf_ok.sum(dim=1).cpu().numpy().astype(np.int64)
            common_count[start:stop] += common.sum(dim=1).cpu().numpy().astype(np.int64)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["trait_id", "trait_name", "protocol", "level", "held_out", "train_n", "n_callrate95", "n_MAF05_nonmonomorphic", "n_COMMON", "COMMON_10K_feasible", "COMMON_50K_feasible", "COMMON_250K_feasible", "COMMON_1M_feasible"]
    output = []
    for partition, call, maf, common, n in zip(partitions, call_count, maf_count, common_count, train_n.astype(int)):
        output.append({**{key: partition[key] for key in ("trait_id", "trait_name", "protocol", "level", "held_out")}, "train_n": int(n),
                       "n_callrate95": int(call), "n_MAF05_nonmonomorphic": int(maf), "n_COMMON": int(common),
                       **{f"COMMON_{label}_feasible": common >= threshold for label, threshold in (("10K", 10_000), ("50K", 50_000), ("250K", 250_000), ("1M", 1_000_000))}})
    with args.out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    primary_ids = {row["trait_id"] for row in trait_rows if row["role"] == "PRIMARY"}
    primary = [row for row in output if row["trait_id"] in primary_ids and row["protocol"] == "blocked"]
    summary = {"schema": "common_feasibility_v1", "status": "PASS" if all(row["COMMON_250K_feasible"] for row in primary) else "PRIMARY_250K_BLOCKER",
               "implementation": "CUDA blockwise genotype-only matrix multiplication", "phenotype_values_read": False,
               "n_partitions": len(output), "primary_trait_ids": sorted(primary_ids), "primary_blocked_partition_n": len(primary),
               "primary_all_250K_feasible": all(row["COMMON_250K_feasible"] for row in primary),
               "minimum_COMMON_by_protocol": {protocol: min(row["n_COMMON"] for row in output if row["protocol"] == protocol) for protocol in sorted({row["protocol"] for row in output})}}
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    if summary["status"] != "PASS":
        raise SystemExit(3)


if __name__ == "__main__":
    main()
