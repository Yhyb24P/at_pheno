"""Exact-refine near-clone edges nominated by the registered 250K d<=0.002 screen.

This is a correction for finite-panel threshold noise: it preserves the
coordinate-hash panel and the scientific d<=0.001 primary cutoff, but uses
the broader pre-registered d<=0.002 sensitivity graph solely as a conservative
candidate screen.  Candidate pairs are then scored on every HDF5 binary state.
"""

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path

import h5py
import numpy as np


def load_nearclone_module():
    path = Path(__file__).with_name("freeze_nearclone_blocks.py")
    spec = importlib.util.spec_from_file_location("nearclone_freeze", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bernoulli_kl(q, p):
    return q * math.log(q / p) + (1 - q) * math.log((1 - q) / (1 - p))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    module = load_nearclone_module()
    with args.samples.open(newline="", encoding="utf-8") as stream:
        ids = [row["accession_id"] for row in csv.DictReader(stream, delimiter="\t")]
    labels = {}
    with args.samples.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            labels[row["accession_id"]] = row["genetic_group"]
    panel = list(csv.DictReader(args.panel.open(newline="", encoding="utf-8"), delimiter="\t"))
    indices = [int(row["hdf5_row"]) for row in panel]
    if len(indices) != 250_000 or indices != sorted(indices) or len(indices) != len(set(indices)):
        raise ValueError("nearclone panel must contain 250K unique physical-order HDF5 rows")
    args.out.mkdir(parents=True)
    with h5py.File(args.hdf5, "r") as handle:
        h5_ids = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["accessions"][:]]
        if h5_ids != ids:
            raise ValueError("HDF5 sample order differs from formal samples")
        estimated = module.hamming_counts(handle, indices)
        screen_limit = int(.002 * len(indices))
        left, right = np.where(np.triu(estimated <= screen_limit, 1))
        if not len(left):
            raise ValueError("0.002 screen has no candidate pairs")
        exact = np.zeros(len(left), dtype=np.int64)
        n_sites = handle["snps"].shape[0]
        for start in range(0, n_sites, 10_000):
            states = handle["snps"][start:start + 10_000, :]
            if np.any((states != 0) & (states != 1)):
                raise ValueError("HDF5 binary-state contract violated")
            exact += np.count_nonzero(states[:, left] != states[:, right], axis=0)
    primary_limit = int(.001 * n_sites)
    exact_edge = exact <= primary_limit
    with (args.out / "nearclone_screen_d002_exact_pairs.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["accession_a", "accession_b", "screen_hamming_distance", "exact_hamming_distance", "exact_primary_edge"], delimiter="\t")
        writer.writeheader()
        for a, b, estimate, observed, edge in zip(left, right, estimated[left, right], exact, exact_edge):
            writer.writerow({"accession_a": ids[a], "accession_b": ids[b], "screen_hamming_distance": estimate / len(indices),
                             "exact_hamming_distance": observed / n_sites, "exact_primary_edge": bool(edge)})
    matrix = np.full((len(ids), len(ids)), primary_limit + 1, dtype=np.int64)
    np.fill_diagonal(matrix, 0)
    matrix[left, right] = exact
    matrix[right, left] = exact
    # Reuse the diagnostic implementation with its denominator explicitly
    # switched from the screening panel to the exact HDF5 universe.
    module.PANEL_SIZE = n_sites
    groups = module.components(matrix, primary_limit)
    diagnostics = module.diagnostics(matrix, groups, primary_limit, ids, labels)
    with (args.out / "nearclone_blocks_exact_d001.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["accession_id", "nearclone_block"], delimiter="\t")
        writer.writeheader()
        for block_id, members in enumerate(groups):
            writer.writerows({"accession_id": ids[index], "nearclone_block": block_id} for index in members)
    with (args.out / "nearclone_exact_component_diagnostics.tsv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(diagnostics)
    # Under independent hash sampling, KL Chernoff bound for a true d=0.001
    # pair exceeding the d=0.002 screen.  Reported as a model bound, not an
    # empirical probability, and family-adjusted across all accession pairs.
    single_bound = math.exp(-len(indices) * bernoulli_kl(.002, .001))
    pair_n = len(ids) * (len(ids) - 1) // 2
    sizes = [len(group) for group in groups]
    summary = {"schema": "nearclone_exact_refinement_v1", "status": "EXACT_PRIMARY_EDGE_REFINEMENT_COMPLETE",
               "primary_cutoff": .001, "primary_max_hamming_count": primary_limit,
               "candidate_screen": {"panel_size": len(indices), "screen_cutoff": .002, "candidate_pairs": int(len(left)),
                                    "hash_sampling_model": "Bernoulli approximation", "single_pair_false_negative_upper_bound": single_bound,
                                    "familywise_union_bound": min(1.0, pair_n * single_bound)},
               "exact_primary": {"n_edges": int(exact_edge.sum()), "n_components": len(groups),
                                 "multi_member_components": sum(size > 1 for size in sizes), "largest_component": max(sizes),
                                 "chain_flag": diagnostics[0]["size"] > 1 and diagnostics[0]["edge_density"] < .25 and diagnostics[0]["max_within_hamming_distance"] > .001,
                                 "status": "PASS" if max(sizes) <= .25 * len(ids) and not (diagnostics[0]["size"] > 1 and diagnostics[0]["edge_density"] < .25 and diagnostics[0]["max_within_hamming_distance"] > .001) else "NEARCLONE_CHAIN_BLOCKER"},
               "noncanonical_note": "Exact refinement corrects screening noise; canonical adoption requires design-manifest update."}
    (args.out / "nearclone_exact_refinement_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["exact_primary"]))


if __name__ == "__main__":
    main()
