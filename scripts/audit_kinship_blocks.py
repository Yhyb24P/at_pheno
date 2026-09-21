"""Audit, but do not freeze, candidate near-kinship blocks from official IBS.

The thresholds are sensitivity analyses.  A threshold becomes a formal split
definition only after it is selected in the registered analysis plan.
"""

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


def components(matrix: np.ndarray, threshold: float) -> list[list[int]]:
    n = matrix.shape[0]
    seen = np.zeros(n, dtype=bool)
    result = []
    for start in range(n):
        if seen[start]:
            continue
        stack, group = [start], []
        seen[start] = True
        while stack:
            node = stack.pop()
            group.append(node)
            neighbours = np.flatnonzero((matrix[node] >= threshold) & ~seen)
            seen[neighbours] = True
            stack.extend(neighbours.tolist())
        result.append(sorted(group))
    return result


def component_diagnostics(matrix: np.ndarray, groups: list[list[int]], threshold: float) -> list[dict]:
    records = []
    for block_id, members in enumerate(groups):
        size = len(members)
        if size == 1:
            minimum, edge_density = None, 0.0
        else:
            values = matrix[np.ix_(members, members)]
            lower = values[np.tril_indices(size, k=-1)]
            minimum = float(np.min(lower))
            edge_density = float(np.mean(lower >= threshold))
        records.append({"candidate_block": block_id, "size": size,
                        "min_pairwise_similarity": minimum, "edge_density": edge_density})
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kinship", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--thresholds", type=float, nargs="+", default=[6.5, 6.8, 7.0, 7.2],
                        help="Similarity thresholds on this file's native scale; calibrate before formal use")
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    with h5py.File(args.kinship, "r") as handle:
        ids = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["accessions"][:]]
        kinship = handle["kinship"][:]
        n_snps = int(handle["n_snps"][()])
    if kinship.shape != (len(ids), len(ids)) or not np.allclose(kinship, kinship.T):
        raise ValueError("Invalid non-symmetric kinship matrix")
    off_diagonal = kinship[np.triu_indices(len(ids), k=1)]
    args.out.mkdir(parents=True)
    summary = {"n_accessions": len(ids), "n_snps": n_snps,
               "off_diagonal_quantiles": {str(q): float(np.quantile(off_diagonal, q))
                                         for q in [0, .5, .9, .95, .99, .999, 1]},
               "thresholds": {}}
    for threshold in args.thresholds:
        groups = components(kinship, threshold)
        groups.sort(key=lambda values: (-len(values), values))
        with (args.out / f"components_ge_{threshold:g}.tsv").open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t")
            writer.writerow(["candidate_block", "accession_id"])
            for block_id, members in enumerate(groups):
                writer.writerows((block_id, ids[index]) for index in members)
        diagnostics = component_diagnostics(kinship, groups, threshold)
        with (args.out / f"components_ge_{threshold:g}_diagnostics.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(diagnostics[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(diagnostics)
        sizes = [len(group) for group in groups]
        summary["thresholds"][str(threshold)] = {"n_components": len(groups),
            "n_multimember_components": sum(size > 1 for size in sizes),
            "largest_component": max(sizes), "accessions_in_multimember_components": sum(size for size in sizes if size > 1),
            "largest_component_min_pairwise_similarity": diagnostics[0]["min_pairwise_similarity"],
            "largest_component_edge_density": diagnostics[0]["edge_density"]}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
