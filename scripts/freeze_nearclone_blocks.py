"""Freeze phenotype-independent HDF5 binary-state near-clone blocks."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np

CHROMS = ("1", "2", "3", "4", "5")
SALT = "at_pheno_nearclone_v1"
CUTOFFS = (0.0005, 0.001, 0.002)


def rank_key(chrom, position, salt=SALT):
    return hashlib.sha256(f"{salt}:{chrom}:{position}".encode()).digest()


def select_panel(positions, regions, panel_size):
    """Hash-rank unique HDF5 coordinates, returning physical-order indices."""
    records = []
    for chrom_index, (start, stop) in enumerate(regions, 1):
        for index, position in enumerate(positions[start:stop], int(start)):
            records.append((rank_key(str(chrom_index), int(position)), chrom_index, int(position), index))
    if len(records) < panel_size:
        raise ValueError("HDF5 has fewer coordinates than requested panel")
    chosen = sorted(records, key=lambda row: row[0])[:panel_size]
    per_chrom = {str(chrom): sum(row[1] == chrom for row in chosen) for chrom in range(1, 6)}
    if not all(per_chrom.values()):
        raise ValueError(f"hash marker panel lacks chromosome coverage: {per_chrom}")
    return sorted(chosen, key=lambda row: (row[1], row[2])), per_chrom


def hamming_counts(handle, indices, block=4096):
    """Exact all-pair binary Hamming counts with bounded HDF5/floating blocks."""
    n = handle["snps"].shape[1]
    ones = np.zeros(n, dtype=np.int64)
    gram = np.zeros((n, n), dtype=np.int64)
    # The selected indices are physically ordered; coalesce contiguous runs
    # only through ordinary bounded slices to keep memory below 25 MiB/block.
    for start in range(0, len(indices), block):
        selected = indices[start:start + block]
        values = handle["snps"][selected, :]
        if np.any((values != 0) & (values != 1)):
            raise ValueError("HDF5 binary-state contract violated in selected panel")
        binary = values.T.astype(np.float32, copy=False)
        ones += np.sum(values, axis=0, dtype=np.int64)
        gram += np.rint(binary @ binary.T).astype(np.int64)
    return ones[:, None] + ones[None, :] - 2 * gram


def components(distance, cutoff):
    n = distance.shape[0]
    active = np.triu(distance <= cutoff, k=1)
    neighbours = [np.flatnonzero(active[i] | active[:, i]) for i in range(n)]
    seen, result = np.zeros(n, dtype=bool), []
    for start in range(n):
        if seen[start]:
            continue
        stack, group = [start], []
        seen[start] = True
        while stack:
            current = stack.pop()
            group.append(current)
            for neighbour in neighbours[current]:
                if not seen[neighbour]:
                    seen[neighbour] = True
                    stack.append(int(neighbour))
        result.append(sorted(group))
    return sorted(result, key=lambda group: (-len(group), group))


def diagnostics(distance, groups, cutoff, ids, group_labels):
    output = []
    for block_id, members in enumerate(groups):
        size = len(members)
        if size == 1:
            max_distance, density, edge_count = 0.0, 1.0, 0
        else:
            values = distance[np.ix_(members, members)]
            lower = values[np.tril_indices(size, -1)]
            edge_count = int(np.count_nonzero(lower <= cutoff))
            density = edge_count / len(lower)
            max_distance = float(np.max(lower))
        composition = {}
        for member in members:
            label = group_labels[ids[member]]
            composition[label] = composition.get(label, 0) + 1
        output.append({"block_id": block_id, "size": size, "edge_count": edge_count,
                       "edge_density": density, "max_within_hamming_count": max_distance,
                       "max_within_hamming_distance": max_distance / PANEL_SIZE,
                       "admixture_composition": json.dumps(dict(sorted(composition.items())))})
    return output


PANEL_SIZE = 250_000


def read_groups(samples):
    with Path(samples).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    ids = [row["accession_id"] for row in rows]
    labels = {row["accession_id"]: row["genetic_group"] for row in rows}
    if len(ids) != len(set(ids)) or any(not labels[item] for item in ids):
        raise ValueError("sample manifest has duplicate IDs or blank groups")
    return ids, labels


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite {args.out}")
    sample_ids, labels = read_groups(args.samples)
    args.out.mkdir(parents=True)
    with h5py.File(args.hdf5, "r") as handle:
        h5_ids = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["accessions"][:]]
        if h5_ids != sample_ids:
            raise ValueError("HDF5 and formal sample manifest have different accession order")
        positions = handle["positions"][:]
        regions = [tuple(map(int, row)) for row in handle["positions"].attrs["chr_regions"]]
        selected, per_chrom = select_panel(positions, regions, PANEL_SIZE)
        with (args.out / "nearclone_marker_panel.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=["chromosome", "position", "hdf5_row", "rank_sha256"], delimiter="\t")
            writer.writeheader()
            for digest, chrom, position, index in selected:
                writer.writerow({"chromosome": chrom, "position": position, "hdf5_row": index, "rank_sha256": digest.hex()})
        distance_count = hamming_counts(handle, [row[3] for row in selected])
    if not np.all(distance_count == distance_count.T) or np.any(np.diag(distance_count)):
        raise ValueError("invalid Hamming matrix")
    summary = {"schema": "nearclone_hamming_v1", "salt": SALT, "panel_size": PANEL_SIZE,
               "panel_per_chromosome": per_chrom, "sample_n": len(sample_ids), "cutoffs": {},
               "sampling_design_exception": "genotype-only coordinate hash panel; no phenotype values used"}
    for cutoff in CUTOFFS:
        max_count = int(np.floor(cutoff * PANEL_SIZE + 1e-12))
        groups = components(distance_count, max_count)
        records = diagnostics(distance_count, groups, max_count, sample_ids, labels)
        if cutoff == 0.001:
            with (args.out / "nearclone_blocks_d001.tsv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["accession_id", "nearclone_block"], delimiter="\t")
                writer.writeheader()
                for block_id, members in enumerate(groups):
                    writer.writerows({"accession_id": sample_ids[index], "nearclone_block": block_id} for index in members)
            with (args.out / "nearclone_component_diagnostics.tsv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(records)
        sizes = [len(group) for group in groups]
        largest = records[0]
        chain = (largest["size"] > 1 and largest["edge_density"] < .25
                 and largest["max_within_hamming_distance"] > .001)
        summary["cutoffs"][str(cutoff)] = {"max_hamming_count": max_count, "n_edges": int(np.count_nonzero(np.triu(distance_count <= max_count, 1))),
            "n_components": len(groups), "multi_member_components": sum(size > 1 for size in sizes),
            "largest_component": max(sizes), "chain_flag": chain,
            "status": "NEARCLONE_CHAIN_BLOCKER" if max(sizes) > .25 * len(sample_ids) or chain else "PASS"}
    (args.out / "nearclone_sensitivity.json").write_text(json.dumps(summary, indent=2) + "\n")
    primary = summary["cutoffs"]["0.001"]
    print(json.dumps({"status": primary["status"], "largest_component": primary["largest_component"], "n_edges": primary["n_edges"]}))
    if primary["status"] != "PASS":
        raise SystemExit(3)


if __name__ == "__main__":
    main()
