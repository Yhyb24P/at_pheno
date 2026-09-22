"""Confirm a selected-panel near-clone chain on all binary HDF5 states."""

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--blocks", type=Path, required=True)
    parser.add_argument("--block-id", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with args.samples.open(newline="", encoding="utf-8") as stream:
        ids = [row["accession_id"] for row in csv.DictReader(stream, delimiter="\t")]
    index = {item: pos for pos, item in enumerate(ids)}
    with args.blocks.open(newline="", encoding="utf-8") as stream:
        members = [row["accession_id"] for row in csv.DictReader(stream, delimiter="\t")
                   if int(row["nearclone_block"]) == args.block_id]
    if len(members) < 2 or any(item not in index for item in members):
        raise ValueError("requested block is not a valid multi-member component")
    member_indices = [index[item] for item in members]
    with h5py.File(args.hdf5, "r") as handle:
        h5_ids = [item.decode() if isinstance(item, bytes) else str(item) for item in handle["accessions"][:]]
        if h5_ids != ids:
            raise ValueError("HDF5 sample order differs from final samples")
        n_sites = handle["snps"].shape[0]
        ones = np.zeros(len(members), dtype=np.int64)
        gram = np.zeros((len(members), len(members)), dtype=np.int64)
        for start in range(0, n_sites, 10000):
            values = handle["snps"][start:start + 10000, member_indices]
            if np.any((values != 0) & (values != 1)):
                raise ValueError("HDF5 binary contract violated")
            binary = values.astype(np.float32, copy=False)
            ones += values.sum(axis=0, dtype=np.int64)
            gram += np.rint(binary.T @ binary).astype(np.int64)
    hamming = ones[:, None] + ones[None, :] - 2 * gram
    lower = hamming[np.tril_indices(len(members), -1)]
    edge_count = int(np.count_nonzero(lower <= .001 * n_sites))
    result = {"status": "FULL_HDF5_CHAIN_CONFIRMATION", "block_id": args.block_id,
              "members": members, "n_sites": n_sites, "cutoff": .001,
              "edge_count": edge_count, "edge_density": edge_count / len(lower),
              "max_hamming_count": int(lower.max()), "max_hamming_distance": float(lower.max() / n_sites),
              "min_hamming_distance": float(lower.min() / n_sites)}
    args.out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
