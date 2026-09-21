"""Full-file QC scan for the 1001G orientation-unknown binary HDF5 product.

It verifies value range and writes chromosome-wise binary-state summaries.
It deliberately does not pretend binary 1 means ALT, and does not run an
O(n²M) all-marker relationship calculation on CPU.
"""

import argparse
import hashlib
import json
from pathlib import Path
import time

import h5py
import numpy as np


def stream_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize_binary_block(block: np.ndarray, n_accessions: int) -> tuple[dict[str, int], np.ndarray]:
    """Return value-range counts and histogram of state-1 carrier counts."""
    counts = {"zero": int(np.sum(block == 0)), "one": int(np.sum(block == 1))}
    counts["other"] = int(block.size - counts["zero"] - counts["one"])
    histogram = np.bincount(np.sum(block == 1, axis=1), minlength=n_accessions + 1)
    return counts, histogram


def minor_state_ge(histogram: np.ndarray, n_accessions: int, minimum: int) -> int:
    state_one_n = np.arange(len(histogram))
    minority_state_n = np.minimum(state_one_n, n_accessions - state_one_n)
    return int(histogram[minority_state_n >= minimum].sum())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hdf5", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--block-markers", type=int, default=10_000)
    args = parser.parse_args()
    if args.out.exists() or args.block_markers < 1:
        raise ValueError("Output exists or invalid block size")
    args.out.mkdir(parents=True)
    start = time.monotonic()
    with h5py.File(args.hdf5, "r") as f:
        snps, positions = f["snps"], f["positions"]
        n_accessions = int(snps.shape[1])
        regions = positions.attrs["chr_regions"]
        chromosomes = positions.attrs["chrs"]
        rows = []
        for chrom, (left, right) in zip(chromosomes, regions):
            chrom = chrom.decode() if isinstance(chrom, bytes) else str(chrom)
            left, right = int(left), int(right)
            counts = {"zero": 0, "one": 0, "other": 0}
            carrier_hist = np.zeros(n_accessions+1, dtype=np.int64)
            for begin in range(left, right, args.block_markers):
                block = np.asarray(snps[begin:min(begin+args.block_markers, right), :])
                block_counts, block_hist = summarize_binary_block(block, n_accessions)
                for key in counts:
                    counts[key] += block_counts[key]
                carrier_hist += block_hist
            np.save(args.out/f"chr{str(chrom)}_one_carrier_histogram.npy", carrier_hist)
            rows.append({"chromosome": chrom, "markers": right-left, **counts,
                         "variable_markers": int(np.sum(carrier_hist[1:-1])),
                         "minor_state_ge_5_markers": minor_state_ge(carrier_hist, n_accessions, 5)})
    source_hash = stream_sha256(args.hdf5)
    report = {"input": str(args.hdf5), "input_sha256": source_hash, "representation": "orientation_unknown_binary",
              "n_accessions": n_accessions, "block_markers": args.block_markers, "chromosomes": rows,
              "elapsed_seconds": time.monotonic()-start,
              "claims_not_supported": ["ALT effect direction", "original missingness", "fold-specific imputation",
                                       "REF/ALT-dependent annotation"],
              "next_step": "Use summaries for feasibility/kinship design; do not label as formal VCF branch."}
    (args.out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
