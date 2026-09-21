"""Full-file QC scan for the 1001G orientation-unknown binary HDF5 product.

It verifies value range and writes chromosome-wise allele/carrier summaries.
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
        regions = positions.attrs["chr_regions"]
        chromosomes = positions.attrs["chrs"]
        rows = []
        for chrom, (left, right) in zip(chromosomes, regions):
            chrom = chrom.decode() if isinstance(chrom, bytes) else str(chrom)
            left, right = int(left), int(right)
            counts = {"zero": 0, "one": 0, "other": 0}
            carrier_hist = np.zeros(snps.shape[1]+1, dtype=np.int64)
            for begin in range(left, right, args.block_markers):
                block = np.asarray(snps[begin:min(begin+args.block_markers, right), :])
                counts["zero"] += int(np.sum(block == 0))
                counts["one"] += int(np.sum(block == 1))
                counts["other"] += int(block.size-np.sum(block == 0)-np.sum(block == 1))
                carrier_hist += np.bincount(np.sum(block == 1, axis=1), minlength=snps.shape[1]+1)
            np.save(args.out/f"chr{str(chrom)}_one_carrier_histogram.npy", carrier_hist)
            rows.append({"chromosome": chrom, "markers": right-left, **counts,
                         "variable_markers": int(np.sum(carrier_hist[1:-1])),
                         "mac_ge_5_markers": int(np.sum(carrier_hist[5:-4]))})
    source_hash = hashlib.sha256(args.hdf5.read_bytes()).hexdigest()
    report = {"input": str(args.hdf5), "input_sha256": source_hash, "representation": "orientation_unknown_binary",
              "n_accessions": 1135, "block_markers": args.block_markers, "chromosomes": rows,
              "elapsed_seconds": time.monotonic()-start,
              "claims_not_supported": ["ALT effect direction", "original missingness", "fold-specific imputation",
                                       "REF/ALT-dependent annotation"],
              "next_step": "Use summaries for feasibility/kinship design; do not label as formal VCF branch."}
    (args.out/"report.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
