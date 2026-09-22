"""Create physical-bin balanced, prefix-nested marker ranks without phenotypes."""

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path

import numpy as np

from at_pheno.formal_provenance import file_sha256

SALTS = ("at_pheno_common_v1", "at_pheno_common_v1_s2", "at_pheno_common_v1_s3")
ALGORITHM = "physical_bin_sha256_round_robin_v1"
BIN_BP = 1_000_000
DTYPE = np.dtype([("chromosome", "u1"), ("position", "<u4"), ("ref_code", "u1"), ("alt_code", "u1"), ("source_record_index", "<u4")])
BASES = "ACGT"


def strata(compact):
    starts = [0]
    previous = (int(compact[0]["chromosome"]), (int(compact[0]["position"]) - 1) // BIN_BP)
    for index in range(1, len(compact)):
        current = (int(compact[index]["chromosome"]), (int(compact[index]["position"]) - 1) // BIN_BP)
        if current != previous:
            starts.append(index)
            previous = current
    starts.append(len(compact))
    return list(zip(starts[:-1], starts[1:]))


def hashed_order(compact, start, stop, salt):
    """Return global column indices in one stratum's salt rank order."""
    digests = np.empty(stop - start, dtype="S32")
    for local, row in enumerate(compact[start:stop]):
        key = f"{int(row['chromosome'])}:{int(row['position'])}:{BASES[int(row['ref_code'])]}:{BASES[int(row['alt_code'])]}"
        digests[local] = hashlib.sha256((salt + "\0" + key).encode()).digest()
    return start + np.argsort(digests, kind="stable")


def physical_rank(compact, salt, round_size=4096):
    """Exact round-robin interleave of hash-sorted chromosome×1-Mb strata."""
    if compact.dtype != DTYPE or not len(compact):
        raise ValueError("invalid variants_compact index")
    queues = [hashed_order(compact, start, stop, salt) for start, stop in strata(compact)]
    rank = np.empty(len(compact), dtype=np.int32)
    cursor, longest = 0, max(map(len, queues))
    for offset in range(0, longest, round_size):
        width = min(round_size, longest - offset)
        table = np.full((width, len(queues)), -1, dtype=np.int32)
        for column, queue in enumerate(queues):
            values = queue[offset:offset + width]
            table[:len(values), column] = values
        values = table.ravel()
        values = values[values >= 0]
        rank[cursor:cursor + len(values)] = values
        cursor += len(values)
    if cursor != len(compact) or len(np.unique(rank)) != len(rank):
        raise ValueError("round-robin rank did not form a permutation")
    return rank, len(queues)


def _hashed_order_worker(arguments):
    compact_path, start, stop, salt = arguments
    compact = np.load(compact_path, mmap_mode="r", allow_pickle=False)
    return hashed_order(compact, start, stop, salt)


def physical_rank_parallel(compact_path, salt, workers, round_size=4096):
    """Parallel stratum hashing with the same deterministic parent interleave."""
    compact = np.load(compact_path, mmap_mode="r", allow_pickle=False)
    bounds = strata(compact)
    if workers <= 1:
        return physical_rank(compact, salt, round_size=round_size)
    with ProcessPoolExecutor(max_workers=min(workers, len(bounds))) as executor:
        queues = list(executor.map(_hashed_order_worker,
                                   [(str(compact_path), start, stop, salt) for start, stop in bounds]))
    rank = np.empty(len(compact), dtype=np.int32)
    cursor, longest = 0, max(map(len, queues))
    for offset in range(0, longest, round_size):
        width = min(round_size, longest - offset)
        table = np.full((width, len(queues)), -1, dtype=np.int32)
        for column, queue in enumerate(queues):
            values = queue[offset:offset + width]
            table[:len(values), column] = values
        values = table.ravel()
        values = values[values >= 0]
        rank[cursor:cursor + len(values)] = values
        cursor += len(values)
    if cursor != len(compact) or len(np.unique(rank)) != len(rank):
        raise ValueError("round-robin rank did not form a permutation")
    return rank, len(queues)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compact", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    compact = np.load(args.compact, mmap_mode="r", allow_pickle=False)
    if compact.dtype != DTYPE:
        raise ValueError(f"variants_compact dtype mismatch: {compact.dtype}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = {"schema": "marker_rank_v1", "algorithm": ALGORITHM, "bin_bp": BIN_BP,
              "source_variants_compact_sha256": file_sha256(args.compact), "n_variants": int(len(compact)), "ranks": []}
    for salt in SALTS:
        rank, stratum_n = physical_rank_parallel(args.compact, salt, args.workers)
        path = args.out_dir / f"marker_rank_{salt}.npy"
        np.save(path, rank, allow_pickle=False)
        output["ranks"].append({"salt": salt, "path": str(path), "sha256": file_sha256(path),
                                "dtype": "int32", "n": int(len(rank)), "n_strata": stratum_n})
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "MARKER_RANKS_FROZEN", "n": output["n_variants"], "salts": list(SALTS)}))


if __name__ == "__main__":
    main()
