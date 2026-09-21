"""Read HDF5 axes, attributes and bounded genotype blocks without recoding alleles."""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def serial(value):
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    if isinstance(value, np.ndarray):
        return [serial(x) for x in value.tolist()]
    if isinstance(value, np.generic):
        return serial(value.item())
    if isinstance(value, (tuple, list)):
        return [serial(x) for x in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("hdf5", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    report = {"input": str(args.hdf5), "bytes": args.hdf5.stat().st_size,
              "schema": {}, "root_attributes": {}, "recoded": False}
    with h5py.File(args.hdf5, "r") as f:
        report["root_attributes"] = {k: serial(v) for k, v in f.attrs.items()}
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                report["schema"][name] = {"shape": list(obj.shape), "dtype": str(obj.dtype),
                    "chunks": obj.chunks, "compression": obj.compression,
                    "attributes": {k: serial(v) for k, v in obj.attrs.items()}}
        f.visititems(visit)
        ids = [str(serial(x)) for x in f["accessions"][:]]
        report["accession_n"] = len(ids)
        report["unique_accession_n"] = len(set(ids))
        report["accessions"] = ids
        positions = f["positions"][:]
        regions = f["positions"].attrs.get("chr_regions")
        report["positions_n"] = len(positions)
        if regions is not None:
            report["chromosome_regions"] = []
            for index, (start, stop) in enumerate(regions):
                start, stop = int(start), int(stop)
                pos = positions[start:stop]
                report["chromosome_regions"].append({"region_index": index, "start": start,
                    "stop": stop, "count": len(pos), "min_position": int(np.min(pos)),
                    "max_position": int(np.max(pos)), "strictly_increasing": bool(np.all(np.diff(pos.astype(np.int64)) > 0))})
        snps = f["snps"]
        if snps.shape != (len(positions), len(ids)):
            raise ValueError("Unexpected SNP matrix axis order; inspect schema")
        values = Counter()
        sampled_rows = []
        for start in np.linspace(0, max(0, snps.shape[0]-256), 16, dtype=int):
            block = snps[start:start+256, :]
            v, count = np.unique(block, return_counts=True)
            values.update({str(serial(a)): int(b) for a, b in zip(v, count)})
            sampled_rows.append([int(start), min(int(start)+256, snps.shape[0])])
        report["genotype_sample_ranges"] = sampled_rows
        report["genotype_sample_value_counts"] = dict(values)
        report["genotype_sample_only"] = True
        report["allele_orientation"] = "not inferred from binary values; inspect metadata/source"
    digest = hashlib.sha256()
    with args.hdf5.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            digest.update(chunk)
    report["sha256"] = digest.hexdigest()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        json.dump(report, stream, indent=2, ensure_ascii=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in report.items() if k not in {"accessions"}}, indent=2))


if __name__ == "__main__":
    main()
