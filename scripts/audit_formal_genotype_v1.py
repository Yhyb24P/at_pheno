"""Independent post-build QA for the formal int8 ALT-dosage dataset.

This reader only accepts a directory carrying BUILD_COMPLETE, then reopens all
final artifacts.  It does not consume phenotypes and never runs prediction.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

import h5py
import numpy as np

from at_pheno.formal_provenance import file_sha256

DNA = "ACGT"
COMPACT_DTYPE = np.dtype([("chromosome", "u1"), ("position", "<u4"),
                          ("ref_code", "u1"), ("alt_code", "u1"),
                          ("source_record_index", "<u4")])


def quantiles(values):
    return {str(q): float(np.quantile(values, q)) for q in (0, .01, .05, .5, .95, .99, 1)}


def top_rows(ids, values, n=20):
    return [{"accession_id": ids[index], "rate": float(values[index])}
            for index in np.argsort(-values, kind="stable")[:n]]


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def vcf_samples(bcftools, vcf):
    result = subprocess.run([str(bcftools), "query", "-l", str(vcf)], check=True,
                            capture_output=True, text=True, encoding="utf-8")
    return result.stdout.splitlines()


def audit_matrix(path, n_samples, n_variants, block=8192):
    matrix = np.load(path, mmap_mode="r", allow_pickle=False)
    if matrix.dtype != np.int8 or matrix.shape != (n_samples, n_variants):
        raise ValueError(f"final matrix contract mismatch: {matrix.dtype} {matrix.shape}")
    values = {str(value): 0 for value in (-1, 0, 1, 2)}
    other = 0
    sample_missing = np.zeros(n_samples, dtype=np.int64)
    sample_het = np.zeros(n_samples, dtype=np.int64)
    variant_missing = np.zeros(n_variants, dtype=np.int32)
    variant_het = np.zeros(n_variants, dtype=np.int16)
    for start in range(0, n_variants, block):
        part = matrix[:, start:start + block]
        for value in (-1, 0, 1, 2):
            values[str(value)] += int(np.count_nonzero(part == value))
        other += int(np.count_nonzero((part < -1) | (part > 2)))
        missing = part == -1
        heterozygous = part == 1
        sample_missing += np.count_nonzero(missing, axis=1)
        sample_het += np.count_nonzero(heterozygous, axis=1)
        variant_missing[start:start + part.shape[1]] = np.count_nonzero(missing, axis=0)
        variant_het[start:start + part.shape[1]] = np.count_nonzero(heterozygous, axis=0)
    if other:
        raise ValueError(f"final matrix contains {other} values outside -1/0/1/2")
    return matrix, values, sample_missing / n_variants, sample_het / n_variants, variant_missing / n_samples, variant_het


def audit_variants(tsv, compact, n_variants):
    compact_array = np.load(compact, mmap_mode="r", allow_pickle=False)
    if compact_array.dtype != COMPACT_DTYPE or compact_array.shape != (n_variants,):
        raise ValueError("variants_compact.npy dtype/shape mismatch")
    counts = {str(chrom): 0 for chrom in range(1, 6)}
    previous = None
    source_previous = 0
    seen_positions = set()
    duplicate_positions = duplicate_quadruples = 0
    with Path(tsv).open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = ["chromosome", "position", "ref", "alt", "source_record_index"]
        if reader.fieldnames != required:
            raise ValueError(f"variants.tsv header mismatch: {reader.fieldnames}")
        for index, row in enumerate(reader):
            if index >= n_variants:
                raise ValueError("variants.tsv has more rows than matrix columns")
            chrom, position, ref, alt, source = (row[key] for key in required)
            key = (int(chrom), int(position), ref, alt)
            if chrom not in counts or int(position) < 1 or ref not in DNA or alt not in DNA or len(ref) != 1 or len(alt) != 1 or ref == alt:
                raise ValueError(f"invalid final variant row {row!r}")
            if previous is not None and key <= previous:
                raise ValueError("variants.tsv is not strictly ordered by chromosome, position, REF, ALT")
            if (int(chrom), int(position)) in seen_positions:
                duplicate_positions += 1
            seen_positions.add((int(chrom), int(position)))
            source_int = int(source)
            if source_int <= source_previous:
                raise ValueError("source_record_index is not strictly increasing")
            expected = compact_array[index]
            if tuple(expected) != (int(chrom), int(position), DNA.index(ref), DNA.index(alt), source_int):
                raise ValueError(f"variants_compact mismatch at matrix column {index}")
            counts[chrom] += 1
            previous, source_previous = key, source_int
        if index + 1 != n_variants:
            raise ValueError(f"variants.tsv rows {index + 1} != matrix columns {n_variants}")
    del compact_array
    return {"n": n_variants, "chromosome_counts": counts, "strict_order": True,
            "duplicate_chromosome_position": duplicate_positions,
            "duplicate_quadruples": duplicate_quadruples, "ref_alt_domain": "ACGT singleton unequal",
            "source_record_index_strictly_increasing": True}


def hdf5_spot_check(matrix, variants_tsv, hdf5_path, prior_spot):
    prior = json.loads(Path(prior_spot).read_text())["spot_check"]
    sites = [(str(row["chrom"]), int(row["pos"])) for row in prior["per_site"]]
    requested = set(sites)
    columns = {}
    with Path(variants_tsv).open(newline="", encoding="utf-8") as stream:
        for index, row in enumerate(csv.DictReader(stream, delimiter="\t")):
            key = (row["chromosome"], int(row["position"]))
            if key in requested:
                columns[key] = index
    if set(columns) != requested:
        raise ValueError("not all deterministic HDF5 spot sites occur in final variants.tsv")
    total = matches = comparable = 0
    with h5py.File(hdf5_path, "r") as handle:
        h5_ids = [value.decode() if isinstance(value, bytes) else str(value) for value in handle["accessions"][:]]
        if h5_ids != [row["accession_id"] for row in read_tsv(Path(variants_tsv).parent / "samples.tsv")]:
            raise ValueError("HDF5 accession order differs from final sample order")
        positions, regions = handle["positions"], handle["positions"].attrs["chr_regions"]
        for chrom, position in sites:
            start, stop = regions[int(chrom) - 1]
            local = int(np.searchsorted(positions[start:stop], position))
            if local >= stop - start or int(positions[start + local]) != position:
                raise ValueError(f"deterministic HDF5 spot coordinate missing: {(chrom, position)}")
            h5 = handle["snps"][start + local, :]
            final = matrix[:, columns[(chrom, position)]]
            mask = final != -1
            # The prior source-derived spot check consistently labels binary
            # 0 as REF.  Heterozygous final calls are non-comparable to binary.
            mask &= final != 1
            comparable += int(mask.sum())
            matches += int(np.count_nonzero(h5[mask] == (final[mask] == 2)))
            total += 1
    return {"sites": total, "comparable_pairs": comparable, "matches": matches,
            "concordance_hdf5_0_is_ref": matches / comparable if comparable else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, required=True)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--bcftools", type=Path, required=True)
    parser.add_argument("--hdf5", type=Path, required=True)
    parser.add_argument("--prior-spot", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    complete = args.dataset / "BUILD_COMPLETE.json"
    if not complete.is_file():
        raise ValueError("refusing QA: BUILD_COMPLETE.json is absent")
    build = json.loads((args.dataset / "genotype_build_manifest.json").read_text())
    n_samples, n_variants = build["shape"]
    samples = read_tsv(args.dataset / "samples.tsv")
    if len(samples) != n_samples or len({row["accession_id"] for row in samples}) != n_samples:
        raise ValueError("sample table is not unique/complete")
    ids = [row["accession_id"] for row in samples]
    vcf_ids = vcf_samples(args.bcftools, args.vcf)
    if ids != vcf_ids:
        raise ValueError("final sample order differs from VCF header")
    matrix, values, miss_s, het_s, miss_v, het_v = audit_matrix(args.dataset / "genotypes.npy", n_samples, n_variants)
    variants = audit_variants(args.dataset / "variants.tsv", args.dataset / "variants_compact.npy", n_variants)
    spot = hdf5_spot_check(matrix, args.dataset / "variants.tsv", args.hdf5, args.prior_spot)
    del matrix
    report = {"schema": "formal_genotype_v1_qc", "status": "FORMAL_GENOTYPE_QC_PASS",
              "matrix": {"shape": [n_samples, n_variants], "dtype": "int8", "file_bytes": (args.dataset / "genotypes.npy").stat().st_size,
                         "sha256": file_sha256(args.dataset / "genotypes.npy"), "value_counts": values, "other_values": 0},
              "samples": {"n": n_samples, "unique": n_samples, "ordered_equality_vcf_header": True,
                          "ordered_equality_hdf5_accessions": True, "genetic_group_mapping_complete": all(bool(row.get("genetic_group")) for row in samples)},
              "variants": variants,
              "missingness": {"per_sample_rate_quantiles": quantiles(miss_s), "top20_per_sample": top_rows(ids, miss_s),
                              "per_variant_call_rate_quantiles": quantiles(1 - miss_v)},
              "heterozygosity": {"per_sample_rate_quantiles": quantiles(het_s), "top20_per_sample": top_rows(ids, het_s),
                                  "per_variant_carrier_count_quantiles": quantiles(het_v)},
              "hdf5_cross_check": spot,
              "artifacts_sha256": {name: file_sha256(args.dataset / name) for name in ("samples.tsv", "variants.tsv", "variants_compact.npy", "genotype_build_manifest.json", "BUILD_COMPLETE.json")}}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "shape": report["matrix"]["shape"], "hdf5_concordance": spot["concordance_hdf5_0_is_ref"]}))


if __name__ == "__main__":
    main()
