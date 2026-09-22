"""Build the formal 1001G v3.1 int8 ALT-dosage matrix without phenotype access.

The builder deliberately consumes two independently audited source products:
the raw VCF (whose bytes are reverified before build) and the ordered
biallelic-SNP catalog (whose key stream is compared line-for-line with
HTSlib's constrained VCF query).  It writes a sample-by-variant int8 memmap
directly; no full float matrix or Python list of variant arrays exists.
"""

import argparse
import csv
import gzip
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

import numpy as np

from at_pheno.formal_provenance import file_sha256, verify_source_vcf

CHROMS = ("1", "2", "3", "4", "5")
DNA_CODE = {"A": 0, "C": 1, "G": 2, "T": 3}
COMPACT_DTYPE = np.dtype([("chromosome", "u1"), ("position", "<u4"),
                          ("ref_code", "u1"), ("alt_code", "u1"),
                          ("source_record_index", "<u4")])
QUERY_FILTER = '(FILTER="PASS" || FILTER=".") && N_ALT=1 && TYPE="snp"'


def catalog_rows(path):
    """Yield validated catalog rows in its documented canonical order."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        expected = ["chromosome", "position", "ref", "alt", "source_record_index"]
        if reader.fieldnames != expected:
            raise ValueError(f"catalog header must be {expected}, got {reader.fieldnames}")
        previous = None
        for row in reader:
            chrom, pos, ref, alt = row["chromosome"], int(row["position"]), row["ref"], row["alt"]
            source_index = int(row["source_record_index"])
            key = (int(chrom), pos, ref, alt)
            if chrom not in CHROMS or pos < 1 or ref not in DNA_CODE or alt not in DNA_CODE or ref == alt:
                raise ValueError(f"invalid catalog row {row!r}")
            if previous is not None and key <= previous:
                raise ValueError("catalog is not strictly ordered by chromosome, position, REF, ALT")
            previous = key
            yield chrom, pos, ref, alt, source_index


def catalog_counts(path):
    counts = {chrom: 0 for chrom in CHROMS}
    for chrom, *_ in catalog_rows(path):
        counts[chrom] += 1
    if not all(counts.values()):
        raise ValueError(f"catalog has an empty nuclear chromosome: {counts}")
    return counts


def parse_query_line(line, expected_samples):
    """Parse one fixed-width `%GT` HTSlib-query line into metadata and int8 dosage.

    The preflight has already proven the unique source FORMAT layout and that
    every accepted GT is a hard diploid biallelic call or fully missing.  This
    byte-level parser rechecks the fixed-width `%GT` representation while
    converting it: `0/0`, `0|0`, `0/1`, `1/0`, `1/1`, `./.` only.
    """
    fields = line.rstrip(b"\n").split(b"\t", 4)
    if len(fields) != 5:
        raise ValueError("HTSlib query row lacks CHROM/POS/REF/ALT/GT fields")
    chrom_b, position_b, ref_b, alt_b, packed = fields
    try:
        chrom, position, ref, alt = chrom_b.decode(), int(position_b), ref_b.decode(), alt_b.decode()
    except UnicodeDecodeError as exc:
        raise ValueError("non-ASCII VCF query metadata") from exc
    if chrom not in CHROMS or position < 1 or ref not in DNA_CODE or alt not in DNA_CODE or ref == alt:
        raise ValueError(f"HTSlib query violated formal structural policy: {(chrom, position, ref, alt)!r}")
    raw = np.frombuffer(packed, dtype=np.uint8)
    # There are 1,134 TAB separators after the first GT, hence 4*n-1 bytes.
    if raw.size != 4 * expected_samples - 1:
        raise ValueError(f"query GT payload has {raw.size} bytes, expected {4 * expected_samples - 1}")
    gt = raw[0::4]
    tail = raw[2::4]
    separators = raw[1::4]
    if np.any((separators != ord("/")) & (separators != ord("|"))):
        raise ValueError("query GT payload has a non-diploid separator")
    if np.any(raw[3::4] != ord("\t")):
        raise ValueError("query GT payload is not fixed-width tab-separated GT")
    valid = ((gt == ord("0")) | (gt == ord("1")) | (gt == ord("."))) & ((tail == ord("0")) | (tail == ord("1")) | (tail == ord(".")))
    if not np.all(valid):
        raise ValueError("query GT payload has an unsupported allele token")
    missing = (gt == ord(".")) & (tail == ord("."))
    partial = ((gt == ord(".")) != (tail == ord(".")))
    if np.any(partial):
        raise ValueError("query GT payload has partial missing GT")
    result = ((gt == ord("1")).astype(np.int8) + (tail == ord("1")).astype(np.int8))
    result[missing] = -1
    return (chrom, position, ref, alt), result


def _catalog_chrom_rows(path, chrom):
    return (row for row in catalog_rows(path) if row[0] == chrom)


def write_chromosome(vcf, bcftools, catalog, matrix_path, chrom, offset, expected_n, expected_samples):
    """Write one non-overlapping memmap column span and verify every key."""
    command = [str(bcftools), "query", "-r", chrom, "-i", QUERY_FILTER,
               "-f", "%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n", str(vcf)]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    matrix = np.load(matrix_path, mmap_mode="r+")
    seen = 0
    catalog_iter = _catalog_chrom_rows(catalog, chrom)
    try:
        for raw_line in process.stdout:
            key, dosage = parse_query_line(raw_line, expected_samples)
            try:
                expected = next(catalog_iter)
            except StopIteration as exc:
                raise ValueError(f"{chrom}: HTSlib query has extra key {key!r}") from exc
            expected_key = expected[:4]
            if key != expected_key:
                raise ValueError(f"{chrom}: query/catalog key mismatch: {key!r} != {expected_key!r}")
            if seen >= expected_n:
                raise ValueError(f"{chrom}: more query rows than allocated ({expected_n})")
            matrix[:, offset + seen] = dosage
            seen += 1
        try:
            extra = next(catalog_iter)
        except StopIteration:
            extra = None
        if extra is not None:
            raise ValueError(f"{chrom}: catalog has extra key {extra[:4]!r}")
    finally:
        matrix.flush()
        del matrix
        if process.stdout is not None:
            process.stdout.close()
        stderr = process.stderr.read().decode("utf-8", errors="replace") if process.stderr else ""
        status = process.wait()
    if status:
        raise ValueError(f"{chrom}: bcftools query exited {status}: {stderr.strip()}")
    if seen != expected_n:
        raise ValueError(f"{chrom}: query rows {seen} != catalog rows {expected_n}")
    return chrom, seen


def read_vcf_samples(bcftools, vcf):
    result = subprocess.run([str(bcftools), "query", "-l", str(vcf)], check=True,
                            capture_output=True, text=True, encoding="utf-8")
    ids = result.stdout.splitlines()
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("VCF has missing or duplicate sample IDs")
    return ids


def admixture_groups(path, sample_ids):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    mapping = {row["id"]: row["group"] for row in rows}
    if len(mapping) != len(rows) or any(not group for group in mapping.values()):
        raise ValueError("admixture metadata has duplicate IDs or blank groups")
    if set(mapping) != set(sample_ids):
        raise ValueError("VCF sample IDs and admixture metadata IDs differ")
    return mapping


def write_metadata(catalog, output, sample_ids, groups):
    samples = output / "samples.tsv.partial"
    with samples.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["accession_id", "genetic_group"], delimiter="\t")
        writer.writeheader()
        writer.writerows({"accession_id": sample, "genetic_group": groups[sample]} for sample in sample_ids)
    rows = sum(1 for _ in catalog_rows(catalog))
    variants = output / "variants.tsv.partial"
    compact = output / "variants_compact.npy.partial"
    array = np.lib.format.open_memmap(compact, mode="w+", dtype=COMPACT_DTYPE, shape=(rows,))
    with variants.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["chromosome", "position", "ref", "alt", "source_record_index"], delimiter="\t")
        writer.writeheader()
        for index, (chrom, position, ref, alt, source_index) in enumerate(catalog_rows(catalog)):
            writer.writerow({"chromosome": chrom, "position": position, "ref": ref, "alt": alt,
                             "source_record_index": source_index})
            array[index] = (int(chrom), position, DNA_CODE[ref], DNA_CODE[alt], source_index)
    array.flush()
    del array
    return rows


def value_counts(matrix_path, block=8192):
    matrix = np.load(matrix_path, mmap_mode="r")
    result = {str(value): 0 for value in (-1, 0, 1, 2)}
    other = 0
    for start in range(0, matrix.shape[1], block):
        part = matrix[:, start:start + block]
        for value in (-1, 0, 1, 2):
            result[str(value)] += int(np.count_nonzero(part == value))
        other += int(np.count_nonzero((part < -1) | (part > 2)))
    return result, other


def finalize(output, expected_samples, expected_variants, source_identity, catalog, catalog_summary):
    matrix_partial = output / "genotypes.npy.partial"
    counts, other = value_counts(matrix_partial)
    matrix = np.load(matrix_partial, mmap_mode="r")
    if matrix.dtype != np.int8 or matrix.shape != (expected_samples, expected_variants) or other:
        raise ValueError("partial matrix failed final dtype/shape/value contract")
    del matrix
    for name in ("genotypes.npy", "samples.tsv", "variants.tsv", "variants_compact.npy"):
        os.replace(output / f"{name}.partial", output / name)
    build = {"schema": "formal_genotype_v1", "status": "BUILD_COMPLETE", "n_samples": expected_samples,
             "n_variants": expected_variants, "dtype": "int8", "shape": [expected_samples, expected_variants],
             "value_counts": counts, "source_identity": source_identity, "catalog_sha256": file_sha256(catalog),
             "catalog_summary_sha256": file_sha256(catalog_summary), "query_filter": QUERY_FILTER,
             "matrix_sha256": file_sha256(output / "genotypes.npy"),
             "samples_sha256": file_sha256(output / "samples.tsv"),
             "variants_sha256": file_sha256(output / "variants.tsv"),
             "variants_compact_sha256": file_sha256(output / "variants_compact.npy")}
    (output / "genotype_build_manifest.json").write_text(json.dumps(build, indent=2) + "\n")
    (output / "genotype_qc.json").write_text(json.dumps({"matrix": {"shape": build["shape"], "dtype": "int8",
        "value_counts": counts, "other_values": other}, "preliminary": True}, indent=2) + "\n")
    (output / "BUILD_COMPLETE.json").write_text(json.dumps({"status": "BUILD_COMPLETE", "build_manifest_sha256": file_sha256(output / "genotype_build_manifest.json")}, indent=2) + "\n")
    return build


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--catalog-summary", type=Path, required=True)
    parser.add_argument("--admixture", type=Path, required=True)
    parser.add_argument("--bcftools", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=5)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError(f"refusing to overwrite existing output directory {args.out}")
    preflight = json.loads(args.preflight.read_text())
    if preflight.get("status") != "FORMAL_VCF_PREFLIGHT_PASS":
        raise ValueError("formal build requires FORMAL_VCF_PREFLIGHT_PASS")
    source_identity = verify_source_vcf(args.source_manifest, args.vcf)
    counts = catalog_counts(args.catalog)
    expected = sum(counts.values())
    if expected != preflight.get("accepted_structural_records"):
        raise ValueError(f"catalog rows {expected} != preflight accepted records {preflight.get('accepted_structural_records')}")
    samples = read_vcf_samples(args.bcftools, args.vcf)
    if len(samples) != preflight.get("n_samples"):
        raise ValueError("VCF sample count differs from preflight")
    groups = admixture_groups(args.admixture, samples)
    args.out.mkdir(parents=True)
    try:
        metadata_rows = write_metadata(args.catalog, args.out, samples, groups)
        if metadata_rows != expected:
            raise ValueError("metadata catalog count changed during build")
        matrix_path = args.out / "genotypes.npy.partial"
        matrix = np.lib.format.open_memmap(matrix_path, mode="w+", dtype=np.int8, shape=(len(samples), expected))
        matrix.flush()
        del matrix
        offsets, cursor = {}, 0
        for chrom in CHROMS:
            offsets[chrom] = cursor
            cursor += counts[chrom]
        if args.workers < 1 or args.workers > len(CHROMS):
            raise ValueError("workers must be 1..5")
        job_args = [(args.vcf, args.bcftools, args.catalog, matrix_path, chrom, offsets[chrom], counts[chrom], len(samples)) for chrom in CHROMS]
        if args.workers == 1:
            results = [write_chromosome(*job) for job in job_args]
        else:
            with multiprocessing.get_context("spawn").Pool(args.workers) as pool:
                results = pool.starmap(write_chromosome, job_args)
        if dict(results) != counts:
            raise ValueError(f"chromosome build counts differ: {results!r} != {counts!r}")
        build = finalize(args.out, len(samples), expected, source_identity, args.catalog, args.catalog_summary)
    except BaseException:
        # Keep partial evidence for diagnosis, but never create BUILD_COMPLETE.
        raise
    print(json.dumps({"status": build["status"], "shape": build["shape"], "matrix_sha256": build["matrix_sha256"]}))


if __name__ == "__main__":
    main()
