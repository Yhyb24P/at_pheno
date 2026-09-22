"""Streaming G0 preflight for a formal 1001G VCF-to-int8 build.

This is intentionally a validation pass, not a converter: it never writes a
genotype matrix and never inspects phenotypes.  It implements the frozen G0
policies for original biallelic nuclear A/C/G/T SNP records, FILTER, dynamic
FORMAT/GT lookup, and accepted-record GT states.
"""

import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path
import subprocess
import sys
from contextlib import contextmanager

from at_pheno.formal_provenance import verify_source_vcf

NUCLEAR = {"1", "2", "3", "4", "5"}
DNA = {"A", "C", "G", "T"}
EXAMPLE_CAP = 4


@contextmanager
def open_vcf(path, pigz_threads=0):
    """Yield decoded VCF text without creating a plaintext temporary VCF.

    `pigz_threads=0` uses Python's portable gzip reader.  A positive value
    uses an explicitly requested pigz subprocess for parallel decompression;
    it is optional so fixture tests and CI do not depend on a system binary.
    """
    with Path(path).open("rb") as probe:
        magic = probe.read(3)
    if magic != b"\x1f\x8b\x08":
        with Path(path).open(encoding="utf-8", errors="strict") as stream:
            yield stream
        return
    if pigz_threads <= 0:
        with gzip.open(path, "rt", encoding="utf-8", errors="strict") as stream:
            yield stream
        return
    process = subprocess.Popen(["pigz", "-dc", "-p", str(pigz_threads), str(path)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                               encoding="utf-8", errors="strict")
    try:
        yield process.stdout
    finally:
        if process.stdout is not None:
            process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        code = process.wait()
        if code:
            raise ValueError(f"PREFLIGHT_BLOCKED: pigz exited {code}: {stderr.strip()}")


def structural_biallelic_snp(fields):
    if len(fields) < 9:
        return False
    chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
    try:
        int(pos)
    except ValueError:
        return False
    return (chrom in NUCLEAR and "," not in alt and len(ref) == len(alt) == 1
            and ref in DNA and alt in DNA and ref != alt)


def accepted_filter(value):
    return value in {"PASS", "."}


def gt_state(token):
    """Classify one dynamically extracted GT token under the G0 contract."""
    if not isinstance(token, str) or not token:
        return "malformed"
    separators = [sep for sep in "/|" if sep in token]
    if not separators and token in {"0", "1", "."}:
        return "haploid"
    if len(separators) != 1:
        return "malformed"
    separator = separators[0]
    parts = token.split(separator)
    if len(parts) == 1:
        return "haploid"
    if len(parts) != 2 or any(part == "" for part in parts):
        return "malformed"
    if parts == [".", "."]:
        return "fully_missing"
    if "." in parts:
        return "partial_missing"
    if not all(part.isdigit() for part in parts):
        return "malformed"
    values = [int(part) for part in parts]
    if any(value > 1 for value in values):
        return "allele_index_gt_1"
    if separator == "|":
        return "phased_heterozygous" if values[0] != values[1] else "phased_homozygous"
    return "unphased_heterozygous" if values[0] != values[1] else "unphased_homozygous"


def new_counts():
    return {"fully_missing": 0, "heterozygous": 0, "phased": 0, "unphased": 0,
            "partial_missing": 0, "haploid": 0, "allele_index_gt_1": 0, "malformed": 0}


def scan_stream(stream, expected_samples):
    filters_all, filters_candidates = Counter(), Counter()
    formats = Counter()
    format_no_gt = format_duplicate_gt = sample_width_mismatch = 0
    gt_counts, unsupported_examples = new_counts(), defaultdict(list)
    total = structural = accepted = 0
    sample_ids = None
    for line in stream:
        if line.startswith("#CHROM"):
            header = line.rstrip("\n").split("\t")
            if len(header) < 10 or header[8] != "FORMAT":
                raise ValueError("PREFLIGHT_BLOCKED: #CHROM lacks the literal FORMAT column")
            sample_ids = header[9:]
            continue
        if not line or line.startswith("#"):
            continue
        total += 1
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        filter_value = fields[6]
        filters_all[filter_value] += 1
        if not structural_biallelic_snp(fields):
            continue
        structural += 1
        filters_candidates[filter_value] += 1
        if not accepted_filter(filter_value):
            continue
        accepted += 1
        layout = fields[8]
        formats[layout] += 1
        keys = layout.split(":") if layout else []
        occurrences = [i for i, key in enumerate(keys) if key == "GT"]
        if not occurrences:
            format_no_gt += 1
            continue
        if len(occurrences) > 1:
            format_duplicate_gt += 1
            continue
        cells = fields[9:]
        if len(cells) != expected_samples:
            sample_width_mismatch += 1
            continue
        index = occurrences[0]
        # Most 1,135 accession rows contain only a handful of repeated
        # hard-call strings (0/0, 1/1, ./., etc.).  Count strings in C
        # first, then parse each distinct token once.  This is exactly
        # equivalent to a per-cell loop for all emitted count fields.
        for cell, occurrences_in_record in Counter(cells).items():
            parts = cell.split(":")
            token = parts[index] if index < len(parts) else ""
            state = gt_state(token)
            if state == "fully_missing":
                gt_counts[state] += occurrences_in_record
            elif state in {"phased_homozygous", "unphased_homozygous", "phased_heterozygous", "unphased_heterozygous"}:
                gt_counts["phased" if state.startswith("phased") else "unphased"] += occurrences_in_record
                if state.endswith("heterozygous"):
                    gt_counts["heterozygous"] += occurrences_in_record
            else:
                gt_counts[state] += occurrences_in_record
                if len(unsupported_examples[state]) < EXAMPLE_CAP:
                    unsupported_examples[state].append({"record": total, "sample_index": "one_or_more",
                                                        "token": token, "format": layout})
    if sample_ids is None:
        raise ValueError("PREFLIGHT_BLOCKED: missing #CHROM header")
    if len(sample_ids) != expected_samples:
        raise ValueError(f"PREFLIGHT_BLOCKED: expected {expected_samples} samples, found {len(sample_ids)}")
    hard = {"records_without_gt": format_no_gt, "records_duplicate_gt": format_duplicate_gt,
            "sample_field_width_mismatch": sample_width_mismatch,
            **{key: gt_counts[key] for key in ("partial_missing", "haploid", "allele_index_gt_1", "malformed")}}
    return {"n_samples": len(sample_ids), "sample_ids": sample_ids, "total_records": total,
            "structural_biallelic_snp_records": structural, "filter_exact_counts": dict(sorted(filters_all.items())),
            "filter_structural_candidate_counts": dict(sorted(filters_candidates.items())),
            "accepted_structural_records": accepted, "named_filter_excluded_records": structural - accepted,
            "format_layout_counts_accepted": dict(sorted(formats.items())), "records_without_gt": format_no_gt,
            "records_duplicate_gt": format_duplicate_gt, "sample_field_width_mismatch": sample_width_mismatch,
            "gt_state_counts_accepted": gt_counts, "unsupported_examples": dict(unsupported_examples),
            "hard_blockers": hard, "status": "FORMAL_VCF_PREFLIGHT_PASS" if not any(hard.values()) else "FILTER_POLICY_BLOCKED"}


def scan(vcf, expected_samples, pigz_threads=0):
    with open_vcf(vcf, pigz_threads=pigz_threads) as stream:
        return scan_stream(stream, expected_samples)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=1135)
    parser.add_argument("--pigz-threads", type=int, default=0,
                        help="positive: use pigz -dc -p N; zero: portable stdlib gzip")
    args = parser.parse_args()
    source = verify_source_vcf(args.source_manifest, args.vcf)
    report = scan(args.vcf, args.expected_samples, pigz_threads=args.pigz_threads)
    report.update({"source_identity": source, "source_manifest": str(args.source_manifest),
                   "filter_policy": {"accepted_exact_values": ["PASS", "."],
                                     "excluded": "any other non-empty named FILTER value"},
                   "format_policy": "dynamic FORMAT.split(':') lookup of unique GT"})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["status"], "accepted": report["accepted_structural_records"]}))
    if report["status"] != "FORMAL_VCF_PREFLIGHT_PASS":
        raise SystemExit(3)


if __name__ == "__main__":
    main()
