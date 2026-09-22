"""Run G0 VCF preflight across indexed chromosomes with HTSlib in parallel."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
from pathlib import Path
import subprocess

from at_pheno.formal_provenance import verify_source_vcf
from formal_vcf_preflight import EXAMPLE_CAP, new_counts, scan_stream

CHROMS = ("1", "2", "3", "4", "5")
# TAIR10 nuclear chromosome lengths.  The indexed VCF is partitioned into
# disjoint closed coordinate intervals, so each original record is streamed
# exactly once while independent Python parsers can run in separate processes.
CHROM_LENGTHS = {"1": 30_427_671, "2": 19_698_289, "3": 23_459_830,
                 "4": 18_585_056, "5": 26_975_502}


def build_regions(workers):
    """Return approximately equal physical-span, non-overlapping VCF regions."""
    if workers < len(CHROMS):
        raise ValueError(f"workers must be at least {len(CHROMS)}")
    total = sum(CHROM_LENGTHS.values())
    pieces = {chrom: max(1, round(workers * CHROM_LENGTHS[chrom] / total)) for chrom in CHROMS}
    # Round-to-nearest can miss the requested process budget.  Add/remove the
    # longest/shortest next chunks deterministically until it matches.
    while sum(pieces.values()) < workers:
        chrom = max(CHROMS, key=lambda key: CHROM_LENGTHS[key] / pieces[key])
        pieces[chrom] += 1
    while sum(pieces.values()) > workers:
        chrom = max((key for key in CHROMS if pieces[key] > 1),
                    key=lambda key: CHROM_LENGTHS[key] / pieces[key])
        pieces[chrom] -= 1
    regions = []
    for chrom in CHROMS:
        width, count = CHROM_LENGTHS[chrom], pieces[chrom]
        for index in range(count):
            start = index * width // count + 1
            end = (index + 1) * width // count
            regions.append((f"{chrom}:{start}-{end}", chrom))
    return regions


def scan_region(bcftools, vcf, region, chrom, expected_samples, threads):
    process = subprocess.Popen([str(bcftools), "view", "--threads", str(threads), "-r", region, str(vcf)],
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, encoding="utf-8")
    try:
        report = scan_stream(process.stdout, expected_samples)
    finally:
        if process.stdout is not None:
            process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        code = process.wait()
    if code:
        raise ValueError(f"PREFLIGHT_BLOCKED: bcftools region {region} exited {code}: {stderr.strip()}")
    return region, chrom, report


def merge_reports(reports):
    # Accept the original chromosome-keyed mapping in unit tests and the
    # region triplets emitted by the multiprocess runner.
    if isinstance(reports, dict):
        report_rows = [(chrom, chrom, reports[chrom]) for chrom in CHROMS]
    else:
        report_rows = list(reports)
    only_reports = [report for _, _, report in report_rows]
    sample_ids = only_reports[0]["sample_ids"]
    if any(report["sample_ids"] != sample_ids for report in only_reports[1:]):
        raise ValueError("PREFLIGHT_BLOCKED: bcftools region streams disagree on sample order")
    result = {"n_samples": len(sample_ids), "sample_ids": sample_ids,
              "total_records": sum(r["total_records"] for r in only_reports),
              "structural_biallelic_snp_records": sum(r["structural_biallelic_snp_records"] for r in only_reports),
              "accepted_structural_records": sum(r["accepted_structural_records"] for r in only_reports),
              "named_filter_excluded_records": sum(r["named_filter_excluded_records"] for r in only_reports)}
    for field in ("filter_exact_counts", "filter_structural_candidate_counts", "format_layout_counts_accepted"):
        merged = Counter()
        for report in only_reports:
            merged.update(report[field])
        result[field] = dict(sorted(merged.items()))
    for field in ("records_without_gt", "records_duplicate_gt", "sample_field_width_mismatch"):
        result[field] = sum(r[field] for r in only_reports)
    result["gt_state_counts_accepted"] = new_counts()
    for report in only_reports:
        for key, value in report["gt_state_counts_accepted"].items():
            result["gt_state_counts_accepted"][key] += value
    examples = defaultdict(list)
    for region, chrom, report in report_rows:
        for kind, rows in report["unsupported_examples"].items():
            for row in rows:
                if len(examples[kind]) < EXAMPLE_CAP:
                    examples[kind].append({"chromosome": chrom, "region": region, **row})
    result["unsupported_examples"] = dict(examples)
    result["hard_blockers"] = {"records_without_gt": result["records_without_gt"],
        "records_duplicate_gt": result["records_duplicate_gt"], "sample_field_width_mismatch": result["sample_field_width_mismatch"],
        **{key: result["gt_state_counts_accepted"][key] for key in ("partial_missing", "haploid", "allele_index_gt_1", "malformed")}}
    result["status"] = "FORMAL_VCF_PREFLIGHT_PASS" if not any(result["hard_blockers"].values()) else "FILTER_POLICY_BLOCKED"
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--bcftools", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--expected-samples", type=int, default=1135)
    parser.add_argument("--workers", type=int, default=12,
                        help="independent indexed regions; each uses one parser process")
    parser.add_argument("--threads-per-region", type=int, default=1)
    args = parser.parse_args()
    source = verify_source_vcf(args.source_manifest, args.vcf)
    if args.threads_per_region < 1:
        raise ValueError("threads-per-region must be positive")
    regions = build_regions(args.workers)
    reports = []
    with ProcessPoolExecutor(max_workers=len(regions)) as executor:
        futures = [executor.submit(scan_region, args.bcftools, args.vcf, region, chrom,
                                   args.expected_samples, args.threads_per_region)
                   for region, chrom in regions]
        for future in as_completed(futures):
            reports.append(future.result())
    report = merge_reports(reports)
    report.update({"source_identity": source, "source_manifest": str(args.source_manifest),
                   "execution": {"engine": "indexed_multiprocess_bcftools", "regions": [region for region, _ in regions],
                                 "workers": len(regions), "threads_per_region": args.threads_per_region},
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
