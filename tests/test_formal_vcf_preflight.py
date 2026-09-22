import gzip
import importlib.util
import os
from pathlib import Path
import sys

import pytest


def module():
    path = Path(__file__).parents[1] / "scripts" / "formal_vcf_preflight.py"
    spec = importlib.util.spec_from_file_location("formal_vcf_preflight", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def parallel_module():
    scripts = Path(__file__).parents[1] / "scripts"
    path = scripts / "run_formal_vcf_preflight_parallel.py"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    spec = importlib.util.spec_from_file_location("run_formal_vcf_preflight_parallel", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def write_vcf(tmp_path, records):
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "tiny.vcf.gz"
    header = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"]
    with gzip.open(path, "wt") as out:
        out.write("\n".join(header + records) + "\n")
    return path


def test_preflight_dynamically_finds_gt_and_applies_filter_policy(tmp_path):
    mod = module()
    vcf = write_vcf(tmp_path, [
        "1\t1\t.\tA\tG\t.\tPASS\t.\tDP:GT:GQ\t10:0/0:20\t9:1|1:20",
        "1\t2\t.\tC\tT\t.\t.\t.\tGT\t./.\t0/1",
        "1\t3\t.\tG\tA\t.\tLowQual\t.\tGT\t0/.\t0",
        "1\t4\t.\tA\tC,G\t.\tPASS\t.\tGT\t0/0\t0/0",
    ])
    result = mod.scan(vcf, 2)
    assert result["status"] == "FORMAL_VCF_PREFLIGHT_PASS"
    assert result["accepted_structural_records"] == 2
    assert result["named_filter_excluded_records"] == 1
    assert result["filter_exact_counts"] == {".": 1, "LowQual": 1, "PASS": 2}
    assert result["format_layout_counts_accepted"] == {"DP:GT:GQ": 1, "GT": 1}
    assert result["gt_state_counts_accepted"] == {"fully_missing": 1, "heterozygous": 1,
        "phased": 1, "unphased": 2, "partial_missing": 0, "haploid": 0,
        "allele_index_gt_1": 0, "malformed": 0}


def test_preflight_blocks_unsupported_gt_or_format(tmp_path):
    mod = module()
    vcf = write_vcf(tmp_path, [
        "1\t1\t.\tA\tG\t.\tPASS\t.\tDP\t10\t9",
        "1\t2\t.\tA\tG\t.\tPASS\t.\tGT:GT\t0/0:0/0\t0/0:0/0",
        "1\t3\t.\tA\tG\t.\tPASS\t.\tGT\t0/.\t2/2",
        "1\t4\t.\tA\tG\t.\tPASS\t.\tGT\t0\t1",
    ])
    result = mod.scan(vcf, 2)
    assert result["status"] == "FILTER_POLICY_BLOCKED"
    assert result["hard_blockers"] == {"records_without_gt": 1, "records_duplicate_gt": 1,
        "sample_field_width_mismatch": 0, "partial_missing": 1, "haploid": 2,
        "allele_index_gt_1": 1, "malformed": 0}


def test_pigz_backend_is_optional_and_reports_nonzero_exit(tmp_path, monkeypatch):
    mod = module()
    vcf = write_vcf(tmp_path, ["1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1"])
    bindir = tmp_path / "bin"
    bindir.mkdir()
    pigz = bindir / "pigz"
    pigz.write_text("#!/bin/sh\ngzip -dc \"$4\"\n")
    pigz.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bindir}:{os.environ['PATH']}")
    assert mod.scan(vcf, 2, pigz_threads=2)["status"] == "FORMAL_VCF_PREFLIGHT_PASS"
    pigz.write_text("#!/bin/sh\nexit 7\n")
    pigz.chmod(0o755)
    with pytest.raises(ValueError, match="pigz exited 7"):
        mod.scan(vcf, 2, pigz_threads=2)


def test_parallel_merger_is_additive_and_rejects_sample_order_divergence(tmp_path):
    serial = module()
    parallel = parallel_module()
    reports = {}
    for chrom in parallel.CHROMS:
        vcf = write_vcf(tmp_path / chrom, [f"{chrom}\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1"])
        reports[chrom] = serial.scan(vcf, 2)
    merged = parallel.merge_reports(reports)
    assert merged["total_records"] == 5
    assert merged["accepted_structural_records"] == 5
    reports["5"]["sample_ids"] = ["S2", "S1"]
    with pytest.raises(ValueError, match="sample order"):
        parallel.merge_reports(reports)
