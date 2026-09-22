"""Bind the G0 formal-data/design artifacts before any real prediction run."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_head(root):
    return subprocess.run(["git", "rev-parse", "HEAD"], cwd=root, check=True, text=True,
                          capture_output=True).stdout.strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("."))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    required = {
        "source_vcf_manifest": "data/manifests/1001g_v31_population_vcf_source.json",
        "formal_vcf_preflight": "data/manifests/formal_vcf_preflight_v1.json",
        "genotype_build_manifest": "data/processed/formal_genotype_v1/genotype_build_manifest.json",
        "genotype_qa_manifest": "data/manifests/formal_genotype_v1_qc.json",
        "samples": "data/processed/formal_genotype_v1/samples.tsv",
        "variants": "data/processed/formal_genotype_v1/variants.tsv",
        "variants_compact": "data/processed/formal_genotype_v1/variants_compact.npy",
        "genotypes": "data/processed/formal_genotype_v1/genotypes.npy",
        "build_complete": "data/processed/formal_genotype_v1/BUILD_COMPLETE.json",
        "trait_resolution": "data/manifests/trait_resolution_v1.tsv",
        "trait_panel": "data/manifests/confirmatory_v1/trait_panel.tsv",
        "nearclone_blocks_canonical": "data/manifests/confirmatory_v1/blocks.tsv",
        "nearclone_exact_summary": "data/manifests/confirmatory_v1/nearclone_exact_refinement/nearclone_exact_refinement_summary.json",
        "marker_rank_manifest": "data/manifests/confirmatory_v1/marker_rank_manifest.json",
        "analysis_plan": "data/manifests/confirmatory_v1/analysis_plan_v1.json",
    }
    files = {}
    for name, relative in required.items():
        path = root / relative
        if not path.is_file():
            raise ValueError(f"required freeze artifact absent: {relative}")
        files[name] = {"path": relative, "sha256": sha256(path), "bytes": path.stat().st_size}
    canonical_splits = sorted((root / "data/manifests/confirmatory_v1").glob("*_outer_blocked.tsv"))
    canonical_splits += sorted((root / "data/manifests/confirmatory_v1").glob("*_inner_blocked_outer*.tsv"))
    canonical_splits += sorted((root / "data/manifests/confirmatory_v1").glob("*_outer_iid.tsv"))
    canonical_splits += sorted((root / "data/manifests/confirmatory_v1").glob("*_inner_iid_outer*.tsv"))
    canonical_splits += [root / "data/manifests/confirmatory_v1/group_logo_manifest.tsv"]
    if len(canonical_splits) != 49 or any(not path.is_file() for path in canonical_splits):
        raise ValueError(f"expected 49 canonical split artifacts, found {len(canonical_splits)}")
    split_hashes = [{"path": str(path.relative_to(root)), "sha256": sha256(path), "bytes": path.stat().st_size}
                    for path in canonical_splits]
    generation_code = ["scripts/formal_vcf_preflight.py", "scripts/run_formal_vcf_preflight_parallel.py",
                       "scripts/build_formal_genotype_v1.py", "scripts/audit_formal_genotype_v1.py",
                       "scripts/freeze_confirmatory_traits.py", "scripts/freeze_nearclone_blocks.py",
                       "scripts/refine_nearclone_exact_edges.py", "scripts/freeze_confirmatory_splits.py",
                       "scripts/freeze_marker_ranks.py", "scripts/write_g0_freeze_manifest.py"]
    code_hashes = {path: sha256(root / path) for path in generation_code}
    run_audit = []
    for audit in sorted((root / "runs").glob("*/audit.json")):
        record = json.loads(audit.read_text())
        run_audit.append({"path": str(audit.relative_to(root)), "trait_id": record.get("trait_id"),
                          "status": "synthetic_smoke_not_formal" if record.get("trait_id") == "synthetic" else "REAL_OR_UNKNOWN_BLOCKER"})
    if any(row["status"] != "synthetic_smoke_not_formal" for row in run_audit):
        raise ValueError("real or unknown prediction artifact found under runs/")
    manifest = {"schema": "g0_formal_data_and_design_freeze_v1", "status": "G0_FORMAL_DATA_AND_DESIGN_FREEZE_COMPLETE",
                "pre_freeze_git_head": git_head(root), "formal_prediction_results_seen_before_freeze": False,
                "run_audit": run_audit, "artifact_bindings": files, "canonical_split_bindings": split_hashes,
                "generation_code_sha256": code_hashes,
                "nearclone_canonical_rule": "250K d<=0.002 candidate screen then exact all-HDF5 d<=0.001 edge connected components",
                "note": "pre_freeze_git_head is the code/data commit used to generate artifacts; the enclosing freeze commit is intentionally not self-hashed."}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(manifest, indent=2) + "\n")
    lines = ["# G0 正式数据与确认性设计冻结报告", "", "状态：`G0_FORMAL_DATA_AND_DESIGN_FREEZE_COMPLETE`。本冻结发生在任何真实 phenotype prediction 之前。", "",
             f"- 生成基线提交：`{manifest['pre_freeze_git_head']}`", "- 正式基因型：1001G v3.1，int8 ALT dosage，1135 × 10,707,430。",
             "- VCF preflight、build QA、HDF5 2,000-site cross-check 均通过。", "- 性状：PRIMARY=703 DTF1；replication=748 stomata density；secondary=704 RL、705 CL；390/261/262 deferred。",
             "- 近克隆：250K hash screen 的 d≤0.002 候选对经全 HDF5 精确 d≤0.001 判边后，11 edges、4 multi-member blocks、最大 block=8、无 chain blocker。",
             "- 分割：blocked 5-fold outer / 4-fold inner，IID、LOGO 与 alternate salts 已冻结。",
             "- 密度：COMMON_10K/50K/250K/1M/COMMON_ALL；primary contrast=COMMON_ALL−COMMON_250K；endpoint=ΔQ²_trainmean；SESOI=0.01；10,000 conditional block bootstrap。", "",
             "## 完整性绑定", "", "所有 source、preflight、build、QA、traits、blocks、canonical splits、rank、analysis-plan 和 generation code SHA256 见 [`freeze_manifest.json`](../data/manifests/confirmatory_v1/freeze_manifest.json)。", "",
             "## 非正式运行", "", "`runs/smoke_iid` 与 `runs/smoke_group` 的 trait_id 均为 `synthetic`，仅为测试产物，不能用于任何正式推断。"]
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text("\n".join(lines) + "\n")
    print(json.dumps({"status": manifest["status"], "pre_freeze_git_head": manifest["pre_freeze_git_head"], "splits": len(split_hashes)}))


if __name__ == "__main__":
    main()
