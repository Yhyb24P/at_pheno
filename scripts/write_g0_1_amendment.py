"""Write the G0.1 confirmatory-design amendment and its hash-bound evidence.

This records only data/design metadata.  It neither reads phenotype values nor
fits a predictive model.  The original G0 v1 freeze is deliberately left
unchanged; this file is an additive amendment made before real prediction.
"""

import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).parents[1]


def digest(relative):
    path = ROOT / relative
    hasher = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(chunk)
    return {"path": relative, "sha256": hasher.hexdigest(), "bytes": path.stat().st_size}


def main():
    source = "data/processed/formal_genotype_v1/genotype_build_manifest.json"
    mirror = "data/manifests/formal_genotype_v1_build.json"
    shutil.copyfile(ROOT / source, ROOT / mirror)
    if (ROOT / source).read_bytes() != (ROOT / mirror).read_bytes():
        raise RuntimeError("tracked build-manifest mirror does not match processed source")
    bindings = [
        "data/manifests/confirmatory_v1/freeze_manifest.json",
        "data/manifests/1001g_hdf5_2026-09-22.json",
        mirror,
        "data/manifests/formal_vcf_preflight_v1.json",
        "data/manifests/formal_genotype_v1_qc.json",
        "data/manifests/confirmatory_v1/trait_panel.tsv",
        "data/manifests/confirmatory_v1/split_summary.json",
        "data/manifests/confirmatory_v1/nearclone_v1/nearclone_marker_panel.tsv",
        "data/manifests/confirmatory_v1/nearclone_v1/nearclone_sensitivity.json",
        "data/manifests/confirmatory_v1/nearclone_exact_refinement/nearclone_screen_d002_exact_pairs.tsv",
        "data/manifests/confirmatory_v1/nearclone_exact_refinement/nearclone_exact_refinement_summary.json",
        "data/manifests/confirmatory_v1/nearclone_exact_refinement/nearclone_blocks_exact_d001.tsv",
        "data/manifests/confirmatory_v1_1/baseline_spec_v1.json",
        "data/manifests/confirmatory_v1_1/common_feasibility.tsv",
        "data/manifests/confirmatory_v1_1/common_feasibility_summary.json",
        "data/manifests/confirmatory_v1/analysis_plan_v1.json",
    ]
    hdf5_manifest_path = ROOT / "data/manifests/1001g_hdf5_2026-09-22.json"
    hdf5_manifest = json.loads(hdf5_manifest_path.read_text(encoding="utf-8"))
    if not hdf5_manifest.get("sha256"):
        raise RuntimeError("HDF5 source manifest is missing its source SHA256")
    manifest = {
        "schema": "g0_formal_data_and_design_freeze_v1_1_amendment",
        "status": "G0_1_CONFIRMATORY_DESIGN_AMENDMENT_COMPLETE",
        "amends": {"tag": "g0-formal-freeze-v1", "purpose": "pre-G1 corrections; original G0 evidence is not rewritten"},
        "formal_prediction_results_seen_before_amendment": False,
        "scientific_corrections": {
            "practical_benefit": "Only claim practical meaningfulness when the conditional 95% block-bootstrap lower confidence bound of Q2(COMMON_ALL)-Q2(COMMON_250K) exceeds SESOI=0.01.",
            "baseline_protocol": "The complete additive-kernel specification is frozen in baseline_spec_v1.json.",
            "primary_feasibility": "Every frozen DTF1 blocked outer and inner fit partition has at least 250,000 fit-partition COMMON markers.",
            "nearclone_interpretation": "All four resolved phenotype cohorts contain no multi-accession d<=0.001 exact near-clone block; primary CV is near-clone-safe/group-balanced, not a relatedness-shift test.",
            "refinement_scope": "Full-HDF5 exact distances refine only d<=0.002 250K-screen nominated candidate pairs, not an exhaustive all-pairs calculation."
        },
        "artifact_bindings": {Path(item).name: digest(item) for item in bindings},
        "nearclone_provenance_binding": {
            "hdf5_source": {
                "path": hdf5_manifest["input"], "sha256": hdf5_manifest["sha256"],
                "bytes": hdf5_manifest["bytes"],
            },
            "hdf5_source_manifest": digest("data/manifests/1001g_hdf5_2026-09-22.json"),
            "coordinate_panel_logical_sha256": digest("data/manifests/confirmatory_v1/nearclone_v1/nearclone_marker_panel.tsv")["sha256"],
            "sensitivity": digest("data/manifests/confirmatory_v1/nearclone_v1/nearclone_sensitivity.json"),
            "screened_exact_pairs": digest("data/manifests/confirmatory_v1/nearclone_exact_refinement/nearclone_screen_d002_exact_pairs.tsv"),
        },
        "build_manifest_mirror": {
            "processed_source": digest(source),
            "tracked_mirror": digest(mirror),
            "byte_identical": True,
        },
    }
    out = ROOT / "data/manifests/confirmatory_v1_1/freeze_amendment_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    report = ROOT / "docs/G0_1_确认性设计修订.md"
    report.write_text(
        "# G0.1 确认性设计修订\n\n"
        "本修订发生在任何真实性状预测之前，不改变 `g0-formal-freeze-v1` 的历史内容。\n\n"
        "- 实用增益仅在 `LCB95[Q²(COMMON_ALL) − Q²(COMMON_250K)] > 0.01` 时声明。\n"
        "- 加性核、调参、数值稳定和失败策略固定于 `baseline_spec_v1.json`。\n"
        "- 237 个冻结训练分区的 COMMON 可行性已作仅基因型 GPU 审计；DTF1 blocked 的 25 个 outer/inner 分区均满足 250K。\n"
        "- 四个正式 cohort 均没有多成员精确 near-clone block，因此 primary CV 是 near-clone-safe/group-balanced，不能解释为 relatedness-shift 泛化。\n"
        "- full-HDF5 距离仅对 250K screen 提名 pair 作 exact refinement；不得误写为 exhaustive all-pairs。\n\n"
        "完整哈希绑定见 `data/manifests/confirmatory_v1_1/freeze_amendment_manifest.json`。\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
