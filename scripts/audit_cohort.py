"""Join public panel metadata, published admixture groups and raw phenotype rows.

Does not aggregate repeats or claim genotype-file verification without --hdf5.
"""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--panel", type=Path, required=True)
    p.add_argument("--admixture", type=Path, required=True)
    p.add_argument("--snapshots", type=Path, nargs="+", required=True)
    p.add_argument("--hdf5", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    panel = json.loads(args.panel.read_text())
    panel_ids = {str(r["Accession_ID"]) for r in panel}
    if len(panel_ids) != len(panel):
        raise ValueError("Duplicate official panel IDs")
    with args.admixture.open() as stream:
        admixture = list(csv.DictReader(stream))
    group = {r["id"]: r["group"] for r in admixture}
    if len(group) != len(admixture) or any(not g for g in group.values()):
        raise ValueError("Invalid admixture IDs/groups")
    if set(group) != panel_ids:
        raise ValueError("Official metadata and admixture ID sets differ")
    genotype_ids = None
    if args.hdf5:
        import h5py
        with h5py.File(args.hdf5, "r") as f:
            raw_ids = f["accessions"][:]
            strings = [s.decode() if isinstance(s, bytes) else str(s) for s in raw_ids]
            genotype_ids = set(strings)
            if len(genotype_ids) != len(strings):
                raise ValueError("Duplicate HDF5 accession IDs")
    effective = panel_ids if genotype_ids is None else panel_ids & genotype_ids
    records, cohorts, duplicates, source_paths = [], [], [], [args.panel, args.admixture]
    if args.hdf5:
        source_paths.append(args.hdf5)
    seen_traits = set()
    for folder in args.snapshots:
        for path in sorted(folder.glob("values_*.json")):
            trait = path.stem.split("_", 1)[1]
            if trait in seen_traits:
                raise ValueError("Duplicate trait snapshots; choose one version")
            seen_traits.add(trait)
            meta_path = folder/f"phenotype_{trait}.json"
            meta = json.loads(meta_path.read_text())
            rows = json.loads(path.read_text())
            source_paths.extend([path, meta_path])
            grouped = defaultdict(list)
            nonfinite = 0
            for row in rows:
                value = row["phenotype_value"]
                if value is None or not math.isfinite(float(value)):
                    nonfinite += 1
                grouped[str(row["accession_id"])].append(row)
            matched = set(grouped) & effective
            unambiguous = set()
            for accession, values in grouped.items():
                valid = [r for r in values if r["phenotype_value"] is not None and math.isfinite(float(r["phenotype_value"]))]
                if len(values) > 1:
                    duplicates.append({"trait_id": trait, "accession_id": accession,
                        "matched": accession in matched, "rows": len(values),
                        "values": [r["phenotype_value"] for r in values],
                        "obs_unit_ids": [r["obs_unit_id"] for r in values]})
                if accession in matched and len(values) == 1 and len(valid) == 1:
                    unambiguous.add(accession)
                if accession in matched:
                    cohorts.append({"trait_id": trait, "accession_id": accession,
                        "genetic_group": group[accession], "raw_row_count": len(values),
                        "status": "unique_finite" if accession in unambiguous else "needs_resolution"})
            records.append({"trait_id": trait, "name": meta["name"], "study": meta["study"],
                "doi": meta["doi"], "unit": meta.get("uo_name"), "growth_conditions": meta.get("growth_conditions"),
                "raw_rows": len(rows), "unique_ids": len(grouped), "nonfinite_rows": nonfinite,
                "official_panel_intersection_n": len(set(grouped)&panel_ids),
                "genotype_intersection_n": None if genotype_ids is None else len(set(grouped)&genotype_ids),
                "unambiguous_matched_n": len(unambiguous),
                "matched_duplicate_ids": sorted(a for a in matched if len(grouped[a]) > 1),
                "unmatched_ids": sorted(set(grouped)-effective),
                "group_counts_matched": dict(sorted(Counter(group[s] for s in matched).items())),
                "group_counts_unambiguous": dict(sorted(Counter(group[s] for s in unambiguous).items())),
                "selection_status": "candidate_not_confirmatory"})
    result = {"panel_n": len(panel_ids), "admixture_group_counts": dict(Counter(group.values())),
        "genotype_file_verified": genotype_ids is not None,
        "genotype_panel_id_difference": None if genotype_ids is None else
            {"panel_only": sorted(panel_ids-genotype_ids), "genotype_only": sorted(genotype_ids-panel_ids)},
        "input_sha256": {str(path): digest(path) for path in source_paths},
        "hdf5_path": None if args.hdf5 is None else str(args.hdf5),
        "script_sha256": digest(Path(__file__)),
        "traits": records}
    (args.out/"audit.json").write_text(json.dumps(result, indent=2, ensure_ascii=False)+"\n")
    (args.out/"duplicates.json").write_text(json.dumps(duplicates, indent=2, ensure_ascii=False)+"\n")
    with (args.out/"cohort_membership.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(cohorts[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(cohorts, key=lambda r: (r["trait_id"], r["accession_id"])))
    print(json.dumps([{k: r[k] for k in ("trait_id", "name", "raw_rows", "unique_ids", "official_panel_intersection_n", "genotype_intersection_n", "unambiguous_matched_n")} for r in records], indent=2))


if __name__ == "__main__":
    main()
