"""Summarize downloaded API rows without inventing genotype intersections."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    catalog = json.loads((args.snapshot/"phenotype_list.json").read_text())
    result = {"catalog_phenotypes": len(catalog),
              "catalog_studies_by_title": len({r["study"] for r in catalog}),
              "warning": "API row counts are not verified genotype intersections or independent endpoints.",
              "sources": json.loads((args.snapshot/"manifest.json").read_text()), "candidates": []}
    for path in sorted(args.snapshot.glob("values_*.json")):
        phenotype = path.stem.split("_", 1)[1]
        values = json.loads(path.read_text())
        meta = json.loads((args.snapshot/f"phenotype_{phenotype}.json").read_text())
        ids = Counter(str(r["accession_id"]) for r in values)
        countries = Counter(r.get("accession_country") or "unknown" for r in values)
        missing = sum(r["phenotype_value"] is None or not math.isfinite(float(r["phenotype_value"])) for r in values)
        result["candidates"].append({"phenotype_id": int(phenotype), "name": meta["name"],
            "doi": meta["doi"], "study": meta["study"], "unit": meta.get("uo_name"),
            "value_rows": len(values), "unique_accessions": len(ids), "nonfinite_or_missing": missing,
            "duplicate_accessions": {k: v for k, v in ids.items() if v > 1},
            "top_country_row_counts": countries.most_common(5),
            "genotype_intersection_n": None, "selection_status": "interface_candidate_not_frozen"})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("x") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({k: v for k, v in result.items() if k != "sources"}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
