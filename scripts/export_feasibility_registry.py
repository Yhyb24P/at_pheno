"""Export the compact, versionable AraPheno feasibility registry."""

import argparse
import csv
from pathlib import Path


KEEP = ["trait_id", "name", "study", "doi", "unit", "catalog_num_values", "raw_row_n",
        "unique_accession_n", "panel_intersection_n", "duplicate_accession_n",
        "max_rows_per_accession", "group_counts", "values_sha256", "retrieval_status"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    with args.registry.open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if any(row["retrieval_status"] != "ok" for row in rows):
        raise ValueError("Registry contains failed retrievals")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=KEEP, delimiter="\t")
        writer.writeheader()
        writer.writerows({field: row[field] for field in KEEP} for row in rows)


if __name__ == "__main__":
    main()
