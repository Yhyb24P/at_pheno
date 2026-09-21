"""Build a reproducible phenotype feasibility registry from the AraPheno API.

The registry describes raw API availability.  It does not resolve replicates,
declare phenotypes comparable, or select a confirmation panel.
"""

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import csv
import hashlib
import json
from pathlib import Path
import time
import urllib.error
import urllib.request


BASE = "https://arapheno.1001genomes.org/rest/phenotype/"


def fetch_json(url, retries=3):
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "at-pheno-research/0.1"})
            with urllib.request.urlopen(request, timeout=45) as response:
                payload = response.read()
            return json.loads(payload), hashlib.sha256(payload).hexdigest()
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            if attempt == retries-1:
                raise
            time.sleep(2**attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--panel", type=Path, required=True)
    parser.add_argument("--admixture", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, help="First N sorted phenotype IDs; for smoke tests")
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        raise ValueError("workers must be 1..8")
    if args.out.exists():
        raise FileExistsError(args.out)
    catalog = json.loads(args.catalog.read_text())
    catalog = sorted(catalog, key=lambda row: int(row["phenotype_id"]))
    if args.limit:
        catalog = catalog[:args.limit]
    panel_ids = {str(row["Accession_ID"]) for row in json.loads(args.panel.read_text())}
    with args.admixture.open() as stream:
        groups = {row["id"]: row["group"] for row in csv.DictReader(stream)}
    if set(groups) != panel_ids:
        raise ValueError("Panel and published admixture IDs differ")
    args.out.mkdir(parents=True)
    raw_dir = args.out/"raw_values"
    raw_dir.mkdir()

    def inspect(meta):
        trait = str(meta["phenotype_id"])
        values, digest = fetch_json(f"{BASE}{trait}/values.json")
        path = raw_dir/f"{trait}.json"
        path.write_text(json.dumps(values, ensure_ascii=False)+"\n")
        ids = [str(row["accession_id"]) for row in values]
        counts = Counter(ids)
        matched = set(ids) & panel_ids
        return {"trait_id": trait, "name": meta.get("name"), "study": meta.get("study"),
                "doi": meta.get("doi"), "unit": meta.get("uo_name"),
                "growth_conditions": meta.get("growth_conditions"), "catalog_num_values": meta.get("num_values"),
                "raw_row_n": len(values), "unique_accession_n": len(counts),
                "panel_intersection_n": len(matched),
                "duplicate_accession_n": sum(n > 1 for n in counts.values()),
                "max_rows_per_accession": max(counts.values(), default=0),
                "group_counts": json.dumps(dict(sorted(Counter(groups[s] for s in matched).items()))),
                "values_sha256": digest, "retrieval_status": "ok", "error": ""}

    records, failures = [], []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        pending = {executor.submit(inspect, meta): meta for meta in catalog}
        for index, future in enumerate(as_completed(pending), 1):
            meta = pending[future]
            try:
                records.append(future.result())
            except Exception as exc:
                failures.append({"trait_id": str(meta["phenotype_id"]), "name": meta.get("name"),
                                 "error": f"{type(exc).__name__}: {exc}"})
            print(f"Completed {index}/{len(catalog)}", flush=True)
    records.sort(key=lambda row: int(row["trait_id"]))
    fields = list(records[0]) if records else ["trait_id"]
    with (args.out/"registry.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    manifest = {"retrieved_at": datetime.now(timezone.utc).isoformat(), "catalog_path": str(args.catalog),
                "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
                "requested": len(catalog), "completed": len(records), "failed": len(failures),
                "panel_n": len(panel_ids), "published_group_source": str(args.admixture),
                "interpretation": "feasibility registry only; duplicate rows need source-grounded resolution"}
    (args.out/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    (args.out/"failures.json").write_text(json.dumps(failures, indent=2)+"\n")
    if failures:
        raise SystemExit("Some API records failed; see failures.json. Do not select a final panel yet.")


if __name__ == "__main__":
    main()
