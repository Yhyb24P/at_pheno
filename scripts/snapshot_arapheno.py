"""Download small public metadata/candidate snapshots, with retrieval provenance.

No genotype download, trait selection or replicate aggregation is performed.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--phenotype", type=int, nargs="*", default=[261, 262])
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    base = "https://arapheno.1001genomes.org/rest/"
    resources = {"phenotype_list.json": base+"phenotype/list.json"}
    for phenotype in args.phenotype:
        resources[f"phenotype_{phenotype}.json"] = base+f"phenotype/{phenotype}.json"
        resources[f"values_{phenotype}.json"] = base+f"phenotype/{phenotype}/values.json"
    records = []
    for filename, url in resources.items():
        record = {"url": url, "filename": filename, "retrieved_at": datetime.now(timezone.utc).isoformat()}
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "at-pheno-research/0.1"})
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = response.read(16*1024*1024+1)
                if len(payload) > 16*1024*1024:
                    raise ValueError("Response exceeds metadata download limit")
                json.loads(payload)
                record["resolved_url"] = response.url
            (args.out/filename).write_bytes(payload)
            record.update(status="ok", bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
        except Exception as exc:
            record.update(status="failed", error=f"{type(exc).__name__}: {exc}")
        records.append(record)
        print(f"{record['status']}: {url}", flush=True)
    (args.out/"manifest.json").write_text(json.dumps(records, indent=2)+"\n")
    if any(r["status"] != "ok" for r in records):
        raise SystemExit("Some resources failed; see manifest.json. Do not treat snapshot as complete.")


if __name__ == "__main__":
    main()
