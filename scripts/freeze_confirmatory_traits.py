"""Freeze G0 trait roles from the SHA256-audited AraPheno resolution registry.

The script never emits phenotype values.  It reads them only to prove each
selected published accession-level target has finite, unique accessions, and
binds that fact to the registry's recorded raw-value SHA256.
"""

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


ROLES = (("PRIMARY", "703", "DTF1"),
         ("INDEPENDENT_REPLICATION", "748", "stomata density"),
         ("SECONDARY_SENSITIVITY", "704", "RL"),
         ("SECONDARY_SENSITIVITY", "705", "CL"),
         ("DEFERRED_DUPLICATE_RESOLUTION", "390", "rosetteDM"),
         ("DEFERRED_DUPLICATE_RESOLUTION", "261", "FT10"),
         ("DEFERRED_DUPLICATE_RESOLUTION", "262", "FT16"))
USABLE = "PUBLISHED_ACCESSION_VALUE_USABLE"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_registry(path):
    with Path(path).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    result = {row["trait_id"]: row for row in rows}
    if len(result) != len(rows):
        raise ValueError("trait-resolution registry duplicates trait_id")
    return result


def raw_summary(path):
    rows = json.loads(Path(path).read_text())
    finite, seen = [], set()
    duplicates = 0
    for row in rows:
        accession = str(row["accession_id"])
        value = row.get("phenotype_value")
        if value is None or not math.isfinite(float(value)):
            continue
        if accession in seen:
            duplicates += 1
        seen.add(accession)
        finite.append(accession)
    return {"finite_row_n": len(finite), "finite_unique_accession_n": len(seen),
            "duplicate_finite_accession_n": duplicates}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--raw-values-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    registry = read_registry(args.registry)
    output = []
    for role, trait_id, name in ROLES:
        if trait_id not in registry:
            raise ValueError(f"trait {trait_id} absent from registry")
        source = args.raw_values_dir / f"{trait_id}.json"
        if not source.is_file():
            raise ValueError(f"trait {trait_id} raw values absent: {source}")
        row = registry[trait_id]
        raw_hash = sha256(source)
        summary = raw_summary(source)
        usable = role != "DEFERRED_DUPLICATE_RESOLUTION"
        if usable:
            if row["target_status"] != USABLE:
                raise ValueError(f"usable trait {trait_id} is not registry-resolved")
            if summary["finite_row_n"] != int(row["finite_panel_n"]) or summary["duplicate_finite_accession_n"]:
                raise ValueError(f"usable trait {trait_id} violates finite/unique registry semantics")
        if row["trait_name"] != name:
            raise ValueError(f"trait {trait_id} registry name differs from frozen role name")
        output.append({"role": role, "trait_id": trait_id, "trait_name": name,
                       "target_status": row["target_status"], "finite_panel_n": row["finite_panel_n"],
                       "raw_finite_unique_accession_n": summary["finite_unique_accession_n"],
                       "raw_values_path": str(source), "raw_values_file_sha256": raw_hash,
                       "registry_api_response_sha256": row["values_sha256"],
                       "hash_domain_relation": "distinct: local JSON reserialization vs historical API response bytes",
                       "local_aggregation": "none", "local_transformation": "none",
                       "target_definition": "AraPheno published accession-level raw numeric value",
                       "source_evidence_status": "registry-bound"})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    print(json.dumps({"status": "CONFIRMATORY_TRAIT_PANEL_FROZEN", "usable": [row["trait_id"] for row in output if row["role"] != "DEFERRED_DUPLICATE_RESOLUTION"]}))


if __name__ == "__main__":
    main()
