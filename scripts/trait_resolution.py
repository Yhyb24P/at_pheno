"""Validation helpers for the trait resolution registry (data-semantics stage).

The registry maps each candidate phenotype to its source, value level,
aggregation/transformation choices and a four-valued target status.  These
helpers hold the machine-checkable rules; the evidence rows themselves live
in data/manifests/trait_resolution_v1.tsv and are provenance (data), not code.
"""

from pathlib import Path

COLUMNS = [
    "trait_id", "trait_name", "study", "source_publication_doi", "arapheno_doi",
    "measurement_definition", "environment", "experimental_unit", "replicate_design",
    "database_value_level", "upstream_aggregation", "local_aggregation",
    "transformation_of_selected_target", "finite_panel_n", "duplicate_accession_n",
    "duplicate_resolution", "target_status", "candidate_role", "values_sha256",
    "primary_evidence", "notes",
]

TARGET_STATUSES = frozenset({
    "PUBLISHED_ACCESSION_VALUE_USABLE",
    "REQUIRES_DUPLICATE_RESOLUTION",
    "UNRESOLVED_SOURCE",
    "EXCLUDED_FOR_SCOPE",
})

_HEX64 = lambda s: isinstance(s, str) and len(s) == 64 and all(c in "0123456789abcdef" for c in s.lower())


def load_tsv(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    if header != COLUMNS:
        raise ValueError(f"Resolution registry header != required {len(COLUMNS)} columns: {header}")
    rows, seen = [], set()
    for n, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        if len(fields) != len(COLUMNS):
            raise ValueError(f"Line {n}: expected {len(COLUMNS)} fields, got {len(fields)}")
        row = dict(zip(COLUMNS, fields))
        if row["trait_id"] in seen:
            raise ValueError(f"Line {n}: duplicate trait_id {row['trait_id']}")
        seen.add(row["trait_id"])
        rows.append(row)
    return rows


def validate(rows):
    """Return a list of rule violations (empty means the registry is coherent)."""
    errors = []
    statuses_seen = set()
    for row in rows:
        tid = row["trait_id"]
        status = row["target_status"]
        if status not in TARGET_STATUSES:
            errors.append(f"{tid}: invalid target_status {status!r}")
            continue
        for counter in ("finite_panel_n", "duplicate_accession_n"):
            if not row[counter].isdigit():
                errors.append(f"{tid}: {counter} must be a non-negative integer")
        if status == "PUBLISHED_ACCESSION_VALUE_USABLE":
            # Preconditions recorded in the delivery spec: the target is the
            # published accession-level value and no conflicting duplicate
            # accession may exist in the database snapshot for this trait.
            if int(row["duplicate_accession_n"] or "0") > 0:
                errors.append(f"{tid}: USABLE with duplicate accessions "
                             f"({row['duplicate_accession_n']}) requires duplicate resolution first")
            if not _HEX64(row["values_sha256"]):
                errors.append(f"{tid}: USABLE requires a 64-hex values_sha256")
        if int(row["duplicate_accession_n"] or "0") > 0 and not row["duplicate_resolution"].strip():
            errors.append(f"{tid}: duplicate accessions without a documented duplicate_resolution")
        if row["target_status"] not in statuses_seen:
            statuses_seen.add(row["target_status"])
    return errors


def usable_rows(rows):
    return [r for r in rows if r["target_status"] == "PUBLISHED_ACCESSION_VALUE_USABLE"]


def duplicate_pending_rows(rows):
    return [r for r in rows if r["target_status"] == "REQUIRES_DUPLICATE_RESOLUTION"]
