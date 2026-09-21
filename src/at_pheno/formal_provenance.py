"""Formal-mode provenance v2 gate: recompute-not-trust bindings.

`PROVENANCE_SCHEMA_V2` names the v2 binding.  v1 records (no
`provenance_schema` key, `input_sha256`-style files) remain readable
for the pilot path, which never validates provenance; the formal gate
(`gate_v2`) rejects them on the missing `provenance_schema` field and
accepts only `at_pheno_formal_v2` records, re-verifying every file
binding against the files on disk at this run.  A recorded hash being
a 64-hex string is explicitly *not* a verification.
"""

import hashlib
import json
from pathlib import Path

PROVENANCE_SCHEMA_V2 = "at_pheno_formal_v2"

# The twelve spec-recorded bindings, in gate order: six SHA256 bindings
# that are recomputed against files on disk, the official source MD5
# (record-only, gate checks 32-hex format), the filter/multiallelic
# policy states, the reference assembly, the allele encoding, and the
# conversion script + commit pair.
SHA256_BINDINGS = (
    ("genotypes_sha256", "genotypes.npy"),
    ("samples_sha256", "samples.tsv"),
    ("variants_sha256", "variants.tsv"),
    ("phenotypes_sha256", "phenotypes.tsv"),
    ("registry_sha256", "the trait-resolution registry TSV"),
    ("source_vcf_sha256", "the source VCF"),
)

REQUIRED_FIELDS = ("provenance_schema", "dataset_id", "genotypes_sha256",
                   "samples_sha256", "variants_sha256", "phenotypes_sha256",
                   "registry_sha256", "source_vcf_sha256", "official_source_md5",
                   "vcf_filter_policy", "multiallelic_policy", "reference_assembly",
                   "allele_encoding", "conversion_script", "conversion_commit",
                   "notes", "phenotype_registry")


def file_sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _ok_sha256(value):
    return (isinstance(value, str) and len(value) == 64
            and all(c in "0123456789abcdef" for c in value))


def _ok_md5(value):
    return (isinstance(value, str) and len(value) == 32
            and all(c in "0123456789abcdef" for c in value))


def bind_dataset_input(dataset_dir, record, source_vcf=None, registry=None):
    """Re-verify the record's six file bindings against the files on disk.

    Every binding requires the target file to exist, the recorded value
    to be a 64-character lowercase hex SHA256, and a re-computed digest
    equal to it.  All mismatches raise a ValueError naming the binding
    field (never a generic "invalid hash").  Called after the gate's
    required-field check; `source_vcf` and `registry` are the recompute
    targets for `source_vcf_sha256` / `registry_sha256` and are required.
    """
    targets = {
        "genotypes_sha256": (dataset_dir / "genotypes.npy", "genotypes.npy"),
        "samples_sha256": (dataset_dir / "samples.tsv", "samples.tsv"),
        "variants_sha256": (dataset_dir / "variants.tsv", "variants.tsv"),
        "phenotypes_sha256": (dataset_dir / "phenotypes.tsv", "phenotypes.tsv"),
        "registry_sha256": (registry, "the trait-resolution registry TSV"),
        "source_vcf_sha256": (source_vcf, "the source VCF"),
    }
    for field, (path, name) in targets.items():
        if path is None or not isinstance(path, Path):
            raise ValueError(
                f"Formal provenance v2: {field} cannot be recomputed; "
                f"the {name} path is required")
        if not path.is_file():
            raise ValueError(f"Formal provenance v2: {field} target {path} is not a file")
        recorded = record.get(field)
        if not _ok_sha256(recorded):
            raise ValueError(
                f"Formal provenance v2: {field} must be a 64-character "
                f"lowercase hex SHA256")
        actual = file_sha256(path)
        if recorded != actual:
            raise ValueError(
                f"Formal provenance v2: {field} does not bind {name} "
                f"(record {recorded}, recomputed {actual} from {path})")


def gate_v2(dataset_dir, trait, source_vcf=None, registry=None):
    """The formal-mode gate for the v2 record.

    Loads `dataset_dir/"provenance.json"`, requires the record to be
    `PROVENANCE_SCHEMA_V2` (a v1 record is rejected here), validates
    every record field, re-verifies all file bindings through
    `bind_dataset_input` (recompute-not-trust), and binds the
    `phenotype_registry` sub-record to `trait` with `resolved is True`.
    Returns the (provenance-path, record) tuple.
    """
    directory = Path(dataset_dir)
    path = directory / "provenance.json"
    if not path.is_file():
        raise ValueError("Formal mode requires data/provenance.json")
    try:
        record = json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise ValueError("Formal provenance is not valid JSON") from exc
    if not isinstance(record, dict):
        raise ValueError("Formal provenance v2: provenance.json must be a JSON object")
    missing = sorted(set(REQUIRED_FIELDS) - set(record))
    if missing:
        raise ValueError(f"Formal provenance v2 missing fields: {', '.join(missing)}")
    if record["provenance_schema"] != PROVENANCE_SCHEMA_V2:
        raise ValueError(
            f"Formal provenance v2: provenance_schema must equal "
            f"'{PROVENANCE_SCHEMA_V2}' (got {record['provenance_schema']!r}; "
            f"a v1 record is not accepted by the formal gate)")
    if not isinstance(record["dataset_id"], str) or not record["dataset_id"].strip():
        raise ValueError("Formal provenance v2: dataset_id must be a nonempty string")
    for field in ("vcf_filter_policy", "multiallelic_policy",
                  "conversion_script", "conversion_commit"):
        if not isinstance(record[field], str) or not record[field].strip():
            raise ValueError(
                f"Formal provenance v2: {field} must be a documented nonempty string")
    if (not isinstance(record["notes"], list)
            or any(not isinstance(note, str) for note in record["notes"])):
        raise ValueError("Formal provenance v2: notes must be a list of strings")
    if record["reference_assembly"] != "TAIR10":
        raise ValueError(
            f"Formal provenance v2: reference_assembly must be 'TAIR10' "
            f"(got {record['reference_assembly']!r})")
    if record["allele_encoding"] != "ALT_dosage_0_1_2":
        raise ValueError(
            f"Formal provenance v2: allele_encoding must be 'ALT_dosage_0_1_2' "
            f"(got {record['allele_encoding']!r})")
    if not _ok_md5(record["official_source_md5"]):
        raise ValueError(
            "Formal provenance v2: official_source_md5 must be a "
            "32-character lowercase hex MD5 (parity with the upstream "
            "official archive is outside the local recompute)")
    sub = record["phenotype_registry"]
    if (not isinstance(sub, dict) or sub.get("trait_id") != trait
            or sub.get("resolved") is not True):
        state = (f"{{'trait_id': {sub.get('trait_id')!r}, "
                 f"'resolved': {sub.get('resolved')!r}}}"
                 if isinstance(sub, dict) else repr(sub))
        raise ValueError(
            f"Formal provenance v2: phenotype_registry must record the "
            f"requested trait {trait!r} with resolved true (got {state})")
    bind_dataset_input(directory, record, source_vcf=source_vcf, registry=registry)
    return path, record
