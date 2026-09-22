"""Formal-mode provenance v2 gate: recompute-not-trust bindings.

`PROVENANCE_SCHEMA_V2` names the v2 binding.  v1 records (no
`provenance_schema` key, `input_sha256`-style files) remain readable
for the pilot path, which never validates provenance; the formal gate
(`gate_v2`) rejects them on the missing `provenance_schema` field and
accepts only `at_pheno_formal_v2` records, re-verifying every file
binding against the files on disk at this run.  A recorded hash being
a 64-hex string is explicitly *not* a verification.

Two binding upgrades (data-semantic phase, 2026-09-22):

* The trait-resolution registry TSV is the truth source for
  `phenotype_registry.resolved`: the gate parses the SHA256-bound
  registry and requires the requested trait to occur exactly once with
  `target_status == PUBLISHED_ACCESSION_VALUE_USABLE`; a record that
  claims `resolved: true` against a registry row in any other state is
  rejected.
* The source-manifest JSON (the VCF-side source audit) is bound by
  SHA256 and cross-checked: its `local_sha256` must equal the record's
  bound `source_vcf_sha256`, its `official_md5` must equal the
  record's `official_source_md5`, and `md5_matches_official` must be
  true.  The gate never re-hashes the source VCF, so official MD5
  parity is carried by the SHA256-bound manifest record, not
  recomputed against the production file.
"""

import hashlib
import json
from pathlib import Path

PROVENANCE_SCHEMA_V2 = "at_pheno_formal_v2"

# The only registry target_status the formal gate accepts as truth.
RESOLVABLE_TARGET_STATUS = "PUBLISHED_ACCESSION_VALUE_USABLE"

# The spec-recorded bindings, in gate order: six SHA256 bindings
# recomputed against files on disk (the four dataset files plus the
# SHA256-bound registry TSV and source-manifest JSON), plus the
# record-only `source_vcf_sha256` that the gate does NOT recompute:
# its identity is carried by the SHA256-bound source-manifest record
# whose `local_sha256`/`official_md5` cross-bind it, and the raw-file
# check happens at the dataset-build stage (`verify_source_vcf`, P0-1
# gate).  The remaining record fields are the filter/multiallelic
# policy states, the reference assembly, the allele encoding, and
# the conversion script + commit pair.
SHA256_BINDINGS = (
    ("genotypes_sha256", "genotypes.npy"),
    ("samples_sha256", "samples.tsv"),
    ("variants_sha256", "variants.tsv"),
    ("phenotypes_sha256", "phenotypes.tsv"),
    ("registry_sha256", "the trait-resolution registry TSV"),
    ("source_vcf_sha256", "the source VCF (record-only; not re-verified by the gate)"),
    ("source_manifest_sha256", "the source-manifest JSON"),
)

REQUIRED_FIELDS = ("provenance_schema", "dataset_id", "genotypes_sha256",
                   "samples_sha256", "variants_sha256", "phenotypes_sha256",
                   "registry_sha256", "source_vcf_sha256", "source_manifest_sha256",
                   "official_source_md5",
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


def bind_dataset_input(dataset_dir, record, registry=None, source_manifest=None):
    """Re-verify the record's SHA256 bindings against the files on disk.

    Every binding requires the target file to exist, the recorded value
    to be a 64-character lowercase hex SHA256, and a re-computed digest
    equal to it.  All mismatches raise a ValueError naming the binding
    field (never a generic "invalid hash").  Called after the gate's
    required-field check; `registry` and `source_manifest` are the
    recompute targets for `registry_sha256` / `source_manifest_sha256`
    and are required.  `source_vcf_sha256` is record-only in the
    experiment-run gate: the gate checks its format and cross-binds it
    to the SHA256-bound source-manifest record (`bind_source_manifest`),
    and never re-reads the raw source VCF; the raw-file check is the
    dataset-build's job (`verify_source_vcf`).
    """
    targets = {
        "genotypes_sha256": (dataset_dir / "genotypes.npy", "genotypes.npy"),
        "samples_sha256": (dataset_dir / "samples.tsv", "samples.tsv"),
        "variants_sha256": (dataset_dir / "variants.tsv", "variants.tsv"),
        "phenotypes_sha256": (dataset_dir / "phenotypes.tsv", "phenotypes.tsv"),
        "registry_sha256": (registry, "the trait-resolution registry TSV"),
        "source_manifest_sha256": (source_manifest, "the source-manifest JSON"),
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


def registry_truth(registry, trait):
    """Parse the SHA256-bound trait-resolution registry TSV into truth.

    Header-driven column lookup for `trait_id` and `target_status`
    (the committed registry is a wide TSV; a minimal two-column
    `trait_id <tab> target_status` fixture parses too).  The requested
    trait must occur exactly once, and that row's `target_status`
    must be `RESOLVABLE_TARGET_STATUS`.  Raises ValueError naming the
    offending field otherwise.
    """
    lines = [line for line in registry.read_text(encoding="utf-8").splitlines()
             if line.strip()]
    if not lines:
        raise ValueError(
            f"Formal provenance v2: trait-resolution registry {registry} is empty")
    header = lines[0].split("\t")
    if "trait_id" not in header or "target_status" not in header:
        raise ValueError(
            "Formal provenance v2: trait-resolution registry header must "
            "carry trait_id and target_status columns "
            f"(got {header!r} in {registry})")
    trait_col = header.index("trait_id")
    status_col = header.index("target_status")
    seen = {}
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) <= max(trait_col, status_col):
            raise ValueError(
                f"Formal provenance v2: malformed registry row "
                f"{line!r} in {registry}")
        key = cells[trait_col]
        if key in seen:
            raise ValueError(
                f"Formal provenance v2: trait-resolution registry {registry} "
                f"lists trait {key!r} more than once; each trait must "
                f"appear exactly once")
        seen[key] = cells[status_col]
    if trait not in seen:
        raise ValueError(
            f"Formal provenance v2: phenotype_registry targets trait "
            f"{trait!r} which is not present in the trait-resolution "
            f"registry {registry} (registry truth is read from the "
            f"SHA256-bound registry, not from the record)")
    status = seen[trait]
    if status != RESOLVABLE_TARGET_STATUS:
        raise ValueError(
            f"Formal provenance v2: target_status for trait {trait!r} "
            f"must be {RESOLVABLE_TARGET_STATUS!r}; the registry {registry} "
            f"records {status!r}, so a record claiming "
            f"phenotype_registry.resolved true is rejected")
    return status


def bind_source_manifest(record, source_manifest):
    """Cross-bind the SHA256-bound source-manifest JSON to the record.

    The manifest (the VCF-side source audit) must agree with the
    record: its `local_sha256` equals the record's bound
    `source_vcf_sha256`, its `official_md5` equals the record's
    `official_source_md5`, and `md5_matches_official` is true.  The
    gate does NOT re-hash the source VCF; official-source MD5
    parity is carried by this record, not recomputed.
    """
    try:
        manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Formal provenance v2: source-manifest {source_manifest} "
            f"is not valid JSON") from exc
    if manifest.get("local_sha256") != record["source_vcf_sha256"]:
        raise ValueError(
            "Formal provenance v2: official source provenance binding: "
            f"source-manifest {source_manifest} records local_sha256 "
            f"{manifest.get('local_sha256')!r} that does not equal the "
            f"bound source_vcf_sha256 {record['source_vcf_sha256']!r}")
    if manifest.get("official_md5") != record["official_source_md5"]:
        raise ValueError(
            "Formal provenance v2: official source MD5: source-manifest "
            f"{source_manifest} records official_md5 "
            f"{manifest.get('official_md5')!r} that does not equal the "
            f"recorded official_source_md5 {record['official_source_md5']!r}")
    if manifest.get("md5_matches_official") is not True:
        raise ValueError(
            "Formal provenance v2: md5_matches_official: the source-manifest "
            f"{source_manifest} must record md5_matches_official true "
            f"(got {manifest.get('md5_matches_official')!r}); the gate does "
            "not re-hash the source VCF, so official-source MD5 parity "
            "is carried by the SHA256-bound manifest record only")


def verify_source_vcf(source_manifest, vcf_path):
    """The dataset-build (P0-1 gate) raw-file identity checker.

    The complement of `bind_source_manifest`: the experiment-run gate
    only cross-binds the record to a *SHA256-bound* source-manifest
    record without re-reading the raw VCF, so this is where the raw
    bytes are actually opened and the manifest's claims are proven
    against the file on disk:

    * the SHA256 of the on-disk file (digested byte-for-byte as stored:
      a .vcf.gz artifact is digested as the compressed artifact, not
      as an unpacked stream it was never unpacked into) equals
      `manifest["local_sha256"]`;
    * its MD5 equals `manifest["official_md5"]`;
    * `manifest["md5_matches_official"]` must be true, so file ==
      official public release digest is cross-bound through the
      manifest + record.

    The raw-file hashing happens here, and only here: a formal run
    must not re-hash the 19 GB source VCF.  Every failure raises a
    ValueError naming the offending field (never a generic
    "invalid hash").
    """
    manifest_path = Path(source_manifest)
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Formal provenance v2: source-manifest {manifest_path} "
            f"is not valid JSON") from exc
    target = Path(vcf_path)
    if not target.is_file():
        raise ValueError(f"Formal provenance v2: source VCF {target} does not exist")
    raw_sha256 = hashlib.sha256()
    raw_md5 = hashlib.md5()
    with target.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b""):
            raw_sha256.update(chunk)
            raw_md5.update(chunk)
    recorded_sha256 = manifest.get("local_sha256")
    recorded_md5 = manifest.get("official_md5")
    if not _ok_sha256(recorded_sha256):
        raise ValueError(
            f"Formal provenance v2: source-manifest {manifest_path} must record "
            f"local_sha256 as a 64-character lowercase hex SHA256 "
            f"(got {recorded_sha256!r})")
    if not _ok_md5(recorded_md5):
        raise ValueError(
            f"Formal provenance v2: source-manifest {manifest_path} must record "
            f"official_md5 as a 32-character lowercase hex MD5 "
            f"(got {recorded_md5!r})")
    if raw_sha256.hexdigest() != recorded_sha256:
        raise ValueError(
            f"Formal provenance v2: source-manifest {manifest_path} records "
            f"local_sha256 {recorded_sha256!r} that does not equal the "
            f"recomputed SHA256 {raw_sha256.hexdigest()!r} of {target}; "
            "the on-disk file is not the audited source VCF")
    if raw_md5.hexdigest() != recorded_md5:
        raise ValueError(
            f"Formal provenance v2: source-manifest {manifest_path} records "
            f"official_md5 {recorded_md5!r} that does not equal the "
            f"recomputed MD5 {raw_md5.hexdigest()!r} of {target}; the file "
            "is not byte-identical to the official release and the "
            "md5_matches_official claim is not reproducible")
    if manifest.get("md5_matches_official") is not True:
        raise ValueError(
            f"Formal provenance v2: source-manifest {manifest_path} must "
            f"record md5_matches_official true so that the official "
            f"release digest cross-binds file and record "
            f"(got {manifest.get('md5_matches_official')!r})")
    return {"sha256": raw_sha256.hexdigest(), "md5": raw_md5.hexdigest()}


def gate_v2(dataset_dir, trait, registry=None, source_manifest=None):
    """The formal-mode gate for the v2 record.

    Loads `dataset_dir/"provenance.json"`, requires the record to be
    `PROVENANCE_SCHEMA_V2` (a v1 record is rejected here), validates
    every record field, re-verifies all SHA256 bindings through
    `bind_dataset_input` (recompute-not-trust), then closes the two
    semantic bindings:

    * `registry_truth` — the requested trait must occur exactly once in
      the SHA256-bound registry TSV with `target_status` equal to
      `PUBLISHED_ACCESSION_VALUE_USABLE`; the record's
      `phenotype_registry.resolved true` is acceptable only when the
      registry truth agrees, and is rejected for every other state.
    * `bind_source_manifest` — the SHA256-bound source-manifest JSON
      must record `local_sha256 == source_vcf_sha256`,
      `official_md5 == official_source_md5` and `md5_matches_official
      true`, so the official-source MD5 claim is carried by a verified
      record rather than re-hashed by the gate.

    Two-layer provenance contract (the dataset-build concept, 2026-09):
    the gate binds the small immutable `source-manifest JSON` and never
    re-reads the 19 GB raw VCF.  Its recorded `local_sha256` /
    `official_md5` cross-bind the record's VCF identity fields, so the
    strength of source-identity evidence is unchanged, only the point
    of verification shifts: `verify_source_vcf` performs the raw-file
    check at the dataset-build stage.

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
            "32-character lowercase hex MD5")
    sub = record["phenotype_registry"]
    if (not isinstance(sub, dict) or sub.get("trait_id") != trait
            or sub.get("resolved") is not True):
        state = (f"{{'trait_id': {sub.get('trait_id')!r}, "
                 f"'resolved': {sub.get('resolved')!r}}}"
                 if isinstance(sub, dict) else repr(sub))
        raise ValueError(
            f"Formal provenance v2: phenotype_registry must record the "
            f"requested trait {trait!r} with resolved true (got {state})")
    # Record-only VCF identity: the experiment-run gate never re-reads
    # the raw VCF.  The value must be a well-formed SHA256, and
    # bind_source_manifest below cross-checks it against the
    # SHA256-bound source manifest; the raw-file check is the
    # dataset-build responsibility (verify_source_vcf).
    if not _ok_sha256(record["source_vcf_sha256"]):
        raise ValueError(
            "Formal provenance v2: source_vcf_sha256 (record-only cross-bound "
            "to the source manifest) must be a 64-character lowercase hex "
            "SHA256")
    bind_dataset_input(directory, record, registry=registry,
                       source_manifest=source_manifest)
    # Semantic closures, both against freshly re-verified on-disk files:
    registry_truth(registry, trait)
    bind_source_manifest(record, source_manifest)
    return path, record
