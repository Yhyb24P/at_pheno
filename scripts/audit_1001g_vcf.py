"""TASK02/TASK03 VCF-side audit of the 1001G v3.1 population VCF (at_pheno project).

Single-file, streaming, deterministic audit.  Pure Python: no numpy, no
network, no bcftools required.  What it does (spec sections in the Chinese
delivery doc tmp/at_pheno_Qwen3.8-27B_数据语义阶段执行交付包.md):

- 2.2 source calibration: official MD5 vs local MD5, local SHA256 (and the
  64-hex check computed, not assumed), URL/download-time provenance, and an
  optional `bcftools view -h` readability probe (only when a path is given).
- 2.3 header + sample audit: contig labels, FORMAT fields (with GT flag),
  filter names, sample IDs + count, and the name-key sample->official-1135
  accession-table mapping.  Three field-name keys are tried in fixed
  order (Accession_ID, Accession_Name, CS_Number); a key resolves the
  mapping when the stringified column values and the VCF sample names
  form the same set of UNIQUE keys (verbatim 1:1, no fuzzy match, no
  phonetic ecotype-name guessing).  File-order equality is reported as
  evidence only (order_match / order_divergences); when the mapping
  stays unresolved the VCF semantic gate fails.
- 2.4 GT-semantics small audit: for a bounded window (default: the first
  400 biallelic-SNP records on each of chrom 1-5, i.e. up to 2,000 records)
  of these biallelic-SNP records: ploidy histogram, phased/unphased,
  state histogram (REF-homo / HET / ALT-homo / missing / partial-missing /
  any-allele-index->1), and haploid single-token GT presence.  No
  interpretation of GT=1 as "HDF5 ALT" and no dosage conversion.
- 3.1/3.2 variant policy + catalog: the seven fixed-precedence buckets, a
  streamed 5-column biallelic-SNP catalog (tsv.gz), duplicate-key counts,
  and a small deterministic sample block so that only the summary JSONs and
  sample block stay under git while the gzipped catalog does not.
- 3.4 step 1 (VCF side only): if a precomputed picked-site list (JSON of
  [chrom, pos] pairs from the mapping script) is supplied, per-site state
  records with the spec-mandated keys.  No HDF5 write-back, no orientation
  written back to the HDF5 file, no assumption that HDF5 state 1 is an ALT
  dosage.

Deliberate limits: no full genotype matrix, no HDF5 binary access, no
bcftools re-decoding of the catalog, no plaintext unpacking of the VCF.
The HDF5-side join (spec 3.3/3.4 step 2) lives in
scripts/audit_hdf5_vcf_mapping.py, a separate agent output.

Importable pure helpers (all covered by tests under a tiny synthetic VCF):
`parse_vcf_line`, `classify_variant_line`, `summary_from_records`,
`gt_window_counts`, `map_samples`, `spot_site_report` /
`spot_site_not_found`, `read_catalog`, `catalog_duplicate_stats`,
`write_summary_outputs`, `run_audit`.
"""

import argparse
import gzip
import hashlib
import json
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

NUCLEAR = ("1", "2", "3", "4", "5")
NUCLEOTIDES = frozenset("ACGT")

# 3.1 policy buckets in fixed precedence; first hit wins and the buckets
# are disjoint.  Multi-allelic records are excluded *before* the N/indel
# checks and are never split (original multi-allelic records are NOT
# decomposed).  `nuclear_biallelic_snp` is the only bucket whose source
# data lines enter the catalog.
CATEGORIES = (
    "malformed", "non_nuclear", "multiallelic", "n_or_noncanonical",
    "indel_or_complex", "ref_equals_alt", "nuclear_biallelic_snp",
)
CAT_SNP = CATEGORIES[-1]

CATALOG_COLUMNS = ("chromosome", "position", "ref", "alt", "source_record_index")
SOURCE_JSON = "1001g_v31_population_vcf_source.json"
HEADER_JSON = "1001g_v31_vcf_header_sample_audit.json"
GT_JSON = "1001g_v31_gt_audit.json"
CATALOG_GZ = "1001g_v31_biallelic_snp_catalog.tsv.gz"
SUMMARY_JSON = "1001g_v31_snp_catalog_summary.json"
SPOT_CHECK_JSON = "1001g_v31_concordance_spot_check.json"

# 2.3 candidate field-name keys, tried in this exact order.  The table ships
# Accession_ID as int, so that column is stringified before the verbatim
# comparison; Accession_Name / CS_Number are compared as-is.
SAMPLE_MAP_KEYS = (
    ("accession_id", "Accession_ID"),
    ("accession_name", "Accession_Name"),
    ("cs_number", "CS_Number"),
)

SPOT_STATE_KEYS = ("homozygous_0", "heterozygous", "homozygous_1",
                   "missing", "out_of_range_or_unexpected")

EXAMPLE_CAP = 4          # bounded excluded-bucket examples per category
SAMPLE_BLOCK = 5        # first/last catalog rows kept in the summary JSON
GZIP_MAGIC = b"\x1f\x8b\x08"

DEFAULT_OUT_DIR = "data/manifests"
DEFAULT_ACC_TABLE = "data/raw/1001g_v3.1/1135_accessions_table.json"


# ---------------------------------------------------------------------------
# streaming helpers
# ---------------------------------------------------------------------------

def stream_kind(path) -> str:
    """'gzip' when the first 3 bytes are the gzip magic (1f 8b 08):
    both plain-gzip and .bgzip files carry it, so either is streamable
    through the stdlib `gzip` reader.  'plain' otherwise (VCF is plain
    text).  The encoder tool behind a 1001G release gz file is unknown:
    plain gzip or htslib-embedded bcftools either way, noted in the JSON
    output instead of being re-decided with bcftools."""
    with Path(path).open("rb") as handle:
        head = handle.read(3)
    return "gzip" if head[:3] == GZIP_MAGIC else "plain"


def open_vcf_stream(path):
    """Context manager yield text lines over the VCF file, whatever
    encoding is on disk (plain text / plain gz / bgzip as one logical
    gzip stream).  Decoding is utf-8 with replacement so any stray byte never
    crashes the pass."""
    encoding = stream_kind(path)
    if encoding == "gzip":
        return gzip.open(Path(path), "rt", encoding="utf-8", errors="replace")
    return open(Path(path), "r", encoding="utf-8", errors="replace")


def file_md5(path) -> str:
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_sha256(path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


# ---------------------------------------------------------------------------
# record parsing / 3.1 classification
# ---------------------------------------------------------------------------

def parse_vcf_line(line, sample_start=None):
    """Parse one physical VCF line.

    Header / blank lines return None.  Data lines with fewer than 7 columns
    (CHROM POS ID REF ALT QUAL FILTER) come back with the `malformed` flag
    (non-integer POS also yields pos None, which routes to the same
    bucket).

    `sample_start` is the 0-based index of the first per-sample column in
    this record.  The default (None -> 8) is correct for headers without
    the VCFv4.2+ FORMAT token in #CHROM.  When the header does carry it
    (the 1001G v3.1 layout: `#CHROM POS ... INFO FORMAT S1 S2 ...`), the
    caller passes 9 so the FORMAT column of the data rows (field index 8,
    e.g. `GT:GQ:DP`) is not misread as the first sample's values."""
    raw = line.rstrip("\n").rstrip("\r")
    if not raw or raw.startswith("#"):
        return None
    fields = raw.split("\t")
    if len(fields) < 7:
        return {"malformed": True, "chrom": None, "pos": None, "ref": None,
                "id": None, "alts": None, "filter": None,
                "sample_fields": None}
    try:
        pos = int(fields[1])
    except ValueError:
        pos = None
    start = 8 if sample_start is None else sample_start
    return {"malformed": False, "chrom": fields[0], "pos": pos,
            "id": fields[2], "ref": fields[3], "alts": fields[4].split(","),
            "info": fields[7] if len(fields) > 7 else ".",
            "filter": fields[6],
            "sample_fields": fields[start:] if len(fields) > start else []}


def classify_record(record) -> str:
    """The 3.1 policy bucket for a parsed record.  Precendence:
    malformed -> non_nuclear -> multiallelic -> n_or_noncanonical ->
    indel_or_complex -> ref_equals_alt -> nuclear_biallelic_snp."""
    if record.get("malformed") or record.get("pos") is None or \
            not record.get("ref") or not record.get("alts"):
        return CATEGORIES[0]
    if record["chrom"] not in NUCLEAR:
        return CATEGORIES[1]
    alts = record["alts"]
    if len(alts) > 1:
        return CATEGORIES[2]
    ref, alt = record["ref"], alts[0]
    if "N" in ref or "N" in alt or \
            any(ch not in NUCLEOTIDES for ch in (ref + alt)):
        return CATEGORIES[3]
    if len(ref) > 1 or len(alt) > 1:
        return CATEGORIES[4]
    if ref == alt:
        return CATEGORIES[5]
    return CAT_SNP


def classify_variant_line(line) -> str:
    """The 3.1 bucket for a physical VCF line (parse then classify).
    Header / blank lines land in `malformed` (they are skipped before
    counting, but classified as such if ever asked)."""
    record = parse_vcf_line(line)
    if record is None:
        return CATEGORIES[0]
    return classify_record(record)


def catalog_line(record, record_index) -> str:
    """The five 3.2 catalog columns, tab-joined in CATALOG_COLUMNS order.
    `source_record_index` is the 1-based data-line index in VCF file order;
    header lines are never counted."""
    return "\t".join((record["chrom"], str(record["pos"]), record["ref"],
                      record["alts"][0], str(record_index)))


# ---------------------------------------------------------------------------
# 3.1/3.2 running summary
# ---------------------------------------------------------------------------

def new_summary() -> dict:
    """Zero-initialized state for one full VCF pass."""
    return {
        "total_records": 0,
        "nuclear_records": 0,
        "per_chrom": {chrom: 0 for chrom in NUCLEAR},
        "type_counts": {cat: 0 for cat in CATEGORIES},
        "catalog_rows": 0,
        "excluded_examples": {cat: [] for cat in CATEGORIES[:-1]},
    }


def update_summary(summary, record) -> str | None:
    """Fold one parsed record (or None for header lines) into the running
    state; returns the category used (None for header).  O(1) per record;
    excluded-bucket examples stop at EXAMPLE_CAP each so state stays small
    through a 10M-line pass."""
    if record is None:
        return None
    summary["total_records"] += 1
    chrom = record.get("chrom")
    if chrom in NUCLEAR:
        summary["nuclear_records"] += 1
        summary["per_chrom"][chrom] += 1
    cat = classify_record(record)
    summary["type_counts"][cat] += 1
    if cat == CAT_SNP:
        summary["catalog_rows"] += 1
    else:
        bucket = summary["excluded_examples"].get(cat)
        if bucket is not None and len(bucket) < EXAMPLE_CAP:
            bucket.append({"chrom": chrom, "pos": record.get("pos"),
                           "ref": record.get("ref"), "alts": record.get("alts"),
                           "filter": record.get("filter")})
    return cat


def summary_from_records(records) -> dict:
    """Single-pass summary over an iterable of parsed records (None entries
    are tolerated, e.g. when header lines are interleaved).  Returns a
    completed `new_summary()` dict; duplicate-key counts and the sample
    block are attached by `write_summary_outputs`."""
    summary = new_summary()
    for record in records:
        update_summary(summary, record)
    return summary


# ---------------------------------------------------------------------------
# 3.2 catalog duplicate + sample block
# ---------------------------------------------------------------------------

def read_catalog(path):
    """Yield (chrom, pos_int, ref, alt, index_int) rows from the streamed
    catalog tsv.gz.  The header line is skipped and never counted.  The
    gzip is read back through `gzip` (never unpacked to a plaintext
    copy on disk)."""
    with gzip.open(Path(path), "rt", encoding="utf-8") as handle:
        header = handle.readline()
        if header.rstrip("\n") != "\t".join(CATALOG_COLUMNS):
            raise ValueError(f"catalog header mismatch: {header!r}")
        for row in handle:
            row = row.rstrip("\n")
            if not row:
                continue
            chrom, pos, ref, alt, index = row.split("\t")
            yield chrom, int(pos), ref, alt, int(index)


def catalog_duplicate_stats(path) -> dict:
    """Re-read the written catalog gz (bounded two-pass over the *catalog*,
    not over the VCF) and report:

    - `catalog_rows` total rows;
    - `duplicate_chrom_position_pairs` = the number of (chrom,pos) keys
      that hold 2 or more catalog rows (duplicates can be different
      ref/alt at the same coordinate);
    - `duplicate_chrom_position_ref_alt_quadruples` = the same for the
      full (chrom,pos,ref,alt) key;
    - `first_rows` / `last_rows`: the first and last SAMPLE_BLOCK rows
      (each as a 5-column list), the deterministic sample block.
    """
    pos_count, full_count = Counter(), Counter()
    first_rows, last_rows = [], []
    rows = 0
    for chrom, pos, ref, alt, index in read_catalog(path):
        rows += 1
        pos_count[(chrom, pos)] += 1
        full_count[(chrom, pos, ref, alt)] += 1
        if len(first_rows) < SAMPLE_BLOCK:
            first_rows.append([chrom, str(pos), ref, alt, str(index)])
        last_rows.append([chrom, str(pos), ref, alt, str(index)])
        if len(last_rows) > SAMPLE_BLOCK:
            last_rows.pop(0)
    return {"catalog_rows": rows,
            "duplicate_chrom_position_pairs":
                sum(v > 1 for v in pos_count.values()),
            "duplicate_chrom_position_ref_alt_quadruples":
                sum(v > 1 for v in full_count.values()),
            "first_rows": first_rows,
            "last_rows": last_rows}


# ---------------------------------------------------------------------------
# 2.3 sample mapping
# ---------------------------------------------------------------------------

def load_accession_table(path):
    """Rows of the official 1135 accession table (JSON: a list of row dicts
    with at least the `Accession_ID` field, or an object whose `rows` key
    holds such a list).  Returns None when missing/unreadable so that the
    mapping gate fails clearly instead of guessing."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, UnicodeDecodeError, ValueError):
        return None
    if isinstance(data, dict):
        data = data.get("rows")
    if not isinstance(data, list) or not all(isinstance(r, dict) for r in data):
        return None
    return data


def map_samples(sample_ids, table_rows) -> dict:
    """2.3 sample mapping by name-key identity.  The candidate columns
    are tried in fixed order (see SAMPLE_MAP_KEYS).  A key *resolves*
    the mapping when the VCF sample names and the stringified column
    values form the same set AND every name on both sides is unique:
    that is a verbatim 1:1 correspondence built from the file's own
    keys (e.g. `Accession_ID`), not an ecotype-name or fuzzy match.
    Delivery doc 2.3 requires a one-to-one correspondence, not a
    position-locked file order: the VCF #CHROM order may diverge from
    the table row order (the 1001G v3.1 VCF lists samples in numeric
    accession order; the accession table keeps a different row order),
    so `order_match` and `order_divergences` are *reported as evidence
    only* and do not gate.  A set mismatch (or a duplicated name
    inside the VCF, or twice in the table) leaves the mapping
    `unresolved` and fails the VCF semantic gate, with the diff sets
    printed for inspection.  No fuzzy / phonetic / ecotype-name
    guessing, per spec 2.3."""
    ids = [str(s) for s in sample_ids]
    note = ("resolved = at least one candidate key whose stringified column "
            "values and the VCF sample names form the same set of UNIQUE "
            "keys (verbatim 1:1 name-key correspondence).  File-order "
            "equality is reported (order_match / order_divergences), it "
            "does not gate.  No fuzzy / phonetic / ecotype-name matching "
            "(spec 2.3).")
    if table_rows is None:
        return {"candidates": [], "resolved_key": None,
                "status": "unresolved", "reason": "accession absent",
                "note": note}
    candidates = []
    resolved_key = None
    for key_name, field in SAMPLE_MAP_KEYS:
        declared = [str(row.get(field, "")) for row in table_rows]
        order_match = ids == declared
        set_match = set(ids) == set(declared)
        unique_keys = (len(set(ids)) == len(ids)
                       and len(set(declared)) == len(declared))
        resolves = set_match and unique_keys
        n_divergences = (0 if order_match else
                         sum(1 for a, b in zip(ids, declared) if a != b))
        candidates.append({
            "key": key_name, "field": field,
            "n_samples": len(ids),
            "n_table_entries": len(table_rows),
            "order_match": order_match,
            "order_divergences": n_divergences,
            "set_match": set_match,
            "unique_keys_on_both_sides": unique_keys,
            "resolves": resolves,
            "in_vcf_not_in_table":
                sorted(set(ids) - set(declared))[:EXAMPLE_CAP],
            "in_table_not_in_vcf":
                sorted(set(declared) - set(ids))[:EXAMPLE_CAP]})
        if resolves and resolved_key is None:
            resolved_key = key_name
    return {"candidates": candidates, "resolved_key": resolved_key,
            "status": "resolved" if resolved_key else "unresolved",
            "note": note}


# ---------------------------------------------------------------------------
# 2.4 GT small audit
# ---------------------------------------------------------------------------

_GT_STATE_KEYS = ("ref_homo", "heterozygous", "alt_homo", "missing",
                  "partial_missing", "allele_index_gt_1", "no_gt")


def _gt_counts(rows) -> dict:
    """Histo counters over lists of per-sample GT cells `rows`
    (each `cells` is one record's GT-column text).  Counts what the VCF
    says; it does *not* map a GT=1 to an HDF5 ALT, and dosage is not
    estimated anywhere."""
    counts = {
        "n_records": len(rows),
        "n_cells": 0,
        "ploidy": {"1": 0, "2": 0, ">2": 0},
        "phased_cells": 0,
        "unphased_cells": 0,
        "n_haploid_gt_cells": 0,
        "haploid_gt_present": False,
        "states": {key: 0 for key in _GT_STATE_KEYS},
    }
    for cells in rows:
        for text in cells:
            counts["n_cells"] += 1
            if not text:
                counts["states"]["no_gt"] += 1
                continue
            phased = "|" in text
            tokens = text.split("|" if phased else "/")
            n = len(tokens)
            counts["ploidy"]["1" if n == 1 else "2" if n == 2 else ">2"] += 1
            if phased:
                counts["phased_cells"] += 1
            else:
                counts["unphased_cells"] += 1
            if n == 1:
                counts["n_haploid_gt_cells"] += 1
                counts["haploid_gt_present"] = True
            present = [tok for tok in tokens if tok not in (".", "")]
            if not present:
                counts["states"]["missing"] += 1
                continue
            values, unparseable = [], False
            for tok in present:
                if tok.isdigit():
                    values.append(int(tok))
                else:
                    unparseable = True
            if unparseable or any(v > 1 for v in values):
                counts["states"]["allele_index_gt_1"] += 1
            elif len(present) < n:
                counts["states"]["partial_missing"] += 1
            elif all(v == 0 for v in values):
                counts["states"]["ref_homo"] += 1
            elif all(v == 1 for v in values):
                counts["states"]["alt_homo"] += 1
            else:
                counts["states"]["heterozygous"] += 1
    return counts


def gt_window_counts(per_chrom_rows) -> dict:
    """2.4 counters for a bounded biallelic-SNP record window.
    `per_chrom_rows` maps chrom label -> a list of per-record lists of
    GT-cell texts (first-annotation column only).  Returns
    {"per_chrom": {...}, "combined": {...}} where each blob carries the
    ploidy / phasing / state histograms and the haploid flag."""
    per_chrom = {chrom: _gt_counts(rows) for chrom, rows in per_chrom_rows.items()}
    combined_rows = []
    for chrom in per_chrom_rows:
        combined_rows.extend(per_chrom_rows[chrom])
    return {"per_chrom": per_chrom, "combined": _gt_counts(combined_rows)}


# ---------------------------------------------------------------------------
# 3.4 spot check, VCF-side step 1
# ---------------------------------------------------------------------------

def spot_state_code(gt_text) -> str:
    """VCF-side state code for one sample's GT column text.  No threshold
    guessing, no backwrite into the HDF5, no dosage assumption; the ploidy
    check belongs to the 2.4 window audit, not to 3.4:

      homozygous_0            all present alleles index 0 (0/0-like)
      heterozygous            mix of 0 and 1
      homozygous_1            all present alleles index 1 (1/1-like)
      missing                 no GT / all-missing / partially missing
      out_of_range_or_unexpected  allele index > 1 or unparseable token
    """
    if not gt_text:
        return "missing"
    phased = "|" in gt_text
    tokens = gt_text.split("|" if phased else "/")
    present = [tok for tok in tokens if tok not in (".", "")]
    if not present:
        return "missing"
    values = []
    for tok in present:
        if not tok.isdigit() or int(tok) > 1:
            return "out_of_range_or_unexpected"
        values.append(int(tok))
    if len({v for v in values}) == 2:
        return "heterozygous"
    return "homozygous_0" if values[0] == 0 else "homozygous_1"


def _spot_counts(gt_texts) -> dict:
    codes = [spot_state_code(text) for text in gt_texts]
    return {key: codes.count(key) for key in SPOT_STATE_KEYS}


def _spot_scores(state_counts) -> tuple:
    """(concordance_if_hdf5_0_is_REF, concordance_if_hdf5_1_is_REF,
    best_orientation, best_concordance) as fractions of `n_comparable`
    (comparable = homozygous and non-missing samples).  The VCF side can
    only record which orientation *the VCF itself* supports; the actual
    HDF5 join happens in the HDF5-side audit.  Ties are surfaced as
    `ambiguous`, not silently broken."""
    n_comparable = state_counts["homozygous_0"] + state_counts["homozygous_1"]
    if n_comparable == 0:
        return None, None, None, None
    h0 = state_counts["homozygous_0"] / n_comparable
    h1 = state_counts["homozygous_1"] / n_comparable
    if h0 > h1:
        return h0, h1, "hdf5_0_is_ref", h0
    if h1 > h0:
        return h0, h1, "hdf5_1_is_ref", h1
    return h0, h1, "ambiguous", h0


def spot_site_report(chrom, pos, gt_texts, per_sample=None) -> dict:
    """Per-site VCF-side spot-check entry (spec-mandated keys included).
    `gt_texts` are the GT-column cells for this record's samples, in
    VCF sample-column order.  `per_sample` -- one state code (the
    value of `spot_state_code` over that sample's GT cell) per VCF
    sample column, same order; when supplied it is passed through
    verbatim under the `per_sample` key so the HDF5-side 3.4 join can
    pair HDF5 accession column i against VCF sample i (the ordering
    agreement is documented on the `--sample-file` side; this key
    records the alignment, it never resolves orientation by itself).
    Omitted when not supplied."""
    state_counts = _spot_counts(gt_texts)
    h0, h1, best_o, best_c = _spot_scores(state_counts)
    entry = {"chrom": chrom, "pos": pos, "record_found": True,
            "concordance_if_hdf5_0_is_REF": h0,
            "concordance_if_hdf5_1_is_REF": h1,
            "best_orientation": best_o, "best_concordance": best_c,
            "n_comparable": state_counts["homozygous_0"] + state_counts["homozygous_1"],
            "n_heterozygous": state_counts["heterozygous"],
            "n_missing_vcf": state_counts["missing"],
            "n_states_out_of_gt": state_counts["out_of_range_or_unexpected"],
            "states": state_counts}
    if per_sample is not None:
        entry["per_sample"] = list(per_sample)
    return entry


def spot_site_not_found(chrom, pos) -> dict:
    """Entry for a picked site that does not match a catalog
    (biallelic-SNP) data line in this VCF: the spec-mandated keys are
    still present, with null scores."""
    zero = {key: 0 for key in SPOT_STATE_KEYS}
    return {"chrom": chrom, "pos": pos, "record_found": False,
            "concordance_if_hdf5_0_is_REF": None,
            "concordance_if_hdf5_1_is_REF": None,
            "best_orientation": None, "best_concordance": None,
            "n_comparable": 0, "n_heterozygous": 0, "n_missing_vcf": 0,
            "n_states_out_of_gt": 0, "states": zero,
            "note": "site not a biallelic-SNP data line of this VCF"}


def load_picked_sites(path):
    """Site-pool JSON from the mapping script: a list of
    `[chrom, pos]` pairs, an object with a top-level `sites` list, or
    the whole `hdf5_vcf_position_mapping_summary` object (pool read from
    `position_mapping.site_pool_sample`).  Used *as-is*, in file order
    -- no re-sorting or resampling here (deterministic selection already
    happened on the mapping side)."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if isinstance(data, dict):
        nested = (data.get("position_mapping") or {})
        data = (data.get("site_pool_sample") or data.get("sites")
                or nested.get("site_pool_sample"))
    if not isinstance(data, list):
        raise ValueError("picked-sites JSON must be a list of [chrom, pos] pairs")
    pairs = []
    for item in data:
        if len(item) != 2:
            raise ValueError(f"picked-sites entries must be [chrom, pos]: {item!r}")
        pairs.append((str(item[0]), int(item[1])))
    return pairs


# ---------------------------------------------------------------------------
# header / source facts
# ---------------------------------------------------------------------------

def _meta(line, tag):
    """(id, length-or-None) from a `##<tag>=<...>` line, if it starts
    with that exact tag; length only matters for contig.  Meta-line
    attributes are separated by ',' (the VCF meta syntax used by the
    1001G release files, `##FORMAT=<ID=GT,Number=1,...>`) and ';'
    (legacy tools); both are handled.  The ID value is taken up to the
    first separator so attribute text never leaks into it.  Without
    the comma split, `##FORMAT=<ID=GT,...>` was parsed as one
    attribute string and the GT-lookup (`"GT" in formats`) silently
    failed, collapsing the whole 2.4 window into the no_gt bucket."""
    prefix = f"##{tag}=<"
    if not line.startswith(prefix):
        return None, None
    body = line[len(prefix):]
    end = body.rfind(">")
    if end != -1:
        body = body[:end]
    _id = length = None
    for part in re.split(r"[;,]", body):
        if part.startswith("ID="):
            _id = part[3:]
        elif part.startswith("length=") and tag == "contig":
            length = part[7:]
    return (_id, length) if _id else (None, None)


def collect_vcf_header(vcf_path) -> dict:
    """Bounded header read: contig labels (name -> length, None when
    absent), FORMAT id list (first-seen order) + GT presence, filter
    ids, and the #CHROM sample header.  Stops at the first #CHROM line."""
    contigs = {}
    filters = []
    formats = []
    sample_ids = []
    format_column_in_chrom = False
    for line in open_vcf_stream(vcf_path):
        if line.startswith("##contig=<"):
            contig_id, contig_len = _meta(line, "contig")
            if contig_id not in contigs:
                contigs[contig_id] = int(contig_len) if (
                    contig_len is not None and contig_len.isdigit()
                ) else None
        elif line.startswith("##FILTER=<"):
            filter_id, _ = _meta(line, "FILTER")
            if filter_id not in filters:
                filters.append(filter_id)
        elif line.startswith("##FORMAT=<"):
            format_id, _ = _meta(line, "FORMAT")
            if format_id not in formats:
                formats.append(format_id)
        elif line.startswith("#CHROM"):
            fields = line.rstrip("\n").split("\t")
            # VCFv4.2+ #CHROM lines carry a literal FORMAT token between
            # INFO and the sample columns (the 1001G v3.1 file does:
            # 1144 columns = fixed 8 + FORMAT + 1135 samples).  When
            # present it must not count as a sample name; the sample
            # list starts one column later.
            format_column_in_chrom = len(fields) > 8 and fields[8] == "FORMAT"
            sample_ids = fields[9:] if format_column_in_chrom else fields[8:]
            break
    # The GT cell position is the assumed 0 (the sample's first field):
    # both released and tool-generated VCFs put GT first.  Declared in
    # the payload; if the GT is not the first sample field in this
    # file's records the missing cells pile up into the `no_gt` bucket
    # of the 2.4 audit, where the assumption fails visibly.
    gt_index = 0 if "GT" in formats else None
    return {"contig_labels": contigs,
            "format_fields": formats,
            "gt_present": gt_index is not None,
            "gt_column_index": gt_index,
            "gt_index_note": "the GT cell position is the assumed 0 (sample's "
                             "first field); not asserted column-by-column "
                             "per record (the VCF spec does not expose the "
                             "per-record sample-field ordering from the "
                             "header)",
            "format_column_in_chrom": format_column_in_chrom,
            "filter_names": filters,
            "n_samples": len(sample_ids),
            "sample_ids": sample_ids,
            "first_samples": sample_ids[:10]}


def probe_bcftools(bcftools, vcf_path) -> dict:
    """Spec 2.2 readability check, only when the caller passes a bcftools
    path: capture `--version` (first line) and whether `view -h <vcf>`
    exits 0.  Never executed when bcftools is None; an unreadable
    binary leaves nulls / False instead of crashing the audit."""
    def run(cmd):
        return subprocess.run(list(cmd), capture_output=True, text=True,
                              timeout=60)
    try:
        version_proc = run([str(bcftools), "--version"])
        lines = [ln for ln in version_proc.stdout.splitlines() if ln.strip()]
        version = lines[0] if lines else None
    except (OSError, subprocess.SubprocessError) as exc:
        return {"bcftools_path": str(bcftools), "bcftools_version": None,
                "view_header_exit_0": False, "probe_error": str(exc)}
    try:
        view_proc = run([str(bcftools), "view", "-h", str(vcf_path)])
        view_ok = view_proc.returncode == 0
        view_error = None
    except (OSError, subprocess.SubprocessError) as exc:
        view_ok, view_error = False, str(exc)
    return {"bcftools_path": str(bcftools), "bcftools_version": version,
            "view_header_exit_0": view_ok, "probe_error": view_error}


def compose_source_json(vcf_path, official_md5=None, source_url=None,
                       download_time=None, download_tool=None,
                       bcftools_probe=None, stream_encoding=None) -> dict:
    """The 2.2 calibration record.  MD5/SHA256 are computed *from the
    file* (never assumed); the 64-hex check is computed, and bcftools is
    always null unless a path was passed.  The note records that the
    catalog was streamed without bcftools and without unpacking the
    gzip to plaintext."""
    local_md5 = file_md5(vcf_path)
    local_sha = file_sha256(vcf_path)
    return {
        "series": "1001g_v31_population_vcf_source",
        "vcf_path": str(vcf_path),
        "file_bytes": Path(vcf_path).stat().st_size,
        "source_url": source_url,
        "download_time": download_time,
        "download_tool": download_tool,
        "official_md5": official_md5,
        "local_md5": local_md5,
        "md5_matches_official": (
            local_md5.lower() == official_md5.lower()
            if official_md5 else None
        ),
        "local_sha256": local_sha,
        "sha256_is_64_hex": (
            isinstance(local_sha, str) and len(local_sha) == 64 and
            all(c in "0123456789abcdef" for c in local_sha.lower())
        ),
        "stream_encoding": stream_encoding,
        "encoder_note": (
            "encoder tool family not identified; plain-gz and .bgzip both "
            "carry the gzip magic and stream through the stdlib gzip "
            "reader; no bcftools re-decode was performed"
        ),
        "catalog_streaming": (
            "pure-Python streaming reader; the catalog was never "
            "re-decoded via bcftools, and the VCF was never unpacked "
            "into a plaintext temp file"
        ),
        "bcftools_header_probe": bcftools_probe,
        "claims_not_supported_by_this_audit": [
            "HDF5 0/1 state semantics (HDF5-side script owns that side)",
            "GT dosage: no GT value is ever rescaled by ploidy",
        ],
    }


# ---------------------------------------------------------------------------
# main audit pass
# ---------------------------------------------------------------------------

def run_audit(vcf, out_dir, accession_table=DEFAULT_ACC_TABLE,
              official_md5=None, source_url=None, download_time=None,
              download_tool=None, bcftools=None, window_chroms=5,
              records_per_chrom=400, picked_sites=None, save_catalog=True) -> dict:
    """One streaming pass over the VCF + one small re-read of the catalog
    gz; writes every JSON artifact under `out_dir` and returns the
    result dict.  Exit-code logic: `exit_code` is 0 only when the sample
    mapping is resolved *and* the catalog carries no duplicate
    (chrom,pos) pair; otherwise 3 (the gate stays UNRESOLVED)."""
    vcf = Path(vcf)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    encoding = stream_kind(vcf)
    header = collect_vcf_header(vcf)
    n_samples = header["n_samples"]
    gt_index = header["gt_column_index"]
    sample_start = (9 if header.get("format_column_in_chrom") else 8)

    # Picked sites: use as-is (deterministic selection already done by
    # the mapping script), indexed for O(1) lookup during the pass.
    targets = set()
    target_records = {}
    if picked_sites:
        for chrom, pos in picked_sites:
            targets.add((chrom, pos))

    # 2.4 window: the first `records_per_chrom` biallelic-SNP records on
    # up to `window_chroms` nuclear chroms (default 400 x 5 = 2,000).
    capture = {chrom: [] for chrom in NUCLEAR[:max(0, window_chroms)]}
    summary = new_summary()

    catalog_handle = None
    if save_catalog:
        catalog_handle = gzip.open(
            out_dir / CATALOG_GZ, "wt", encoding="utf-8")
        catalog_handle.write("\t".join(CATALOG_COLUMNS) + "\n")

    record_index = 0
    sample_row_width_mismatch = None
    with open_vcf_stream(vcf) as stream:
        for line in stream:
            record = parse_vcf_line(line, sample_start)
            if record is None:
                continue
            record_index += 1  # data-line index, 1-based, header not counted
            if sample_row_width_mismatch is None:
                wf = len(record.get("sample_fields") or [])
                if wf != n_samples:
                    sample_row_width_mismatch = wf
            cat = update_summary(summary, record)
            if record.get("malformed") or record.get("pos") is None:
                continue
            if cat != CAT_SNP:
                continue
            chrom, pos = record["chrom"], record["pos"]
            if catalog_handle is not None:
                catalog_handle.write(catalog_line(record, record_index) + "\n")
            if chrom in capture and len(capture[chrom]) < records_per_chrom:
                capture[chrom].append(record)
            if (chrom, pos) in targets:
                target_records.setdefault((chrom, pos), record)
    if catalog_handle is not None:
        catalog_handle.close()

    # Catalog duplicate keys: second (bounded) pass over the small
    # catalog gz.  The two-pass view is what keeps memory flat; the
    # VCF stream itself is never unpacked into a plaintext file on
    # disk.
    catalog_path = out_dir / CATALOG_GZ
    dup_stats = catalog_duplicate_stats(catalog_path) if save_catalog else None
    if dup_stats is not None:
        summary["duplicate_chrom_position_pairs"] = \
            dup_stats["duplicate_chrom_position_pairs"]
        summary["duplicate_chrom_position_ref_alt_quadruples"] = \
            dup_stats["duplicate_chrom_position_ref_alt_quadruples"]

    # 2.3 mapping (the VCF semantic gate)
    table_rows = load_accession_table(accession_table)
    mapping = map_samples(header["sample_ids"], table_rows)

    # 2.4 window completion: pad/truncate each record's cells to the
    # header sample count and turn "GT:DP:AD"-style cells into GT-only
    # text.
    def gt_rows_for(records):
        rows = []
        for record in records:
            cells = [cell_text(f, gt_index)
                     for f in (record.get("sample_fields") or [])]
            if len(cells) < n_samples:
                cells += [""] * (n_samples - len(cells))
            rows.append(cells[:n_samples])
        return rows

    gt_audit = gt_window_counts({chrom: gt_rows_for(recs)
                                 for chrom, recs in capture.items()})
    gt_payload = {"window": {
                      "chroms": sorted(capture, key=lambda c: int(c)),
                      "records_per_chrom_cap": records_per_chrom,
                      "n_records_audited": gt_audit["combined"]["n_records"],
                      "format_column_assumed": sample_start == 9,
                      "sample_fields_width_mismatch": sample_row_width_mismatch,
                      "note": ("bounded window; no full genotype matrix; "
                               "GT states are reported as-is -- no "
                               "GT=1==HDF5-ALT and no dosage assumption")},
                  "per_chrom": gt_audit["per_chrom"],
                  "combined": gt_audit["combined"]}

    # 3.4 step 1 (VCF side)
    spot_entries = []
    if picked_sites:
        for chrom, pos in picked_sites:
            record = target_records.get((chrom, pos))
            if record is None:
                spot_entries.append(spot_site_not_found(chrom, pos))
            else:
                # `cells` is one GT text per sample column (VCF order);
                # the per-sample VCF-side state-code row records the
                # alignment against the HDF5 accession row so the HDF5
                # side can pair sample i with sample i.
                cells = gt_rows_for([record])[0]
                spot_entries.append(spot_site_report(
                    chrom, pos, cells,
                    per_sample=[spot_state_code(c) for c in cells]))
    spot_payload = None
    if picked_sites:
        spot_payload = {
            "n_sites": len(spot_entries),
            "n_samples": n_samples,
            "site_selection": ("the supplied picked-site list is used "
                               "as-is, in order; no resampling or "
                               "re-sorting in this script"),
            "state_code_legend": {
                "homozygous_0": "all present GT alleles indexed 0 (0/0-class)",
                "heterozygous": "mix of 0 and 1",
                "homozygous_1": "all present GT alleles indexed 1 (1/1-class)",
                "missing": "no GT / all-missing / partially-missing",
                "out_of_range_or_unexpected":
                    "allele index > 1, or an unparseable token"},
            "low_concordance_policy": ("low-concordance sites are listed "
                                       "for inspection: release/coordinate, "
                                       "multiallelic, imputation, "
                                       "heterozygosity, sample order, "
                                       "variant representation -- no "
                                       "threshold-based pass (spec 3.4)"),
            "per_site": spot_entries}

    # Outputs -------------------------------------------------------------
    source_payload = compose_source_json(
        vcf_path=vcf, official_md5=official_md5, source_url=source_url,
        download_time=download_time, download_tool=download_tool,
        bcftools_probe=(probe_bcftools(bcftools, vcf) if bcftools else None),
        stream_encoding=encoding)
    header_payload = {**header, "vcf_path": str(vcf),
                      "sample_id_map": mapping,
                      "sample_gate_note": (
                          "an unresolved sample_id_map leaves the VCF "
                          "semantic gate open (VCF_SEMANTICS_UNRESOLVED), "
                          "per spec 2.3"
                      )}
    _write_json(out_dir / SOURCE_JSON, source_payload)
    _write_json(out_dir / HEADER_JSON, header_payload)
    _write_json(out_dir / GT_JSON, gt_payload)
    write_summary_outputs(out_dir, summary, dup_stats,
                          catalog_path if save_catalog else None)
    if spot_payload is not None:
        _write_json(out_dir / SPOT_CHECK_JSON, spot_payload)

    # Gate ------------------------------------------------------------------
    catalog_problem = (save_catalog and
                       dup_stats["duplicate_chrom_position_pairs"] > 0)
    gate_ok = (mapping["status"] == "resolved" and not catalog_problem)
    result = {
        "status": "pass" if gate_ok else "fail",
        "gate": ("VCF_SEMANTICS_READY_FOR_CONVERSION" if gate_ok
                  else "VCF_SEMANTICS_UNRESOLVED"),
        "exit_code": 0 if gate_ok else 3,
        "vcf_path": str(vcf),
        "stream_encoding": encoding,
        "mapping_status": mapping["status"],
        "mapping": mapping,
        "header": header,
        "summary": summary,
        "gt_audit": gt_payload,
        "spot_check": spot_payload,
        "source": source_payload,
        "outputs": {name: str(out_dir / name) for name in
                    (SOURCE_JSON, HEADER_JSON, GT_JSON, SUMMARY_JSON, CATALOG_GZ,
                     SPOT_CHECK_JSON)
                    if name != SPOT_CHECK_JSON or spot_payload is not None
                    if not (name == CATALOG_GZ and not save_catalog)},
    }
    return result


def cell_text(sample_field, gt_index) -> str:
    """Return the GT-column text of one sample cell.  `gt_index` is the
    0-based position of the GT cell in the per-sample FORMAT cells -- the
    assumed convention is index 0 (GT-first sample cells, the release-file
    layout); None when the header declares no GT.  A missing index or a
    truncated cell returns "" and is counted as the `no_gt` cell
    downstream."""
    if gt_index is None:
        return ""
    parts = sample_field.split(":") if ":" in sample_field else [sample_field]
    return parts[gt_index] if len(parts) > gt_index else ""


def write_summary_outputs(out_dir, summary, dup_stats=None, catalog_path=None) -> Path:
    """Write the 3.2 `SUMMARY_JSON`.  When `dup_stats` (the return of
    `catalog_duplicate_stats`) is supplied the duplicate-key counts and
    the deterministic sample block (first/last SAMPLE_BLOCK catalog rows +
    catalog file SHA256) are embedded so *only* the small summary JSON
    ever lands under -- with the catalog gz itself kept out of git.  The
    written path is returned for chaining."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    dup_pairs = dup_stats["duplicate_chrom_position_pairs"] if dup_stats else None
    dup_quads = dup_stats["duplicate_chrom_position_ref_alt_quadruples"] \
        if dup_stats else None
    sample_block = None
    if catalog_path is not None:
        block_rows = dup_stats["first_rows"] if dup_stats else []
        last_rows = dup_stats["last_rows"] if dup_stats else []
        sample_block = {
            "first_rows": block_rows[:SAMPLE_BLOCK],
            "last_rows": last_rows[-SAMPLE_BLOCK:],
            "catalog_file": Path(catalog_path).name,
            "catalog_file_sha256": file_sha256(catalog_path)
                if Path(catalog_path).exists() else None,
            "note": ("deterministic sample block + catalog gz SHA256: the "
                     "small summary stays under git, the catalog gz does not")}
    payload = {
        "total_records": summary["total_records"],
        "nuclear_records": summary["nuclear_records"],
        "per_chrom": dict(summary["per_chrom"]),
        "type_counts": dict(summary["type_counts"]),
        "catalog_rows": summary["catalog_rows"],
        "excluded_examples": {k: v for k, v in
                              summary["excluded_examples"].items()},
        "duplicate_chrom_position_pairs": dup_pairs,
        "duplicate_chrom_position_ref_alt_quadruples": dup_quads,
        "catalog_sample_block": sample_block,
    }
    out_path = out_dir / SUMMARY_JSON
    _write_json(out_path, payload)
    return out_path


def _write_json(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
                    encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="1001G v3.1 population VCF audit (TASK02/TASK03, VCF "
                    "side only).")
    parser.add_argument("--vcf", required=True,
                        help="path to the population VCF (.vcf or .vcf.gz)")
    parser.add_argument("--out", default=DEFAULT_OUT_DIR,
                        help="output directory (default: %(default)s)")
    parser.add_argument("--accession-table", default=DEFAULT_ACC_TABLE,
                        help="official 1135 accession table JSON "
                             "(default: %(default)s)")
    parser.add_argument("--official-md5", default=None,
                        help="official MD5 of the release file, if known")
    parser.add_argument("--source-url", default=None,
                        help="release URL the file was downloaded from")
    parser.add_argument("--download-time", default=None,
                        help="ISO timestamp of the download")
    parser.add_argument("--download-tool", default=None,
                        help="aria2c / curl -C - / wget -c ...")
    parser.add_argument("--bcftools", default=None,
                        help="path to a bcftools binary for the 2.2 "
                             "header-read probe (optional)")
    parser.add_argument("--gt-window-chrom", type=int, default=5,
                        help="how many leading nuclear chroms to include "
                             "in the 2.4 window (default: %(default)s)")
    parser.add_argument("--gt-window-records-per-chrom", type=int,
                        default=400,
                        help="per-chrom biallelic-SNP record cap for the "
                             "2.4 window (default: %(default)s) -> up to "
                             "5*400=2,000 records on chrom 1-5")
    parser.add_argument("--picked-sites", default=None,
                        help="JSON of [chrom, pos] pairs from the mapping "
                             "script (3.4 step 1); used as-is")
    parser.add_argument("--no-catalog", action="store_true",
                        help="skip the catalog gz (use for small tests)")
    args = parser.parse_args(argv)

    picked = load_picked_sites(args.picked_sites) if args.picked_sites else None
    result = run_audit(vcf=args.vcf, out_dir=args.out,
                       accession_table=args.accession_table,
                       official_md5=args.official_md5,
                       source_url=args.source_url,
                       download_time=args.download_time,
                       download_tool=args.download_tool,
                       bcftools=args.bcftools,
                       window_chroms=args.gt_window_chrom,
                       records_per_chrom=args.gt_window_records_per_chrom,
                       picked_sites=picked,
                       save_catalog=not args.no_catalog)

    # Machine summary: one-line JSON on stdout.
    summary_line = {"module": "scripts.audit_1001g_vcf",
                    "run": result["gate"],
                    "exit_code": result["exit_code"],
                    "status": result["status"],
                    "mapping": result["mapping_status"],
                    "catalog_rows": result["summary"]["catalog_rows"],
                    "dup_chr_pos": result["summary"].get(
                        "duplicate_chrom_position_pairs"),
                    "n_picked": result["spot_check"]["n_sites"]
                        if result["spot_check"] else 0,
                    "out_dir": args.out}
    print(json.dumps(summary_line, sort_keys=True))
    return result["exit_code"]


if __name__ == "__main__":
    sys.exit(main())
