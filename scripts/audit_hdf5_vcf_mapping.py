"""TASK03 HDF5-side half of the 1001G v3.1 population-VCF semantic audit.

Implements spec section 3.3 (HDF5 position mapping) and 3.4 step 2
(genotype concordance spot-check, HDF5 side) of TASK03 of the delivery
doc tmp/at_pheno_Qwen3.8-27B_数据语义阶段执行交付包.md.  The VCF-side half
(2.x, 3.1, 3.2, 3.4 step 1: checksums, sample mapping, GT audit, the
biallelic-SNP catalog and the VCF-side per-site spot JSON) is finished
in scripts/audit_1001g_vcf.py; this module is its HDF5-side join
partner and the only place of the two that touches the imputation HDF5.

Style: one plain module in the scripts/ house style
(scripts/trait_resolution.py), no classes, no h5py and no numpy
requirement in the pure core (fixtures fake the positions axis and the
per-site state lists with plain sequences or numpy arrays).  h5py, the
`data` extra, is needed only by the CLI-adjacent helper
`load_hdf5_positions` for real HDF5 files.  The 18 GB population VCF
under data/raw/1001g_v3.1/ is deliberately never touched here; tests
run on tiny self-made fixtures under tmp_path only.

Pure helpers (all fixture-tested, none of them opens HDF5/VCF files):

- `read_catalog`            streaming (chrom, pos, ref, alt, idx) rows
                            from the 3.2 biallelic-SNP catalog tsv.gz
                            (re-implements the audit_1001g_vcf contract
                            so this module is importable standalone);
- `validate_region_spans`   contiguity / tiling validation of the
                            [start, stop) region spans;
- `mapping_from_catalog`    spec 3.3 streaming coordinate join +
                            per-coordinate statistics + the fixed
                            hash/salt site pool of spec 3.4;
- `selected_sites`          deterministic [chrom, pos] pool ordering;
- `index_vcf_sites`         (chrom, pos) -> VCF-side per-site entry;
- `index_site_values`       (chrom, pos) -> HDF5-side 0/1 state list;
- `reconcile_spot`          spec 3.4 step 2: per-site + overall
                            reading under the two orientation
                            hypotheses, and the `suspect_sites`
                            enumeration (no threshold-forced pass);
- `fuse_site_mapping`       join-level fusion of mapping summary +
                            VCF-side per-site JSON + sample overlap;
- `spot_gate`               PASS / UNRESOLVED verdict mapping;
- `synthetic_hdf5_fixture`  the 5-region mini positions fixture.

CLI: `main(argv=None)` writes exactly two contract artifacts into
`--out`:

    hdf5_vcf_position_mapping_summary.json
    hdf5_vcf_concordance_spot_check.json

and prints a single-line machine summary JSON
(`{"module","status","exit_code", ...}`).  Exit codes: 0 = PASS,
3 = gate UNRESOLVED, 4 = input unreadable / malformed (one-line JSON
error, no traceback).
"""

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

# The seven spec-mandated per-site / overall keys, in spec order.
SPEC_KEYS = (
    "concordance_if_hdf5_0_is_REF",
    "concordance_if_hdf5_1_is_REF",
    "best_orientation",
    "best_concordance",
    "n_comparable",
    "n_heterozygous",
    "n_missing_vcf",
)
KEY_VP0 = SPEC_KEYS[0]
KEY_VP1 = SPEC_KEYS[1]
KEY_BO = SPEC_KEYS[2]
KEY_BC = SPEC_KEYS[3]
KEY_NC = SPEC_KEYS[4]
KEY_NH = SPEC_KEYS[5]
KEY_NM = SPEC_KEYS[6]

NUCLEAR = ("1", "2", "3", "4", "5")
NUCLEOTIDES = frozenset("ACGT")

# The 3.4 state-code names mirror audit_1001g_vcf.SPOT_STATE_KEYS; the
# two orientation hypotheses follow the spec ordering.
SPOT_STATE_KEYS = ("homozygous_0", "heterozygous", "homozygous_1",
                   "missing", "out_of_range_or_unexpected")
ORIENT_LABELS = ("hdf5_0_is_ref", "hdf5_1_is_ref")
# label -> (the concordance key measuring THAT label, the other one)
_PKEYS = {"hdf5_0_is_ref": (KEY_VP0, KEY_VP1),
          "hdf5_1_is_ref": (KEY_VP1, KEY_VP0)}
HOMO_STATE_CODES = ("homozygous_0", "homozygous_1")

# Spec 3.4: fixed hash/salt, 1,000-2,000 sites -> cap at 2,000.
SALT = "at_pheno_task03_v1"
SITE_POOL_CAP = 2000
EXAMPLE_CAP = 20

# Sub-check items named by spec 3.4, carried alongside any suspected
# site enumeration below.
SUSPECT_CHECKLIST = (
    "release/coordinate (build release / coordinate system, incl. "
    "any page/coordination offset)",
    "multiallelic & imputation (multiallelic / imputation decomposition "
    "residue in the variant representation)",
    "heterozygosity (an HDF5 0/1 pair cannot carry a true heterozygote)",
    "SAMPLE-ORDERING (shared-sample alignment: VCF sample-header order "
    "vs HDF5 accession order)",
    "version/representation (REF/ALT direction, flavor, ALT-list "
    "order between release builds)",
)

NO_DOSAGE_NOTE = ("HDF5 0/1 values are never converted into ALT-"
                  "dosages; the two orientation hypotheses are label-"
                  "only, and no orientation is ever written back to the "
                  "HDF5 file (spec 3.4).")


# ---------------------------------------------------------------------------
# streaming helpers
# ---------------------------------------------------------------------------

def file_sha256(path) -> str:
    """SHA-256 of `path` in 1 MiB chunks (provenance record; the 18 GB
    VCF is never re-hashed here, only small files like the catalog gz)."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


_CATALOG_COLUMNS = ("chromosome", "position", "ref", "alt",
                    "source_record_index")


def read_catalog(path):
    """Yield (chrom, pos_int, ref, alt, idx_int) rows from the 3.2
    biallelic-SNP catalog tsv.gz written by audit_1001g_vcf.  The
    single header line is skipped and never counted; the gzip is read
    through stdlib gzip and is never unpacked into a plaintext copy
    on disk (3.2 contract)."""
    with gzip.open(Path(path), "rt", encoding="utf-8") as handle:
        header = handle.readline()
        if header.rstrip("\n") != "\t".join(_CATALOG_COLUMNS):
            raise ValueError(f"catalog header mismatch: {header!r}")
        for row in handle:
            row = row.rstrip("\n")
            if not row:
                continue
            chrom, pos, ref, alt, idx = row.split("\t")
            yield chrom, int(pos), ref, alt, int(idx)


def catalog_content_md5(path):
    """MD5 of the whole UNPACKED TSV byte stream (header line included).

    This is the logical-content hash: identical across re-gzips of the
    same TSV text (gzip header mtime is not part of it).  Contrast with
    the .tsv.gz file SHA256 (`compressed_artifact_sha256`): that one
    covers the exact on-disk .gz bytes and changes whenever the file is
    re-gzipped even with identical content, so it verifies transfer of
    a specific artifact, NOT content identity.  Use
    logical_content_md5 to compare audit runs; use compressed
    artifact sha for file transfer only."""
    with gzip.open(Path(path), "rb") as handle:
        h = hashlib.md5()
        for chunk in iter(lambda: handle.read(1024*1024), b""):
            h.update(chunk)
        return h.hexdigest()


# ---------------------------------------------------------------------------
# 3.3 position mapping
# ---------------------------------------------------------------------------

def _row_policy_ok(row) -> bool:
    """3.1 policy predicate: nuclear chrom, biallelic, ref/alt a single
    A/C/G/T.  Whatever else a catalog row carries (multi-allelic ALT
    lists, N, non-nuclear, indel-shaped) never enters the matching
    sets -- such rows are tagged and counted, not split (spec 3.3/3.4:
    no fuzzy matching, no decomposition)."""
    chrom, _pos, ref, alt, _idx = row
    return (chrom in NUCLEAR
            and len(ref) == 1 and ref in NUCLEOTIDES
            and len(alt) == 1 and alt in NUCLEOTIDES)


def validate_region_spans(region_spans, chroms, n_positions):
    """The [start, stop) region spans join the positions axis to the
    per-region chrom labels; they must be consecutive (each start
    equals the previous stop), ordered, and tile the prefix of the
    axis.  Returns the spans as (start, stop) tuples; raises ValueError
    on a manifest/span provenance problem (never an ambiguous state)."""
    if len(chroms) != len(region_spans):
        raise ValueError(
            f"chrom label count ({len(chroms)}) differs from region-span "
            f"count ({len(region_spans)})")
    out = []
    prev_stop = 0
    for index, span in enumerate(region_spans):
        start, stop = int(span[0]), int(span[1])
        if start != prev_stop:
            raise ValueError(
                f"region span {index} starts at {start} but the previous "
                f"span stops at {prev_stop} (spans must be consecutive)")
        if stop <= start:
            raise ValueError(
                f"region span {index} [{start}, {stop}) is empty or "
                f"inverted")
        out.append((start, stop))
        prev_stop = stop
    if prev_stop > n_positions:
        raise ValueError(
            f"region spans end at {prev_stop} but the positions axis has "
            f"{n_positions} rows")
    return out


def selected_sites(coords, cap=SITE_POOL_CAP, salt=SALT):
    """The spec 3.4 site pool: deterministic `md5(chr:pos:salt)` order
    over the coordinate rows where a hold count is 1, capped at `cap`
    (the spec range of 1,000-2000 sites).  Returns [[chrom, pos], ...]
    -- directly usable as the VCF-side --picked-sites input (a list of
    [chrom, pos] pairs with string/int-tolerant pos)."""
    ordered = sorted(
        coords,
        key=lambda key: hashlib.md5(
            f"{key[0]}:{key[1]}:{salt}".encode("utf-8")).hexdigest())
    return [[str(chrom), str(pos)] for chrom, pos in ordered[:cap]]


def mapping_from_catalog(cat_rows, h5_positions, region_spans,
                         chroms=None) -> dict:
    """Spec 3.3: a streaming-merge / sorted-join coordinate join of the
    yielded catalog rows against the positions axis.

    cat_rows -- iterable of (chrom, pos_int, ref, alt, idx) tuples as
        yielded by `read_catalog` (the biallelic-SNP policy set; rows
        that fail the 3.1 predicate (multi-allelic etc.) are counted
        in `policy_skipped_rows` and never participate in a match);
    h5_positions -- the positions axis as a numpy-style sequence: the
        5 consecutive nuclear regions in the order of `region_spans`
        (fixture arrays in tests; `load_hdf5_positions` in
        production);
    region_spans -- [start, stop) index pairs into `h5_positions` (see
        validate_region_spans);
    chroms -- one label per region (default "1".."N", N =
        n_regions; pass explicit labels when the spans are non-nuclear
        or beyond 5 regions).

    Returns a dict of the 3.3 quantities.  All per-coord statistics
    live in `(chrom, pos)` counters, so `coord_hits` can be stated
    without re-storing any row (bounded memory) and without ever
    inferring a 0/1 REF/ALT orientation from the match (spec 3.3; see
    orientation_note).
    """
    if chroms is None:
        chroms = [NUCLEAR[i] if i < len(NUCLEAR) else str(i + 1)
                  for i in range(len(region_spans))]
    chroms = [str(c) for c in chroms]
    spans = validate_region_spans(region_spans, chroms, len(h5_positions))

    # HDF5 side: (chrom, pos) -> row count
    h5_counts = {}
    for chrom, (start, stop) in zip(chroms, spans):
        for value in h5_positions[start:stop]:
            key = (chrom, int(value))
            h5_counts[key] = h5_counts.get(key, 0) + 1
    h5_keys = set(h5_counts)

    # VCF side: candidate counts per (chrom, pos) and per full
    # quadruple (quadruple = (chrom, pos, ref, alt))
    cat_pairs = {}
    cat_quads = {}
    catalog_rows = 0
    skipped_examples = []
    for row in cat_rows:
        chrom, pos, ref, alt, _idx = row
        if not _row_policy_ok(row):
            skipped_examples.append([str(chrom), str(pos), ref, alt])
            continue
        key = (str(chrom), int(pos))
        cat_pairs[key] = cat_pairs.get(key, 0) + 1
        quad = (str(chrom), int(pos), str(ref), str(alt))
        cat_quads[quad] = cat_quads.get(quad, 0) + 1
        catalog_rows += 1

    # Three-way split inside the coordinate space
    unique_matches = [k for k in sorted(cat_pairs)
                      if k in h5_keys and cat_pairs[k] == 1]
    multi_candidates = [k for k in sorted(cat_pairs)
                        if k in h5_keys and cat_pairs[k] > 1]
    hdf5_only = [k for k in sorted(h5_keys) if k not in cat_pairs]
    vcf_only = [k for k in sorted(cat_pairs) if k not in h5_keys]

    dup_h5 = {k: c for k, c in h5_counts.items() if c > 1}
    dup_h5_rows = sum(c - 1 for c in dup_h5.values())
    dup_h5_per_chrom = {chrom: sum(1 for (cc, _p) in dup_h5 if cc == chrom)
                        for chrom in chroms}

    pool = selected_sites(unique_matches)

    return {
        "n_regions": len(chroms),
        "chroms": list(chroms),
        "region_spans": [[s, e] for s, e in spans],
        "hdf5_coordinate_rows": sum(h5_counts.values()),
        "unique_coordinates": len(h5_keys),
        "duplicated_positions": len(dup_h5),
        "duplicated_position_rows": dup_h5_rows,
        "duplicated_positions_per_chrom": dup_h5_per_chrom,
        "duplicated_coordinate_examples":
            [[k[0], str(k[1])] for k in sorted(dup_h5)[:EXAMPLE_CAP]],
        "catalog_rows": catalog_rows,
        "unique_matches": len(unique_matches),
        "unique_matches_sample": [list(k) for k in pool[:EXAMPLE_CAP]],
        "multiple_candidates": len(multi_candidates),
        "multiple_candidates_sample": [
            {"chrom": k[0], "pos": str(k[1]),
             "n_catalog_candidates": int(cat_pairs[k])}
            for k in multi_candidates[:EXAMPLE_CAP]],
        "hdf5_only_coordinates": len(hdf5_only),
        "hdf5_only_sample": [[k[0], str(k[1])]
                             for k in hdf5_only[:EXAMPLE_CAP]],
        "vcf_only_coordinates": len(vcf_only),
        "vcf_only_sample": [[k[0], str(k[1])]
                            for k in vcf_only[:EXAMPLE_CAP]],
        "duplicate_chrom_position_pairs":
            sum(1 for c in cat_pairs.values() if c > 1),
        "duplicate_chrom_position_ref_alt_quadruples":
            sum(1 for c in cat_quads.values() if c > 1),
        "policy_skipped_rows": len(skipped_examples),
        "policy_skipped_examples": skipped_examples[:EXAMPLE_CAP],
        "coord_hits": len(unique_matches) + len(multi_candidates),
        "coord_hits_note": ("HDF5 coordinate rows holding at least one "
                            "biallelic-catalog candidate; stated from "
                            "the per-coordinate dict, not from re-stored "
                            "rows (bounded memory)"),
        "site_pool": {
            "n_pool_capped": len(pool),
            "cap": SITE_POOL_CAP,
            "salt": SALT,
            "selection": ("generation-order-hash md5(chrom:pos:salt) "
                          "over the coordinates with a 1-instance "
                          "catalog count"),
        },
        "site_pool_sample": pool,
        "no_orientation_claimed": True,
        "orientation_note": ("spec 3.3: coordinate matching alone never "
                              "implies a 0/1 REF/ALT direction; "
                              + NO_DOSAGE_NOTE),
    }


# ---------------------------------------------------------------------------
# 3.4 step 2: per-site concordance join
# ---------------------------------------------------------------------------

_STATE_TOKENS = {
    "0": "homozygous_0",
    "homozygous_0": "homozygous_0",
    "1": "homozygous_1",
    "homozygous_1": "homozygous_1",
    "H": "heterozygous",
    ".": "missing",
    "X": "out_of_range_or_unexpected",
    "heterozygous": "heterozygous",
    "missing": "missing",
    "out_of_range_or_unexpected": "out_of_range_or_unexpected",
}


def _vcf_state_code(token):
    """One shared-sample VCF-side state token -> a SPOT_STATE_KEYS code.
    Accepts 0/1 as int or str, H/heterozygous, ./missing/blank, the
    audit_1001g_vcf state-code names, and 'X'; anything else is
    out_of_range_or_unexpected (state classification never raises)."""
    if token in (0, "0"):
        return "homozygous_0"
    if token in (1, "1"):
        return "homozygous_1"
    text = str(token).strip().lower()
    if text in ("h", "hetero", "heterozygous"):
        return "heterozygous"
    if text in ("", ".", "missing", "none"):
        return "missing"
    return _STATE_TOKENS.get(text, "out_of_range_or_unexpected")


def _h5_state_bucket(value):
    """One HDF5-side shared-sample state -> 'zero' / 'one' /
    'missing' / 'invalid'.  Out-of-range values (anything that is not
    a 0/1 integer, 0/1 string, or a missing marker) are reported,
    never coerced into a state or a dosage."""
    if value is None:
        return "missing"
    if isinstance(value, bool):
        return "one" if value else "zero"
    if isinstance(value, (int, float)):
        iv = int(value)
        if iv == 0:
            return "zero"
        return "one" if iv == 1 else "invalid"
    text = str(value).strip()
    if text in ("", ".", "missing", "None"):
        return "missing"
    if text == "0":
        return "zero"
    if text == "1":
        return "one"
    return "invalid"


def index_vcf_sites(vcf_side):
    """(chrom, pos_int) -> VCF-side per-site entry.  Accepts the whole
    VCF-side JSON dict as emitted by audit_1001g_vcf.py (its
    `per_site` list is honored), a bare per-site list, or a nested
    {chrom: {pos: entry}} / {chrom: [entry...]} mapping.  Every entry
    carries the SPEC_KEYS keys; an entry without a matching record on
    the VCF side (spot_site_not_found shape) still carries them, as
    null/zero counts."""
    if isinstance(vcf_side, dict):
        if isinstance(vcf_side.get("per_site"), list):
            entries = vcf_side["per_site"]
        else:
            entries = []
            for value in vcf_side.values():
                if isinstance(value, dict):
                    entries.extend(value.values())
                elif isinstance(value, (list, tuple)):
                    entries.extend(value)
    else:
        entries = list(vcf_side or [])
    index = {}
    for entry in entries:
        index[(str(entry.get("chrom", "")), int(entry.get("pos", -1)))] = entry
    return index


def index_site_values(values):
    """(chrom, pos_int) -> HDF5-side per-shared-sample state list.
    `values` is the datasource-join input: in production a join result
    over the shared-sample accessions; in tests a tiny JSON file.
    Accepted shapes: nested {chrom: {pos: [states...]}}, flat
    {(chrom, pos): [states...]} (tuple keys), or a list of
    {"chrom", "pos", "values"} dicts.  An absent key says the site is
    missing on the HDF5 side; None entries inside a list say a shared
    sample is missing (counted, never imputed)."""
    index = {}

    def absorb(chrom, pos, states):
        index[(str(chrom), int(pos))] = list(states or [])

    if isinstance(values, dict):
        if all(isinstance(v, list) for v in values.values()):
            for key, states in values.items():
                if isinstance(key, (tuple, list)) and len(key) == 2:
                    absorb(key[0], key[1], states)
                else:
                    absorb(key, "0", states)
        else:
            for chrom, inner in values.items():
                for pos, states in inner.items():
                    absorb(chrom, pos, states)
    else:
        for entry in values or []:
            absorb(entry["chrom"], entry["pos"], entry.get("values", []))
    return index


def reconcile_spot(hdf5_vals, vcf_per_site, orientation_pass=None,
                   sample_overlap=None) -> dict:
    """Spec 3.4 step 2: the per-site AND the overall reading under the
    two orientation hypotheses.

    Arguments:
      `hdf5_vals` -- the HDF5-side per-shared-sample 0/1 states for the
        same sites, indexed via `index_site_values` (a dict keys, a
        flat {(chrom,pos): [...]}, or a list of {"chrom","pos",
        "values"}; production: a datasource join, tests: a JSON
        fixture).
      `vcf_per_site` -- the VCF-side per-site JSON (SPEC_KEYS entries,
        as emitted by audit_1001g_vcf.py), indexed via
        `index_vcf_sites`; an entry may carry a `per_sample` list of
        per-shared-sample state tokens (0/1, H, ., X, or the
        state-code names) which is the alignment reference for the
        pair-wise comparison.
      `orientation_pass` -- optional hypothesized-orientation label
        ("hdf5_0_is_ref" / "hdf5_1_is_ref"); when given, a site whose
        passed orientation underperforms its alternative is
        enumerated as a suspect.  The label is diagnostic only: it is
        not written back to the HDF5 file and never threshold-passed.
      `sample_overlap` -- shared-sample context dict
        (shared_accessions / shared_n / n_samples_vcf /
        n_samples_hdf5 / note), passed through as-is; the
        SAMPLE-ORDERING fact stays a manual cross-check.

    Per site the result reports the SPEC_KEYS (exact spec spellings;
    `n_comparable` / `n_heterozygous` / `n_missing_vcf` are re-derived
    from the shared-sample `per_sample` alignment reference when it is
    present, else mirror the VCF-side entry counts; 0 when the site
    has no VCF-side entry), plus:
      concordance readings under both orientation hypotheses,
      best_orientation ("ambiguous" on a tie, None when fewer than
      one comparable pair), best_concordance, the pair-set size, the
      HDF5-side missing / invalid counts, and a `suspect` flag with
      reasons.

    Suspect enumeration (spec 3.4: a LIST, not a boolean verdict, no
    threshold-forced pass) marks a site when:
      (a) no comparable pair although the VCF side reports comparable
          samples (imputation / missing state / sample-ordering
          surface);
      (b) the two orientation concordances tie on a non-empty pair
          set (ambiguous orientation signal);
      (c) all shared-sample HDF5 states at the site are missing or
          invalid (datasource-join / imputation coverage);
      (d) the site is absent from the HDF5-side values entirely
          (release/coordinate offset);
      (e) an orientation_pass was given and it underperforms its
          alternative on the measured pair set.
    """
    if orientation_pass not in (None,) + tuple(_PKEYS):
        raise ValueError(
            f"orientation_pass must be one of {list(_PKEYS)} or None, "
            f"got {orientation_pass!r}")
    h5_index = index_site_values(hdf5_vals)
    vcf_index = index_vcf_sites(vcf_per_site) if vcf_per_site else {}
    site_keys = list(dict.fromkeys(list(h5_index) + list(vcf_index)))

    per_site = []
    tot_match0 = 0
    tot_pairs = 0
    tot_nc = 0
    tot_nh = 0
    tot_nm = 0
    for chrom, pos in site_keys:
        in_h5 = (chrom, pos) in h5_index
        in_vcf = (chrom, pos) in vcf_index
        h5_states = h5_index.get((chrom, pos), [])
        vcf_entry = vcf_index.get((chrom, pos)) or {}
        per_sample = vcf_entry.get("per_sample") if in_vcf else None

        # VCF-side spec keys (0 when absent).  The entry counts cover
        # the whole VCF population; the 3.4 spot check reports counts
        # over the SHARED aligned samples, so when the `per_sample`
        # alignment reference is present the three counts are
        # re-derived from it (homo codes -> comparable, H ->
        # heterozygous, '.' -> missing); without it the entry mirror
        # is used.
        states = vcf_entry.get("states") or {}
        if per_sample is not None:
            codes = [_vcf_state_code(t) for t in per_sample]
            n_comparable = sum(1 for c in codes if c in HOMO_STATE_CODES)
            n_heterozygous = sum(1 for c in codes if c == "heterozygous")
            n_missing_vcf = sum(1 for c in codes if c == "missing")
        else:
            n_comparable = int(
                vcf_entry.get("n_comparable")
                or (states.get("homozygous_0", 0)
                    + states.get("homozygous_1", 0)))
            n_heterozygous = int(vcf_entry.get("n_heterozygous")
                                 or states.get("heterozygous", 0))
            n_missing_vcf = int(vcf_entry.get("n_missing_vcf")
                                or states.get("missing", 0))

        buckets = {"zero": 0, "one": 0, "missing": 0, "invalid": 0}
        for value in h5_states:
            buckets[_h5_state_bucket(value)] += 1

        # Pair-wise concordance over the shared samples (0/1 space only;
        # the two hypothesis readings are complements by construction,
        # so no dosage conversion is possible here)
        c0 = c1 = None
        best_orientation = None
        best_concordance = None
        n_aligned = 0
        n_pairs = 0
        match0 = 0
        if in_h5 and in_vcf and per_sample is not None:
            n_aligned = min(len(per_sample), len(h5_states))
            for i in range(n_aligned):
                code = _vcf_state_code(per_sample[i])
                if code not in HOMO_STATE_CODES:
                    continue
                bucket = _h5_state_bucket(h5_states[i])
                if bucket not in ("zero", "one"):
                    continue
                n_pairs += 1
                if (bucket == "zero") == (code == "homozygous_0"):
                    match0 += 1
            tot_match0 += match0
            tot_pairs += n_pairs
            if n_pairs:
                c0 = match0 / n_pairs
                c1 = (n_pairs - match0) / n_pairs
                if c0 > c1:
                    best_orientation = "hdf5_0_is_ref"
                    best_concordance = c0
                elif c1 > c0:
                    best_orientation = "hdf5_1_is_ref"
                    best_concordance = c1
                else:
                    best_orientation = "ambiguous"
                    best_concordance = c0

        suspect_reasons = []
        if (not in_h5 and in_vcf):
            suspect_reasons.append(
                "site is on the VCF side but has no HDF5-side value: "
                "release/coordinate offset or imputation-coverage "
                "candidate (spec 3.4 cross-check: release/coordinate, "
                "version/representation)")
        if in_h5 and len(h5_states) and \
                (buckets["missing"] + buckets["invalid"]) == len(h5_states):
            suspect_reasons.append(
                "all shared-sample HDF5 states at this site are "
                "missing/invalid: datasource-join or imputation "
                "coverage issue")
        if in_h5 and in_vcf and per_sample is not None \
                and n_pairs == 0 and n_comparable:
            suspect_reasons.append(
                "no comparable pair although the VCF side reports "
                "comparable samples: inspect missing-states and "
                "SAMPLE-ORDERING (spec 3.4: heterozygosity, "
                "imputation, sample order)")
        if best_orientation == "ambiguous":
            suspect_reasons.append(
                "orientation tie on a non-empty pair set: the 0/1 "
                "reference direction is not resolvable from this site "
                "alone; inspect the suspect-checklist factors")
        if orientation_pass in _PKEYS and n_pairs:
            pass_key, other_key = _PKEYS[orientation_pass]
            val_pass = c0 if pass_key == KEY_VP0 else c1
            val_other = c1 if pass_key == KEY_VP0 else c0
            if val_pass is not None and val_pass < val_other:
                suspect_reasons.append(
                    "orientation_pass=%r underperforms its alternative "
                    "(%.4f < %.4f): the passed hypothesis contradicts "
                    "the measured concordance; no threshold-forced "
                    "pass (spec 3.4)" % (orientation_pass, val_pass,
                                          val_other))

        per_site.append({
            KEY_VP0: c0,
            KEY_VP1: c1,
            KEY_BO: best_orientation,
            KEY_BC: best_concordance,
            KEY_NC: n_comparable,
            KEY_NH: n_heterozygous,
            KEY_NM: n_missing_vcf,
            "chrom": chrom,
            "pos": pos,
            "record_found": bool(vcf_entry.get("record_found", True))
                if in_vcf else False,
            "vcf_side_present": in_vcf,
            "hdf5_side_present": in_h5,
            "per_sample_available": per_sample is not None,
            "n_aligned": n_aligned,
            "n_pairs_comparable": n_pairs,
            "n_hdf5_states": len(h5_states) if in_h5 else 0,
            "n_hdf5_missing": buckets["missing"],
            "n_hdf5_invalid": buckets["invalid"],
            "suspect": bool(suspect_reasons),
            "suspect_reasons": suspect_reasons,
            "suspect_checklist": SUSPECT_CHECKLIST if suspect_reasons else [],
        })
        tot_nc += n_comparable
        tot_nh += n_heterozygous
        tot_nm += n_missing_vcf

    overall_c0 = (tot_match0 / tot_pairs) if tot_pairs else None
    overall_c1 = ((tot_pairs - tot_match0) / tot_pairs) if tot_pairs else None
    if overall_c0 is not None and overall_c0 > overall_c1:
        overall_orientation = "hdf5_0_is_ref"
        overall_best = overall_c0
    elif overall_c0 is not None and overall_c1 > overall_c0:
        overall_orientation = "hdf5_1_is_ref"
        overall_best = overall_c1
    elif overall_c0 is not None:
        overall_orientation = "ambiguous"
        overall_best = overall_c0
    else:
        overall_orientation = None
        overall_best = None

    overall = {
        KEY_VP0: overall_c0,
        KEY_VP1: overall_c1,
        KEY_BO: overall_orientation,
        KEY_BC: overall_best,
        KEY_NC: tot_nc,
        KEY_NH: tot_nh,
        KEY_NM: tot_nm,
        "n_pairs_comparable": tot_pairs,
        "n_sites_with_pairs": sum(1 for e in per_site
                                  if e["n_pairs_comparable"]),
    }
    return {
        "n_sites": len(per_site),
        "per_site": per_site,
        "overall": overall,
        "suspect_sites": [e for e in per_site if e["suspect"]],
        "suspect_checklist": SUSPECT_CHECKLIST,
        "orientation_pass": orientation_pass,
        "sample_overlap": sample_overlap or {},
        "no_dosage_note": NO_DOSAGE_NOTE,
        "note": ("concordance over the shared samples, per pair "
                 "(homozygous VCF state vs same-sample HDF5 0/1 "
                 "state); the two hypothesis readings are complements "
                 "on this 0/1 space -- no GT re-encoding, no dosing, "
                 "no orientation write-back (spec 3.4)"),
    }


# ---------------------------------------------------------------------------
# 3.3 + 3.4 fusion: mapping pool x VCF-side per-site x shared samples
# ---------------------------------------------------------------------------

def fuse_site_mapping(mapping_summary, vcf_side_per_site, sample_overlap):
    """Join-level fusion of the 3.3 mapping summary, the VCF-side
    3.4-step-1 per-site JSON, and the shared-sample overlap context
    (the documented manual input of `--sample-file`).  Returns a dict
    with:

    `site_pool` -- how many of the fixed hash/salt pool entries the
      VCF side resolved (`pool_sites_with_vcf` / `..._missing` lists),
    `sample_overlap` -- shared_n / n_samples_vcf / n_samples_hdf5 /
      shared_fraction (+ the note that SAMPLE-ORDERING remains a
      manual cross-check, never inferred),
    `orientation_status` -- "unresolved" until the spot join
      (`reconcile_spot`) resolves it; plus the no-dosage note.
    """
    pool_set = set()
    for item in (mapping_summary or {}).get("site_pool_sample", []):
        if isinstance(item, (list, tuple)) and len(item) == 2:
            # pool entries are [chrom, pos] pairs; match the (chrom-str,
            # pos-int) keying of index_vcf_sites
            pool_set.add((str(item[0]), int(item[1])))

    vcf_index = index_vcf_sites(vcf_side_per_site) \
        if vcf_side_per_site is not None else None
    found_in_vcf = []
    missing_in_vcf = []
    if vcf_index is not None:
        for key in sorted(pool_set):
            if key in vcf_index:
                found_in_vcf.append([key[0], str(key[1])])
            else:
                missing_in_vcf.append([key[0], str(key[1])])

    so = sample_overlap or {}
    shared_n = so.get("shared_n")
    if shared_n is None and so.get("shared_accessions"):
        shared_n = len(so["shared_accessions"])
    n_vcf = so.get("n_samples_vcf")
    n_h5 = so.get("n_samples_hdf5")
    if None not in (shared_n, n_vcf, n_h5) and min(n_vcf, n_h5) > 0:
        shared_fraction = float(shared_n) / float(min(n_vcf, n_h5))
    else:
        shared_fraction = None

    return {
        "site_pool": {
            "n": len(pool_set),
            "salt": (mapping_summary or {})
                    .get("site_pool", {}).get("salt"),
            "pool_sites_with_vcf": found_in_vcf if vcf_index else None,
            "pool_sites_missing_in_vcf": missing_in_vcf if vcf_index
                                       else None,
            "n_missing_in_vcf": (len(missing_in_vcf)
                                 if vcf_index else None),
            "note": ("pool entry `missing` = the VCF-side step-1 did "
                     "not resolve that site to a biallelic-SNP data "
                     "line (or the record was not in the picked "
                     "site-pool): a release/coordinate misalignment "
                     "candidate, not a silent pass"),
        },
        "sample_overlap": {
            "shared_n": shared_n,
            "n_samples_vcf": n_vcf,
            "n_samples_hdf5": n_h5,
            "shared_fraction": shared_fraction,
            "note": so.get("note"),
            "order_cross_check": ("SAMPLE-ORDERING remains a manual "
                                  "cross-check: the shared-sample order "
                                  "must equal VCF sample-header order "
                                  "intersected with the HDF5 accession "
                                  "order (spec 3.4); this audit records "
                                  "the overlap but claims no ordering "
                                  "fact"),
        },
        "orientation_status": "unresolved",
        "orientation_note": NO_DOSAGE_NOTE,
        "verdict": None,
    }


# ---------------------------------------------------------------------------
# gate: PASS vs UNRESOLVED for this HDF5-side half
# ---------------------------------------------------------------------------

_MANDATORY_3_3_KEYS = (
    "hdf5_coordinate_rows",
    "unique_coordinates",
    "unique_matches",
    "multiple_candidates",
    "hdf5_only_coordinates",
    "vcf_only_coordinates",
)


def spot_gate(mapping_summary, fused, spot_result):
    """PASS / UNRESOLVED verdict mapping for the 3.3 + 3.4 step-2
    facts owned by this half.

    Returns {"status": "PASS" | "UNRESOLVED", "reasons": [..],
    "gate": "VCF_SEMANTICS_READY_FOR_CONVERSION" |
    "VCF_SEMANTICS_UNRESOLVED"}.

    PASS requires ALL of:
      - every spec 3.3 quantity computed on `mapping_summary`;
      - NO duplicated `(chr,pos)` rows and NO duplicated
        `(chr,pos,ref,alt)` quadruples in the biallelic catalog
        (a duplicated catalog row alone keeps the gate UNRESOLVED;
        this is by design, not a threshold);
      - the 3.4 spot check was performed (`spot_result` not None);
      - no `suspect_sites` remain after the enumeration;
      - the overall orientation resolves to a label (not `ambiguous`,
        not None);
      - the shared-sample overlap is either fully stated
        (shared_fraction >= 1.0) or undeclared (None -- nothing to
        resolve here); a partial overlap (< 1.0) is an unexplained
        systematic-sample-misalignment candidate -> UNRESOLVED.

    The VCF-side prerequisite facts (source checksum ok, sample
    mapping unambiguous, GT ploidy / ALT-dosage semantics clear,
    reproducible filter policy) belong to audit_1001g_vcf.py and are
    check-required on the VCF gate; they are echoed from that JSON
    when available but not re-decided here.
    """
    reasons = []
    m = mapping_summary or {}
    if not m:
        reasons.append("3.3 mapping summary was not computed")
    else:
        for key in _MANDATORY_3_3_KEYS:
            if m.get(key) is None:
                reasons.append(f"3.3 quantity {key} was not computed")
        if int(m.get("duplicate_chrom_position_pairs") or 0) > 0:
            reasons.append(
                f"duplicated (chr,pos) catalog rows: "
                f"{m.get('duplicate_chrom_position_pairs')} "
                "(duplicated rows -> gate UNRESOLVED by design)")
        if int(m.get("duplicate_chrom_position_ref_alt_quadruples")
               or 0) > 0:
            reasons.append(
                f"duplicated (chr,pos,ref,alt) quadruples: "
                f"{m.get('duplicate_chrom_position_ref_alt_quadruples')}"
                " (duplicated quadruples -> gate UNRESOLVED by design)")

    if spot_result is None:
        reasons.append(
            "the spec-3.4 spot check was not performed "
            "(--per-site-vcf / --hdf5-values join missing): spot "
            "not run -> no coverage of systematic misalignment -> "
            "UNRESOLVED")
    else:
        suspects = spot_result.get("suspect_sites") or []
        if suspects:
            reasons.append(
                f"{len(suspects)} suspect site(s) enumerated by the "
                "3.4 join; systematic sample/coordinate misalignment "
                "not explained -- no threshold-forced pass is "
                "applied (spec 3.4)")
        overall = spot_result.get("overall") or {}
        if int(overall.get("n_pairs_comparable") or 0) == 0:
            reasons.append(
                "no comparable pair in the spot check: the coordinate/"
                "state mapping cannot be certified by 3.4")
        elif overall.get(KEY_BO) == "ambiguous":
            reasons.append(
                "overall orientation is ambiguous: the two concordance "
                "readings tie, so the HDF5 0/1 reference direction is "
                "not resolvable by the spot check alone "
                "(spec 3.4: enumerate, do not threshold, do not "
                "write back)")

    if fused:
        so_block = fused.get("sample_overlap") or {}
        shared_fraction = so_block.get("shared_fraction")
        if shared_fraction is not None and shared_fraction < 1.0:
            reasons.append(
                f"shared-sample overlap declared partial "
                f"({shared_fraction:.4f} < 1.0): the shared-sample "
                "alignment cannot be fully certified; UNRESOLVED "
                "until the shared-sample ordering is independently "
                "verified (SAMPLE-ORDERING cross-check)")

    status = "UNRESOLVED" if reasons else "PASS"
    verdict = (
        "VCF_SEMANTICS_READY_FOR_CONVERSION"
        if not reasons else "VCF_SEMANTICS_UNRESOLVED")
    return {"status": status, "reasons": reasons, "gate": verdict}


# ---------------------------------------------------------------------------
# fixtures / real-file helpers
# ---------------------------------------------------------------------------

def synthetic_hdf5_fixture():
    """The deterministic 5-region mini fixture that the CLI's
    `--hdf5 synthetic` uses in tests: 26 positions across 5
    consecutive regions (one chrom each, ~30 rows as the delivery
    spec requests), with a single duplicated position per
    duplicated region to try the 3.3 duplicate stats.  No h5py,
    no HDF5 binary, no VCF file involved."""
    import numpy as np
    regions = [
        [100, 100, 150, 200, 250, 300, 350, 400],      # chrom "1"
        [50, 51, 60, 70, 80, 85],                      # chrom "2"
        [333, 300, 300, 700],                          # chrom "3"
        [80, 81, 82, 83],                              # chrom "4"
        [900, 901, 902, 903],                          # chrom "5"
    ]
    positions = np.array(
        [value for region in regions for value in region],
        dtype=np.int32)
    spans = []
    start = 0
    for region in regions:
        stop = start + len(region)
        spans.append([start, stop])
        start = stop
    chroms = [str(i + 1) for i in range(len(regions))]
    return {"positions": positions,
            "region_spans": spans,
            "chroms": chroms,
            "n_regions": len(regions),
            "n_positions": int(len(positions)),
            "synthetic": True}


def load_hdf5_positions(path):
    """Read the positions axis + the `chr_regions` spans + `chrs`
    labels from a real HDF5 file; requires the `data` extra (h5py).
    Opaque error when h5py is not installed (documents the `uv sync
    --extra data` requirement)."""
    try:
        import h5py
    except ImportError as exc:
        raise RuntimeError(
            "real-HDF5 loading needs the `data` extra (h5py); install "
            "with `uv sync --extra data`, or use --hdf5 synthetic for "
            "the test fixture") from exc
    with h5py.File(str(path), "r") as handle:
        positions = handle["positions"][:]
        attrs = dict(handle["positions"].attrs)
    # The attribute values come off h5py as numpy objects; the loader
    # must convert to plain Python without ever truth-testing them
    # (`attrs.get("chr_regions") or []` raises ValueError on a
    # multi-element ndarray).  chrs elements arrive as fixed-width
    # bytes (|S5 in the release file) and must be decoded, not
    # str()-ed (str(b"1") is "b'1'", which would corrupt the labels).
    spans, chrs = _spans_chrs_from_attrs(attrs)
    return {"positions": positions, "region_spans": spans,
            "chroms": chrs, "attrs": attrs}


def _spans_chrs_from_attrs(attrs):
    """The region spans + chrom labels off the positions-axis attrs
    dict, in plain Python.  `attrs` is the raw `dict(handle[...].attrs)`
    (numpy values in production, plain lists in tests).  Both keys may
    be absent: spans fall back to empty, chrom labels to the leading
    NUCLEAR labels up to the region count."""
    regions = attrs.get("chr_regions")
    chrs_attr = attrs.get("chrs")
    spans = ([list(int(v) for v in pair) for pair in regions]
             if regions is not None else [])
    if chrs_attr is not None:
        chroms = [label.decode("utf-8") if isinstance(label, (bytes, bytearray))
                  else str(label)
                  for label in chrs_attr]
    else:
        chroms = [NUCLEAR[i] if i < len(NUCLEAR) else str(i + 1)
                   for i in range(max(1, len(spans)))]
    return spans, chroms
# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _write_json(path, payload) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False)
                    + "\n", encoding="utf-8")
    return path


def main(argv=None) -> int:
    """The 3.3 + 3.4-step-2 run over the chosen inputs; writes the two
    contract JSONs under `--out` and prints a single-line machine
    summary.  Exit codes: 0 PASS / 3 UNRESOLVED / 4 bad-or-missing
    input (printed as a one-line JSON error, no traceback)."""
    parser = argparse.ArgumentParser(
        description="TASK03 HDF5-side 3.3 position mapping + 3.4 "
                    "spec-concordance spot-check. No HDF5 write-back, "
                    "no GT-dose conversion, the 18 GB population VCF "
                    "is not read here (see audit_1001g_vcf.py for that "
                    "half).")
    parser.add_argument("--hdf5", required=True,
                        help="Path to a HDF5 file carrying the "
                             "positions axis, or the literal "
                             "`synthetic` for the built-in 5-region "
                             "mini fixture (test path only; no h5py, "
                             "no file access)")
    parser.add_argument("--catalog", required=True,
                        help="The task-02 3.2 biallelic-SNP catalog "
                             "tsv.gz (audit_1001g_vcf output)")
    parser.add_argument("--per-site-vcf", default=None,
                        help="The task-01 VCF-side per-site JSON "
                             "1001g_v31_concordance_spot_check.json "
                             "(audit_1001g_vcf's 3.4 step-1 output)")
    parser.add_argument("--hdf5-values", default=None,
                        help="HDF5-side per-site 0/1 state JSON "
                             "(in production this is a datasource "
                             "join; in tests a tiny fixture file). "
                             "Needed to run the 3.4 step-2 join "
                             "together with --per-site-vcf")
    parser.add_argument("--sample-file", default=None,
                        help="Documented manual shared-sample "
                             "context JSON: {shared_accessions?, "
                             "shared_n?, n_samples_vcf?, "
                             "n_samples_hdf5?, note?}. The pure "
                             "fuse join accepts it as a "
                             "sample-overlap input; it is "
                             "declared, not inferred (the "
                             "SAMPLE-ORDERING cross-check "
                             "remains manual)")
    parser.add_argument("--out", required=True,
                        help="Output directory for the two JSON "
                             "artifacts (created if absent)")
    parser.add_argument("--region-spans", default=None,
                        help="Optional JSON list-of-pairs "
                             "[start, stop] region-span override "
                             "(takes precedence over the manifest/"
                             "HDF5 attrs)")
    parser.add_argument("--hdf5-manifest", default=None,
                        help="Optional manifest JSON "
                             "(data/manifests/1001g_hdf5_*.json) "
                             "holding the region spans + chrs labels "
                             "for the HDF5 positions axis")
    parser.add_argument("--orientation-pass", default=None,
                        help="Optional orientation_pass label: "
                             "hdf5_0_is_ref / hdf5_1_is_ref. "
                             "Passed to the 3.4 join as a "
                             "diagnostic hypothesis only; nothing "
                             "is ever written back to the HDF5 "
                             "file")
    args = parser.parse_args(argv)
    try:
        return _run_cli(args)
    except (OSError, ValueError, RuntimeError) as exc:
        # Unreadable / malformed / h5py-missing case: one-line JSON
        # error, exit 4 (distinct from the gate's 0/3).
        print(json.dumps({
            "module": "at_pheno_audit_hdf5_vcf_mapping",
            "error": f"{type(exc).__name__}: {exc}",
            "exit_code": 4,
        }, sort_keys=True))
        return 4


def _run_cli(args) -> int:
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- positions axis + region spans -----------------------
    manifest_provenance = None
    if args.hdf5 == "synthetic":
        fixture = synthetic_hdf5_fixture()
        positions = fixture["positions"]
        spans = fixture["region_spans"]
        chroms = fixture["chroms"]
        source_note = "synthetic 5-region mini fixture (test path)"
    else:
        loaded = load_hdf5_positions(args.hdf5)
        positions = loaded["positions"]
        spans = loaded["region_spans"]
        chroms = loaded["chroms"]
        source_note = str(args.hdf5)

    # --hdf5-manifest: record provenance (sha256 / byte size) and
    # take the region spans + chrom labels from it if no
    # --region-spans override was supplied.
    if args.hdf5_manifest:
        manifest = _load_json(args.hdf5_manifest)
        manifest_provenance = {
            "path": str(args.hdf5_manifest),
            "sha256_recorded": manifest.get("sha256"),
            "bytes_recorded": manifest.get("bytes"),
            "positions_n_recorded": manifest.get("positions_n"),
        }
        pos_attrs = (manifest.get("schema", {})
                    .get("positions", {}).get("attributes", {}))
        if args.region_spans is None:
            if pos_attrs.get("chr_regions"):
                spans = [list(p) for p in pos_attrs["chr_regions"]]
            if pos_attrs.get("chrs"):
                chroms = [str(c) for c in pos_attrs["chrs"]]
    # --region-spans always beats everything else
    if args.region_spans:
        spans = [list(p) for p in json.loads(args.region_spans)]
        chroms = [str(i + 1) for i in range(len(spans))]
    if not spans:
        raise ValueError(
            "no region spans: pass --region-spans or --hdf5-manifest "
            "(or a real HDF5 file carrying the positions attrs)")
    if not chroms:
        chroms = [str(i + 1) for i in range(len(spans))]
    if len(chroms) != len(spans):
        # Re-derive chunk labels over the spans on mismatch.
        chroms = [str(i + 1) for i in range(len(spans))]
    # The pure core never trusts a span that out-runs the axis
    spans = [[int(s), int(e)] for s, e in spans]
    checked_spans = validate_region_spans(spans, chroms, int(len(positions)))

    # ---- 3.3 mapping ----------------------------------------
    mapping = mapping_from_catalog(read_catalog(args.catalog),
                                   positions, checked_spans,
                                   chroms=chroms)
    catalog_provenance = {
        "path": str(args.catalog),
        "compressed_artifact_sha256": file_sha256(args.catalog),
        "logical_content_md5": catalog_content_md5(args.catalog),
        "catalog_rows": mapping.get("catalog_rows"),
        "duplicate_chrom_position_pairs":
            mapping.get("duplicate_chrom_position_pairs"),
        "duplicate_chrom_position_ref_alt_quadruples":
            mapping.get("duplicate_chrom_position_ref_alt_quadruples"),
        "note": ("provenance recorded, not re-verified against the "
                 "VCF-side delivery notes (3.3/3.4 cross-check "
                 "residues); the biallelic-SNP policy set only. "
                 "Hash distinction: compressed_artifact_sha256 is the "
                 "byte-identity hash of one .tsv.gz file (gzip header "
                 "carries mtime; verifies transfer of that artifact, "
                 "NOT content identity). logical_content_md5 is the "
                 "MD5 of the unpacked TSV stream (stable across "
                 "re-gzips of the same content; the hash to compare "
                 "audit runs against)."),
    }

    # ---- 3.4 step-2 join ------------------------------------
    vcf_side = _load_json(args.per_site_vcf) if args.per_site_vcf else None
    h5_vals = _load_json(args.hdf5_values) if args.hdf5_values else None
    sample_overlap = _load_json(args.sample_file) if args.sample_file else None
    fused = fuse_site_mapping(mapping, vcf_side, sample_overlap)
    spot = None
    if vcf_side is not None and h5_vals is not None:
        spot = reconcile_spot(h5_vals, vcf_side,
                              orientation_pass=args.orientation_pass,
                              sample_overlap=sample_overlap)
    gated = spot_gate(mapping, fused, spot)

    # ---- write the two contract JSONs -----------------------
    mapping_payload = {
        "series": "hdf5_vcf_position_mapping_summary",
        "inputs": {
            "hdf5": source_note,
            "catalog": str(args.catalog),
            "per_site_vcf": args.per_site_vcf,
            "hdf5_values": args.hdf5_values,
            "sample_file": args.sample_file,
            "orientation_pass": args.orientation_pass,
        },
        "hdf5_provenance": manifest_provenance,
        "catalog_provenance": catalog_provenance,
        "position_mapping": mapping,
        "fused": fused,
        "gate": gated,
    }
    _write_json(out_dir / "hdf5_vcf_position_mapping_summary.json",
                mapping_payload)

    if spot is not None:
        spot_payload = {
            "series": "hdf5_vcf_concordance_spot_check",
            "inputs": {
                "per_site_vcf": args.per_site_vcf,
                "hdf5_values": args.hdf5_values,
                "sample_file": args.sample_file,
                "orientation_pass": args.orientation_pass,
            },
            "hdf5_provenance": manifest_provenance,
            "spot_check": spot,
            "fused": fused,
            "gate": gated,
        }
    else:
        spot_payload = {
            "series": "hdf5_vcf_concordance_spot_check",
            "inputs": {
                "per_site_vcf": args.per_site_vcf,
                "hdf5_values": args.hdf5_values,
                "sample_file": args.sample_file,
                "orientation_pass": args.orientation_pass,
            },
            "hdf5_provenance": manifest_provenance,
            "spot_check": None,
            "spot_check_unrun_reason": ("3.4 step-2 join needs both "
                                        "--per-site-vcf and "
                                        "--hdf5-values; provide the "
                                        "datasource-join per-site 0/1 "
                                        "values to run it"),
            "fused": fused,
            "gate": gated,
        }
    _write_json(out_dir / "hdf5_vcf_concordance_spot_check.json",
                spot_payload)

    # ---- one-line machine summary + exit code ----------------
    one_line = {
        "module": "at_pheno_audit_hdf5_vcf_mapping",
        "status": gated["status"],
        "gate": gated["gate"],
        "exit_code": 0 if gated["status"] == "PASS" else 3,
        "hdf5_coordinate_rows": mapping.get("hdf5_coordinate_rows"),
        "unique_matches": mapping.get("unique_matches"),
        "multiple_candidates": mapping.get("multiple_candidates"),
        "hdf5_only": mapping.get("hdf5_only_coordinates"),
        "vcf_only": mapping.get("vcf_only_coordinates"),
        "coord_hits": mapping.get("coord_hits"),
        "dup_pairs": mapping.get("duplicate_chrom_position_pairs"),
        "dup_quads":
            mapping.get("duplicate_chrom_position_ref_alt_quadruples"),
        "site_pool_capped": mapping.get("site_pool", {}).get(
            "n_pool_capped"),
        "n_sites": spot["n_sites"] if spot else None,
        "n_suspects": (len(spot.get("suspect_sites") or [])
                       if spot else None),
        "overall_orientation": (spot.get("overall", {})
                                .get(KEY_BO) if spot else None),
        "overall_best_concordance": (spot.get("overall", {})
                                     .get(KEY_BC) if spot else None),
        "reasons": gated["reasons"],
        "out_dir": str(out_dir),
    }
    print(json.dumps(one_line, sort_keys=True))
    return 0 if gated["status"] == "PASS" else 3


if __name__ == "__main__":
    sys.exit(main())
