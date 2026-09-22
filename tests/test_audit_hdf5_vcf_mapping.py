"""Fixture tests for the TASK03 HDF5-side 3.3/3.4 audit
(scripts/audit_hdf5_vcf_mapping.py).

Self-made fixtures under tmp_path only (a tiny catalog gz, the
module's built-in 5-region mini positions fixture, small per-site
JSON files): the 10.7M-row HDF5 positions axis and the 18 GB
population VCF are never touched.  Coverage:

- 3.3 quantities on the mini fixture: total rows / unique /
  duplicated positions (per-chrom), unique matches, multi-candidate
  coordinates, HDF5-only vs VCF-only coordinates, duplicate
  (chr,pos) pairs and (chr,pos,ref,alt) quads on the catalog, the
  policy-skipped rows that prove a multiallelic candidate never
  participates, and `coord_hits`;
- region-span validation errors;
- the deterministic fixed-hash/salt site pool (cap + salt);
- 3.4 step-2 join: the two-orientation concordance readings (clearly
  different for the known-better orientation), the tie case ->
  "ambiguous", suspect-site enumeration as a LIST (no threshold
  verdict), counts-only sites, and the every-site key presence;
- fuse (pool coverage + shared-sample overlap) and the gate verdict
  mapping (PASS; duplicated rows alone -> UNRESOLVED; no-spot ->
  UNRESOLVED; ambiguous/suspects -> UNRESOLVED; partial overlap ->
  UNRESOLVED);
- CLI round-trip through `main()`: the exit-0 / exit-3 / exit-4
  paths, the two JSON artifacts, and the exact spec keys in the
  emitted per-site / overall dicts (one-line machine summary
  captured on stdout).
"""

import gzip
import importlib.util
import json
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
_MOD_PATH = _REPO / "scripts" / "audit_hdf5_vcf_mapping.py"
_MOD = None


def _module():
    global _MOD
    if _MOD is None:
        spec = importlib.util.spec_from_file_location(
            "audit_hdf5_vcf_mapping", _MOD_PATH)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _MOD = module
    return _MOD


# ---------------------------------------------------------------- fixtures

FIXTURE_SPANS = [[0, 8], [8, 14], [14, 18], [18, 22], [22, 26]]
FIXTURE_CHROMS = ["1", "2", "3", "4", "5"]
# 26 positions: chrom-1 has a duplicated 100, chrom-3 a duplicated 300
FIXTURE_POSITIONS = [
    # region 0 (chrom 1, 8 rows)
    100, 100, 150, 200, 250, 300, 350, 400,
    # region 1 (chrom 2, 6 rows)
    50, 51, 60, 70, 80, 85,
    # region 2 (chrom 3, 4 rows)
    333, 300, 300, 700,
    # region 3 (chrom 4, 4 rows)
    80, 81, 82, 83,
    # region 4 (chrom 5, 4 rows)
    900, 901, 902, 903,
]

# The clean catalog: 6 valid biallelic rows + 1 multiallelic-style
# row that the 3.1 policy must exclude.  (2,51)/(4,81)/(5,900) are
# unique matches; (3,300) is only reachable through ONE valid row
# (the skipped row must not add a second candidate); (5,4242) is
# VCF-only (not in the positions axis).
CLEAN_ROWS = [
    ("1", 100, "A", "C", 1),
    ("3", 300, "C", "G", 2),
    ("5", 4242, "T", "C", 3),
    ("3", 300, "A", "C,G", 4),   # policy-skipped: 2-allele ALT
    ("2", 51, "G", "C", 5),
    ("4", 81, "C", "A", 6),
    ("5", 900, "G", "C", 7),
]
# The duplicated-row variant (spec 3.2 duplicated-(chr,pos) proof):
# two additional records at the same coordinate -> a duplicated pair
# AND a duplicated (ref,alt) quadruple.
DUP_ROWS = CLEAN_ROWS + [
    ("1", 100, "A", "C", 8),
    ("1", 100, "G", "T", 9),
]

CLEAN_EXPECT = {
    "hdf5_coordinate_rows": 26,
    "unique_coordinates": 24,
    "duplicated_positions": 2,          # (1,100) and (3,300)
    "duplicated_position_rows": 2,
    "catalog_rows": 6,
    "unique_matches": 5,
    "multiple_candidates": 0,
    "hdf5_only_coordinates": 19,
    "vcf_only_coordinates": 1,          # (5,4242)
    "duplicate_chrom_position_pairs": 0,
    "duplicate_chrom_position_ref_alt_quadruples": 0,
    "policy_skipped_rows": 1,
    "coord_hits": 5,
}
DUP_EXPECT = {
    **CLEAN_EXPECT,
    "catalog_rows": 8,
    "unique_matches": 4,                # (1,100) is no longer unique
    "multiple_candidates": 1,
    "duplicate_chrom_position_pairs": 1,
    "duplicate_chrom_position_ref_alt_quadruples": 1,
    "coord_hits": 5,                   # hits unchanged
}


def write_catalog(tmp, rows, name="catalog.gz"):
    path = tmp / name
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("chromosome\tposition\tref\talt\t"
                     "source_record_index\n")
        for chrom, pos, ref, alt, idx in rows:
            handle.write(f"{chrom}\t{pos}\t{ref}\t{alt}\t{idx}\n")
    return path


def spot_entry(chrom, pos, n_hom0, n_het, n_hom1, n_missing,
               per_sample, record_found=True):
    """A VCF-side 3.4-step-1 per-site entry: the spec-mandated keys,
    the `states` block, and the (datasource-join) `per_sample`
    re-search alignment codes."""
    return {
        "chrom": chrom,
        "pos": pos,
        "record_found": record_found,
        "concordance_if_hdf5_0_is_REF": None,
        "concordance_if_hdf5_1_is_REF": None,
        "best_orientation": None,
        "best_concordance": None,
        "n_comparable": n_hom0 + n_hom1,
        "n_heterozygous": n_het,
        "n_missing_vcf": n_missing,
        "n_states_out_of_gt": 0,
        "states": {"homozygous_0": n_hom0, "heterozygous": n_het,
                   "homozygous_1": n_hom1, "missing": n_missing,
                   "out_of_range_or_unexpected": 0},
        "per_sample": per_sample,
    }


def make_vcf_side_json(entries):
    return {"n_sites": len(entries), "n_samples": 8, "per_site": entries}


# ------------------------------------------------- 3.3 mapping ---

def test_mapping_clean_3_3_quantities():
    mod = _module()
    fixture = mod.synthetic_hdf5_fixture()
    mapping = mod.mapping_from_catalog(iter(CLEAN_ROWS),
                                       fixture["positions"],
                                       fixture["region_spans"])
    for key, want in CLEAN_EXPECT.items():
        assert mapping[key] == want, key
    assert mapping["chroms"] == FIXTURE_CHROMS
    assert mapping["n_regions"] == 5
    assert mapping["duplicated_positions_per_chrom"] == \
        {"1": 1, "2": 0, "3": 1, "4": 0, "5": 0}
    assert mapping["duplicated_coordinate_examples"] == \
        [["1", "100"], ["3", "300"]]
    # bounded-memory claim: coord_hits is stated from the dict,
    # not from re-stored rows
    assert mapping["coord_hits"] == mapping["unique_matches"] + \
        mapping["multiple_candidates"]
    assert mapping["no_orientation_claimed"] is True
    assert mapping["site_pool"]["n_pool_capped"] == 5
    # the multiallelic-shaped row proves non-participation (3.1):
    # it adds no candidate anywhere; (3,300) remains a UNIQUE match
    # and multiple_candidates stays at 0
    assert mapping["policy_skipped_examples"] == [["3", "300", "A", "C,G"]]
    assert ["3", "300"] in mapping["unique_matches_sample"]
    assert mapping["multiple_candidates_sample"] == []
    assert mapping["vcf_only_sample"] == [["5", "4242"]]


def test_mapping_duplicated_rows_surface_in_3_3():
    mod = _module()
    fixture = mod.synthetic_hdf5_fixture()
    mapping = mod.mapping_from_catalog(iter(DUP_ROWS),
                                       fixture["positions"],
                                       fixture["region_spans"])
    for key, want in DUP_EXPECT.items():
        assert mapping[key] == want, key
    assert mapping["multiple_candidates_sample"] == [
        {"chrom": "1", "pos": "100", "n_catalog_candidates": 3}]


def test_region_span_validation():
    mod = _module()
    fixture = mod.synthetic_hdf5_fixture()
    n = len(fixture["positions"])
    good = mod.validate_region_spans(FIXTURE_SPANS, FIXTURE_CHROMS, n)
    assert good == [(0, 8), (8, 14), (14, 18), (18, 22), (22, 26)]
    # two consecutive runs do not tile
    try:
        mod.validate_region_spans(
            [[0, 8], [9, 14], [14, 18], [18, 22], [22, 26]],
            FIXTURE_CHROMS, n)
        raise AssertionError("expected a ValueError on a gap")
    except ValueError as exc:
        assert "consecutive" in str(exc)
    # a span runs past the end
    try:
        mod.validate_region_spans(
            [[0, 8], [8, 14], [14, 18], [18, 22], [22, 30]],
            FIXTURE_CHROMS, n)
        raise AssertionError("expected a ValueError on overshoot")
    except ValueError as exc:
        assert "positions axis" in str(exc)
    # an empty run is inverted
    try:
        mod.validate_region_spans(
            [[0, 8], [8, 8], [8, 18], [18, 22], [22, 26]],
            FIXTURE_CHROMS, n)
        raise AssertionError("expected a ValueError on an empty run")
    except ValueError as exc:
        assert "empty or inverted" in str(exc)
    # the chrom-label count must match the run count
    try:
        mod.validate_region_spans(FIXTURE_SPANS, ["1", "2"], n)
        raise AssertionError("expected a ValueError on label count")
    except ValueError as exc:
        assert "differs" in str(exc)


def test_synthetic_fixture_shape():
    mod = _module()
    fixture = mod.synthetic_hdf5_fixture()
    assert fixture["n_regions"] == 5 == len(fixture["region_spans"])
    assert len(fixture["positions"]) == 26
    assert fixture["chroms"] == FIXTURE_CHROMS
    # the duplicated spots: chrom-1 region and chrom-3 region each
    # carry one duplicated position
    seg1 = [int(v) for v in fixture["positions"][0:8]]
    seg3 = [int(v) for v in fixture["positions"][14:18]]
    assert seg1.count(100) == 2
    assert seg3.count(300) == 2
    # the spans tile the axis
    start = 0
    for s, e in fixture["region_spans"]:
        assert s == start
        start = e
    assert start == len(fixture["positions"])


def test_selected_sites_deterministic_capped():
    mod = _module()
    input_coords = [("1", 100), ("3", 300), ("5", 900), ("4", 81),
                    ("2", 51)]
    m1 = mod.selected_sites(input_coords)
    m2 = mod.selected_sites(input_coords)
    assert m1 == m2                      # deterministic across calls
    assert len(m1) == 5
    for pair in m1:
        assert len(pair) == 2 and isinstance(pair[0], str)
    assert sorted(tuple(p) for p in m1) == \
        [("1", "100"), ("2", "51"), ("3", "300"), ("4", "81"),
         ("5", "900")]
    capped = mod.selected_sites(input_coords, cap=3)
    assert len(capped) == 3
    assert set(map(tuple, capped)).issubset(set(map(tuple, m1)))
    other_salt = mod.selected_sites(input_coords, cap=3,
                                     salt="second-salt")
    again = mod.selected_sites(input_coords, cap=3,
                               salt="second-salt")
    assert other_salt == again           # stable for a fixed salt
    assert len(other_salt) == 3


# ------------------------------------------------- 3.4 step 2 ---

def _all_site_values():
    """Three-site clean fixture: the two orientation hypotheses are
    clearly different at (2,51) and (4,81); (5,900) is the tie case."""
    h5_vals = {
        "2": {"51": [0, 0, 1, None, 0, 1, 0, 1]},
        "4": {"81": [0, 1, 1, 0, 0, 1, 0, 0]},
        "5": {"900": [0, 0, 0, 0, 0, 0, 0, 0]},
    }
    vcf_side = [
        spot_entry("2", 51, 4, 1, 1, 1,
                   ["0", "0", "1", ".", "0", "1", "0", "H"]),
        spot_entry("4", 81, 3, 0, 3, 0,
                   ["0", "0", "1", "1", "0", "1"]),
        spot_entry("5", 900, 3, 0, 3, 0,
                   ["0", "1", "0", "1", "0", "1"]),
    ]
    return h5_vals, vcf_side


def test_reconcile_spot_two_orientations_and_tie():
    mod = _module()
    h5_vals, vcf_side = _all_site_values()
    spot = mod.reconcile_spot(h5_vals, vcf_side, orientation_pass=None,
                              sample_overlap={"shared_n": 8,
                                              "n_samples_vcf": 8,
                                              "n_samples_hdf5": 8})
    assert spot["n_sites"] == 3
    assert spot["overall"]["n_sites_with_pairs"] == 3
    by_site = {(e["chrom"], e["pos"]): e for e in spot["per_site"]}
    e1, e2, e3 = by_site[("2", 51)], by_site[("4", 81)], \
        by_site[("5", 900)]
    # every per-site dict carries all the spec-mandated keys
    required = set(mod.SPEC_KEYS)
    for entry in spot["per_site"]:
        assert required.issubset(entry.keys()), entry
    # (2,51): the state-0-is-REF orientation clearly wins
    assert e1["n_pairs_comparable"] == 6
    assert e1["concordance_if_hdf5_0_is_REF"] == 1.0
    assert e1["concordance_if_hdf5_1_is_REF"] == 0.0
    assert e1["best_orientation"] == "hdf5_0_is_ref"
    assert e1["best_concordance"] == 1.0
    assert e1["n_comparable"] == 6
    assert e1["n_heterozygous"] == 1
    assert e1["n_missing_vcf"] == 1
    assert e1["n_hdf5_missing"] == 1      # the shared sample with no h5
    assert e1["n_hdf5_states"] == 8
    assert e1["suspect"] is False
    # (4,81): state-0-is-REF wins 4/6 vs 2/6
    assert e2["n_pairs_comparable"] == 6
    assert e2["concordance_if_hdf5_0_is_REF"] == 4 / 6
    assert e2["concordance_if_hdf5_1_is_REF"] == 2 / 6
    assert e2["best_orientation"] == "hdf5_0_is_ref"
    assert e2["best_concordance"] == 4 / 6
    # (5,900): tie -> the orientation is ambiguous, NOT a conclusion
    assert e3["concordance_if_hdf5_0_is_REF"] == 0.5
    assert e3["concordance_if_hdf5_1_is_REF"] == 0.5
    assert e3["best_orientation"] == "ambiguous"
    assert e3["best_concordance"] == 0.5
    # the tie site is the only enumerated suspect -- a LIST, not an
    # arbitrary verdict; there is no threshold-based forced pass
    assert len(spot["suspect_sites"]) == 1
    assert spot["suspect_sites"][0]["pos"] == 900
    assert spot["suspect_sites"][0]["suspect"] is True
    assert spot["suspect_checklist"]  # the cross-check items are present
    # overall reading: 13 of 18 comparable pairs in state 0 -> ref
    assert spot["overall"] == {
        "concordance_if_hdf5_0_is_REF": 13 / 18,
        "concordance_if_hdf5_1_is_REF": 5 / 18,
        "best_orientation": "hdf5_0_is_ref",
        "best_concordance": 13 / 18,
        "n_comparable": 18,
        "n_heterozygous": 1,
        "n_missing_vcf": 1,
        "n_pairs_comparable": 18,
        "n_sites_with_pairs": 3,
    }
    # the explicit no-dosage note must be present
    assert "ALT" in spot["no_dosage_note"]
    assert "written back" in spot["no_dosage_note"]


def test_state_code_named_tokens_equivalent_to_symbol_tokens():
    """The VCF side emits the STATE-CODE NAMES (homozygous_0,
    heterozygous, missing, ...) as `per_sample` entries.  The 3.4
    join must treat those names exactly like the symbol tokens --
    otherwise the whole shared-sample pair set empties out (the
    exact regression observed on the 1135 x 2000 pool scan, which
    reported n_pairs_comparable = 0 while both sides were fully
    populated)."""
    mod = _module()
    # Token-level behaviour: every state-code NAME in the module's
    # token table must round-trip through _vcf_state_code, and a
    # foreign token must land out_of_range_or_unexpected.
    state_names = sorted(set(mod._STATE_TOKENS.values()))
    for name in state_names:
        assert mod._vcf_state_code(name) == name, name
    assert mod._vcf_state_code("zebra") == "out_of_range_or_unexpected"

    h5_vals, vcf_side = _all_site_values()
    # symbol tokens -> their state-code names, from the module's own
    # table (no hand-transcribed literals)
    symbol_tokens = ("0", "1", "H", ".")
    to_name = {sym: mod._STATE_TOKENS[sym] for sym in symbol_tokens}
    named_side = [dict(entry, per_sample=[
                        to_name.get(t, t)
                        for t in (entry["per_sample"] or [])])
                  for entry in vcf_side]
    from_symbols = mod.reconcile_spot(h5_vals, vcf_side)
    from_named = mod.reconcile_spot(h5_vals, named_side)
    # the named-code list must pair EXACTLY the same sample pairs --
    # identical overall block, per-site orientations, and the
    # counts re-derived from `per_sample` itself
    assert from_symbols["overall"] == from_named["overall"]
    assert from_named["overall"]["n_pairs_comparable"] == 18
    assert from_named["overall"]["best_orientation"] == "hdf5_0_is_ref"
    assert [e["chrom"] + "/" + str(e["pos"]) for e in from_named["per_site"]] \
        == [e["chrom"] + "/" + str(e["pos"]) for e in from_symbols["per_site"]]
    assert [e["suspect"] for e in from_named["per_site"]] == \
        [e["suspect"] for e in from_symbols["per_site"]]


def test_reconcile_spot_suspect_enumeration_and_orientation_pass():
    mod = _module()
    h5_vals, vcf_side = _all_site_values()
    # add three more sites: D (hdf5 side is missing), E (no hdf5
    # side at all), F (out-of-range h5 states)
    h5_vals["5"] = dict(h5_vals["5"],
                        **{"33": [None, None],
                           "442": [0, "x", 7]})
    vcf_side += [
        spot_entry("5", 33, 2, 0, 0, 0, ["0", "0"]),
        spot_entry("5", 44, 2, 0, 0, 0, ["0", "0"]),
        spot_entry("5", 442, 4, 0, 0, 0, ["0", "0", "0", "0"]),
    ]
    # without an orientation_pass: D (no comparable pair), E (no h5
    # side), 900 (tie); F has a pair -> not a suspect
    spot = mod.reconcile_spot(h5_vals, vcf_side, orientation_pass=None)
    suspects = {(e["chrom"], e["pos"]) for e in spot["suspect_sites"]}
    assert suspects == {("5", 33), ("5", 44), ("5", 900)}
    by_pos = {e["pos"]: e for e in spot["per_site"]}
    e_d, e_e, e_f = by_pos[33], by_pos[44], by_pos[442]
    assert "no comparable pair" in " | ".join(e_d["suspect_reasons"])
    assert "no HDF5-side value" in " | ".join(e_e["suspect_reasons"])
    assert e_f["suspect"] is False
    assert e_f["n_hdf5_invalid"] == 2     # "x" and 7 are not coerced
    assert e_f["n_hdf5_missing"] == 0
    assert e_f["n_pairs_comparable"] == 1  # only the "0" line is usable
    assert e_f["concordance_if_hdf5_0_is_REF"] == 1.0
    # overall still resolves even with the suspects
    overall = spot["overall"]
    assert overall["best_orientation"] == "hdf5_0_is_ref"
    assert overall["n_pairs_comparable"] > 0
    # with orientation_pass=hdf5_1_is_ref the mismatch sites
    # (2,51)/(4,81) are enumerated on top of the existing suspects
    spot_pass = mod.reconcile_spot(h5_vals, vcf_side,
                                   orientation_pass="hdf5_1_is_ref",
                                   sample_overlap={})
    suspects_pass = {(e["chrom"], e["pos"])
                     for e in spot_pass["suspect_sites"]}
    assert {("2", 51), ("4", 81)} <= suspects_pass
    e1 = next(e for e in spot_pass["per_site"] if e["pos"] == 51)
    assert any("orientation_pass" in reason
               for reason in e1["suspect_reasons"])
    assert spot_pass["orientation_pass"] == "hdf5_1_is_ref"
    # an unknown orientation_pass label is rejected
    try:
        mod.reconcile_spot(h5_vals, vcf_side,
                           orientation_pass="state_2_is_ref")
        raise AssertionError("expected a ValueError for that label")
    except ValueError as exc:
        assert "orientation_pass" in str(exc)


def test_reconcile_spot_counts_only_site():
    mod = _module()
    vcf_only = [spot_entry("5", 44, 2, 0, 0, 0, ["0", "0"])]
    spot = mod.reconcile_spot({}, vcf_only, sample_overlap=None)
    assert spot["n_sites"] == 1
    entry = spot["per_site"][0]
    # a VCF-only site: no 0/1 map, all spec keys are still present
    # with null/zero values
    assert entry["hdf5_side_present"] is False
    assert entry["vcf_side_present"] is True
    assert entry["record_found"] is True
    assert entry["concordance_if_hdf5_0_is_REF"] is None
    assert entry["concordance_if_hdf5_1_is_REF"] is None
    assert entry["best_orientation"] is None
    assert entry["best_concordance"] is None
    assert entry["n_comparable"] == 2
    assert entry["n_heterozygous"] == 0
    assert entry["n_missing_vcf"] == 0
    assert entry["suspect"] is True
    assert spot["suspect_sites"][0]["pos"] == 44
    # overall: no comparable pair -> the orientation stays None
    assert spot["overall"]["best_orientation"] is None
    assert spot["overall"]["n_pairs_comparable"] == 0


def test_index_vcf_sites_shapes():
    mod = _module()
    entry_a = {"chrom": "2", "pos": 51, "record_found": True}
    entry_b = {"chrom": "4", "pos": 81, "record_found": True}
    whole_json = {"n_sites": 2, "per_site": [entry_a, entry_b]}
    assert set(mod.index_vcf_sites(whole_json)) == {("2", 51), ("4", 81)}
    assert set(mod.index_vcf_sites([entry_a, entry_b])) == \
        {("2", 51), ("4", 81)}
    assert set(mod.index_vcf_sites({"2": {"51": entry_a},
                                   "4": {"81": entry_b}})) == \
        {("2", 51), ("4", 81)}
    assert set(mod.index_vcf_sites({"2": [entry_a], "4": [entry_b]})) == \
        {("2", 51), ("4", 81)}


def test_index_site_values_shapes():
    mod = _module()
    nested = {"2": {"51": [0, 1]}, "4": {"81": [None]}}
    flat = {("2", 51): [0], ("4", 81): [1, 1]}
    entries = [{"chrom": "5", "pos": 90, "values": [1, 1]}]
    assert mod.index_site_values(nested) == {("2", 51): [0, 1],
                                            ("4", 81): [None]}
    assert mod.index_site_values(flat) == {("2", 51): [0],
                                          ("4", 81): [1, 1]}
    assert mod.index_site_values(entries) == {("5", 90): [1, 1]}
    assert mod.index_site_values({}) == {}
    assert mod.index_site_values(None) == {}


# ------------------------------------------------- fuse + gate ---

def test_fuse_site_mapping_pool_and_overlap():
    mod = _module()
    fixture = mod.synthetic_hdf5_fixture()
    mapping = mod.mapping_from_catalog(iter(CLEAN_ROWS),
                                       fixture["positions"],
                                       fixture["region_spans"])
    vcf_side = [
        spot_entry("2", 51, 4, 1, 1, 1,
                   ["0", "0", "1", ".", "0", "1", "0", "H"]),
        spot_entry("5", 900, 1, 0, 2, 0, ["1", "1", "0"],
                   record_found=False),
    ]
    fused = mod.fuse_site_mapping(
        mapping, vcf_side,
        {"shared_n": 8, "n_samples_vcf": 16, "n_samples_hdf5": 16})
    assert fused["site_pool"]["n"] == 5
    assert fused["site_pool"]["pool_sites_with_vcf"] == \
        [["2", "51"], ["5", "900"]]
    assert sorted(fused["site_pool"]["pool_sites_missing_in_vcf"]) == \
        [["1", "100"], ["3", "300"], ["4", "81"]]
    assert fused["sample_overlap"]["shared_fraction"] == 0.5
    # a 12/13 overlap stays a declared partial overlap
    fused_partial = mod.fuse_site_mapping(
        mapping, vcf_side,
        {"shared_n": 12, "n_samples_vcf": 13, "n_samples_hdf5": 16})
    assert fused_partial["sample_overlap"]["shared_fraction"] == 12 / 13
    # no shared context -> the fraction stays None
    fused_none = mod.fuse_site_mapping(mapping, None, None)
    assert fused_none["sample_overlap"]["shared_fraction"] is None
    assert fused_none["orientation_status"] == "unresolved"
    assert fused_none["verdict"] is None
    # no vcf side -> the pool search is not runnable
    assert fused_none["site_pool"]["pool_sites_with_vcf"] is None
    assert fused_none["site_pool"]["pool_sites_missing_in_vcf"] is None


def _clean_mapping_dict():
    return {
        "hdf5_coordinate_rows": 26,
        "unique_coordinates": 24,
        "unique_matches": 5,
        "multiple_candidates": 0,
        "hdf5_only_coordinates": 19,
        "vcf_only_coordinates": 1,
        "duplicate_chrom_position_pairs": 0,
        "duplicate_chrom_position_ref_alt_quadruples": 0,
    }


def _clean_fused_block(shared_fraction=None):
    return {"sample_overlap": {"shared_fraction": shared_fraction,
                               "shared_n": None,
                               "n_samples_vcf": None,
                               "n_samples_hdf5": None}}


def _clean_spot():
    return {"suspect_sites": [],
            "overall": {"concordance_if_hdf5_0_is_REF": 0.8,
                        "concordance_if_hdf5_1_is_REF": 0.2,
                        "best_orientation": "hdf5_0_is_ref",
                        "best_concordance": 0.8,
                        "n_comparable": 10,
                        "n_heterozygous": 1,
                        "n_missing_vcf": 1,
                        "n_pairs_comparable": 10}}


def test_spot_gate_pass():
    mod = _module()
    gate = mod.spot_gate(_clean_mapping_dict(),
                         _clean_fused_block(), _clean_spot())
    assert gate["status"] == "PASS"
    assert gate["gate"] == "VCF_SEMANTICS_READY_FOR_CONVERSION"
    assert gate["reasons"] == []


def test_spot_gate_duplicated_rows_alone_unresolved():
    mod = _module()
    mapping = dict(_clean_mapping_dict(),
                   duplicate_chrom_position_pairs=1,
                   duplicate_chrom_position_ref_alt_quadruples=1)
    gate = mod.spot_gate(mapping, _clean_fused_block(), _clean_spot())
    assert gate["status"] == "UNRESOLVED"
    assert len(gate["reasons"]) == 2
    assert "duplicated (chr,pos) catalog rows" in gate["reasons"][0]
    assert "quadruples" in gate["reasons"][1]
    assert gate["gate"] == "VCF_SEMANTICS_UNRESOLVED"


def test_spot_gate_no_spot_check_unresolved():
    mod = _module()
    gate = mod.spot_gate(_clean_mapping_dict(),
                         _clean_fused_block(), None)
    assert gate["status"] == "UNRESOLVED"
    assert len(gate["reasons"]) == 1
    assert "spot check was not performed" in gate["reasons"][0]


def test_spot_gate_ambiguous_and_suspects_unresolved():
    mod = _module()
    ambiguous = {"suspect_sites": [],
                 "overall": {"concordance_if_hdf5_0_is_REF": 0.5,
                             "concordance_if_hdf5_1_is_REF": 0.5,
                             "best_orientation": "ambiguous",
                             "best_concordance": 0.5,
                             "n_pairs_comparable": 4}}
    gate = mod.spot_gate(_clean_mapping_dict(),
                         _clean_fused_block(), ambiguous)
    assert gate["status"] == "UNRESOLVED"
    assert any("ambiguous" in reason for reason in gate["reasons"])
    suspects = {"suspect_sites": [
        {"chrom": "2", "pos": 10}, {"chrom": "3", "pos": 20}],
        "overall": {"concordance_if_hdf5_0_is_REF": 0.7,
                    "concordance_if_hdf5_1_is_REF": 0.3,
                    "best_orientation": "hdf5_0_is_ref",
                    "best_concordance": 0.7,
                    "n_pairs_comparable": 9}}
    gate = mod.spot_gate(_clean_mapping_dict(),
                         _clean_fused_block(), suspects)
    assert gate["status"] == "UNRESOLVED"
    assert any("suspect site" in reason for reason in gate["reasons"])


def test_spot_gate_partial_shared_overlap_unresolved():
    mod = _module()
    gate = mod.spot_gate(_clean_mapping_dict(),
                         _clean_fused_block(0.75), _clean_spot())
    assert gate["status"] == "UNRESOLVED"
    assert any("partial" in reason for reason in gate["reasons"])
    assert any("0.7500" in reason for reason in gate["reasons"])


# ------------------------------------------------- CLI ---

def _cli_input_files(tmp_path, tie=False, dup_catalog=False):
    catalog = write_catalog(tmp_path,
                            DUP_ROWS if dup_catalog else CLEAN_ROWS)
    entries = [
        spot_entry("2", 51, 4, 1, 1, 1,
                   ["0", "0", "1", ".", "0", "1", "0", "H"]),
        spot_entry("4", 81, 3, 0, 3, 0,
                   ["0", "0", "1", "1", "0", "1"]),
    ]
    if tie:
        entries.append(spot_entry("5", 900, 3, 0, 3, 0,
                                  ["0", "1", "0", "1", "0", "1"]))
        h5_values = {"2": {"51": [0, 0, 1, None, 0, 1, 0, 1]},
                     "4": {"81": [0, 1, 1, 0, 0, 1, 0, 0]},
                     "5": {"900": [0] * 8}}
    else:
        entries.append(spot_entry("5", 900, 1, 0, 2, 0, ["1", "1", "0"]))
        h5_values = {"2": {"51": [0, 0, 1, None, 0, 1, 0, 1]},
                     "4": {"81": [0, 1, 1, 0, 0, 1, 0, 0]},
                     "5": {"900": [1, 1, 1, 0, 0, 0, 0, 0]}}
    vcf_path = tmp_path / "per_site_vcf.json"
    vcf_path.write_text(json.dumps(make_vcf_side_json(entries)),
                        encoding="utf-8")
    h5_path = tmp_path / "h5_values.json"
    h5_path.write_text(json.dumps(h5_values), encoding="utf-8")
    return {"catalog": catalog, "per_site_vcf": vcf_path,
            "h5_values": h5_path}


def _run_cli(tmp_path, args):
    """Run main() with stdout captured; return (exit_code, json line)."""
    import contextlib
    import io
    mod = _module()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = mod.main(args)
    return code, json.loads(buf.getvalue())


def test_cli_round_trip_pass(tmp_path):
    out = tmp_path / "out_pass"
    fi = _cli_input_files(tmp_path)
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(fi["catalog"]),
         "--per-site-vcf", str(fi["per_site_vcf"]),
         "--hdf5-values", str(fi["h5_values"]),
         "--out", str(out)])
    assert code == 0
    assert line["status"] == "PASS"
    assert line["gate"] == "VCF_SEMANTICS_READY_FOR_CONVERSION"
    assert line["exit_code"] == 0
    # 3.3 quantities reach the machine line
    assert line["hdf5_coordinate_rows"] == 26
    assert line["unique_matches"] == 5
    assert line["multiple_candidates"] == 0
    assert line["vcf_only"] == 1
    assert line["dup_pairs"] == 0
    assert line["dup_quads"] == 0
    assert line["coord_hits"] == 5
    # spot checks also reach the machine line
    assert line["n_sites"] == 3
    assert line["n_suspects"] == 0
    assert line["overall_orientation"] == "hdf5_0_is_ref"
    # both JSON artifacts exist and carry the exact spec keys
    mapping_payload = json.loads(
        (out / "hdf5_vcf_position_mapping_summary.json").read_text())
    position_mapping = mapping_payload["position_mapping"]
    for key, want in CLEAN_EXPECT.items():
        assert position_mapping[key] == want, key
    assert mapping_payload["gate"]["status"] == "PASS"
    # Goal 5: catalog provenance keeps the two hashes distinct and the
    # bare `sha256` key must not reappear as a "content" claim
    cat = mapping_payload["catalog_provenance"]
    assert "sha256" not in cat, "no bare sha256 key: artifact vs content split"
    content = ("chromosome\tposition\tref\talt\tsource_record_index\n"
               + "".join(f"{c}\t{p}\t{r}\t{a}\t{i}\n"
                         for c, p, r, a, i in CLEAN_ROWS))
    assert cat["logical_content_md5"] == _md5(content)
    assert len(cat["compressed_artifact_sha256"]) == 64
    spot_payload = json.loads(
        (out / "hdf5_vcf_concordance_spot_check.json").read_text())
    spot_check = spot_payload["spot_check"]
    required = set(_module().SPEC_KEYS)
    assert spot_check["n_sites"] == 3
    for entry in spot_check["per_site"]:
        assert required.issubset(entry.keys())
    assert required.issubset(spot_check["overall"].keys())
    assert spot_check["overall"]["n_sites_with_pairs"] == 3
    assert spot_check["suspect_sites"] == []
    assert spot_check["overall"]["best_orientation"] == "hdf5_0_is_ref"
    assert spot_payload["gate"]["gate"] == \
        "VCF_SEMANTICS_READY_FOR_CONVERSION"


def _md5(text):
    import hashlib as _hashlib
    return _hashlib.md5(text.encode("utf-8")).hexdigest()


def test_catalog_regzip_stable_md5_unstable_artifact_sha(tmp_path):
    """The same TSV text re-gzipped yields identical
    logical_content_md5 but different compressed-artifact SHA256
    (gzip mtime header) — the two hashes are inherently different
    claims and must never be conflated."""
    mod = _module()
    payload = ("chromosome\tposition\tref\talt\tsource_record_index\n"
               "1\t100\tA\tC\t1\n").encode("utf-8")
    with open(tmp_path / "a.tsv.gz", "wb") as f1, open(tmp_path / "b.tsv.gz", "wb") as f2:
        with gzip.GzipFile(fileobj=f1, mode="wb", mtime=0) as g1:
            g1.write(payload)
        with gzip.GzipFile(fileobj=f2, mode="wb", mtime=999999) as g2:
            g2.write(payload)
    a, b = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    assert mod.catalog_content_md5(a) == mod.catalog_content_md5(b)
    assert mod.catalog_content_md5(a) == _md5(payload.decode("utf-8"))
    assert mod.file_sha256(a) != mod.file_sha256(b), \
        "re-gzip changes the artifact hash: artifact claims are transport, not content"


def test_cli_round_trip_unresolved_paths(tmp_path):
    # tie path (clean catalog): the suspect enumeration -> UNRESOLVED
    fi_tie = _cli_input_files(tmp_path, tie=True)
    out_tie = tmp_path / "out_tie"
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(fi_tie["catalog"]),
         "--per-site-vcf", str(fi_tie["per_site_vcf"]),
         "--hdf5-values", str(fi_tie["h5_values"]),
         "--out", str(out_tie)])
    assert code == 3
    assert line["status"] == "UNRESOLVED"
    tie_spot = json.loads(
        (out_tie / "hdf5_vcf_concordance_spot_check.json").read_text())
    assert tie_spot["gate"]["status"] == "UNRESOLVED"
    # the overall orientation still resolves despite the suspects
    assert tie_spot["spot_check"]["overall"]["best_orientation"] == \
        "hdf5_0_is_ref"
    # the tie site (5,900) alone is the enumerated suspect
    assert len(tie_spot["spot_check"]["suspect_sites"]) == 1
    assert tie_spot["spot_check"]["suspect_sites"][0]["pos"] == 900
    assert any("suspect site" in reason
               for reason in tie_spot["gate"]["reasons"])
    # the machine line still carries the counts
    assert line["exit_code"] == 3 and line["n_sites"] == 3

    # the duplicated-row path (clean spot): duplicated rows alone
    # -> UNRESOLVED
    fi_dup = _cli_input_files(tmp_path, dup_catalog=True)
    out_dup = tmp_path / "out_dup"
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(fi_dup["catalog"]),
         "--per-site-vcf", str(fi_dup["per_site_vcf"]),
         "--hdf5-values", str(fi_dup["h5_values"]),
         "--out", str(out_dup)])
    assert code == 3
    dup_mapping = json.loads(
        (out_dup / "hdf5_vcf_position_mapping_summary.json").read_text())
    assert dup_mapping["position_mapping"][
        "duplicate_chrom_position_pairs"] == 1
    assert dup_mapping["position_mapping"][
        "duplicate_chrom_position_ref_alt_quadruples"] == 1
    assert dup_mapping["gate"]["status"] == "UNRESOLVED"
    assert any("duplicated (chr,pos) catalog rows" in reason
               for reason in dup_mapping["gate"]["reasons"])

    # the no-spot path: no per-site-vcf / hdf5-values -> UNRESOLVED
    fi_plain = _cli_input_files(tmp_path)
    out_nospot = tmp_path / "out_nospot"
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(fi_plain["catalog"]),
         "--out", str(out_nospot)])
    assert code == 3
    assert line["n_sites"] is None and line["n_suspects"] is None
    nospot = json.loads(
        (out_nospot / "hdf5_vcf_concordance_spot_check.json").read_text())
    assert nospot["spot_check"] is None
    assert any("spot check was not performed" in reason
               for reason in nospot["gate"]["reasons"])


def test_cli_orientation_pass_and_bad_input(tmp_path):
    # an orientation_pass that underperforms -> the mismatched sites
    fi = _cli_input_files(tmp_path)
    out = tmp_path / "out_passfail"
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(fi["catalog"]),
         "--per-site-vcf", str(fi["per_site_vcf"]),
         "--hdf5-values", str(fi["h5_values"]),
         "--orientation-pass", "hdf5_1_is_ref",
         "--out", str(out)])
    assert code == 3
    spot = json.loads(
        (out / "hdf5_vcf_concordance_spot_check.json").read_text())
    assert any("orientation_pass" in reason
               for site in spot["spot_check"]["suspect_sites"]
               for reason in site["suspect_reasons"])

    # a missing catalog: non-zero exit code, one-line JSON error, no
    # traceback
    fi2 = _cli_input_files(tmp_path)
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(tmp_path / "missing_catalog.gz"),
         "--out", str(tmp_path / "out_bad")])
    assert code == 4
    assert line["module"] == "at_pheno_audit_hdf5_vcf_mapping"
    assert "error" in line and "exit_code" in line
    # a region-span override that breaks the reassembling -> ValueError
    # -> exit 4
    code, line = _run_cli(
        tmp_path,
        ["--hdf5", "synthetic",
         "--catalog", str(fi2["catalog"]),
         "--region-spans",
         "[[0, 8], [8, 14], [14, 18], [18, 22], [23, 26]]",
         "--out", str(tmp_path / "out_spans")])
    assert code == 4
    assert "region span" in line["error"]


# ---------------------------------------------------------------- spans/attrs conversion (h5py loader)

def test_spans_chrs_from_attrs_plain_and_bytes():
    mod = _module()
    # h5py stores chrs as fixed-width bytes (dtype |S5 in the release
    # file); decoding is required -- str(b"1") would give "b'1'"
    # labels and corrupt the chrom join.
    spans, chrs = mod._spans_chrs_from_attrs(
        {"chr_regions": [[0, 2], [2, 4]], "chrs": [b"1", b"2"]})
    assert spans == [[0, 2], [2, 4]]
    assert chrs == ["1", "2"]


def test_spans_chrs_from_attrs_provided_lists_not_truth_tested():
    mod = _module()
    # The loader must key off `is not None`, never a truth test: an
    # explicitly provided-but-empty list of regions has to survive the
    # read, not be swapped for a default.
    spans, chrs = mod._spans_chrs_from_attrs({"chr_regions": [], "chrs": []})
    assert spans == [] and chrs == []
    # Absent keys: no spans, chrom labels fall back to the leading
    # NUCULAR labels for one region.
    spans2, chrs2 = mod._spans_chrs_from_attrs({})
    assert spans2 == [] and chrs2 == ["1"]


def test_spans_chrs_from_attrs_numeric_pairs_and_label_cast():
    mod = _module()
    # numpy-style scalars reach the loader as indexable numbers; the
    # spans come back as plain ints and non-bytes labels as str(text).
    spans, chrs = mod._spans_chrs_from_attrs(
        {"chr_regions": [[0, 2.0], [2, 4]], "chrs": [1, 2]})
    assert spans == [[0, 2], [2, 4]]
    assert chrs == ["1", "2"]
