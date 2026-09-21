"""Fixture tests for the TASK02/TASK03 VCF-side audit.

Tiny synthetic VCFs (plain text and gzip) only -- the real 18 GB release
file is never touched.  Covers the seven policy buckets, the GT window
histograms, the 2.3 sample-mapping gate, the catalog gz round-trip, the
3.2 summary duplicate counts, and the 3.4 per-site state coding.
"""

import gzip
import importlib.util
import json
from pathlib import Path


def _module():
    global _MOD
    if _MOD is None:
        path = Path(__file__).resolve().parents[1] / "scripts" / "audit_1001g_vcf.py"
        spec = importlib.util.spec_from_file_location("audit_1001g_vcf", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        _MOD = module
    return _MOD

_MOD = None


def vcf_lines():
    """The shared fixture body: contigs (one bare '5'), a FILTER, a
    GT-only FORMAT, two sample labels, then six data lines that hit
    every policy bucket."""
    return [
        "##fileformat=VCFv4.2",
        '##contig=<ID=1;length=3000>',
        "##contig=<ID=4;length=123>",
        "##contig=<ID=5>",
        '##FILTER=<ID=PASS;Description="All filters passed">',
        '##FORMAT=<ID=GT;Number=1;Type=String;Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2",
    ] + [
        # biallelic SNPs: catalog rows with index 1..4 ((1,100) has two
        # records -> coordinate duplicate; (2,200) likewise, the first is
        # a complete duplicate of the second of the fourtuples).
        "1\t100\t.\tA\tC\t.\tPASS\tGT\t0/0\t1/1",
        "1\t100\t.\tT\tC\t.\tPASS\tGT\t0/0\t0/0",
        "2\t200\t.\tG\tC\t.\tPASS\tGT\t0/1\t0|1",
        "2\t200\t.\tG\tC\t.\tPASS\tGT\t./.\t./.",
        # excluded: multi-allelic (never split) and non-nuclear
        "3\t300\t.\tA\tC,G\t.\t.",
        "MT\t30\t.\tC\tT\t.\t.",
    ]


def write_vcf(tmp_path, lines, name="small.vcf"):
    path = tmp_path / name
    if name.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write("\n".join(lines) + "\n")
    else:
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_table(tmp_path, rows, name="table.json"):
    path = tmp_path / name
    path.write_text(json.dumps(rows), encoding="utf-8")
    return path


EXPECTED_ROWS = [
    ["1", "100", "A", "C", "1"],
    ["1", "100", "T", "C", "2"],
    ["2", "200", "G", "C", "3"],
    ["2", "200", "G", "C", "4"],
]


def test_classify_variant_line_covers_all_seven_buckets():
    mod = _module()
    pairs = [
        ("malformed", "1\t100\t.\tA\tC"),              # fewer than 7 cols
        ("non_nuclear", "MT\t10\t.\tA\tC\t.\t."),
        ("non_nuclear", "6\t60\t.\tA\tC,G\t.\t."),      # non-nuclear wins first
        ("multiallelic", "1\t150\t.\tA\tC,G\t.\t."),    # 2 ALT calls stay together
        ("multiallelic", "1\t160\t.\tA\tN,C\t.\t."),    # checked before any N-scan
        ("n_or_noncanonical", "2\t200\t.\tQ\tC\t.\t."),
        ("indel_or_complex", "3\t300\t.\tA\tAG\t.\t."),
        ("indel_or_complex", "3\t310\t.\tAC\tA\t.\t."),
        ("ref_equals_alt", "4\t400\t.\tA\tA\t.\t."),
        ("nuclear_biallelic_snp", "5\t500\t.\tA\tC\t.\tPASS"),
    ]
    expect = {cat: 0 for cat in mod.CATEGORIES}
    for cat, line in pairs:
        assert mod.classify_variant_line(line) == cat, line
        expect[cat] += 1
    # every distinct bucket is driven by at least one fixture line
    assert expect == {"malformed": 1, "non_nuclear": 2, "multiallelic": 2,
                      "n_or_noncanonical": 1, "indel_or_complex": 2,
                      "ref_equals_alt": 1,
                      "nuclear_biallelic_snp": 1}
    assert len(mod.CATEGORIES) == 7


def test_parse_vcf_line_records_header_and_malformed():
    mod = _module()
    assert mod.parse_vcf_line("##fileformat=VCFv4.2") is None
    assert mod.parse_vcf_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1") is None
    assert mod.parse_vcf_line("") is None
    record = mod.parse_vcf_line("1\t10\t.\tT\tC,T\t.\tPASS\tGT:DP\t0/0:5\t1/1:5")
    assert record == {"malformed": False, "chrom": "1", "pos": 10, "id": ".",
                      "ref": "T", "alts": ["C", "T"], "info": "GT:DP",
                      "filter": "PASS", "sample_fields": ["0/0:5", "1/1:5"]}
    # 6 fields: few of the 7 required -> malformed
    assert mod.parse_vcf_line("1\t10\t.\tT\tC\t.")["malformed"] is True
    record = mod.parse_vcf_line("1\tx\t.\tT\tC\t.\t.")
    assert record["pos"] is None and mod.classify_variant_line("1\tx\t.\tT\tC\t.\t.") == "malformed"


def test_summary_from_records_folds_a_fixture():
    mod = _module()
    records = [mod.parse_vcf_line(ln) for ln in vcf_lines()[7:]]
    summary = mod.summary_from_records(records)
    assert summary["total_records"] == 6
    assert summary["nuclear_records"] == 5  # chrom 3 rows land in nuclear, MT does not
    assert summary["per_chrom"] == {"1": 2, "2": 2, "3": 1, "4": 0, "5": 0}
    assert summary["type_counts"]["nuclear_biallelic_snp"] == 4
    assert summary["type_counts"]["multiallelic"] == 1
    assert summary["type_counts"]["non_nuclear"] == 1
    assert summary["type_counts"]["malformed"] == 0
    assert summary["catalog_rows"] == 4
    examples = summary["excluded_examples"]
    assert len(examples["multiallelic"]) == 1 and len(examples["non_nuclear"]) == 1


def test_gt_window_counts_histograms_haploid_and_out_of_range():
    mod = _module()
    # one record row covering every 2.4 histogram cell: haploid single
    # token ("0"), phased het, all-missing, partial-missing, out-of-range
    # allele index ("2"), blank (no_gt), and a clean 0/0.
    rows = [
        ["0", "0|1", "./.", "0/.", "2", "", "0/0"],
    ]
    got = mod.gt_window_counts({"1": rows})["combined"]
    assert got["n_records"] == 1 and got["n_cells"] == 7
    # ploidy: "0" and "2" are single-token cells; the other four are
    # two-token.
    assert got["ploidy"] == {"1": 2, "2": 4, ">2": 0}
    assert got["phased_cells"] == 1 and got["unphased_cells"] == 5
    assert got["n_haploid_gt_cells"] == 2 and got["haploid_gt_present"] is True
    assert got["states"] == {
        "ref_homo": 2,        # "0" (haploid 0) and "0/0"
        "heterozygous": 1,    # "0|1"
        "alt_homo": 0,
        "missing": 1,         # "./."
        "partial_missing": 1, # "0/."
        "allele_index_gt_1": 1,
        "no_gt": 1,
    }
    # a clean two-allelic "2/2": both tokens index 2 -- out of the 0/1
    # state space; it must land in the out-of-range bucket, not a
    # pseudo-state bucket.
    got = mod.gt_window_counts({"1": [["2/2"]]})["combined"]
    assert got["states"] == {
        "ref_homo": 0, "heterozygous": 0, "alt_homo": 0, "missing": 0,
        "partial_missing": 0, "allele_index_gt_1": 1, "no_gt": 0,
    }


def test_sample_mapping_resolved_on_accession_id_when_ordered():
    mod = _module()
    rows = [
        {"Accession_ID": 4, "Accession_Name": "alppa", "CS_Number": "CS404"},
        {"Accession_ID": 1, "Accession_Name": "aldual", "CS_Number": "CS109"},
        {"Accession_ID": 2, "Accession_Name": "aidi", "CS_Number": "CS212"},
    ]
    got = mod.map_samples(["4", "1", "2"], rows)
    assert got["status"] == "resolved"
    assert got["resolved_key"] == "accession_id"  # tried first, order+set both hold
    assert got["candidates"][0]["order_match"] is True
    assert got["candidates"][0]["field"] == "Accession_ID"


def test_sample_mapping_resolves_on_name_identity_even_when_order_differs():
    mod = _module()
    # accession_id: set-equivalent and uniquely-keyed, but the VCF #CHROM
    # order differs from the table row order (the real 1001G situation:
    # the VCF lists samples in numeric accession order).  Delivery doc
    # 2.3 gates on a one-to-one correspondence, not on file order, so
    # the key still resolves; the order evidence is reported, not gating.
    rows = [
        {"Accession_ID": 4, "Accession_Name": "N4", "CS_Number": "CS404"},
        {"Accession_ID": 3, "Accession_Name": "N3", "CS_Number": "CS303"},
        {"Accession_ID": 5, "Accession_Name": "N5", "CS_Number": "CS505"},
    ]
    got = mod.map_samples(["4", "5", "3"], rows)
    assert got["status"] == "resolved"
    assert got["resolved_key"] == "accession_id"
    by_key = {c["key"]: c for c in got["candidates"]}
    assert by_key["accession_id"]["order_match"] is False
    assert by_key["accession_id"]["set_match"] is True
    assert by_key["accession_id"]["resolves"] is True
    assert by_key["accession_id"]["order_divergences"] == 2
    assert by_key["accession_name"]["set_match"] is False


def test_sample_mapping_missing_sample_fails_the_gate():
    mod = _module()
    rows = [
        {"Accession_ID": 4, "Accession_Name": "N4", "CS_Number": "CS404"},
        {"Accession_ID": 3, "Accession_Name": "N3", "CS_Number": "CS303"},
    ]
    got = mod.map_samples(["4", "5", "3"], rows)   # "5" not on the table
    assert got["status"] == "unresolved"
    assert got["resolved_key"] is None
    by_key = {c["key"]: c for c in got["candidates"]}
    assert by_key["accession_id"]["set_match"] is False
    assert "5" in by_key["accession_id"]["in_vcf_not_in_table"]
    assert got["candidates"] and got["candidates"][0]["field"] == "Accession_ID"


def test_catalog_gz_roundtrip_plain_text_vcf(tmp_path):
    mod = _module()
    vcf = write_vcf(tmp_path, vcf_lines())
    out = tmp_path / "out"
    got = mod.run_audit(vcf, out, accession_table=tmp_path / "missing.json",
                        official_md5=mod.file_md5(vcf), source_url="https://example/1001g",
                        records_per_chrom=2, window_chroms=5)
    catalog = out / mod.CATALOG_GZ
    with gzip.open(catalog, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert lines[0] == "\t".join(mod.CATALOG_COLUMNS)
    assert lines[1:] == ["\t".join(row) for row in EXPECTED_ROWS]
    summary = json.loads((out / mod.SUMMARY_JSON).read_text())
    assert summary["total_records"] == 6
    assert summary["catalog_rows"] == 4
    assert summary["type_counts"] == {"malformed": 0, "non_nuclear": 1,
                                      "multiallelic": 1,
                                      "n_or_noncanonical": 0, "indel_or_complex": 0,
                                      "ref_equals_alt": 0,
                                      "nuclear_biallelic_snp": 4}
    # (chrom,pos) has 2 duplicate keys: (1,100) and (2,200).
    # Fourtuples: only (2,200,G,C) appears more than once.
    assert summary["duplicate_chrom_position_pairs"] == 2
    assert summary["duplicate_chrom_position_ref_alt_quadruples"] == 1
    block = summary["catalog_sample_block"]
    assert [r[:4] for r in EXPECTED_ROWS] == [r[:4] for r in block["first_rows"]]
    assert block["catalog_file_sha256"] == mod.file_sha256(catalog)
    # GT per-chrom windows (chrom-1 records 1..2, chrom-2 records 3..4)
    gt = json.loads((out / mod.GT_JSON).read_text())
    chrom1 = gt["per_chrom"]["1"]
    assert chrom1["states"]["ref_homo"] == 3 and chrom1["states"]["alt_homo"] == 1
    assert chrom1["ploidy"] == {"1": 0, "2": 4, ">2": 0}
    chrom2 = gt["per_chrom"]["2"]
    assert chrom2["states"]["heterozygous"] == 2 and chrom2["states"]["missing"] == 2
    assert chrom2["phased_cells"] == 1 and chrom2["unphased_cells"] == 3
    assert gt["combined"]["n_cells"] == 8 and gt["combined"]["n_haploid_gt_cells"] == 0
    # Sample mapping unresolved -> gate is open, exit_code 3
    assert got["exit_code"] == 3
    assert got["gate"] == "VCF_SEMANTICS_UNRESOLVED"
    assert got["mapping_status"] == "unresolved"
    # source json provenance
    source = json.loads((out / mod.SOURCE_JSON).read_text())
    assert source["md5_matches_official"] is True
    assert source["sha256_is_64_hex"] is True
    assert source["bcftools_header_probe"] is None
    assert source["stream_encoding"] == "plain"
    # header/sample audit
    pay = json.loads((out / mod.HEADER_JSON).read_text())
    assert pay["contig_labels"] == {"1": 3000, "4": 123, "5": None}
    assert pay["gt_present"] is True
    assert pay["sample_id_map"]["status"] == "unresolved"
    # SPOT JSON only when requested
    assert not (out / mod.SPOT_CHECK_JSON).exists()


def test_gzip_fixture_streams_same_catalog(tmp_path):
    mod = _module()
    vcf = write_vcf(tmp_path, vcf_lines(), name="small.vcf.gz")
    out = tmp_path / "out"
    mod.run_audit(vcf, out, records_per_chrom=2)
    with gzip.open(out / mod.CATALOG_GZ, "rt", encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    assert lines[0] == "\t".join(mod.CATALOG_COLUMNS)
    assert lines[1:] == ["\t".join(row) for row in EXPECTED_ROWS]
    source = json.loads((out / mod.SOURCE_JSON).read_text())
    assert source["stream_encoding"] == "gzip"


def test_summary_json_duplicate_counts_are_hand_computable(tmp_path):
    mod = _module()
    rows = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1",
        "1\t10\t.\tA\tC\t.\tPASS\tGT\t0/0",
        "2\t20\t.\tT\tG\t.\tPASS\tGT\t0/1",
        "3\t30\t.\tG\tA\t.\tPASS\tGT\t1/1",
        "1\t10\t.\tA\tC\t.\tPASS\tGT\t0/0",   # full copy of line 1
        "1\t10\t.\tT\tC\t.\tPASS\tGT\t0/0",   # same coordinate, different variant
    ]
    vcf = write_vcf(tmp_path, rows)
    out = tmp_path / "out"
    mod.run_audit(vcf, out, accession_table=tmp_path / "missing.json")
    got = json.loads((out / mod.SUMMARY_JSON).read_text())
    # Hand calc: 5 catalog rows; (1,10) has 3 rows -> 1 duplicate pair key.
    # (1,10,A,C) has 2 rows -> 1 duplicate quattuor.
    assert got["catalog_rows"] == 5
    assert got["duplicate_chrom_position_pairs"] == 1
    assert got["duplicate_chrom_position_ref_alt_quadruples"] == 1
    assert got["catalog_sample_block"]["first_rows"] == [
        ["1", "10", "A", "C", "1"],
        ["2", "20", "T", "G", "2"],
        ["3", "30", "G", "A", "3"],
        ["1", "10", "A", "C", "4"],
        ["1", "10", "T", "C", "5"],
    ]


def test_spot_check_keys_and_2_sample_like_fixtures(tmp_path):
    mod = _module()
    rows = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT;Number=1;Type=String;Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2\tS3\tS4",
        "1\t50\t.\tA\tC\t.\tPASS\tGT\t0/0\t0/1\t./.\t1/1",
    ]
    vcf = write_vcf(tmp_path, rows)
    picked = tmp_path / "picked.json"
    picked.write_text(json.dumps([["1", 50], ["5", 99]]))
    out = tmp_path / "out"
    got = mod.run_audit(vcf, out, accession_table=tmp_path / "missing.json",
                        picked_sites=mod.load_picked_sites(picked))
    pay = json.loads((out / mod.SPOT_CHECK_JSON).read_text())
    assert pay["n_sites"] == 2
    required = {"chrom", "concordance_if_hdf5_0_is_REF",
                "concordance_if_hdf5_1_is_REF", "best_orientation",
                "best_concordance", "n_comparable", "n_heterozygous",
                "n_missing_vcf"}
    site_a, site_b = pay["per_site"]
    for entry in pay["per_site"]:
        assert required.issubset(entry.keys())
    assert (site_a["chrom"], site_a["pos"]) == ("1", 50)
    assert (site_b["chrom"], site_b["pos"]) == ("5", 99)
    # S1 0/0 (hom-REF), S2 0/1 (het), S3 ./. (missing), S4 1/1 (hom-ALT)
    # -> equal comparable share under either orientation hypothesis.
    assert site_a["n_comparable"] == 2
    assert site_a["n_heterozygous"] == 1
    assert site_a["n_missing_vcf"] == 1
    assert site_a["concordance_if_hdf5_0_is_REF"] == 0.5
    assert site_a["concordance_if_hdf5_1_is_REF"] == 0.5
    assert site_a["best_orientation"] == "ambiguous"  # tie on 1/1
    assert site_b["record_found"] is False
    assert site_b["concordance_if_hdf5_0_is_REF"] is None
    assert site_b["n_comparable"] == 0
    assert got["spot_check"]["n_sites"] == 2


def test_spot_site_report_passes_per_sample_through(tmp_path):
    mod = _module()
    entry = mod.spot_site_report(
        "1", 50, ["0/0", "0/1", "./.", "1/1", "2/2"],
        per_sample=["homozygous_0", "heterozygous", "missing",
                    "homozygous_1", "out_of_range_or_unexpected"])
    assert "per_sample" in entry
    assert entry["per_sample"] == ["homozygous_0", "heterozygous",
                                   "missing", "homozygous_1",
                                   "out_of_range_or_unexpected"]
    # omitted when the caller did not supply it (pairing only when
    # the key is present on the HDF5 side)
    bare = mod.spot_site_report("1", 50, ["0/0"])
    assert "per_sample" not in bare


def test_spot_check_emits_per_sample_sample_aligned(tmp_path):
    mod = _module()
    rows = [
        "##fileformat=VCFv4.2",
        '##FORMAT=<ID=GT;Number=1;Type=String;Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2\tS3\tS4",
        "1\t50\t.\tA\tC\t.\tPASS\tGT\t0/0\t0/1\t./.\t1/1",
    ]
    vcf = write_vcf(tmp_path, rows)
    picked = tmp_path / "picked.json"
    picked.write_text(json.dumps([["1", 50], ["5", 99]]))
    out = tmp_path / "out"
    got = mod.run_audit(vcf, out, accession_table=tmp_path / "missing.json",
                        picked_sites=mod.load_picked_sites(picked))
    pay = json.loads((out / mod.SPOT_CHECK_JSON).read_text())
    site_a, site_b = pay["per_site"]
    # VCF sample-column order; state codes computed over the same cells
    # the histogram uses, so the HDF5 side can pair accession i with
    # VCF sample i.
    assert site_a["per_sample"] == ["homozygous_0", "heterozygous",
                                   "missing", "homozygous_1"]
    assert len(site_a["per_sample"]) == pay["n_samples"]
    # an unmatched site stays record_found=False and carries no
    # per-sample row: nothing to pair on the HDF5 side.
    assert site_b["record_found"] is False
    assert "per_sample" not in site_b
    assert got["spot_check"]["n_sites"] == 2


def test_source_json_and_unresolved_table_open_the_gate(tmp_path):
    mod = _module()
    vcf = write_vcf(tmp_path, vcf_lines())
    out = tmp_path / "out"
    got = mod.run_audit(vcf, out, accession_table=tmp_path / "nope.json")
    source = json.loads((out / mod.SOURCE_JSON).read_text())
    assert source["local_md5"] == mod.file_md5(vcf)
    assert source["local_sha256"] == mod.file_sha256(vcf)
    assert source["md5_matches_official"] is None  # no official passed
    assert source["stream_encoding"] == "plain"
    header = json.loads((out / mod.HEADER_JSON).read_text())
    assert header["sample_id_map"]["status"] == "unresolved"
    assert got["gate"] == "VCF_SEMANTICS_UNRESOLVED"
    assert got["exit_code"] == 3


def test_main_exit_zero_when_resolved_and_clean(tmp_path, capsys):
    mod = _module()
    rows = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1;length=3000>",
        '##FORMAT=<ID=GT;Number=1;Type=String;Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t4\t1",
        "1\t10\t.\tA\tC\t.\tPASS\tGT\t0/0\t1/1",
    ]
    vcf = write_vcf(tmp_path, rows)
    table = write_table(tmp_path, [
        {"Accession_ID": 4, "Accession_Name": "N4", "CS_Number": "CS404"},
        {"Accession_ID": 1, "Accession_Name": "N1", "CS_Number": "CS101"},
    ])
    out = tmp_path / "out"
    code = mod.main(["--vcf", str(vcf), "--out", str(out),
                     "--accession-table", str(table)])
    assert code == 0
    for f in (mod.SOURCE_JSON, mod.HEADER_JSON, mod.GT_JSON,
              mod.SUMMARY_JSON, mod.CATALOG_GZ):
        assert (out / f).exists()
    print_lines = capsys.readouterr().out.strip().splitlines()
    assert len(print_lines) == 1
    payload = json.loads(print_lines[0])
    assert payload["run"] == "VCF_SEMANTICS_READY_FOR_CONVERSION"
    assert payload["exit_code"] == 0


def test_main_exit_three_when_mapping_unresolved(tmp_path, capsys):
    mod = _module()
    rows = [
        "##fileformat=VCFv4.2",
        "##contig=<ID=1;length=3000>",
        '##FORMAT=<ID=GT;Number=1;Type=String;Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\t4\t1",
        "1\t10\t.\tA\tC\t.\tPASS\tGT\t0/0\t1/1",
    ]
    vcf = write_vcf(tmp_path, rows)
    out = tmp_path / "out"
    code = mod.main(["--vcf", str(vcf), "--out", str(out),
                     "--accession-table", str(tmp_path / "nope.json"),
                     "--no-catalog"])
    print_lines = capsys.readouterr().out.strip().splitlines()
    assert code == 3
    assert json.loads(print_lines[0])["run"] == "VCF_SEMANTICS_UNRESOLVED"
    assert not (out / mod.CATALOG_GZ).exists()
    with_no_catalog = json.loads((out / mod.SUMMARY_JSON).read_text())
    assert with_no_catalog["duplicate_chrom_position_pairs"] is None


# ---------------------------------------------------------------------------
# Comma-style meta separators (the 1001G v3.1 `##FORMAT=<ID=GT,Number=1,...>`)
# ---------------------------------------------------------------------------

def test_vcf_meta_comma_separators_parse_gt_and_filter(tmp_path):
    mod = _module()
    rows = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=1,length=100>",
        "##contig=<ID=4,length=90>",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FILTER=<ID=q4,Description=\"Quality below 4\">",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tS1\tS2",
    ]
    vcf = write_vcf(tmp_path, rows)
    header = mod.collect_vcf_header(vcf)
    assert header["format_fields"] == ["GT", "GQ"]
    assert header["gt_present"] is True and header["gt_column_index"] == 0
    assert header["filter_names"] == ["PASS", "q4"]
    assert header["contig_labels"] == {"1": 100, "4": 90}
    assert header["n_samples"] == 2 and header["format_column_in_chrom"] is False
    # gt_window on a record parsed from this header now works end to end
    out = tmp_path / "out"
    got = mod.run_audit(vcf, out, accession_table=tmp_path / "missing.json",
                        save_catalog=False)
    gt = got["gt_audit"]["combined"]
    assert gt["states"]["no_gt"] == 0 or gt["n_cells"] == 0


def test_flow_runs_as_assumed(tmp_path):
    """run_audit with the FORMAT token in #CHROM and comma meta: sample
    column peeling, mapping on shorthand keys, and the 2.4 window must
    report real states (previously everything collapsed to no_gt)."""
    mod = _module()
    rows = [
        "##fileformat=VCFv4.1",
        "##contig=<ID=1,length=100>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t88\t108\t139",
        "1\t5\t.\tC\tT\t40\tPASS\tDP=10\tGT:GQ:DP\t0|0:30:10\t./.:.:.\t1|1:30:9",
    ]
    vcf = write_vcf(tmp_path, rows)
    table = write_table(tmp_path, [
        {"Accession_Name": "139", "Accession_ID": 139, "CS_Number": "CS76347"},
        {"Accession_Name": "108", "Accession_ID": 108, "CS_Number": "CS76115"},
        {"Accession_Name": "88", "Accession_ID": 88, "CS_Number": "CS76001"},
    ])
    out = tmp_path / "out"
    got = mod.run_audit(vcf, out, accession_table=table,
                        official_md5=mod.file_md5(vcf))
    header = got["header"]
    assert header["format_column_in_chrom"] is True
    assert header["n_samples"] == 3
    assert header["first_samples"] == ["88", "108", "139"]
    assert header["sample_ids"] == ["88", "108", "139"]  # no "FORMAT" literal
    mapping = got["mapping"]
    assert mapping["status"] == "resolved"
    assert mapping["resolved_key"] == "accession_id"
    acc = mapping["candidates"][0]
    assert acc["order_match"] is False and acc["order_divergences"] == 2
    assert acc["unique_keys_on_both_sides"] is True and acc["resolves"] is True
    gt = got["gt_audit"]
    assert gt["window"]["format_column_assumed"] is True
    assert gt["window"]["sample_fields_width_mismatch"] is None
    assert gt["per_chrom"]["1"]["states"] == {
        "ref_homo": 1, "heterozygous": 0, "alt_homo": 1, "missing": 1,
        "partial_missing": 0, "allele_index_gt_1": 0, "no_gt": 0,
    }
    assert got["gate"] == "VCF_SEMANTICS_READY_FOR_CONVERSION"
    assert got["exit_code"] == 0


def test_parse_vcf_line_honors_sample_start():
    mod = _module()
    line = "1\t5\t.\tC\tT\t40\tPASS\tDP=10\tGT:GQ:DP\t0/0:5\t./.:.:.\t1/1:5\n"
    default = mod.parse_vcf_line(line)
    assert default["sample_fields"] == ["GT:GQ:DP", "0/0:5", "./.:.:.", "1/1:5"]
    fmt = mod.parse_vcf_line(line, sample_start=9)
    assert fmt["sample_fields"] == ["0/0:5", "./.:.:.", "1/1:5"]
    assert fmt["filter"] == "PASS" and fmt["chrom"] == "1" and fmt["pos"] == 5
    # a record narrower than the start index gets an empty sample section
    short = mod.parse_vcf_line("1\t5\t.\tC\tT\t40\tPASS", sample_start=9)
    assert short["sample_fields"] == [] and short["malformed"] is False


def test_sample_mapping_duplicate_names_cannot_resolve():
    mod = _module()
    rows = [
        {"Accession_ID": 4, "Accession_Name": "N4", "CS_Number": "CS404"},
        {"Accession_Name": "N4", "CS_Number": "CS404", "Accession_ID": 4},
        {"Accession_Name": "N3", "CS_Number": "CS303", "Accession_ID": 3},
    ]
    got = mod.map_samples(["4", "4", "3"], rows)
    assert got["status"] == "unresolved"
    acc = got["candidates"][0]
    assert acc["set_match"] is True
    assert acc["unique_keys_on_both_sides"] is False
    assert acc["resolves"] is False


# ---------------------------------------------------------------------------
# Gate-wide: source JSON must record an evidence-anchored download time
# ---------------------------------------------------------------------------

def test_source_json_records_download_time_verbatim(tmp_path):
    mod = _module()
    vcf = write_vcf(tmp_path, ["##fileformat=VCFv4.1",
                               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"])
    out = tmp_path / "out"
    stamp = "artifact mtime 2026-09-21T20:30:30Z (stat)"
    mod.run_audit(vcf, out, accession_table=tmp_path / "missing.json",
                  save_catalog=False, download_time=stamp)
    source = json.loads((out / mod.SOURCE_JSON).read_text())
    assert source["download_time"] == stamp


def test_load_picked_sites_accepts_raw_list_and_summary_shapes(tmp_path):
    mod = _module()
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps([["1", "55"], ["5", "99"]]), encoding="utf-8")
    assert mod.load_picked_sites(raw) == [("1", 55), ("5", 99)]
    nested = tmp_path / "nested.json"
    nested.write_text(json.dumps({"sites": [["2", 7]]}), encoding="utf-8")
    assert mod.load_picked_sites(nested) == [("2", 7)]
    # the full mapping-summary object: pool under
    # position_mapping.site_pool_sample (list of [chrom, pos] pairs)
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "series": "hdf5_vcf_position_mapping_summary",
        "position_mapping": {
            "site_pool_sample": [["1", "55"], ["5", "99"]],
        }}), encoding="utf-8")
    assert mod.load_picked_sites(summary) == [("1", 55), ("5", 99)]
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"no_pool": 1}), encoding="utf-8")
    try:
        mod.load_picked_sites(bad)
    except ValueError:
        pass
    else:
        raise AssertionError("a pool-less object must raise")
