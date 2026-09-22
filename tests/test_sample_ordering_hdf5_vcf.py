"""Committed-file check: the SAMPLE-ORDERING flag between the two
1001G sides.

The audit stage of 2026-09 recorded `hdf5_0_is_ref: label-only` and
left SAMPLE-ORDERING (the manual cross-check of which accession lines
up with which VCF sample column) as an open item at the gate.  This
test closes it with a deterministic, committable equivalence check
against the two committed manifest JSONs (no 12 GB re-read):

* `1001g_hdf5_2026-09-22.json` -> `accessions`: the 1135 sample labels
  of the HDF5 genotype matrix, in HDF5 sample order;
* `1001g_v31_vcf_header_sample_audit.json` -> `sample_ids`: the 1135
  sample names from the VCF header, in VCF sample-column order.

Ordered equality is the claim (not set equality): position `i` refers
to the same accession on both sides, 1135/1135, 0 divergences.
"""

import json
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]

N_EXPECTED = 1135


def _manifest_key(base_name, key):
    return json.loads(
        (_REPO / "data/manifests" / base_name).read_text(encoding="utf-8"))[key]


def test_sample_ordering_hdf5_vcf_is_one_to_one_and_ordered():
    accessions = _manifest_key("1001g_hdf5_2026-09-22.json", "accessions")
    sample_ids = _manifest_key(
        "1001g_v31_vcf_header_sample_audit.json", "sample_ids")
    # Defensive: keep both sides comparable (the committed files use
    # strings, but the check must not silently break if a future audit
    # emits ints or mixed types).
    accessions = [str(x) for x in accessions]
    sample_ids = [str(x) for x in sample_ids]

    assert len(accessions) == N_EXPECTED, f"HDF5 side carries {len(accessions)}, expected {N_EXPECTED}"
    assert len(sample_ids) == N_EXPECTED, f"VCF side carries {len(sample_ids)}, expected {N_EXPECTED}"
    assert len(set(accessions)) == N_EXPECTED, "HDF5 accessions must be 1:1 unique"
    assert len(set(sample_ids)) == N_EXPECTED, "VCF sample ids must be 1:1 unique"
    assert accessions == sample_ids, (
        "SAMPLE-ORDERING check: HDF5 accessions and VCF sample_ids must be "
        "identical in both set and order (0 divergences)")
