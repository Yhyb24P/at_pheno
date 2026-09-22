"""Tiny VCF/GT -> int8 formal storage fixtures (P0-1, minimal).

The formal genotype contract of P0-1 (frozen 2026-09): int8
ALT-dosage matrix, shape n_samples x n_variants, value set
{-1, 0, 1, 2} (-1 missing, 0/1/2 ALT dosage); the tiny GT fixture
must round-trip 0/0, 0/1, 1/1, ./. -> 0, 1, 2, -1, and -1 must not
leak into any dosage statistic.
"""

import numpy as np
import pytest

from at_pheno.cli import load_dataset, write_table
from at_pheno.core import additive_kernel, fit_qc
from at_pheno.vcf_dosage import gt_cells_to_int8_dosage, gt_dosage_row

GT_COLUMN = ["0/0", "0/1", "1/1", "./."]


def test_gt_column_roundtrips_contract_int8():
    row = gt_dosage_row(GT_COLUMN)
    assert row.dtype == np.int8
    assert row.tolist() == [0, 1, 2, -1]


def test_gt_matrix_is_int8_samples_by_variants():
    cells = [
        ["0/0", "0/1", "1/1", "./."],
        ["1/1", "0/0", "./.", "0/1"],
    ]
    matrix = gt_cells_to_int8_dosage(cells)
    assert matrix.dtype == np.int8
    # (n_samples, n_variants) — samples along rows, variants along columns
    assert matrix.shape == (4, 2)
    assert matrix[0].tolist() == [0, 2]
    assert matrix[3].tolist() == [-1, 1]
    assert set(np.unique(matrix)) <= {-1, 0, 1, 2}


def test_phased_cells_share_dosage():
    assert gt_dosage_row(["0|0", "0|1", "1|0", "1|1"]).tolist() == [0, 1, 1, 2]


@pytest.mark.parametrize("bad", ["2/3", "./1", "0/0/1", "0", "", "0/."])
def test_unsupported_gt_cells_rejected(bad):
    with pytest.raises(ValueError, match="int8"):
        gt_dosage_row(["0/0", bad, "1/1", "./."])


def test_fit_qc_int8_matrix_excludes_minus_one_from_dosage_math():
    # (4 samples, 2 variants) int8; sample 2 is missing in variant 1
    x = np.array([[0, 0], [1, 1], [2, -1], [0, 2]], dtype=np.int8)
    qc = fit_qc(x, [0, 1, 2, 3], min_call_rate=0.0, min_maf=0.0, min_mac=0)
    assert qc.columns.tolist() == [0, 1]
    # If -1 were treated as a real dosage, variant 1 would have
    # mean(0+1-1+2)/4 = 0.5 and p = 0.25 instead:
    assert qc.means.tolist() == [0.75, 1.0]        # v1 mean over called [0,1,2]
    assert qc.frequencies.tolist() == [0.375, 0.5]
    assert qc.call_rates.tolist() == [1.0, 0.75]
    assert qc.mac.tolist() == [3, 3]               # v1 called: min(3, 2*3-3) = 3
    # the default 0.95 call-rate threshold drops the part-missing marker
    assert fit_qc(x, [0, 1, 2, 3]).columns.tolist() == [0]


def test_additive_kernel_int8_missing_deviation_is_zero():
    x = np.array([[0, 0], [1, 1], [2, -1], [0, 2]], dtype=np.int8)
    qc = fit_qc(x, [0, 1, 2, 3], min_call_rate=0.0, min_maf=0.0, min_mac=0)
    k, cross = additive_kernel(x, [0, 1, 2, 3], [0, 1, 2, 3], qc, block_size=4)
    # Hand-check: variant 0 values [0,1,2,0] mean 0.75;
    # variant 1 values [0,1,-1,2] called-mean 1.0; the missing entry
    # contributes zero deviation, NOT (-1 - mean).
    # (n_samples, n_variants) deviations; k is (train x train) and
    # cross is (test x train) — with train == test here both equal
    # dev @ dev.T / denom.
    dev = np.array([[-0.75, -1.0], [0.25, 0.0], [1.25, 0.0], [-0.75, 1.0]])
    denom = 2 * (0.375 * 0.625) + 2 * (0.5 * 0.5)
    exp = dev @ dev.T / denom
    assert np.allclose(k, exp)
    assert np.allclose(cross, exp)


def _make_int8_dataset(tmp_path, out_of_range=False):
    data = tmp_path/"data"
    data.mkdir()
    x = np.full((6, 8), -1, dtype=np.int8)
    x[0] = [0, 0, 1, 1, 2, -1, 0, 1]
    x[1] = [0, 1, 1, 2, 0, 0, -1, 2]
    x[2] = [2, -1, 0, 0, 1, 2, 1, 0]
    x[3] = [1, 2, 2, -1, 0, 1, 2, -1]
    x[4] = [0, 0, -1, 1, 2, 0, 0, 1]
    x[5] = [1, 1, 0, 2, -1, 2, -1, 0]
    if out_of_range:
        x[0, 0] = -2
    np.save(data/"genotypes.npy", x)
    write_table(data/"samples.tsv",
                [{"accession_id": f"s{i}", "genetic_group": ("G0" if i % 2 else "G1")}
                 for i in range(6)],
                ["accession_id", "genetic_group"])
    write_table(data/"variants.tsv",
                [{"chromosome": str(j // 4 + 1), "position": str((j % 4 + 1) * 1000),
                  "ref": "A", "alt": "G"} for j in range(8)],
                ["chromosome", "position", "ref", "alt"])
    write_table(data/"phenotypes.tsv",
                [{"accession_id": f"s{i}", "trait_id": "t1", "value": str(float(i))}
                 for i in range(6)],
                ["accession_id", "trait_id", "value"])
    return data


def test_load_dataset_accepts_int8_contract(tmp_path):
    data = _make_int8_dataset(tmp_path)
    x, rows, ids, variants, y, groups, audit = load_dataset(data, "t1")
    assert x.dtype == np.int8
    assert audit["genotyped_n"] == 6 and audit["marker_n"] == 8
    assert len(y) == 6


def test_formal_load_dataset_accepts_only_int8_contract(tmp_path):
    data = _make_int8_dataset(tmp_path)
    assert load_dataset(data, "t1", formal=True)[0].dtype == np.int8


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.int16, np.uint8])
def test_formal_load_dataset_rejects_non_int8_even_when_values_are_valid(tmp_path, dtype):
    data = _make_int8_dataset(tmp_path)
    x = np.load(data/"genotypes.npy")
    if np.issubdtype(dtype, np.unsignedinteger):
        x = np.where(x == -1, 0, x)
    np.save(data/"genotypes.npy", x.astype(dtype))
    with pytest.raises(ValueError, match=rf"dtype int8; got dtype {np.dtype(dtype)}"):
        load_dataset(data, "t1", formal=True)


def test_load_dataset_rejects_out_of_contract_int8_values(tmp_path):
    data = _make_int8_dataset(tmp_path, out_of_range=True)
    with pytest.raises(ValueError, match="int8 ALT-dosage"):
        load_dataset(data, "t1")


def test_formal_load_dataset_rejects_out_of_contract_int8_values(tmp_path):
    data = _make_int8_dataset(tmp_path, out_of_range=True)
    with pytest.raises(ValueError, match="int8 ALT-dosage"):
        load_dataset(data, "t1", formal=True)
