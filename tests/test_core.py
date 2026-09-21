import numpy as np
import pytest

from at_pheno.core import (additive_kernel, binary_additive_kernel, fit_binary_qc,
                           fit_qc, marker_order, predict, random_folds, scores, select)


def test_test_only_variant_is_excluded_and_test_values_cannot_change_qc():
    x = np.array([[0, 0, 0], [2, 0, np.nan], [0, 0, 2], [2, 2, 0]], dtype=float)
    qc = fit_qc(x, [0, 1, 2], min_call_rate=0.5)
    assert qc.columns.tolist() == [0, 2]
    x[3] = [np.nan, 0, 2]
    other = fit_qc(x, [0, 1, 2], min_call_rate=0.5)
    for name in vars(qc):
        np.testing.assert_array_equal(getattr(qc, name), getattr(other, name))


def test_block_kernel_matches_dense_primal_with_missing_and_intercept():
    rng = np.random.default_rng(2)
    x = rng.integers(0, 3, (12, 15)).astype(float)
    x[0, 0] = np.nan
    train, test = np.arange(9), np.arange(9, 12)
    qc = fit_qc(x, train, min_call_rate=0.5)
    k, cross = additive_kernel(x, train, test, qc, block_size=2)
    z = np.nan_to_num(x[:, qc.columns]-qc.means)
    scale = np.sum(2*qc.frequencies*(1-qc.frequencies))
    np.testing.assert_allclose(k, z[train]@z[train].T/scale)
    y, alpha = rng.normal(size=9)+5, 0.7
    weights = np.linalg.solve(z[train].T@z[train]+alpha*scale*np.eye(z.shape[1]), z[train].T@(y-y.mean()))
    np.testing.assert_allclose(predict(k, cross, y, alpha), y.mean()+z[test]@weights)


def test_density_is_nested_and_hash_order_stable_under_column_permutation():
    ids = [f"1:{j}:A:G" for j in range(1, 11)]
    order = marker_order(ids, "a")
    reverse = ids[::-1]
    assert [ids[i] for i in order] == [reverse[i] for i in marker_order(reverse, "a")]
    x = np.array([[0]*10, [2]*10], dtype=float)
    qc = fit_qc(x, [0, 1])
    assert set(select(qc, order, 3).columns) < set(select(qc, order, 7).columns)


def test_metric_degenerate_values_and_mean_centering_identity():
    y, p, c = np.array([1., 2, 4]), np.array([2., 1, 3]), np.array([1., 1, 2])
    assert scores(y, p, c)["rmse"] == scores(y-c, p-c, c-c)["rmse"]
    assert scores(y, np.ones(3), np.ones(3))["pcc"] is None
    with pytest.raises(ValueError):
        scores(y, np.array([np.nan]*3), c)


def test_random_split_is_stable_under_sample_reordering():
    ids = ["z", "a", "q", "b", "v", "n"]
    original = dict(zip(ids, random_folds(ids, 3, 9)))
    reverse = dict(zip(ids[::-1], random_folds(ids[::-1], 3, 9)))
    assert original == reverse


def test_inner_qc_excludes_variant_seen_only_in_inner_validation():
    x = np.array([[0, 0], [2, 0], [0, 0], [2, 2]], dtype=float)
    assert fit_qc(x, [0, 1, 2, 3]).columns.tolist() == [0, 1]
    assert fit_qc(x, [0, 1, 2]).columns.tolist() == [0]


def test_orientation_unknown_binary_kernel_is_allele_flip_invariant():
    rng = np.random.default_rng(11)
    x = rng.integers(0, 2, (10, 12), dtype=np.int8)
    train, test = np.arange(7), np.arange(7, 10)
    qc = fit_binary_qc(x, train, min_minor_state_frequency=0)
    k, cross = binary_additive_kernel(x, train, test, qc, block_size=3)
    flipped = x.copy()
    flipped[:, [1, 4, 8]] = 1-flipped[:, [1, 4, 8]]
    other = fit_binary_qc(flipped, train, min_minor_state_frequency=0)
    fk, fcross = binary_additive_kernel(flipped, train, test, other, block_size=2)
    np.testing.assert_allclose(k, fk)
    np.testing.assert_allclose(cross, fcross)
