"""Independent reference solves and leakage-sensitive numerical contracts."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import torch
from experiment import solve_grid
from preprocessing.panels import fit_qc, centered, ld_order, physical_order, group_folds


def test_cg_matches_primal_ridge_with_missing_dosages():
    rng = np.random.default_rng(42)
    x = rng.integers(0, 3, (47, 31)).astype(float)
    x[0, :5] = np.nan
    xv = rng.integers(0, 3, (13, 31)).astype(float)
    y = rng.normal(size=47)
    cfg = dict(call_rate_min=.9, maf_min=0, relative_residual_max=1e-8, cg_max_iter=1000)
    mean, keep = fit_qc(x, cfg)
    z = centered(x, mean, np.flatnonzero(keep))
    v = centered(xv, mean, np.flatnonzero(keep))
    scale = np.sum(z * z) / len(z)
    k, cross = z @ z.T / scale, v @ z.T / scale
    rhos = [1e-4, .01, 1., 100.]
    prediction, residual, _ = solve_grid(torch.tensor(k), torch.tensor(cross), y, rhos, cfg)
    for j, rho in enumerate(rhos):
        beta = np.linalg.solve(z.T @ z + len(z) * rho * scale * np.eye(z.shape[1]), z.T @ (y - y.mean()))
        np.testing.assert_allclose(prediction[:, j], y.mean() + v @ beta, atol=1e-7, rtol=1e-7)
    assert max(residual) < 1e-8


def test_ld_crosses_bin_boundaries_and_qc_is_training_only():
    rng = np.random.default_rng(73)
    x = rng.integers(0, 3, (100, 4)).astype(float)
    x[:, 1] = x[:, 0]  # linked markers straddle the physical-bin boundary
    x[:, 3] = 0
    cfg = dict(call_rate_min=.98, maf_min=0, ld_window_bp=100, ld_r2_threshold=.8, device='cpu')
    mean, keep = fit_qc(x, cfg)
    c = np.array(['c1'] * 4); p = np.array([990, 1010, 5000, 6000])
    physical = physical_order(c, p, 1, 1000)
    order, graph, _, n = ld_order(x, mean, keep, c, p, physical, cfg)
    assert graph[0] == {1} and graph[1] == {0}
    assert n == 2 and set(order[:2]) in [{0, 2}, {1, 2}]
    assert not keep[3]
    before = mean.copy()
    centered(np.full((3, 4), 99.), mean, np.flatnonzero(keep))
    np.testing.assert_array_equal(mean, before)


def test_group_folds_keep_members_together():
    groups = np.array(['a'] * 70 + ['b'] * 12 + [f'g{i}' for i in range(500)])
    a = group_folds(groups, 4, 101, 'full')
    np.testing.assert_array_equal(a, group_folds(groups, 4, 101, 'full'))
    for g in set(groups):
        assert len(set(a[groups == g])) == 1


def test_legacy_scaling_matches_training_only_standard_scaler():
    sys.path.append(str(Path(__file__).resolve().parents[2] / 'fprior/.deps'))
    from sklearn.preprocessing import StandardScaler
    from preprocessing.legacy import fit_transform, transform
    x = np.array([[1, 3, 2], [2, 3, 1], [3, 3, 3], [99, -20, 5]], dtype=np.float32)
    train = np.arange(3)
    mean, scale = fit_transform(x, train, True)
    reference = StandardScaler().fit(x[train])
    np.testing.assert_allclose(mean, reference.mean_)
    np.testing.assert_allclose(scale, reference.scale_)
    np.testing.assert_allclose(transform(x, np.arange(4), np.arange(3), mean, scale),
                               reference.transform(x), atol=1e-6)
