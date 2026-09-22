"""Blockwise additive kernel and strict train-only preprocessing.

The kernel is centered dosage ZZ'/sum_j 2p_j(1-p_j), without HWE claims.
All callers must fit a fresh QC object for each inner or outer training set.
"""

from dataclasses import dataclass
import hashlib

import numpy as np


@dataclass
class QC:
    columns: np.ndarray
    means: np.ndarray
    frequencies: np.ndarray
    call_rates: np.ndarray
    mac: np.ndarray


@dataclass
class BinaryQC:
    """QC for orientation-unknown, complete 0/1 marker calls.

    The allele labelled 1 is deliberately not interpreted as ALT.  Flipping a
    whole marker from 0↔1 leaves the centered relationship kernel unchanged.
    """
    columns: np.ndarray
    means: np.ndarray
    state_one_frequencies: np.ndarray
    minor_state_count: np.ndarray


def blocks(columns, size):
    if size < 1:
        raise ValueError("block_size must be positive")
    for start in range(0, len(columns), size):
        yield start, columns[start:start + size]


def fit_qc(x, train, min_call_rate=0.95, min_maf=0.05, min_mac=0, block_size=4096):
    """Diploid ALT dosage 0/1/2; NaN means missing (never zero).

    For int8 input the formal storage convention is `-1` = missing,
    `0/1/2` = ALT dosage: missing entries are excluded from call
    counts, dosage sums and means (never treated as a dosage value).
    """
    if not 0 <= min_call_rate <= 1 or not 0 <= min_maf <= 0.5 or min_mac < 0:
        raise ValueError("Invalid QC thresholds")
    train = np.asarray(train, dtype=int)
    if len(train) < 2 or len(np.unique(train)) != len(train):
        raise ValueError("At least two unique training samples are required")
    kept, means, freqs, rates, macs = [], [], [], [], []
    for _, cols in blocks(np.arange(x.shape[1]), block_size):
        g = np.asarray(x[np.ix_(train, cols)])
        if g.dtype == np.int8:
            # Formal storage convention: -1 = missing (excluded from
            # every dosage statistic, never a dosage value),
            # 0/1/2 = ALT dosage.  int64 accumulation avoids int8
            # overflow in the per-column dose sums.
            g = g.astype(np.int64)
            if np.any((g < -1) | (g > 2)):
                raise ValueError(
                    "Expected int8 ALT-dosage values -1/0/1/2 "
                    "(-1 missing, 0/1/2 ALT dosage)")
            n = np.sum(g != -1, axis=0)
            total = np.sum(np.where(g == -1, 0, g), axis=0)
        else:
            g = g.astype(float)
            if np.any(~(np.isnan(g) | (g == 0) | (g == 1) | (g == 2))):
                raise ValueError("Expected diploid hard calls 0/1/2 or NaN")
            n = np.sum(~np.isnan(g), axis=0)
            total = np.nansum(g, axis=0)
        mean = np.divide(total, n, out=np.zeros(len(cols)), where=n > 0)
        p = mean / 2
        mac = np.minimum(total, 2 * n - total)
        rate = n / len(train)
        eligible = ((n > 0) & (rate >= min_call_rate) & (mac > 0)
                    & (np.minimum(p, 1-p) >= min_maf) & (mac >= min_mac))
        kept.extend(cols[eligible])
        means.extend(mean[eligible])
        freqs.extend(p[eligible])
        rates.extend(rate[eligible])
        macs.extend(mac[eligible])
    if not kept:
        raise ValueError("No eligible markers in this training partition")
    return QC(np.array(kept, dtype=int), np.array(means), np.array(freqs),
              np.array(rates), np.array(macs))


def fit_binary_qc(x, train, min_minor_state_frequency=0.05, min_minor_state_count=0, block_size=4096):
    """Fit QC for complete binary calls with no allele-orientation claim."""
    if not 0 <= min_minor_state_frequency <= 0.5 or min_minor_state_count < 0:
        raise ValueError("Invalid binary QC thresholds")
    train = np.asarray(train, dtype=int)
    if len(train) < 2 or len(np.unique(train)) != len(train):
        raise ValueError("At least two unique training samples are required")
    kept, means, freqs, macs = [], [], [], []
    for _, cols in blocks(np.arange(x.shape[1]), block_size):
        g = np.asarray(x[np.ix_(train, cols)], dtype=float)
        if np.any(~((g == 0) | (g == 1))):
            raise ValueError("Expected complete binary calls 0/1")
        mean = np.mean(g, axis=0)
        mac = np.minimum(np.sum(g, axis=0), len(train)-np.sum(g, axis=0))
        eligible = ((mac > 0) & (np.minimum(mean, 1-mean) >= min_minor_state_frequency)
                    & (mac >= min_minor_state_count))
        kept.extend(cols[eligible])
        means.extend(mean[eligible])
        freqs.extend(mean[eligible])
        macs.extend(mac[eligible])
    if not kept:
        raise ValueError("No eligible binary markers in this training partition")
    return BinaryQC(np.array(kept, dtype=int), np.array(means),
                    np.array(freqs), np.array(macs))


def marker_order(ids, salt):
    """Stable global hash ranking: nesting is exact, chromosome balance approximate."""
    if len(set(ids)) != len(ids):
        raise ValueError("Duplicate variant identifiers")
    digests = [hashlib.sha256(f"{salt}\0{s}".encode()).digest() for s in ids]
    return np.array(sorted(range(len(ids)), key=lambda i: (digests[i], ids[i])))


def select(qc, order, density):
    lookup = {int(c): i for i, c in enumerate(qc.columns)}
    eligible = [lookup[int(c)] for c in order if int(c) in lookup]
    if density != "all":
        if not isinstance(density, int) or density < 1:
            raise ValueError("density must be a positive integer or 'all'")
        eligible = eligible[:density]
    # Restore physical/manifest order for reproducible matrix accumulation.
    chosen = np.sort(eligible)
    return type(qc)(*(getattr(qc, field)[chosen] for field in vars(qc)))


def additive_kernel(x, train, test, qc, block_size=4096):
    train, test = np.asarray(train), np.asarray(test)
    k = np.zeros((len(train), len(train)))
    cross = np.zeros((len(test), len(train)))
    denominator = float(np.sum(2 * qc.frequencies * (1 - qc.frequencies)))
    if denominator <= 0:
        raise ValueError("Degenerate kernel")
    int8 = x.dtype == np.int8
    for start, cols in blocks(qc.columns, block_size):
        mean = qc.means[start:start + len(cols)]
        a_raw = np.asarray(x[np.ix_(train, cols)], dtype=float)
        b_raw = np.asarray(x[np.ix_(test, cols)], dtype=float)
        a = a_raw - mean
        b = b_raw - mean
        if int8:
            # Formal int8 convention: -1 missing entries contribute a
            # zero deviation to the kernel product, never (-1 - mean)
            # (that would count a missing call as a negative dosage,
            # which the frozen contract forbids).
            a = np.where(a_raw == -1.0, 0.0, a)
            b = np.where(b_raw == -1.0, 0.0, b)
        else:
            a = np.nan_to_num(a, nan=0.0)
            b = np.nan_to_num(b, nan=0.0)
        k += a @ a.T
        cross += b @ a.T
    return k / denominator, cross / denominator


def binary_additive_kernel(x, train, test, qc, block_size=4096):
    """Orientation-invariant additive kernel for a BinaryQC fitted on train."""
    train, test = np.asarray(train), np.asarray(test)
    k = np.zeros((len(train), len(train)))
    cross = np.zeros((len(test), len(train)))
    # This is a binary-state relationship kernel, not diploid GBLUP:
    # Var(state)=p(1-p), and p is orientation-unknown.
    denominator = float(np.sum(qc.state_one_frequencies * (1 - qc.state_one_frequencies)))
    if denominator <= 0:
        raise ValueError("Degenerate binary kernel")
    for start, cols in blocks(qc.columns, block_size):
        mean = qc.means[start:start + len(cols)]
        a = np.asarray(x[np.ix_(train, cols)], dtype=float) - mean
        b = np.asarray(x[np.ix_(test, cols)], dtype=float) - mean
        k += a @ a.T
        cross += b @ a.T
    return k / denominator, cross / denominator


def predict(k, cross, y, alpha):
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive and finite")
    # Columns were centered on this exact training set; intercept is y.mean().
    center = float(np.mean(y))
    return center + cross @ np.linalg.solve(k + alpha * np.eye(len(y)), y-center)


def scores(y, prediction, baseline):
    y, prediction, baseline = map(np.asarray, (y, prediction, baseline))
    if y.shape != prediction.shape or y.shape != baseline.shape or y.size == 0:
        raise ValueError("Aligned nonempty metric arrays required")
    if not all(np.isfinite(a).all() for a in (y, prediction, baseline)):
        raise ValueError("Metrics cannot contain missing or nonfinite values")
    sse = float(np.sum((y-prediction)**2))
    denom = float(np.sum((y-baseline)**2))
    pcc = None
    if y.size > 1 and np.std(y) > 0 and np.std(prediction) > 0:
        pcc = float(np.corrcoef(y, prediction)[0, 1])
    return {"n": int(y.size), "pcc": pcc, "rmse": float(np.sqrt(sse/y.size)),
            "mae": float(np.mean(np.abs(y-prediction))),
            "q2_train_mean": None if denom == 0 else 1-sse/denom}


def random_folds(ids, folds, seed):
    """IID accession CV pilot. Not an ancestry/kinship-controlled splitter."""
    if folds < 2 or folds > len(ids):
        raise ValueError("Invalid fold count")
    canonical = np.argsort(np.asarray(ids, dtype=str))
    rng = np.random.default_rng(seed)
    shuffled = rng.permutation(canonical)
    result = np.empty(len(ids), dtype=int)
    result[shuffled] = np.arange(len(ids)) % folds
    return result


def grouped_folds(groups):
    """Leave one declared genetic group out; labels are externally audited."""
    if any(not str(g).strip() for g in groups) or len(set(groups)) < 2:
        raise ValueError("At least two nonempty groups are required")
    mapping = {g: i for i, g in enumerate(sorted(set(groups)))}
    return np.array([mapping[g] for g in groups])
