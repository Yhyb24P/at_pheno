"""Training-only marker processing and phenotype-free nested marker orders."""
from __future__ import annotations
import hashlib
import numpy as np


def group_folds(groups, n_folds, seed, context):
    unique, counts = np.unique(groups, return_counts=True)
    order = sorted(range(len(unique)), key=lambda i: (-counts[i], hashlib.sha256(
        f'{seed}\t{context}\t{unique[i]}'.encode()).hexdigest()))
    loads = np.zeros(n_folds, dtype=int)
    assignment = {}
    for i in order:
        fold = int(np.argmin(loads))
        assignment[unique[i]] = fold
        loads[fold] += counts[i]
    return np.array([assignment[g] for g in groups])


def fit_qc(x, config):
    called = np.isfinite(x)
    count = called.sum(axis=0)
    mean = np.divide(np.nansum(x, axis=0, dtype=np.float64), count,
                     out=np.zeros(x.shape[1]), where=count > 0)
    maf = np.minimum(mean / 2, 1 - mean / 2)
    variance = np.divide(np.nansum((x - mean) ** 2, axis=0), count,
                         out=np.zeros_like(mean), where=count > 0)
    keep = ((count / len(x) >= config['call_rate_min']) &
            (maf >= config['maf_min']) & (variance > 0))
    return mean, keep


def centered(x, mean, indices):
    z = np.asarray(x[:, indices], dtype=np.float64) - mean[indices]
    z[~np.isfinite(z)] = 0
    return z


def physical_order(chrom, pos, seed, bin_bp):
    bins = {}
    for i, (c, p) in enumerate(zip(chrom, pos)):
        bins.setdefault((c, int(p) // bin_bp), []).append(i)
    for key, values in bins.items():
        values.sort(key=lambda i: hashlib.sha256(f'{seed}\t{chrom[i]}\t{pos[i]}\t{i}'.encode()).digest())
    return interleave(bins)


def interleave(bins):
    keys = sorted(bins)
    return np.array([bins[k][j] for j in range(max(map(len, bins.values())))
                     for k in keys if j < len(bins[k])], dtype=int)


def ld_order(x_train, mean, keep, chrom, pos, physical, config):
    """Greedy local LD representatives, followed by suppressed markers.

    Edges span <= configured distance, including across physical-bin boundaries.
    A representative covers the most still-unrepresented r2>=threshold neighbours;
    ties use the fixed physical order. No validation genotypes enter this function.
    """
    import torch
    eligible = np.flatnonzero(keep)
    neighbours = {int(i): set() for i in eligible}
    rank = np.empty(len(physical), dtype=int)
    rank[physical] = np.arange(len(physical))
    summaries = []
    device = config['device']
    window = config['ld_window_bp']
    for c in np.unique(chrom):
        ids = eligible[chrom[eligible] == c]
        ids = ids[np.argsort(pos[ids], kind='stable')]
        pairs, total, high = 0, 0.0, 0
        # Blocks cap memory while including all local pairs, including bin edges.
        for start in range(0, len(ids), 128):
            left = ids[start:start + 128]
            stop = np.searchsorted(pos[ids], pos[left[-1]] + window, side='right')
            right = ids[start:stop]
            zl = centered(x_train, mean, left)
            zr = centered(x_train, mean, right)
            zl /= np.linalg.norm(zl, axis=0)
            zr /= np.linalg.norm(zr, axis=0)
            with torch.no_grad():
                corr = (torch.as_tensor(zl, device=device).T @ torch.as_tensor(zr, device=device)).cpu().numpy()
            r2 = np.clip(corr ** 2, 0, 1)
            mask = (pos[right][None, :] > pos[left][:, None]) & (pos[right][None, :] - pos[left][:, None] <= window)
            pairs += int(mask.sum())
            total += float(r2[mask].sum())
            ii, jj = np.where(mask & (r2 >= config['ld_r2_threshold']))
            high += len(ii)
            for a, b in zip(left[ii], right[jj]):
                neighbours[int(a)].add(int(b))
                neighbours[int(b)].add(int(a))
        summaries.append(dict(chromosome=str(c), local_pairs=pairs,
                              mean_r2=total / pairs if pairs else None, high_ld_pairs=high))
    remaining = set(map(int, eligible))
    representatives = set()
    # Heap updates only changed neighbourhoods; deterministic ties.
    import heapq
    heap = [(-len(neighbours[i]), int(rank[i]), i) for i in remaining]
    heapq.heapify(heap)
    while remaining:
        neg, tie, i = heapq.heappop(heap)
        if i not in remaining:
            continue
        degree = len(neighbours[i] & remaining)
        if -neg != degree:
            heapq.heappush(heap, (-degree, tie, i))
            continue
        representatives.add(i)
        removed = {i} | (neighbours[i] & remaining)
        remaining.difference_update(removed)
    # Stable physical interleaving retains broad chromosome coverage in each tier.
    order = np.array([i for i in physical if i in representatives] +
                     [i for i in physical if keep[i] and i not in representatives] +
                     [i for i in physical if not keep[i]], dtype=int)
    return order, neighbours, summaries, len(representatives)


def panel_stats(indices, chrom, pos, neighbours, bin_bp):
    chosen = set(map(int, indices))
    rows = []
    for c in np.unique(chrom):
        p = np.sort(pos[indices[chrom[indices] == c]])
        spacing = np.diff(p)
        catalog = pos[chrom == c]
        rows.append(dict(chromosome=str(c), markers=len(p),
                         occupied_bins=len(set((p // bin_bp).tolist())),
                         catalog_bins=len(set((catalog // bin_bp).tolist())),
                         span_bp=int(p[-1] - p[0]) if len(p) else 0,
                         median_spacing_bp=float(np.median(spacing)) if len(spacing) else None,
                         max_spacing_bp=int(spacing.max()) if len(spacing) else None,
                         high_ld_pairs=sum(len(neighbours.get(int(i), set()) & chosen)
                                           for i in indices if chrom[i] == c) // 2))
    return rows
