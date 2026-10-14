"""Independent soybean density experiment; fixed grouped outer evaluation."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import sys
import time
import numpy as np
import torch
from preprocessing.panels import (group_folds, fit_qc, centered, physical_order,
                                  ld_order, panel_stats)


def read_tsv(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False))
    tmp.replace(path)


def write_tsv(path, rows):
    with Path(path).open('w') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def solve_grid(k, cross, y, rhos, config):
    """Batched preconditioned CG, checked against the original ridge equation."""
    yt = torch.as_tensor(y, dtype=torch.float64, device=k.device)
    target = yt - yt.mean()
    shifts = torch.as_tensor(rhos, dtype=torch.float64, device=k.device) * len(y)
    b = target[:, None].expand(-1, len(rhos))
    a = torch.zeros_like(b)
    residual = b.clone()
    diagonal = torch.diag(k)[:, None] + shifts
    z = residual / diagonal
    direction = z.clone()
    rz = (residual * z).sum(0)
    norm = torch.linalg.vector_norm(target)
    tolerance = config['relative_residual_max'] * .1
    active = torch.ones(len(rhos), dtype=torch.bool, device=k.device)
    iterations = torch.zeros(len(rhos), dtype=torch.int64, device=k.device)
    for step in range(config['cg_max_iter']):
        product = k @ direction + direction * shifts
        denominator = (direction * product).sum(0)
        if bool((denominator[active] <= 0).any()):
            raise RuntimeError('CG non-positive curvature')
        ratio = torch.where(active, rz / torch.where(active, denominator, 1), 0)
        a += direction * ratio
        residual -= product * ratio
        converged = torch.linalg.vector_norm(residual, dim=0) <= tolerance * norm
        iterations[active] += 1
        active &= ~converged
        if not bool(active.any()):
            break
        new_z = residual / diagonal
        new_rz = (residual * new_z).sum(0)
        beta = torch.where(active, new_rz / torch.where(active, rz, 1), 0)
        direction = new_z + direction * beta
        direction[:, ~active] = 0
        rz = new_rz
    actual = torch.linalg.vector_norm(k @ a + a * shifts - b, dim=0) / norm
    if not bool(torch.isfinite(actual).all()) or float(actual.max()) > config['relative_residual_max']:
        raise RuntimeError(f'CG residual failed: {actual.tolist()}')
    pred = yt.mean() + cross @ a
    if not bool(torch.isfinite(pred).all()):
        raise RuntimeError('Nonfinite predictions')
    return pred.cpu().numpy(), actual.cpu().tolist(), iterations.cpu().tolist()


def load_inputs(config, root):
    data = (root / config['data']).resolve()
    manifest = json.loads((data / 'partition_manifest.json').read_text())
    # Validate permitted non-label artifacts against existing source identity.
    source = json.loads((data / 'source_manifest.json').read_text())
    records = {Path(a['path']).name: a for a in source['artifacts']}
    for name in ['samples.tsv', 'variants.tsv', 'splits.tsv']:
        if hashlib.sha256((data / name).read_bytes()).hexdigest() != records[name]['sha256']:
            raise ValueError(f'Input identity mismatch: {name}')
    rows = read_tsv(data / 'splits.tsv')
    dev = [r for r in rows if r['partition'] == 'development']
    sealed = {r['accession_id'] for r in rows if r['partition'] == 'sealed'}
    ids = [r['accession_id'] for r in dev]
    if len(ids) != 14020 or len(set(ids)) != 14020 or set(ids) & sealed or len(sealed) != 3506:
        raise ValueError('Cohort mismatch')
    groups = np.array([r['group_id'] for r in dev])
    folds = np.array([int(r['outer_fold']) for r in dev])
    if set(groups) & {r['group_id'] for r in rows if r['partition'] == 'sealed'}:
        raise ValueError('Sealed group overlap')
    for g in set(groups):
        if len(set(folds[groups == g])) != 1:
            raise ValueError('Outer group overlap')
    # Reuse the existing protected per-row development reader, no teacher imports.
    sys.path.insert(0, str(root.parent / 'fprior' / 'src'))
    from fprior.io import load_development_phenotypes
    y = load_development_phenotypes(data / 'partition_manifest.json', ids)
    samples = read_tsv(data / 'samples.tsv')
    mapping = {r['normalized_id']: int(r['sample_index']) for r in samples}
    sample_indices = np.array([mapping[i] for i in ids])
    variants = read_tsv(data / 'variants.tsv')
    chrom = np.array([r['CHROM'] for r in variants])
    pos = np.array([int(r['POS']) for r in variants])
    dosage_path = data / 'dosages.npy'
    stat = dosage_path.stat()
    dosage_identity = dict(sha256=config['dosages_sha256'], bytes=stat.st_size, mtime_ns=stat.st_mtime_ns)
    prior_path = (root / config['output'] / 'identity.json').resolve()
    prior = json.loads(prior_path.read_text()) if prior_path.exists() else {}
    if prior.get('inputs', {}).get('dosages.npy') != dosage_identity:
        digest = hashlib.sha256()
        with dosage_path.open('rb') as stream:
            while block := stream.read(8 * 1024 * 1024):
                digest.update(block)
        if digest.hexdigest() != config['dosages_sha256']:
            raise ValueError('Genotype identity mismatch')
    raw = np.load(dosage_path, mmap_mode='r')
    if raw.shape != (41346, 20087) or len(variants) != 41346:
        raise ValueError('Genotype shape mismatch')
    # Only permitted development columns are selected from the read-only source.
    x = np.empty((len(ids), len(variants)), dtype=np.float32)
    for start in range(0, len(variants), 1024):
        x[:, start:start + 1024] = raw[start:start + 1024, sample_indices].T
    identity = {name: records[name]['sha256'] for name in ['samples.tsv', 'variants.tsv', 'splits.tsv']}
    identity['dosages.npy'] = dosage_identity
    identity['partition_manifest.json'] = hashlib.sha256((data / 'partition_manifest.json').read_bytes()).hexdigest()
    return x, y, ids, groups, folds, chrom, pos, variants, identity


def context(x, y, train, test, chrom, pos, physical, config, directory, tag, rhos):
    started = time.perf_counter()
    xt, xv = x[train], x[test]
    mean, keep = fit_qc(xt, config)
    ld, neighbours, ld_summary, representatives = ld_order(xt, mean, keep, chrom, pos, physical, config)
    orders = {'physical': physical, 'ld': ld}
    write_json(directory / f'{tag}_qc_ld.json', dict(training_samples=len(train),
        validation_samples=len(test), qc_retained=int(keep.sum()),
        ld_representatives=representatives, ld=ld_summary))
    np.savez(directory / f'{tag}_preprocessing.npz', mean=mean, keep=keep,
             physical_order=physical, ld_order=ld, train=train, test=test)
    results = {}
    snapshots = {}
    ref_full = None
    for strategy, order in orders.items():
        n = len(train)
        kernel = torch.zeros((n, n), dtype=torch.float64, device=config['device'])
        cross = torch.zeros((len(test), n), dtype=torch.float64, device=config['device'])
        previous = set()
        for budget in config['budgets']:
            # Budget is the pre-QC prefix. No post-QC replacement in physical panels.
            selected = order[:budget]
            selected = selected[keep[selected]]
            if not len(selected):
                raise ValueError('Empty marker panel')
            if budget == config['budgets'][-1] and strategy == 'ld':
                if set(selected) != set(ref_full['indices']):
                    raise ValueError('Full panel strategy mismatch')
                results[(strategy, budget)] = ref_full
                continue
            added = np.array([i for i in selected if int(i) not in previous])
            if not previous <= set(map(int, selected)):
                raise ValueError('Non-nested marker panel')
            # Chunked multiplication bounds memory and accumulates unnormalised K.
            for start in range(0, len(added), 4096):
                chunk = added[start:start + 4096]
                zt = torch.as_tensor(centered(xt, mean, chunk), device=config['device'])
                zv = torch.as_tensor(centered(xv, mean, chunk), device=config['device'])
                kernel.addmm_(zt, zt.T)
                cross.addmm_(zv, zt.T)
                del zt, zv
            previous = set(map(int, selected))
            scale = float(torch.trace(kernel) / n)
            if scale <= 0:
                raise ValueError('Zero kernel scale')
            # Unit mean training diagonal and n*rho keep penalty scale comparable.
            normalized = kernel / scale
            prediction, residuals, iterations = solve_grid(normalized, cross / scale, y[train], rhos, config)
            record = dict(prediction=prediction, indices=selected.copy(), scale=scale,
                          residuals=residuals, iterations=iterations)
            results[(strategy, budget)] = record
            if budget == config['budgets'][-1]:
                ref_full = record
            if tag == 'outer':
                stats = panel_stats(selected, chrom, pos, neighbours, config['physical_bin_bp'])
                write_json(directory / f'{strategy}_{budget}_coverage.json', stats)
                np.save(directory / f'{strategy}_{budget}_indices.npy', selected)
                # Exact trace participation ratio; sampled entries for full-panel comparison.
                sample = np.linspace(0, n - 1, min(512, n), dtype=int)
                idx = torch.as_tensor(sample, device=kernel.device)
                snapshot = normalized[idx[:, None], idx[None, :]].cpu().numpy()
                snapshots[(strategy, budget)] = (snapshot, float(n * n / torch.sum(normalized ** 2)))
            log(f'{directory.name} {tag} {strategy} {budget}: retained={len(selected)}, '
                f'CG={max(iterations)}, residual={max(residuals):.2g}')
            del normalized
        del kernel, cross
    if tag == 'outer':
        snapshots[('ld', config['budgets'][-1])] = snapshots[('physical', config['budgets'][-1])]
        full = snapshots[('physical', config['budgets'][-1])][0]
        mask = ~np.eye(len(full), dtype=bool)
        diagnostics = []
        for (strategy, budget), (a, rank) in snapshots.items():
            diagnostics.append(dict(strategy=strategy, budget=budget,
                trace_participation_ratio=rank, diagnostic_training_samples=len(full),
                sampled_relative_frobenius=float(np.linalg.norm(a - full) / np.linalg.norm(full)),
                sampled_offdiagonal_correlation=float(np.corrcoef(a[mask], full[mask])[0, 1])))
        write_json(directory / 'relationship_diagnostics.json', diagnostics)
        # Emit an explicit full-panel alias; both strategies share identical fit.
        np.save(directory / f'ld_{config["budgets"][-1]}_indices.npy', ref_full['indices'])
        write_json(directory / f'ld_{config["budgets"][-1]}_coverage.json',
                   panel_stats(ref_full['indices'], chrom, pos, neighbours, config['physical_bin_bp']))
    write_json(directory / f'{tag}_cost.json', dict(seconds=time.perf_counter() - started,
        peak_host_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0))
    return results


def run_fold(fold, inputs, physical, config, output):
    x, y, ids, groups, folds, chrom, pos, variants, identity = inputs
    directory = output / f'fold_{fold}'
    directory.mkdir(exist_ok=True)
    train, test = np.flatnonzero(folds != fold), np.flatnonzero(folds == fold)
    inner = group_folds(groups[train], config['inner_folds'], config['seed'], f'inner:{fold}:full')
    write_tsv(directory / 'inner_splits.tsv', [dict(accession_id=ids[i], group_id=groups[i], inner_fold=int(f))
                                               for i, f in zip(train, inner)])
    started = time.perf_counter()
    torch.cuda.reset_peak_memory_stats() if torch.cuda.is_available() else None
    rhos = list(config['rho_grid'])
    contexts = {}
    while True:
        scores = {(s, b): np.zeros(len(rhos)) for s in ['physical', 'ld'] for b in config['budgets']}
        for k in range(config['inner_folds']):
            ti, vi = train[inner != k], train[inner == k]
            if set(groups[ti]) & set(groups[vi]):
                raise ValueError('Inner group overlap')
            result = context(x, y, ti, vi, chrom, pos, physical, config, directory, f'inner_{k}', rhos)
            for key, record in result.items():
                scores[key] += ((record['prediction'] - y[vi, None]) ** 2).sum(0)
        selected = {key: int(np.argmin(values)) for key, values in scores.items()}
        extend_low = any(i == 0 for i in selected.values()) and rhos[0] > config['rho_bounds'][0]
        extend_high = any(i == len(rhos) - 1 for i in selected.values()) and rhos[-1] < config['rho_bounds'][1]
        if not extend_low and not extend_high:
            break
        if extend_low:
            rhos.insert(0, max(config['rho_bounds'][0], rhos[0] / 10))
        if extend_high:
            rhos.append(min(config['rho_bounds'][1], rhos[-1] * 10))
        log(f'fold_{fold}: extending shared rho grid to {rhos[0]}..{rhos[-1]}')
    fitted_rhos = sorted({rhos[i] for i in selected.values()})
    result = context(x, y, train, test, chrom, pos, physical, config, directory, 'outer', fitted_rhos)
    records, predictions = [], []
    for (strategy, budget), record in result.items():
        selected_rho = rhos[selected[(strategy, budget)]]
        j = fitted_rhos.index(selected_rho)
        records.append(dict(strategy=strategy, budget=budget, actual_markers=len(record['indices']),
            rho=selected_rho, search_rhos=rhos, inner_mse=(scores[(strategy, budget)] / len(train)).tolist(),
            boundary_selected=selected[(strategy, budget)] in [0, len(rhos) - 1],
            scale=record['scale'], relative_residual=record['residuals'][j], iterations=record['iterations'][j]))
        predictions.extend(dict(accession_id=ids[i], group_id=groups[i], outer_fold=fold,
            strategy=strategy, budget=budget, actual_markers=len(record['indices']),
            observed=float(y[i]), prediction=float(pred), train_mean=float(y[train].mean()))
            for i, pred in zip(test, record['prediction'][:, j]))
    write_tsv(directory / 'predictions.tsv', predictions)
    write_json(directory / 'search.json', records)
    write_json(directory / 'complete.json', dict(status='COMPLETED', fold=fold,
        train_samples=len(train), test_samples=len(test), seconds=time.perf_counter() - started,
        peak_host_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        peak_cuda_allocated_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0))
    log(f'fold_{fold} completed in {time.perf_counter() - started:.1f}s')


def report(config, output):
    rows = []
    for fold in range(10):
        if not (output / f'fold_{fold}' / 'complete.json').exists():
            raise ValueError(f'Incomplete outer fold {fold}')
        rows += read_tsv(output / f'fold_{fold}' / 'predictions.tsv')
    write_tsv(output / 'predictions.tsv', rows)
    metrics, arrays = [], {}
    for strategy in ['physical', 'ld']:
        for budget in config['budgets']:
            selected = sorted([r for r in rows if r['strategy'] == strategy and int(r['budget']) == budget],
                              key=lambda r: r['accession_id'])
            if len(selected) != 14020 or len({r['accession_id'] for r in selected}) != 14020:
                raise ValueError('Prediction coverage mismatch')
            y = np.array([float(r['observed']) for r in selected])
            pred = np.array([float(r['prediction']) for r in selected])
            error = (y - pred) ** 2
            arrays[(strategy, budget)] = (selected, error)
            metrics.append(dict(strategy=strategy, budget=budget, pooled_mse=float(error.mean()),
                pcc=float(np.corrcoef(y, pred)[0, 1]),
                r2=float(1 - error.sum() / ((y - y.mean()) ** 2).sum()),
                mae=float(np.abs(y - pred).mean()),
                actual_markers_min=min(int(r['actual_markers']) for r in selected),
                actual_markers_max=max(int(r['actual_markers']) for r in selected)))
    reference = arrays[('physical', config['budgets'][0])][0]
    groups = np.array([r['group_id'] for r in reference])
    folds = np.array([int(r['outer_fold']) for r in reference])
    unique, inverse = np.unique(groups, return_inverse=True)
    counts = np.bincount(inverse)
    rng = np.random.default_rng(config['seed'])
    comparisons = []
    pairs = []
    for strategy in ['physical', 'ld']:
        pairs += [((strategy, a), (strategy, b), 'adjacent') for a, b in zip(config['budgets'][:-1], config['budgets'][1:])]
        pairs.append(((strategy, config['budgets'][0]), (strategy, config['budgets'][-1]), 'full_vs_sparse'))
    pairs += [(('physical', b), ('ld', b), 'strategy') for b in config['budgets'][:-1]]
    for left, right, kind in pairs:
        a, ea = arrays[left]
        b, eb = arrays[right]
        if [r['accession_id'] for r in a] != [r['accession_id'] for r in b]:
            raise ValueError('Paired identity mismatch')
        diff = eb - ea
        group_delta = np.bincount(inverse, weights=diff)
        group_base = np.bincount(inverse, weights=ea)
        bootstrap, gains = [], []
        # Stratified group resampling keeps fold proportions and all group members.
        strata = [np.unique(inverse[folds == f]) for f in range(10)]
        for _ in range(config['bootstrap_replicates']):
            draw = np.concatenate([rng.choice(s, len(s), replace=True) for s in strata])
            bootstrap.append(float(group_delta[draw].sum() / counts[draw].sum()))
            gains.append(float(-group_delta[draw].sum() / group_base[draw].sum()))
        comparisons.append(dict(kind=kind, reference_strategy=left[0], reference_budget=left[1],
            candidate_strategy=right[0], candidate_budget=right[1], delta_mse=float(diff.mean()),
            relative_mse_gain=float(-diff.sum() / ea.sum()),
            delta_mse_ci95=np.quantile(bootstrap, [.025, .975]).tolist(),
            relative_gain_ci95=np.quantile(gains, [.025, .975]).tolist(),
            improved_outer_folds=sum(float(diff[folds == f].mean()) < 0 for f in range(10))))
    write_tsv(output / 'metrics.tsv', metrics)
    write_json(output / 'comparisons.json', comparisons)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for strategy in ['physical', 'ld']:
        m = [a for a in metrics if a['strategy'] == strategy]
        ax.plot([a['budget'] for a in m], [a['pooled_mse'] for a in m], 'o-', label=strategy)
    ax.set(xlabel='Pre-QC SNP budget', ylabel='Pooled out-of-fold MSE (% dry weight squared)')
    ax.legend(); fig.tight_layout(); fig.savefig(output / 'density_performance.png', dpi=180)
    fig.savefig(output / 'density_performance.pdf'); plt.close(fig)
    costs = [json.loads((output / f'fold_{f}' / 'complete.json').read_text()) for f in range(10)]
    write_json(output / 'costs.json', costs)
    write_json(output / 'complete.json', dict(status='COMPLETED', development_samples=14020,
        outer_folds=10, workflows=10, unique_full_panel_fit=True,
        sealed_phenotypes_read=False, sealed_genotypes_used_for_ld=False,
        bootstrap_scope='fixed predictions; fold-stratified whole groups; no refitting'))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['profile', 'evaluate', 'report'])
    parser.add_argument('--config', default=str(Path(__file__).with_name('config.json')))
    args = parser.parse_args()
    path = Path(args.config).resolve()
    config = json.loads(path.read_text())
    root = path.parent
    output = (root / config['output']).resolve()
    from threadpoolctl import threadpool_limits
    threadpool_limits(config['blas_threads'])
    torch.set_num_threads(config['blas_threads'])
    if config['device'] == 'cuda' and not torch.cuda.is_available():
        raise RuntimeError('Configured CUDA unavailable')
    if args.action == 'report':
        report(config, output); return
    inputs = load_inputs(config, root)
    source_files = [Path(__file__), root / 'preprocessing' / 'panels.py', root.parent / 'fprior' / 'src' / 'fprior' / 'io.py']
    identity = dict(config=config, inputs=inputs[-1],
        implementation={str(p.relative_to(root.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files})
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'identity.json').exists():
        if json.loads((output / 'identity.json').read_text()) != identity:
            raise ValueError('Cannot resume changed experiment identity')
    else:
        write_json(output / 'identity.json', identity)
        write_json(output / 'environment.json', dict(python=sys.version, numpy=np.__version__,
            torch=torch.__version__, gpu=torch.cuda.get_device_name() if torch.cuda.is_available() else None))
    x, y, ids, groups, folds, chrom, pos, variants, _ = inputs
    physical = physical_order(chrom, pos, config['seed'], config['physical_bin_bp'])
    write_tsv(output / 'marker_catalog.tsv', [dict(index=i, **v) for i, v in enumerate(variants)])
    for fold in ([0] if args.action == 'profile' else range(10)):
        if (output / f'fold_{fold}' / 'complete.json').exists():
            log(f'fold_{fold}: reusing completed fit with matching identity'); continue
        log(f'fold_{fold}: starting full-training grouped nested evaluation')
        run_fold(fold, inputs, physical, config, output)
    if args.action == 'evaluate':
        report(config, output)


if __name__ == '__main__':
    main()
