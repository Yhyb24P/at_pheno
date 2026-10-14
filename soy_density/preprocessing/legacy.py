"""Read the existing BiMSGP tables without changing its cohort or marker order."""
from pathlib import Path
import csv
import hashlib
import numpy as np


def load_legacy(config, root):
    import pyarrow.parquet as pq
    paths = {name: (root / config[name]).resolve() for name in
             ['genotypes', 'phenotypes', 'folds', 'bimsgp_results', 'upstream_gblup_results']}
    source = pq.ParquetFile(paths['genotypes'])
    markers = source.schema_arrow.names[1:]
    table = source.read()
    ids = table.column(0).to_pylist()
    x = np.empty((len(ids), len(markers)), dtype=np.float32)
    for j in range(len(markers)):
        x[:, j] = table.column(j + 1).to_numpy()
    del table
    if x.shape != (config['expected_samples'], config['expected_markers']) or len(set(ids)) != len(ids):
        raise ValueError('Legacy genotype shape or identity mismatch')
    # Only the legacy phenotype table is used; no official-cohort labels are opened.
    table = pq.read_table(paths['phenotypes'], columns=['acid', config['trait']],
                          filters=[('acid', 'in', ids)])
    label_ids = table.column('acid').to_pylist()
    labels = table.column(config['trait']).to_pylist()
    if len(set(label_ids)) != len(label_ids):
        raise ValueError('Duplicate legacy phenotype IDs')
    mapping = dict(zip(label_ids, labels))
    y = np.array([float(mapping[i]) for i in ids], dtype=np.float32).astype(np.float64)
    if not np.isfinite(y).all() or not np.isfinite(x).all():
        raise ValueError('Legacy cohort differs from the finite complete-case paper cohort')
    with paths['folds'].open() as stream:
        rows = list(csv.DictReader(stream))
    assignment = {r['sample_id']: int(r['fold']) - 1 for r in rows}
    if set(assignment) != set(ids) or len(rows) != len(ids) or set(assignment.values()) != set(range(10)):
        raise ValueError('Legacy fold identity mismatch')
    folds = np.array([assignment[i] for i in ids])
    # Hash only these small existing input files; no checkpoint loading or data audit.
    identity = {name: dict(path=config[name], bytes=p.stat().st_size,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for name, p in paths.items()}
    return x, y, ids, folds, markers, identity


def fit_transform(x, train, standardize):
    mean = x[train].mean(axis=0, dtype=np.float64)
    if standardize:
        # Population variance agrees with StandardScaler on complete data.
        variance = np.mean((x[train] - mean) ** 2, axis=0)
        scale = np.sqrt(variance)
        scale[scale == 0] = 1
    else:
        scale = np.ones(x.shape[1])
    return mean, scale


def transform(x, rows, columns, mean, scale):
    return (x[np.ix_(rows, columns)].astype(np.float64) - mean[columns]) / scale[columns]
