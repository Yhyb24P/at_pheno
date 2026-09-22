import importlib.util
from pathlib import Path

import h5py
import numpy as np


def module():
    path = Path(__file__).parents[1] / "scripts" / "freeze_nearclone_blocks.py"
    spec = importlib.util.spec_from_file_location("nearclone", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_hash_panel_is_deterministic_physically_sorted_and_has_chrom_coverage():
    mod = module()
    positions = np.array([1, 4, 9, 2, 5, 3, 7, 11, 6, 8], dtype=np.int32)
    regions = [(0, 3), (3, 5), (5, 7), (7, 8), (8, 10)]
    a, per_chrom = mod.select_panel(positions, regions, 10)
    b, _ = mod.select_panel(positions, regions, 10)
    assert a == b
    assert all(per_chrom.values())
    assert [(row[1], row[2]) for row in a] == sorted((row[1], row[2]) for row in a)


def test_exact_binary_hamming_and_connected_components(tmp_path):
    mod = module()
    path = tmp_path / "binary.h5"
    # samples 0/1 differ at one selected marker; sample 2 differs at all.
    with h5py.File(path, "w") as handle:
        handle.create_dataset("snps", data=np.array([[0, 0, 1], [1, 1, 0], [0, 1, 1], [1, 1, 0]], dtype=np.int8))
    with h5py.File(path, "r") as handle:
        distance = mod.hamming_counts(handle, [0, 1, 2, 3], block=2)
    assert distance.tolist() == [[0, 1, 4], [1, 0, 3], [4, 3, 0]]
    assert mod.components(distance, 1) == [[0, 1], [2]]
