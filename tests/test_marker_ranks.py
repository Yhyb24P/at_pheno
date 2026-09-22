import importlib.util
from pathlib import Path
import sys

import numpy as np


def module():
    path = Path(__file__).parents[1] / "scripts" / "freeze_marker_ranks.py"
    spec = importlib.util.spec_from_file_location("ranks", path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = result
    spec.loader.exec_module(result)
    return result


def compact(mod):
    rows = []
    for chrom in range(1, 6):
        for position in (100, 200, 300):
            rows.append((chrom, position, 0, 1, len(rows) + 1))
    return np.array(rows, dtype=mod.DTYPE)


def test_rank_is_deterministic_permutation_and_physically_spread(tmp_path):
    mod = module()
    values = compact(mod)
    first, strata = mod.physical_rank(values, "salt", round_size=2)
    second, _ = mod.physical_rank(values, "salt", round_size=2)
    third, _ = mod.physical_rank(values, "other", round_size=2)
    path = tmp_path / "compact.npy"
    np.save(path, values)
    parallel, _ = mod.physical_rank_parallel(path, "salt", workers=2, round_size=2)
    assert first.tolist() == second.tolist()
    assert sorted(first.tolist()) == list(range(len(values)))
    assert strata == 5
    assert len(set(first[:strata])) == strata
    assert first.tolist() != third.tolist()
    assert parallel.tolist() == first.tolist()
