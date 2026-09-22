import importlib.util
from pathlib import Path


def module():
    path = Path(__file__).parents[1] / "scripts" / "refine_nearclone_exact_edges.py"
    spec = importlib.util.spec_from_file_location("nearclone_refine", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_conservative_screen_kl_bound_is_small_at_registered_cutoffs():
    mod = module()
    bound = mod.bernoulli_kl(.002, .001)
    assert bound > 0
    # 250K hash-selected states make a true d=.001 pair crossing .002
    # extraordinarily unlikely under the declared Bernoulli approximation.
    assert 1135 * 1134 // 2 * __import__("math").exp(-250_000 * bound) < 1e-30
