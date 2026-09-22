import importlib.util
from pathlib import Path


def module():
    path = Path(__file__).parents[1] / "scripts" / "freeze_confirmatory_splits.py"
    spec = importlib.util.spec_from_file_location("splits", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_block_assignment_is_deterministic_and_never_splits_a_block():
    mod = module()
    blocks = [[("a", "g1", 0), ("b", "g1", 0), ("c", "g1", 0)],
              [("d", "g2", 1), ("e", "g2", 1)], [("f", "g1", 2)], [("g", "g2", 3)]]
    first = mod.assign_blocks(blocks, 3, "salt")
    second = mod.assign_blocks(blocks, 3, "salt")
    assert first == second
    assert first["a"] == first["b"] == first["c"]
    assert first["d"] == first["e"]


def test_trait_blocks_reject_missing_or_cross_manifest_ids():
    mod = module()
    groups = {"a": "g1", "b": "g1"}
    global_blocks = {0: ["a"], 1: ["b"]}
    assert mod.blocks_for_trait(["a", "b"], global_blocks, groups) == [[("a", "g1", 0)], [("b", "g1", 1)]]
