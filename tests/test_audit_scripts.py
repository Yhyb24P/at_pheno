import importlib.util
from pathlib import Path

import numpy as np


def script_module(name):
    path = Path(__file__).parents[1] / "scripts" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_binary_scan_uses_explicit_minority_state_rule():
    scan = script_module("scan_hdf5_binary")
    histogram = np.zeros(12, dtype=int)  # n=11; states 5 and 6 pass, 7 fails.
    histogram[5], histogram[6], histogram[7] = 3, 4, 9
    assert scan.minor_state_ge(histogram, 11, 5) == 7
    counts, carriers = scan.summarize_binary_block(np.array([[0, 1, 1], [0, 2, 1]]), 3)
    assert counts == {"zero": 2, "one": 3, "other": 1}
    assert carriers[2] == 1 and carriers[1] == 1


def test_feasibility_counts_exclude_nonfinite_values():
    registry = script_module("scan_arapheno_feasibility")
    values = [{"accession_id": 1, "phenotype_value": 2.0},
              {"accession_id": 1, "phenotype_value": float("nan")},
              {"accession_id": 2, "phenotype_value": None},
              {"accession_id": 3, "phenotype_value": 1.0}]
    counts, _, _, _ = registry.phenotype_counts(values, {"1", "2"})
    assert counts["raw_row_n"] == 4 and counts["finite_row_n"] == 2
    assert counts["panel_intersection_n"] == 2 and counts["finite_panel_intersection_n"] == 1


def test_kinship_component_diagnostics_expose_single_linkage_chains():
    audit = script_module("audit_kinship_blocks")
    matrix = np.array([[1, .9, .1], [.9, 1, .9], [.1, .9, 1]])
    groups = audit.components(matrix, .8)
    stats = audit.component_diagnostics(matrix, groups, .8)
    assert groups == [[0, 1, 2]]
    assert stats[0]["min_pairwise_similarity"] == .1
    assert stats[0]["edge_density"] == 2/3
