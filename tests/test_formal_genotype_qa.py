import csv
import importlib.util
from pathlib import Path

import numpy as np
import pytest


def qa_module():
    path = Path(__file__).parents[1] / "scripts" / "audit_formal_genotype_v1.py"
    spec = importlib.util.spec_from_file_location("formal_qa", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def write_variants(tmp_path, rows, compact_rows):
    tsv = tmp_path / "variants.tsv"
    with tsv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["chromosome", "position", "ref", "alt", "source_record_index"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    mod = qa_module()
    compact = tmp_path / "variants_compact.npy"
    array = np.lib.format.open_memmap(compact, mode="w+", dtype=mod.COMPACT_DTYPE, shape=(len(compact_rows),))
    array[:] = compact_rows
    array.flush()
    return tsv, compact


def test_independent_matrix_qa_rejects_out_of_contract_int8(tmp_path):
    mod = qa_module()
    matrix = tmp_path / "genotypes.npy"
    np.save(matrix, np.array([[0, 2], [-1, -2]], dtype=np.int8))
    with pytest.raises(ValueError, match="outside -1/0/1/2"):
        mod.audit_matrix(matrix, 2, 2, block=1)


def test_independent_variant_qa_checks_compact_and_source_order(tmp_path):
    mod = qa_module()
    rows = [{"chromosome": "1", "position": "55", "ref": "C", "alt": "T", "source_record_index": "1"},
            {"chromosome": "1", "position": "56", "ref": "T", "alt": "A", "source_record_index": "2"}]
    tsv, compact = write_variants(tmp_path, rows, [(1, 55, 1, 3, 1), (1, 56, 3, 0, 2)])
    report = mod.audit_variants(tsv, compact, 2)
    assert report["strict_order"] is True and report["duplicate_chromosome_position"] == 0
    array = np.load(compact, mmap_mode="r+")
    array[1]["position"] = 99
    array.flush()
    with pytest.raises(ValueError, match="variants_compact mismatch"):
        mod.audit_variants(tsv, compact, 2)
