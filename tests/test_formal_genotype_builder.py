import gzip
import importlib.util
from pathlib import Path

import numpy as np
import pytest


def builder_module():
    path = Path(__file__).parents[1] / "scripts" / "build_formal_genotype_v1.py"
    spec = importlib.util.spec_from_file_location("formal_builder", path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def test_fixed_width_query_gt_parser_converts_all_contract_states():
    mod = builder_module()
    # %CHROM/%POS/%REF/%ALT then [%GT] exactly as emitted by bcftools query.
    line = b"1\t55\tC\tT\t0/0\t0|1\t1/1\t./.\n"
    key, dosage = mod.parse_query_line(line, 4)
    assert key == ("1", 55, "C", "T")
    assert dosage.dtype == np.int8
    assert dosage.tolist() == [0, 1, 2, -1]


@pytest.mark.parametrize("line, message", [
    (b"1\t55\tC\tT\t0/0\t0/.\n", "partial missing"),
    (b"1\t55\tC\tT\t0/0\t2/2\n", "unsupported allele"),
    (b"1\t55\tC\tT\t0/0\t0/0/0\n", "payload has"),
    (b"1\t55\tC\tT\t0/0\t0:0\n", "non-diploid separator"),
])
def test_fixed_width_query_gt_parser_rejects_semantic_drift(line, message):
    mod = builder_module()
    with pytest.raises(ValueError, match=message):
        mod.parse_query_line(line, 2)


def test_catalog_rows_rejects_ordering_and_preserves_source_indices(tmp_path):
    mod = builder_module()
    catalog = tmp_path / "catalog.tsv.gz"
    with gzip.open(catalog, "wt") as stream:
        stream.write("chromosome\tposition\tref\talt\tsource_record_index\n")
        stream.write("1\t55\tC\tT\t1\n1\t56\tT\tA\t2\n")
    assert list(mod.catalog_rows(catalog)) == [("1", 55, "C", "T", 1), ("1", 56, "T", "A", 2)]
    with gzip.open(catalog, "wt") as stream:
        stream.write("chromosome\tposition\tref\talt\tsource_record_index\n")
        stream.write("1\t56\tT\tA\t2\n1\t55\tC\tT\t1\n")
    with pytest.raises(ValueError, match="strictly ordered"):
        list(mod.catalog_rows(catalog))


def test_write_chromosome_rejects_catalog_query_key_mismatch(tmp_path, monkeypatch):
    mod = builder_module()
    catalog = tmp_path / "catalog.tsv.gz"
    with gzip.open(catalog, "wt") as stream:
        stream.write("chromosome\tposition\tref\talt\tsource_record_index\n1\t55\tC\tT\t1\n")
    matrix_path = tmp_path / "genotypes.npy.partial"
    np.lib.format.open_memmap(matrix_path, mode="w+", dtype=np.int8, shape=(2, 1)).flush()

    class FakeProcess:
        def __init__(self):
            self.stdout = FakeStream([b"1\t56\tC\tT\t0/0\t1/1\n"])
            self.stderr = None
        def wait(self):
            return 0

    class FakeStream:
        def __init__(self, lines):
            self._lines = lines
        def __iter__(self):
            return iter(self._lines)
        def close(self):
            pass

    monkeypatch.setattr(mod.subprocess, "Popen", lambda *args, **kwargs: FakeProcess())
    with pytest.raises(ValueError, match="key mismatch"):
        mod.write_chromosome("source.vcf.gz", "bcftools", catalog, matrix_path, "1", 0, 1, 2)
