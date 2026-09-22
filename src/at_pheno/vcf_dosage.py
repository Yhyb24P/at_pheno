"""VCF GT hard-call conversion into the formal int8 ALT-dosage matrix.

Formal genotype storage contract (frozen 2026-09): the dataset matrix is
stored as ``int8`` with shape ``n_samples x n_variants`` and value set
``{-1, 0, 1, 2}`` — ``-1`` missing, ``0/1/2`` ALT dosage.  No
float32 full matrix is ever materialized: the converter processes one
variant row (``n_samples`` GT cells) at a time and writes int8
directly; the caller stacks rows into the storage layout.

Only hard biallelic calls are supported.  Phase separators (``|``)
are accepted since they carry no dosage difference; every other GT
token (multi-allelic ``i/1`` states, half-calls, single-token cells,
empty cells) is rejected by name.
"""

import numpy as np

_INT8_DOSAGE = np.int8


def _cell_to_dosage(cell, sample_index):
    """One GT cell -> an int8 ALT dosage in {-1, 0, 1, 2}.

    ``"/" and "|" cells with state ``0``/``1`` map by ALT state count
    (``0/0`` -> 0, ``0/1`` -> 1, ``1/1`` -> 2); a fully missing
    genotype (``.`` in both states) maps to ``-1``.  Raises
    ValueError naming the offending cell and sample index otherwise.
    """
    if not isinstance(cell, str) or not cell:
        raise ValueError(
            f"Unsupported GT token for ALT-dosage int8 conversion: cell "
            f"of sample #{sample_index} is {cell!r}; every sample must "
            "carry a hard diploid cell")
    if "/" in cell:
        tokens = cell.split("/")
    else:
        tokens = cell.split("|")
    if len(tokens) != 2:
        raise ValueError(
            f"Unsupported GT token for ALT-dosage int8 conversion: cell "
            f"{cell!r} of sample #{sample_index} is not a two-state "
            "diploid cell (got " + ", ".join(repr(tokens)) + ")")
    left, right = tokens
    for token in (left, right):
        if token not in {"0", "1", "."}:
            raise ValueError(
                f"Unsupported GT token for ALT-dosage int8 conversion: "
                f"state {token!r} in cell {cell!r} of sample "
                f"#{sample_index} is neither hard call 0/1 nor missing . "
                "(multi-allelic states are out of the -1/0/1/2 contract)")
    if left == "." and right == ".":
        return -1
    if left == "." or right == ".":
        raise ValueError(
            f"Half-call GT token {cell!r} of sample #{sample_index} is "
            "unsupported by the int8 hard-call contract: -1 is reserved "
            "for a fully missing genotype")
    return int(left != "0") + int(right != "0")


def gt_dosage_row(cells):
    """One variant row of ``n_samples`` VCF-ordered GT cells -> int8.

    Returns an ``np.int8`` array of shape ``(n_samples,)``.
    """
    row = np.empty(len(cells), dtype=_INT8_DOSAGE)
    for index, cell in enumerate(cells):
        row[index] = _cell_to_dosage(cell, index)
    return row


def gt_cells_to_int8_dosage(cells_by_variant):
    """GT rows -> the formal storage layout int8 ``(n_samples, n_variants)``.

    ``cells_by_variant`` is a sequence of variant rows, each a
    sequence of GT cell strings in exact VCF header sample order.
    The output materializes directly as the int8 storage matrix —
    no wider float dtype ever appears in the path.
    """
    sample_by_sample = [gt_dosage_row(cells) for cells in cells_by_variant]
    if not sample_by_sample:
        return np.empty((0, 0), dtype=_INT8_DOSAGE)
    n_samples = len(sample_by_sample[0])
    for row in sample_by_sample:
        if len(row) != n_samples:
            raise ValueError("GT rows disagree on the number of samples")
    matrix = np.stack(sample_by_sample).T
    return np.ascontiguousarray(matrix, dtype=_INT8_DOSAGE)
