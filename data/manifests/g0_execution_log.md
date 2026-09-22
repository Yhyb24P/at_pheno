# G0 execution log

## TASK00

- Current evidence: `HEAD` and `origin/master` are both `d8c219b0913d8930d48477551508bd214df96db8`; the worktree is clean; `uv run pytest -q` reports 106 passed; source VCF/HDF5 and source manifests exist; 338 GiB disk and 83 GiB available RAM were observed.
- Current hypothesis: the fixed G0 baseline is reproducible and the host has sufficient capacity for the formal int8 build.
- Minimal verification: record reproducible runtime/storage facts in `stage_g0_execution_start.json`; do not inspect phenotype values or run prediction.

## TASK01

- Current evidence: `load_dataset()` accepts any numeric dtype with values 0/1/2 (and NaN for floating dtypes), even after formal mode has passed provenance validation; this does not enforce the frozen formal int8 storage contract.
- Current hypothesis: passing an explicit `formal=True` storage context from `run()` to `load_dataset()` can enforce `np.int8` and `{-1,0,1,2}` without changing the pilot float-plus-NaN path.
- Minimal verification: add exact dtype/value fixture tests for formal int8, float32, int16, uint8, and invalid int8; retain the existing pilot float fixture.

### TASK01 test-fixture correction

- Current evidence: the new formal dtype gate correctly rejects two pre-existing formal CLI fixtures because `demo()` deliberately produces float32 pilot storage.
- Current hypothesis: the formal CLI tests must explicitly convert only their synthetic fixture matrices to `int8` with `-1` missing before constructing provenance hashes; production pilot behavior remains unchanged.
- Minimal verification: preserve the same provenance/gate assertions and re-run the complete test suite.
