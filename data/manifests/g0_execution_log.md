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

## TASK02

- Current evidence: the historical VCF audit validates structural records and only a 2,000-site GT window; it does not freeze record-level FILTER/FORMAT/GT behavior over the full VCF.
- Current hypothesis: a separate streaming preflight can verify source identity, count exact FILTER and FORMAT layouts, dynamically locate GT, and fail only when an accepted structural record contains an unsupported GT state.
- Minimal verification: exercise PASS/`.` acceptance, named FILTER exclusion, dynamic `DP:GT:GQ`, missing/duplicate GT, malformed widths, and every unsupported GT class on a synthetic gzipped VCF before the real full pass.

### TASK02 performance correction

- Current evidence: the first real preflight was genuine (28 minutes, 23.9 GB cumulative source reads, one CPU core saturated) but its per-cell Python loop made its projected runtime unsuitable for a 12.15-billion-cell source.
- Current hypothesis: collapsing identical per-record sample-field strings with `Counter` preserves every GT state count and blocker while eliminating repeated Python token parsing for the overwhelmingly repeated hard-call cells.
- Minimal verification: retain synthetic output equality, then restart exactly one full preflight; no phenotype or prediction operation is introduced.

### TASK02 input-backend correction

- Current evidence: this host has 24 CPU cores and `pigz`, but no `bcftools`, `tabix`, or VCF region index; a GPU is not an appropriate accelerator for gzipped tabular genotype parsing.
- Current hypothesis: an explicit `pigz` streaming backend can parallelize release decompression without creating a plaintext VCF, while the existing Python semantic state machine remains the single source of policy truth.
- Minimal verification: preserve the stdlib-gzip default for fixtures/CI; add a controlled backend test with a fake `pigz` executable and require nonzero decompressor exits to fail clearly.

### TASK02 HTSlib correction

- Current evidence: the local `audit-1001g` conda environment already contains bcftools/HTSlib/tabix 1.24, but the VCF has no `.csi`/`.tbi` index; the prior shell PATH check was incomplete.
- Current hypothesis: a CSI created by HTSlib will permit chromosome-region reads and make the full semantic preflight and later build CPU-parallel without modifying the source VCF.
- Minimal verification: validate the compressed artifact with `bgzip -t`, build a CSI with bcftools, then require `bcftools index --stats` to enumerate all five nuclear contigs before using any region stream.

### TASK02 indexed parallelization

- Current evidence: HTSlib built a 93 KiB CSI and reports Chr1–Chr5 record spans; region reads are now authoritative and available.
- Current hypothesis: source verification once, followed by five indexed `bcftools view -r` streams with bounded per-worker statistics, preserves the same policy result while using the available CPU threads.
- Minimal verification: the worker must reuse the same stream state machine as the serial fixture path; merger tests must reject mismatched sample order and preserve additive counters.

### TASK02 CPU-saturation correction

- Current evidence: the five chromosome workers ran correctly but their Python parsing shared one interpreter GIL, leaving approximately 6 of 24 CPU cores busy despite a multithreaded HTSlib invocation.
- Current hypothesis: splitting the CSI-indexed TAIR10 coordinate ranges into 12 disjoint regions and using independent processes will preserve exact source coverage while allowing the parser and HTSlib reader stages to use the available CPU capacity.
- Minimal verification: generate contiguous, gap-free intervals covering each chromosome exactly once; retain the same stream state machine and additive/sample-order merger test before restarting the full scan.

## TASK03

- Current evidence: the full VCF preflight now passes with 10,707,430 accepted records, and the existing source-derived biallelic catalog supplies the required ordered `(chromosome, position, REF, ALT, source_record_index)` rows with a recorded SHA256.
- Current hypothesis: a five-chromosome HTSlib query can write non-overlapping column ranges of one preallocated `int8` `.npy` memmap, while each worker stream-checks its query keys against the catalog; this avoids a Python list/float matrix and retains the original VCF record index without a second semantic VCF parser.
- Minimal verification: test fixed-width GT-byte conversion, catalog/query key mismatch rejection, FILTER/structural query expression, atomic partial-to-final finalization, and a real VCF interval against independent direct extraction before any full build.
