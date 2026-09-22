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

## TASK05

- Current evidence: `BUILD_COMPLETE.json` exists only after the formal builder's value-domain/shape/hash checks; however its preliminary QC is intentionally insufficient because it was produced in the build process rather than by an independent reread of final artifacts.
- Current hypothesis: a post-build QA reader can stream final matrix blocks and `variants.tsv`, compare the compact runtime index, verify VCF/HDF5 sample order, and recompute the deterministic 2,000-site HDF5 orientation check without reading phenotypes.
- Minimal verification: fixture-test matrix-domain and variant-order failures, then run the QA only against the final `BUILD_COMPLETE` artifact and emit its own SHA256-bound report.

### TASK06 provenance-contract correction

- Current evidence: selected local `raw_values/*.json` bytes do not equal `trait_resolution_v1.tsv.values_sha256`; inspection of the acquisition code proves that field is the remote HTTP response digest, whereas the local files were intentionally parsed and re-serialized with `json.dumps(...)+"\\n"`. Remote re-fetch currently fails with a connection reset, so byte identity cannot be re-established by a new download.
- Current hypothesis: these are distinct hash domains, not evidence of a changed phenotype table. The freeze record must retain the registry's remote-response hash, add a separately recomputed local-file SHA256, and require semantic invariants (finite unique accession count, no selected duplicates, trait name/status) before it can use the local snapshot.
- Minimal verification: do not overwrite the registry hash; emit both hashes and their explicit non-equality, validate every selected raw snapshot against its registry counts/status, and make the remote re-fetch failure an auditable limitation rather than a silent pass.

## TASK07

- Current evidence: the official IBS kinship matrix has an uncalibrated native scale around 6–7 and cannot define a Hamming cutoff. The official binary HDF5 has complete 0/1 states, unique positions, and the exact 1,135 final sample order.
- Current hypothesis: a hash-selected 250K physical-coordinate panel plus blockwise binary Gram accumulation yields exact all-pair Hamming distances without REF/ALT orientation or phenotype values; connected components and their internal edge diagnostics can then enforce the frozen chain blocker rule.
- Minimal verification: validate deterministic panel selection and hand-calculated binary Hamming distances on fixtures; require every chromosome to appear and compare final HDF5 sample order to the formal genotype sample manifest before scanning the full panel.

### TASK07 chain-evidence check

- Current evidence: the selected-panel primary component has 13 members, 15/78 internal edges (density 0.1923), and a maximum internal Hamming distance 0.001292 versus the 0.001 edge cutoff; this meets the specified low-density chain condition, but the prose `>>` merits a direct full-HDF5 confirmation before treating it as a durable design blocker.
- Current hypothesis: recomputing only the affected component's pairwise binary Hamming distances across all 10,709,949 HDF5 states can confirm or falsify the chain pattern without modifying the frozen 250K primary definition.
- Minimal verification: load only the 13 affected columns in bounded HDF5 chunks; record full-panel maximum distance and edge density at the same 0.001 numerical cutoff, with no phenotype access.

### TASK07 chain-topology characterization

- Current evidence: full-HDF5 confirmation establishes a sparse 13-member component but summary metrics alone do not show whether clique-based blocking, pairwise edge constraints, or simple component blocking would yield materially different sample units.
- Current hypothesis: recording the exact full-HDF5 threshold-edge topology and maximum clique sizes will make the methodological choice explicit; it is diagnostic only and cannot silently replace the frozen connected-components rule.
- Minimal verification: re-use the same full-HDF5 bounded computation, emit only accession IDs and pair distances below the existing cutoff, and label all alternative topology summaries non-canonical.

### TASK07 sampling-correction investigation

- Current evidence: the 250K rule yields 16 primary edges and 2,686 edges at the already registered 0.002 sensitivity cutoff, while full-HDF5 examination of the apparent 13-member component removes several 250K-only edges. At a true distance of 0.001, a 250K binary sample has standard error about 0.000063; using 0.002 as a candidate screen is therefore more than fifteen standard errors above the primary boundary.
- Current hypothesis: retain the frozen 250K coordinate panel and 0.001 primary *edge definition*, but use its registered 0.002 sensitivity graph only to nominate candidate pairs, then compute exact all-HDF5 Hamming distance for every candidate before the 0.001 connected-component operation. This removes threshold sampling noise without changing the biological threshold or using phenotype information.
- Minimal verification: materialize the 0.002 candidate pairs, exact-check all of them in bounded HDF5 chunks, record the conservative screen's model-based false-negative bound, and verify that the full exact primary graph is a subset of candidates before considering any replacement of the provisional blocker.

## TASK08

- Current evidence: exact-refined full-HDF5 blocks pass the chain rule; the selected trait snapshots have unique published accession IDs, but the prior cohort-membership artifact does not include 703/704/705/748.
- Current hypothesis: the split generator can read only accession IDs from each registry-bound raw snapshot, intersect them with exact blocks and published group labels, and validate the pre-frozen finite-panel count without inspecting any `phenotype_value` field.
- Minimal verification: fixture-test deterministic block assignment, outer/inner block non-crossing, group-only balancing, and failure on a trait ID count that disagrees with the frozen panel before emitting canonical, alternate-salt, IID, or LOGO manifests.

## TASK09

- Current evidence: the final variants compact index is physically ordered and SHA256-bound; 10.7M Python marker strings plus one global hash sort would be needlessly memory-heavy and would not guarantee early physical spread.
- Current hypothesis: sorting SHA256 ranks inside each contiguous chromosome × 1-Mb stratum, then emitting fixed-size rounds across strata, creates an exact prefix-nested rank with physical balance and only bounded per-stratum hash arrays.
- Minimal verification: fixture-test salt determinism, prefix nesting, and first-round multi-stratum coverage; bind each generated int32 rank artifact to the final compact-index SHA256 and freeze all QC/density/endpoint/bootstrap rules without executing predictions.

### TASK09 rank-throughput correction

- Current evidence: three 10.7M SHA256 rank artifacts require roughly 32M independent short hashes; a serial Python loop would underuse the available CPU although each chromosome-bin stratum is independent before deterministic interleaving.
- Current hypothesis: worker processes can rank independent contiguous compact-index strata and return them in fixed stratum order; the parent retains the prescribed round-robin interleave, so concurrency cannot affect rank identity.
- Minimal verification: retain a serial fixture oracle, compare parallel and serial rank arrays byte-for-byte on a small compact array, then use bounded worker count for the formal artifacts.

## TASK10

- Current evidence: formal build, QA, trait/split/near-clone/rank/analysis artifacts exist and no real phenotype prediction has been run; `runs/smoke_iid` and `runs/smoke_group` contain explicitly synthetic pilot outputs.
- Current hypothesis: a manifest generated from final on-disk artifact hashes, canonical split set, code hashes and the pre-freeze Git HEAD will provide an auditable handoff without falsely binding a self-referential commit hash.
- Minimal verification: recompute required artifact SHA256 values, reject any run whose audit trait is not `synthetic`, ensure all expected canonical files exist, run the complete test suite and a clean-tree/remote/CI verification after the freeze commit.
