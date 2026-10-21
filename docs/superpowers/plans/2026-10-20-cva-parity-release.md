# CVA MATLAB and Python Parity Release Implementation Plan

Status: delivered in the draft PR. Numerical outcomes and the MATLAB batch-shutdown limitation are recorded in `docs/parity-report.md`.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkboxes for tracking.

**Goal:** Deliver legacy-preserving MATLAB MTEX 6.1.1 and MTEX for Python 1.0.0 CVA/PGA implementations with measured cross-language parity and a reviewable draft PR.

**Architecture:** The unchanged root MATLAB files are the legacy reference. New MATLAB functions in `matlab/mtex61/` retain the same names and property outputs, and a Python package in `python/cva_mtex/` uses MTEX Python means and specimen-frame left logarithms with NumPy tensor reduction. Shared MATLAB-generated fixtures and bounded mixed-phase EBSD checks define parity; speed work follows parity.

**Tech Stack:** MATLAB R2025b, MTEX 6.1.1, Python 3.13, MTEX for Python 1.0.0, NumPy, pytest.

## Global Constraints

- Do not edit root `PGA.m`, `grainsCVA.m`, or `gridCVA.m`, or the separate dirty ODA CVA checkout.
- Default grid CVA uses the mean of each center-phase-masked 3-by-3 window; center reference is a separately named method.
- Preserve `T = sum(t_i t_i')` in radians squared, MTEX symmetry handling, specimen-frame left tangents, descending eigenvalues, and axial direction comparisons.
- Preserve per-phase symmetry, grain/window membership, filters, output properties, and the 10-degree-halfwidth antipodal bulk-vector calculation unless a measured exception is documented.
- No MTEX 7 compatibility claim: the Python version tested here is MTEX for Python 1.0.0.
- Write behavior tests before each implementation, see them fail for the missing behavior, and rerun after implementation.

---

### Task 1: MATLAB parity fixtures and reference checks

**Files:**
- Create: `tests/matlab/testCVAParity.m`
- Create: `tests/matlab/exportCVAFixtures.m`
- Create: `tests/fixtures/README.md`
- Create: `tests/fixtures/pga_rotations.csv`, `tests/fixtures/pga_expected.json`

**Interfaces:** Produces fixed orientation matrices, expected mean/tensor/eigensystem, and a MATLAB check callable from `matlab -batch`; later tasks use the same fixtures.

- [x] Write MATLAB tests that compare the root `PGA` with a proposed `matlab/mtex61/PGA` on nonidentity, `222`-symmetry, near-zero, and near-degenerate orientation sets; compare tensors and axial/subspace results, not quaternion signs. Confirm the new-path test fails because that function is absent.
- [x] Export matrices, crystal symmetry, reference mean, tensor, eigenvalues, and canonicalized axial eigenvectors to JSON/CSV from root MATLAB. Include radians and coordinate-frame metadata. Verify regenerated fixture values against the committed files.
- [x] Run `/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('/Users/Z/Documents/MATLAB/mtex-6.1.1'); startup_mtex('noMenu'); addpath(fullfile(pwd,'matlab','mtex61'),'-begin'); r=runtests('tests/matlab/testCVAParity.m'); assertSuccess(r);"` from the CVA worktree once the new functions exist; record the exit code and failures. The test must switch paths explicitly when it calls root legacy functions.
- [x] Commit fixtures and tests with a message describing the reference contract.

### Task 2: MATLAB PGA and per-grain CVA

**Files:**
- Create: `matlab/mtex61/PGA.m`, `matlab/mtex61/grainsCVA.m`
- Create: `tests/matlab/testGrainsCVAMtex61.m`

**Interfaces:** `PGA(ori, oriRef)` returns `[eV,mags,T]`; `grainsCVA(grains,ebsd,Name,Value)` returns legacy-compatible `[gCVA,bv]` and an optional diagnostics struct. `method='batched'` is default; `method='loop'` is the reference. `computeBulkVector=false`, `storeTensor=false`, and `showProgress=false` may be used to isolate numerical tests and benchmarks.

- [x] Write a failing test for `PGA(ori,oriRef)`, full tensor agreement with the root `PGA(ori)` when `oriRef=mean(ori)`, and descending eigenvalues.
- [x] Implement the MTEX 6.1.1-specific left-vector PGA without legacy version parsing; rerun the PGA tests.
- [x] Write failing tests for loop/batched agreement on selected grain IDs, phase IDs, all three magnitudes, axial axes, and optional ODT; include the MTEX forsterite three-phase data and a bounded selection of grains with GOS above 0.5 degrees.
- [x] Implement the grain path with phase-specific orientation construction, batched left-log evaluation, six-component `accumarray` tensor reduction, and the exact legacy selection and density defaults. Reuse `grains.meanOrientation` only if parity limits hold; otherwise calculate the legacy mean per grain and batch the remaining work.
- [x] Run the complete MATLAB grain tests and a repeated same-input loop/batched timing; record number of grains, maximum relative tensor/eigenvalue errors, maximum axial error, and timings. Commit the implementation and tests.

### Task 3: MATLAB mean-referenced grid CVA with bounded batching

**Files:**
- Create: `matlab/mtex61/gridCVA.m`
- Create: `tests/matlab/testGridCVAMtex61.m`

**Interfaces:** `gridCVA(ebsd, [grains], Name,Value)` returns legacy-compatible `[eCVA,bv]` and optional diagnostics. `reference='mean'` is default; `reference='center'` is a named alternative. Optional `sameGrain`, `maxKernels`, `computeBulkVector`, `storeTensor`, and `showProgress` allow bounded checks. Do not alter root `gridCVA`.

- [x] Write failing tests comparing root and new mean-reference outputs for a 20-by-20 forsterite crop: retained IDs, center phase, KOS/KAX, means, tensors, magnitudes, and axial axes. Add mixed-phase and phase-boundary cases: remove other-phase neighbors while retaining valid center-phase windows.
- [x] Implement exact legacy window selection/order and per-window MTEX `mean(o)`; batch post-mean left logs, six tensor components, and eigensystems in bounded phase-specific chunks. Preserve optional boundary exclusion, lower-hemisphere convention, KDE, and output properties. Remove global `warning off` from the new path.
- [x] Run the MATLAB grid tests and compare repeated timings to the root grid function on identical bounded input. If batching does not improve mean-reference timing materially, report that honestly and keep correctness; pagewise mean optimization requires a separate passing test for MTEX symmetry projection and conditional second mean.
- [x] Commit the grid implementation and tests.

### Task 4: Python PGA parity core

**Files:**
- Create: `python/pyproject.toml`, `python/cva_mtex/__init__.py`, `python/cva_mtex/pga.py`, `python/cva_mtex/results.py`
- Create: `python/tests/test_pga.py`

**Interfaces:** `pga(orientations, reference=None) -> PGAResult`, where `PGAResult` records mean/reference, tangent vectors, 3-by-3 unnormalized tensor, descending eigenvalues, and 3-by-3 specimen-frame axial eigenvectors. MTEX Python orientation inputs are required for the parity route.

- [x] Write pytest cases loading both committed PGA fixtures through `mtex.orientation.byMatrix` and asserting mean, tensor, eigenvalue, and axial results; see them fail because `cva_mtex` is absent.
- [x] Implement with `orientations.mean()`, `orientations.log(reference, side='left').data`, `T=t.T @ t`, and `numpy.linalg.eigh`; validate input sizes, finite values, and radians metadata. Do not call ODA `pyoda.pga` as a substitute.
- [x] Run `PYTHONPATH=python '/Users/Z/Documents/EBSD processing python app/.local_artifacts/mtex-python/.venv/bin/python' -m pytest python/tests/test_pga.py -q`, then the whole Python suite. Commit core and tests.

### Task 5: Python grain and grid adapters

**Files:**
- Create: `python/cva_mtex/grains.py`, `python/cva_mtex/grid.py`, `python/cva_mtex/density.py`
- Create: `python/tests/test_grains.py`, `python/tests/test_grid.py`
- Create: `tests/fixtures/grid_expected.json`, `tests/fixtures/grid_center_expected.json`, and the three-phase grain fixture `tests/fixtures/grains_expected.json`

**Interfaces:** `grains_cva(ebsd, grains, ...)` and `grid_cva(ebsd, ...)` return result objects keyed by source grain or center-pixel IDs, with phase IDs, tensor/eigensystem, means, KOS/KAX where applicable, and optional bulk axis. The Python adapters do not silently mutate caller EBSD data.

- [x] Write failing Python tests for three-phase grain grouping, center-phase 3-by-3 masks, optional same-grain masking, minimum count and spread, retained IDs, and root-MATLAB fixture values. Include an explicit test that a phase-boundary window can remain valid after other-phase neighbors are removed.
- [x] Implement arrays and MTEX Python orientation objects grouped by phase, using each phase's crystal symmetry. Match MATLAB's grain selection and grid order; call the PGA core for the strict mean-reference mode. Add bounded chunking only where it preserves outputs.
- [x] Implement the optional antipodal density/bulk axis using MTEX Python and validate against MATLAB on a fixed axis set. If a density API difference prevents parity, surface a documented limitation instead of silently changing the estimator.
- [x] Run all Python tests and a bounded real three-phase comparison against MATLAB fixtures; commit adapters and tests.

### Task 6: Documentation, license, benchmark, and PR

**Files:**
- Modify: `README.md`
- Create: `LICENSE`, `examples/matlab_mtex61.m`, `examples/python_cva.py`, `docs/parity-report.md`, `docs/benchmarks.md`

**Interfaces:** Users can select root legacy MATLAB, enhanced MATLAB, or Python without shadowing ambiguity; report gives source versions, commands, input counts, tolerances, measured differences, and speed distributions.

- [x] Write the compatibility table, installation/path guidance, full multi-phase examples, method/units explanation, and warnings that center reference and ODA pyoda are not direct parity modes. Add MIT license text to match the existing README claim.
- [x] Run the full MATLAB and Python test commands, MATLAB and Python examples, `git diff --check`, package metadata validation, and a fresh same-input benchmark. Inspect outputs and record all passing counts and limitations in the parity report.
- [x] Review the diff for untouched root files and untouched ODA checkout; have a separate reviewer check scientific contract and API/doc clarity. Commit docs and reports.
- [x] Push the `codex/cva-matlab-python-parity` branch, create a draft PR against the public CVA repository, and attach the PR to this chat.
