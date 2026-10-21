# CVA MATLAB and Python parity release design

## Purpose and scope

Release a maintained CVA/PGA path for MATLAB with MTEX 6.1.1 and for MTEX for Python 1.0.0 while retaining the public MATLAB functions as a usable legacy reference. The initial release covers `PGA`, `grainsCVA`, and `gridCVA`. The unpublished `reversePGA.m` in the separate ODA checkout is outside this release.

The repository root's `PGA.m`, `grainsCVA.m`, and `gridCVA.m` remain byte-for-byte unchanged so existing MATLAB paths and calls continue to resolve. The new MATLAB implementation lives in `matlab/mtex61/`; users select that directory explicitly on their MATLAB path. The Python implementation lives in `python/` as an installable package. A compatibility table in the README states which path and MTEX version each implementation was tested with. The README currently claims MIT licensing but the repository has no license file; add the matching MIT license text as release hygiene.

## Scientific contract

MATLAB MTEX 6.1.1 is the numerical reference. `PGA` uses MTEX's symmetry-projected quaternion mean when no reference is supplied, MTEX's specimen-frame `SO3TangentSpace.leftVector` logarithms relative to that mean, the **unnormalized** tensor `T = sum(t_i * t_i')` in radians squared, and a descending symmetric eigensystem. Eigenvectors are axial; comparisons use the smaller of the angles to `v` and `-v`. A near-degenerate eigenvalue pair does not define stable individual axes, so validation compares its invariant subspace or tensor instead.

Grain CVA retains indexed-pixel selection, grain membership, minimum three pixels, GOS threshold of 0.01 degrees, per-phase crystal symmetry, tensor/eigenvalue conventions, lower-hemisphere sign, and the existing optional antipodal density-based bulk vector. Grid CVA retains the 3-by-3 neighborhood and source ordering; for an indexed center it keeps only indexed pixels of the **center phase** (and, when requested, center grain), requires at least three retained pixels and more than 0.01 degrees orientation spread, and uses the MTEX mean of those pixels by default. A boundary-crossing window can remain valid after its other-phase pixels are removed. KOS, KAX, mean orientation, tensors, eigenaxes, and density output retain their documented meanings.

Display conventions do not rotate the stored CVA vectors. Axial signs and plot presentation must be distinguished from specimen-frame coordinates.

## Implementations and speed

The MATLAB `matlab/mtex61/` functions keep the public signatures and output-property names. Grain processing groups selected pixels by phase, obtains symmetry-correct reference orientations, evaluates left tangent vectors in batches, and accumulates six independent tensor components by grain. Reusing grain-object means is acceptable only if recorded parity tolerances hold across synthetic, single-phase, and mixed-phase checks; retain the loop implementation as the reference and offer an exact-mean route if needed. Existing ODA research code is a candidate, not the release source of truth.

For mean-referenced grid CVA, first retain MTEX's `mean(o)` for each valid window and batch the subsequent tangent, tensor, and eigen calculations in phase-specific bounded chunks. Benchmark that against the legacy loop. A later pagewise quaternion-mean implementation is warranted only if the mean calls dominate and it reproduces MTEX's symmetry projection and conditional second pass for every tested window. The center-reference grid mode may be exposed as a separately named alternative; it is **not** a speed-equivalent version of mean-referenced CVA. On one forsterite probe its leading axis differed by a median 14.0 degrees across 490 common kernels.

The Python package uses MTEX for Python 1.0.0 `orientation.mean()` and `orientation.log(ref, side="left")` for the parity path, with NumPy reduction and eigendecomposition. It provides `pga`, `grains_cva`, and `grid_cva` interfaces accepting MTEX Python orientation/EBSD objects, plus documented results with units, references, phase IDs, masks, and axial directions. It does not silently use the current ODA `pyoda.pga`: that implementation uses an intrinsic mean and body-frame tangent vectors, which gave a 42.19-degree principal-axis difference on a controlled nonidentity fixture. An independent SciPy kernel is a future optimization only after it matches the MTEX path.

## Validation and release evidence

Shared, bounded MATLAB-generated fixtures record orientation matrices, phase/grain IDs, grid positions, crystal/specimen symmetry, mean/reference orientation, tangent frame, tensor, eigenvalues, axial eigenvectors, KOS/KAX, and retained output IDs. Include no-symmetry and nontrivial symmetry cases, rotated means, phase-boundary windows, nonindexed pixels, same-grain selection, low spread, and near-degenerate tensors. Compare MATLAB legacy to enhanced MATLAB and both to Python; compare means, tensors, spectra, axes/subspaces, grain/kernel selection, and density/bulk-vector behavior separately. Any intentionally changed numerical or selection behavior is reported with its trigger and measured effect rather than hidden under a parity claim.

Initial acceptance targets for well-separated eigensystems: identical selected IDs and phase assignments; mean and principal-axis differences no greater than 0.01 degrees; relative eigenvalue/tensor differences no greater than 1e-4 for enhanced MATLAB versus legacy; tighter controlled Python/MATLAB PGA tensor agreement on the order of 1e-10 relative. Set absolute floors for nearly zero tensors and compare degenerate eigenspaces rather than individual eigenvectors. The final test report records actual maxima, fixture scope, MTEX/MATLAB/Python versions, and any exceptions. Timing uses identical bounded datasets and options, reports sample size and repeated wall-time distributions, and never substitutes the center-reference grid mode for a parity benchmark.

Read-only pre-release evidence: the ODA MATLAB CVA suite passed 22 tests on MATLAB R2025b with MTEX 6.1.1. A three-phase forsterite probe returned the same 109 grains for its loop and batched methods, with maximum axial leading-axis difference 0.00080 degrees and maximum relative leading-eigenvalue difference 5.18e-6. Controlled Python MTEX PGA tensors agreed with MATLAB within about 5e-16 absolute Frobenius norm. Those probes are not a release pass; the public root grain/grid functions completed a bounded smoke test before MATLAB crashed during batch-process shutdown, so that run did not exit cleanly.

## Delivery boundary

Deliver source, fixtures, repeatable checks, examples, versioned documentation, and a reviewable draft PR from an isolated `codex/` branch. Preserve the dirty ODA CVA checkout untouched. Do not claim MTEX 7 compatibility from the tested Python package's `1.0.0` version. Keep `reversePGA` and a dedicated UI/workbench out of this release.
