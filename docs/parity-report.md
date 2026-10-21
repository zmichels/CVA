# CVA parity report

## Reference and scope

The unchanged root MATLAB functions are the historical reference. The maintained MATLAB functions in `matlab/mtex61/` and the `python/cva_mtex/` package were compared with MATLAB R2025b Update 5, MTEX 6.1.1, Python 3.13.14, and MTEX for Python 1.0.0 on macOS ARM64. This is a bounded numerical validation, not a claim that every EBSD dataset, crystal frame, symmetry boundary, density tie, or later MTEX release produces identical results.

The shared definition is the symmetry-aware MTEX orientation mean, specimen-frame **left** tangent vectors in radians, unnormalized `T = Σ tᵢ tᵢᵀ` in radians squared, descending eigenvalues, and axial eigenvectors. Grain and grid outputs use the lower specimen hemisphere. Grid membership is the indexed center phase within a 3-by-3 window, optionally restricted to the center grain. `reference='mean'` is the default. `reference='center'` is a separate calculation and must be named in analyses.

## Numerical results

| Comparison | Sample | Observed maximum discrepancy |
| --- | --- | --- |
| MATLAB maintained grain loop versus batched | 24 grains, 1,488 pixels, Forsterite/Enstatite/Diopside | Tensor relative Frobenius 1.59×10⁻¹⁶; all-eigenvalue relative 3.83×10⁻¹⁴; all-axis antipodal angle 1.21×10⁻⁶° |
| MATLAB maintained versus root legacy grain | Same 24 grains | Selected grain and phase IDs identical; tensor and all three eigenvalues within 10⁻¹² absolute; sampled KDE bulk vector identical within 10⁻¹² per component |
| MTEX Python PGA versus MATLAB PGA | Rotated no-symmetry and `222` cases, plus near-zero and near-degenerate cases | Tensor absolute Frobenius at most 3.95×10⁻¹⁶; ordinary tangent components at most 3.89×10⁻¹⁵ rad |
| Python versus MATLAB grain CVA | Same 24 three-phase grains | IDs/phases identical; tensor relative 4.36×10⁻¹³; eigenvalue relative 2.43×10⁻¹³; leading axial angle 1.21×10⁻⁶° |
| MATLAB maintained versus root mean grid | 324 kernels in a 20-by-20 forsterite crop | IDs identical; tensor and leading magnitude 2.71×10⁻²⁰ absolute; leading axial angle 1.48×10⁻⁶°; mean 4.83×10⁻⁶°; bulk axial angle 0° |
| Python versus MATLAB mean grid | 15 default and 15 same-grain mixed-phase synthetic kernels | Ordered center/member IDs identical; tensor relative 2.27×10⁻¹²; KOS 2.16×10⁻¹³ rad; KAX component 1.26×10⁻¹¹; leading axial angle ≤1.48×10⁻⁶° |
| Python versus MATLAB center grid | 15 default and 15 same-grain kernels from the same input | IDs/members identical; tensor relative ≤1.47×10⁻¹²; KOS 1.29×10⁻¹³ rad; KAX component 5.74×10⁻¹⁵ |
| Python versus MATLAB antipodal density peak | Four fixed lower-hemisphere axes, 0.25° grid, 10° halfwidth | Same sampled bulk axis within 3.61×10⁻¹⁶ per component; peak-density relative 2.56×10⁻¹⁰ |

Near-zero rotations need a separate numerical floor: MATLAB MTEX 6.1.1 clamps tangent components around 10⁻⁹ rad to zero, while MTEX for Python retains values up to 4.30×10⁻⁹ rad in the fixture. Its tensor difference is only 6.23×10⁻¹⁷ rad². Tests use an explicit 5×10⁻⁹ rad absolute tangent tolerance for that case. Individual axes in the near-degenerate fixture are not compared; the tensor and eigenvalues define the stable result.

For exact zero misorientation, a rotation axis is mathematically undefined. MATLAB uses a reference-rotated crystal-x fallback in its `axis` operation; the Python adapter reproduces that convention when calculating raw arithmetic KAX, including a center pixel compared with itself. The tensor and KOS do not depend on this fallback.

## Checks and reproducibility

The MATLAB tests are in `tests/matlab/`; the Python tests are in `python/tests/`. `tests/matlab/exportCVAFixtures.m` regenerates PGA and three-phase grain references; `tests/fixtures/exportGridFixturesMtex61.m` independently rebuilds the compact mean- and center-grid fixtures from physical XY neighbors and direct per-window MTEX operations; `tests/matlab/exportDensityFixture.m` regenerates the fixed-axis density reference. The grid exporter does not call the maintained `gridCVA`, so its masks and tensors are an independent oracle. The committed fixtures record matrices, symmetry labels, IDs, membership, means/references, tensors, spectra, axes, KOS/KAX, and units where applicable.

From the repository root, with MTEX 6.1.1 available (replace the local executable
and MTEX/Python environment paths for your machine):

```sh
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('/Users/Z/Documents/MATLAB/mtex-6.1.1'); startup_mtex('noMenu'); r=runtests('tests/matlab/testCVAParity.m'); disp(r); assertSuccess(r);"
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('/Users/Z/Documents/MATLAB/mtex-6.1.1'); startup_mtex('noMenu'); r=runtests('tests/matlab/testGrainsCVANoDensity.m'); disp(r); assertSuccess(r);"
/Applications/MATLAB_R2025b.app/bin/matlab -batch "addpath('/Users/Z/Documents/MATLAB/mtex-6.1.1'); startup_mtex('noMenu'); s=matlab.unittest.TestSuite.fromFile('tests/matlab/testGridCVAMtex61.m'); s=s(~contains({s.Name},'Forsterite')); r=run(s); disp(r); assertSuccess(r);"
PYTHONPATH=python '/Users/Z/Documents/EBSD processing python app/.local_artifacts/mtex-python/.venv/bin/python' -m pytest python/tests -q
```

The current MATLAB PGA suite has 2 passing tests with clean process exit; the separate no-density three-phase grain test and no-density grid subset have 1 and 9 passing tests, respectively, with clean exits. The Python suite has 24 passing tests with clean exit, including native and array-mapping cases where MTEX's unindexed phase has ID 1 and six neighboring pixels. Earlier KDE-containing grain and grid comparisons passed their numerical assertions, but those MATLAB `-batch` processes subsequently ended with a native segmentation violation during shutdown (exit 137). The crash recurred with the unchanged root functions and with `-singleCompThread`; it must not be counted as a clean full-suite pass. The bounded MATLAB and Python examples completed with exit 0 while bulk KDE was disabled. See [benchmarks](benchmarks.md) for the timing scope.

The real three-phase grain fixture does not substitute for a broad real-map grid audit. The grid fixture is a controlled 6-by-6 two-phase map with nonindexed and phase-boundary cases. Density parity was checked on one fixed axis set; arbitrary multimodal peaks and exact ties remain outside the evidence. Full lattice/reference-axis conventions should be supplied through native crystal-frame objects when analyzing other symmetries.

## MATLAB MTEX 7.2 compatibility check

MATLAB R2025b Update 5 with MTEX 7.2.0 rejects multiplication of MTEX tangent-vector objects. The maintained `PGA` now converts the left tangent vectors to their 3-by-N Cartesian matrix before forming the same unnormalized outer-product sum. The root legacy `PGA` is unchanged. MTEX 7.2 also warns when a 3-by-3 coordinate matrix is passed to `vector3d`; the maintained grain path now supplies explicit x/y/z component arrays.

The separate `tests/matlab/testCVA_Mtex72.m` run passed four tests and exited cleanly. Its four PGA fixtures matched 6.1.1 tensors within 3.86×10⁻¹⁶ Frobenius norm and eigenvalues within 2.27×10⁻¹⁶ absolute. A synthetic three-phase map yielded three selected grains in both loop and batched modes and 49 valid mean-referenced grid windows. The automated suite also compares 24 fixed grain orientation sets and 60 mean/center grid cases to committed 6.1.1 fixtures, including selected IDs, phase IDs, grid members/masks, tensors, and spectra. The observed grain tensor relative discrepancy was at most 4.36×10⁻¹³ and grid tensor absolute discrepancy at most 6.31×10⁻¹⁶.

SO(3) orientation density and S² direction density both returned finite values under MTEX 7.2.0. The S² check used `v=[xvector,yvector,zvector,xvector]; f=calcDensity(v,'halfwidth',10*degree); eval(f,xvector)` and returned 45.9274630581. An isolated S² `calcDensity` call and even `nfsftmex('precompute',1,1000,1,0)` then caused MATLAB `-batch` to segfault during process shutdown (exit 137) on this macOS installation; the numerical calls completed before the crash. Startup, CTF import, `calcGrains`, and SO(3) density batches exited cleanly. The automated 7.2 suite disables optional CVA bulk-axis density for a reliable exit code. The S² calculation has **not** yet been checked in the interactive MATLAB GUI; that is the next environment check before relying on an interactive bulk-axis workflow. The isolated fast-transform MEX invocation reproduces the shutdown trigger independently of CVA code, so exit 137 is not evidence of a CVA tensor discrepancy.
