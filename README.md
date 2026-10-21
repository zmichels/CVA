# CVA and PGA

Crystallographic Vorticity Axis (CVA) analysis uses Principal Geodesic Analysis (PGA) of EBSD orientations to describe orientation dispersion within grains or local 3-by-3 windows. The leading tensor eigenvector is an **axis**: its sign is interchangeable. Grain and grid outputs place that axis in the lower specimen hemisphere for consistent plotting.

This repository has three paths:

| Path | Use | Tested environment |
| --- | --- | --- |
| [`PGA.m`](PGA.m), [`grainsCVA.m`](grainsCVA.m), [`gridCVA.m`](gridCVA.m) at the root | Unchanged historical MATLAB entry points | MATLAB R2025b, MTEX 6.1.1 bounded checks |
| [`matlab/mtex61/`](matlab/mtex61/) | Maintained MATLAB functions with the same names, phase-aware batching, and explicit options | MATLAB R2025b, MTEX 6.1.1 and 7.2.0; 7.2 automated checks exclude S² bulk density because of a native batch-shutdown crash |
| [`python/`](python/) | Python PGA, grain CVA, and grid CVA using MTEX for Python | Python 3.13, **MTEX for Python 1.0.0** |

The tested Python package is called MTEX for Python 1.0.0; its version number is independent of MATLAB MTEX 7. The separate unpublished `reversePGA.m` in the ODA project is not part of this release.

## Choose the MATLAB path deliberately

MATLAB resolves a function in its current folder before paths added with `addpath`. Run a project from outside this repository root, add exactly the implementation directory you want, and check `which('PGA')`, `which('grainsCVA')`, and `which('gridCVA')` before a long analysis.

```matlab
repo = '/absolute/path/to/CVA';
addpath('/absolute/path/to/mtex-6.1.1');
startup_mtex('noMenu');
addpath(fullfile(repo,'matlab','mtex61'),'-begin');
assert(strcmp(which('PGA'),fullfile(repo,'matlab','mtex61','PGA.m')))

ebsd = mtexdata('forsterite','silent');
[grains,ebsd.grainId] = calcGrains(ebsd('indexed'),'angle',10*degree);
selected = grains(grains.GOS > 0.5*degree);
selected = selected(1:min(24,length(selected)));
[gCVA,~] = grainsCVA(selected,ebsd,'method','batched','computeBulkVector',false);
[eCVA,~] = gridCVA(ebsd,'reference','mean','maxKernels',500,'computeBulkVector',false);
```

The maintained MATLAB functions default to the historical antipodal density-based bulk axis. `computeBulkVector=false` makes exploratory runs cheaper and avoids plotting-grid initialization. `grainsCVA(...,'method','loop')` is a reference for checking the batched result. `gridCVA` batches work **after** calculating each window's MTEX mean, preserving the default mean-referenced analysis. Its explicit `'reference','center'` option changes the analysis and must be reported as such. Other options and returned diagnostics are documented in each function header and the [MATLAB example](examples/matlab_mtex61.m).

For MTEX 7.2.0, initialize that installation and run the bounded compatibility suite from outside the repository root:

```matlab
restoredefaultpath
addpath('/absolute/path/to/mtex-7.2.0')
startup_mtex('noMenu')
r = runtests('/absolute/path/to/CVA/tests/matlab/testCVA_Mtex72.m');
assertSuccess(r)
```

This checks the maintained PGA against the committed 6.1.1 fixtures and compares grain loop/batched results on a synthetic three-phase map, along with mean- and center-referenced grid fixtures. The historical root functions remain the 6.1.1 edition. The [parity report](docs/parity-report.md) records the separate S² density calculation and batch-shutdown limitation.

## Python

Install the local package in an environment with MTEX for Python 1.0.0:

```sh
python -m pip install './python[test]'
PYTHONPATH=python python -m pytest python/tests -q
```

`from cva_mtex import pga, grains_cva, grid_cva` provides the three analyses. `pga` accepts MTEX Python orientations; the adapters accept MTEX Python EBSD and grain objects or documented array mappings. They retain source grain/center IDs, phase IDs, member IDs, tensors, eigenvalues, and axial directions. The grid adapter keeps only neighbors of the center pixel's phase; a phase-boundary window can remain valid with three or more same-phase pixels. See the [Python package guide](python/README.md) and [example](examples/python_cva.py) for inputs and outputs.

## What parity means here

The reference tensor is the **unnormalized** `T = Σ tᵢ tᵢᵀ` in radians squared, where each `tᵢ` is MTEX's specimen-frame left logarithm relative to the symmetry-aware mean orientation. Its eigenvalues therefore depend on the number of included orientations; do not compare magnitudes across unequal grain or window sizes as if they were normalized variances. Eigenvectors are axial, and individual axes of a nearly degenerate eigenspace are not stable.

The maintained MATLAB and Python paths were checked against fixed rotations, crystal-symmetry equivalents, a bounded three-phase forsterite grain fixture, phase-boundary grids, and a fixed antipodal density peak. The [parity report](docs/parity-report.md) gives versions, commands, errors, speed measurements, and limitations. The existing ODA `pyoda.pga` is useful for ODA work but has a different default mean and tangent frame, so its axes are not interchangeable with this CVA parity path.

This is a numerical specimen-frame analysis. A `plottingConvention` changes display, not the physical CVA vectors.

## Reference and license

Z. D. Michels, S. C. Kruckenberg, J. R. Davis, and B. Tikoff (2015), “Determining vorticity axes from grain-scale dispersion of crystallographic orientations,” *Geology*, doi:[10.1130/G36868.1](https://doi.org/10.1130/G36868.1).

The code is provided under the [MIT license](LICENSE). For questions about the original method, contact Zachary Michels at zachary.michels@gmail.com.
