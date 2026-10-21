# CVA with MTEX for Python

Requires Python >=3.12. Tested with Python 3.13.14 and **MTEX for Python 1.0.0**. This package does not
claim MTEX 7 compatibility. Install from this directory with `python -m pip
install '.[test]'`, or from its parent using `python -m pip install './python[test]'`.

```python
import mtex
from cva_mtex import pga, grains_cva, grid_cva

# matrices: (N, 3, 3), crystal-to-specimen orientation matrices
orientations = mtex.orientation.byMatrix(matrices, mtex.crystalFrame('mmm'))
result = pga(orientations)
# MTEX EBSD and grain2d objects:
grain_result = grains_cva(ebsd, grains)
grid_result = grid_cva(ebsd, reference='mean', same_grain=True)
```

`pga` uses the MTEX quaternion mean and `log(reference, side='left')` with an
**unnormalized** tensor in radians squared. `reference` is optional for PGA.
`PGAResult.eigenvectors` stores descending axial eigenvectors in columns.
PGA preserves the eigensolver's arbitrary axial sign; adapters project each
axis to z <= 0. Display conventions do not rotate coordinates.

Adapters also accept a dictionary with `rotations` (N×3×3), `phase_ids`,
`grain_ids`, `id`, `x`, `y`, and `phases` mapping integer phase IDs to MTEX
crystal frames or symmetry strings. Phase IDs <= 0 are unindexed; a positive
phase ID is indexed only when its `phases` entry names a crystal frame or
symmetry. A missing entry, `None`, or a `notIndexedFrame` marks that phase
unindexed, including phase ID 1 in MATLAB-exported arrays. For native MTEX
EBSD, the adapter uses `ebsd.isIndexed` and excludes padding. An indexed
crystal phase may still have ID 1 when its frame is present.
For grain selection a dictionary may supply `id`, `phase_ids`, `num_pixel`,
and `gos`; GOS is in radians. With no grain metadata the adapter derives
count and GOS from indexed samples. MTEX EBSD inputs require precomputed
`prop['grainId']` for grains or `same_grain=True`; the package never invokes
`calcGrains` or mutates the input. Use full native crystal frames for lattice
and reference-axis conventions, especially for monoclinic/triclinic data.

`CVAResult.records` is an ordered mapping from source grain/center IDs to
`CVARecord`, which includes phase ID, exact retained member IDs and PGA result.
Grains are ordered by grain ID. Grid centers use increasing x then increasing y
(MATLAB column-major order); member IDs follow MATLAB's 3×3 row sequence.
Grid inputs must use a regular square coordinate grid; missing coordinates are
omitted from windows. Only interior, indexed centers with at least three same
phase pixels and maximum reference disorientation > 0.01 degrees are retained.
`same_grain=True` additionally filters each window by its center grain ID.
`reference='center'` is a distinct scientific analysis; it does not substitute
for the default mean analysis. Selection, KOS, KAX and the tensor all use the
selected reference, matching maintained MATLAB center mode.
KOS is in radians. KAX preserves MATLAB's arithmetic mean of specimen axes,
including its magnitude and NaN propagation; it is not normalized. A zero
misorientation has no physical axis. To match MATLAB, its conventional axis
contribution is reference-rotated crystal x; MATLAB's multi-orientation no-symmetry
branch uses the negative of this vector. The self member stays in center-mode
KAX; it is not silently dropped.

`bulk=True` additionally evaluates MTEX's antipodal de la Vallee Poussin
density, 10 degree halfwidth, on a 0.25 degree grid. The bulk axis is the
lower-hemisphere sampled maximum; `density` contains the sampled grid array.
A fixed four-axis fixture matches MATLAB's peak axis and peak value; the full
density field and arbitrary multimodal peak/tie parity has not been
exhaustively established. Empty samples return no
bulk axis/density. Density computation is optional because the fine grid adds cost.

The bounded parity suite uses committed MATLAB fixtures: four PGA cases,
24 real forsterite grains across three phases, 15 kernels for each combination of mean/center reference and default/same-grain
masking, and a fixed density peak. For the near-zero fixture MATLAB clamps
~1e-9 rad logarithms to zero while Python retains them (maximum 4.30e-9 rad,
tensor difference 6.23e-17 rad²); tests apply an explicit absolute floor.
Near-degenerate eigensystems are checked through tensors and eigenvalues.
Individual axes in degenerate eigenspaces are not stable scientific outputs.

Run from the repository root:

```sh
PYTHONPATH=python python -m pytest python/tests -q
```
