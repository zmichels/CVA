# MATLAB numerical reference fixtures

`pga_rotations.csv` has no header: each row is a numeric case ID followed by
nine rotation-matrix entries in **row-major** order (`m11,m12,...,m33`).
`pga_expected.json` records case names, crystal symmetry, MTEX mean matrices,
unnormalised dispersion tensors, descending eigenvalues, canonicalised axial
eigenvectors (columns), and specimen-frame tangent vectors (3 by N).
The four cases cover a rotated triclinic mean, orthorhombic `222` symmetry,
near-zero dispersion, and a degenerate spectrum. Compare the tensor or the
whole degenerate subspace for that last case; individual axes are arbitrary.

`grains_expected.json` records a bounded selection of 24 grains from all three
indexed phases of MTEX's bundled forsterite dataset. Grain reconstruction uses
10 degree segmentation; eligible grains have at least three pixels and GOS
above 0.5 degrees. Select the first eight eligible grain IDs in each phase.
Each grain stores its phase/ID, point group, its member orientations as an
N by 3 by 3 array, exact MTEX mean, tensor, descending eigenvalues and axial
vectors (columns). Pixel membership is represented by the containing grain.
This fixture tests grain PGA numerics; full EBSD selection is checked live by
`testGrainsCVAMtex61` against the legacy grain routine.

All tangents use MTEX `SO3TangentSpace.leftVector` (specimen coordinates), in
radians. Tensors and eigenvalues have radians-squared units and use a **sum**,
without division by sample count. Crystal symmetry is per case/grain;
specimen symmetry is triclinic. Stored axis signs have no physical meaning.

Regenerate from repository root using MATLAB R2025b and MTEX 6.1.1:

```matlab
addpath('/path/to/mtex-6.1.1'); startup_mtex('noMenu');
addpath('tests/matlab'); exportCVAFixtures;
```

The exporter explicitly calls the untouched repository-root `PGA.m`, regardless
of an active maintained MATLAB path. Tests read the committed files without
regenerating them. To check reproducibility, export into a temporary existing
directory and compare the resulting files byte for byte.

`grid_expected.json` and `grid_center_expected.json` cover every eligible indexed
interior center in a regular 6-by-6, two-phase synthetic map, both with and without
same-grain masking. `exportGridFixturesMtex61` rebuilds the selected member IDs
from physical XY offsets and uses direct MTEX means, logs, and eigensystems;
it does not call the maintained grid function. `density_expected.json` records
one fixed four-axis antipodal density peak. Regenerate these independently from
the repository root with:

```matlab
addpath('/path/to/mtex-6.1.1'); startup_mtex('noMenu');
addpath('tests/fixtures'); exportGridFixturesMtex61;
addpath('tests/matlab'); exportDensityFixture;
```

The density exporter can finish its numerical work but trigger a native MATLAB
shutdown crash in the tested batch environment. Confirm its JSON values and
the process exit status separately.
