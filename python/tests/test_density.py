import json
from pathlib import Path
import numpy as np
from cva_mtex.density import bulk_axis

def test_fixed_matlab_antipodal_density_peak():
    fixture=json.loads((Path(__file__).resolve().parents[2]/'tests/fixtures/density_expected.json').read_text())
    axis,density=bulk_axis(np.array(fixture['axis_xyz']))
    np.testing.assert_allclose(axis,fixture['peak_xyz'],atol=1e-12)
    np.testing.assert_allclose(np.max(density),fixture['peak_density'],rtol=1e-8)
