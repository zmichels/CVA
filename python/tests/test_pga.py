"""MATLAB 6.1.1 reference fixtures; eigenvectors are compared axially."""
import json
from pathlib import Path
import numpy as np
import pytest
import mtex
from cva_mtex import pga

FIXTURES = Path(__file__).resolve().parents[2] / 'tests' / 'fixtures'

@pytest.mark.parametrize('case_id,expected', list(enumerate(json.loads((FIXTURES/'pga_expected.json').read_text())['cases'],1)), ids=lambda item: item['name'] if isinstance(item,dict) else str(item))
def test_matlab_pga_reference(case_id, expected):
    rows = np.loadtxt(FIXTURES/'pga_rotations.csv', delimiter=',')
    matrices = rows[rows[:,0]==case_id,1:].reshape(-1,3,3)
    o = mtex.orientation.byMatrix(matrices, mtex.crystalFrame(expected['crystal_symmetry']))
    result = pga(o)
    np.testing.assert_allclose(result.reference.matrix(), expected['mean'], atol=2e-14)
    np.testing.assert_allclose(result.tensor, expected['tensor'], atol=2e-14)
    np.testing.assert_allclose(result.tangent_vectors, np.array(expected['tangent_vectors']).T, atol=5e-9 if expected['name']=='near_zero' else 2e-14)
    np.testing.assert_allclose(result.eigenvalues, expected['eigenvalues'], atol=2e-14)
    if expected['name'] not in ('near_degenerate','near_zero'):
        np.testing.assert_allclose(abs(np.sum(result.eigenvectors*np.array(expected['eigenvectors']),axis=0)), 1, atol=2e-12)
    assert result.tangent_frame == 'specimen' and result.tensor_units == 'rad^2'

def test_explicit_reference_and_unnormalized_tensor():
    o = mtex.orientation.byMatrix(np.loadtxt(FIXTURES/'pga_rotations.csv',delimiter=',')[:6,1:].reshape(-1,3,3),mtex.crystalFrame('1'))
    ref = o[0]
    result = pga(o, reference=ref)
    t = o.log(ref, side='left').data
    np.testing.assert_allclose(result.tensor, t.T@t)
    assert result.reference is ref

@pytest.mark.parametrize('value', [np.eye(3), [], None])
def test_requires_mtex_orientations(value):
    with pytest.raises(TypeError):
        pga(value)

def test_rejects_empty_and_nonfinite():
    for matrices in (np.empty((0,3,3)), np.full((3,3,3),np.nan)):
        with np.errstate(invalid='ignore'), pytest.raises(ValueError):
            pga(mtex.orientation.byMatrix(matrices,mtex.crystalFrame('1')))
