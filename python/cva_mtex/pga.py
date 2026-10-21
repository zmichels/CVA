"""MTEX mean and left logarithm; no intrinsic/body-frame substitution."""
import numpy as np
from mtex import orientation
from .results import PGAResult

def pga(orientations, reference=None):
    """Return the unnormalized specimen-frame PGA tensor in radians squared.

    ``reference=None`` uses the MTEX symmetry-projected quaternion mean.
    Eigenvectors are axial and stored as columns in descending spectral order.
    """
    if not isinstance(orientations, orientation):
        raise TypeError('orientations must be an MTEX orientation object')
    if np.asarray(orientations.matrix()).size == 0:
        raise ValueError('at least one orientation is required')
    if not np.isfinite(orientations.matrix()).all():
        raise ValueError('orientations must be finite')
    if reference is not None:
        if not isinstance(reference, orientation):
            raise TypeError('reference must be an MTEX orientation object')
        if np.asarray(reference.matrix()).shape != (3,3) or not np.isfinite(reference.matrix()).all():
            raise ValueError('reference must be one finite orientation')
    else:
        reference = orientations.mean()
    tangent = np.asarray(orientations.log(reference, side='left').data, dtype=float).reshape(-1,3)
    if not np.isfinite(tangent).all():
        raise ValueError('MTEX returned nonfinite tangent vectors')
    tensor = tangent.T @ tangent
    values, vectors = np.linalg.eigh(tensor)
    return PGAResult(reference, tangent, tensor, values[::-1], vectors[:,::-1])
