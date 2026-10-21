"""Optional antipodal de la Vallee Poussin density using MTEX Python 1.0.0.

Fixed-axis fixture agrees with MATLAB at the sampled peak. This bounded
validation does not imply every possible multimodal maximum is identical.
"""
import numpy as np
import mtex

def bulk_axis(axes):
    """Return (lower-hemisphere peak, sampled density) on the 0.25 degree grid.

    Empty samples return (None, None). The estimator uses MTEX's default
    harmonic cutoff with a 10 degree kernel and both axial representatives.
    """
    axes=np.asarray(axes,dtype=float).reshape(-1,3)
    if not len(axes):return None,None
    if not np.isfinite(axes).all() or np.any(np.linalg.norm(axes,axis=1)==0):
        raise ValueError('density axes must be finite and nonzero')
    axes=axes/np.linalg.norm(axes,axis=1)[:,None]
    sample=mtex.vector3d(np.concatenate([axes,-axes]),antipodal=True)
    grid=mtex.plotS2Grid(resolution=.25*mtex.degree,antipodal=True)
    density=np.asarray(mtex.calcDensity(sample,halfwidth=10*mtex.degree,antipodal=True).eval(grid))
    peak=grid.data.reshape(-1,3)[np.argmax(density)].copy()
    if peak[2]>0:peak*=-1
    return peak,density
