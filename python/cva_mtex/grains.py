"""Indexed grain CVA with explicit source IDs and phase-specific MTEX frames."""
import numpy as np
from .pga import pga
from .results import CVARecord
from ._data import arrays,orientations,lower,finish

def grains_cva(ebsd,grains=None,*,min_pixels=3,min_gos=np.deg2rad(.01),bulk=False):
    """Return CVA keyed by sorted source grain IDs, without mutating EBSD.

    A supplied MTEX grain2d or mapping (id, phase_ids, num_pixel, gos) supplies
    MATLAB's grain selection metadata. Without it, compute count and mean
    angular spread from indexed pixels. GOS and thresholds are radians.
    """
    d=arrays(ebsd)
    if d.get('grain_ids') is None:raise ValueError('EBSD grain IDs are required')
    selected=None
    if grains is not None:
        if isinstance(grains,dict):
            ids,phase,count,gos=(np.asarray(grains[k]) for k in ('id','phase_ids','num_pixel','gos'))
        else:ids,phase,count,gos=grains.id,grains.phaseId,grains.numPixel,grains.GOS
        selected=set(np.asarray(ids)[(np.asarray(phase)>0)&(np.asarray(count)>=min_pixels)&(np.asarray(gos)>min_gos)].tolist())
    records={}
    for gid in np.unique(d['grain_ids'][d['_indexed']]):
        if selected is not None and gid not in selected:continue
        indices=np.flatnonzero((d['grain_ids']==gid)&d['_indexed'])
        if len(indices)<min_pixels:continue
        phases=np.unique(d['phase_ids'][indices])
        if len(phases)!=1:raise ValueError(f'grain {gid} contains multiple indexed phases')
        phase=int(phases[0]);o=orientations(d,indices,phase);result=pga(o)
        if selected is None and np.mean(o.angle(result.reference))<=min_gos:continue
        records[int(gid)]=CVARecord(int(gid),phase,d['id'][indices].copy(),lower(result))
    return finish(records,bulk)
