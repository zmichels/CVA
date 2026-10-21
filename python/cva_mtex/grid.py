"""Mean-referenced 3-by-3 CVA windows on complete square coordinate grids."""
import numpy as np
from .pga import pga
from .results import CVARecord
from ._data import arrays,orientations,lower,finish

def grid_cva(ebsd,*,same_grain=False,reference='mean',min_pixels=3,min_spread=np.deg2rad(.01),bulk=False):
    """Analyze interior centers in MATLAB column-major grid order.

    Retain only indexed pixels of the center phase, and optionally its grain.
    IDs are source IDs; member_ids exposes the exact retained window mask.
    reference="center" is a distinct center-referenced analysis; selection, KOS,
    KAX and tensor all use the selected reference.
    KOS is the maximum disorientation from that reference, in radians; KAX
    is MTEX's mean specimen-frame disorientation axis. Invalid windows vanish.
    """
    if reference not in ('mean','center'):raise ValueError('reference must be mean or center')
    d=arrays(ebsd)
    if same_grain and d.get('grain_ids') is None:raise ValueError('same_grain requires grain IDs')
    xs=np.unique(d['x']);ys=np.unique(d['y'])
    if len(xs)>2 and not np.allclose(np.diff(xs),np.diff(xs)[0]):raise ValueError('regular square grid required')
    if len(ys)>2 and not np.allclose(np.diff(ys),np.diff(ys)[0]):raise ValueError('regular square grid required')
    if len(xs)>1 and len(ys)>1 and not np.isclose(xs[1]-xs[0],ys[1]-ys[0]):raise ValueError('square grid required')
    lookup={(x,y):i for i,(x,y) in enumerate(zip(d['x'],d['y']))}
    if len(lookup)!=len(d['id']):raise ValueError('duplicate grid coordinates')
    records={}
    for ix in range(1,len(xs)-1):
        for iy in range(1,len(ys)-1):
            center=lookup.get((xs[ix],ys[iy]))
            if center is None or not d['_indexed'][center]:continue
            phase=int(d['phase_ids'][center]);indices=[]
            for jy in range(iy-1,iy+2):
                for jx in range(ix-1,ix+2):
                    k=lookup.get((xs[jx],ys[jy]))
                    if k is not None and d['_indexed'][k] and d['phase_ids'][k]==phase and (not same_grain or d['grain_ids'][k]==d['grain_ids'][center]):indices.append(k)
            indices=np.asarray(indices,dtype=int)
            if len(indices)<min_pixels:continue
            o=orientations(d,indices,phase)
            ref=o.mean() if reference=='mean' else orientations(d,[center],phase)[0]
            result=pga(o,reference=ref)
            kos=float(np.max(o.angle(ref)))
            if kos<=min_spread:continue
            axes=np.asarray(o.axis(ref).data).reshape(-1,3).copy()
            # MATLAB quaternion.axis uses crystal x at zero, then orientation.axis
            # rotates it to the specimen frame. Python pairAxis defaults to specimen x.
            zero=np.asarray(o.angle(ref)).reshape(-1)<=2e-10
            if np.any(zero):
                fallback=np.asarray(ref.matrix())[:,0]
                if o.frameA.numSym()==1 and o.frameB.numSym()==1 and len(indices)>1:
                    fallback=-fallback
                axes[zero]=fallback
            kax=np.mean(axes,axis=0)
            if not np.isfinite(kax).all():continue
            id=int(d['id'][center]);records[id]=CVARecord(id,phase,d['id'][indices].copy(),lower(result),kos,kax)
    return finish(records,bulk)
