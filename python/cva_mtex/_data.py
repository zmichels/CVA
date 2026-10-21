"""Read-only normalization of MTEX EBSD or explicit per-phase array inputs."""
import numpy as np
import mtex
from dataclasses import replace

def arrays(ebsd):
    if isinstance(ebsd,dict):
        d=dict(ebsd)
        native_indexed=None
    elif isinstance(ebsd,mtex.EBSD):
        d=dict(rotations=ebsd.rotations.matrix(),x=ebsd.pos.x,y=ebsd.pos.y,
               phase_ids=ebsd.phaseId,id=ebsd.id,phases=ebsd.CSList,
               grain_ids=ebsd.prop.get('grainId'),frame=ebsd.frame)
        native_indexed=np.asarray(ebsd.isIndexed & ~ebsd.isPadding).reshape(-1)
    else:
        raise TypeError('ebsd must be MTEX EBSD or an array mapping')
    for key in ('id','phase_ids','x','y'):
        d[key]=np.asarray(d[key]).reshape(-1)
    n=len(d['id'])
    d['rotations']=np.asarray(d['rotations'],dtype=float).reshape(-1,3,3)
    if any(len(d[k])!=n for k in ('rotations','phase_ids','x','y')):
        raise ValueError('EBSD array lengths must agree')
    if len(np.unique(d['id']))!=n:
        raise ValueError('source IDs must be unique')
    if d.get('grain_ids') is not None:
        d['grain_ids']=np.asarray(d['grain_ids']).reshape(-1)
        if len(d['grain_ids'])!=n: raise ValueError('grain IDs length must match')
    if native_indexed is not None:
        d['_indexed']=native_indexed
    else:
        phases=d['phases']
        d['_indexed']=np.zeros(n,dtype=bool)
        for phase in np.unique(d['phase_ids']):
            if phase<=0:continue
            frame=_phase_frame(phases,phase)
            if isinstance(frame,str):
                is_indexed=frame.strip().lower()!='notindexed'
            else:
                is_indexed=bool(getattr(frame,'isIndexed',False))
            if is_indexed:d['_indexed'][d['phase_ids']==phase]=True
    return d

def _phase_frame(phases,phase):
    phase=int(phase)
    if isinstance(phases,dict):return phases.get(phase,phases.get(str(phase)))
    try:return phases[phase]
    except (IndexError,KeyError):return None

def orientations(d,indices,phase):
    cs=_phase_frame(d['phases'],phase)
    if cs is None: raise ValueError(f'crystal symmetry missing for phase {phase}')
    if isinstance(cs,str): cs=mtex.crystalFrame(cs)
    frames=[cs]
    if d.get('frame') is not None: frames.append(d['frame'])
    return mtex.orientation.byMatrix(d['rotations'][indices],*frames)

def lower(result):
    axes=result.eigenvectors.copy()
    axes[:,axes[2]>0]*=-1
    return replace(result,eigenvectors=axes)

def finish(records,bulk):
    from .results import CVAResult
    if not bulk:return CVAResult(records)
    from .density import bulk_axis
    axis,density=bulk_axis(np.array([r.pga.eigenvectors[:,0] for r in records.values()]))
    return CVAResult(records,axis,density)
