import copy
import numpy as np
from cva_mtex import grains_cva,pga
import mtex

def test_three_phases_ids_and_symmetry(multiphase):
    before=copy.deepcopy(multiphase)
    result=grains_cva(multiphase)
    assert list(result.records)==[10,20,30]
    for gid,record in result.records.items():
        mask=multiphase['grain_ids']==gid
        phase=int(multiphase['phase_ids'][mask][0])
        expected=pga(mtex.orientation.byMatrix(multiphase['rotations'][mask],mtex.crystalFrame(multiphase['phases'][phase])))
        assert record.phase_id==phase
        np.testing.assert_allclose(record.pga.tensor,expected.tensor)
        assert (record.pga.eigenvectors[2]<=0).all()
    for key in ('rotations','phase_ids','grain_ids'):
        np.testing.assert_array_equal(multiphase[key],before[key])

def test_selection_minimum_count_and_gos(multiphase):
    grains=dict(id=[30,10,20,99],phase_ids=[3,1,2,1],num_pixel=[2,12,12,3],gos=[1,0,1,1])
    assert list(grains_cva(multiphase,grains).records)==[20]

def test_native_mtex_ebsd_and_grain_objects_are_accepted(multiphase):
    d=multiphase
    e=mtex.EBSD(mtex.vector3d(np.column_stack([d['x'],d['y'],np.zeros(len(d['id']))])),
                mtex.rotation.byMatrix(d['rotations']),d['phase_ids'],
                [mtex.notIndexedFrame()]+[mtex.crystalFrame(d['phases'][p]) for p in (1,2,3)],
                prop={'grainId':d['grain_ids'].copy()},id=d['id'])
    before=e.rotations.data.copy()
    result=grains_cva(e)
    assert list(result.records)==[10,20,30]
    np.testing.assert_array_equal(e.rotations.data,before)
    grains=mtex.grain2d([10,20,30],[1,2,3],e.CSList,[12,12,12],
                        mtex.rotation.byMatrix(np.repeat(np.eye(3)[None],3,axis=0)),[1.,1.,1.],e.frame)
    assert list(grains_cva(e,grains).records)==[10,20,30]

def test_real_three_phase_matlab_fixture():
    import json
    from pathlib import Path
    fixture=json.loads((Path(__file__).resolve().parents[2]/'tests/fixtures/grains_expected.json').read_text())
    cases=fixture['grains'];mats=[];phases=[];gids=[];symmetries={}
    for case in cases:
        n=len(case['rotations']);mats.extend(case['rotations']);phases.extend([case['phase_id']]*n);gids.extend([case['grain_id']]*n)
        symmetries[case['phase_id']]=case['crystal_symmetry']
    n=len(mats)
    result=grains_cva(dict(rotations=mats,phase_ids=phases,grain_ids=gids,id=np.arange(n),x=np.arange(n),y=np.zeros(n),phases=symmetries))
    assert list(result.records)==sorted(c['grain_id'] for c in cases)
    assert len(set(r.phase_id for r in result.records.values()))==3
    for case in cases:
        r=result.records[case['grain_id']].pga
        np.testing.assert_allclose(r.tensor,case['tensor'],atol=1e-12,rtol=1e-10)
        np.testing.assert_allclose(r.eigenvalues,case['eigenvalues'],atol=1e-12,rtol=1e-10)
        np.testing.assert_allclose(abs(np.sum(r.eigenvectors*np.asarray(case['eigenvectors']),axis=0)),1,atol=1e-8)
        expected=mtex.orientation.byMatrix(np.asarray(case['mean']),mtex.crystalFrame(case['crystal_symmetry']))
        assert float(r.reference.angle(expected))<np.deg2rad(.01)
