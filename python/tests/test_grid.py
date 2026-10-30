import numpy as np
import pytest
from cva_mtex import grid_cva

def test_phase_boundary_window_survives_and_masks_exactly(multiphase):
    result=grid_cva(multiphase)
    assert list(result.records)==[108,109,110,111,114,115,116,117,120,121,122,123,126,127,128,129]
    for id,record in result.records.items():
        center=np.flatnonzero(multiphase['id']==id)[0]
        members=np.isin(multiphase['id'],record.member_ids)
        assert np.all(multiphase['phase_ids'][members]==multiphase['phase_ids'][center])
        assert len(record.member_ids)==6
        assert record.kos>np.deg2rad(.01)
        assert np.isfinite(record.kax).all()

def test_same_grain_mask_and_invalid_centers(multiphase):
    multiphase['grain_ids']=np.arange(36)
    assert len(grid_cva(multiphase,same_grain=True).records)==0
    multiphase['phase_ids'][7]=0
    assert 108 not in grid_cva(multiphase).records

def test_selected_center_ids_match_full_grid_in_spatial_order(multiphase):
    full=grid_cva(multiphase)
    selected=grid_cva(multiphase,center_ids=[129,108,117,129])
    assert list(selected.records)==[108,117,129]
    for source_id,record in selected.records.items():
        original=full.records[source_id]
        np.testing.assert_array_equal(record.member_ids,original.member_ids)
        np.testing.assert_allclose(record.pga.tensor,original.pga.tensor)
    assert not grid_cva(multiphase,center_ids=[]).records
    with pytest.raises(ValueError,match='source ID not present'):
        grid_cva(multiphase,center_ids=[999])

def test_selected_grid_covaries_with_specimen_rotation_keepxy(multiphase):
    import mtex
    d=multiphase
    ebsd=mtex.EBSD(mtex.vector3d(np.column_stack([d['x'],d['y'],np.zeros(len(d['id']))])),
                   mtex.rotation.byMatrix(d['rotations']),d['phase_ids'],
                   [mtex.notIndexedFrame()]+[mtex.crystalFrame(d['phases'][p]) for p in (1,2,3)],
                   prop={'grainId':d['grain_ids']},id=d['id'])
    rotation=mtex.rotation.byAxisAngle(mtex.vector3d([1,2,3]),37*mtex.degree)
    moved=ebsd.rotate(rotation,keepXY=True)
    np.testing.assert_array_equal(moved.pos.data,ebsd.pos.data)
    before=grid_cva(ebsd,center_ids=[108,117,129])
    after=grid_cva(moved,center_ids=[108,117,129])
    assert before.records.keys()==after.records.keys()
    matrix=np.asarray(rotation.matrix())
    for source_id,record in before.records.items():
        changed=after.records[source_id]
        np.testing.assert_array_equal(changed.member_ids,record.member_ids)
        np.testing.assert_allclose(changed.pga.tensor,matrix@record.pga.tensor@matrix.T,atol=1e-12)
        np.testing.assert_allclose(changed.pga.eigenvalues,record.pga.eigenvalues,atol=1e-12)

def test_low_spread_is_excluded(multiphase):
    multiphase['rotations'][:]=multiphase['rotations'][0]
    assert len(grid_cva(multiphase).records)==0

@pytest.mark.parametrize('reference', ['mean','center'])
def test_matlab_grid_fixture(reference):
    import json
    from pathlib import Path
    filename='grid_expected.json' if reference=='mean' else 'grid_center_expected.json'
    fixture=json.loads((Path(__file__).resolve().parents[2]/'tests/fixtures'/filename).read_text())
    d=fixture['input'];d['phases']={p['id']:p['symmetry'] for p in d['phases']}
    for mode in ('default','same_grain'):
        result=grid_cva(d,same_grain=mode=='same_grain',reference=reference)
        assert list(result.records)==[r['id'] for r in fixture[mode]]
        for expected in fixture[mode]:
            r=result.records[expected['id']]
            np.testing.assert_array_equal(r.member_ids,expected['member_ids'])
            np.testing.assert_allclose(r.pga.reference.matrix(),expected['mean'],atol=1e-12)
            np.testing.assert_allclose(abs(np.sum(r.pga.eigenvectors*np.asarray(expected['eigenvectors']),axis=0)),1,atol=1e-8)
            np.testing.assert_allclose(r.pga.tensor,expected['tensor'],atol=1e-12,rtol=1e-9)
            np.testing.assert_allclose(r.pga.eigenvalues,expected['eigenvalues'],atol=1e-12,rtol=1e-9)
            np.testing.assert_allclose(r.kos,expected['kos'],atol=1e-12)
            np.testing.assert_allclose(r.kax,expected['kax'],atol=1e-9)

def test_explicit_center_reference_is_distinct(multiphase):
    mean=grid_cva(multiphase,reference='mean')
    center=grid_cva(multiphase,reference='center')
    common=mean.records.keys() & center.records.keys()
    assert common
    assert any(np.linalg.norm(mean.records[i].pga.tensor-center.records[i].pga.tensor)>1e-8 for i in common)
    import pytest
    with pytest.raises(ValueError,match='reference'):
        grid_cva(multiphase,reference='intrinsic')

def test_kax_keeps_matlab_arithmetic_mean_magnitude(multiphase):
    import mtex
    result=grid_cva(multiphase)
    r=next(iter(result.records.values()))
    mask=np.isin(multiphase['id'],r.member_ids)
    o=mtex.orientation.byMatrix(multiphase['rotations'][mask],mtex.crystalFrame(multiphase['phases'][r.phase_id]))
    np.testing.assert_allclose(r.kax,np.mean(o.axis(o.mean()).data,axis=0),atol=1e-12)
    assert np.linalg.norm(r.kax)<1

def test_center_reference_controls_selection_kos_and_kax(multiphase):
    import mtex
    d=multiphase
    result=grid_cva(d,reference='center')
    assert result.records
    for source_id,r in result.records.items():
        indices=np.flatnonzero(np.isin(d['id'],r.member_ids))
        o=mtex.orientation.byMatrix(d['rotations'][indices],mtex.crystalFrame(d['phases'][r.phase_id]))
        center=np.flatnonzero(d['id']==source_id)[0]
        ref=mtex.orientation.byMatrix(d['rotations'][center],mtex.crystalFrame(d['phases'][r.phase_id]))
        np.testing.assert_allclose(r.kos,np.max(o.angle(ref)),atol=1e-12)
        axes=o.axis(ref).data.copy()
        zero=o.angle(ref)<=2e-10
        fallback=ref.matrix()[:,0]
        if o.frameA.numSym()==1 and o.frameB.numSym()==1: fallback=-fallback
        axes[zero]=fallback
        np.testing.assert_allclose(r.kax,np.mean(axes,axis=0),atol=1e-12)
    mean=grid_cva(d,reference='mean')
    id=next(i for i in result.records if result.records[i].kos>mean.records[i].kos+1e-6)
    threshold=(result.records[id].kos+mean.records[id].kos)/2
    assert id in grid_cva(d,reference='center',min_spread=threshold).records
    assert id not in grid_cva(d,reference='mean',min_spread=threshold).records
