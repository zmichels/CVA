import copy

import mtex
import numpy as np

from cva_mtex import grains_cva, grid_cva, pga


def test_unindexed_phase_one_is_skipped_in_mapping_with_six_neighbors(multiphase):
    data = copy.deepcopy(multiphase)
    # This is the MATLAB fixture convention: phase 1 is notIndexed, and its
    # missing crystal frame must not be mistaken for an indexed phase.
    del data['phases'][1]
    center = np.flatnonzero((data['x'] == 1) & (data['y'] == 1))[0]
    neighborhood = (np.abs(data['x'] - 1) <= 1) & (np.abs(data['y'] - 1) <= 1)
    assert data['phase_ids'][center] == 1
    assert np.count_nonzero(neighborhood & (data['phase_ids'] == 1)) == 6

    grain_result = grains_cva(data)
    grid_result = grid_cva(data)
    assert list(grain_result.records) == [20, 30]
    assert grid_result.records
    assert int(data['id'][center]) not in grid_result.records
    assert all(record.phase_id in (2, 3) for record in grid_result.records.values())


def test_native_mtex_notindexed_frame_is_respected_at_phase_one(multiphase):
    data = multiphase
    position = mtex.vector3d(np.column_stack([data['x'], data['y'], np.zeros(len(data['id']))]))
    frames = [mtex.notIndexedFrame(), mtex.notIndexedFrame(),
              mtex.crystalFrame(data['phases'][2]), mtex.crystalFrame(data['phases'][3])]
    ebsd = mtex.EBSD(position, mtex.rotation.byMatrix(data['rotations']),
                     data['phase_ids'], frames,
                     prop={'grainId': data['grain_ids'].copy()}, id=data['id'])
    assert not np.any(ebsd.isIndexed[data['phase_ids'] == 1])
    grain_result = grains_cva(ebsd)
    grid_result = grid_cva(ebsd)
    assert list(grain_result.records) == [20, 30]
    assert grid_result.records
    assert all(record.phase_id in (2, 3) for record in grid_result.records.values())


def absorbed_pixel_map():
    """A measured grain enclosing a pixel with no acquired orientation."""
    x, y = np.meshgrid(np.arange(5), np.arange(5))
    angles = np.linspace(-0.02, 0.02, 25).reshape(5, 5)
    rotations = np.broadcast_to(np.eye(3), (5, 5, 3, 3)).copy()
    rotations[..., 0, 0] = np.cos(angles)
    rotations[..., 0, 1] = -np.sin(angles)
    rotations[..., 1, 0] = np.sin(angles)
    rotations[..., 1, 1] = np.cos(angles)
    phases = np.ones((5, 5), dtype=int)
    phases[2, 2] = 0
    ebsd = mtex.EBSD(
        mtex.vector3d(np.stack([x, y, np.zeros_like(x)], axis=-1)),
        mtex.rotation.byMatrix(rotations), phases,
        [mtex.notIndexedFrame(), mtex.crystalFrame('1')],
        id=np.arange(101, 126).reshape(5, 5),
    )
    ebsd.rotations.data[2, 2] = np.nan
    grains = mtex.calcGrains(ebsd, angle=np.deg2rad(10))
    return ebsd, grains


def test_absorbed_unmeasured_pixel_is_excluded_from_grain_cva():
    ebsd, grains = absorbed_pixel_map()
    center_id = int(ebsd.id[2, 2])
    assert ebsd.grainId[2, 2] > 0
    assert ebsd.isIndexed[2, 2]
    assert np.isnan(ebsd.rotations.data[2, 2]).all()

    measured = np.isfinite(ebsd.rotations.data).all(axis=-1)
    expected = pga(ebsd[measured].orientations)
    for result in (grains_cva(ebsd, grains), grains_cva(ebsd)):
        assert len(result.records) == 1
        record = next(iter(result.records.values()))
        assert center_id not in record.member_ids
        assert len(record.member_ids) == 24
        np.testing.assert_allclose(record.pga.tensor, expected.tensor)


def test_absorbed_unmeasured_pixel_is_excluded_from_grid_cva():
    ebsd, _ = absorbed_pixel_map()
    center_id = int(ebsd.id[2, 2])
    assert not grid_cva(ebsd, center_ids=[center_id]).records
    for same_grain in (False, True):
        result = grid_cva(ebsd, same_grain=same_grain)
        assert result.records
        assert center_id not in result.records
        assert all(center_id not in record.member_ids for record in result.records.values())


def test_phase_assigned_missing_orientation_is_excluded_in_array_mapping():
    ebsd, _ = absorbed_pixel_map()
    data = dict(
        rotations=ebsd.rotations.matrix(), x=ebsd.pos.x, y=ebsd.pos.y,
        phase_ids=ebsd.phaseId, grain_ids=ebsd.grainId, id=ebsd.id,
        phases={1: '1'}, frame=ebsd.frame,
    )
    center_id = int(ebsd.id[2, 2])
    grain_result = grains_cva(data)
    grid_result = grid_cva(data)
    assert len(grain_result.records) == 1
    assert center_id not in next(iter(grain_result.records.values())).member_ids
    assert grid_result.records
    assert center_id not in grid_result.records
    assert all(center_id not in record.member_ids for record in grid_result.records.values())
