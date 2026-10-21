import copy

import mtex
import numpy as np

from cva_mtex import grains_cva, grid_cva


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
