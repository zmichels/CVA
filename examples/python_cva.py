"""Runnable bounded CVA example using the committed cross-language fixtures.

Run from the repository root: PYTHONPATH=python python examples/python_cva.py
"""
import json
from pathlib import Path

import mtex
import numpy as np

from cva_mtex import pga, grains_cva, grid_cva


FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures"


def main():
    grain_cases = json.loads((FIXTURES / "grains_expected.json").read_text())["grains"]
    first = grain_cases[0]
    orientations = mtex.orientation.byMatrix(
        np.asarray(first["rotations"]), mtex.crystalFrame(first["crystal_symmetry"])
    )
    analysis = pga(orientations)

    matrices, phase_ids, grain_ids, symmetries = [], [], [], {}
    for case in grain_cases:
        count = len(case["rotations"])
        matrices.extend(case["rotations"])
        phase_ids.extend([case["phase_id"]] * count)
        grain_ids.extend([case["grain_id"]] * count)
        symmetries[case["phase_id"]] = case["crystal_symmetry"]
    count = len(matrices)
    grain_input = dict(
        rotations=matrices, phase_ids=phase_ids, grain_ids=grain_ids,
        id=np.arange(count), x=np.arange(count), y=np.zeros(count), phases=symmetries,
    )
    grains = grains_cva(grain_input)

    grid_fixture = json.loads((FIXTURES / "grid_expected.json").read_text())
    grid_input = grid_fixture["input"]
    grid_input["phases"] = {p["id"]: p["symmetry"] for p in grid_input["phases"]}
    windows = grid_cva(grid_input, reference="mean")

    print(f"PGA leading eigenvalue: {analysis.eigenvalues[0]:.8g} rad^2")
    print(f"Grain CVA: {len(grains.records)} grains across {len(symmetries)} phases")
    print(f"Grid CVA: {len(windows.records)} retained center pixels")


if __name__ == "__main__":
    main()
