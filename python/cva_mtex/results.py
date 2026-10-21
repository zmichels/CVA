"""Scientific results: matrices use columns for descending axial eigenvectors."""
from dataclasses import dataclass
from typing import Any
import numpy as np

@dataclass(frozen=True)
class PGAResult:
    reference: Any
    tangent_vectors: np.ndarray
    tensor: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    tangent_frame: str = 'specimen'
    tangent_units: str = 'rad'
    tensor_units: str = 'rad^2'

@dataclass(frozen=True)
class CVARecord:
    source_id: int
    phase_id: int
    member_ids: np.ndarray
    pga: PGAResult
    kos: float | None = None
    kax: np.ndarray | None = None

@dataclass(frozen=True)
class CVAResult:
    records: dict[int, CVARecord]
    bulk_axis: np.ndarray | None = None
    density: Any = None
