from __future__ import annotations

import numpy as np
from pymatgen.core import Structure

from .frames import LocalFrame


def _unit(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm < 1e-12:
        raise ValueError("Cannot normalize a zero-length vector")
    return vector / norm


def parse_vector(text: str) -> np.ndarray:
    cleaned = text.strip().replace("(", "").replace(")", "")
    tokens = cleaned.replace(",", " ").split()
    if len(tokens) != 3:
        raise ValueError(f"Vector must have three components, got: {text!r}")
    try:
        return np.array([float(token) for token in tokens], dtype=float)
    except ValueError as exc:
        raise ValueError(f"Vector must contain numbers, got: {text!r}") from exc


def build_vector_frame(
    structure: Structure,
    center_index: int,
    x_vector: np.ndarray,
    z_vector: np.ndarray,
) -> LocalFrame:
    """Build a local frame from explicit Cartesian x and z directions.

    The supplied z direction is normalized first. x is projected perpendicular
    to z, y=z×x, and x is recomputed as y×z to enforce an exactly orthonormal,
    right-handed frame.
    """
    if not 0 <= center_index < len(structure):
        raise ValueError(f"center_index {center_index} is outside the structure")

    z = _unit(z_vector)
    x_seed = np.asarray(x_vector, dtype=float) - float(np.dot(x_vector, z)) * z
    x = _unit(x_seed)
    y = _unit(np.cross(z, x))
    x = _unit(np.cross(y, z))

    R = np.column_stack((x, y, z))
    if not np.allclose(R.T @ R, np.eye(3), atol=1e-10):
        raise ValueError("Constructed vector frame is not orthonormal")
    if not np.isclose(np.linalg.det(R), 1.0, atol=1e-10):
        raise ValueError("Constructed vector frame is not a proper right-handed rotation")

    return LocalFrame(
        site_index=center_index,
        center_symbol=structure[center_index].specie.symbol,
        ligand_indices=(),
        ligand_symbols=(),
        x=x,
        y=y,
        z=z,
        rotation_local_to_global=R,
        pair_indices=(),
        pair_mean_distances=(),
        provider="explicit-vectors",
        mode="full-3d",
        metadata={
            "input_x_vector": np.asarray(x_vector, dtype=float).tolist(),
            "input_z_vector": np.asarray(z_vector, dtype=float).tolist(),
        },
    )
