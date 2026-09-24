from __future__ import annotations

import numpy as np
from pymatgen.core import Structure


def lattice_cosines(lattice_matrix: np.ndarray) -> tuple[float, float, float]:
    """Return cos(alpha), cos(beta), cos(gamma) for row-vector lattice matrices."""
    lattice = np.asarray(lattice_matrix, dtype=float)
    a, b, c = lattice
    la, lb, lc = np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)
    cos_alpha = float(np.dot(b, c) / (lb * lc))
    cos_beta = float(np.dot(a, c) / (la * lc))
    cos_gamma = float(np.dot(a, b) / (la * lb))
    return cos_alpha, cos_beta, cos_gamma


def is_vesta_rhombohedral_like(lattice_matrix: np.ndarray) -> bool:
    """Detect the rhombohedral/trigonal setting used by the legacy VESTA helper.

    This is intentionally a narrow compatibility heuristic, not a general model
    of VESTA's display coordinate system.
    """
    lattice = np.asarray(lattice_matrix, dtype=float)
    lengths = np.linalg.norm(lattice, axis=1)
    cos_alpha, cos_beta, cos_gamma = lattice_cosines(lattice)
    equal_lengths = np.allclose(lengths, lengths[0], rtol=1e-4, atol=1e-4)
    two_obtuse_one_acute = cos_alpha < -0.5 and cos_beta < -0.5 and cos_gamma > 0.5
    return bool(equal_lengths and two_obtuse_one_acute)


def vesta_rhombohedral_transform(lattice_matrix: np.ndarray) -> np.ndarray:
    """Return an orthogonal row-vector transform matching the legacy VESTA view.

    Row vectors transform as ``v_display = v_poscar @ C``. The returned matrix C
    is a proper rotation with det(C)=+1.
    """
    a, b, _ = np.asarray(lattice_matrix, dtype=float)

    z_axis = a + b
    z_axis = z_axis / np.linalg.norm(z_axis)

    diagonal_axis = a - b
    diagonal_axis = diagonal_axis - np.dot(diagonal_axis, z_axis) * z_axis
    diagonal_axis = diagonal_axis / np.linalg.norm(diagonal_axis)

    perpendicular_axis = np.cross(z_axis, diagonal_axis)
    perpendicular_axis = perpendicular_axis / np.linalg.norm(perpendicular_axis)

    x_axis = -(perpendicular_axis + diagonal_axis) / np.sqrt(2.0)
    y_axis = -(perpendicular_axis - diagonal_axis) / np.sqrt(2.0)
    C = np.column_stack([x_axis, y_axis, z_axis])

    if not np.allclose(C.T @ C, np.eye(3), atol=1e-10):
        raise ValueError("VESTA compatibility transform is not orthogonal")
    if not np.isclose(np.linalg.det(C), 1.0, atol=1e-10):
        raise ValueError("VESTA compatibility transform is not a proper rotation")
    return C


def resolve_cartesian_transform(structure: Structure, mode: str = "poscar") -> tuple[np.ndarray, str]:
    """Resolve the working Cartesian frame used by manual frame construction.

    Parameters
    ----------
    mode
        ``poscar``: use the physical Cartesian frame of the POSCAR.
        ``auto``: use the legacy VESTA rhombohedral transform only when the
        narrow compatibility heuristic matches.
        ``vesta``: request that transform explicitly; for non-matching cells an
        informative error is raised rather than silently guessing.
    """
    if mode not in {"poscar", "auto", "vesta"}:
        raise ValueError("cartesian frame must be 'poscar', 'auto', or 'vesta'")

    identity = np.eye(3)
    lattice = np.asarray(structure.lattice.matrix, dtype=float)
    if mode == "poscar":
        return identity, "poscar"

    detected = is_vesta_rhombohedral_like(lattice)
    if mode == "auto":
        if detected:
            return vesta_rhombohedral_transform(lattice), "vesta(auto)"
        return identity, "poscar(auto)"

    if not detected:
        raise ValueError(
            "--cart-frame vesta currently implements only the legacy "
            "rhombohedral/trigonal VESTA compatibility transform; this cell "
            "does not match that heuristic. Use --cart-frame poscar or define "
            "the desired axes explicitly."
        )
    return vesta_rhombohedral_transform(lattice), "vesta"
