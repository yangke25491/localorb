from __future__ import annotations

import numpy as np

# Real d-orbital order used throughout localorb and by common VASP post-processing:
D_ORBITALS = ("dxy", "dyz", "dz2", "dxz", "dx2-y2")


def _normalized_traceless_quadratics() -> np.ndarray:
    """Return an orthonormal Cartesian-tensor representation of real d orbitals.

    A real d angular function can be represented by a symmetric traceless
    quadratic form r^T Q r. Frobenius-normalized Q tensors form an orthonormal
    representation of the l=2 subspace, making spatial rotations straightforward.
    """
    mats = []

    q = np.zeros((3, 3)); q[0, 1] = q[1, 0] = 0.5; mats.append(q)          # xy
    q = np.zeros((3, 3)); q[1, 2] = q[2, 1] = 0.5; mats.append(q)          # yz
    q = np.diag([-0.5, -0.5, 1.0]); mats.append(q)                         # z2
    q = np.zeros((3, 3)); q[0, 2] = q[2, 0] = 0.5; mats.append(q)          # xz
    q = np.diag([0.5, -0.5, 0.0]); mats.append(q)                          # x2-y2

    out = []
    for q in mats:
        out.append(q / np.sqrt(np.sum(q * q)))
    return np.asarray(out)


_D_BASIS = _normalized_traceless_quadratics()


def d_rotation_matrix(rotation_local_to_global: np.ndarray) -> np.ndarray:
    """Return the 5x5 real-d basis transformation for a local frame.

    Parameters
    ----------
    rotation_local_to_global
        3x3 matrix whose columns are local x, y, z unit vectors expressed in
        global Cartesian coordinates.

    Returns
    -------
    T : (5, 5) ndarray
        Rows are local real-d orbitals expanded in the global real-d basis:
        |d_local[a]> = sum_m T[a,m] |d_global[m]>.

    The orbital order is ``D_ORBITALS``.
    """
    R = np.asarray(rotation_local_to_global, dtype=float)
    if R.shape != (3, 3):
        raise ValueError("rotation_local_to_global must be 3x3")
    if not np.allclose(R.T @ R, np.eye(3), atol=1e-8):
        raise ValueError("rotation matrix is not orthonormal")
    if np.linalg.det(R) < 0.0:
        raise ValueError("rotation matrix must be right-handed")

    T = np.empty((5, 5), dtype=float)
    for a, q_local in enumerate(_D_BASIS):
        q_global = R @ q_local @ R.T
        for m, q_basis in enumerate(_D_BASIS):
            T[a, m] = float(np.sum(q_global * q_basis))

    # Numerical sanity: l=2 representation of a proper 3D rotation is orthogonal.
    if not np.allclose(T @ T.T, np.eye(5), atol=1e-8):
        raise RuntimeError("constructed d-orbital rotation is not orthogonal")
    return T


def rotate_d_coefficients(coefficients: np.ndarray, rotation_local_to_global: np.ndarray) -> np.ndarray:
    """Transform complex global d-orbital projection amplitudes to a local frame.

    The final axis must have length 5 in ``D_ORBITALS`` order. Rotation is done
    on amplitudes, not squared weights, so phase/interference information is kept.
    """
    c = np.asarray(coefficients)
    if c.shape[-1] != 5:
        raise ValueError("last coefficient dimension must have length 5")
    T = d_rotation_matrix(rotation_local_to_global)
    return np.einsum("am,...m->...a", T, c)


def rotated_d_weights(coefficients: np.ndarray, rotation_local_to_global: np.ndarray) -> np.ndarray:
    c_local = rotate_d_coefficients(coefficients, rotation_local_to_global)
    return np.real(c_local * np.conjugate(c_local))
