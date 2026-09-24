import numpy as np

from localorb.orbitals import D_ORBITALS, d_rotation_matrix, rotate_d_coefficients


def rot_z(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def test_identity_rotation():
    T = d_rotation_matrix(np.eye(3))
    assert np.allclose(T, np.eye(5), atol=1e-12)


def test_d_rotation_is_orthogonal():
    T = d_rotation_matrix(rot_z(0.37))
    assert np.allclose(T @ T.T, np.eye(5), atol=1e-10)
    assert np.isclose(np.linalg.det(T), 1.0, atol=1e-10)


def test_45_degree_z_rotation_maps_x2y2_to_xy_subspace():
    T = d_rotation_matrix(rot_z(np.pi / 4.0))
    idx_xy = D_ORBITALS.index("dxy")
    idx_x2y2 = D_ORBITALS.index("dx2-y2")
    row = T[idx_x2y2]
    assert np.isclose(abs(row[idx_xy]), 1.0, atol=1e-10)
    mask = np.ones(5, dtype=bool)
    mask[idx_xy] = False
    assert np.allclose(row[mask], 0.0, atol=1e-10)


def test_complex_amplitudes_preserve_total_d_weight():
    rng = np.random.default_rng(7)
    c = rng.normal(size=(4, 3, 5)) + 1j * rng.normal(size=(4, 3, 5))
    out = rotate_d_coefficients(c, rot_z(0.61))
    assert np.allclose(np.sum(abs(c) ** 2, axis=-1), np.sum(abs(out) ** 2, axis=-1))
