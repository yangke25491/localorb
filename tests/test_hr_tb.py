import json

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from localorb.hr import Wannier90HR, read_hr, write_hr
from localorb.tb import build_basis_transform, rotate_hr_basis
from localorb.vectors import build_vector_frame


def make_structure():
    lattice = Lattice.cubic(10.0)
    return Structure(lattice, ["Ni"], [[0.5, 0.5, 0.5]])


def make_hr(num_wann=5, nrpts=2):
    rng = np.random.default_rng(1234)
    hamiltonians = np.zeros((nrpts, num_wann, num_wann), dtype=complex)

    # R=0 block is Hermitian so its spectrum can be compared directly.
    a = rng.normal(size=(num_wann, num_wann)) + 1j * rng.normal(
        size=(num_wann, num_wann)
    )
    hamiltonians[0] = 0.5 * (a + a.conjugate().T)

    # A generic hopping block for a nonzero R.
    if nrpts > 1:
        hamiltonians[1] = rng.normal(size=(num_wann, num_wann)) + 1j * rng.normal(
            size=(num_wann, num_wann)
        )

    return Wannier90HR(
        comment="localorb test hr",
        num_wann=num_wann,
        degeneracies=np.ones(nrpts, dtype=int),
        r_vectors=np.array([[0, 0, 0], [1, 0, 0]][:nrpts], dtype=int),
        hamiltonians=hamiltonians,
    )


def test_hr_roundtrip_preserves_data(tmp_path):
    original = make_hr()
    path = tmp_path / "wannier90_hr.dat"
    write_hr(path, original)
    parsed = read_hr(path)

    assert parsed.comment == original.comment
    assert parsed.num_wann == original.num_wann
    assert np.array_equal(parsed.degeneracies, original.degeneracies)
    assert np.array_equal(parsed.r_vectors, original.r_vectors)
    assert np.allclose(parsed.hamiltonians, original.hamiltonians, atol=1e-9)


def test_full_d_basis_rotation_is_unitary_and_preserves_spectrum():
    structure = make_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 1.0, 0.0]),
        z_vector=np.array([0.0, 0.0, 1.0]),
    )
    basis_map = {
        "num_wann": 5,
        "index_base": 1,
        "groups": [
            {
                "center_atom": "Ni1",
                "indices": [1, 2, 3, 4, 5],
                "orbitals": ["dxy", "dyz", "dz2", "dxz", "dx2-y2"],
            }
        ],
    }

    B, metadata = build_basis_transform(structure, [frame], basis_map, num_wann=5)
    assert np.allclose(B @ B.conjugate().T, np.eye(5), atol=1e-12)
    assert metadata["unitarity_error"] < 1e-12

    hr = make_hr()
    rotated, diagnostics = rotate_hr_basis(hr, B)
    evals_old = np.linalg.eigvalsh(hr.hamiltonians[0])
    evals_new = np.linalg.eigvalsh(rotated.hamiltonians[0])
    assert np.allclose(evals_old, evals_new, atol=1e-10)
    assert diagnostics["unitarity_error"] < 1e-12
    assert diagnostics["max_frobenius_norm_change"] < 1e-10


def test_orbital_order_permutation_is_supported():
    structure = make_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 0.0, 0.0]),
        z_vector=np.array([0.0, 0.0, 1.0]),
    )
    order = ["dz2", "dx2-y2", "dxy", "dxz", "dyz"]
    basis_map = {
        "num_wann": 5,
        "groups": [
            {
                "center_atom": "Ni1",
                "indices": [1, 2, 3, 4, 5],
                "orbitals": order,
            }
        ],
    }
    B, _ = build_basis_transform(structure, [frame], basis_map, num_wann=5)
    assert np.allclose(B, np.eye(5), atol=1e-12)


def test_reduced_eg_basis_is_rejected_for_general_rotation():
    structure = make_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 1.0, 0.0]),
        z_vector=np.array([0.0, 0.0, 1.0]),
    )
    basis_map = {
        "num_wann": 2,
        "groups": [
            {
                "center_atom": "Ni1",
                "indices": [1, 2],
                "orbitals": ["dz2", "dx2-y2"],
            }
        ],
    }
    with pytest.raises(ValueError, match="complete five-d-orbital"):
        build_basis_transform(structure, [frame], basis_map, num_wann=2)


def test_basis_map_num_wann_must_match_hr_size():
    structure = make_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 0.0, 0.0]),
        z_vector=np.array([0.0, 0.0, 1.0]),
    )
    basis_map = {
        "num_wann": 10,
        "groups": [
            {
                "center_atom": "Ni1",
                "indices": [1, 2, 3, 4, 5],
            }
        ],
    }
    with pytest.raises(ValueError, match="declares num_wann=10"):
        build_basis_transform(structure, [frame], basis_map, num_wann=5)
