import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from localorb.compat_vesta import (
    is_vesta_rhombohedral_like,
    resolve_cartesian_transform,
    vesta_rhombohedral_transform,
)
from localorb.manual import build_manual_frame, resolve_site_selector
from localorb.rotate import rotate_structure_to_frame
from localorb.vectors import build_vector_frame


def make_manual_structure():
    lattice = Lattice.cubic(10.0)
    species = ["Ni", "O", "O", "O"]
    cart = [
        [5.0, 5.0, 5.0],
        [7.0, 5.0, 5.0],
        [5.0, 7.0, 5.0],
        [5.0, 5.0, 7.0],
    ]
    return Structure(lattice, species, cart, coords_are_cartesian=True)


def make_boundary_structure():
    lattice = Lattice.cubic(10.0)
    species = ["Ni", "O", "O"]
    frac = [
        [0.95, 0.50, 0.50],
        [0.05, 0.50, 0.50],
        [0.95, 0.70, 0.50],
    ]
    return Structure(lattice, species, frac)


def make_rhombohedral_like_structure():
    # Equal-length lattice vectors with the same narrow angular pattern used by
    # the legacy VESTA compatibility heuristic:
    # cos(alpha), cos(beta) < -0.5 and cos(gamma) > 0.5.
    lattice_matrix = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.8, 0.6, 0.0],
            [-0.8, -0.2, np.sqrt(0.32)],
        ]
    ) * 6.0
    return Structure(Lattice(lattice_matrix), ["Ni"], [[0.0, 0.0, 0.0]])


def test_manual_full_3d_frame_matches_selected_bonds():
    structure = make_manual_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1",
        plane_atom="O2",
        mode="full-3d",
    )
    assert np.allclose(frame.x, [1.0, 0.0, 0.0])
    assert np.allclose(frame.y, [0.0, 1.0, 0.0])
    assert np.allclose(frame.z, [0.0, 0.0, 1.0])
    assert np.isclose(np.linalg.det(frame.rotation_local_to_global), 1.0)


def test_manual_fixed_z_keeps_poscar_cartesian_z():
    structure = make_manual_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1",
        plane_atom="O2",
        mode="fixed-z",
    )
    assert np.allclose(frame.z, [0.0, 0.0, 1.0])
    assert np.allclose(frame.x, [1.0, 0.0, 0.0])


def test_generated_labels_and_one_based_integer_selectors_agree():
    structure = make_manual_structure()
    ni_label = resolve_site_selector(structure, "Ni1")
    ni_integer = resolve_site_selector(structure, "1", index_base=1)
    o_label = resolve_site_selector(structure, "O2")
    o_integer = resolve_site_selector(structure, "3", index_base=1)
    assert ni_label.site_index == ni_integer.site_index == 0
    assert o_label.site_index == o_integer.site_index == 2


def test_nearest_periodic_image_is_used_for_manual_bond():
    structure = make_boundary_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1",
        plane_atom="O2",
        mode="fixed-z",
        nearest_image=True,
    )
    # O1 at x=0.05 should be interpreted as the image at x=1.05 relative to
    # Ni at x=0.95, so local +x points along Cartesian +x.
    assert np.allclose(frame.x, [1.0, 0.0, 0.0])
    assert frame.metadata["x_image"] == [1, 0, 0]


def test_explicit_image_overrides_nearest_image_choice():
    structure = make_boundary_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1@0,0,0",
        plane_atom="O2",
        mode="fixed-z",
        nearest_image=True,
    )
    assert frame.metadata["x_image"] == [0, 0, 0]
    assert frame.x[0] < 0.0


def test_vector_provider_orthonormalizes_user_axes():
    structure = make_manual_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 1.0, 0.25]),
        z_vector=np.array([0.0, 0.0, 2.0]),
    )
    R = frame.rotation_local_to_global
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)
    assert np.allclose(frame.z, [0.0, 0.0, 1.0])
    assert abs(float(np.dot(frame.x, frame.z))) < 1e-12


def test_vesta_compatibility_transform_is_proper_rotation():
    structure = make_rhombohedral_like_structure()
    assert is_vesta_rhombohedral_like(structure.lattice.matrix)
    C = vesta_rhombohedral_transform(structure.lattice.matrix)
    assert np.allclose(C.T @ C, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(C), 1.0, atol=1e-12)
    resolved, name = resolve_cartesian_transform(structure, "vesta")
    assert name == "vesta"
    assert np.allclose(resolved, C)


def test_vesta_mode_rejects_unrelated_cells_instead_of_guessing():
    structure = make_manual_structure()
    with pytest.raises(ValueError, match="legacy rhombohedral/trigonal"):
        resolve_cartesian_transform(structure, "vesta")


def test_structure_rotation_preserves_fractional_coordinates_and_volume():
    structure = make_manual_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 1.0, 0.0]),
        z_vector=np.array([0.0, 0.0, 1.0]),
    )
    rotated = rotate_structure_to_frame(structure, frame)
    assert np.allclose(rotated.frac_coords, structure.frac_coords, atol=1e-12)
    assert np.isclose(rotated.volume, structure.volume, rtol=1e-12, atol=1e-12)
