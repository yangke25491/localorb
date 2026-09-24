import numpy as np
from pymatgen.core import Lattice, Structure

from localorb.compat_vesta import (
    is_vesta_rhombohedral_like,
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
        [6.0, 5.0, 7.0],
    ]
    return Structure(lattice, species, cart, coords_are_cartesian=True)


def test_manual_full_3d_identity_axes():
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


def test_manual_fixed_z_projects_tilted_bond():
    structure = make_manual_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O3",
        plane_atom="O2",
        mode="fixed-z",
    )
    assert np.allclose(frame.z, [0.0, 0.0, 1.0], atol=1e-12)
    assert np.allclose(frame.x, [1.0, 0.0, 0.0], atol=1e-12)
    assert np.allclose(frame.y, [0.0, 1.0, 0.0], atol=1e-12)


def test_vesta_style_selectors_and_explicit_image():
    structure = make_manual_structure()
    ni = resolve_site_selector(structure, "1", index_base=1)
    oxygen = resolve_site_selector(structure, "O1@0,0,1", index_base=1)
    assert ni.site_index == 0
    assert oxygen.site_index == 1
    assert oxygen.explicit_image
    assert np.array_equal(oxygen.image_shift, [0, 0, 1])


def test_vector_provider_orthonormalizes_axes():
    structure = make_manual_structure()
    frame = build_vector_frame(
        structure,
        center_index=0,
        x_vector=np.array([1.0, 1.0, 0.2]),
        z_vector=np.array([0.0, 0.0, 1.0]),
    )
    R = frame.rotation_local_to_global
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)
    assert np.allclose(frame.z, [0.0, 0.0, 1.0], atol=1e-12)
    assert np.isclose(np.dot(frame.x, frame.z), 0.0, atol=1e-12)


def test_legacy_vesta_rhombohedral_transform_is_proper_rotation():
    # Equal-length lattice with cos(gamma)=+0.6 and cos(alpha)=cos(beta)=-0.6.
    lattice = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.6, 0.8, 0.0],
            [-0.6, -0.3, np.sqrt(0.55)],
        ]
    )
    assert is_vesta_rhombohedral_like(lattice)
    C = vesta_rhombohedral_transform(lattice)
    assert np.allclose(C.T @ C, np.eye(3), atol=1e-12)
    assert np.isclose(np.linalg.det(C), 1.0, atol=1e-12)


def test_rotation_preserves_cartesian_input_geometry():
    structure = make_manual_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O3",
        plane_atom="O2",
        mode="full-3d",
    )
    rotated = rotate_structure_to_frame(structure, frame)
    before = [structure.get_distance(0, i) for i in range(1, len(structure))]
    after = [rotated.get_distance(0, i) for i in range(1, len(rotated))]
    assert np.allclose(before, after, atol=1e-12)
    assert np.isclose(structure.volume, rotated.volume, atol=1e-10)
