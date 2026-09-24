import numpy as np
from pymatgen.core import Lattice, Structure

from localorb.manual import build_manual_frame, resolve_site_selector
from localorb.rotate import rotate_structure_to_frame


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


def test_vesta_style_selector_and_one_based_integer():
    structure = make_manual_structure()
    assert resolve_site_selector(structure, "Ni1").site_index == 0
    assert resolve_site_selector(structure, "O2").site_index == 2
    assert resolve_site_selector(structure, "1").site_index == 0
    assert resolve_site_selector(structure, "3").site_index == 2


def test_manual_full_3d_frame():
    structure = make_manual_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1",
        plane_atom="O2",
        mode="full-3d",
    )
    assert frame.provider == "manual-atoms"
    assert frame.mode == "full-3d"
    assert np.allclose(frame.x, [1.0, 0.0, 0.0], atol=1e-12)
    assert np.allclose(frame.y, [0.0, 1.0, 0.0], atol=1e-12)
    assert np.allclose(frame.z, [0.0, 0.0, 1.0], atol=1e-12)
    assert np.isclose(np.linalg.det(frame.rotation_local_to_global), 1.0)


def test_manual_fixed_z_projects_x_bond():
    lattice = Lattice.cubic(12.0)
    structure = Structure(
        lattice,
        ["Ni", "O", "O"],
        [[6, 6, 6], [8, 8, 8], [4, 8, 6]],
        coords_are_cartesian=True,
    )
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1",
        plane_atom="O2",
        mode="fixed-z",
    )
    root2 = np.sqrt(0.5)
    assert np.allclose(frame.z, [0.0, 0.0, 1.0], atol=1e-12)
    assert np.allclose(frame.x, [root2, root2, 0.0], atol=1e-12)
    assert np.allclose(frame.y, [-root2, root2, 0.0], atol=1e-12)


def test_nearest_periodic_image_is_recorded():
    lattice = Lattice.cubic(10.0)
    structure = Structure(
        lattice,
        ["Ni", "O", "O"],
        [[0.95, 0.5, 0.5], [0.05, 0.5, 0.5], [0.95, 0.7, 0.5]],
        coords_are_cartesian=False,
    )
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1",
        plane_atom="O2",
        mode="fixed-z",
        nearest_image=True,
    )
    assert frame.metadata["x_image"] == [1, 0, 0]
    assert np.allclose(frame.x, [1.0, 0.0, 0.0], atol=1e-12)


def test_explicit_image_overrides_nearest_image():
    lattice = Lattice.cubic(10.0)
    structure = Structure(
        lattice,
        ["Ni", "O", "O"],
        [[0.95, 0.5, 0.5], [0.05, 0.5, 0.5], [0.95, 0.7, 0.5]],
        coords_are_cartesian=False,
    )
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O1@0,0,0",
        plane_atom="O2",
        mode="fixed-z",
        nearest_image=True,
    )
    assert frame.metadata["x_image"] == [0, 0, 0]
    assert np.allclose(frame.x, [-1.0, 0.0, 0.0], atol=1e-12)


def test_manual_rotation_preserves_pair_distances():
    structure = make_manual_structure()
    frame = build_manual_frame(
        structure,
        center="Ni1",
        x_atom="O2",
        plane_atom="O3",
        mode="full-3d",
    )
    rotated = rotate_structure_to_frame(structure, frame)
    for i in range(1, len(structure)):
        assert np.isclose(structure.get_distance(0, i), rotated.get_distance(0, i), atol=1e-12)
