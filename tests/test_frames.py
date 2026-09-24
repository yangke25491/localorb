import numpy as np
from pymatgen.core import Lattice, Structure

from localorb.frames import build_local_frame
from localorb.rotate import rotate_structure_to_frame
from localorb.wannier import render_projection_block


def make_elongated_octahedron():
    lattice = Lattice.cubic(12.0)
    center = np.array([6.0, 6.0, 6.0])
    species = ["Ni"] + ["O"] * 6
    cart = [
        center,
        center + [2.0, 0.0, 0.0],
        center + [-2.0, 0.0, 0.0],
        center + [0.0, 2.0, 0.0],
        center + [0.0, -2.0, 0.0],
        center + [0.0, 0.0, 2.3],
        center + [0.0, 0.0, -2.3],
    ]
    return Structure(lattice, species, cart, coords_are_cartesian=True)


def test_frame_is_orthonormal_and_right_handed():
    s = make_elongated_octahedron()
    f = build_local_frame(s, 0, ligand="O", coordination=6, cutoff=3.0)
    R = f.rotation_local_to_global
    assert np.allclose(R.T @ R, np.eye(3), atol=1e-10)
    assert np.isclose(np.linalg.det(R), 1.0, atol=1e-10)
    assert np.allclose(np.abs(f.z), [0.0, 0.0, 1.0], atol=1e-10)


def test_rotated_structure_preserves_distances():
    s = make_elongated_octahedron()
    f = build_local_frame(s, 0, ligand="O", coordination=6, cutoff=3.0)
    r = rotate_structure_to_frame(s, f)
    before = sorted(s.get_distance(0, i) for i in range(1, 7))
    after = sorted(r.get_distance(0, i) for i in range(1, 7))
    assert np.allclose(before, after, atol=1e-10)


def test_wannier_block_uses_site_specific_axes():
    s = make_elongated_octahedron()
    f = build_local_frame(s, 0, ligand="O", coordination=6, cutoff=3.0)
    text = render_projection_block(s, [f], ["dz2", "dx2-y2"])
    assert "begin projections" in text
    assert ":dz2:z=" in text
    assert ":dx2-y2:z=" in text
    assert "end projections" in text
