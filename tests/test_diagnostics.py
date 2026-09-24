import numpy as np
from pymatgen.core import Lattice, Structure

from localorb.diagnostics import diagnose_frame
from localorb.frames import build_local_frame


def make_octahedron(z_distance=2.3, x_distance=2.0, y_distance=2.0):
    lattice = Lattice.cubic(12.0)
    center = np.array([6.0, 6.0, 6.0])
    species = ["Ni"] + ["O"] * 6
    cart = [
        center,
        center + [x_distance, 0.0, 0.0],
        center + [-x_distance, 0.0, 0.0],
        center + [0.0, y_distance, 0.0],
        center + [0.0, -y_distance, 0.0],
        center + [0.0, 0.0, z_distance],
        center + [0.0, 0.0, -z_distance],
    ]
    return Structure(lattice, species, cart, coords_are_cartesian=True)


def test_elongated_octahedron_has_clean_geometry_and_unique_z():
    structure = make_octahedron(z_distance=2.4)
    frame = build_local_frame(structure, 0, ligand="O", cutoff=3.0)
    quality = frame.metadata["quality"]
    assert quality["max_opposition_error"] < 1e-12
    assert quality["raw_axis_orthogonality_error"] < 1e-12
    assert quality["z_gap_fraction"] > 0.02

    diag = diagnose_frame(frame)
    assert diag.status == "PASS"
    assert not diag.errors


def test_cubic_octahedron_warns_that_z_choice_is_degenerate():
    structure = make_octahedron(z_distance=2.0, x_distance=2.0, y_distance=2.0)
    frame = build_local_frame(structure, 0, ligand="O", cutoff=3.0)
    diag = diagnose_frame(frame)
    assert diag.status == "WARN"
    assert any("z-axis assignment is nearly degenerate" in msg for msg in diag.warnings)


def test_distorted_opposite_pair_can_be_flagged_without_invalidating_rotation():
    structure = make_octahedron(z_distance=2.4)
    # Bend one +x ligand away from exact opposition while keeping sixfold shell.
    structure.translate_sites([1], [0.0, 0.7, 0.0], frac_coords=False)
    frame = build_local_frame(structure, 0, ligand="O", cutoff=3.0)
    diag = diagnose_frame(frame, max_opposition_error_warn=0.01)

    assert np.allclose(
        frame.rotation_local_to_global.T @ frame.rotation_local_to_global,
        np.eye(3),
        atol=1e-10,
    )
    assert diag.status in {"PASS", "WARN"}
    if diag.metrics.get("max_opposition_error", 0.0) > 0.01:
        assert any("not strongly opposite" in msg for msg in diag.warnings)
