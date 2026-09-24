import json

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar

from localorb.cli import main
from localorb.specfile import build_frames_from_spec_file


def make_two_site_structure():
    lattice = Lattice.cubic(20.0)
    species = ["Ni", "Ni", "O", "O", "O", "O"]
    cart = [
        [5.0, 5.0, 5.0],
        [15.0, 15.0, 15.0],
        [7.0, 5.0, 5.0],
        [5.0, 7.0, 5.0],
        [17.0, 15.0, 15.0],
        [15.0, 17.0, 15.0],
    ]
    return Structure(lattice, species, cart, coords_are_cartesian=True)


def write_spec(path):
    payload = {
        "defaults": {
            "index_base": 1,
            "nearest_image": True,
            "cart_frame": "poscar",
        },
        "frames": [
            {
                "provider": "manual",
                "center_atom": "Ni1",
                "x_atom": "O1",
                "plane_atom": "O2",
                "mode": "full-3d",
            },
            {
                "provider": "vectors",
                "center_atom": "Ni2",
                "x_vector": [1.0, 0.0, 0.0],
                "z_vector": [0.0, 0.0, 1.0],
            },
        ],
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def test_mixed_multi_site_spec_builds_one_frame_per_center(tmp_path):
    structure = make_two_site_structure()
    spec = tmp_path / "frames.json"
    write_spec(spec)

    frames = build_frames_from_spec_file(structure, spec)
    assert [frame.site_index for frame in frames] == [0, 1]
    assert [frame.provider for frame in frames] == ["manual-atoms", "explicit-vectors"]
    for frame in frames:
        R = frame.rotation_local_to_global
        assert np.allclose(R.T @ R, np.eye(3), atol=1e-12)
        assert np.isclose(np.linalg.det(R), 1.0, atol=1e-12)


def test_array_form_periodic_images_are_accepted(tmp_path):
    structure = make_two_site_structure()
    spec = tmp_path / "images.json"
    spec.write_text(
        json.dumps(
            {
                "frames": [
                    {
                        "provider": "manual",
                        "center_atom": "Ni1",
                        "x_atom": "O1",
                        "plane_atom": "O2",
                        "center_image": [0, 0, 0],
                        "x_image": [0, 0, 0],
                        "plane_image": [0, 0, 0],
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    frame = build_frames_from_spec_file(structure, spec)[0]
    assert frame.metadata["center_image"] == [0, 0, 0]
    assert frame.metadata["x_image"] == [0, 0, 0]
    assert frame.metadata["plane_image"] == [0, 0, 0]


def test_duplicate_center_sites_are_rejected(tmp_path):
    structure = make_two_site_structure()
    spec = tmp_path / "duplicate.json"
    spec.write_text(
        json.dumps(
            {
                "frames": [
                    {
                        "provider": "vectors",
                        "center_atom": "Ni1",
                        "x_vector": [1, 0, 0],
                        "z_vector": [0, 0, 1],
                    },
                    {
                        "provider": "vectors",
                        "center_atom": "Ni1",
                        "x_vector": [0, 1, 0],
                        "z_vector": [0, 0, 1],
                    },
                ]
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="more than once"):
        build_frames_from_spec_file(structure, spec)


def test_frames_file_drives_wannier_cli_for_multiple_sites(tmp_path):
    structure = make_two_site_structure()
    poscar = tmp_path / "POSCAR"
    Poscar(structure).write_file(poscar)
    spec = tmp_path / "frames.json"
    write_spec(spec)
    output = tmp_path / "projections.win"

    rc = main(
        [
            "wannier",
            str(poscar),
            "--frames-file",
            str(spec),
            "--orbitals",
            "dz2,dx2-y2",
            "-o",
            str(output),
        ]
    )
    assert rc == 0
    text = output.read_text(encoding="utf-8")
    assert text.count(":dz2:z=") == 2
    assert text.count(":dx2-y2:z=") == 2
