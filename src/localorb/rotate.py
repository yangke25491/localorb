from __future__ import annotations

import numpy as np
from pymatgen.core import Lattice, Structure

from .frames import LocalFrame


def rotate_structure_to_frame(structure: Structure, frame: LocalFrame) -> Structure:
    """Rigidly rotate lattice and Cartesian positions so frame becomes global xyz.

    If R maps local -> global, then R.T maps global -> local. Applying R.T to
    every Cartesian vector produces a representation in which x',y',z' coincide
    with global x,y,z. Fractional coordinates are preserved under a simultaneous
    rotation of lattice and Cartesian positions.
    """
    Q = frame.rotation_local_to_global.T
    old_lat = np.asarray(structure.lattice.matrix, dtype=float)
    # pymatgen lattice vectors are rows; active Cartesian rotation acts on each row.
    new_lat = old_lat @ Q.T
    new_structure = Structure(
        Lattice(new_lat),
        [site.specie for site in structure],
        structure.frac_coords,
        coords_are_cartesian=False,
        site_properties=structure.site_properties,
    )
    return new_structure


def frame_alignment_error(frame: LocalFrame) -> float:
    Q = frame.rotation_local_to_global.T
    transformed = Q @ frame.rotation_local_to_global
    return float(np.max(np.abs(transformed - np.eye(3))))
