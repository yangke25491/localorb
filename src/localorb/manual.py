from __future__ import annotations

from dataclasses import dataclass
import re

import numpy as np
from pymatgen.core import Structure

from .compat_vesta import resolve_cartesian_transform
from .frames import LocalFrame


_SELECTOR_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*?)(\d+)$")


@dataclass(frozen=True)
class SiteSelection:
    site_index: int
    image_shift: np.ndarray
    explicit_image: bool
    selector: str


def _unit(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm < 1e-12:
        raise ValueError("Cannot define a local axis from a zero-length vector")
    return vector / norm


def parse_image_shift(text: str) -> np.ndarray:
    cleaned = text.strip().replace("(", "").replace(")", "")
    tokens = cleaned.replace(",", " ").split()
    if len(tokens) != 3:
        raise ValueError(f"Image shift must have three integers, got: {text!r}")
    try:
        return np.array([int(token) for token in tokens], dtype=int)
    except ValueError as exc:
        raise ValueError(f"Image shift must contain integers, got: {text!r}") from exc


def _resolve_generated_label(structure: Structure, label: str) -> int:
    match = _SELECTOR_RE.match(label)
    if match is None:
        raise ValueError(f"Unrecognized atom selector: {label!r}")
    symbol, ordinal_text = match.groups()
    ordinal = int(ordinal_text)
    if ordinal < 1:
        raise ValueError("Element-local atom numbering starts from 1")

    matches = [
        index
        for index, site in enumerate(structure)
        if site.specie.symbol.lower() == symbol.lower()
    ]
    if ordinal > len(matches):
        raise ValueError(
            f"Selector {label!r} requests {symbol}{ordinal}, but the structure contains "
            f"only {len(matches)} {symbol} sites"
        )
    return matches[ordinal - 1]


def resolve_site_selector(structure: Structure, selector: str, index_base: int = 1) -> SiteSelection:
    """Resolve a VASP/VESTA-style atom selector.

    Supported forms are global integer indices and generated labels such as
    ``Ni2`` or ``O7``. An optional periodic image can be appended as
    ``@i,j,k``. By default, bare integer indices are interpreted as 1-based,
    matching POSCAR/VESTA atom numbering.
    """
    selector = selector.strip()
    if not selector:
        raise ValueError("Atom selector cannot be empty")

    if "@" in selector:
        atom_text, image_text = selector.split("@", 1)
        image = parse_image_shift(image_text)
        explicit_image = True
    else:
        atom_text = selector
        image = np.zeros(3, dtype=int)
        explicit_image = False

    atom_text = atom_text.strip()
    if atom_text.lstrip("+-").isdigit():
        raw_index = int(atom_text)
        site_index = raw_index - index_base
        if not 0 <= site_index < len(structure):
            raise ValueError(
                f"Atom index {raw_index} is outside the structure for index_base={index_base}"
            )
    else:
        site_index = _resolve_generated_label(structure, atom_text)

    return SiteSelection(
        site_index=site_index,
        image_shift=image,
        explicit_image=explicit_image,
        selector=selector,
    )


def _nearest_image_shift(
    structure: Structure,
    center_frac: np.ndarray,
    target_site_index: int,
) -> np.ndarray:
    target_frac = np.asarray(structure[target_site_index].frac_coords, dtype=float)
    _, image = structure.lattice.get_distance_and_image(center_frac, target_frac)
    return np.asarray(image, dtype=int)


def _materialize_selection(
    structure: Structure,
    selection: SiteSelection,
    center_frac: np.ndarray | None,
    nearest_image: bool,
) -> tuple[np.ndarray, np.ndarray]:
    shift = np.asarray(selection.image_shift, dtype=int)
    if nearest_image and center_frac is not None and not selection.explicit_image:
        shift = _nearest_image_shift(structure, center_frac, selection.site_index)
    frac = np.asarray(structure[selection.site_index].frac_coords, dtype=float) + shift
    cart = np.asarray(structure.lattice.get_cartesian_coords(frac), dtype=float)
    return shift, cart


def build_manual_frame(
    structure: Structure,
    center: str,
    x_atom: str,
    plane_atom: str,
    mode: str = "full-3d",
    nearest_image: bool = True,
    index_base: int = 1,
    cartesian_frame: str = "poscar",
) -> LocalFrame:
    """Build a local frame from three explicitly selected atoms.

    ``cartesian_frame`` controls only the working Cartesian coordinates used to
    interpret the selected bonds. ``poscar`` is the safe default. ``auto`` and
    ``vesta`` enable the narrow rhombohedral/trigonal VESTA compatibility layer
    ported from the user's legacy rotation script. The final LocalFrame is always
    returned in the physical POSCAR Cartesian frame.

    In ``fixed-z`` mode the physically fixed direction is always the POSCAR
    Cartesian z axis. If a VESTA working frame is active, that physical z axis is
    first expressed in the working frame, matching the legacy script's behavior.
    """
    if mode not in {"full-3d", "fixed-z"}:
        raise ValueError("mode must be 'full-3d' or 'fixed-z'")

    center_sel = resolve_site_selector(structure, center, index_base=index_base)
    x_sel = resolve_site_selector(structure, x_atom, index_base=index_base)
    plane_sel = resolve_site_selector(structure, plane_atom, index_base=index_base)

    center_shift, center_cart = _materialize_selection(
        structure, center_sel, center_frac=None, nearest_image=False
    )
    center_frac = np.asarray(structure[center_sel.site_index].frac_coords, dtype=float) + center_shift
    x_shift, x_cart = _materialize_selection(
        structure, x_sel, center_frac=center_frac, nearest_image=nearest_image
    )
    plane_shift, plane_cart = _materialize_selection(
        structure, plane_sel, center_frac=center_frac, nearest_image=nearest_image
    )

    C, resolved_cartesian_frame = resolve_cartesian_transform(structure, cartesian_frame)
    center_work = center_cart @ C
    x_work = x_cart @ C
    plane_work = plane_cart @ C

    x_vector = x_work - center_work
    plane_vector = plane_work - center_work

    if mode == "full-3d":
        x_w = _unit(x_vector)
        y_seed = plane_vector - float(np.dot(plane_vector, x_w)) * x_w
        y_w = _unit(y_seed)
        z_w = _unit(np.cross(x_w, y_w))
        y_w = _unit(np.cross(z_w, x_w))
    else:
        z_poscar = np.array([0.0, 0.0, 1.0], dtype=float)
        z_w = _unit(z_poscar @ C)
        x_seed = x_vector - float(np.dot(x_vector, z_w)) * z_w
        x_w = _unit(x_seed)
        y_w = _unit(np.cross(z_w, x_w))

    # Row-vector convention: v_work = v_poscar @ C, therefore an axis expressed
    # in the working frame maps back as axis_poscar = axis_work @ C.T.
    x = _unit(x_w @ C.T)
    y = _unit(y_w @ C.T)
    z = _unit(z_w @ C.T)
    R = np.column_stack((x, y, z))

    if not np.allclose(R.T @ R, np.eye(3), atol=1e-10):
        raise ValueError("Constructed manual frame is not orthonormal")
    if not np.isclose(np.linalg.det(R), 1.0, atol=1e-10):
        raise ValueError("Constructed manual frame is not a proper right-handed rotation")

    plane_projected = plane_vector - float(np.dot(plane_vector, z_w)) * z_w
    if np.linalg.norm(plane_projected) > 1e-12:
        plane_alignment = float(np.dot(_unit(plane_projected), y_w))
    else:
        plane_alignment = float("nan")

    ligand_indices = (x_sel.site_index, plane_sel.site_index)
    ligand_symbols = (
        structure[x_sel.site_index].specie.symbol,
        structure[plane_sel.site_index].specie.symbol,
    )

    return LocalFrame(
        site_index=center_sel.site_index,
        center_symbol=structure[center_sel.site_index].specie.symbol,
        ligand_indices=ligand_indices,
        ligand_symbols=ligand_symbols,
        x=x,
        y=y,
        z=z,
        rotation_local_to_global=R,
        pair_indices=(),
        pair_mean_distances=(),
        provider="manual-atoms",
        mode=mode,
        metadata={
            "center_selector": center,
            "x_selector": x_atom,
            "plane_selector": plane_atom,
            "center_image": center_shift.tolist(),
            "x_image": x_shift.tolist(),
            "plane_image": plane_shift.tolist(),
            "nearest_image": bool(nearest_image),
            "index_base": int(index_base),
            "cartesian_frame": resolved_cartesian_frame,
            "cartesian_transform": C.tolist(),
            "fixed_z_reference": "poscar-cartesian-z" if mode == "fixed-z" else None,
            "plane_y_alignment": plane_alignment,
        },
    )
