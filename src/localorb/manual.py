from __future__ import annotations

from dataclasses import dataclass
import re

import numpy as np
from pymatgen.core import Structure

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


def _selection_frac_coords(structure: Structure, selection: SiteSelection) -> np.ndarray:
    return np.asarray(structure[selection.site_index].frac_coords, dtype=float) + selection.image_shift


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
) -> LocalFrame:
    """Build a local frame from three explicitly selected atoms.

    Parameters
    ----------
    center
        Center atom selector, e.g. ``Ni2`` or ``6``.
    x_atom
        Atom whose center->atom bond defines +x (full-3d), or whose in-plane
        projection defines +x (fixed-z).
    plane_atom
        Second atom used to define the xy plane and the sign/orientation of y.
    mode
        ``full-3d`` constructs x from the selected bond, y by Gram-Schmidt from
        the plane atom, and z=x×y. ``fixed-z`` keeps Cartesian z=(0,0,1),
        projects the selected x bond into the xy plane, then sets y=z×x.
    nearest_image
        If true, selectors without an explicit ``@i,j,k`` image are moved to
        the nearest periodic image relative to the selected center.
    index_base
        Base for bare integer selectors. Default 1 matches POSCAR/VESTA.
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

    x_vector = x_cart - center_cart
    plane_vector = plane_cart - center_cart

    if mode == "full-3d":
        x = _unit(x_vector)
        y_seed = plane_vector - float(np.dot(plane_vector, x)) * x
        y = _unit(y_seed)
        z = _unit(np.cross(x, y))
        # Recompute y to make the frame exactly orthonormal and right handed.
        y = _unit(np.cross(z, x))
    else:
        z = np.array([0.0, 0.0, 1.0], dtype=float)
        x_seed = x_vector - float(np.dot(x_vector, z)) * z
        x = _unit(x_seed)
        y = _unit(np.cross(z, x))

    R = np.column_stack((x, y, z))
    if not np.allclose(R.T @ R, np.eye(3), atol=1e-10):
        raise ValueError("Constructed manual frame is not orthonormal")
    if not np.isclose(np.linalg.det(R), 1.0, atol=1e-10):
        raise ValueError("Constructed manual frame is not a proper right-handed rotation")

    plane_projected = plane_vector - float(np.dot(plane_vector, z)) * z
    if np.linalg.norm(plane_projected) > 1e-12:
        plane_alignment = float(np.dot(_unit(plane_projected), y))
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
            "plane_y_alignment": plane_alignment,
        },
    )
