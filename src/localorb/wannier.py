from __future__ import annotations

from typing import Iterable

from pymatgen.core import Structure

from .frames import LocalFrame


_VALID_ORBITALS = {
    "s",
    "px", "py", "pz",
    "dxy", "dyz", "dz2", "dxz", "dx2-y2",
}


def _fmt(v) -> str:
    return ",".join(f"{float(x):.10f}" for x in v)


def validate_orbitals(orbitals: Iterable[str]) -> list[str]:
    out = [o.strip() for o in orbitals if o.strip()]
    bad = [o for o in out if o not in _VALID_ORBITALS]
    if bad:
        raise ValueError(f"Unsupported orbital(s): {', '.join(bad)}")
    if not out:
        raise ValueError("At least one orbital must be selected")
    return out


def projection_line(structure: Structure, frame: LocalFrame, orbital: str) -> str:
    """Return one Wannier90 projection using a site-specific local frame.

    We use an explicit fractional-coordinate center so that symmetry-equivalent
    sites may still carry different local axes in distorted/supercell systems.
    """
    site = structure[frame.site_index]
    frac = _fmt(site.frac_coords)
    zaxis = _fmt(frame.z)
    xaxis = _fmt(frame.x)
    return f"f={frac}:{orbital}:z={zaxis}:x={xaxis}"


def render_projection_block(
    structure: Structure,
    frames: Iterable[LocalFrame],
    orbitals: Iterable[str],
    include_comments: bool = True,
) -> str:
    orbitals = validate_orbitals(orbitals)
    lines = ["begin projections"]
    for frame in frames:
        if include_comments:
            lines.append(
                f"! site {frame.site_index} {frame.center_symbol}; "
                f"x={_fmt(frame.x)} z={_fmt(frame.z)}"
            )
        for orbital in orbitals:
            lines.append(projection_line(structure, frame, orbital))
    lines.append("end projections")
    return "\n".join(lines) + "\n"
