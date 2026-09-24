from __future__ import annotations

from pathlib import Path

import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.io.vasp import Procar

from .frames import LocalFrame
from .orbitals import D_ORBITALS, rotate_d_coefficients


# VASP/Pymatgen conventional orbital order for phase-factor-resolved PROCAR:
# s, py, pz, px, dxy, dyz, dz2, dxz, dx2-y2, ...
_D_SLICE = slice(4, 9)


def rotate_procar_d(
    procar_path: str | Path,
    frames: list[LocalFrame],
) -> dict:
    """Rotate phase-resolved VASP d projections into site-local frames.

    Requires a PROCAR containing complex projection amplitudes/phase factors
    readable by ``pymatgen.io.vasp.Procar.phase_factors`` (typically generated
    with a phase-resolved LORBIT mode such as 12/14, depending on VASP version).

    Returns a dictionary suitable for ``numpy.savez_compressed``.
    """
    procar = Procar(str(procar_path))
    if not getattr(procar, "phase_factors", None):
        raise ValueError(
            "PROCAR does not expose phase_factors. localorb must rotate complex "
            "projection amplitudes, not already-squared orbital weights."
        )

    spins = [s for s in (Spin.up, Spin.down) if s in procar.phase_factors]
    if not spins:
        raise ValueError("No spin channel with phase-resolved projections found")

    site_indices = np.asarray([f.site_index for f in frames], dtype=int)
    rotated_by_spin = []

    for spin in spins:
        phase = np.asarray(procar.phase_factors[spin])
        if phase.ndim != 4:
            raise ValueError(f"Unexpected PROCAR phase_factors shape: {phase.shape}")
        if phase.shape[-1] < 9:
            raise ValueError(
                "PROCAR does not contain the five resolved d channels expected "
                "in columns dxy,dyz,dz2,dxz,dx2-y2"
            )

        per_site = []
        for frame in frames:
            if frame.site_index >= phase.shape[2]:
                raise IndexError(
                    f"Site index {frame.site_index} exceeds PROCAR ion count {phase.shape[2]}"
                )
            c_global = phase[:, :, frame.site_index, _D_SLICE]
            per_site.append(rotate_d_coefficients(c_global, frame.rotation_local_to_global))
        # nk, nb, nsite, 5
        rotated_by_spin.append(np.stack(per_site, axis=2))

    coeff = np.stack(rotated_by_spin, axis=0)
    weight = np.real(coeff * np.conjugate(coeff))
    return {
        "coefficients": coeff,
        "weights": weight,
        "site_indices": site_indices,
        "orbital_names": np.asarray(D_ORBITALS),
        "spin_values": np.asarray([int(s) for s in spins], dtype=int),
    }


def save_rotated_procar_npz(
    output: str | Path,
    procar_path: str | Path,
    frames: list[LocalFrame],
) -> None:
    payload = rotate_procar_d(procar_path, frames)
    np.savez_compressed(output, **payload)
