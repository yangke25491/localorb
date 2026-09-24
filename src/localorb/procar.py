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

    The returned dictionary contains both the rotated numerical data and the
    exact local-frame matrices needed to reproduce the transformation.
    """
    if not frames:
        raise ValueError("At least one local frame is required")
    site_indices_list = [frame.site_index for frame in frames]
    if len(set(site_indices_list)) != len(site_indices_list):
        raise ValueError("PROCAR rotation received duplicate local frames for the same site")

    procar = Procar(str(procar_path))
    if not getattr(procar, "phase_factors", None):
        raise ValueError(
            "PROCAR does not expose phase_factors. localorb must rotate complex "
            "projection amplitudes, not already-squared orbital weights."
        )

    spins = [s for s in (Spin.up, Spin.down) if s in procar.phase_factors]
    if not spins:
        raise ValueError("No spin channel with phase-resolved projections found")

    site_indices = np.asarray(site_indices_list, dtype=int)
    rotated_by_spin = []
    conservation_errors = []

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
        per_site_conservation = []
        for frame in frames:
            if frame.site_index >= phase.shape[2]:
                raise IndexError(
                    f"Site index {frame.site_index} exceeds PROCAR ion count {phase.shape[2]}"
                )
            c_global = phase[:, :, frame.site_index, _D_SLICE]
            c_local = rotate_d_coefficients(c_global, frame.rotation_local_to_global)
            per_site.append(c_local)

            global_sum = np.sum(np.real(c_global * np.conjugate(c_global)), axis=-1)
            local_sum = np.sum(np.real(c_local * np.conjugate(c_local)), axis=-1)
            per_site_conservation.append(float(np.max(np.abs(local_sum - global_sum))))

        # nk, nb, nsite, 5
        rotated_by_spin.append(np.stack(per_site, axis=2))
        conservation_errors.append(per_site_conservation)

    coeff = np.stack(rotated_by_spin, axis=0)
    weight = np.real(coeff * np.conjugate(coeff))
    rotations = np.stack([frame.rotation_local_to_global for frame in frames], axis=0)

    return {
        "coefficients": coeff,
        "weights": weight,
        "site_indices": site_indices,
        "center_symbols": np.asarray([frame.center_symbol for frame in frames]),
        "providers": np.asarray([frame.provider for frame in frames]),
        "modes": np.asarray([frame.mode for frame in frames]),
        "rotation_local_to_global": rotations,
        "rotation_global_to_local": np.transpose(rotations, (0, 2, 1)),
        "orbital_names": np.asarray(D_ORBITALS),
        "spin_values": np.asarray([int(s) for s in spins], dtype=int),
        "weight_conservation_max_abs": np.asarray(conservation_errors, dtype=float),
        "source_procar": np.asarray(str(procar_path)),
    }


def save_rotated_procar_npz(
    output: str | Path,
    procar_path: str | Path,
    frames: list[LocalFrame],
) -> None:
    payload = rotate_procar_d(procar_path, frames)
    np.savez_compressed(output, **payload)
