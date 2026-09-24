from __future__ import annotations

from dataclasses import dataclass, asdict

import numpy as np

from .frames import LocalFrame


@dataclass(frozen=True)
class FrameDiagnostics:
    site_index: int
    provider: str
    mode: str
    orthonormal_error: float
    determinant: float
    status: str
    warnings: tuple[str, ...]
    errors: tuple[str, ...]
    metrics: dict[str, float]

    def as_dict(self) -> dict:
        return asdict(self)


def diagnose_frame(
    frame: LocalFrame,
    *,
    orthonormal_tol: float = 1e-8,
    determinant_tol: float = 1e-8,
    max_opposition_error_warn: float = 0.10,
    max_raw_axis_dot_warn: float = 0.15,
    min_z_gap_fraction_warn: float = 0.02,
) -> FrameDiagnostics:
    """Evaluate numerical validity and geometric ambiguity of a local frame.

    Hard failures are reserved for an invalid rotation matrix. Geometry-dependent
    quantities are warnings because distorted structures may be physically valid.

    For the automatic octahedral provider, three metadata metrics are used when
    available:

    ``max_opposition_error``
        ``max(1 + cos(theta_pair))`` for the three chosen opposite ligand pairs.
        Zero is ideal.

    ``raw_axis_orthogonality_error``
        Maximum absolute dot product between the three ligand-defined axes before
        Gram-Schmidt orthogonalization. Zero is ideal.

    ``z_gap_fraction``
        Relative separation between the selected z-pair mean bond length and the
        nearest competing pair. A small value means the geometric choice of which
        axis should be called z is nearly degenerate, not necessarily incorrect.
    """
    R = np.asarray(frame.rotation_local_to_global, dtype=float)
    gram = R.T @ R
    orth_error = float(np.max(np.abs(gram - np.eye(3))))
    determinant = float(np.linalg.det(R))

    warnings: list[str] = []
    errors: list[str] = []
    metrics: dict[str, float] = {}

    if orth_error > orthonormal_tol:
        errors.append(
            f"rotation matrix is not orthonormal: max|R^T R-I|={orth_error:.3e}"
        )
    if abs(determinant - 1.0) > determinant_tol:
        errors.append(f"rotation matrix is not a proper rotation: det(R)={determinant:.12f}")

    quality = frame.metadata.get("quality") if isinstance(frame.metadata, dict) else None
    if isinstance(quality, dict):
        for name in (
            "max_opposition_error",
            "raw_axis_orthogonality_error",
            "z_gap_fraction",
        ):
            value = quality.get(name)
            if value is not None:
                metrics[name] = float(value)

        opposition = metrics.get("max_opposition_error")
        if opposition is not None and opposition > max_opposition_error_warn:
            warnings.append(
                "selected ligand pairs are not strongly opposite "
                f"(max 1+cos(theta)={opposition:.3f})"
            )

        axis_error = metrics.get("raw_axis_orthogonality_error")
        if axis_error is not None and axis_error > max_raw_axis_dot_warn:
            warnings.append(
                "ligand-defined axes are substantially non-orthogonal before "
                f"orthogonalization (max |dot|={axis_error:.3f})"
            )

        z_gap = metrics.get("z_gap_fraction")
        if z_gap is not None and z_gap < min_z_gap_fraction_warn:
            warnings.append(
                "automatic z-axis assignment is nearly degenerate in bond length "
                f"(relative gap={z_gap:.3e}); inspect or use a manual/vector frame"
            )

    if errors:
        status = "FAIL"
    elif warnings:
        status = "WARN"
    else:
        status = "PASS"

    return FrameDiagnostics(
        site_index=frame.site_index,
        provider=frame.provider,
        mode=frame.mode,
        orthonormal_error=orth_error,
        determinant=determinant,
        status=status,
        warnings=tuple(warnings),
        errors=tuple(errors),
        metrics=metrics,
    )
