from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

import numpy as np
from pymatgen.core import Structure


@dataclass(frozen=True)
class LocalFrame:
    site_index: int
    center_symbol: str
    ligand_indices: tuple[int, ...]
    ligand_symbols: tuple[str, ...]
    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    rotation_local_to_global: np.ndarray
    pair_indices: tuple[tuple[int, int], ...]
    pair_mean_distances: tuple[float, ...]
    provider: str = "auto-octahedral"
    mode: str = "full-3d"
    metadata: dict[str, object] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "site_index": self.site_index,
            "center_symbol": self.center_symbol,
            "provider": self.provider,
            "mode": self.mode,
            "ligand_indices": list(self.ligand_indices),
            "ligand_symbols": list(self.ligand_symbols),
            "x": self.x.tolist(),
            "y": self.y.tolist(),
            "z": self.z.tolist(),
            "rotation_local_to_global": self.rotation_local_to_global.tolist(),
            "rotation_global_to_local": self.rotation_local_to_global.T.tolist(),
            "opposite_pairs": [list(p) for p in self.pair_indices],
            "pair_mean_distances": list(self.pair_mean_distances),
            "metadata": self.metadata,
        }


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError("Cannot normalize a zero-length vector")
    return v / n


def _canonical_sign(v: np.ndarray) -> np.ndarray:
    """Choose a deterministic sign for an axis without changing its line."""
    v = np.asarray(v, dtype=float)
    k = int(np.argmax(np.abs(v)))
    return v if v[k] >= 0 else -v


def find_center_sites(structure: Structure, element: str) -> list[int]:
    return [i for i, site in enumerate(structure) if site.specie.symbol == element]


def _candidate_ligands(
    structure: Structure,
    center_index: int,
    ligand: str | None,
    coordination: int,
    cutoff: float | None,
) -> list[tuple[int, np.ndarray, float, str]]:
    center = structure[center_index]
    if cutoff is None:
        # Use a generous shell, then keep the nearest requested neighbors.
        cutoff = 5.0
    neigh = structure.get_neighbors(center, cutoff)
    rows: list[tuple[int, np.ndarray, float, str]] = []
    for nn in neigh:
        symbol = nn.specie.symbol
        if ligand is not None and symbol != ligand:
            continue
        vec = np.asarray(nn.coords - center.coords, dtype=float)
        rows.append((int(nn.index), vec, float(nn.nn_distance), symbol))
    rows.sort(key=lambda r: r[2])
    if len(rows) < coordination:
        raise ValueError(
            f"Site {center_index} has only {len(rows)} matching neighbors within "
            f"{cutoff:.3f} Å; need {coordination}. Increase --cutoff or check --ligand."
        )
    return rows[:coordination]


def _best_opposite_pairing(vectors: list[np.ndarray]) -> tuple[tuple[int, int], ...]:
    """Pair six ligand directions into three maximally opposite pairs.

    The score minimizes sum(1 + cos(theta)); exact opposite vectors score zero.
    """
    n = len(vectors)
    if n != 6:
        raise ValueError("Automatic octahedral opposite-pairing currently requires 6 ligands")
    u = [_unit(v) for v in vectors]

    def rec(remaining: tuple[int, ...]):
        if not remaining:
            yield ()
            return
        i = remaining[0]
        for pos in range(1, len(remaining)):
            j = remaining[pos]
            rest = remaining[1:pos] + remaining[pos + 1 :]
            for tail in rec(rest):
                yield ((i, j),) + tail

    best = None
    best_score = np.inf
    for pairing in rec(tuple(range(n))):
        score = sum(1.0 + float(np.dot(u[i], u[j])) for i, j in pairing)
        if score < best_score:
            best_score = score
            best = pairing
    assert best is not None
    return tuple(best)


def _axis_from_pair(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    # Opposite ligands define an unoriented line. Difference is robust even if
    # center-ligand distances are unequal.
    return _canonical_sign(_unit(_unit(v1) - _unit(v2)))


def _automatic_quality_metrics(
    vectors: list[np.ndarray],
    pairing_local: tuple[tuple[int, int], ...],
    axes: list[np.ndarray],
    mean_dist: list[float],
    z_idx: int,
) -> dict[str, object]:
    unit_vectors = [_unit(v) for v in vectors]
    pair_cosines = [
        float(np.dot(unit_vectors[i], unit_vectors[j])) for i, j in pairing_local
    ]
    opposition_errors = [1.0 + value for value in pair_cosines]

    raw_axis_dot_products = []
    for i in range(3):
        for j in range(i + 1, 3):
            raw_axis_dot_products.append(float(np.dot(axes[i], axes[j])))

    competing = [float(mean_dist[k]) for k in range(3) if k != z_idx]
    nearest_competing_distance = min(
        competing,
        key=lambda value: abs(value - float(mean_dist[z_idx])),
    )
    scale = max(float(np.mean(mean_dist)), 1e-12)
    z_gap_fraction = abs(float(mean_dist[z_idx]) - nearest_competing_distance) / scale

    return {
        "pair_cosines": pair_cosines,
        "pair_opposition_errors": opposition_errors,
        "max_opposition_error": float(max(opposition_errors)),
        "raw_axis_dot_products": raw_axis_dot_products,
        "raw_axis_orthogonality_error": float(
            max(abs(value) for value in raw_axis_dot_products)
        ),
        "selected_z_pair_local_index": int(z_idx),
        "selected_z_mean_distance": float(mean_dist[z_idx]),
        "nearest_competing_mean_distance": float(nearest_competing_distance),
        "z_gap_fraction": float(z_gap_fraction),
    }


def build_local_frame(
    structure: Structure,
    center_index: int,
    ligand: str | None = None,
    coordination: int = 6,
    cutoff: float | None = None,
    z_policy: str = "longest",
) -> LocalFrame:
    """Construct a right-handed local orbital frame around one site.

    For coordination=6, opposite ligand directions are paired geometrically.
    The pair with the largest mean bond length is local z by default, which is
    useful for tetragonally elongated octahedra. The remaining pair most nearly
    orthogonal to z seeds x; Gram-Schmidt and y=z×x enforce orthonormality.

    The frame metadata also records pre-orthogonalization geometry diagnostics so
    users can detect ambiguous or strongly distorted automatic assignments.
    """
    rows = _candidate_ligands(structure, center_index, ligand, coordination, cutoff)
    if coordination != 6:
        raise NotImplementedError(
            "Initial release implements automatic frame construction for coordination=6."
        )

    vectors = [r[1] for r in rows]
    pairing_local = _best_opposite_pairing(vectors)

    axes = []
    mean_dist = []
    for i, j in pairing_local:
        axes.append(_axis_from_pair(vectors[i], vectors[j]))
        mean_dist.append(0.5 * (rows[i][2] + rows[j][2]))

    if z_policy == "longest":
        z_idx = int(np.argmax(mean_dist))
    elif z_policy == "shortest":
        z_idx = int(np.argmin(mean_dist))
    else:
        raise ValueError("z_policy must be 'longest' or 'shortest'")

    quality = _automatic_quality_metrics(vectors, pairing_local, axes, mean_dist, z_idx)

    z = _unit(axes[z_idx])
    remaining = [k for k in range(3) if k != z_idx]
    # Choose the line closest to 90 degrees from z as the x seed.
    x_idx = min(remaining, key=lambda k: abs(float(np.dot(axes[k], z))))
    x0 = axes[x_idx] - float(np.dot(axes[x_idx], z)) * z
    x = _canonical_sign(_unit(x0))
    y = _unit(np.cross(z, x))
    # Recompute x to eliminate numerical non-orthogonality while preserving RH frame.
    x = _unit(np.cross(y, z))

    R = np.column_stack((x, y, z))
    if np.linalg.det(R) < 0:
        y = -y
        R = np.column_stack((x, y, z))

    ligand_indices = tuple(int(r[0]) for r in rows)
    ligand_symbols = tuple(r[3] for r in rows)
    global_pairs = tuple((ligand_indices[i], ligand_indices[j]) for i, j in pairing_local)

    return LocalFrame(
        site_index=center_index,
        center_symbol=structure[center_index].specie.symbol,
        ligand_indices=ligand_indices,
        ligand_symbols=ligand_symbols,
        x=x,
        y=y,
        z=z,
        rotation_local_to_global=R,
        pair_indices=global_pairs,
        pair_mean_distances=tuple(float(v) for v in mean_dist),
        provider="auto-octahedral",
        mode="full-3d",
        metadata={
            "coordination": int(coordination),
            "ligand_filter": ligand,
            "cutoff": cutoff,
            "z_policy": z_policy,
            "quality": quality,
        },
    )


def build_frames_for_element(
    structure: Structure,
    center: str,
    ligand: str | None = None,
    coordination: int = 6,
    cutoff: float | None = None,
    z_policy: str = "longest",
    sites: Iterable[int] | None = None,
) -> list[LocalFrame]:
    indices = list(sites) if sites is not None else find_center_sites(structure, center)
    return [
        build_local_frame(
            structure,
            i,
            ligand=ligand,
            coordination=coordination,
            cutoff=cutoff,
            z_policy=z_policy,
        )
        for i in indices
    ]
