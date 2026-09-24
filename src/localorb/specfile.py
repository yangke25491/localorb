from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
from pymatgen.core import Structure

from .frames import LocalFrame, build_local_frame
from .manual import build_manual_frame, resolve_site_selector
from .vectors import build_vector_frame


def _as_vector(value: Any, name: str) -> np.ndarray:
    if isinstance(value, str):
        tokens = value.replace(",", " ").split()
        if len(tokens) != 3:
            raise ValueError(f"{name} must have three components")
        return np.array([float(token) for token in tokens], dtype=float)
    array = np.asarray(value, dtype=float)
    if array.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {array.shape}")
    return array


def _image_text(value: Any | None, name: str) -> str | None:
    if value is None:
        return None
    if isinstance(value, str):
        return value
    array = np.asarray(value)
    if array.shape != (3,):
        raise ValueError(f"{name} must have three integer components")
    try:
        integers = [int(item) for item in array]
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain integers") from exc
    if not np.allclose(np.asarray(array, dtype=float), integers, atol=0.0):
        raise ValueError(f"{name} must contain integers")
    return ",".join(str(value) for value in integers)


def _selector_with_image(
    selector: str,
    *,
    specific_image: Any | None,
    shared_image: Any | None,
) -> str:
    if "@" in selector:
        return selector
    image_value = specific_image if specific_image is not None else shared_image
    image = _image_text(image_value, "periodic image")
    return selector if image is None else f"{selector}@{image}"


def _merged_defaults(defaults: dict[str, Any], entry: dict[str, Any]) -> dict[str, Any]:
    merged = dict(defaults)
    merged.update(entry)
    return merged


def _build_manual_from_entry(
    structure: Structure,
    entry: dict[str, Any],
) -> LocalFrame:
    required = ["center_atom", "x_atom", "plane_atom"]
    missing = [name for name in required if not entry.get(name)]
    if missing:
        raise ValueError("manual frame entry missing: " + ", ".join(missing))

    shared_image = entry.get("image")
    center = _selector_with_image(
        str(entry["center_atom"]),
        specific_image=entry.get("center_image"),
        shared_image=shared_image,
    )
    x_atom = _selector_with_image(
        str(entry["x_atom"]),
        specific_image=entry.get("x_image"),
        shared_image=shared_image,
    )
    plane_atom = _selector_with_image(
        str(entry["plane_atom"]),
        specific_image=entry.get("plane_image", entry.get("y_image")),
        shared_image=shared_image,
    )

    return build_manual_frame(
        structure,
        center=center,
        x_atom=x_atom,
        plane_atom=plane_atom,
        mode=str(entry.get("mode", entry.get("manual_mode", "full-3d"))),
        nearest_image=bool(entry.get("nearest_image", True)),
        index_base=int(entry.get("index_base", 1)),
        cartesian_frame=str(entry.get("cart_frame", entry.get("cartesian_frame", "poscar"))),
    )


def _build_vectors_from_entry(
    structure: Structure,
    entry: dict[str, Any],
) -> LocalFrame:
    if not entry.get("center_atom"):
        raise ValueError("vectors frame entry requires center_atom")
    if "x_vector" not in entry or "z_vector" not in entry:
        raise ValueError("vectors frame entry requires x_vector and z_vector")

    selection = resolve_site_selector(
        structure,
        str(entry["center_atom"]),
        index_base=int(entry.get("index_base", 1)),
    )
    return build_vector_frame(
        structure,
        center_index=selection.site_index,
        x_vector=_as_vector(entry["x_vector"], "x_vector"),
        z_vector=_as_vector(entry["z_vector"], "z_vector"),
    )


def _build_auto_from_entry(
    structure: Structure,
    entry: dict[str, Any],
) -> LocalFrame:
    if "site" in entry:
        site_index = int(entry["site"])
    elif entry.get("center_atom"):
        selection = resolve_site_selector(
            structure,
            str(entry["center_atom"]),
            index_base=int(entry.get("index_base", 1)),
        )
        site_index = selection.site_index
    else:
        raise ValueError("auto frame entry requires site or center_atom")

    cutoff = entry.get("cutoff")
    return build_local_frame(
        structure,
        site_index,
        ligand=entry.get("ligand"),
        coordination=int(entry.get("coordination", 6)),
        cutoff=None if cutoff is None else float(cutoff),
        z_policy=str(entry.get("z_policy", "longest")),
    )


def build_frame_from_spec_entry(
    structure: Structure,
    entry: dict[str, Any],
) -> LocalFrame:
    provider = str(entry.get("provider", "manual")).lower()
    if provider in {"manual", "manual-atoms"}:
        return _build_manual_from_entry(structure, entry)
    if provider in {"vectors", "explicit-vectors"}:
        return _build_vectors_from_entry(structure, entry)
    if provider in {"auto", "auto-octahedral"}:
        return _build_auto_from_entry(structure, entry)
    raise ValueError(f"Unknown frame provider in specification: {provider!r}")


def load_frame_spec_file(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON frame specification {path}: {exc}") from exc

    if not isinstance(payload, dict):
        raise ValueError("Frame specification root must be a JSON object")
    if "frames" not in payload:
        raise ValueError("Frame specification must contain a 'frames' array")
    if not isinstance(payload["frames"], list) or not payload["frames"]:
        raise ValueError("Frame specification 'frames' must be a non-empty array")
    if "defaults" in payload and not isinstance(payload["defaults"], dict):
        raise ValueError("Frame specification 'defaults' must be an object")
    return payload


def build_frames_from_spec_file(
    structure: Structure,
    path: str | Path,
) -> list[LocalFrame]:
    payload = load_frame_spec_file(path)
    defaults = dict(payload.get("defaults", {}))
    frames: list[LocalFrame] = []
    seen_sites: set[int] = set()

    for position, raw_entry in enumerate(payload["frames"]):
        if not isinstance(raw_entry, dict):
            raise ValueError(f"frames[{position}] must be a JSON object")
        entry = _merged_defaults(defaults, raw_entry)
        frame = build_frame_from_spec_entry(structure, entry)
        if frame.site_index in seen_sites:
            raise ValueError(
                f"Frame specification defines site {frame.site_index} more than once; "
                "use exactly one orbital frame per center site"
            )
        seen_sites.add(frame.site_index)
        frames.append(frame)

    return frames
