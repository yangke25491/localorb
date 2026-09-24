from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from .diagnostics import diagnose_frame
from .frames import build_frames_for_element, build_local_frame
from .manual import build_manual_frame, resolve_site_selector
from .orbitals import D_ORBITALS, d_rotation_matrix
from .procar import save_rotated_procar_npz
from .rotate import frame_alignment_error, rotate_structure_to_frame
from .vectors import build_vector_frame, parse_vector
from .wannier import render_projection_block


def _add_manual_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--center-atom", help="Manual/vector center selector, e.g. Ni2, 6, or Ni2@0,0,1")
    parser.add_argument("--x-atom", help="Manual atom defining local +x, e.g. O4 or 12@0,0,1")
    parser.add_argument("--plane-atom", help="Manual second atom defining the local xy plane, e.g. O7")
    parser.add_argument(
        "--manual-mode",
        choices=["fixed-z", "full-3d"],
        default="full-3d",
        help="Manual frame convention: keep working-frame z fixed, or construct a full 3D frame",
    )
    parser.add_argument(
        "--cart-frame",
        choices=["poscar", "auto", "vesta"],
        default="poscar",
        help=(
            "Working Cartesian frame for manual atoms. 'poscar' is the safe default; "
            "'auto'/'vesta' enable the narrow legacy rhombohedral VESTA compatibility layer."
        ),
    )
    parser.add_argument(
        "--image",
        default=None,
        help="Shared periodic image for manual selectors, e.g. 0,0,1; explicit ATOM@i,j,k wins",
    )
    parser.add_argument("--center-image", default=None, help="Periodic image override for --center-atom")
    parser.add_argument("--x-image", default=None, help="Periodic image override for --x-atom")
    parser.add_argument(
        "--plane-image",
        "--y-image",
        dest="plane_image",
        default=None,
        help="Periodic image override for --plane-atom (alias: --y-image)",
    )
    parser.add_argument(
        "--no-nearest-image",
        dest="nearest_image",
        action="store_false",
        default=True,
        help="Do not move manual atoms without an explicit image to their nearest periodic image",
    )
    parser.add_argument(
        "--index-base",
        type=int,
        choices=[0, 1],
        default=1,
        help="Index base for bare manual integer selectors; default 1 matches POSCAR/VESTA",
    )


def _add_vector_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--x-vector", help="Explicit local x direction in POSCAR Cartesian coordinates, e.g. 1,1,0")
    parser.add_argument("--z-vector", help="Explicit local z direction in POSCAR Cartesian coordinates, e.g. 0,0,1")


def _common_frame_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", help="Input structure readable by pymatgen, e.g. POSCAR")
    parser.add_argument(
        "--provider",
        choices=["auto", "manual", "vectors"],
        default="auto",
        help="Frame source: automatic octahedral geometry, explicit atoms, or explicit vectors",
    )
    parser.add_argument("--center", help="Auto provider: center element, e.g. Ni")
    parser.add_argument("--ligand", default=None, help="Auto provider: ligand element filter, e.g. O")
    parser.add_argument("--coordination", type=int, default=6, help="Auto provider: nearest ligands (default: 6)")
    parser.add_argument("--cutoff", type=float, default=None, help="Auto provider: neighbor search cutoff in Å")
    parser.add_argument(
        "--z-policy",
        choices=["longest", "shortest"],
        default="longest",
        help="Auto provider: which opposite ligand pair defines local z",
    )
    parser.add_argument(
        "--sites",
        default=None,
        help="Auto provider: comma-separated 0-based site indices; default all matching --center",
    )
    _add_manual_args(parser)
    _add_vector_args(parser)


def _parse_sites(text: str | None):
    if text is None:
        return None
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def _require_manual_args(args) -> None:
    missing = [
        name
        for name, value in [
            ("--center-atom", args.center_atom),
            ("--x-atom", args.x_atom),
            ("--plane-atom", args.plane_atom),
        ]
        if not value
    ]
    if missing:
        raise SystemExit("Manual provider requires " + ", ".join(missing))


def _selector_with_image(selector: str, specific_image: str | None, shared_image: str | None) -> str:
    if "@" in selector:
        return selector
    image = specific_image if specific_image is not None else shared_image
    return selector if image is None else f"{selector}@{image}"


def _manual_frame_from_args(structure: Structure, args, mode: str | None = None):
    _require_manual_args(args)
    center = _selector_with_image(args.center_atom, args.center_image, args.image)
    x_atom = _selector_with_image(args.x_atom, args.x_image, args.image)
    plane_atom = _selector_with_image(args.plane_atom, args.plane_image, args.image)
    return build_manual_frame(
        structure,
        center=center,
        x_atom=x_atom,
        plane_atom=plane_atom,
        mode=mode or args.manual_mode,
        nearest_image=args.nearest_image,
        index_base=args.index_base,
        cartesian_frame=args.cart_frame,
    )


def _vector_frame_from_args(structure: Structure, args):
    if not args.center_atom or not args.x_vector or not args.z_vector:
        raise SystemExit("Vector provider requires --center-atom, --x-vector, and --z-vector")
    selection = resolve_site_selector(structure, args.center_atom, index_base=args.index_base)
    return build_vector_frame(
        structure,
        center_index=selection.site_index,
        x_vector=parse_vector(args.x_vector),
        z_vector=parse_vector(args.z_vector),
    )


def _frames_from_args(args):
    structure = Structure.from_file(args.structure)
    if args.provider == "manual":
        return structure, [_manual_frame_from_args(structure, args)]
    if args.provider == "vectors":
        return structure, [_vector_frame_from_args(structure, args)]

    if not args.center:
        raise SystemExit("Auto provider requires --center, e.g. --center Ni")
    frames = build_frames_for_element(
        structure,
        center=args.center,
        ligand=args.ligand,
        coordination=args.coordination,
        cutoff=args.cutoff,
        z_policy=args.z_policy,
        sites=_parse_sites(args.sites),
    )
    if not frames:
        raise SystemExit(f"No {args.center} sites found")
    return structure, frames


def _single_frame_from_args(structure: Structure, args):
    if args.provider == "manual":
        return _manual_frame_from_args(structure, args)
    if args.provider == "vectors":
        return _vector_frame_from_args(structure, args)
    if args.site is None:
        raise SystemExit("Auto provider requires --site for this command")
    return build_local_frame(
        structure,
        args.site,
        ligand=args.ligand,
        coordination=args.coordination,
        cutoff=args.cutoff,
        z_policy=args.z_policy,
    )


def _print_frame(frame) -> None:
    print(f"site {frame.site_index:4d} {frame.center_symbol}")
    print(f"  provider = {frame.provider}")
    print(f"  mode     = {frame.mode}")
    if frame.ligand_indices:
        print("  ligands  =", ", ".join(map(str, frame.ligand_indices)))
    print("  x =", " ".join(f"{v:+.8f}" for v in frame.x))
    print("  y =", " ".join(f"{v:+.8f}" for v in frame.y))
    print("  z =", " ".join(f"{v:+.8f}" for v in frame.z))
    print(f"  det(R) = {float(np.linalg.det(frame.rotation_local_to_global)):.10f}")
    quality = frame.metadata.get("quality") if isinstance(frame.metadata, dict) else None
    if isinstance(quality, dict):
        print(
            "  auto quality: "
            f"opposition={float(quality['max_opposition_error']):.4f}, "
            f"raw-axis-dot={float(quality['raw_axis_orthogonality_error']):.4f}, "
            f"z-gap={float(quality['z_gap_fraction']):.4f}"
        )
    if frame.metadata:
        print("  metadata =", json.dumps(frame.metadata, ensure_ascii=False))


def cmd_inspect(args) -> int:
    _, frames = _frames_from_args(args)
    for frame in frames:
        _print_frame(frame)
    return 0


def cmd_report(args) -> int:
    _, frames = _frames_from_args(args)
    payload = {
        "convention": {
            "rotation_local_to_global": "columns are local x,y,z unit vectors in POSCAR Cartesian coordinates",
            "rotation_global_to_local": "transpose/inverse of rotation_local_to_global",
            "auto_indices": "0-based pymatgen/VASP ordering",
            "manual_integer_selectors": f"index_base={args.index_base}",
        },
        "frames": [frame.as_dict() for frame in frames],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


def cmd_validate(args) -> int:
    _, frames = _frames_from_args(args)
    diagnostics = [
        diagnose_frame(
            frame,
            max_opposition_error_warn=args.max_opposition_error,
            max_raw_axis_dot_warn=args.max_raw_axis_dot,
            min_z_gap_fraction_warn=args.min_z_gap,
        )
        for frame in frames
    ]

    for frame, diag in zip(frames, diagnostics):
        print(
            f"[{diag.status}] site {frame.site_index} {frame.center_symbol} "
            f"provider={frame.provider} mode={frame.mode} "
            f"orth_err={diag.orthonormal_error:.3e} det={diag.determinant:.10f}"
        )
        for name, value in diag.metrics.items():
            print(f"  {name} = {value:.6g}")
        for message in diag.warnings:
            print(f"  WARNING: {message}")
        for message in diag.errors:
            print(f"  ERROR: {message}")

    if args.output:
        payload = {
            "thresholds": {
                "max_opposition_error_warn": args.max_opposition_error,
                "max_raw_axis_dot_warn": args.max_raw_axis_dot,
                "min_z_gap_fraction_warn": args.min_z_gap,
            },
            "diagnostics": [diag.as_dict() for diag in diagnostics],
        }
        Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(f"Wrote {args.output}")

    if any(diag.status == "FAIL" for diag in diagnostics):
        return 1
    if args.strict and any(diag.status == "WARN" for diag in diagnostics):
        return 2
    return 0


def cmd_wannier(args) -> int:
    structure, frames = _frames_from_args(args)
    orbitals = [x.strip() for x in args.orbitals.split(",") if x.strip()]
    text = render_projection_block(structure, frames, orbitals)
    Path(args.output).write_text(text, encoding="utf-8")
    print(f"Wrote {args.output} with {len(frames)} site-dependent local frame(s)")
    return 0


def _write_rotated_poscar(structure: Structure, frame, output: Path) -> None:
    rotated = rotate_structure_to_frame(structure, frame)
    # Always write a self-consistent Direct-coordinate POSCAR. This deliberately
    # avoids the legacy bug where Cartesian input could rotate the lattice while
    # leaving atomic Cartesian coordinates unchanged.
    Poscar(rotated).write_file(str(output), direct=True)
    print(f"Wrote {output}")
    print(f"Provider: {frame.provider}; mode: {frame.mode}")
    print(f"Frame alignment error: {frame_alignment_error(frame):.3e}")


def _mode_output_path(path: str, mode: str) -> Path:
    p = Path(path)
    suffix = mode.replace("-", "_")
    if p.suffix:
        return p.with_name(f"{p.stem}_{suffix}{p.suffix}")
    return p.with_name(f"{p.name}_{suffix}")


def cmd_rotate_poscar(args) -> int:
    structure = Structure.from_file(args.structure)
    if args.provider == "manual" and args.both_manual_modes:
        for mode in ("fixed-z", "full-3d"):
            frame = _manual_frame_from_args(structure, args, mode=mode)
            _write_rotated_poscar(structure, frame, _mode_output_path(args.output, mode))
        return 0

    frame = _single_frame_from_args(structure, args)
    _write_rotated_poscar(structure, frame, Path(args.output))
    return 0


def cmd_d_matrix(args) -> int:
    structure = Structure.from_file(args.structure)
    frame = _single_frame_from_args(structure, args)
    T = d_rotation_matrix(frame.rotation_local_to_global)
    print("orbital order:", " ".join(D_ORBITALS))
    print(f"provider: {frame.provider}; mode: {frame.mode}")
    np.set_printoptions(precision=10, suppress=True)
    print(T)
    return 0


def cmd_procar(args) -> int:
    _, frames = _frames_from_args(args)
    save_rotated_procar_npz(args.output, args.procar, frames)
    print(f"Wrote {args.output}")
    print("arrays: coefficients[spin,k,band,site,orbital], weights, site_indices, orbital_names, spin_values")
    print("orbital order:", " ".join(D_ORBITALS))
    return 0


def _add_single_frame_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure")
    parser.add_argument("--provider", choices=["auto", "manual", "vectors"], default="auto")
    parser.add_argument("--site", type=int, default=None, help="Auto provider: 0-based center site index")
    parser.add_argument("--ligand", default=None)
    parser.add_argument("--coordination", type=int, default=6)
    parser.add_argument("--cutoff", type=float, default=None)
    parser.add_argument("--z-policy", choices=["longest", "shortest"], default="longest")
    _add_manual_args(parser)
    _add_vector_args(parser)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="localorb",
        description="Build local orbital frames from crystal structures.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("inspect", help="Print local frames")
    _common_frame_args(p)
    p.set_defaults(func=cmd_inspect)

    p = sub.add_parser("report", help="Write local frames and rotation matrices to JSON")
    _common_frame_args(p)
    p.add_argument("-o", "--output", default="local_frames.json")
    p.set_defaults(func=cmd_report)

    p = sub.add_parser("validate", help="Validate rotation matrices and flag ambiguous automatic frames")
    _common_frame_args(p)
    p.add_argument(
        "--max-opposition-error",
        type=float,
        default=0.10,
        help="Warn when max(1+cos(theta)) of opposite pairs exceeds this value",
    )
    p.add_argument(
        "--max-raw-axis-dot",
        type=float,
        default=0.15,
        help="Warn when pre-orthogonalization max |axis_i dot axis_j| exceeds this value",
    )
    p.add_argument(
        "--min-z-gap",
        type=float,
        default=0.02,
        help="Warn when automatic z bond-length separation is below this relative gap",
    )
    p.add_argument("--strict", action="store_true", help="Return exit code 2 when warnings are present")
    p.add_argument("-o", "--output", default=None, help="Optional JSON diagnostics output")
    p.set_defaults(func=cmd_validate)

    p = sub.add_parser("wannier", help="Generate a Wannier90 projections block")
    _common_frame_args(p)
    p.add_argument(
        "--orbitals",
        default="dz2,dx2-y2",
        help="Comma-separated orbitals (default: dz2,dx2-y2)",
    )
    p.add_argument("-o", "--output", default="projections.win")
    p.set_defaults(func=cmd_wannier)

    p = sub.add_parser("rotate-poscar", help="Rigidly rotate a structure to one local frame")
    _add_single_frame_args(p)
    p.add_argument("--both-manual-modes", action="store_true", help="For provider=manual, write fixed-z and full-3d outputs")
    p.add_argument("-o", "--output", default="POSCAR.rotated")
    p.set_defaults(func=cmd_rotate_poscar)

    p = sub.add_parser("d-matrix", help="Print the 5x5 real-d orbital rotation matrix for one frame")
    _add_single_frame_args(p)
    p.set_defaults(func=cmd_d_matrix)

    p = sub.add_parser("procar", help="Rotate phase-resolved PROCAR d amplitudes into local frames")
    _common_frame_args(p)
    p.add_argument("--procar", default="PROCAR", help="Path to phase-resolved PROCAR")
    p.add_argument("-o", "--output", default="PROCAR_LOCAL.npz")
    p.set_defaults(func=cmd_procar)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
