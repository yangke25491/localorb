from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from .frames import build_frames_for_element, build_local_frame
from .manual import build_manual_frame
from .orbitals import D_ORBITALS, d_rotation_matrix
from .procar import save_rotated_procar_npz
from .rotate import frame_alignment_error, rotate_structure_to_frame
from .wannier import render_projection_block


def _add_manual_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--center-atom", help="Manual center selector, e.g. Ni2, 6, or Ni2@0,0,1")
    parser.add_argument("--x-atom", help="Manual atom defining local +x, e.g. O4 or 12@0,0,1")
    parser.add_argument("--plane-atom", help="Manual second atom defining the local xy plane, e.g. O7")
    parser.add_argument(
        "--manual-mode",
        choices=["fixed-z", "full-3d"],
        default="full-3d",
        help="Manual frame convention: keep Cartesian z fixed, or construct a full 3D frame",
    )
    parser.add_argument(
        "--no-nearest-image",
        dest="nearest_image",
        action="store_false",
        default=True,
        help="Do not move manual atoms without @i,j,k to their nearest periodic image",
    )
    parser.add_argument(
        "--index-base",
        type=int,
        choices=[0, 1],
        default=1,
        help="Index base for bare manual integer selectors; default 1 matches POSCAR/VESTA",
    )


def _common_frame_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", help="Input structure readable by pymatgen, e.g. POSCAR")
    parser.add_argument(
        "--provider",
        choices=["auto", "manual"],
        default="auto",
        help="How to construct local frames: automatic octahedral geometry or explicit atoms",
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


def _manual_frame_from_args(structure: Structure, args):
    _require_manual_args(args)
    return build_manual_frame(
        structure,
        center=args.center_atom,
        x_atom=args.x_atom,
        plane_atom=args.plane_atom,
        mode=args.manual_mode,
        nearest_image=args.nearest_image,
        index_base=args.index_base,
    )


def _frames_from_args(args):
    structure = Structure.from_file(args.structure)
    if args.provider == "manual":
        return structure, [_manual_frame_from_args(structure, args)]

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
            "rotation_local_to_global": "columns are local x,y,z unit vectors in global Cartesian coordinates",
            "rotation_global_to_local": "transpose/inverse of rotation_local_to_global",
            "auto_indices": "0-based pymatgen/VASP ordering",
            "manual_integer_selectors": f"index_base={args.index_base}",
        },
        "frames": [frame.as_dict() for frame in frames],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


def cmd_wannier(args) -> int:
    structure, frames = _frames_from_args(args)
    orbitals = [x.strip() for x in args.orbitals.split(",") if x.strip()]
    text = render_projection_block(structure, frames, orbitals)
    Path(args.output).write_text(text, encoding="utf-8")
    print(f"Wrote {args.output} with {len(frames)} site-dependent local frame(s)")
    return 0


def cmd_rotate_poscar(args) -> int:
    structure = Structure.from_file(args.structure)
    frame = _single_frame_from_args(structure, args)
    rotated = rotate_structure_to_frame(structure, frame)
    # Poscar writes a consistent structure from fractional coordinates, so both
    # Direct and Cartesian input structures are safely handled.
    Poscar(rotated).write_file(args.output, direct=True)
    print(f"Wrote {args.output}")
    print(f"Provider: {frame.provider}; mode: {frame.mode}")
    print(f"Frame alignment error: {frame_alignment_error(frame):.3e}")
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
    parser.add_argument("--provider", choices=["auto", "manual"], default="auto")
    parser.add_argument("--site", type=int, default=None, help="Auto provider: 0-based center site index")
    parser.add_argument("--ligand", default=None)
    parser.add_argument("--coordination", type=int, default=6)
    parser.add_argument("--cutoff", type=float, default=None)
    parser.add_argument("--z-policy", choices=["longest", "shortest"], default="longest")
    _add_manual_args(parser)


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
