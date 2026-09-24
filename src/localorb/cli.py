from __future__ import annotations

import argparse
import json
from pathlib import Path

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar

from .frames import build_frames_for_element, build_local_frame
from .rotate import frame_alignment_error, rotate_structure_to_frame
from .wannier import render_projection_block


def _common_frame_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("structure", help="Input structure readable by pymatgen, e.g. POSCAR")
    parser.add_argument("--center", required=True, help="Center element, e.g. Ni")
    parser.add_argument("--ligand", default=None, help="Ligand element filter, e.g. O")
    parser.add_argument("--coordination", type=int, default=6, help="Number of nearest ligands (default: 6)")
    parser.add_argument("--cutoff", type=float, default=None, help="Neighbor search cutoff in Å")
    parser.add_argument("--z-policy", choices=["longest", "shortest"], default="longest",
                        help="Which opposite ligand pair defines local z")
    parser.add_argument("--sites", default=None,
                        help="Comma-separated 0-based site indices; default: all atoms matching --center")


def _parse_sites(text: str | None):
    if text is None:
        return None
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def _frames_from_args(args):
    structure = Structure.from_file(args.structure)
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


def cmd_inspect(args) -> int:
    _, frames = _frames_from_args(args)
    for f in frames:
        print(f"site {f.site_index:4d} {f.center_symbol}")
        print("  ligands:", ", ".join(map(str, f.ligand_indices)))
        print("  x =", " ".join(f"{v:+.8f}" for v in f.x))
        print("  y =", " ".join(f"{v:+.8f}" for v in f.y))
        print("  z =", " ".join(f"{v:+.8f}" for v in f.z))
        print(f"  det(R) = {float(__import__('numpy').linalg.det(f.rotation_local_to_global)):.10f}")
    return 0


def cmd_report(args) -> int:
    _, frames = _frames_from_args(args)
    payload = {
        "convention": {
            "rotation_local_to_global": "columns are local x,y,z unit vectors in global Cartesian coordinates",
            "indices": "0-based pymatgen/VASP atom ordering",
        },
        "frames": [f.as_dict() for f in frames],
    }
    Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"Wrote {args.output}")
    return 0


def cmd_wannier(args) -> int:
    structure, frames = _frames_from_args(args)
    orbitals = [x.strip() for x in args.orbitals.split(",") if x.strip()]
    text = render_projection_block(structure, frames, orbitals)
    Path(args.output).write_text(text, encoding="utf-8")
    print(f"Wrote {args.output} with {len(frames)} site-dependent local frames")
    return 0


def cmd_rotate_poscar(args) -> int:
    structure = Structure.from_file(args.structure)
    frame = build_local_frame(
        structure,
        args.site,
        ligand=args.ligand,
        coordination=args.coordination,
        cutoff=args.cutoff,
        z_policy=args.z_policy,
    )
    rotated = rotate_structure_to_frame(structure, frame)
    Poscar(rotated).write_file(args.output)
    print(f"Wrote {args.output}")
    print(f"Frame alignment error: {frame_alignment_error(frame):.3e}")
    return 0


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
    p.add_argument("--orbitals", default="dz2,dx2-y2",
                   help="Comma-separated orbitals (default: dz2,dx2-y2)")
    p.add_argument("-o", "--output", default="projections.win")
    p.set_defaults(func=cmd_wannier)

    p = sub.add_parser("rotate-poscar", help="Rigidly rotate a structure to one site's local frame")
    p.add_argument("structure")
    p.add_argument("--site", type=int, required=True, help="0-based center site index")
    p.add_argument("--ligand", default=None)
    p.add_argument("--coordination", type=int, default=6)
    p.add_argument("--cutoff", type=float, default=None)
    p.add_argument("--z-policy", choices=["longest", "shortest"], default="longest")
    p.add_argument("-o", "--output", default="POSCAR.rotated")
    p.set_defaults(func=cmd_rotate_poscar)

    return parser


def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
