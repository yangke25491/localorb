# localorb

`localorb` builds physically meaningful local orbital frames directly from crystal structures and connects those frames to Wannier90, VASP-oriented workflows, and future tight-binding/post-processing backends.

The central idea is simple: a global Cartesian frame is often not the natural orbital frame of a transition-metal coordination polyhedron. Instead of forcing one correction method onto every calculation, `localorb` provides multiple strategies and lets the user choose the one appropriate to the structure and workflow.

## What localorb does

- Detect local coordination environments around selected center atoms.
- Construct a right-handed local orthonormal frame `(x', y', z')` from ligand geometry.
- Generate site-dependent Wannier90 projections such as local `dz2` and `dx2-y2` orbitals.
- Rotate an entire POSCAR when a single common local frame is appropriate.
- Export reproducible frame reports containing local axes and rotation matrices.
- Provide a backend-oriented architecture for future PROCAR, Hamiltonian, and orbital-rotation tools.

## Why this exists

VASP orbital labels such as `dxy`, `dz2`, and `dx2-y2` are defined relative to a global Cartesian frame. In tilted, rotated, distorted, supercell, or disordered coordination environments, those labels may not coincide with the physically meaningful local crystal-field orbitals. `localorb` addresses this at the structure level and can feed the resulting local frames into the representation most appropriate for the calculation.

## Installation

```bash
python -m pip install -e .
```

Requirements are intentionally minimal: `numpy` and `pymatgen`.

## Quick start

Inspect all Ni-centered oxygen octahedra in a POSCAR:

```bash
localorb inspect POSCAR --center Ni --ligand O --coordination 6
```

Write a JSON report:

```bash
localorb report POSCAR --center Ni --ligand O --coordination 6 -o local_frames.json
```

Generate Wannier90 local `eg` projections:

```bash
localorb wannier POSCAR --center Ni --ligand O --coordination 6 \
  --orbitals dz2,dx2-y2 -o projections.win
```

Generate all five local d orbitals:

```bash
localorb wannier POSCAR --center Ni --ligand O --coordination 6 \
  --orbitals dxy,dyz,dz2,dxz,dx2-y2 -o projections.win
```

Rotate the entire POSCAR so the local frame of a chosen site becomes the global Cartesian frame:

```bash
localorb rotate-poscar POSCAR --site 10 --ligand O --coordination 6 \
  -o POSCAR.rotated
```

This rigid rotation is appropriate only when one common frame is meaningful. For structures with site-dependent octahedral tilts, use site-dependent Wannier projections instead.

## Local-frame convention

For an approximately octahedral environment, `localorb` pairs opposite ligand directions and chooses three approximate axes. The axis associated with the most axial pair (by default the pair with the largest mean center-ligand distance) is taken as local `z`. A second pair defines local `x`; `y = z x x`, followed by orthonormalization. Axis signs are canonicalized for reproducibility.

The local-to-global frame matrix is stored as

```text
R = [ x'  y'  z' ]
```

with the three local unit vectors as columns in global Cartesian coordinates. Thus a local coordinate vector transforms as

```text
v_global = R @ v_local
```

and

```text
v_local = R.T @ v_global
```

## Wannier90 strategy

Wannier90 allows projection functions to carry their own local `z` and `x` axes. `localorb wannier` therefore writes one projection per selected site using that site's local frame, rather than rotating the whole structure when different sites require different frames.

## Roadmap

1. Geometry + Wannier90 backend (implemented in the initial release).
2. Optional whole-structure POSCAR rotation (implemented in the initial release).
3. VASP PROCAR complex-coefficient rotation backend.
4. Rotation of Wannier/TB Hamiltonians between global and local orbital bases.
5. Additional coordination templates: square-planar, tetrahedral, trigonal-prismatic, user-defined ligand groups.
6. Validation utilities for comparing orbital-resolved bands before/after basis changes.

## Scientific caution

Automatic local-axis construction is a geometric convention, not a substitute for physical judgment. For strongly distorted or ambiguous environments, inspect the reported ligand assignment and use explicit site/axis overrides when necessary.

## License

MIT
