# localorb

`localorb` builds physically meaningful local orbital frames directly from crystal structures and connects those frames to Wannier90, VASP-oriented workflows, and orbital-basis transformations.

The central idea is simple: a global Cartesian frame is often not the natural orbital frame of a transition-metal coordination polyhedron. Instead of forcing one correction method onto every calculation, `localorb` provides multiple strategies and lets the user choose the one appropriate to the structure and workflow.

## Implemented workflows

1. **Structure → local frames**: identify ligand environments and construct a right-handed `(x', y', z')` frame for each selected site.
2. **Structure → Wannier90**: generate site-dependent local projections such as `dz2` and `dx2-y2` using Wannier90's native `z=` and `x=` projection axes.
3. **Structure → rotated POSCAR**: rigidly rotate the whole structure when a single common local frame is physically appropriate.
4. **Structure → 5×5 d-orbital matrix**: construct the real-d representation of the 3D frame rotation without depending on RotSph.
5. **POSCAR + phase-resolved PROCAR → local d projections**: rotate complex d-orbital amplitudes for each site and write a compressed analysis file.

This means the project does not assume that one method is always best. For a simple crystal, rotating the structure can be cleanest. For tilted/distorted supercells, site-dependent Wannier projections are usually better. For already-finished VASP calculations, the PROCAR backend provides a post-processing route.

## Why this exists

VASP orbital labels such as `dxy`, `dz2`, and `dx2-y2` are defined relative to a global Cartesian frame. In tilted, rotated, distorted, supercell, or disordered coordination environments, those labels may not coincide with the physically meaningful local crystal-field orbitals. `localorb` derives local frames from the structure and applies them at the most appropriate stage of the workflow.

## Installation

```bash
python -m pip install -e .
```

For development/tests:

```bash
python -m pip install -e '.[test]'
pytest -q
```

Requirements are intentionally small: `numpy` and `pymatgen`.

## Quick start

Inspect all Ni-centered oxygen octahedra:

```bash
localorb inspect POSCAR --center Ni --ligand O --coordination 6
```

Write a reusable JSON report containing axes and rotation matrices:

```bash
localorb report POSCAR --center Ni --ligand O --coordination 6 \
  -o local_frames.json
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

Rotate the entire POSCAR so one site's local frame becomes the global frame:

```bash
localorb rotate-poscar POSCAR --site 10 --ligand O --coordination 6 \
  -o POSCAR.rotated
```

Print the corresponding real-d 5×5 basis transformation:

```bash
localorb d-matrix POSCAR --site 10 --ligand O --coordination 6
```

Rotate phase-resolved VASP d projections for all Ni sites:

```bash
localorb procar POSCAR --center Ni --ligand O --coordination 6 \
  --procar PROCAR -o PROCAR_LOCAL.npz
```

`PROCAR_LOCAL.npz` contains:

```text
coefficients[spin, k, band, site, orbital]   complex amplitudes
weights[spin, k, band, site, orbital]        |coefficient|^2
site_indices                                  original 0-based atom indices
orbital_names                                 dxy dyz dz2 dxz dx2-y2
spin_values                                   pymatgen spin labels
```

The PROCAR route requires phase-resolved complex projection amplitudes. `localorb` deliberately does **not** rotate already-squared orbital weights, because that would lose interference/phase information.

## Local-frame convention

For an approximately octahedral environment, `localorb` geometrically pairs the six ligands into three opposite pairs. By default the pair with the largest mean center-ligand distance is chosen as local `z` (`--z-policy longest`), which is useful for tetragonally elongated octahedra. The remaining pair closest to perpendicular to `z` seeds local `x`; then `y = z × x`, followed by orthonormalization and deterministic sign canonicalization.

The local-to-global frame matrix is

```text
R = [ x'  y'  z' ]
```

with local unit vectors as columns in global Cartesian coordinates:

```text
v_global = R @ v_local
v_local  = R.T @ v_global
```

## Wannier90 strategy

Wannier90 allows each projection to carry its own local `z` and `x` axes. `localorb wannier` writes explicit fractional-coordinate projection centers so symmetry-equivalent sites may still use different local frames in distorted or disordered supercells.

This is usually the preferred route for low-energy models because the Wannier basis is local from the beginning rather than corrected only after the calculation.

## Choosing a method

| Situation | Recommended localorb route |
|---|---|
| All relevant polyhedra share one orientation | `rotate-poscar` or `wannier` |
| Different sites have different octahedral tilts | `wannier` |
| Need a site-local tight-binding/Wannier basis | `wannier` |
| VASP calculation is already finished | `procar` |
| Need the explicit orbital basis transformation | `d-matrix` / Python API |
| Strongly ambiguous/distorted coordination | inspect first; use explicit/manual conventions before production |

## Python API

```python
from pymatgen.core import Structure
from localorb.frames import build_local_frame
from localorb.orbitals import d_rotation_matrix

s = Structure.from_file("POSCAR")
frame = build_local_frame(s, 10, ligand="O", coordination=6)
T_d = d_rotation_matrix(frame.rotation_local_to_global)
```

## Current scope and roadmap

Implemented in the initial release:

- octahedral (coordination-6) geometry inference;
- site-dependent Wannier90 projection generation;
- rigid POSCAR rotation;
- native real-d 5×5 rotation matrices;
- phase-resolved PROCAR d-amplitude rotation;
- automated tests on Python 3.10–3.12.

Next targets:

- rotation of Wannier/TB Hamiltonians between global and site-local orbital bases;
- explicit user axis overrides and ligand-pair overrides;
- square-planar, tetrahedral, trigonal-prismatic, and user-defined coordination templates;
- orbital-resolved band/DOS plotting helpers;
- validation utilities comparing equivalent representations before/after rotation.

## Scientific caution

Automatic local-axis construction is a geometric convention, not a substitute for physical judgment. In nearly regular octahedra there may be no unique physical choice of which equivalent axis should be called `z`. In strongly distorted or low-symmetry environments, always inspect the ligand pairing and local-frame report before using the generated basis in production calculations.

## License

MIT
