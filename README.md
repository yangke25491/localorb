# localorb

`localorb` builds physically meaningful local orbital frames directly from crystal structures and connects those frames to Wannier90, VASP-oriented workflows, and tight-binding/post-processing backends.

The central idea is simple: a global Cartesian frame is often not the natural orbital frame of a transition-metal coordination polyhedron. There is no single best correction strategy for every material. `localorb` therefore separates **how a local frame is defined** from **what that frame is used for**.

## Core model

Every frame provider produces the same object:

```text
POSCAR / structure
       |
       v
  frame provider
       |
       v
R = [x' y' z']   (local -> global)
       |
       +--> Wannier90 local projections
       +--> rigid POSCAR rotation
       +--> real-d 5x5 rotation matrix
       +--> phase-resolved PROCAR rotation
       +--> future Wannier/TB Hamiltonian rotation
```

This means an automatically detected NiO6 frame and a manually chosen `Ni2 -> O4` frame can feed exactly the same downstream workflow.

## Frame providers

### 1. Automatic octahedral provider

Use this when many symmetry-related or supercell sites should be processed automatically.

```bash
localorb inspect POSCAR \
  --provider auto \
  --center Ni \
  --ligand O \
  --coordination 6
```

For coordination 6, `localorb` pairs approximately opposite ligands into three axes. By default the opposite pair with the largest mean center-ligand distance defines local `z`; the remaining pair most nearly orthogonal to `z` seeds local `x`.

```bash
--z-policy longest
--z-policy shortest
```

Automatic geometry is convenient for large supercells, but it is a convention rather than a substitute for physical judgment.

### 2. Manual atom-defined provider

Use this when you know which bonds should define the physically meaningful local axes.

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d
```

The three selected atoms mean:

```text
center atom
center -> x-atom       defines local +x
center/x-atom/plane-atom define the local xy plane
```

For `full-3d`, the frame is built as

```text
x = normalize(R_x - R_center)
y = Gram-Schmidt(R_plane - R_center, x)
z = x cross y
```

For systems where the physical `z` direction must remain the global Cartesian `z` axis, use:

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode fixed-z
```

Then

```text
z = (0,0,1)
x = in-plane projection of center -> x-atom
y = z cross x
```

This is useful for octahedral rotations around `c` where the main problem is an in-plane `dx2-y2 <-> dxy` basis rotation rather than a full 3D tilt.

## VASP/VESTA-style atom selectors

The manual provider accepts generated element-local labels and global integer indices:

```text
Ni2
O4
O7
6
12
15
```

Bare integers are **1-based by default**, matching POSCAR/VESTA-style atom numbering. Use `--index-base 0` if you explicitly want Python/pymatgen indexing.

Periodic images can be selected explicitly:

```text
O4@0,0,1
O7@-1,0,0
```

If no explicit image is supplied, manual ligand atoms are moved to their nearest periodic image relative to the selected center by default. Disable this with:

```bash
--no-nearest-image
```

An explicit `@i,j,k` always takes precedence over automatic nearest-image selection.

## Installation

```bash
python -m pip install -e .
```

Requirements are intentionally minimal: `numpy` and `pymatgen`.

## Inspect and export frames

Automatic:

```bash
localorb inspect POSCAR --center Ni --ligand O
```

Manual:

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  --manual-mode full-3d
```

Write a reproducible JSON report:

```bash
localorb report POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  -o local_frames.json
```

The report records the frame provider, mode, local axes, rotation matrices, periodic images, and selector metadata.

## Wannier90: define the local orbitals at the source

For most Wannier/TB workflows this is the preferred route because the initial projection functions themselves are local orbitals.

Automatic, all Ni sites:

```bash
localorb wannier POSCAR \
  --provider auto \
  --center Ni --ligand O \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Manual, one explicitly defined Ni site:

```bash
localorb wannier POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  --manual-mode full-3d \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Generate all five local d orbitals:

```bash
--orbitals dxy,dyz,dz2,dxz,dx2-y2
```

The generated Wannier90 projections include site-dependent local `z` and `x` axes, so a tilted supercell does not need one impossible global frame shared by all sites.

## Rotate an entire POSCAR

Use whole-structure rotation only when one common local frame is meaningful.

Automatic frame at one 0-based pymatgen site:

```bash
localorb rotate-poscar POSCAR \
  --provider auto \
  --site 10 --ligand O \
  -o POSCAR.rotated
```

Manual/VESTA-style frame:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode fixed-z \
  -o POSCAR.rotated
```

`localorb` rebuilds the rotated structure through pymatgen and writes a consistent Direct-coordinate POSCAR. This avoids the common error of rotating lattice vectors while leaving Cartesian atomic coordinates unrotated.

## Real-d orbital rotation matrix

Print the `5 x 5` real-d basis transformation for a frame:

```bash
localorb d-matrix POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7
```

Orbital order is

```text
dxy dyz dz2 dxz dx2-y2
```

The implementation uses the symmetric-traceless rank-2 tensor representation, avoiding dependence on Euler-angle conventions.

## PROCAR local-orbital backend

For a phase-resolved VASP PROCAR, rotate the complex d-projection amplitudes rather than the already-squared orbital weights:

```bash
localorb procar POSCAR \
  --provider auto \
  --center Ni --ligand O \
  --procar PROCAR \
  -o PROCAR_LOCAL.npz
```

or for one manually defined site:

```bash
localorb procar POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  --procar PROCAR \
  -o PROCAR_LOCAL.npz
```

The NPZ contains local-frame complex coefficients and weights without pretending to be a native VASP PROCAR file.

## Coordinate convention

The local-to-global frame matrix is

```text
R = [ x'  y'  z' ]
```

with local unit vectors stored as columns in global Cartesian coordinates:

```text
v_global = R @ v_local
v_local  = R.T @ v_global
```

Keeping this convention explicit is important when the same frame is reused for structure rotations, spherical-harmonic rotations, Wannier projectors, and Hamiltonian basis transformations.

## Choosing a strategy

| Situation | Recommended route |
|---|---|
| One common crystal orientation | `rotate-poscar` |
| Octahedra rotated differently at different sites | site-dependent `wannier` |
| Large supercell / disorder | `provider auto` + validation |
| Publication-quality frame definition known from chemistry | `provider manual` |
| Existing VASP calculation, need local fatbands/PDOS | `procar` |
| Only in-plane octahedral rotation matters | `manual-mode fixed-z` |
| Genuine 3D octahedral tilt | `manual-mode full-3d` or auto geometry |

## Scientific caution

Automatic local-axis construction is geometric. Strong Jahn-Teller distortion, ligand vacancies, unusual coordination, or near-degenerate axis choices can make the physical convention ambiguous. In those cases, use the manual provider and record the selected atoms/images in the JSON report.

## Roadmap

- Multiple manual frame specifications in one command for arbitrary supercells.
- Optional VESTA display-frame compatibility layer kept separate from the physical frame core.
- Square-planar, tetrahedral, trigonal-prismatic, and user-defined coordination templates.
- Rotation of Wannier/TB Hamiltonians between global and site-local orbital bases.
- Validation plots comparing orbital-resolved bands before and after basis changes.

## License

MIT
