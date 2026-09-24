# localorb

`localorb` builds physically meaningful local orbital frames directly from crystal structures and connects those frames to Wannier90, VASP-oriented workflows, and tight-binding/post-processing backends.

The project deliberately separates **how a local frame is defined** from **what the frame is used for**. A manually selected `Ni2 -> O4` frame, an automatically detected NiO6 frame, and an explicitly supplied pair of Cartesian axes all produce the same `LocalFrame` object and can feed the same downstream tools.

For complete end-to-end examples, see [`docs/workflows.md`](docs/workflows.md). Release notes are in [`CHANGELOG.md`](CHANGELOG.md).

## Core model

```text
POSCAR / structure
       |
       v
  frame provider
       |
       v
R = [x' y' z']   (local -> POSCAR Cartesian)
       |
       +--> validation / reproducibility report
       +--> Wannier90 local projections
       +--> rigid POSCAR rotation
       +--> real-d 5x5 rotation matrix
       +--> phase-resolved PROCAR rotation
       +--> future Wannier/TB Hamiltonian rotation
```

The local-to-global convention is

```text
v_global = R @ v_local
v_local  = R.T @ v_global
```

with `x'`, `y'`, `z'` stored as columns of `R`.

## Installation

```bash
python -m pip install -e .
```

Requirements are intentionally small: `numpy` and `pymatgen`.

---

# Frame providers

## 1. Automatic octahedral provider

Use this for large supercells or many equivalent transition-metal sites:

```bash
localorb inspect POSCAR \
  --provider auto \
  --center Ni \
  --ligand O \
  --coordination 6
```

For sixfold coordination, `localorb` pairs approximately opposite ligands into three axes. By default the pair with the largest mean center-ligand distance defines local `z`:

```bash
--z-policy longest
```

or use

```bash
--z-policy shortest
```

Automatic geometry is convenient, but it is still a convention. For publication-quality work where the physical axes are known, prefer the manual or explicit-vector providers.

## 2. Manual atom-defined provider

This is the preferred route when you know which bonds should define the local orbital axes.

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d
```

The three atoms mean:

```text
center atom
center -> x-atom                  defines local +x
center/x-atom/plane-atom          define the local xy plane
```

For `full-3d`:

```text
x = normalize(R_x - R_center)
y = Gram-Schmidt(R_plane - R_center, x)
z = x cross y
```

For an in-plane octahedral rotation where the physical POSCAR Cartesian `z` direction should remain fixed:

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
physical z = POSCAR Cartesian (0,0,1)
x = projection of center -> x-atom perpendicular to that z
y = z cross x
```

If a VESTA working frame is enabled, the same physical POSCAR `z` is first expressed in that working frame and converted back afterwards. This matches the behavior of the legacy rotation script.

This mode is particularly useful when the main issue is an in-plane `dx2-y2 <-> dxy` basis rotation rather than a genuine 3D tilt.

### VASP/VESTA-style selectors

The manual provider understands generated element-local labels and global integer indices:

```text
Ni2
O4
O7
6
12
15
```

Bare integer selectors are 1-based by default, matching common POSCAR/VESTA usage. Use

```bash
--index-base 0
```

for Python/pymatgen indexing.

Explicit periodic images can be written inline:

```text
O4@0,0,1
O7@-1,0,0
```

or with the legacy-script style options:

```bash
--image 0,0,1
--center-image 0,0,0
--x-image 0,0,1
--y-image 0,0,1
```

`--plane-image` and `--y-image` are aliases. Inline `ATOM@i,j,k` has highest priority, then the atom-specific image option, then shared `--image`.

If no explicit image is supplied, ligand atoms are moved to the nearest periodic image relative to the center. Disable this with:

```bash
--no-nearest-image
```

### Optional VESTA display-frame compatibility

The safe default is always:

```bash
--cart-frame poscar
```

For the rhombohedral/trigonal display convention used by the legacy `rotate_poscar_by_vesta_atom.py` workflow, `localorb` also provides:

```bash
--cart-frame auto
--cart-frame vesta
```

This is deliberately implemented as a narrow compatibility layer, not as a claim to reproduce every possible VESTA display orientation. `vesta` raises an error if the cell does not match the supported rhombohedral-like heuristic; `auto` falls back to the POSCAR Cartesian frame.

## 3. Explicit-vector provider

Use this when the physical axes are already known analytically or from another code:

```bash
localorb inspect POSCAR \
  --provider vectors \
  --center-atom Ni2 \
  --x-vector 1,1,0 \
  --z-vector 0,0,1
```

`z` is normalized first, `x` is projected perpendicular to `z`, and the right-handed frame is completed by `y = z x x`.

This provider is useful for known crystallographic directions, externally fitted local frames, or debugging the automatic/manual providers.

---

# Backends

## Inspect a frame

```bash
localorb inspect POSCAR --center Ni --ligand O
```

or

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7
```

The output includes local axes, provider, mode, ligand information, `det(R)`, and automatic geometry-quality metrics when available.

## Validate frames and flag ambiguous automatic choices

For a large supercell, validate before generating a production Wannier projection block:

```bash
localorb validate POSCAR \
  --provider auto \
  --center Ni --ligand O \
  -o frame_validation.json
```

Automatic octahedral validation reports:

```text
max_opposition_error
raw_axis_orthogonality_error
z_gap_fraction
```

The first two quantify how close the selected ligand geometry is to three opposite orthogonal axes. `z_gap_fraction` measures how clearly the chosen z pair is distinguished from the nearest competing pair by mean bond length.

A `WARN` does not automatically mean the frame is physically wrong. For example, a cubic octahedron naturally has nearly degenerate choices for which equivalent axis should be called z. Use `--strict` when you want warnings to produce a non-zero exit code for automated workflows.

## Export a reproducible JSON report

```bash
localorb report POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  -o local_frames.json
```

The report records local/global rotation matrices, selected atoms, periodic images, frame convention, and compatibility metadata.

## Wannier90: define local orbitals at the source

For Wannier/TB work this is generally the preferred route.

Automatic all-Ni example:

```bash
localorb wannier POSCAR \
  --provider auto \
  --center Ni --ligand O \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Manual one-site example:

```bash
localorb wannier POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  --manual-mode full-3d \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Explicit-vector example:

```bash
localorb wannier POSCAR \
  --provider vectors \
  --center-atom Ni2 \
  --x-vector 1,0,0 \
  --z-vector 0,0,1 \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

The generated Wannier90 projection block contains site-dependent local `z` and `x` axes, so each transition-metal site can have its own orbital frame.

## Rotate an entire POSCAR

Use whole-structure rotation only when one common frame is meaningful.

Automatic:

```bash
localorb rotate-poscar POSCAR \
  --provider auto \
  --site 10 --ligand O \
  -o POSCAR.rotated
```

Manual:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode fixed-z \
  -o POSCAR.rotated
```

Generate both manual conventions in one run:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  --both-manual-modes \
  -o POSCAR_rotated
```

which writes files with `_fixed_z` and `_full_3d` suffixes.

Rotated POSCAR output is always rebuilt from the rotated lattice plus unchanged fractional coordinates and written in Direct form. This avoids the classic Cartesian-input bug where the lattice is rotated but Cartesian atom coordinates are accidentally left unchanged.

## Real-d orbital rotation matrix

```bash
localorb d-matrix POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7
```

Orbital order is

```text
dxy dyz dz2 dxz dx2-y2
```

The implementation uses the symmetric-traceless rank-2 tensor representation and therefore does not depend on Euler-angle conventions.

## PROCAR local-orbital backend

For a phase-resolved VASP PROCAR, rotate the complex d-projection amplitudes rather than already-squared weights:

```bash
localorb procar POSCAR \
  --provider auto \
  --center Ni --ligand O \
  --procar PROCAR \
  -o PROCAR_LOCAL.npz
```

or use the manual/vector provider in exactly the same way.

The NPZ output stores local-frame complex coefficients and weights without pretending to be a native VASP PROCAR file.

---

# Choosing a strategy

| Situation | Recommended route |
|---|---|
| One common crystal orientation | `rotate-poscar` |
| Different octahedral tilts at different sites | site-dependent `wannier` |
| Large supercell / disorder | `provider auto` + `validate` |
| Physical ligand directions known | `provider manual` |
| Crystallographic axes known analytically | `provider vectors` |
| Existing VASP calculation, need local fatbands/PDOS | `procar` |
| Only in-plane rotation matters | `manual-mode fixed-z` |
| Genuine 3D tilt | `manual-mode full-3d` or auto geometry |
| Need legacy VESTA rhombohedral display convention | manual + `--cart-frame vesta` |

# Scientific caution

Automatic local-axis construction is geometric. Strong Jahn-Teller distortion, ligand vacancies, unusual coordination, or nearly degenerate axis choices can make the physical convention ambiguous. In those cases, use the manual or explicit-vector provider and save a JSON report alongside the calculation.

The VESTA compatibility transform is intentionally narrow and should not be treated as a universal VESTA coordinate converter.

# Roadmap

- Multiple manual frame specifications in one command for arbitrary supercells.
- Square-planar, tetrahedral, trigonal-prismatic, and user-defined coordination templates.
- Rotation of Wannier/TB Hamiltonians between global and site-local orbital bases.
- Validation plots comparing orbital-resolved bands before and after basis changes.
- Export helpers for downstream PyProcar/fatband workflows.

## License

MIT
