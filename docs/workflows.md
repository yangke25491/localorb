# Practical workflows

This document collects recommended `localorb` workflows for real calculations. The key rule is to decide the **physical local frame first**, then choose the backend that best matches the calculation.

## 1. Publication-quality single-site frame: manual atoms

Use this when you know which metal-ligand bonds define the physically meaningful axes.

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d
```

Validate the resulting rotation:

```bash
localorb validate POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d
```

Save the exact convention used:

```bash
localorb report POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d \
  -o Ni2_local_frame.json
```

The JSON report should be archived with the calculation so the orbital convention is reproducible.

## 2. Wannier90: define local orbitals before Wannierization

This is usually the preferred route for low-energy tight-binding work.

```bash
localorb wannier POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Copy the generated `begin projections ... end projections` block into the production `.win` file.

For a multi-site tilted structure, use the automatic provider when the same geometric rule is meaningful for all sites:

```bash
localorb validate POSCAR \
  --provider auto \
  --center Ni \
  --ligand O \
  --coordination 6 \
  -o frame_validation.json

localorb wannier POSCAR \
  --provider auto \
  --center Ni \
  --ligand O \
  --coordination 6 \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Inspect any site reported as `WARN` before using the automatic projection block in a production calculation.

### Interpreting the automatic diagnostics

`localorb validate` reports three geometry-sensitive metrics for automatic octahedral frames:

- `max_opposition_error = max(1 + cos(theta_pair))`
  - `0` is perfectly opposite.
- `raw_axis_orthogonality_error = max |axis_i dot axis_j|`
  - `0` is perfectly orthogonal before Gram-Schmidt.
- `z_gap_fraction`
  - measures how clearly one opposite ligand pair is distinguished from the nearest competing pair by mean bond length.

A small `z_gap_fraction` is not automatically wrong. In a nearly cubic octahedron all three axes may be physically equivalent, so geometry alone cannot uniquely decide which one should be called `z`.

For low-energy models where orbital identity matters strongly, use a manual or explicit-vector convention whenever the automatic z assignment is ambiguous.

## 3. In-plane octahedral rotation: fixed-z

If the physically meaningful c axis is already the POSCAR Cartesian z direction and only the in-plane ligand axes are rotated, use:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode fixed-z \
  -o POSCAR_fixed_z
```

The same frame can be used without rotating the structure:

```bash
localorb wannier POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode fixed-z \
  --orbitals dxy,dyz,dz2,dxz,dx2-y2
```

This is useful when the main orbital mixing is an in-plane `dx2-y2 <-> dxy` rotation and the local z axis should remain fixed.

## 4. Genuine 3D octahedral tilt: full-3d

For a real tilt where the local apical direction is not the global z direction:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d \
  -o POSCAR_full_3d
```

To compare the two conventions in one run:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --both-manual-modes \
  -o POSCAR_rotated
```

This writes `POSCAR_rotated_fixed_z` and `POSCAR_rotated_full_3d`.

## 5. VESTA periodic images

Inline image syntax is the clearest and most reproducible:

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4@0,0,1 \
  --plane-atom O7@0,0,1
```

Legacy-style shared or per-atom image options are also supported:

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --image 0,0,1 \
  --center-image 0,0,0
```

Priority is:

1. inline `ATOM@i,j,k`,
2. atom-specific image option,
3. shared `--image`,
4. nearest periodic image if no image is explicit.

## 6. Explicit crystallographic axes

When the local directions are known analytically, bypass ligand detection entirely:

```bash
localorb inspect POSCAR \
  --provider vectors \
  --center-atom Ni2 \
  --x-vector 1,1,0 \
  --z-vector 0,0,1
```

The same provider works with `wannier`, `rotate-poscar`, `d-matrix`, and `procar`.

## 7. Existing VASP calculation: rotate PROCAR projections

For a phase-resolved PROCAR:

```bash
localorb procar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --procar PROCAR \
  -o PROCAR_LOCAL.npz
```

The transformation is applied to the complex d-orbital amplitudes before squaring. Do not rotate already-squared `LORBIT=11` weights as if they were linear orbital coefficients.

The NPZ backend intentionally does not pretend to be a native VASP PROCAR file. It stores a clean numerical representation for downstream plotting or analysis.

## 8. Suggested workflow for large supercells

For a supercell with many transition-metal sites:

```bash
# 1. Build and validate automatic frames
localorb validate POSCAR \
  --provider auto --center Ni --ligand O \
  -o all_Ni_frame_validation.json

# 2. Save the actual frames
localorb report POSCAR \
  --provider auto --center Ni --ligand O \
  -o all_Ni_frames.json

# 3. Generate site-dependent Wannier90 projections
localorb wannier POSCAR \
  --provider auto --center Ni --ligand O \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

Then manually inspect sites flagged by `validate`. If only a few sites are problematic, define those frames explicitly rather than forcing one automatic rule onto every local environment.

## 9. Recommended comparison study

When introducing local frames into an existing workflow, compare at least three calculations:

1. original global-axis Wannier projections,
2. whole-structure rotated POSCAR when one global rotation is meaningful,
3. site-dependent local Wannier90 projections.

Compare:

- Wannier spreads,
- interpolation error,
- orbital-resolved band character,
- local Wannier-function shape,
- onsite energies,
- dominant hopping amplitudes,
- symmetry-equivalent hopping relations.

The purpose is not to force every representation to look identical. The purpose is to establish which representation gives the cleanest and most physically reproducible low-energy basis for the problem being studied.
