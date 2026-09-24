# localorb

`localorb` is a structure-aware toolkit for defining physically meaningful local orbital frames and carrying the same convention consistently through VASP, Wannier90, and tight-binding workflows.

The central design rule is:

> **separate how a local frame is defined from what that frame is used for.**

A manually selected `Ni2 -> O4` frame, an automatically detected NiO6 frame, and explicitly supplied crystallographic axes all produce the same `LocalFrame` object and can drive the same downstream backends.

Current version: **0.4.0**

## Documentation

- [`docs/workflows.md`](docs/workflows.md): practical VASP/Wannier90 workflows
- [`docs/frame_spec.md`](docs/frame_spec.md): multi-site frame JSON files for supercells
- [`docs/tb_rotation.md`](docs/tb_rotation.md): exact full-5d `wannier90_hr.dat` basis rotation
- [`CHANGELOG.md`](CHANGELOG.md): release history
- [`examples/frames.example.json`](examples/frames.example.json): mixed frame-provider example
- [`examples/basis_map.example.json`](examples/basis_map.example.json): complete-d Wannier basis map

---

## Physical model

Every provider produces

```text
R = [ x'  y'  z' ]
```

with the local Cartesian axes stored as columns in POSCAR Cartesian coordinates:

```text
v_global = R @ v_local
v_local  = R.T @ v_global
```

The same `R` then feeds multiple backends:

```text
POSCAR / structure
       |
       v
  frame provider
       |
       v
R = [x' y' z']
       |
       +--> inspect / validate / JSON report
       +--> Wannier90 local projections
       +--> rigid POSCAR rotation
       +--> real-d 5x5 rotation matrix
       +--> phase-resolved PROCAR rotation
       +--> exact full-d wannier90_hr.dat basis rotation
```

---

## Installation

```bash
git clone https://github.com/yangke25491/localorb.git
cd localorb
python -m pip install -e .
```

For development/tests:

```bash
python -m pip install -e '.[test]'
pytest -q
```

Requirements are intentionally small: `numpy` and `pymatgen`.

---

# 1. Frame providers

## Automatic octahedral geometry

Useful for large supercells where the same geometric convention should be applied repeatedly:

```bash
localorb inspect POSCAR \
  --provider auto \
  --center Ni \
  --ligand O \
  --coordination 6
```

The six nearest selected ligands are paired into three approximately opposite directions. The local `z` pair is selected with

```bash
--z-policy longest
```

or

```bash
--z-policy shortest
```

Automatic frames carry geometry-quality metadata and can be checked with `localorb validate`.

## Manual atom-defined frame

Preferred when the physically meaningful ligand directions are known:

```bash
localorb inspect POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d
```

The convention is

```text
center -> x-atom       defines local +x
center/x/plane atoms   define the local xy plane
```

`full-3d` constructs the complete local frame. `fixed-z` keeps the physical POSCAR Cartesian z axis fixed and only determines the in-plane axes:

```bash
--manual-mode fixed-z
```

This is useful when the relevant issue is an in-plane octahedral rotation rather than a genuine 3D tilt.

### VASP/VESTA-style atom selectors

Supported selectors include

```text
Ni2
O4
6
12
O4@0,0,1
```

Bare integers are 1-based by default. Periodic images may also be supplied with

```bash
--image 0,0,1
--center-image 0,0,0
--x-image 0,0,1
--y-image 0,0,1
```

If no image is explicit, the ligand is moved to its nearest periodic image relative to the center. Disable that behavior with

```bash
--no-nearest-image
```

The optional legacy VESTA rhombohedral/trigonal display compatibility layer is available through

```bash
--cart-frame vesta
```

but the safe default is always the physical POSCAR Cartesian frame:

```bash
--cart-frame poscar
```

## Explicit Cartesian vectors

When the desired axes are already known analytically:

```bash
localorb inspect POSCAR \
  --provider vectors \
  --center-atom Ni2 \
  --x-vector 1,1,0 \
  --z-vector 0,0,1
```

---

# 2. Validate and record the frame convention

For automatic supercell processing:

```bash
localorb validate POSCAR \
  --provider auto \
  --center Ni --ligand O \
  -o frame_validation.json
```

Automatic diagnostics include:

```text
max_opposition_error
raw_axis_orthogonality_error
z_gap_fraction
```

A warning is not automatically a failure. For example, a nearly cubic octahedron naturally has nearly degenerate choices for which equivalent axis should be called `z`.

Save the actual numerical frame matrices used in production:

```bash
localorb report POSCAR \
  --provider auto \
  --center Ni --ligand O \
  -o local_frames.json
```

---

# 3. Multi-site supercell frame files

For inequivalent sites, use one JSON file instead of one command per atom:

```bash
localorb inspect POSCAR --frames-file frames.json
localorb validate POSCAR --frames-file frames.json -o frame_validation.json
```

A single file may mix providers:

```json
{
  "frames": [
    {
      "provider": "manual",
      "center_atom": "Ni1",
      "x_atom": "O1",
      "plane_atom": "O2",
      "mode": "full-3d"
    },
    {
      "provider": "vectors",
      "center_atom": "Ni2",
      "x_vector": [1, 0, 0],
      "z_vector": [0, 0, 1]
    }
  ]
}
```

See [`docs/frame_spec.md`](docs/frame_spec.md) for the full schema.

---

# 4. Wannier90: define local orbitals at the source

For new Wannier/TB calculations this is generally the preferred route.

```bash
localorb wannier POSCAR \
  --provider auto \
  --center Ni --ligand O \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

For an inequivalent supercell:

```bash
localorb wannier POSCAR \
  --frames-file frames.json \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

The output uses Wannier90 site-specific local `z=` and `x=` projection axes, so different sites may carry different physically meaningful orbital frames without rotating the whole crystal.

---

# 5. Rotate an entire POSCAR

When one common local frame is meaningful for the entire structure:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 \
  --x-atom O4 \
  --plane-atom O7 \
  --manual-mode full-3d \
  -o POSCAR.rotated
```

Compare fixed-z and full-3d conventions in one run:

```bash
localorb rotate-poscar POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7 \
  --both-manual-modes \
  -o POSCAR_rotated
```

Rotated structures are rebuilt through pymatgen and written in Direct coordinates. This avoids the common Cartesian-input failure mode where only the lattice is rotated while Cartesian atom coordinates are left unchanged.

---

# 6. Real d-orbital rotation matrix

```bash
localorb d-matrix POSCAR \
  --provider manual \
  --center-atom Ni2 --x-atom O4 --plane-atom O7
```

The fixed orbital order is

```text
dxy dyz dz2 dxz dx2-y2
```

The implementation uses the symmetric-traceless rank-2 tensor representation, so it does not depend on an Euler-angle convention.

---

# 7. Existing VASP calculation: rotate PROCAR amplitudes

For a phase-resolved PROCAR:

```bash
localorb procar POSCAR \
  --frames-file frames.json \
  --procar PROCAR \
  -o PROCAR_LOCAL.npz
```

The transformation is applied to the complex d-orbital projection amplitudes **before squaring**. Already-squared orbital weights are not treated as linear coefficients.

The NPZ stores:

- local complex coefficients and weights,
- site indices,
- frame provider/mode,
- local/global rotation matrices,
- orbital names and spin channels,
- total-d-weight conservation diagnostics.

It deliberately does not pretend to be a native VASP PROCAR file.

---

# 8. Existing complete-d Wannier model: rotate `wannier90_hr.dat`

For an existing Wannier model containing a complete five-d block for every transformed site:

```bash
localorb tb-rotate POSCAR \
  --frames-file frames.json \
  --hr wannier90_hr.dat \
  --basis-map basis_map.json \
  -o wannier90_hr_local.dat \
  --transform-output basis_transform.npz
```

The exact basis convention is

```text
|new_a> = sum_m B[a,m] |old_m>
H_local(R) = B* H_global(R) B^T
```

where `B*` is element-wise complex conjugation. The current spatial d-orbital matrices are real, so for them this is equivalently `B H B^T`.

`localorb` verifies that `B` is unitary and that the Frobenius norm of every real-space Hamiltonian block is preserved.

## Why complete five-d blocks are required

A general spatial rotation of d orbitals acts in the full `l=2` space:

```text
dxy dyz dz2 dxz dx2-y2
```

A reduced `eg = {dz2, dx2-y2}` model is **not generally closed** under an arbitrary 3D rotation because the rotated orbitals can contain `t2g` components.

Therefore `localorb` deliberately rejects post-hoc arbitrary 3D rotation of a two-orbital `eg` or three-orbital `t2g` block. For reduced low-energy models, define the desired local axes during Wannier90 projection generation instead.

See [`docs/tb_rotation.md`](docs/tb_rotation.md) for the basis-map format and mathematical details.

---

## Strategy guide

| Situation | Recommended route |
|---|---|
| One common crystal orientation | `rotate-poscar` |
| Different local tilts at different sites | site-dependent `wannier` |
| Large supercell / disorder | `provider auto` + `validate` |
| Physical ligand directions known | `provider manual` |
| Crystallographic axes known analytically | `provider vectors` |
| Existing VASP result, need local orbital weights | `procar` |
| Existing complete five-d Wannier model | `tb-rotate` |
| Reduced `eg` Wannier model | define local axes during Wannierization |
| Only in-plane octahedral rotation matters | `manual-mode fixed-z` |
| Genuine 3D tilt | `manual-mode full-3d` or validated auto frame |

## Scientific caution

A local orbital frame is a physical convention, not just a plotting choice. Automatic geometry is useful for scale, but strong Jahn-Teller distortions, vacancies, unusual coordination, or nearly degenerate axes may make the physically appropriate convention ambiguous. Use `validate`, inspect flagged sites, and archive both the human-readable frame specification and the resulting numerical frame matrices with production calculations.

## License

MIT
