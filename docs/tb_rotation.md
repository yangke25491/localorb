# Rotating Wannier90 Hamiltonians into local d-orbital frames

`localorb tb-rotate` changes the **orbital basis** of an existing `wannier90_hr.dat` while preserving the real-space Hamiltonian information.

The central convention is:

```text
|new_a> = sum_m B[a,m] |old_m>
```

so for a general complex basis matrix every real-space hopping block transforms as

```text
H_local(R) = B* H_global(R) B^T
```

where `B*` denotes element-wise complex conjugation and `B^T` the transpose.

For the present spatial real-d rotations, `B` is real, so this reduces to

```text
H_local(R) = B H_global(R) B^T = B H_global(R) B^dagger
```

The transformation is unitary when each transformed orbital block is complete.

## Important physical restriction: use the complete five-d subspace

A general three-dimensional spatial rotation acts in the full `l=2` space:

```text
dxy dyz dz2 dxz dx2-y2
```

and is represented by a `5 x 5` matrix. In general, a rotated `dz2` or `dx2-y2` orbital contains components from the nominal `t2g` channels as well.

Therefore `localorb` deliberately **rejects** a request to rotate only a reduced two-orbital `eg` block or three-orbital `t2g` block under an arbitrary spatial rotation.

This is a physical closure condition, not a software limitation.

If your Wannier model contains only `dz2` and `dx2-y2`, the safest choices are:

1. define the desired local axes during Wannier90 projection generation, or
2. rebuild a complete five-d Wannier basis before applying an arbitrary post-hoc 3D rotation.

Special symmetry-restricted rotations may preserve a smaller subspace, but `localorb` does not silently assume such a restriction.

---

## Inputs

A typical exact post-processing workflow needs four pieces of information:

```text
POSCAR
frames.json
wannier90_hr.dat
basis_map.json
```

`frames.json` says what local Cartesian frame belongs to each physical site.

`basis_map.json` says which Wannier indices form the five d orbitals of that site and in what orbital order they appear in `wannier90_hr.dat`.

### Example frame file

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

### Example basis map

```json
{
  "index_base": 1,
  "num_wann": 10,
  "groups": [
    {
      "center_atom": "Ni1",
      "indices": [1, 2, 3, 4, 5],
      "orbitals": ["dxy", "dyz", "dz2", "dxz", "dx2-y2"]
    },
    {
      "center_atom": "Ni2",
      "indices": [6, 7, 8, 9, 10],
      "orbitals": ["dxy", "dyz", "dz2", "dxz", "dx2-y2"]
    }
  ]
}
```

The orbital list may use another ordering, provided it contains each of the five real d orbitals exactly once. `localorb` permutes the `5 x 5` rotation matrix into that Wannier ordering.

Unmapped Wannier functions are left unchanged.

---

## Command

```bash
localorb tb-rotate POSCAR \
  --frames-file frames.json \
  --hr wannier90_hr.dat \
  --basis-map basis_map.json \
  -o wannier90_hr_local.dat \
  --transform-output basis_transform.npz
```

Outputs:

```text
wannier90_hr_local.dat
basis_transform.npz
```

The transformed `hr.dat` keeps:

- the same `R` vectors,
- the same degeneracy list,
- the same number of Wannier functions,
- all unmapped orbitals unchanged.

The NPZ file stores the exact full basis matrix `B` plus JSON metadata describing the applied site/orbital groups.

---

## Numerical checks

`localorb` checks that the row-wise basis coefficients are orthonormal:

```text
B B^dagger = I
```

and reports the maximum unitarity error.

Because the operation is unitary, the Frobenius norm of each real-space Hamiltonian block should also be invariant:

```text
||H_local(R)||_F = ||H_global(R)||_F
```

The command prints the maximum change across all `R` blocks.

For the `R=0` onsite block, the eigenvalues are also invariant under this basis transformation.

These checks verify the numerical transformation. They do not by themselves prove that the chosen local frame is the physically best convention; use `localorb validate` and inspect the frame definition separately.

---

## Relation to Wannier90 projections

There are two conceptually different workflows:

### Preferred when starting a new Wannier calculation

```text
POSCAR
  -> localorb frame
  -> site-dependent Wannier90 projection axes
  -> Wannierization
  -> wannier90_hr.dat already expressed in the desired local-orbital convention
```

This is usually the cleanest route.

### Useful for an existing complete-d Wannier model

```text
existing global-axis five-d wannier90_hr.dat
  + localorb frame
  -> exact unitary basis transformation
  -> local-axis wannier90_hr_local.dat
```

The second route is especially useful for controlled comparisons between global and local orbital conventions without repeating the DFT calculation.

---

## Spinful models

The current basis-map implementation treats the mapped indices literally. If a spinful Wannier Hamiltonian contains separate five-d blocks for each spin sector, list each five-orbital block as a separate group with non-overlapping indices. Reusing the same physical `center_atom` in multiple basis-map groups is allowed because the orbital blocks may belong to different spin sectors; the `frames.json` file itself still defines only one spatial frame per physical site.

Spinor rotations associated with changing the **spin quantization axis** are a separate operation and are not performed by the real-space orbital rotation described here.
