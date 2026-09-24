# Multi-site frame specification files

For supercells with several inequivalent local orbital frames, command-line arguments become cumbersome. `localorb` therefore accepts a JSON frame specification file through:

```bash
--frames-file frames.json
```

The same file can be used with:

```text
inspect
report
validate
wannier
procar
```

For single-frame commands (`rotate-poscar`, `d-matrix`), choose one entry with:

```bash
--frame-index N
```

where `N` is the zero-based position in the JSON `frames` array. If the file contains exactly one frame, `--frame-index` is optional.

See [`../examples/frames.example.json`](../examples/frames.example.json) for a complete mixed-provider example.

## File structure

```json
{
  "defaults": {
    "index_base": 1,
    "nearest_image": true,
    "cart_frame": "poscar"
  },
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
      "x_vector": [1, 1, 0],
      "z_vector": [0, 0, 1]
    }
  ]
}
```

`defaults` is optional. Values defined inside an individual frame override the top-level defaults.

Each center site may appear only once. This prevents silently generating two incompatible orbital conventions for the same atom.

---

## Manual frame entries

Required fields:

```json
{
  "provider": "manual",
  "center_atom": "Ni1",
  "x_atom": "O1",
  "plane_atom": "O2"
}
```

Optional fields:

```text
mode                full-3d | fixed-z
nearest_image       true | false
index_base           0 | 1
cart_frame           poscar | auto | vesta
image                shared periodic image
center_image         center-specific periodic image
x_image              x-atom-specific periodic image
plane_image          plane-atom-specific periodic image
y_image              alias of plane_image
```

Inline images remain supported and have highest priority:

```json
{
  "provider": "manual",
  "center_atom": "Ni2",
  "x_atom": "O4@0,0,1",
  "plane_atom": "O7@0,0,1",
  "mode": "fixed-z"
}
```

---

## Explicit-vector entries

```json
{
  "provider": "vectors",
  "center_atom": "Ni3",
  "x_vector": [1.0, 1.0, 0.0],
  "z_vector": [0.0, 0.0, 1.0]
}
```

Vectors may also be written as strings:

```json
{
  "provider": "vectors",
  "center_atom": "Ni3",
  "x_vector": "1,1,0",
  "z_vector": "0,0,1"
}
```

`z_vector` is normalized first. `x_vector` is projected perpendicular to z before the right-handed orthonormal frame is completed.

---

## Automatic entries

An automatic entry can use a zero-based site index:

```json
{
  "provider": "auto",
  "site": 10,
  "ligand": "O",
  "coordination": 6,
  "z_policy": "longest"
}
```

or a VESTA-style atom selector:

```json
{
  "provider": "auto",
  "center_atom": "Ni4",
  "ligand": "O",
  "coordination": 6,
  "z_policy": "longest"
}
```

Optional fields are:

```text
ligand
coordination
cutoff
z_policy             longest | shortest
index_base            for center_atom integer selectors
```

---

# Typical supercell workflow

## 1. Define the frames

Create `frames.json` manually, mixing providers where useful. For example, most sites can use automatic geometry while a few strongly distorted sites are defined manually.

## 2. Inspect all frames

```bash
localorb inspect POSCAR --frames-file frames.json
```

## 3. Validate

```bash
localorb validate POSCAR \
  --frames-file frames.json \
  -o frame_validation.json
```

## 4. Save the exact frame matrices

```bash
localorb report POSCAR \
  --frames-file frames.json \
  -o local_frames.json
```

## 5. Generate one Wannier90 projection block

```bash
localorb wannier POSCAR \
  --frames-file frames.json \
  --orbitals dz2,dx2-y2 \
  -o projections.win
```

## 6. Rotate an existing phase-resolved PROCAR in the same local convention

```bash
localorb procar POSCAR \
  --frames-file frames.json \
  --procar PROCAR \
  -o PROCAR_LOCAL.npz
```

This keeps the local-frame convention used for Wannier90 and VASP post-processing synchronized.

---

# Reproducibility recommendation

For a production calculation, keep these files together:

```text
POSCAR
frames.json
local_frames.json
frame_validation.json
wannier90.win
```

`frames.json` records the human/algorithmic choices. `local_frames.json` records the actual numerical rotation matrices generated from those choices. Keeping both makes the orbital convention reproducible even if local-frame construction algorithms evolve later.
