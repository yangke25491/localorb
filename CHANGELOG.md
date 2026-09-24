# Changelog

All notable user-facing changes to `localorb` are recorded here.

## 0.4.0 - 2026-09-24

### Added

- Multi-site JSON frame specifications through `--frames-file`.
  - Mix automatic, manual-atom, and explicit-vector providers in one supercell.
  - Optional shared defaults.
  - String or array-form periodic-image specifications.
  - Duplicate physical center sites are rejected.
- `localorb validate` geometry diagnostics for automatic octahedral frames.
- `docs/frame_spec.md` and `examples/frames.example.json` for reproducible supercell frame definitions.
- Standard `wannier90_hr.dat` parser/writer.
- Exact full-five-d Wannier Hamiltonian basis rotation.
- `localorb tb-rotate` command.
- Explicit `basis_map.json` describing which Wannier indices form each site's five real d orbitals.
- Exact basis-transform NPZ output containing the unitary transformation matrix and metadata.
- `examples/basis_map.example.json` and `docs/tb_rotation.md`.
- Numerical TB checks:
  - basis unitarity,
  - real-space Hamiltonian Frobenius-norm conservation,
  - regression tests for onsite-block eigenvalue conservation.
- Enhanced PROCAR NPZ output with frame matrices, providers/modes, and total-d-weight conservation diagnostics.

### Physical safeguards

- Arbitrary three-dimensional post-hoc orbital rotation is allowed only for a complete five-d Wannier block.
- Reduced `eg`-only or `t2g`-only blocks are rejected for a general spatial rotation because those subspaces are not generally closed under the full `l=2` rotation representation.
- For reduced low-energy models, the recommended workflow is to define the local axes during Wannier90 projection generation rather than force an incomplete post-hoc basis transformation.

### Notes

- A frame specification file and a Wannier basis map solve two different bookkeeping problems: the former maps physical sites to local Cartesian axes, while the latter maps Wannier indices to orbital blocks.
- Unmapped Wannier functions remain unchanged during `tb-rotate`.

## 0.3.0 - 2026-09-24

### Added

- Three interchangeable local-frame providers:
  - automatic octahedral geometry,
  - manual atom-defined frames,
  - explicit Cartesian vectors.
- VASP/VESTA-style atom selectors such as `Ni2`, `O4`, global integer indices, and explicit periodic images such as `O4@0,0,1`.
- Shared and per-atom image options: `--image`, `--center-image`, `--x-image`, `--plane-image`/`--y-image`.
- Manual `fixed-z` and `full-3d` frame conventions.
- Optional narrow VESTA rhombohedral/trigonal display-frame compatibility layer.
- `--both-manual-modes` for writing fixed-z and full-3d rotated POSCAR files in one run.
- Explicit-vector provider using `--x-vector` and `--z-vector`.
- Real-d `5 x 5` orbital rotation matrices in the order `dxy dyz dz2 dxz dx2-y2`.
- Phase-resolved PROCAR local-orbital rotation backend with NPZ output.
- `localorb validate` command for numerical frame validation and automatic-frame ambiguity diagnostics.
- Automatic octahedral quality metrics:
  - opposite-pair quality,
  - raw-axis orthogonality,
  - z-axis bond-length degeneracy.
- Regression tests for manual frames, periodic images, VESTA compatibility, explicit vectors, diagnostics, orbital rotations, and structure rotations.

### Changed

- Rotated POSCAR files are rebuilt through pymatgen and written in Direct coordinates. This fixes the legacy failure mode where a Cartesian-input POSCAR could have its lattice rotated without rotating the Cartesian atom coordinates.
- The safe default Cartesian frame for manual definitions is now the physical POSCAR Cartesian frame. VESTA compatibility must be requested explicitly or through the narrow `auto` compatibility mode.
- Local-frame convention is explicit throughout the codebase:

  `R = [x' y' z']` maps local Cartesian coordinates to POSCAR Cartesian coordinates.

### Notes

- The automatic octahedral provider is intended for scalable geometry-based processing, not as a substitute for a physically chosen orbital convention.
- A nearly degenerate z-axis warning can be completely legitimate for cubic or nearly cubic octahedra; it means the geometry alone does not uniquely distinguish which equivalent axis should be named z.

## 0.2.0

- Initial modular structure/Wannier90/POSCAR/PROCAR implementation.

## 0.1.0

- Initial project skeleton.
