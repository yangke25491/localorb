from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Wannier90HR:
    comment: str
    num_wann: int
    degeneracies: np.ndarray
    r_vectors: np.ndarray
    hamiltonians: np.ndarray

    @property
    def nrpts(self) -> int:
        return int(len(self.degeneracies))


def read_hr(path: str | Path) -> Wannier90HR:
    """Read a standard ``wannier90_hr.dat`` file.

    Hamiltonian matrices are returned as ``hamiltonians[ir, m, n]`` with Python
    zero-based indices. The R-vector order and degeneracy list are preserved.
    """
    path = Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 4:
        raise ValueError(f"{path} is too short to be a Wannier90 hr.dat file")

    comment = lines[0]
    try:
        num_wann = int(lines[1].split()[0])
        nrpts = int(lines[2].split()[0])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid Wannier90 hr.dat header in {path}") from exc
    if num_wann <= 0 or nrpts <= 0:
        raise ValueError("num_wann and nrpts must be positive")

    line_index = 3
    degeneracies: list[int] = []
    while len(degeneracies) < nrpts:
        if line_index >= len(lines):
            raise ValueError("hr.dat ended while reading degeneracies")
        try:
            degeneracies.extend(int(token) for token in lines[line_index].split())
        except ValueError as exc:
            raise ValueError(
                f"Invalid degeneracy line {line_index + 1} in {path}"
            ) from exc
        line_index += 1
    if len(degeneracies) != nrpts:
        raise ValueError(
            f"Expected {nrpts} degeneracies but parsed {len(degeneracies)}; "
            "check the hr.dat header"
        )

    expected_records = nrpts * num_wann * num_wann
    remaining = [line for line in lines[line_index:] if line.strip()]
    if len(remaining) != expected_records:
        raise ValueError(
            f"Expected {expected_records} Hamiltonian records, found {len(remaining)}"
        )

    hamiltonians = np.zeros((nrpts, num_wann, num_wann), dtype=complex)
    r_vectors = np.zeros((nrpts, 3), dtype=int)

    record_index = 0
    for ir in range(nrpts):
        seen: set[tuple[int, int]] = set()
        current_r = None
        for _ in range(num_wann * num_wann):
            tokens = remaining[record_index].split()
            record_index += 1
            if len(tokens) < 7:
                raise ValueError("Each hr.dat Hamiltonian record must contain at least 7 fields")
            try:
                r = tuple(int(tokens[i]) for i in range(3))
                m = int(tokens[3]) - 1
                n = int(tokens[4]) - 1
                value = float(tokens[5]) + 1j * float(tokens[6])
            except ValueError as exc:
                raise ValueError(f"Invalid Hamiltonian record: {' '.join(tokens)}") from exc

            if not (0 <= m < num_wann and 0 <= n < num_wann):
                raise ValueError(f"Wannier indices out of range in record m={m+1}, n={n+1}")
            if current_r is None:
                current_r = r
                r_vectors[ir] = r
            elif r != current_r:
                raise ValueError(
                    "Unexpected R-vector change inside a num_wann^2 hr.dat block; "
                    "file ordering is not the standard Wannier90 ordering"
                )
            if (m, n) in seen:
                raise ValueError(f"Duplicate matrix element for R={r}, m={m+1}, n={n+1}")
            seen.add((m, n))
            hamiltonians[ir, m, n] = value

        if len(seen) != num_wann * num_wann:
            raise ValueError(f"Incomplete Hamiltonian block for R={current_r}")

    return Wannier90HR(
        comment=comment,
        num_wann=num_wann,
        degeneracies=np.asarray(degeneracies, dtype=int),
        r_vectors=r_vectors,
        hamiltonians=hamiltonians,
    )


def write_hr(path: str | Path, hr: Wannier90HR) -> None:
    """Write ``Wannier90HR`` data in standard ``wannier90_hr.dat`` ordering."""
    path = Path(path)
    degeneracies = np.asarray(hr.degeneracies, dtype=int)
    r_vectors = np.asarray(hr.r_vectors, dtype=int)
    hamiltonians = np.asarray(hr.hamiltonians, dtype=complex)

    if degeneracies.ndim != 1:
        raise ValueError("degeneracies must be one-dimensional")
    nrpts = len(degeneracies)
    if r_vectors.shape != (nrpts, 3):
        raise ValueError(f"r_vectors must have shape ({nrpts}, 3)")
    if hamiltonians.shape != (nrpts, hr.num_wann, hr.num_wann):
        raise ValueError(
            "hamiltonians must have shape "
            f"({nrpts}, {hr.num_wann}, {hr.num_wann})"
        )

    with path.open("w", encoding="utf-8") as handle:
        handle.write(f"{hr.comment}\n")
        handle.write(f"{hr.num_wann:12d}\n")
        handle.write(f"{nrpts:12d}\n")
        for start in range(0, nrpts, 15):
            chunk = degeneracies[start : start + 15]
            handle.write("".join(f"{int(value):5d}" for value in chunk) + "\n")

        for ir in range(nrpts):
            r1, r2, r3 = (int(value) for value in r_vectors[ir])
            for m in range(hr.num_wann):
                for n in range(hr.num_wann):
                    value = hamiltonians[ir, m, n]
                    handle.write(
                        f"{r1:5d}{r2:5d}{r3:5d}{m+1:5d}{n+1:5d}"
                        f"{value.real:18.10f}{value.imag:18.10f}\n"
                    )
