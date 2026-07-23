"""End-to-end MPI smoke test for the independent two-level CMS command."""

from __future__ import annotations

import math

import openseesmp as ops


def fixed_free_chain(number_of_free_dofs: int) -> None:
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    for tag in range(1, number_of_free_dofs + 2):
        ops.node(tag, float(tag - 1))
    ops.fix(1, 1)
    ops.uniaxialMaterial("Elastic", 1, 1.0)
    for tag in range(1, number_of_free_dofs + 1):
        ops.element("truss", tag, tag, tag + 1, 1.0, 1)
    for tag in range(2, number_of_free_dofs + 2):
        ops.mass(tag, 1.0)


def exact_eigenvalues(number_of_free_dofs: int, number_of_modes: int) -> list[float]:
    return [
        2.0
        - 2.0
        * math.cos((2 * mode - 1) * math.pi / (2 * number_of_free_dofs + 1))
        for mode in range(1, number_of_modes + 1)
    ]


def solve(number_of_modes: int) -> list[float]:
    values = ops.eigen(
        "-ladrunoCMS",
        "-hierarchy",
        "logical",
        "-level1",
        2,
        "-level2",
        2,
        "-modesL2",
        4,
        "-modesL1",
        8,
        "-tol",
        1.0e-8,
        "-maxEnrich",
        2,
        "-denseMax",
        128,
        "-verifyAssembly",
        "full",
        "-verifyFullMaxBytes",
        1048576,
        number_of_modes,
    )
    return [float(value) for value in values]


def main() -> None:
    if ops.getNP() != 4:
        raise RuntimeError("this smoke test requires exactly four MPI ranks")
    number_of_free_dofs = 16
    number_of_modes = 4
    fixed_free_chain(number_of_free_dofs)
    expected = exact_eigenvalues(number_of_free_dofs, number_of_modes)

    first = solve(number_of_modes)
    second = solve(number_of_modes)  # exercises safe same-class SOE replacement
    for label, computed in (("first", first), ("second", second)):
        if len(computed) != number_of_modes:
            raise RuntimeError(f"{label} solve returned {len(computed)} modes")
        for mode, (actual, reference) in enumerate(zip(computed, expected), start=1):
            relative_error = abs(actual - reference) / max(abs(reference), 1.0e-30)
            if relative_error > 1.0e-8:
                raise RuntimeError(
                    f"{label} mode {mode}: relative eigenvalue error {relative_error:.3e}"
                )

    if ops.getPID() == 0:
        print("ladruno CMS OpenSeesPyMP two-level smoke passed")
    ops.wipe()


if __name__ == "__main__":
    main()
