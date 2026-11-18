"""Truss scratch-matrix leak across `wipe` (WP-228, apeGmsh #2859 desk triage).

`Truss`, `TrussSection` (and `Truss2`) point `theMatrix` at a CLASS-STATIC
scratch matrix chosen by numDOF (`trussM6` for every 6-DOF truss).
`getTangentStiff()` / `getInitialStiff()` wrote only the `dimension x dimension`
translational blocks and zeroed the matrix only in the `L == 0` branch. So a
2-D `ndf = 3` truss (6 DOF, rotations at rows/cols 2 and 5) handed the
assembler whatever the previous user of `trussM6` left in its rotation slots.

The poisoner is a 3-D `ndf = 3` truss with `-rho`: its `getMass()` writes the
lumped mass `0.5*rho*L` at (2,2) and (5,5) -- for a 3-D truss those are the
z translations. After `wipe`, the first stiffness assembly of a 2-D frame
with ndf-3 truss braces then carries a fake rotational spring of that value at
every brace end, until some truss `getMass()` zeroes the matrix again. The
first eigen solve is silently wrong. The desk triage measured lambda_1 =
311.600977296998 clean and 311.601700975889 poisoned on a 3-storey ENT-braced
frame, the shift linear in rho, on stock openseespy 3.7.1 and the fork alike.

Oracle: the same frame solved twice in one process, once right after the
poison and once after a second `wipe` (by then the frame's own braces have
rewritten the scratch matrix cleanly). The eigenvalues must agree to round-off.
Before the fix the poisoned lambda_1 is high by ~1e-3 relative (rho chosen
large so the gap is unmistakable); after it the two runs are bit-for-bit
comparable. `corotTruss` / `corotTrussSection` zero `K` before filling and are
kept as controls so a regression in the sibling elements is caught too.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E = 2.0e11
RHO_POISON = 8000.0          # mass per unit length -> 0.5*rho*L = 4000 at (2,2),(5,5)
N_MODES = 4

ELEMENTS = ["Truss", "TrussSection", "corotTruss", "corotTrussSection"]


def _truss(kind, tag, i, j, area, mat, sec, *extra):
    if kind in ("Truss", "corotTruss"):
        ops.element(kind, tag, i, j, area, mat, *extra)
    else:
        ops.element(kind, tag, i, j, sec, *extra)


def _poison(kind):
    """3-D ndf-3 truss with -rho, mass assembled by an eigen solve."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 1.0, 0.0, 0.0)
    ops.node(3, 2.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(3, 0, 1, 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.section("Uniaxial", 1, 1, "P")
    _truss(kind, 1, 1, 2, 0.01, 1, 1, "-rho", RHO_POISON)
    _truss(kind, 2, 2, 3, 0.01, 1, 1, "-rho", RHO_POISON)
    ops.eigen("-fullGenLapack", 1)
    ops.wipe()


def _frame_eigen(kind):
    """3-storey, 1-bay 2-D frame (ndf 3) with X truss braces in every storey."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    bay, h = 6.0, 3.0
    for lvl in range(4):
        ops.node(10 * lvl + 1, 0.0, lvl * h)
        ops.node(10 * lvl + 2, bay, lvl * h)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 1, 1, 1)
    for lvl in range(1, 4):
        for n in (10 * lvl + 1, 10 * lvl + 2):
            ops.mass(n, 1.0e4, 1.0e4, 1.0)
    ops.geomTransf("Linear", 1)
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.section("Uniaxial", 1, 1, "P")
    tag = 100
    for lvl in range(3):
        a, b = 10 * lvl, 10 * (lvl + 1)
        for i, j in ((a + 1, b + 1), (a + 2, b + 2)):                 # columns
            tag += 1
            ops.element("elasticBeamColumn", tag, i, j, 0.02, E, 2.0e-4, 1)
        tag += 1                                                       # beam
        ops.element("elasticBeamColumn", tag, b + 1, b + 2, 0.02, E, 3.0e-4, 1)
        for i, j in ((a + 1, b + 2), (a + 2, b + 1)):                 # X braces
            tag += 1
            _truss(kind, tag, i, j, 0.002, 1, 1)
    lam = list(ops.eigen("-fullGenLapack", N_MODES))
    ops.wipe()
    return lam


@pytest.mark.parametrize("kind", ELEMENTS)
def test_truss_scratch_matrix_does_not_leak_across_wipe(kind):
    _poison(kind)
    lam_after_poison = _frame_eigen(kind)
    lam_clean = _frame_eigen(kind)

    assert all(math.isfinite(x) and x > 0.0 for x in lam_clean), lam_clean
    for k, (p, c) in enumerate(zip(lam_after_poison, lam_clean)):
        assert p == pytest.approx(c, rel=1e-12, abs=0.0), (
            "%s: mode %d lambda %.15g after a -rho truss poison vs %.15g clean "
            "(rel %.3e) -- the class-static truss matrix leaked rotation "
            "entries into the 2-D ndf-3 stiffness" % (kind, k + 1, p, c, (p - c) / c))
