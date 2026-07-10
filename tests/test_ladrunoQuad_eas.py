"""LadrunoQuad ``-formulation eas`` — Q1/E4 Simo-Rifai enhanced assumed strain
(ADR 25, Phase 3, classTag 33007).

4 enhanced parameters (2 natural bubbles xi/eta x 2 dofs), the 2D sibling of
LadrunoBrick's E9 -- ported from the in-tree ``EnhancedQuad`` mode map
(Simo & Rifai 1990) but wired through LadrunoQuad's own architecture (formB,
shapeFunction, inner-Newton + static condensation shared with -formulation eas
on the brick). Small-strain only; no artificial stabilization (ADR 20's
beta-Tikhonov regularization was refuted for the brick and is NOT ported here
-- bare EAS).

Validation gates (per the ADR 25 P3 plan):
  * factory: 'eas' now builds (P2's refusal test flips to a build test)
  * rank / zero-energy: exactly 3 rigid-body modes, no spurious mechanism
    from the enhanced field
  * constant-stress patch: on an UNDISTORTED unit square (where the standard
    field alone already satisfies int M^T sigma = 0, i.e. alpha -> 0 at the
    solution) eas reproduces the same closed-form stress as std/bbar --
    the "reduce-to-std" gate
  * distorted-mesh relief: on the same distorted quad used throughout this
    suite, eas is strictly more flexible than std (the enhanced field can
    only relax the stiffness, never stiffen it -- same inequality pattern as
    the existing bbar-vs-std locking-relief test)
  * bending relief vs Euler-Bernoulli: a single, severely elongated
    (aspect-ratio 10) Q1 element in tip-shear bending is the classic
    Wilson/EAS benchmark -- std locks (parasitic shear), eas comes close to
    a converged (fine-mesh) reference. Comparative, not closed-form, to stay
    robust to the exact numerical ratio.

Plan: Ladruno_implementation/25_ladruno_plane_elements_adr.md.
"""
import math

import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy, check_constant_stress

pytestmark = [pytest.mark.zone_a]

# Same deliberately distorted (positive-Jacobian) quad used throughout the
# plane-element suite. CCW node order.
_QNODES = {
    1: (0.00, 0.00),
    2: (1.00, 0.10),
    3: (1.10, 1.00),
    4: (0.05, 0.95),
}
_QCONN = [1, 2, 3, 4]
_QBASE = [1, 2]            # fully fixed
_QLOADS = {                # mixed so the full 8x8 K is exercised
    3: (10.0, -4.0),
    4: (3.0, 7.0),
}
_THK = 0.7


def _static_solve(nodes, conn, base, loads, place_fn, E=1000.0, nu=0.3):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in nodes.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for n in base:
        ops.fix(n, 1, 1)

    place_fn(1)  # matTag = 1

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy) in loads.items():
        ops.load(n, fx, fy)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    disps = []
    for n in conn:
        disps.extend(ops.nodeDisp(n))
    forces = list(ops.eleResponse(1, "forces"))
    stresses = list(ops.eleResponse(1, "stresses"))
    return disps, stresses, forces


def _assert_close(a, b, label, rtol=1e-9, atol=1e-9):
    assert len(a) == len(b), f"{label}: length mismatch {len(a)} vs {len(b)}"
    for i, (x, y) in enumerate(zip(a, b)):
        tol = atol + rtol * max(abs(x), abs(y))
        assert abs(x - y) <= tol, (
            f"{label}[{i}]: {x!r} vs {y!r} (|d|={abs(x - y):.3e} > {tol:.3e})"
        )


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


# --------------------------------------------------------------------------
# factory: eas now builds
# --------------------------------------------------------------------------
def test_formulation_eas_builds():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QNODES.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoQuad", 1, *_QCONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-formulation", "eas")
    tags = ops.getEleTags() or []
    if isinstance(tags, int):
        tags = [tags]
    assert 1 in tags, "eas must build (ADR 25 P3)"


# --------------------------------------------------------------------------
# rank / zero-energy (no spurious mechanisms from the enhanced field)
# --------------------------------------------------------------------------
def test_quad_eas_no_spurious_modes():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QNODES.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.element("LadrunoQuad", 1, *_QCONN, 1, "-thick", _THK,
                "-type", "PlaneStrain", "-formulation", "eas")
    assert_zero_energy(list(_QNODES.keys()), ndf=2, ndm=2)


# --------------------------------------------------------------------------
# constant-stress patch / reduce-to-std: on an undistorted unit square the
# standard field alone already satisfies the enhanced orthogonality
# condition, so eas must reproduce the same closed-form stress as std.
# --------------------------------------------------------------------------
def test_quad_eas_constant_stress_patch():
    """Same patch test as test_ladrunoPlane_element.test_quad_constant_stress_patch
    (constant uniaxial stress on a unit square), run at -formulation eas."""
    E, nu, t, P = 1000.0, 0.25, 0.7, 5.0
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in sq.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-thick", t,
                "-type", "PlaneStrain", "-formulation", "eas")
    ops.fix(1, 1, 1)   # left-bottom: u=v=0
    ops.fix(4, 1, 0)   # left-top:    u=0
    ops.fix(2, 0, 1)   # right-bottom: v=0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 0.5 * P, 0.0)
    ops.load(3, 0.5 * P, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    sxx = P / (1.0 * t)
    check_constant_stress(1, 4, [sxx, 0.0, 0.0], E, sxx / E)


def test_quad_eas_matches_std_on_compatible_field():
    """On the same unit-square patch test, eas and std must agree to ~1e-8:
    when the standard displacement field already satisfies
    int M^T sigma = 0, the inner Newton drives alpha -> 0 and eas reduces
    to std exactly (up to Newton tolerance)."""
    E, nu, t, P = 1000.0, 0.25, 0.7, 5.0
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}

    def _run(formulation):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for tag, (x, y) in sq.items():
            ops.node(tag, x, y)
        ops.nDMaterial("ElasticIsotropic", 1, E, nu)
        ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-thick", t,
                    "-type", "PlaneStrain", "-formulation", formulation)
        ops.fix(1, 1, 1)
        ops.fix(4, 1, 0)
        ops.fix(2, 0, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 0.5 * P, 0.0)
        ops.load(3, 0.5 * P, 0.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 1.0)
        ops.algorithm("Linear")
        ops.analysis("Static")
        assert ops.analyze(1) == 0
        disps = []
        for n in (1, 2, 3, 4):
            disps.extend(ops.nodeDisp(n))
        return disps

    std = _run("std")
    eas = _run("eas")
    _assert_close(eas, std, "disp (eas vs std, compatible field)", rtol=1e-7, atol=1e-9)


# --------------------------------------------------------------------------
# distorted-mesh relief: eas is strictly more flexible than std (same
# inequality pattern as test_bbar_relieves_volumetric_locking).
# --------------------------------------------------------------------------
def test_quad_eas_at_least_as_flexible_as_std():
    std = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-type", "PlaneStrain",
                                              "-formulation", "std"))
    eas = _static_solve(_QNODES, _QCONN, _QBASE, _QLOADS,
                        lambda m: ops.element("LadrunoQuad", 1, *_QCONN, m,
                                              "-thick", _THK, "-type", "PlaneStrain",
                                              "-formulation", "eas"))
    assert _norm(eas[0]) >= _norm(std[0]) * (1.0 - 1e-9), (
        f"eas |u|={_norm(eas[0]):.6e} < std |u|={_norm(std[0]):.6e} "
        "(the enhanced field can only relax stiffness, never stiffen it)"
    )


# --------------------------------------------------------------------------
# bending relief: severely elongated single-element cantilever (aspect
# ratio 10) under tip shear -- the classic Wilson/EAS incompatible-mode
# benchmark. std locks (parasitic shear), eas should land much closer to a
# fine-mesh converged reference. Comparative (not closed-form) on purpose.
# --------------------------------------------------------------------------
def _cantilever_tip_deflection(formulation, nx, L=10.0, h=1.0, t=1.0,
                                E=1000.0, nu=0.0, Q=1.0):
    """nx elements along the length, 1 element through the height. Fixed at
    x=0 (both node layers, all dof); tip shear Q split evenly between the
    two tip nodes (y-direction). Returns the average tip uy."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    def bot(i): return 2 * i + 1
    def top(i): return 2 * i + 2

    for i in range(nx + 1):
        x = L * i / nx
        ops.node(bot(i), x, 0.0)
        ops.node(top(i), x, h)

    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    for i in range(nx):
        ops.element("LadrunoQuad", i + 1, bot(i), bot(i + 1), top(i + 1), top(i),
                    1, "-thick", t, "-type", "PlaneStress", "-formulation", formulation)

    ops.fix(bot(0), 1, 1)
    ops.fix(top(0), 1, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(bot(nx), 0.0, 0.5 * Q)
    ops.load(top(nx), 0.0, 0.5 * Q)

    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0

    return 0.5 * (ops.nodeDisp(bot(nx), 2) + ops.nodeDisp(top(nx), 2))


def test_quad_eas_relieves_bending_locking_vs_std():
    L, h = 10.0, 1.0   # aspect ratio 10 per element in the coarse (nx=1) mesh
    std_coarse = _cantilever_tip_deflection("std", nx=1, L=L, h=h)
    eas_coarse = _cantilever_tip_deflection("eas", nx=1, L=L, h=h)
    reference  = _cantilever_tip_deflection("std", nx=40, L=L, h=h)  # converged, aspect ratio 0.25/element

    assert abs(std_coarse) < abs(reference), (
        "sanity: the single coarse std element should lock (under-predict "
        f"the converged deflection) -- std={std_coarse:.6e}, ref={reference:.6e}"
    )
    err_std = abs(reference - std_coarse)
    err_eas = abs(reference - eas_coarse)
    assert err_eas < 0.5 * err_std, (
        f"eas should land much closer to the converged reference than std: "
        f"eas={eas_coarse:.6e} (err={err_eas:.3e}), "
        f"std={std_coarse:.6e} (err={err_std:.3e}), ref={reference:.6e}"
    )
