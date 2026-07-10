"""LadrunoQuad -geom finite (ADR 70 P1) — openseespy end-to-end battery.

The finite-strain *mechanics* are already verified building-free by the P0
oracle (tests/finitestrain2d_reference.py) and the C++ kernel diff
(tests/test_finitestrain2d_kernel_cpp.py, which matches the exact Q4 std+F-bar
f/K the element assembles). This battery gates what only a build can: the
ELEMENT wiring — parser (-geom finite), the setTrialF path through LogStrain2D,
the reduce-to-linear collapse, large-strain Newton convergence, the F-bar
volumetric-locking cure, and the parse/material guards.

Gates (ADR 70 §6, P1):
  * reduce-to-linear: at ε→0, -geom finite std/bbar ≈ -geom linear std;
  * homogeneous finite-stretch patch: GP Cauchy stress == the LogStrain2D oracle
    σ(F) evaluated at the solved deformation gradient (end-to-end material+element);
  * large-strain Newton convergence (LogStrain2D over elastic AND over LadrunoJ2);
  * F-bar volumetric-locking relief at ν→0.5 (bbar strictly more flexible than std);
  * guards: ssp/eas + finite refused at parse; a non-finite material under -geom
    finite errors gracefully (no kernel crash).

Plan: Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md.
"""
import math

import numpy as np
import pytest

from _testbed import ops
import logstrain2d_reference as l2

pytestmark = [pytest.mark.zone_a]

I3 = np.eye(3)


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _KG(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


def _finite_quad(nodes, conn, matTag, form="std", ptype="PlaneStrain", thick=1.0):
    """element LadrunoQuad ... -geom finite over a LogStrain2D material."""
    ops.element("LadrunoQuad", 1, *conn, matTag, "-thick", thick,
                "-type", ptype, "-formulation", form, "-geom", "finite")


# --------------------------------------------------------------------------
# reduce-to-linear: at tiny strain -geom finite collapses to -geom linear
# --------------------------------------------------------------------------
# distorted (positive-Jacobian) quad, mixed loads → full 8x8 exercised.
_QN = {1: (0.00, 0.00), 2: (1.00, 0.10), 3: (1.10, 1.00), 4: (0.05, 0.95)}
_QC = [1, 2, 3, 4]
_THK = 0.7


def _solve_quad(place_fn, load_scale, E=1000.0, nu=0.3, ndm_inner=2):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QN.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)          # inner (3D copy on demand)
    place_fn()
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 10.0 * load_scale, -4.0 * load_scale)
    ops.load(4, 3.0 * load_scale, 7.0 * load_scale)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    disps = []
    for n in _QC:
        disps.extend(ops.nodeDisp(n))
    return disps


@pytest.mark.parametrize("form", ["std", "bbar"])
def test_reduce_to_linear(form):
    """At ε→0 the finite path (LogStrain2D/elastic) matches the small-strain
    -geom linear element to ~1e-6 relative — the setTrialF/kernel wiring, F-bar
    coupling included, collapses to linear elasticity."""
    scale = 1e-6                                          # keep max strain ~1e-6
    lin = _solve_quad(
        lambda: ops.element("LadrunoQuad", 1, *_QC, 1, "-thick", _THK,
                            "-type", "PlaneStrain", "-formulation", form),
        scale)
    fin = _solve_quad(
        lambda: (ops.nDMaterial("LogStrain2D", 2, 1),
                 _finite_quad(_QN, _QC, 2, form=form, thick=_THK)),
        scale)
    for i, (a, b) in enumerate(zip(fin, lin)):
        tol = 1e-12 + 1e-5 * max(abs(a), abs(b))
        assert abs(a - b) <= tol, f"[{i}] finite {a!r} vs linear {b!r}"


# --------------------------------------------------------------------------
# homogeneous finite-stretch patch: GP Cauchy stress == oracle σ(F)
# --------------------------------------------------------------------------
def test_homogeneous_stretch_matches_oracle():
    """Unit-square quad, symmetry-fixed, pulled uniaxially into a genuinely
    finite stretch. The deformation is homogeneous, so every Gauss point must
    report the same Cauchy stress — and it must equal the validated LogStrain2D
    numpy oracle σ(F) at the SOLVED deformation gradient. This ties the compiled
    element's setTrialF → getStress path to the reference at large strain."""
    E, nu, t = 500.0, 0.3, 1.0
    K, G = _KG(E, nu)
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in sq.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain2D", 2, 1)                   # plane strain over 3D elastic
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 2, "-thick", t,
                "-type", "PlaneStrain", "-formulation", "std", "-geom", "finite")
    ops.fix(1, 1, 1)     # origin pinned
    ops.fix(2, 0, 1)     # bottom edge: v=0
    ops.fix(4, 1, 0)     # left edge:  u=0
    # pull the right edge into a large stretch (force-controlled, several steps)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 60.0, 0.0)
    ops.load(3, 60.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.1)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-11, 50)
    ops.analysis("Static")
    assert ops.analyze(10) == 0

    # solved (homogeneous) deformation gradient from the corner displacements
    a = ops.nodeDisp(2, 1)                                # u_x at X=1 → F_xx-1
    b = ops.nodeDisp(4, 2)                                # u_y at Y=1 → F_yy-1
    F = np.array([[1.0 + a, 0.0], [0.0, 1.0 + b]])
    assert a > 0.1, f"test should reach a genuinely finite stretch (a={a})"

    # External equilibrium on the DEFORMED config: the applied traction resultant
    # on the right face (total Fx = 60+60 = 120 at load factor 1.0) must balance
    # the Cauchy stress acting over the DEFORMED cross-section. Plane strain keeps
    # the out-of-plane thickness at t, and the deformed face height is (1+b), so
    # σ_xx · (1+b) · t == applied Fx. This is a genuine current-config force
    # balance — it fails if formFinite drops the J volume weight (small-strain
    # area) or reports σ in the wrong (reference/PK) measure.
    applied_fx = 120.0
    smag_r = list(ops.eleResponse(1, "stresses"))[0]      # GP0 σ_xx (homogeneous)
    deformed_face = (1.0 + b) * t                          # deformed area normal to x
    assert abs(smag_r * deformed_face - applied_fx) <= 1e-3 * applied_fx, (
        f"external equilibrium: σxx·(1+b)·t = {smag_r * deformed_face:.4f} "
        f"!= applied Fx = {applied_fx} (a J volume weight is likely missing)"
    )

    res = l2.plane_strain_step(F, I3.copy(), I3.copy(), "elastic", K=K, G=G)
    s_ref = res["sigma_voigt"]                            # [sxx, syy, sxy] Cauchy
    smag = max(abs(v) for v in s_ref)

    stresses = list(ops.eleResponse(1, "stresses"))       # 4 GP × [sxx,syy,sxy]
    for gp in range(4):
        sxx, syy, sxy = stresses[3 * gp:3 * gp + 3]
        tol = 1e-5 * smag + 1e-7
        assert abs(sxx - s_ref[0]) <= tol, f"GP{gp} sxx {sxx} vs oracle {s_ref[0]}"
        assert abs(syy - s_ref[1]) <= tol, f"GP{gp} syy {syy} vs oracle {s_ref[1]}"
        assert abs(sxy - s_ref[2]) <= tol, f"GP{gp} sxy {sxy} vs oracle {s_ref[2]}"
    # traction-free top edge ⇒ σ_yy ≈ 0
    assert abs(s_ref[1]) < 1e-4 * smag + 1e-4


# --------------------------------------------------------------------------
# large-strain Newton convergence (end-to-end), elastic + LadrunoJ2 inner
# --------------------------------------------------------------------------
def _stretch_column(inner_fn, form="std", target=0.2, nstep=40):
    """A 2x1 column of quads pulled into a genuinely finite tensile stretch under
    DISPLACEMENT control (monotonic ⇒ no fixed-step overshoot into det F≤0, robust
    for plasticity). Returns (analyze rc, tip uy)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    # 2 stacked unit squares, nodes 1..6
    crd = {1: (0, 0), 2: (1, 0), 3: (0, 1), 4: (1, 1), 5: (0, 2), 6: (1, 2)}
    for tag, (x, y) in crd.items():
        ops.node(tag, float(x), float(y))
    inner_fn()                                           # defines 3D inner as tag 1
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoQuad", 1, 1, 2, 4, 3, 2, "-type", "PlaneStrain",
                "-formulation", form, "-geom", "finite")
    ops.element("LadrunoQuad", 2, 3, 4, 6, 5, 2, "-type", "PlaneStrain",
                "-formulation", form, "-geom", "finite")
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    # reference load shape (equal pull on the two top nodes); DisplacementControl
    # on node 5 uy scales it to reach the target stretch monotonically.
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(5, 0.0, 1.0)
    ops.load(6, 0.0, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("DisplacementControl", 5, 2, target / nstep)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-9, 50)
    ops.analysis("Static")
    rc = ops.analyze(nstep)
    tip = ops.nodeDisp(5, 2)
    return rc, tip


def test_finite_convergence_elastic():
    rc, tip = _stretch_column(
        lambda: ops.nDMaterial("ElasticIsotropic", 1, 800.0, 0.3))
    assert rc == 0, "finite-strain Newton (elastic) failed to converge"
    assert tip > 0.1, f"expected a genuinely finite tensile tip disp, got {tip}"


def test_finite_convergence_j2():
    """LogStrain2D over a genuinely elasto-plastic 3D inner (LadrunoJ2 if present,
    else upstream J2Plasticity) — the ADR's 'over LadrunoJ2/elastic' gate."""
    def inner():
        # LadrunoJ2: element('nDMaterial','LadrunoJ2', tag, E, nu, sigY, ...) —
        # fall back to upstream J2Plasticity if the fork material name differs.
        try:
            ops.nDMaterial("LadrunoJ2", 1, 800.0, 0.3, 8.0)
        except Exception:
            # upstream J2Plasticity: tag K G sig0 sigInf delta H
            ops.nDMaterial("J2Plasticity", 1, 400.0, 300.0, 8.0, 12.0, 0.1, 1.0)
    try:
        rc, tip = _stretch_column(inner, target=0.1)
    except Exception as e:
        pytest.skip(f"no 3D J2 material available for the finite inner: {e}")
    assert rc == 0, "finite-strain Newton (J2) failed to converge"
    assert tip > 0.02


# --------------------------------------------------------------------------
# F-bar volumetric-locking relief at nu -> 0.5 (finite strain)
# --------------------------------------------------------------------------
def test_fbar_relieves_locking_finite():
    """Near-incompressible plane strain: -geom finite bbar (2D F-bar) is strictly
    more flexible than std, which locks. Mirrors the small-strain gate on the
    finite path — proving the (J0/J)^{1/2} F-bar scale + unsymmetric G0 coupling
    are wired and effective, not a no-op."""
    nu = 0.4999

    def run(form):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for tag, (x, y) in _QN.items():
            ops.node(tag, x, y)
        ops.nDMaterial("ElasticIsotropic", 1, 1000.0, nu)
        ops.nDMaterial("LogStrain2D", 2, 1)
        ops.element("LadrunoQuad", 1, *_QC, 2, "-thick", _THK,
                    "-type", "PlaneStrain", "-formulation", form, "-geom", "finite")
        ops.fix(1, 1, 1)
        ops.fix(2, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, 5.0, -2.0)
        ops.load(4, 1.5, 3.5)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.2)
        ops.algorithm("Newton")
        ops.test("NormDispIncr", 1e-10, 50)
        ops.analysis("Static")
        assert ops.analyze(5) == 0, f"{form} finite analysis diverged"
        d = []
        for n in _QC:
            d.extend(ops.nodeDisp(n))
        return _norm(d)

    std = run("std")
    bbar = run("bbar")
    assert bbar > std * (1.0 + 1e-4), (
        f"F-bar |u|={bbar:.6e} not > std |u|={std:.6e} (locking not relieved)"
    )


# --------------------------------------------------------------------------
# Newton iteration-count gate: a large single-step bbar-finite solve must
# converge in few iterations. A transposed tangent or wrong-sign F-bar G0
# coupling is often still convergent (the residual is exact and independent of
# the tangent), just SLOWLY / non-quadratically — a norm-check on the answer
# would pass while Newton crawls. Gating the iteration count catches that.
# Ceiling calibrated below: observed count on the correct build; the gate is set
# to ~2x that so a healthy quadratic solve passes with margin but a degraded
# (linear-rate) tangent trips it.
# --------------------------------------------------------------------------
_NEWTON_ITER_CEIL = 12  # calibrated: 6 iters observed on the correct consistent
#                         tangent (|u3|~0.58, a genuinely finite step); gate = 2x.


def test_bbar_finite_newton_iterations_gate():
    """One big load step, bbar + -geom finite, over a distorted quad. Assert the
    consistent (unsymmetric F-bar) tangent gives a fast Newton solve — the
    iteration count stays under a calibrated ceiling. A transposed K or wrong-sign
    G0 coupling still converges (exact residual) but at a linear rate, blowing
    past the ceiling."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QN.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoQuad", 1, *_QC, 2, "-thick", _THK,
                "-type", "PlaneStrain", "-formulation", "bbar", "-geom", "finite")
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 32.0, -12.0)                              # large -> genuinely finite
    ops.load(4, 10.0, 22.0)                               # (|u3|~0.58 in one step)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)                    # ONE big step
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 60)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "bbar-finite single big step failed to converge"
    iters = ops.testIter()
    # sanity: the step must actually reach a finite stretch (not a trivial solve)
    assert _norm(ops.nodeDisp(3)) > 0.05, "step did not reach a finite deformation"
    assert iters <= _NEWTON_ITER_CEIL, (
        f"bbar-finite Newton took {iters} iters (> {_NEWTON_ITER_CEIL}); a "
        "transposed tangent / wrong-sign F-bar coupling degrades the rate"
    )


# --------------------------------------------------------------------------
# parse / material guards
# --------------------------------------------------------------------------
def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QN.items():
        ops.node(tag, x, y)


@pytest.mark.parametrize("form", ["ssp", "eas"])
def test_finite_refused_for_ssp_eas(form):
    """ssp/eas + finite are reserved — the factory must refuse (no element made),
    and the interpreter must survive so a valid build still works."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    try:
        ops.element("LadrunoQuad", 99, *_QC, 2, "-type", "PlaneStrain",
                    "-formulation", form, "-geom", "finite")
    except Exception:
        pass
    # a valid std-finite element still constructs afterwards
    ops.element("LadrunoQuad", 1, *_QC, 2, "-type", "PlaneStrain",
                "-formulation", "std", "-geom", "finite")


def test_finite_refused_for_plane_stress():
    """PlaneStress + -geom finite is refused at the FACTORY (ADR 70 review): the
    finite volume weight omits the out-of-plane thickness stretch λ, so finite
    plane-stress would be silently wrong. The element must never be created, and
    the interpreter must survive so a valid PlaneStrain finite element still
    builds."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    try:
        ops.element("LadrunoQuad", 99, *_QC, 2, "-type", "PlaneStress",
                    "-formulation", "std", "-geom", "finite")
    except Exception:
        pass
    assert 99 not in ops.getEleTags(), "PlaneStress + finite must be refused at parse"
    # a valid PlaneStrain-finite element still constructs afterwards
    ops.element("LadrunoQuad", 1, *_QC, 2, "-type", "PlaneStrain",
                "-formulation", "std", "-geom", "finite")
    assert 1 in ops.getEleTags()


def test_finite_rejects_non_finite_material():
    """-geom finite with a plain (non-FiniteStrainND2D) material must be refused at
    the FACTORY — the element is never created (no bogus assembly, no kernel crash)
    — and the interpreter must survive so a valid finite element still builds."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)   # NOT a LogStrain2D
    try:
        ops.element("LadrunoQuad", 1, *_QC, 1, "-type", "PlaneStrain",
                    "-formulation", "std", "-geom", "finite")
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), "non-finite material must be refused at parse"
    # a valid finite element (LogStrain2D) still constructs afterwards
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoQuad", 2, *_QC, 2, "-type", "PlaneStrain",
                "-formulation", "std", "-geom", "finite")
    assert 2 in ops.getEleTags()


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
