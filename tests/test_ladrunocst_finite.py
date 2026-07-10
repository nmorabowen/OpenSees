"""LadrunoCST -geom finite (ADR 70 P2) — openseespy end-to-end battery.

The finite-strain *mechanics* are already verified building-free by the P0
oracle (tests/finitestrain2d_reference.py, which includes the T3 std lane AND
pins the T3 F-bar no-op) and the C++ kernel diff
(tests/test_finitestrain2d_kernel_cpp.py). This battery gates what only a build
can: the ELEMENT wiring — parser (-geom finite), the setTrialF path through
LogStrain2D at the single centroid GP, the reduce-to-linear collapse,
large-strain Newton convergence, the det-F step-cut, and the honest-baseline
volumetric over-stiffness (ADR 70 §3.2: element-level F-bar is a NO-OP on the
constant-strain T3, so CST-finite locks — LadrunoLST (P3) is the usable
finite-strain triangle).

Gates (ADR 70 §6, P2):
  * homogeneous finite-stretch patch: GP Cauchy stress == the LogStrain2D oracle
    σ(F) at the solved F, on BOTH triangles, + deformed-config equilibrium;
  * reduce-to-linear: at ε→0, -geom finite ≈ -geom linear;
  * det F ≤ 0 step-cut: an inverting step fails gracefully (analyze rc != 0,
    interpreter survives);
  * honest baseline pinned: near-incompressible CST-finite locks hard relative
    to the F-bar quad (the ratio collapses as ν→0.5);
  * guards: PlaneStress+finite and non-finite materials refused at parse.

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


# distorted (positive-Jacobian) mesh shared with the quad battery, split into
# two CCW triangles along the 1-3 diagonal; mixed loads → all free DOFs exercised.
_QN = {1: (0.00, 0.00), 2: (1.00, 0.10), 3: (1.10, 1.00), 4: (0.05, 0.95)}
_TRIS = [(1, 2, 3), (1, 3, 4)]
_THK = 0.7


def _place_tris(matTag, tris=_TRIS, thick=_THK, geom=None, ptype="PlaneStrain"):
    for k, conn in enumerate(tris):
        extra = ("-geom", geom) if geom else ()
        ops.element("LadrunoCST", k + 1, *conn, matTag, "-thick", thick,
                    "-type", ptype, *extra)


# --------------------------------------------------------------------------
# reduce-to-linear: at tiny strain -geom finite collapses to -geom linear
# --------------------------------------------------------------------------
def _solve_tris(place_fn, load_scale):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QN.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)    # inner (3D copy on demand)
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
    for n in sorted(_QN):
        disps.extend(ops.nodeDisp(n))
    return disps


def test_reduce_to_linear():
    """At ε→0 the finite path (LogStrain2D/elastic) matches the small-strain
    -geom linear element to ~1e-5 relative — the setTrialF/kernel wiring
    collapses to linear elasticity."""
    scale = 1e-6                                          # keep max strain ~1e-6
    lin = _solve_tris(lambda: _place_tris(1), scale)
    fin = _solve_tris(
        lambda: (ops.nDMaterial("LogStrain2D", 2, 1),
                 _place_tris(2, geom="finite")),
        scale)
    for i, (a, b) in enumerate(zip(fin, lin)):
        tol = 1e-12 + 1e-5 * max(abs(a), abs(b))
        assert abs(a - b) <= tol, f"[{i}] finite {a!r} vs linear {b!r}"


# --------------------------------------------------------------------------
# homogeneous finite-stretch patch: GP Cauchy stress == oracle σ(F)
# --------------------------------------------------------------------------
def test_homogeneous_stretch_matches_oracle():
    """Unit square tiled by two CSTs, symmetry-fixed, pulled uniaxially into a
    genuinely finite stretch. The homogeneous deformation is exactly
    representable by the linear T3 field, so BOTH triangles must report the same
    Cauchy stress — and it must equal the validated LogStrain2D numpy oracle
    σ(F) at the SOLVED deformation gradient. This ties the compiled element's
    setTrialF → getStress path to the reference at large strain."""
    E, nu, t = 500.0, 0.3, 1.0
    K, G = _KG(E, nu)
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    tris = [(1, 2, 3), (1, 3, 4)]
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in sq.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain2D", 2, 1)                   # plane strain over 3D elastic
    for k, conn in enumerate(tris):
        ops.element("LadrunoCST", k + 1, *conn, 2, "-thick", t,
                    "-type", "PlaneStrain", "-geom", "finite")
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

    # External equilibrium on the DEFORMED config: applied Fx = 120 must balance
    # the Cauchy stress over the DEFORMED cross-section σ_xx·(1+b)·t. Fails if
    # formFinite drops the J volume weight or reports σ in a reference measure.
    applied_fx = 120.0
    sxx_e1 = list(ops.eleResponse(1, "stresses"))[0]      # centroid GP σ_xx
    deformed_face = (1.0 + b) * t
    assert abs(sxx_e1 * deformed_face - applied_fx) <= 1e-3 * applied_fx, (
        f"external equilibrium: σxx·(1+b)·t = {sxx_e1 * deformed_face:.4f} "
        f"!= applied Fx = {applied_fx} (a J volume weight is likely missing)"
    )

    res = l2.plane_strain_step(F, I3.copy(), I3.copy(), "elastic", K=K, G=G)
    s_ref = res["sigma_voigt"]                            # [sxx, syy, sxy] Cauchy
    smag = max(abs(v) for v in s_ref)

    for ele in (1, 2):                                    # both triangles, 1 GP each
        sxx, syy, sxy = list(ops.eleResponse(ele, "stresses"))
        tol = 1e-5 * smag + 1e-7
        assert abs(sxx - s_ref[0]) <= tol, f"ele{ele} sxx {sxx} vs oracle {s_ref[0]}"
        assert abs(syy - s_ref[1]) <= tol, f"ele{ele} syy {syy} vs oracle {s_ref[1]}"
        assert abs(sxy - s_ref[2]) <= tol, f"ele{ele} sxy {sxy} vs oracle {s_ref[2]}"
    # traction-free top edge ⇒ σ_yy ≈ 0
    assert abs(s_ref[1]) < 1e-4 * smag + 1e-4


# --------------------------------------------------------------------------
# large-strain Newton convergence (end-to-end), elastic + J2 inner
# --------------------------------------------------------------------------
def _stretch_column(inner_fn, target=0.2, nstep=40):
    """A 2x1 column of unit squares, each split into two CSTs, pulled into a
    genuinely finite tensile stretch under DISPLACEMENT control (monotonic ⇒ no
    fixed-step overshoot into det F≤0, robust for plasticity)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    crd = {1: (0, 0), 2: (1, 0), 3: (0, 1), 4: (1, 1), 5: (0, 2), 6: (1, 2)}
    for tag, (x, y) in crd.items():
        ops.node(tag, float(x), float(y))
    inner_fn()                                           # defines 3D inner as tag 1
    ops.nDMaterial("LogStrain2D", 2, 1)
    tris = [(1, 2, 4), (1, 4, 3), (3, 4, 6), (3, 6, 5)]  # all CCW
    for k, conn in enumerate(tris):
        ops.element("LadrunoCST", k + 1, *conn, 2, "-type", "PlaneStrain",
                    "-geom", "finite")
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
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
    else upstream J2Plasticity) — the end-to-end plasticity lane."""
    def inner():
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
# det F ≤ 0 step-cut: an inverting trial step fails GRACEFULLY
# --------------------------------------------------------------------------
def test_detF_step_cut():
    """One huge displacement-controlled compressive step drives a triangle
    inverted (det F ≤ 0). setTrialF must reject the trial (LogStrain2D guard →
    updateFinite returns < 0), so analyze reports failure — no crash, no
    silently-accepted inverted state — and the interpreter survives to run a
    valid model afterwards."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
    for tag, (x, y) in sq.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    for k, conn in enumerate([(1, 2, 3), (1, 3, 4)]):
        ops.element("LadrunoCST", k + 1, *conn, 2, "-type", "PlaneStrain",
                    "-geom", "finite")
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, -1.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    # push the loaded corner 1.5 to the LEFT in ONE step: the first Newton trial
    # is the linear prediction u_x(3) = -1.5, which turns the unit square inside
    # out (det F < 0 in both triangles).
    ops.integrator("DisplacementControl", 3, 1, -1.5)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.analysis("Static")
    rc = ops.analyze(1)
    assert rc != 0, "an element-inverting step must FAIL (det F step-cut), not pass"

    # the interpreter survives: a fresh valid finite model still solves
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in sq.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    for k, conn in enumerate([(1, 2, 3), (1, 3, 4)]):
        ops.element("LadrunoCST", k + 1, *conn, 2, "-type", "PlaneStrain",
                    "-geom", "finite")
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.fix(4, 1, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# --------------------------------------------------------------------------
# the honest baseline, pinned: near-incompressible CST-finite LOCKS
# --------------------------------------------------------------------------
def test_cst_finite_locks_near_incompressible():
    """ADR 70 §3.2 / P2 gate: element-level F-bar is a no-op on the T3 (J is
    element-constant — nothing to average), so CST -geom finite is the honest,
    volumetrically over-stiff baseline. Pin it on a bending-dominated
    plane-strain cantilever (8x2 cells; CST = each cell split in two triangles)
    against the locking-free F-bar quad (bbar + finite) on the same grid: the
    CST/quad tip-deflection ratio COLLAPSES as ν→0.5 (calibrated on the correct
    build: 0.50 at ν=0.3 → 0.069 at ν=0.4999, a 7.3x collapse; gate at 0.35x
    with margin). If someone later bolts a fake per-element 'F-bar' onto the T3
    (or the quad's F-bar regresses), this gate trips."""
    NX, NY, L, H = 8, 2, 4.0, 1.0

    def tip(nu, cst):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)

        def tag(i, j):
            return j * (NX + 1) + i + 1

        for j in range(NY + 1):
            for i in range(NX + 1):
                ops.node(tag(i, j), L * i / NX, H * j / NY)
        ops.nDMaterial("ElasticIsotropic", 1, 1000.0, nu)
        ops.nDMaterial("LogStrain2D", 2, 1)
        e = 1
        for j in range(NY):
            for i in range(NX):
                n1, n2, n3, n4 = tag(i, j), tag(i + 1, j), tag(i + 1, j + 1), tag(i, j + 1)
                if cst:
                    ops.element("LadrunoCST", e, n1, n2, n3, 2,
                                "-type", "PlaneStrain", "-geom", "finite")
                    e += 1
                    ops.element("LadrunoCST", e, n1, n3, n4, 2,
                                "-type", "PlaneStrain", "-geom", "finite")
                    e += 1
                else:
                    ops.element("LadrunoQuad", e, n1, n2, n3, n4, 2,
                                "-type", "PlaneStrain", "-formulation", "bbar",
                                "-geom", "finite")
                    e += 1
        for j in range(NY + 1):
            ops.fix(tag(0, j), 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for j in range(NY + 1):
            ops.load(tag(NX, j), 0.0, -0.5 / (NY + 1))   # tip shear, finite-but-mild
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.2)
        ops.algorithm("Newton")
        ops.test("NormDispIncr", 1e-9, 50)
        ops.analysis("Static")
        assert ops.analyze(5) == 0, f"cantilever diverged (nu={nu}, cst={cst})"
        return abs(ops.nodeDisp(tag(NX, NY // 2), 2))

    ratio_mild = tip(0.3, cst=True) / tip(0.3, cst=False)
    ratio_inc = tip(0.4999, cst=True) / tip(0.4999, cst=False)
    # compressible: the CST mesh is stiffer but same order (~0.5); near-
    # incompressible: the T3 locks volumetrically and the ratio collapses (~0.07).
    assert ratio_mild > 0.3, (
        f"sanity: compressible CST/quad ratio should be same-order, got {ratio_mild:.4f}"
    )
    assert ratio_inc < 0.35 * ratio_mild, (
        f"CST-finite should LOCK near incompressibility (honest baseline): "
        f"ratio(nu=0.4999)={ratio_inc:.4f} vs ratio(nu=0.3)={ratio_mild:.4f}"
    )


# --------------------------------------------------------------------------
# parse / material guards
# --------------------------------------------------------------------------
def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QN.items():
        ops.node(tag, x, y)


def test_finite_refused_for_plane_stress():
    """PlaneStress + -geom finite is refused at the FACTORY (the finite volume
    weight omits the out-of-plane thickness stretch λ). The element must never
    be created, and the interpreter must survive."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    try:
        ops.element("LadrunoCST", 99, 1, 2, 3, 2, "-type", "PlaneStress",
                    "-geom", "finite")
    except Exception:
        pass
    assert 99 not in ops.getEleTags(), "PlaneStress + finite must be refused at parse"
    # a valid PlaneStrain-finite element still constructs afterwards
    ops.element("LadrunoCST", 1, 1, 2, 3, 2, "-type", "PlaneStrain",
                "-geom", "finite")
    assert 1 in ops.getEleTags()


def test_finite_rejects_non_finite_material():
    """-geom finite with a plain (non-FiniteStrainND2D) material must be refused
    at the FACTORY — the element is never created (no bogus assembly, no kernel
    crash) — and the interpreter must survive."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)   # NOT a LogStrain2D
    try:
        ops.element("LadrunoCST", 1, 1, 2, 3, 1, "-type", "PlaneStrain",
                    "-geom", "finite")
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), "non-finite material must be refused at parse"
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoCST", 2, 1, 2, 3, 2, "-type", "PlaneStrain",
                "-geom", "finite")
    assert 2 in ops.getEleTags()


def test_bogus_geom_refused():
    """'-geom corot' (or any non linear|finite token) is refused at parse."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    try:
        ops.element("LadrunoCST", 1, 1, 2, 3, 1, "-type", "PlaneStrain",
                    "-geom", "corot")
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), "-geom corot must be refused at parse"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
