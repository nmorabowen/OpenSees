"""LadrunoLST -geom finite (ADR 70 P3) — openseespy end-to-end battery.

The finite-strain *mechanics* are verified building-free by the P0 oracle
(tests/finitestrain2d_reference.py — it carries a full T6 std lane: FD
consistent tangent, homogeneous patch, rigid-rotation objectivity) and the C++
kernel diff. This battery gates the ELEMENT wiring: parser, the setTrialF path
at the 3 interior points, reduce-to-linear, and large-strain Newton
convergence. std ONLY — the centroid F-bar was REFUTED on the T6 at P3
(rank-deficient: 2 conformal zero-energy modes; see
test_ladrunolst_element.py's refutation pin and ADR 70 §9), so bbar (and with
it any F-bar-finite lane) is parser-refused.

Gates (ADR 70 §6, P3 finite half, amended by the refutation):
  * homogeneous finite-stretch patch on a 2-LST unit square: every integration
    point's Cauchy stress == the LogStrain2D oracle σ(F) at the solved F, plus
    the deformed-config equilibrium check σxx·(1+b)·t = applied Fx;
  * reduce-to-linear at ε→0;
  * large-strain Newton convergence over elastic and a 3D J2 inner;
  * std-finite Newton iteration-count ceiling (a transposed K still converges —
    linearly; the ceiling catches it);
  * det F ≤ 0 step-cut (graceful failure, interpreter survives);
  * parse guards (PlaneStress+finite, non-finite material, bbar+finite).

Plan: Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md.
"""
import math

import numpy as np
import pytest

from _testbed import ops
import logstrain2d_reference as l2

pytestmark = [pytest.mark.zone_a]

I3 = np.eye(3)


def _KG(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


def _norm(v):
    return math.sqrt(sum(x * x for x in v))


# --------------------------------------------------------------------------
# meshes
# --------------------------------------------------------------------------
# distorted single T6 (straight edges, midsides at midpoints)
_C = {1: (0.00, 0.00), 2: (1.00, 0.10), 3: (0.15, 0.95)}
_T6N = {
    1: _C[1], 2: _C[2], 3: _C[3],
    4: ((_C[1][0] + _C[2][0]) / 2, (_C[1][1] + _C[2][1]) / 2),
    5: ((_C[2][0] + _C[3][0]) / 2, (_C[2][1] + _C[3][1]) / 2),
    6: ((_C[3][0] + _C[1][0]) / 2, (_C[3][1] + _C[1][1]) / 2),
}
_CONN = [1, 2, 3, 4, 5, 6]
_BASE = [1, 2, 4]                 # bottom edge fixed
_LOADS = {3: (6.0, -5.0), 5: (2.0, 3.0), 6: (-1.5, 4.0)}
_THK = 0.7

# unit square tiled by TWO T6s (diagonal 1-3); nodes 5..9 = edge/diagonal mids
_SQN = {
    1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0),
    5: (0.5, 0.0),                     # mid bottom (1-2)
    6: (1.0, 0.5),                     # mid right  (2-3)
    7: (0.5, 0.5),                     # mid diagonal (1-3)
    8: (0.5, 1.0),                     # mid top    (3-4)
    9: (0.0, 0.5),                     # mid left   (4-1)
}
_SQE = [
    (1, 2, 3, 5, 6, 7),                # corners CCW; mids (1-2),(2-3),(3-1)
    (1, 3, 4, 7, 8, 9),                # corners CCW; mids (1-3),(3-4),(4-1)
]


def _place_sq(matTag, geom="finite", thick=1.0):
    for k, conn in enumerate(_SQE):
        ops.element("LadrunoLST", k + 1, *conn, matTag, "-thick", thick,
                    "-type", "PlaneStrain", "-geom", geom)


# --------------------------------------------------------------------------
# reduce-to-linear: at tiny strain -geom finite collapses to -geom linear
# --------------------------------------------------------------------------
def _solve_single(place_fn, load_scale):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    place_fn()
    for n in _BASE:
        ops.fix(n, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n, (fx, fy) in _LOADS.items():
        ops.load(n, fx * load_scale, fy * load_scale)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-12, 30)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    disps = []
    for n in _CONN:
        disps.extend(ops.nodeDisp(n))
    return disps


def test_reduce_to_linear():
    """At ε→0 the finite path (LogStrain2D/elastic) matches the small-strain
    -geom linear element to ~1e-5 relative — the setTrialF/kernel wiring
    collapses to linear elasticity."""
    scale = 1e-6
    lin = _solve_single(
        lambda: ops.element("LadrunoLST", 1, *_CONN, 1, "-thick", _THK,
                            "-type", "PlaneStrain"),
        scale)
    fin = _solve_single(
        lambda: (ops.nDMaterial("LogStrain2D", 2, 1),
                 ops.element("LadrunoLST", 1, *_CONN, 2, "-thick", _THK,
                             "-type", "PlaneStrain", "-geom", "finite")),
        scale)
    for i, (a, b) in enumerate(zip(fin, lin)):
        tol = 1e-12 + 1e-5 * max(abs(a), abs(b))
        assert abs(a - b) <= tol, f"[{i}] finite {a!r} vs linear {b!r}"


# --------------------------------------------------------------------------
# homogeneous finite-stretch patch: every IP's Cauchy stress == oracle σ(F)
# --------------------------------------------------------------------------
def test_homogeneous_stretch_matches_oracle():
    """Two-LST unit square, symmetry-fixed, pulled uniaxially into a genuinely
    finite stretch by the CONSISTENT quadratic-edge loads (corner 1/6, midside
    2/3 of the resultant). The homogeneous deformation is exactly representable
    by the quadratic T6 field, so all 6 integration points must report the same
    Cauchy stress == the LogStrain2D numpy oracle σ(F) at the solved F."""
    E, nu, t = 500.0, 0.3, 1.0
    K, G = _KG(E, nu)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _SQN.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain2D", 2, 1)
    _place_sq(2, thick=t)
    ops.fix(1, 1, 1)                   # origin pinned
    ops.fix(2, 0, 1)                   # bottom edge: v=0
    ops.fix(5, 0, 1)
    ops.fix(4, 1, 0)                   # left edge: u=0
    ops.fix(9, 1, 0)
    # uniform traction resultant Fx=120 on the quadratic right edge (2-6-3):
    # consistent loads 1/6, 2/3, 1/6
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 20.0, 0.0)
    ops.load(6, 80.0, 0.0)
    ops.load(3, 20.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 0.1)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-11, 50)
    ops.analysis("Static")
    assert ops.analyze(10) == 0

    a = ops.nodeDisp(2, 1)             # u_x at X=1 → F_xx − 1
    b = ops.nodeDisp(4, 2)             # u_y at Y=1 → F_yy − 1
    F = np.array([[1.0 + a, 0.0], [0.0, 1.0 + b]])
    assert a > 0.1, f"test should reach a genuinely finite stretch (a={a})"

    # deformed-config equilibrium: σxx·(1+b)·t must balance the applied Fx=120
    applied_fx = 120.0
    sxx_e1 = list(ops.eleResponse(1, "stresses"))[0]
    deformed_face = (1.0 + b) * t
    assert abs(sxx_e1 * deformed_face - applied_fx) <= 1e-3 * applied_fx, (
        f"external equilibrium: σxx·(1+b)·t = {sxx_e1 * deformed_face:.4f} "
        f"!= applied Fx = {applied_fx} (a J volume weight is likely missing)"
    )

    res = l2.plane_strain_step(F, I3.copy(), I3.copy(), "elastic", K=K, G=G)
    s_ref = res["sigma_voigt"]
    smag = max(abs(v) for v in s_ref)

    for ele in (1, 2):
        stresses = list(ops.eleResponse(ele, "stresses"))   # 3 IP × [sxx,syy,sxy]
        for gp in range(3):
            sxx, syy, sxy = stresses[3 * gp:3 * gp + 3]
            tol = 1e-5 * smag + 1e-7
            assert abs(sxx - s_ref[0]) <= tol, f"ele{ele} GP{gp} sxx {sxx} vs {s_ref[0]}"
            assert abs(syy - s_ref[1]) <= tol, f"ele{ele} GP{gp} syy {syy} vs {s_ref[1]}"
            assert abs(sxy - s_ref[2]) <= tol, f"ele{ele} GP{gp} sxy {sxy} vs {s_ref[2]}"
    assert abs(s_ref[1]) < 1e-4 * smag + 1e-4    # traction-free top ⇒ σyy ≈ 0


# --------------------------------------------------------------------------
# large-strain Newton convergence (end-to-end), elastic + J2 inner
# --------------------------------------------------------------------------
def _stretch_square(inner_fn, target=0.2, nstep=40):
    """The 2-LST unit square pulled vertically under DISPLACEMENT control on
    the top-mid node (monotonic ⇒ robust for plasticity)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _SQN.items():
        ops.node(tag, float(x), float(y))
    inner_fn()                                   # defines 3D inner as tag 1
    ops.nDMaterial("LogStrain2D", 2, 1)
    _place_sq(2)
    for n in (1, 5, 2):                          # bottom edge fixed
        ops.fix(n, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 0.0, 1.0 / 6.0)                  # consistent quadratic-edge shape
    ops.load(8, 0.0, 2.0 / 3.0)
    ops.load(3, 0.0, 1.0 / 6.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.integrator("DisplacementControl", 8, 2, target / nstep)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-9, 50)
    ops.analysis("Static")
    rc = ops.analyze(nstep)
    tip = ops.nodeDisp(8, 2)
    return rc, tip


def test_finite_convergence_elastic():
    rc, tip = _stretch_square(
        lambda: ops.nDMaterial("ElasticIsotropic", 1, 800.0, 0.3))
    assert rc == 0, "finite-strain Newton (elastic) failed to converge"
    assert tip > 0.1, f"expected a genuinely finite tensile tip disp, got {tip}"


def test_finite_convergence_j2():
    """LogStrain2D over a genuinely elasto-plastic 3D inner (LadrunoJ2 if
    present, else upstream J2Plasticity)."""
    def inner():
        try:
            ops.nDMaterial("LadrunoJ2", 1, 800.0, 0.3, 8.0)
        except Exception:
            ops.nDMaterial("J2Plasticity", 1, 400.0, 300.0, 8.0, 12.0, 0.1, 1.0)
    try:
        rc, tip = _stretch_square(inner, target=0.1)
    except Exception as e:
        pytest.skip(f"no 3D J2 material available for the finite inner: {e}")
    assert rc == 0, "finite-strain Newton (J2) failed to converge"
    assert tip > 0.02


# --------------------------------------------------------------------------
# Newton iteration-count gate (consistent spatial tangent)
# --------------------------------------------------------------------------
_NEWTON_ITER_CEIL = 12  # calibrated after build: healthy quadratic solve is
#                         a handful of iters on this step; gate = ~2x.


def test_std_finite_newton_iterations_gate():
    """One big load step, std + -geom finite. A transposed or inconsistent K
    still converges (the residual is exact) but at a LINEAR rate — the
    iteration ceiling catches what a converged-answer check can't."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoLST", 1, *_CONN, 2, "-thick", _THK,
                "-type", "PlaneStrain", "-geom", "finite")
    for n in _BASE:
        ops.fix(n, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 20.0, -8.0)                       # large -> genuinely finite
    ops.load(5, 7.0, 11.0)
    ops.load(6, -5.0, 13.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)            # ONE big step
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 60)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "std-finite single big step failed to converge"
    iters = ops.testIter()
    assert _norm(ops.nodeDisp(3)) > 0.05, "step did not reach a finite deformation"
    assert iters <= _NEWTON_ITER_CEIL, (
        f"std-finite Newton took {iters} iters (> {_NEWTON_ITER_CEIL}); an "
        "inconsistent/transposed tangent degrades the rate"
    )


# --------------------------------------------------------------------------
# det F ≤ 0 step-cut: an inverting trial step fails GRACEFULLY
# --------------------------------------------------------------------------
def test_detF_step_cut():
    """One huge displacement-controlled step turns the square inside out
    (det F ≤ 0). setTrialF must reject the trial → analyze reports failure —
    no crash, no silently-accepted inverted state — and the interpreter
    survives to run a valid model afterwards."""
    def build(load_dir):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for tag, (x, y) in _SQN.items():
            ops.node(tag, x, y)
        ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
        ops.nDMaterial("LogStrain2D", 2, 1)
        _place_sq(2)
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.fix(5, 0, 1)
        ops.fix(4, 1, 0)
        ops.fix(9, 1, 0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(3, load_dir, 0.0)
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Transformation")

    build(-1.0)
    ops.integrator("DisplacementControl", 3, 1, -1.5)   # corner across the mesh
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 10)
    ops.analysis("Static")
    rc = ops.analyze(1)
    assert rc != 0, "an element-inverting step must FAIL (det F step-cut), not pass"

    build(1.0)                                          # interpreter survives
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-10, 25)
    ops.analysis("Static")
    assert ops.analyze(1) == 0


# --------------------------------------------------------------------------
# serialization: database round-trip preserves the FINITE element
# --------------------------------------------------------------------------
def test_finite_database_roundtrip():
    """sendSelf/recvSelf + broker reconstruction: solve a finite-stretch state,
    round-trip through the FE_Datastore, and probe the element-owned resisting
    force — under -geom finite it is assembled through the spatial/kernel path,
    so a recvSelf that silently dropped the geom int (reverting to linear)
    changes it at this deformed state."""
    from _testbed.roundtrip import database_roundtrip

    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for tag, (x, y) in _T6N.items():
            ops.node(tag, x, y)
        ops.nDMaterial("ElasticIsotropic", 1, 800.0, 0.3)
        ops.nDMaterial("LogStrain2D", 2, 1)
        ops.element("LadrunoLST", 1, *_CONN, 2, "-thick", _THK,
                    "-type", "PlaneStrain", "-geom", "finite")
        for n in _BASE:
            ops.fix(n, 1, 1)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for n, (fx, fy) in _LOADS.items():
            ops.load(n, 2.0 * fx, 2.0 * fy)          # genuinely finite state
        ops.system("FullGeneral")
        ops.numberer("Plain")
        ops.constraints("Plain")
        ops.integrator("LoadControl", 0.25)
        ops.algorithm("Newton")
        ops.test("NormDispIncr", 1e-10, 50)
        ops.analysis("Static")
        assert ops.analyze(4) == 0

    build()
    assert _norm(ops.nodeDisp(3)) > 1e-2, "roundtrip state would be vacuous"
    database_roundtrip(build, probe_nodes=[3, 5, 6], ndf=2,
                       probe_fn=lambda: ops.eleResponse(1, "forces"))


# --------------------------------------------------------------------------
# parse / material guards
# --------------------------------------------------------------------------
def _fresh():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _T6N.items():
        ops.node(tag, x, y)


def test_finite_refused_for_plane_stress():
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    try:
        ops.element("LadrunoLST", 99, *_CONN, 2, "-type", "PlaneStress",
                    "-geom", "finite")
    except Exception:
        pass
    assert 99 not in ops.getEleTags(), "PlaneStress + finite must be refused at parse"
    ops.element("LadrunoLST", 1, *_CONN, 2, "-type", "PlaneStrain",
                "-geom", "finite")
    assert 1 in ops.getEleTags()


def test_finite_rejects_non_finite_material():
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)   # NOT a LogStrain2D
    try:
        ops.element("LadrunoLST", 1, *_CONN, 1, "-type", "PlaneStrain",
                    "-geom", "finite")
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), "non-finite material must be refused at parse"
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoLST", 2, *_CONN, 2, "-type", "PlaneStrain",
                "-geom", "finite")
    assert 2 in ops.getEleTags()


def test_bbar_finite_refused():
    """bbar (and with it F-bar-finite) stays parser-refused on the LST — the
    P3 rank refutation (see test_ladrunolst_element.py)."""
    _fresh()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    try:
        ops.element("LadrunoLST", 1, *_CONN, 2, "-type", "PlaneStrain",
                    "-formulation", "bbar", "-geom", "finite")
    except Exception:
        pass
    assert 1 not in ops.getEleTags(), "bbar+finite must stay refused (rank)"


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
