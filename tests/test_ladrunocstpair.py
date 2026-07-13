"""LadrunoCSTPair (ADR 70 P4a) — openseespy end-to-end battery.

The patch MATH is verified building-free by tests/cstpair_reference.py
(FD-exact tangent at finite stress, row-form == dSNPO 15.37/15.38, reduce-to-
2-CSTs, 3-RBM rank). This battery gates what only a build can prove:

  * the ASSEMBLED tangent equals the oracle's — printA vs assemble_pair at a
    genuinely finite, INHOMOGENEOUS state (the crown gate: it catches a
    transposed copy, a symmetrized K, missing/wrong cross blocks, wrong
    v-weights — everything the Newton-converges test alone cannot);
  * homogeneous finite-stretch patch: J-bar = J so the pair reduces exactly to
    two CST-finite elements; GP stresses match the LogStrain2D oracle sigma(F)
    + deformed-config equilibrium;
  * the LOCKING-RELIEF FLIP (the element's reason to exist): on the P2
    battery's near-incompressible cantilever the plain-CST/F-bar-quad tip
    ratio collapses (pinned at ~0.07); the pair must hold O(1);
  * free-pair zero-energy = exactly 3 RBM (T0 gate);
  * large-strain Newton convergence (elastic + J2 inner) with an iteration
    ceiling; det F <= 0 step-cut; Jbar response; parse guards.

Plan: Ladruno_implementation/70_ladruno_plane_finite_triangles_adr.md (§9 P4).
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.fem_checks import assert_zero_energy
import cstpair_reference as cp
import logstrain2d_reference as l2

pytestmark = [pytest.mark.zone_a]

I3 = np.eye(3)


def _KG(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


# mildly distorted CCW quad (same shape as the oracle's REF4)
_QN = {1: (0.0, 0.0), 2: (2.0, 0.1), 3: (2.2, 1.8), 4: (0.1, 1.9)}
_THK = 0.7


def _fresh(E=1000.0, nu=0.3, thick=_THK, coords=_QN):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in coords.items():
        ops.node(tag, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain2D", 2, 1)
    ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, 2, "-thick", thick,
                "-type", "PlaneStrain")


def _static(tol=1e-10, iters=50, dlam=1.0):
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.integrator("LoadControl", dlam)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", tol, iters)
    ops.analysis("Static")


# --------------------------------------------------------------------------
# T0: free pair has exactly 3 rigid-body zero-energy modes
# --------------------------------------------------------------------------
def test_free_pair_3_rbm():
    _fresh()
    assert_zero_energy([1, 2, 3, 4], ndf=2, ndm=2)


# --------------------------------------------------------------------------
# crown gate: assembled tangent == oracle (incl. unsymmetric cross blocks)
# --------------------------------------------------------------------------
def test_assembled_tangent_matches_oracle_at_finite_stress():
    """Drive the pair to a genuinely finite, INHOMOGENEOUS state (unequal nodal
    loads so J1 != J2 and the stress-proportional cross blocks are live), then
    diff the ASSEMBLED free-DOF tangent (printA, FullGeneral) against the numpy
    oracle's exact 8x8 at the solved displacement state. Unsym-aware by
    construction — a transposed row-major copy, a symmetrized K, or wrong
    v_s/v_p weights all fail here."""
    E, nu = 1000.0, 0.3
    _fresh(E, nu)
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 30.0, -12.0)
    ops.load(4, -8.0, 25.0)
    _static(dlam=0.2)
    assert ops.analyze(5) == 0      # -> load factor 1.0

    # full 8-dof displacement vector in node order 1..4
    u = np.zeros(8)
    for n in range(1, 5):
        d = ops.nodeDisp(n)
        u[2 * (n - 1)], u[2 * (n - 1) + 1] = d[0], d[1]
    # confirm the state is inhomogeneous + finite (J1 != J2, strain ~ O(1e-2+))
    assert np.max(np.abs(u)) > 1e-3

    K_, G_ = _KG(E, nu)
    X = np.array([_QN[n] for n in range(1, 5)], dtype=float)
    _, K8 = cp.assemble_pair(u, K_, G_, thickness=_THK, X=X)

    # equation-number map (-1 = constrained)
    eqn = {}
    for n in range(1, 5):
        for comp, q in enumerate(ops.nodeDOFs(n)):
            if q >= 0:
                eqn[q] = 2 * (n - 1) + comp
    nfree = len(eqn)
    # FullGenLinSOE stores A COLUMN-major (addA: A[col*size+row] += m(j,i)) —
    # a C-order reshape here would silently transpose the matrix and pass a
    # transposed element tangent / fail a correct one. order="F" is load-bearing.
    Kasm = np.array(ops.printA("-ret")).reshape(nfree, nfree, order="F")
    Kora = np.zeros((nfree, nfree))
    for r in range(nfree):
        for c in range(nfree):
            Kora[r, c] = K8[eqn[r], eqn[c]]

    scale = np.max(np.abs(Kora))
    err = np.max(np.abs(Kasm - Kora)) / scale
    assert err < 1e-8, f"assembled tangent != oracle (rel {err:.3e})"
    # and the state must genuinely exercise the unsymmetric coupling
    asym = np.max(np.abs(Kora - Kora.T)) / scale
    assert asym > 1e-8, "test state too symmetric to gate the cross blocks"


# --------------------------------------------------------------------------
# homogeneous stretch: J-bar = J -> pair == 2 CST-finite; stress == oracle
# --------------------------------------------------------------------------
def test_homogeneous_stretch_matches_oracle_and_cst():
    E, nu, t = 500.0, 0.3, 1.0
    K_, G_ = _KG(E, nu)
    sq = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}

    def run(place):
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for tag, (x, y) in sq.items():
            ops.node(tag, float(x), float(y))
        ops.nDMaterial("ElasticIsotropic", 1, E, nu)
        ops.nDMaterial("LogStrain2D", 2, 1)
        place()
        ops.fix(1, 1, 1)
        ops.fix(2, 0, 1)
        ops.fix(4, 1, 0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 60.0, 0.0)
        ops.load(3, 60.0, 0.0)
        _static(dlam=0.1)
        assert ops.analyze(10) == 0      # -> load factor 1.0
        return [c for n in sorted(sq) for c in ops.nodeDisp(n)]

    u_pair = run(lambda: ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, 2,
                                     "-thick", t, "-type", "PlaneStrain"))
    a = ops.nodeDisp(2, 1)          # stretch: x -> (1+a) x
    b = ops.nodeDisp(3, 2)          # lateral contraction
    F = np.array([[1.0 + a, 0.0], [0.0, 1.0 + b]])

    # both GP stresses equal the LogStrain2D oracle sigma(F) at the solved F
    ref = l2.plane_strain_step(F, I3.copy(), I3.copy(), "elastic", K=K_, G=G_)
    s_ref = ref["sigma_voigt"]
    s_pair = ops.eleResponse(1, "stresses")
    assert len(s_pair) == 6
    for t_ in range(2):
        for c in range(3):
            assert abs(s_pair[3 * t_ + c] - s_ref[c]) <= 1e-6 * max(1.0, abs(s_ref[c])), \
                f"tri {t_+1} comp {c}: {s_pair[3*t_+c]} vs oracle {s_ref[c]}"

    # deformed-config equilibrium: sigma_xx * (deformed height) * t = total load
    total = s_ref[0] * (1.0 + b) * t
    assert abs(total - 120.0) <= 1e-6 * 120.0

    # Jbar response = det F under the homogeneous state
    jbar = ops.eleResponse(1, "Jbar")[0]
    assert abs(jbar - np.linalg.det(F)) < 1e-9

    # identical to two plain CST-finite elements (J-bar = J -> F-bar = F)
    def place_csts():
        ops.element("LadrunoCST", 1, 1, 2, 3, 2, "-thick", t,
                    "-type", "PlaneStrain", "-geom", "finite")
        ops.element("LadrunoCST", 2, 1, 3, 4, 2, "-thick", t,
                    "-type", "PlaneStrain", "-geom", "finite")
    u_cst = run(place_csts)
    for i, (p, c) in enumerate(zip(u_pair, u_cst)):
        assert abs(p - c) <= 1e-9 + 1e-8 * max(abs(p), abs(c)), \
            f"[{i}] pair {p!r} vs 2xCST {c!r}"


# --------------------------------------------------------------------------
# the reason this element exists: near-incompressible locking-relief FLIP
# --------------------------------------------------------------------------
def _cantilever_tip(nu, kind, nx=8, ny=2):
    """P2-battery cantilever (nx x ny unit cells), tip-loaded; returns the tip
    deflection. kind: 'pair' (one LadrunoCSTPair per cell), 'cst' (two plain
    CST-finite per cell), 'quad' (one F-bar LadrunoQuad per cell)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    def tag(i, j):
        return 1 + j * (nx + 1) + i

    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(tag(i, j), float(i), float(j))
    E = 1000.0
    ops.nDMaterial("ElasticIsotropic", 1, E, nu)
    ops.nDMaterial("LogStrain2D", 2, 1)
    k = 0
    for j in range(ny):
        for i in range(nx):
            n1, n2 = tag(i, j), tag(i + 1, j)
            n3, n4 = tag(i + 1, j + 1), tag(i, j + 1)
            k += 1
            if kind == "pair":
                ops.element("LadrunoCSTPair", k, n1, n2, n3, n4, 2,
                            "-thick", 1.0, "-type", "PlaneStrain")
            elif kind == "cst":
                ops.element("LadrunoCST", 2 * k - 1, n1, n2, n3, 2, "-thick", 1.0,
                            "-type", "PlaneStrain", "-geom", "finite")
                ops.element("LadrunoCST", 2 * k, n1, n3, n4, 2, "-thick", 1.0,
                            "-type", "PlaneStrain", "-geom", "finite")
            elif kind == "quad":
                ops.element("LadrunoQuad", k, n1, n2, n3, n4, 2, "-thick", 1.0,
                            "-type", "PlaneStrain", "-formulation", "bbar",
                            "-geom", "finite")
    for j in range(ny + 1):
        ops.fix(tag(0, j), 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(tag(nx, ny), 0.0, -1.0)
    _static(tol=1e-9)
    assert ops.analyze(1) == 0
    return abs(ops.nodeDisp(tag(nx, ny), 2))


def test_pair_relieves_locking_flip():
    """At nu -> 0.5 the plain-CST/F-bar-quad tip ratio collapses (the P2
    honest-baseline pin sits at ~0.07). The pair must (a) hold an O(1) ratio
    to the F-bar quad (dSNPO Fig 15.8: 'virtually coincide') and (b) beat the
    plain CST by a wide factor — the FLIP of the P2 pin."""
    nu = 0.4999
    tip_pair = _cantilever_tip(nu, "pair")
    tip_cst = _cantilever_tip(nu, "cst")
    tip_quad = _cantilever_tip(nu, "quad")

    ratio_quad = tip_pair / tip_quad
    ratio_cst = tip_pair / tip_cst
    assert ratio_quad > 0.7, \
        f"pair still locks vs F-bar quad at nu={nu}: ratio {ratio_quad:.3f}"
    assert ratio_quad < 1.5, \
        f"pair suspiciously SOFTER than the F-bar quad: {ratio_quad:.3f}"
    assert ratio_cst > 5.0, \
        f"pair barely beats plain CST ({ratio_cst:.2f}x) — patch inactive?"

    # sanity at mild nu: everything agrees loosely
    mild_pair = _cantilever_tip(0.3, "pair")
    mild_quad = _cantilever_tip(0.3, "quad")
    assert 0.7 < mild_pair / mild_quad < 1.3


# --------------------------------------------------------------------------
# Newton convergence at large strain (elastic + J2), with iteration ceiling
# --------------------------------------------------------------------------
def _stretch_column(inner_fn, target=0.2, nstep=40):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    crd = {1: (0, 0), 2: (1, 0), 3: (0, 1), 4: (1, 1), 5: (0, 2), 6: (1, 2)}
    for tag, (x, y) in crd.items():
        ops.node(tag, float(x), float(y))
    inner_fn()                                           # defines 3D inner as tag 1
    ops.nDMaterial("LogStrain2D", 2, 1)
    # two stacked pair cells (CCW quads)
    ops.element("LadrunoCSTPair", 1, 1, 2, 4, 3, 2, "-type", "PlaneStrain")
    ops.element("LadrunoCSTPair", 2, 3, 4, 6, 5, 2, "-type", "PlaneStrain")
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
    ops.test("NormDispIncr", 1e-9, 12)   # iteration CEILING: consistent-tangent gate
    ops.analysis("Static")
    rc = ops.analyze(nstep)
    return rc, ops.nodeDisp(5, 2)


def test_finite_convergence_elastic():
    rc, tip = _stretch_column(
        lambda: ops.nDMaterial("ElasticIsotropic", 1, 800.0, 0.3))
    assert rc == 0, "pair Newton (elastic) failed within the iteration ceiling"
    assert tip > 0.1


def test_finite_convergence_j2():
    def inner():
        try:
            ops.nDMaterial("LadrunoJ2", 1, 800.0, 0.3, 8.0)
        except Exception:
            ops.nDMaterial("J2Plasticity", 1, 400.0, 300.0, 8.0, 12.0, 0.1, 1.0)
    try:
        rc, tip = _stretch_column(inner, target=0.1)
    except Exception as e:
        pytest.skip(f"no 3D J2 material available for the finite inner: {e}")
    assert rc == 0, "pair Newton (J2) failed within the iteration ceiling"
    assert tip > 0.02


# --------------------------------------------------------------------------
# det F <= 0 step-cut: an inverting trial step fails GRACEFULLY
# --------------------------------------------------------------------------
def test_detF_step_cut():
    _fresh()
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, -1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    # one huge downward step drives triangle inversion
    ops.integrator("DisplacementControl", 3, 2, -5.0)
    ops.algorithm("Newton")
    ops.test("NormDispIncr", 1e-9, 10)
    ops.analysis("Static")
    rc = ops.analyze(1)
    assert rc != 0, "an element-inverting step must fail, not 'converge'"
    # interpreter survives: a fresh small model still runs
    _fresh()
    ops.fix(1, 1, 1)
    ops.fix(2, 0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 1.0, 0.0)
    _static()
    assert ops.analyze(1) == 0


# --------------------------------------------------------------------------
# parse guards
# --------------------------------------------------------------------------
def _guard_model():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for tag, (x, y) in _QN.items():
        ops.node(tag, float(x), float(y))


def test_plane_stress_refused():
    _guard_model()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    with pytest.raises(Exception):
        ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, 2, "-type", "PlaneStress")
    assert ops.getEleTags() in ([], None)


def test_non_finite_material_refused():
    _guard_model()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    with pytest.raises(Exception):
        ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, 1)
    assert ops.getEleTags() in ([], None)


def test_geom_linear_refused():
    _guard_model()
    ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.3)
    ops.nDMaterial("LogStrain2D", 2, 1)
    with pytest.raises(Exception):
        ops.element("LadrunoCSTPair", 1, 1, 2, 3, 4, 2, "-geom", "linear")
    assert ops.getEleTags() in ([], None)


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))
