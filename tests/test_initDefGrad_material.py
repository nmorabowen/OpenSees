"""InitDefGradNDMaterial — multiplicative staged-activation wrapper battery.

InitDefGrad is the finite-strain sibling of InitStrain: it makes a finite-strain
continuum element BORN STRESS-FREE at the deformed geometry when appended
mid-stage. At activation it captures the birth deformation gradient F0 and
thereafter feeds the inner finite-strain material the RELATIVE gradient
F_rel = F·F0⁻¹ (multiplicative split F = F_rel·F0), so it is objective under a
post-birth rigid rotation where the additive InitStrain offset is not.

It wraps an inner FiniteStrainNDMaterial (LogStrain) and IS-A FiniteStrainNDMaterial,
so LadrunoBrick -geom finite consumes it with zero element edits. The material has
no Python setTrialF entry point, so every test drives it through a single
LadrunoBrick -geom finite element by imposing the full 24-dof displacement field.

Coverage:
  * -noInitF (pass-through)  ≡ bare LogStrain  (stress AND force);
  * explicit -F0 = G, impose F = G            -> stress-free birth;
  * explicit -F0 = G, impose F = Q·G          -> objectivity (zero stress under
                                                  post-birth rigid rotation);
  * InitDefGrad(-F0=G) at F = M  Cauchy σ  ==  bare LogStrain at F = M·G⁻¹
                                                  (reduce-to-relative);
  * consistent tangent with -F0 = G == FD of the resisting force (the wrapper must
    not break Newton — the J-in-spatial-tangent consistency arbiter);
  * genuine TWO-PHASE staged construction: deform+commit a holder, THEN append an
    auto-capture InitDefGrad brick on the same nodes -> it is born with ~zero force.

Design note: Ladruno_implementation/staged_deformation_gradiend.md.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# A distorted positive-Jacobian reference hex (same shape as the finite battery).
_NODES = {
    1: (0.00, 0.00, 0.00),
    2: (1.00, 0.10, 0.05),
    3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10),
    5: (0.00, 0.05, 1.00),
    6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10),
    8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_E = 200.0
_NU = 0.3

# A representative finite birth gradient G (stretch + shear, det > 0).
_G = np.array([[1.12, 0.06, -0.04],
               [0.05, 0.94, 0.07],
               [-0.02, 0.03, 1.08]])


def _f0_args(G):
    """Row-major 9-tuple for the -F0 option."""
    return ["-F0", *[float(v) for v in np.asarray(G).reshape(9)]]


def _build(mat_spec, eltag=1, mbase=100):
    """Add a LadrunoBrick -geom finite with the material described by mat_spec.

    mat_spec is one of:
      ("bare",)                 -> nDMaterial LogStrain over ElasticIsotropic
      ("noinit",)               -> InitDefGrad -noInitF over LogStrain (pass-through)
      ("f0", G)                 -> InitDefGrad -F0 G over LogStrain
      ("auto",)                 -> InitDefGrad (auto-capture) over LogStrain
    eltag is the element tag the test will query; mbase offsets the material tags
    so two coexisting elements (the staged-construction test) never collide.
    Returns the element tag.
    """
    ela = mbase + 1
    log = mbase + 2
    wrap = mbase + 3
    ops.nDMaterial("ElasticIsotropic", ela, _E, _NU)
    ops.nDMaterial("LogStrain", log, ela)
    kind = mat_spec[0]
    if kind == "bare":
        mtag = log
    elif kind == "noinit":
        ops.nDMaterial("InitDefGrad", wrap, log, "-noInitF")
        mtag = wrap
    elif kind == "f0":
        ops.nDMaterial("InitDefGrad", wrap, log, *_f0_args(mat_spec[1]))
        mtag = wrap
    elif kind == "auto":
        ops.nDMaterial("InitDefGrad", wrap, log)
        mtag = wrap
    else:
        raise ValueError(kind)
    ops.element("LadrunoBrick", eltag, *_CONN, mtag, "-formulation", "std",
                "-geom", "finite")
    return eltag


def _affine_disp(F):
    """Nodal displacements u_a = (F - I) X_a for a uniform deformation gradient."""
    u = np.zeros(24)
    I = np.eye(3)
    F = np.asarray(F, dtype=float)
    for tag, (x, y, z) in _NODES.items():
        X = np.array([x, y, z])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = (F - I) @ X
    return u


def _setup_nodes():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)


def _impose_and_solve(u, mat_spec):
    """Single-element model: impose the full 24-dof field u and solve one step."""
    _setup_nodes()
    _build(mat_spec)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        b = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[b + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _stresses(u, mat_spec):
    assert _impose_and_solve(u, mat_spec) == 0, "imposed-displacement solve failed"
    s = ops.eleResponse(1, "stresses")
    assert s and len(s) == 48, "LadrunoBrick 'stresses' response missing (8 GP x 6)"
    return np.array(s, dtype=float).reshape(8, 6)


def _force(u, mat_spec):
    assert _impose_and_solve(u, mat_spec) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _stiff(u, mat_spec):
    assert _impose_and_solve(u, mat_spec) == 0, "imposed-displacement solve failed"
    K = ops.eleResponse(1, "stiff")
    assert K and len(K) == 24 * 24, "LadrunoBrick must answer eleResponse('stiff') 24x24"
    return np.array(K, dtype=float).reshape(24, 24)


# --------------------------------------------------------------------------- #
#  -noInitF is a pure pass-through: identical to the bare inner.              #
# --------------------------------------------------------------------------- #
def test_noinit_is_passthrough():
    M = [[1.10, 0.05, -0.02],
         [0.03, 0.96, 0.04],
         [-0.01, 0.02, 1.07]]
    u = _affine_disp(M)
    s_bare = _stresses(u, ("bare",))
    s_wrap = _stresses(u, ("noinit",))
    f_bare = _force(u, ("bare",))
    f_wrap = _force(u, ("noinit",))
    sc = max(np.abs(s_bare).max(), 1.0)
    fc = max(np.abs(f_bare).max(), 1.0)
    assert np.abs(s_wrap - s_bare).max() <= 1.0e-10 * sc, "-noInitF changed the stress"
    assert np.abs(f_wrap - f_bare).max() <= 1.0e-10 * fc, "-noInitF changed the force"


# --------------------------------------------------------------------------- #
#  Stress-free birth: -F0 = G, impose F = G  ->  F_rel = I  ->  zero stress.   #
# --------------------------------------------------------------------------- #
def test_explicit_F0_births_stress_free():
    u = _affine_disp(_G)                       # F = G at every GP
    s = _stresses(u, ("f0", _G))
    assert np.abs(s).max() <= 1.0e-7 * _E, (
        f"born stress-free expected; got Cauchy stress {np.abs(s).max():.3e}"
    )
    f = _force(u, ("f0", _G))
    assert np.abs(f).max() <= 1.0e-7 * _E, f"nonzero birth force {np.abs(f).max():.3e}"


# --------------------------------------------------------------------------- #
#  Objectivity: -F0 = G, impose F = Q·G (post-birth rigid rotation) -> 0.      #
#  This is exactly what the additive InitStrain offset CANNOT do.             #
# --------------------------------------------------------------------------- #
def test_post_birth_rigid_rotation_is_stress_free():
    th = 0.4
    cz, sz = math.cos(th), math.sin(th)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(0.3), math.sin(0.3)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    Q = Ry @ Rz
    u = _affine_disp(Q @ _G)                   # F = Q·G  ->  F_rel = Q  ->  C_rel = I
    s = _stresses(u, ("f0", _G))
    assert np.abs(s).max() <= 1.0e-6 * _E, (
        f"post-birth rotation induced stress {np.abs(s).max():.3e} (not objective)"
    )


# --------------------------------------------------------------------------- #
#  Reduce-to-relative: InitDefGrad(-F0=G) at F = M gives the same Cauchy σ as  #
#  the bare inner at F = M·G⁻¹ (both have F_rel = M·G⁻¹).                       #
# --------------------------------------------------------------------------- #
def test_reduces_to_relative_gradient():
    M = np.array([[1.18, 0.09, -0.05],
                  [0.04, 0.90, 0.08],
                  [-0.03, 0.06, 1.13]])
    s_wrap = _stresses(_affine_disp(M), ("f0", _G))
    s_bare = _stresses(_affine_disp(M @ np.linalg.inv(_G)), ("bare",))
    sc = max(np.abs(s_bare).max(), 1.0)
    assert np.abs(s_wrap - s_bare).max() <= 1.0e-6 * sc, (
        "InitDefGrad(-F0=G) at F=M must match bare inner at F=M·G⁻¹ "
        f"(max dev {np.abs(s_wrap - s_bare).max():.3e})"
    )


# --------------------------------------------------------------------------- #
#  The wrapper must preserve the consistent element tangent (Newton).         #
# --------------------------------------------------------------------------- #
def test_wrapped_consistent_tangent_matches_finite_difference():
    M = [[1.13, 0.07, -0.04],
         [0.05, 0.92, 0.06],
         [-0.02, 0.04, 1.09]]
    u = _affine_disp(M)
    spec = ("f0", _G)

    K = _stiff(u, spec)
    assert np.isfinite(K).all()

    h = 1.0e-6
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_force(up, spec) - _force(um, spec)) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    assert err <= 1.0e-4 * scale, (
        f"wrapped tangent != FD of force: max abs err {err:.3e} vs scale "
        f"{scale:.3e} (the wrapper broke the consistent tangent)"
    )


# --------------------------------------------------------------------------- #
#  Genuine staged construction: deform+commit a holder element, THEN append    #
#  an auto-capture InitDefGrad brick on the same nodes. It captures F0 at the   #
#  deformed birth state and is born carrying ~zero force, while the holder      #
#  still carries the full stress.                                              #
# --------------------------------------------------------------------------- #
def test_two_phase_auto_capture_births_stress_free():
    Fdef = [[1.10, 0.04, -0.03],
            [0.03, 0.95, 0.05],
            [-0.02, 0.03, 1.06]]
    u = _affine_disp(Fdef)

    _setup_nodes()
    # Phase 1 — holder brick (bare LogStrain), imposed deformation, commit.
    _build(("bare",), eltag=1, mbase=110)      # element tag 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        b = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[b + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "phase-1 holder solve failed"

    f_holder = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f_holder).max() > 1.0e-3 * _E, "holder should be stressed after phase 1"

    # Phase 2 — append the auto-capture InitDefGrad brick on the SAME nodes; add
    # zero further load so the nodes stay at the committed deformed state. The new
    # element's first setTrialF is at F = Fdef -> it auto-captures F0 = Fdef.
    eltag2 = _build(("auto",), eltag=2, mbase=120)   # element tag 2, on the same _CONN
    ops.integrator("LoadControl", 0.0)         # no extra load -> loadFactor stays 1
    assert ops.analyze(1) == 0, "phase-2 staged-append solve failed"

    f_new = np.array(ops.eleForce(eltag2), dtype=float)
    assert np.abs(f_new).max() <= 1.0e-6 * _E, (
        f"appended element not born stress-free: force {np.abs(f_new).max():.3e}"
    )
    # the holder is unchanged (still carrying the load)
    f_holder2 = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f_holder2 - f_holder).max() <= 1.0e-6 * _E, (
        "holder force changed when the stress-free element was appended"
    )

    # And the captured birth gradient is reported back as Fdef (row-major).
    F0resp = ops.eleResponse(eltag2, "material", "1", "F0")
    if F0resp and len(F0resp) == 9:
        assert np.allclose(np.array(F0resp, dtype=float),
                           np.array(Fdef).reshape(9), atol=1.0e-8), (
            "reported F0 != the deformed birth gradient"
        )
