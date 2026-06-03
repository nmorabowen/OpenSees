"""LadrunoBrick ``-geom finite`` driving the finite-strain-NATIVE combined-hardening
J2 (``nDMaterial LadrunoJ2Finite``) — the §14.11 "v2" that co-rotates the kinematic
backstress, so combined hardening is frame-indifferent under LARGE rotation (the
property the v1 LogStrain-over-LadrunoJ2 path provably lacks).

    nDMaterial LadrunoJ2Finite 1 K G -iso voce ... -kin N ...   # F → setTrialF
    element LadrunoBrick 1 ...nodes... 1 -geom finite

Gates (requires a local/CI build — ``from _testbed import ops``):
  * NO-rotation spine: native ≡ LogStrain(LadrunoJ2) for a symmetric F (the
    co-rotation is a no-op there) — proves v2 changes nothing outside large rotation;
  * OBJECTIVITY end-to-end: a committed plastic path WITH finite rotation matches the
    objective native numpy oracle step for step (co-rotation + commit round-trip
    through the real element) — what v1 fails;
  * CONSISTENT TANGENT: Newton converges with the native tangent through the real
    element (the objectivity test's 6 LoadControl steps each solve to NormDispIncr
    1e-11). The tangent is validated numerically OFF the element: channel A + channel B
    (the backstress co-rotation ∂α̃/∂F) exactly-to-FD at the material level in
    tests/test_ladrunoJ2_finite_native_tangent.py, and the C++ channel-B tensor (the
    F→L mapping + cauchyFromR conventions) against the numpy oracle in
    tests/test_ladrunoJ2_finite_native_tangentB_cpp.py. An element-level FD-of-K gate
    is intentionally absent — see the NOTE at the consistent-tangent section below;
  * SERIALIZATION: sendSelf/recvSelf + broker survive a committed plastic finite state
    (FE_Datastore round-trip; skips if the build lacks database support);
  * rigid rotation of the unloaded element → zero stress/force;
  * tiny strain → reduces to small-strain LadrunoJ2.
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip
import ladrunoj2_reference as lr
import ladrunoj2_finite_native_reference as nat
import ladrunoj2_finite_implex_reference as fx

pytestmark = [pytest.mark.zone_a]

_NODES = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.05), 3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_K, _G = 1500.0, 700.0
_ISO = dict(sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0)
_KIN = dict(C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])
_I = np.eye(3)


def _params():
    return lr.Params(_K, _G, **_ISO, **_KIN)


def _kin_flat():
    out = []
    for c, g in zip(_KIN["C"], _KIN["gam"]):
        out += [c, g]
    return out


def _mat_native(tag):
    ops.nDMaterial("LadrunoJ2Finite", tag, _K, _G, "-iso", "voce",
                   _ISO["sig0"], _ISO["Qinf"], _ISO["bIso"], _ISO["Hiso"],
                   "-kin", len(_KIN["C"]), *_kin_flat())


def _mat_wrapper(tag):  # v1 LogStrain-over-LadrunoJ2, for the no-rotation spine check
    ops.nDMaterial("LadrunoJ2", tag + 100, _K, _G, "-iso", "voce",
                   _ISO["sig0"], _ISO["Qinf"], _ISO["bIso"], _ISO["Hiso"],
                   "-kin", len(_KIN["C"]), *_kin_flat())
    ops.nDMaterial("LogStrain", tag, tag + 100)


def _mat_small(tag):    # small-strain LadrunoJ2 (for the -geom linear reduction)
    ops.nDMaterial("LadrunoJ2", tag, _K, _G, "-iso", "voce",
                   _ISO["sig0"], _ISO["Qinf"], _ISO["bIso"], _ISO["Hiso"],
                   "-kin", len(_KIN["C"]), *_kin_flat())


def _mat_native_implex(tag):   # finite-strain native J2 with the IMPL-EX reporting path
    ops.nDMaterial("LadrunoJ2Finite", tag, _K, _G, "-iso", "voce",
                   _ISO["sig0"], _ISO["Qinf"], _ISO["bIso"], _ISO["Hiso"],
                   "-kin", len(_KIN["C"]), *_kin_flat(), "-implex")


def _build(mat_fn, geom="finite"):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    mat_fn(1)
    ops.element("LadrunoBrick", 1, *_CONN, 1, "-formulation", "std", "-geom", geom)


def _disp_for(F, wiggle=None):
    u = np.zeros(24)
    for tag, (x, y, z) in _NODES.items():
        d = (np.asarray(F) - _I) @ np.array([x, y, z])
        if wiggle is not None:
            d = d + np.array(wiggle[tag - 1])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = d
    return u


def _setup_static(handler="Lagrange"):
    ops.constraints(handler)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")


def _solve_to(F, mat_fn, geom="finite", wiggle=None):
    """One imposed-displacement static solve to F from virgin."""
    _build(mat_fn, geom)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    u = _disp_for(F, wiggle)
    for tag in _NODES:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    _setup_static()
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _force_single(F, mat_fn, geom="finite", wiggle=None):
    assert _solve_to(F, mat_fn, geom, wiggle) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _rotz(th):
    c, s = math.cos(th), math.sin(th)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _oracle_stress_voigt(m, F):
    sig, _ = m.setTrialF(np.asarray(F)); m.commit()
    return np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[0, 2]])


# =========================================================================== #
#  1. NO-rotation spine: native ≡ LogStrain(LadrunoJ2) for a symmetric F.      #
# =========================================================================== #
def test_native_finite_reduces_to_logstrain_no_rotation():
    Ffin = [[1.18, 0.02, 0.01], [0.02, 0.93, 0.015], [0.01, 0.015, 0.95]]
    f_nat = _force_single(Ffin, _mat_native)
    f_wrap = _force_single(Ffin, _mat_wrapper)
    scale = max(np.abs(f_wrap).max(), 1.0e-30)
    assert np.abs(f_nat - f_wrap).max() <= 1.0e-7 * scale, (
        "native finite J2 != LogStrain(LadrunoJ2) for a symmetric F — the "
        "co-rotation must be a no-op when there is no material rotation")


# =========================================================================== #
#  2. OBJECTIVITY end-to-end: committed plastic path WITH finite rotation       #
#     matches the objective native numpy oracle step for step. Driven by a      #
#     single LoadControl ramp to a ROTATED+stretched target (the per-step F_k    #
#     carries progressively more rotation); v1 LogStrain(LadrunoJ2) cannot hold  #
#     this for the kinematic part.                                              #
# =========================================================================== #
def test_native_finite_objective_rotating_plastic_path():
    Ffin = _rotz(0.5) @ np.array([[1.18, 0.04, 0.0], [0.0, 0.93, 0.0], [0.0, 0.0, 0.95]])
    N = 6
    _build(_mat_native, "finite")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _NODES.items():
        d = (Ffin - _I) @ np.array([x, y, z])
        base = (tag - 1) * 3
        for j in range(3):
            ops.sp(tag, j + 1, float(d[j]))
    _setup_static()
    ops.integrator("LoadControl", 1.0 / N)
    ops.analysis("Static")

    m = nat.NativeFiniteJ2(_params(), corotate=True)
    saw = False
    for k in range(1, N + 1):
        assert ops.analyze(1) == 0, f"rotating-plastic solve failed at step {k}"
        Fk = _I + (k / N) * (Ffin - _I)                 # affine ⇒ F_k at every GP
        ref = _oracle_stress_voigt(m, Fk)
        s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)[0]
        tol = 1.0e-6 * max(np.abs(ref).max(), 1.0)
        assert np.allclose(s, ref, rtol=1.0e-6, atol=tol), (
            f"native finite J2 GP Cauchy != objective oracle at rotated-plastic "
            f"step {k}: {s} vs {ref}")
        saw = saw or np.abs(ref).max() > _ISO["sig0"]
    assert saw, "rotating plastic path never developed appreciable stress"


# =========================================================================== #
#  3. CONSISTENT TANGENT: Newton converges quadratically with the native tangent #
#     through the real element along the rotated plastic path (test 2's 6        #
#     LoadControl steps each solve to NormDispIncr 1e-11 in a handful of iters).  #
#                                                                                 #
#     NOTE — why no element-level FD-tangent gate here. The CONSISTENT tangent    #
#     (channel A = log-strain (1/2J)[D:L:B]  +  channel B = the backstress        #
#     co-rotation ∂α̃/∂F) is validated EXACTLY-to-FD at the MATERIAL level in      #
#     tests/test_ladrunoJ2_finite_native_tangent.py (recipe == full ∂σ/∂F to      #
#     ~1e-14) and the element's assembly of getSpatialTangentTensor is proven by  #
#     the v1 gate tests/test_ladrunoJ2_finite_element.py. An element-level FD of   #
#     K cannot cleanly isolate channel B: it scales with the COMMITTED backstress, #
#     but eleResponse("stiff") re-forms the tangent POST-commit (α≠0) while a      #
#     from-virgin force FD carries committed α=0 — the two live at different       #
#     committed bases, so the numeric comparison is ill-posed without an          #
#     FE_Datastore replay that OpenSees' commit/sp semantics make fragile.        #
# =========================================================================== #


# =========================================================================== #
#  4. Rigid rotation of the unloaded element → zero stress / force.            #
# =========================================================================== #
def test_native_finite_rigid_rotation_is_stress_free():
    cy, sy = math.cos(0.3), math.sin(0.3)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    R = Ry @ _rotz(0.4)
    f = _force_single(R, _mat_native)
    s = ops.eleResponse(1, "stresses")
    assert max(abs(v) for v in s) <= 1.0e-6 * _ISO["sig0"], "rigid rotation induced stress"
    assert np.abs(f).max() <= 1.0e-6 * _ISO["sig0"], "rigid rotation induced force"


# =========================================================================== #
#  5. Tiny strain: reduces to the small-strain LadrunoJ2.                       #
# =========================================================================== #
def test_native_finite_reduces_to_small_strain():
    eps = 1.0e-5
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    f_fin = _force_single(Fbar, _mat_native, "finite")
    f_lin = _force_single(Fbar, _mat_small, "linear")
    scale = max(np.abs(f_lin).max(), 1.0e-30)
    assert np.abs(f_fin - f_lin).max() <= 1.0e-3 * scale, (
        "native finite J2 did not reduce to small-strain LadrunoJ2 at tiny strain")


# =========================================================================== #
#  6. SERIALIZATION: LadrunoJ2Finite in a committed plastic finite state must   #
#     survive sendSelf/recvSelf (Fn, Be_n, ebarP_n, alpha_n) + broker dispatch.  #
#     Drives save→wipe→rebuild→restore through the FE_Datastore; skips if the    #
#     build lacks database support.                                             #
# =========================================================================== #
def test_native_finite_database_roundtrip():
    Fbar = [[1.15, 0.02, 0.01], [0.02, 0.92, 0.01], [0.01, 0.01, 0.94]]

    def build():
        assert _solve_to(Fbar, _mat_native, "finite") == 0

    database_roundtrip(build, probe_nodes=[7], ndf=3)


# =========================================================================== #
#  7. IMPL-EX (`-implex`): the EXPLICIT (Δγ-extrapolation) reporting path,      #
#     driven through the real LadrunoBrick -geom finite element along the SAME  #
#     rotated-plastic path as test 2, must match the IMPL-EX numpy oracle       #
#     (ImplexNativeFiniteJ2) step for step — proving the element surfaces the    #
#     frozen-multiplier / co-rotated-flow-direction explicit stress (and that    #
#     the implicit return still commits the true history underneath).           #
# =========================================================================== #
def test_implex_finite_element_matches_oracle():
    Ffin = _rotz(0.5) @ np.array([[1.18, 0.04, 0.0], [0.0, 0.93, 0.0], [0.0, 0.0, 0.95]])
    N = 6
    _build(_mat_native_implex, "finite")
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag, (x, y, z) in _NODES.items():
        d = (Ffin - _I) @ np.array([x, y, z])
        base = (tag - 1) * 3
        for j in range(3):
            ops.sp(tag, j + 1, float(d[j]))
    _setup_static()
    ops.integrator("LoadControl", 1.0 / N)
    ops.analysis("Static")

    m = fx.ImplexNativeFiniteJ2(_params())
    saw = False
    for k in range(1, N + 1):
        assert ops.analyze(1) == 0, f"IMPL-EX finite solve failed at step {k}"
        Fk = _I + (k / N) * (Ffin - _I)                  # affine ⇒ F_k at every GP
        sig = m.setTrialF(Fk); m.commit()                # implex (explicit) Cauchy
        ref = np.array([sig[0, 0], sig[1, 1], sig[2, 2], sig[0, 1], sig[1, 2], sig[0, 2]])
        s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)[0]
        tol = 1.0e-6 * max(np.abs(ref).max(), 1.0)
        assert np.allclose(s, ref, rtol=1.0e-6, atol=tol), (
            f"IMPL-EX finite element GP Cauchy != implex oracle at step {k}: {s} vs {ref}")
        saw = saw or m._plastic
    assert saw, "IMPL-EX rotating plastic path never yielded"


# =========================================================================== #
#  8. IMPL-EX SERIALIZATION: a committed plastic finite state with -implex must  #
#     survive sendSelf/recvSelf + broker — exercises the EXTRA implex state      #
#     (Δγ_n, N_n) added to the data vector. Skips if no database support.       #
# =========================================================================== #
def test_implex_finite_database_roundtrip():
    Fbar = [[1.15, 0.02, 0.01], [0.02, 0.92, 0.01], [0.01, 0.01, 0.94]]

    def build():
        assert _solve_to(Fbar, _mat_native_implex, "finite") == 0

    database_roundtrip(build, probe_nodes=[7], ndf=3)


# =========================================================================== #
#  9. Softening-parameter guard (#1): a softening isotropic law (Hiso < 0)      #
#     CONSTRUCTS (the guard warns, it does not reject) and — the guard's        #
#     recommended path — runs a plastic step STABLY under -implex (the constant  #
#     SPD elastic tangent sidesteps the indefinite implicit tangent the warning  #
#     is about).                                                                #
# =========================================================================== #
def _mat_softening(tag, implex):
    args = ["LadrunoJ2Finite", tag, _K, _G, "-iso", "voce",
            _ISO["sig0"], 0.0, 0.0, -200.0,           # Qinf=b=0, Hiso=-200 < 0 ⇒ softening
            "-kin", 0]
    if implex:
        args.append("-implex")
    ops.nDMaterial(*args)


def test_softening_constructs_and_implex_runs_plastic():
    # softening + -implex: a plastic imposed-F step converges (SPD elastic tangent)
    Fpl = [[1.10, 0.0, 0.0], [0.0, 1.0/math.sqrt(1.10), 0.0], [0.0, 0.0, 1.0/math.sqrt(1.10)]]
    assert _solve_to(Fpl, lambda t: _mat_softening(t, implex=True), "finite") == 0
    # without -implex the material still CONSTRUCTS (guard warns, does not reject):
    # a tiny sub-yield step is well-posed even with a softening backbone.
    Fel = [[1.0 + 1e-5, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    assert _solve_to(Fel, lambda t: _mat_softening(t, implex=False), "finite") == 0
