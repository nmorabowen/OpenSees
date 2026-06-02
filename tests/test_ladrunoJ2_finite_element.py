"""LadrunoBrick ``-geom finite`` driving finite-strain combined-hardening J2 —
the end-to-end element acceptance for LogStrain-over-LadrunoJ2.

    nDMaterial LadrunoJ2 1 K G -iso voce ... -kin N ...   # small-strain inner
    nDMaterial LogStrain 2 1                              # Hencky finite-strain lift
    element LadrunoBrick 1 ...nodes... 2 -geom finite     # F → setTrialF

LadrunoJ2 presents the exact J2ThreeDimensional interface the LogStrain plastic-
inner protocol requires (engineering-shear strain in / true stress out, linear
elastic, stateful εᵖ subtraction, constant getInitialTangent), so it wraps with
NO change to LogStrain — the verified combined-hardening return map
(LadrunoJ2Kernel.h) runs inside.

Gates (requires a local/CI build — ``from _testbed import ops``):
  * consistent tangent K == FD of the resisting force, in the plastic range
    (combined hardening) — the headline arbiter;
  * a uniaxial homogeneous finite stretch into plasticity → constant Cauchy
    stress matching the numpy direct Box-14.4 chain (tests/ladrunoj2_reference.py);
  * pure rigid rotation of the (unloaded) element → zero stress / force;
  * tiny strain → finite reduces to the small-strain LadrunoJ2.

Objectivity of the *kinematic* part under LARGE rotation is the de Souza Neto
§14.11 v2 boundary (backstress does not co-rotate in the v1 wrapper) and is pinned
in tests/test_ladrunoJ2_finite.py; tests here that combine plasticity with finite
rotation use ISOTROPIC hardening (exact), while the no-rotation tests exercise the
full Chaboche combined hardening.
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip
import ladrunoj2_reference as lr
import logstrain_reference as ls

pytestmark = [pytest.mark.zone_a]

_NODES = {
    1: (0.00, 0.00, 0.00), 2: (1.00, 0.10, 0.05), 3: (1.10, 1.00, 0.00),
    4: (0.05, 0.95, 0.10), 5: (0.00, 0.05, 1.00), 6: (1.00, 0.00, 1.05),
    7: (1.05, 1.00, 1.10), 8: (0.00, 1.00, 0.95),
}
_CONN = [1, 2, 3, 4, 5, 6, 7, 8]
_K, _G = 1500.0, 700.0
# full Chaboche combined hardening (used where there is no finite rotation)
_ISO = dict(sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0)
_KIN = dict(C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])


def _mat_combined(tag):
    kin = []
    for c, g in zip(_KIN["C"], _KIN["gam"]):
        kin += [c, g]
    ops.nDMaterial("LadrunoJ2", tag, _K, _G, "-iso", "voce",
                   _ISO["sig0"], _ISO["Qinf"], _ISO["bIso"], _ISO["Hiso"],
                   "-kin", len(_KIN["C"]), *kin)


def _mat_iso(tag):
    ops.nDMaterial("LadrunoJ2", tag, _K, _G, "-iso", "voce",
                   _ISO["sig0"], _ISO["Qinf"], _ISO["bIso"], _ISO["Hiso"], "-kin", 0)


def _build(geom, mat_fn, combined_inner=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for tag, (x, y, z) in _NODES.items():
        ops.node(tag, x, y, z)
    mat_fn(1)
    if geom == "finite":
        ops.nDMaterial("LogStrain", 2, 1)
        mtag = 2
    else:
        mtag = 1
    ops.element("LadrunoBrick", 1, *_CONN, mtag, "-formulation", "std", "-geom", geom)
    return mtag


def _impose_and_solve(u, geom, mat_fn):
    _build(geom, mat_fn)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for tag in _NODES:
        base = (tag - 1) * 3
        for d in range(3):
            ops.sp(tag, d + 1, float(u[base + d]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    return ops.analyze(1)


def _resisting_force(u, mat_fn, geom="finite"):
    assert _impose_and_solve(u, geom, mat_fn) == 0, "imposed-displacement solve failed"
    return np.array(ops.eleForce(1), dtype=float)


def _affine_disp(Fbar, wiggle=None):
    u = np.zeros(24)
    I = np.eye(3)
    for tag, (x, y, z) in _NODES.items():
        X = np.array([x, y, z])
        d = (np.asarray(Fbar) - I) @ X
        if wiggle is not None:
            d = d + np.array(wiggle[tag - 1])
        u[(tag - 1) * 3:(tag - 1) * 3 + 3] = d
    return u


_WIGGLE = [
    (0.000, 0.000, 0.000), (0.010, -0.007, 0.005), (-0.006, 0.009, -0.008),
    (0.007, 0.004, 0.010), (-0.009, 0.006, -0.005), (0.005, -0.010, 0.007),
    (-0.008, 0.006, 0.009), (0.010, -0.005, -0.006),
]


# --------------------------------------------------------------------------- #
#  Headline: consistent tangent == FD of resisting force, in the plastic range #
#  (full Chaboche combined hardening; a stretch with NO large rigid rotation). #
# --------------------------------------------------------------------------- #
def test_finite_j2_consistent_tangent_matches_fd():
    # A symmetric-ish stretch well into yield (deviatoric stretch ~20%): the
    # polar rotation is ~I, so combined hardening is exercised without invoking
    # the §14.11 large-rotation backstress-corotation boundary.
    Fbar = [[1.20, 0.02, 0.01],
            [0.02, 0.90, 0.015],
            [0.01, 0.015, 0.93]]
    u = _affine_disp(Fbar, wiggle=_WIGGLE)

    assert _impose_and_solve(u, "finite", _mat_combined) == 0
    # guard against a vacuous test: the state MUST be plastic, else K==Kfd==elastic
    # tangent and the combined-hardening elastoplastic tangent is never exercised.
    ebar = list(ops.eleResponse(1, "material", 1, "equivalentPlasticStrain"))
    assert ebar and ebar[0] > 0.0, "tangent test is vacuous — element did not yield"
    K = np.array(ops.eleResponse(1, "stiff"), dtype=float).reshape(24, 24)
    f0 = ops.eleForce(1)
    assert np.isfinite(K).all()

    h = 1.0e-6
    Kfd = np.zeros((24, 24))
    for d in range(24):
        up = u.copy(); up[d] += h
        um = u.copy(); um[d] -= h
        Kfd[:, d] = (_resisting_force(up, _mat_combined) -
                     _resisting_force(um, _mat_combined)) / (2.0 * h)

    scale = np.abs(K).max()
    err = np.abs(K - Kfd).max()
    assert err <= 2.0e-4 * scale, (
        f"finite J2 tangent != FD of resisting force: max abs err {err:.3e} "
        f"vs scale {scale:.3e}")


# --------------------------------------------------------------------------- #
#  Uniaxial homogeneous finite stretch into plasticity → constant Cauchy       #
#  stress matching the numpy direct Box-14.4 chain (R = I ⇒ combined is exact). #
# --------------------------------------------------------------------------- #
def test_finite_j2_uniaxial_patch_matches_oracle():
    lam = 1.18                      # ~16.6% log axial strain — solidly plastic
    lat = 1.0 / math.sqrt(lam)      # isochoric-ish lateral (just a prescribed F)
    Fm = np.diag([lam, lat, lat])
    u = _affine_disp(Fm.tolist())
    assert _impose_and_solve(u, "finite", _mat_combined) == 0

    s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)
    s0 = s[0]
    assert np.abs(s - s0).max() <= 1.0e-6 * max(np.abs(s0).max(), 1.0), (
        "uniaxial finite patch: Cauchy stress not constant across Gauss points")

    # numpy direct chain (single monotonic step from the reference config)
    P = lr.Params(_K, _G, **_ISO, **_KIN)
    Be_tr = Fm @ Fm.T
    eps_e_tr = ls.hencky_strain(Be_tr)
    res = lr.return_map(P, eps_e_tr, np.zeros((3, 3)), 0.0,
                        [np.zeros((3, 3))] * P.nBack)
    assert res["plastic"], "oracle did not yield — increase lam"
    J = float(np.linalg.det(Fm))
    sig = res["stress"] / J
    sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2],
                      sig[0, 1], sig[1, 2], sig[2, 0]])
    tol = 1.0e-6 * max(np.abs(sig_v).max(), 1.0)
    assert np.allclose(s0, sig_v, rtol=1.0e-6, atol=tol), (
        f"finite J2 patch: GP Cauchy {s0} != oracle {sig_v}")


# --------------------------------------------------------------------------- #
#  MULTI-STEP committed finite-plastic path: the only test that exercises the   #
#  COMPILED commitState -> next setTrialF round-trip (Be_n / epsFeed_n carry,    #
#  be_from_hencky_voigt recovery) across several committed plastic steps. A      #
#  sign/index slip that is identity on step 1 (Be_n=I) but bites on step 2+      #
#  would pass every single-step gate but fail here.                             #
# --------------------------------------------------------------------------- #
def test_finite_j2_multistep_committed_matches_oracle():
    # Symmetric F (polar R = I exactly) so combined hardening stays in the
    # exactly-objective regime; imposed in N increments F_k = I + (k/N)(Ffinal-I)
    # via LoadControl, committing each step. Compare every committed GP Cauchy
    # stress to the numpy direct Box-14.4 chain stepped over the SAME F_k.
    Ffinal = np.array([[1.18, 0.02, 0.01],
                       [0.02, 0.93, 0.015],
                       [0.01, 0.015, 0.95]])
    N = 6
    _build("finite", _mat_combined)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    I = np.eye(3)
    for tag, (x, y, z) in _NODES.items():
        d = (Ffinal - I) @ np.array([x, y, z])
        base = (tag - 1) * 3
        for j in range(3):
            ops.sp(tag, j + 1, float(d[j]))
    ops.constraints("Lagrange")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-11, 50, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / N)
    ops.analysis("Static")

    P = lr.Params(_K, _G, **_ISO, **_KIN)
    Be_n, Fn, ebarP, alpha = I.copy(), I.copy(), 0.0, [np.zeros((3, 3))] * P.nBack
    saw_plastic = False
    for k in range(1, N + 1):
        assert ops.analyze(1) == 0, f"multistep finite solve failed at step {k}"
        Fk = I + (k / N) * (Ffinal - I)             # affine => F_k at every GP
        Fd = Fk @ np.linalg.inv(Fn)
        eps_e_tr = ls.hencky_strain(Fd @ Be_n @ Fd.T)
        res = lr.return_map(P, eps_e_tr, np.zeros((3, 3)), ebarP, alpha)
        saw_plastic = saw_plastic or res["plastic"]
        J = float(np.linalg.det(Fk))
        sig = res["stress"] / J
        sig_v = np.array([sig[0, 0], sig[1, 1], sig[2, 2],
                          sig[0, 1], sig[1, 2], sig[2, 0]])
        # advance the oracle's committed state exactly as the C++ adaptor does
        Be_n = ls.iso_function(2.0 * (eps_e_tr - res["epsP"]), np.exp)
        Fn, ebarP, alpha = Fk, res["ebarP"], res["alpha"]

        s = np.array(ops.eleResponse(1, "stresses"), dtype=float).reshape(8, 6)[0]
        tol = 1.0e-6 * max(np.abs(sig_v).max(), 1.0)
        assert np.allclose(s, sig_v, rtol=1.0e-6, atol=tol), (
            f"multistep finite step {k}: GP Cauchy {s} != oracle {sig_v}")
    assert saw_plastic, "multistep finite path never yielded"


# --------------------------------------------------------------------------- #
#  Pure rigid rotation of the unloaded element → zero stress / force.          #
# --------------------------------------------------------------------------- #
def test_finite_j2_rigid_rotation_is_stress_free():
    th = 0.4
    cz, sz = math.cos(th), math.sin(th)
    Rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    cy, sy = math.cos(0.3), math.sin(0.3)
    Ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    R = Ry @ Rz
    u = _affine_disp(R)
    assert _impose_and_solve(u, "finite", _mat_combined) == 0
    s = ops.eleResponse(1, "stresses")
    smax = max(abs(v) for v in s)
    assert smax <= 1.0e-6 * _ISO["sig0"], f"rigid rotation induced stress {smax:.3e}"
    f = np.array(ops.eleForce(1), dtype=float)
    assert np.abs(f).max() <= 1.0e-6 * _ISO["sig0"], "rigid rotation induced force"


# --------------------------------------------------------------------------- #
#  Tiny strain: finite reduces to the small-strain LadrunoJ2 (same std kernel).#
# --------------------------------------------------------------------------- #
def test_finite_j2_reduces_to_small_strain():
    eps = 1.0e-5
    Fbar = [[1.0 + eps, 0.5 * eps, 0.0],
            [0.5 * eps, 1.0 - 0.3 * eps, 0.0],
            [0.0, 0.0, 1.0 + 0.2 * eps]]
    u = _affine_disp(Fbar)
    f_fin = _resisting_force(u, _mat_combined, "finite")
    f_lin = _resisting_force(u, _mat_combined, "linear")   # small-strain LadrunoJ2
    scale = max(np.abs(f_lin).max(), 1.0e-30)
    assert np.abs(f_fin - f_lin).max() <= 1.0e-3 * scale, (
        "finite J2 did not reduce to small-strain LadrunoJ2 at tiny strain")


# --------------------------------------------------------------------------- #
#  Serialization: LogStrain(LadrunoJ2) in a committed plastic finite state must #
#  survive sendSelf/recvSelf (the new protocol state Fn/Be_n/epsFeed_n is       #
#  serialized in LogStrainNDMaterial). This drives save->wipe->rebuild->restore #
#  through the FE_Datastore, exercising the element + both materials' send/recv #
#  and the broker reconstruction — guarding against a dropped or mis-sized      #
#  field (the classic bolt-on-state regression). Skips if the build lacks DB.   #
# --------------------------------------------------------------------------- #
def test_finite_j2_database_roundtrip():
    u = _affine_disp([[1.15, 0.02, 0.01],
                      [0.02, 0.92, 0.01],
                      [0.01, 0.01, 0.94]])

    def build():
        assert _impose_and_solve(u, "finite", _mat_combined) == 0

    database_roundtrip(build, probe_nodes=[7], ndf=3)
