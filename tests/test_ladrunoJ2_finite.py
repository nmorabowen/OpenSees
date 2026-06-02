"""Finite-strain combined-hardening J2 — the LogStrain-over-LadrunoJ2 composition,
proven at the math level WITHOUT building OpenSees.

This is the numpy counterpart of the deliverable "LadrunoBrick -geom finite gets
large-strain combined-hardening J2". The element drives a FiniteStrainNDMaterial
(LogStrainNDMaterial) by setTrialF(F); LogStrain wraps an UNCHANGED small-strain
3D NDMaterial via the plastic-inner state protocol. Because LadrunoJ2 presents
the exact J2ThreeDimensional interface (engineering-shear strain in / true stress
out, linear elastic, stateful εᵖ subtraction, constant getInitialTangent), it is
a valid inner — so the finite-strain path is literally `LogStrain` over
`LadrunoJ2`, with the verified combined-hardening return map (LadrunoJ2Kernel.h)
running inside.

Here we drive the protocol adaptor over the numpy oracle StatefulLadrunoJ2
(tests/ladrunoj2_reference.py) and verify it reproduces the DIRECT Box-14.4
elastic-trial chain step for step along a finite path with stretch + shear +
finite rotation + unloading, with three-term Chaboche kinematic hardening. The
C++ kernel itself is checked separately (test_ladrunoJ2_kernel_cpp.py); the
log-strain kinematics/tangent are checked in the LogStrain suite; this file
proves the COMPOSITION (the εᵖ-double-count-free protocol) for combined hardening.
"""
import numpy as np
import pytest

import ladrunoj2_reference as lr
import logstrain_reference as ls

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

I3 = np.eye(3)
# A genuine metal with three-term Chaboche AF kinematic + Voce+linear isotropic.
PARAMS = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                   C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])
# Pure isotropic (Voce+linear) — the part the log-strain wrapper lifts EXACTLY.
PARAMS_ISO = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0)


# --------------------------------------------------------------------------- #
#  Direct Box-14.4 chain: elastic-trial Hencky strain in, no εᵖ carry.        #
#  (Uses the verified return map with epsP_n = 0, i.e. the trial IS elastic.) #
# --------------------------------------------------------------------------- #
def _direct_step(F, Be_n, Fn, ebarP_n, alpha_n, params=PARAMS):
    Fd = F @ np.linalg.inv(Fn)
    Be_tr = Fd @ Be_n @ Fd.T
    eps_e_tr = ls.hencky_strain(Be_tr)
    res = lr.return_map(params, eps_e_tr, np.zeros((3, 3)), ebarP_n, alpha_n)
    tau = res["stress"]                         # elastic-trial in ⇒ τ = Dᵉ:εᵉ
    eps_e = eps_e_tr - res["epsP"]              # εᵉ_{n+1} = εᵉᵗʳ − Δεᵖ
    J = float(np.linalg.det(F))
    Be_new = ls.iso_function(2.0 * eps_e, np.exp)
    return dict(sigma=tau / J, Be=Be_new, alpha=res["alpha"], ebarP=res["ebarP"],
                eps_e=eps_e, J=J, plastic=res["plastic"])


# --------------------------------------------------------------------------- #
#  The plastic-inner protocol adaptor (mirrors LogStrainNDMaterial.cpp).       #
# --------------------------------------------------------------------------- #
class LogStrainAdaptor:
    def __init__(self, inner: lr.StatefulLadrunoJ2):
        self.inner = inner
        self.Be_n = I3.copy()
        self.eps_feed_n = np.zeros((3, 3))
        self.Fn = I3.copy()

    def setTrialF(self, F):
        Fd = F @ np.linalg.inv(self.Fn)
        Be_tr = Fd @ self.Be_n @ Fd.T
        eps_e_trial = ls.hencky_strain(Be_tr)
        eps_e_n = ls.hencky_strain(self.Be_n)
        eps_feed = self.eps_feed_n + (eps_e_trial - eps_e_n)   # neutralize inner εᵖ
        self.inner.setTrialStrain(eps_feed)
        tau = self.inner.getStress()
        J = float(np.linalg.det(F))
        eps_e_np1 = self.inner.getInitialTangent_compliance(tau)   # Cᵉ:τ recovery
        self._stage = (ls.iso_function(2.0 * eps_e_np1, np.exp), eps_feed.copy(), F.copy())
        return tau / J, eps_e_np1, J

    def commitState(self):
        self.Be_n, self.eps_feed_n, self.Fn = self._stage
        self.inner.commitState()


def _path():
    """Finite deformation path: growing axial stretch + shear + rotation, then
    partial unloading."""
    Fs = []
    for k in range(1, 13):
        s = 1.0 + 0.03 * k
        g = 0.02 * k
        th = 0.05 * k
        U = np.array([[s, g, 0.0],
                      [0.0, 1.0 / np.sqrt(s), 0.0],
                      [0.0, 0.0, 1.0 / np.sqrt(s)]])
        c, sn = np.cos(th), np.sin(th)
        R = np.array([[c, -sn, 0.0], [sn, c, 0.0], [0.0, 0.0, 1.0]])
        Fs.append(R @ U)
    for k in range(5):
        Fs.append(Fs[len(Fs) - 2 - k])     # unload back down the path
    return Fs


def test_finite_combined_hardening_matches_direct_chain():
    """Adaptor-over-StatefulLadrunoJ2 reproduces the direct elastic-trial chain
    step for step — proving the protocol neutralizes the εᵖ subtraction for
    Chaboche combined hardening (backstress lives inside the inner)."""
    adaptor = LogStrainAdaptor(lr.StatefulLadrunoJ2(PARAMS))
    Be_n, Fn, ebarP_n, alpha_n = I3.copy(), I3.copy(), 0.0, [np.zeros((3, 3))] * PARAMS.nBack

    saw_plastic = False
    for F in _path():
        sig_proto, _, _ = adaptor.setTrialF(F)
        direct = _direct_step(F, Be_n, Fn, ebarP_n, alpha_n)
        saw_plastic = saw_plastic or direct["plastic"]
        assert np.allclose(sig_proto, direct["sigma"], rtol=1e-9, atol=1e-9), (
            f"protocol σ disagrees with direct chain:\n{sig_proto}\n{direct['sigma']}")
        adaptor.commitState()
        Be_n, Fn, ebarP_n, alpha_n = direct["Be"], F, direct["ebarP"], direct["alpha"]
    assert saw_plastic, "path never yielded — strengthen the loading"


def test_finite_plastic_incompressibility():
    """Deviatoric (J2) plastic flow ⇒ det Fᵖ = 1 ⇒ det bᵉ = J², tr εᵉ = ln J at
    every step, even with combined hardening and through the protocol."""
    adaptor = LogStrainAdaptor(lr.StatefulLadrunoJ2(PARAMS))
    for F in _path():
        _, eps_e, J = adaptor.setTrialF(F)
        assert np.trace(eps_e) == pytest.approx(np.log(J), abs=1e-10)
        Be_new = adaptor._stage[0]
        assert np.linalg.det(Be_new) == pytest.approx(J * J, rel=1e-9)
        adaptor.commitState()


def _load_then_superpose_rotation(params):
    """Load into the plastic range along a fixed stretch (committing each step),
    then superpose a finite rigid rotation Q on the committed state. Returns
    (sigma_before, Q, sigma_after) where frame indifference requires
    sigma_after == Q sigma_before Qᵀ."""
    adaptor = LogStrainAdaptor(lr.StatefulLadrunoJ2(params))
    F_def = np.array([[1.25, 0.10, 0.0], [0.0, 1.0 / np.sqrt(1.25), 0.0],
                      [0.0, 0.0, 1.0 / np.sqrt(1.25)]])
    sig_ref = None
    for a in np.linspace(0.2, 1.0, 8):
        Fa = I3 + a * (F_def - I3)
        sig_ref, _, _ = adaptor.setTrialF(Fa)
        adaptor.commitState()
    th = 0.6
    c, s = np.cos(th), np.sin(th)
    Q = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    sig_rot, _, _ = adaptor.setTrialF(Q @ F_def)
    return sig_ref, Q, sig_rot


def test_finite_isotropic_objectivity_is_exact():
    """The log-strain lift is EXACTLY frame-indifferent for the isotropic spine:
    the elastic state co-rotates via Bᵉᵗʳ = F_Δ Bᵉ_n F_Δᵀ and isotropic yield
    depends only on ‖s‖ and ε̄ᵖ — so superposing a finite rotation rotates the
    Cauchy stress rigidly, σ_rot = Q σ_ref Qᵀ."""
    sig_ref, Q, sig_rot = _load_then_superpose_rotation(PARAMS_ISO)
    assert np.allclose(sig_rot, Q @ sig_ref @ Q.T, rtol=1e-9, atol=1e-9), (
        "isotropic finite-strain J2 is not objective (it must be)")


@pytest.mark.xfail(strict=True, reason=(
    "Kinematic backstress does not co-rotate in the v1 log-strain wrapper: the "
    "inner material stores α in a fixed frame, so a large superposed rotation is "
    "NOT frame-indifferent for combined hardening. This is the de Souza Neto §14.11 "
    "boundary — scoped to v2 (a finite-strain-NATIVE J2 reusing LadrunoJ2Kernel "
    "with a co-rotated backstress; the kernel extraction is what enables it). "
    "Isotropic hardening and moderate rotations are unaffected. When v2 lands and "
    "co-rotates α, this test flips to PASS and the strict-xfail flags it."))
def test_finite_combined_hardening_large_rotation_objectivity_is_v2():
    sig_ref, Q, sig_rot = _load_then_superpose_rotation(PARAMS)
    assert np.allclose(sig_rot, Q @ sig_ref @ Q.T, rtol=1e-7, atol=1e-7)
