"""Prototype + verification of the PLASTIC-INNER state protocol for the
log-strain finite-strain adaptor — the piece that lets it wrap a *stateful*
small-strain plastic material (J2ThreeDimensional, DruckerPrager3D, …) without
double-counting the plastic strain.

The problem: a stock OpenSees small-strain material treats setTrialStrain() input
as TOTAL strain and subtracts its OWN committed plastic strain εᵖ. But the
multiplicative adaptor already carries the elastic history in bᵉ_n and wants to
feed the *trial elastic* Hencky strain εᵉᵗʳ — so εᵖ gets double-counted.

The protocol (rigorous multiplicative bᵉ carry, reusing an UNMODIFIED stateful
material by NEUTRALIZING its εᵖ subtraction):
  adaptor state: bᵉ_n  and  ε_feed_n (the total strain last fed to the inner).
  invariant (preserved by construction): the inner's committed elastic strain
      equals εᵉ_n = ½ln bᵉ_n, hence its committed εᵖ_n = ε_feed_n − εᵉ_n.
  per step:
      εᵉᵗʳ = ½ln(F_Δ bᵉ_n F_Δᵀ);  εᵉ_n = ½ln bᵉ_n
      ε_feed = ε_feed_n + (εᵉᵗʳ − εᵉ_n)
      → inner sees ε_feed − εᵖ_n = εᵉᵗʳ  (exactly the log elastic trial) ✓
      τ = inner stress;  σ = τ/J
      εᵉ_{n+1} = Cᵉ : τ  (recover from stress via the elastic compliance)
      bᵉ_{n+1} = exp[2 εᵉ_{n+1}]
  This reproduces de Souza Neto Box 14.3/14.4 exactly; hardening α persists
  naturally inside the inner.

This file verifies the protocol numerically against the DIRECT (stateless,
elastic-strain-in) Box 14.4 chain already in logstrain_reference.matisu_step.
Pure numpy, no openseespy. Mirrors what the C++ LogStrainNDMaterial will do.
"""
import numpy as np
import pytest

import logstrain_reference as lr

I3 = np.eye(3)
J2 = dict(K=1.5e3, G=7.0e2, sig_y0=10.0, H=50.0)


# --------------------------------------------------------------------------- #
#  A STATEFUL small-strain J2, mimicking an OpenSees NDMaterial: it stores εᵖ  #
#  and α, and setTrialStrain() input is the TOTAL strain (εᵉ = ε − εᵖ_n).      #
# --------------------------------------------------------------------------- #
class StatefulSmallStrainJ2:
    def __init__(self, K, G, sig_y0, H):
        self.p = dict(K=K, G=G, sig_y0=sig_y0, H=H)
        self.eps_p = np.zeros((3, 3)); self.alpha = 0.0           # committed
        self._eps_p_tr = self.eps_p.copy(); self._alpha_tr = 0.0  # trial
        self._tau = np.zeros((3, 3))

    def setTrialStrain(self, eps_total):                 # eps_total: 3×3 tensor
        eps_e_tr = eps_total - self.eps_p                # subtract COMMITTED εᵖ
        tau, eps_e, alpha, D, _ = lr.j2_return_map(eps_e_tr, self.alpha, **self.p)
        self._tau = tau
        self._eps_p_tr = eps_total - eps_e               # trial plastic strain
        self._alpha_tr = alpha

    def getStress(self):          return self._tau
    def getInitialTangent(self):  return lr.elastic_tangent_D(self.p["K"], self.p["G"])
    def commitState(self):
        self.eps_p = self._eps_p_tr.copy(); self.alpha = self._alpha_tr


# --------------------------------------------------------------------------- #
#  The PLASTIC-INNER protocol adaptor (what the C++ will implement)           #
# --------------------------------------------------------------------------- #
class LogStrainAdaptor:
    def __init__(self, inner, Ke, Ge):
        self.inner = inner
        self.Ke, self.Ge = Ke, Ge          # inner elastic constants → Cᵉ recovery
        self.Be_n = I3.copy()
        self.eps_feed_n = np.zeros((3, 3))
        self.Fn = I3.copy()

    def _recover_eps_e(self, tau):         # εᵉ = Cᵉ : τ  (isotropic elastic compliance)
        tr = np.trace(tau)
        return (tau - (tr/3.0)*I3)/(2.0*self.Ge) + (tr/(9.0*self.Ke))*I3

    def setTrialF(self, F):
        Fd = F @ np.linalg.inv(self.Fn)
        Be_tr = Fd @ self.Be_n @ Fd.T
        eps_e_trial = lr.hencky_strain(Be_tr)
        eps_e_n = lr.hencky_strain(self.Be_n)
        eps_feed = self.eps_feed_n + (eps_e_trial - eps_e_n)   # neutralize inner εᵖ
        self.inner.setTrialStrain(eps_feed)
        tau = self.inner.getStress()
        J = float(np.linalg.det(F))
        sigma = tau / J
        eps_e_np1 = self._recover_eps_e(tau)
        self._stage = (lr.iso_function(2.0*eps_e_np1, np.exp), eps_feed.copy(), F.copy())
        return sigma, tau, J, eps_e_np1

    def commitState(self):
        self.Be_n, self.eps_feed_n, self.Fn = self._stage
        self.inner.commitState()                       # ← the inner owns εᵖ, α


# --------------------------------------------------------------------------- #
#  A deformation path with plastic loading, shear, rotation and unloading.    #
# --------------------------------------------------------------------------- #
def _path():
    Fs = []
    th = 0.0
    for k in range(1, 13):
        s = 1.0 + 0.02*k                       # growing axial stretch (into yield)
        g = 0.015*k                            # growing shear
        th = 0.04*k                            # growing rotation about z
        U = np.array([[s, g, 0.0],
                      [0.0, 1.0/np.sqrt(s), 0.0],
                      [0.0, 0.0, 1.0/np.sqrt(s)]])
        c, sn = np.cos(th), np.sin(th)
        R = np.array([[c, -sn, 0.0], [sn, c, 0.0], [0.0, 0.0, 1.0]])
        Fs.append(R @ U)
    for k in range(5):                         # partial unloading
        Fs.append(Fs[-1] * 0.0 + Fs[len(Fs)-2-k])
    return Fs


@pytest.mark.zone_a
@pytest.mark.t1
def test_plastic_protocol_matches_direct_chain():
    """Adaptor-over-stateful-J2 must reproduce the direct stateless Box-14.4
    chain (matisu_step 'j2') step-for-step along a finite-strain plastic path."""
    inner = StatefulSmallStrainJ2(**J2)
    adaptor = LogStrainAdaptor(inner, J2["K"], J2["G"])

    Be_n, alpha_n, Fn = I3.copy(), 0.0, I3.copy()      # direct-chain state
    saw_plastic = False
    for F in _path():
        s_proto, _, J, _ = adaptor.setTrialF(F)

        direct = lr.matisu_step(F, Be_n, Fn, "j2", alpha_n=alpha_n, **J2)
        saw_plastic = saw_plastic or direct["plastic"]

        assert np.allclose(s_proto, direct["sigma"], rtol=1e-9, atol=1e-9), (
            f"protocol σ disagrees with direct chain:\n{s_proto}\n{direct['sigma']}")

        adaptor.commitState()
        Be_n, alpha_n, Fn = direct["Be"], direct["alpha"], F
    assert saw_plastic, "path never yielded — strengthen the loading"


@pytest.mark.zone_a
@pytest.mark.t1
def test_plastic_protocol_exact_incompressibility():
    """Deviatoric (J2) plastic flow ⇒ det Fᵖ = 1 ⇒ det bᵉ = J² and tr εᵉ = ln J
    at every step (Remark 14.4), even through the stateful-inner protocol."""
    inner = StatefulSmallStrainJ2(**J2)
    adaptor = LogStrainAdaptor(inner, J2["K"], J2["G"])
    for F in _path():
        _, _, J, eps_e = adaptor.setTrialF(F)
        assert np.trace(eps_e) == pytest.approx(np.log(J), abs=1e-10)
        Be_new = adaptor._stage[0]
        assert np.linalg.det(Be_new) == pytest.approx(J*J, rel=1e-9)
        adaptor.commitState()
