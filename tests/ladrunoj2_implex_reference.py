"""IMPL-EX (Oliver–Huespe–Cante 2008) oracle for combined-hardening J2 — the numpy
reference the C++ `nDMaterial LadrunoJ2[Finite] -implex` Δγ-extrapolation path must
match. Written BEFORE the C++ (the fork's oracle-first pattern).

WHAT IMPL-EX IS (here, for plastic flow — NOT the damage variable). The classic
implicit/explicit split decouples the GLOBAL Newton tangent from the LOCAL return map:

  * IMPLICIT pass (runs every step, in the background): the UNCHANGED return map
    (`ladrunoj2_reference.return_map`) gives the true Δγ, εᵖ, ε̄ᵖ, α and the flow
    direction N = (εᵖ_{n+1} − εᵖ_n)/Δγ. THIS is what gets committed — the committed
    plastic history is therefore byte-identical to a fully-implicit run on the same
    prescribed strain path.

  * EXPLICIT (implex) stress reported to the solver: freeze an EXTRAPOLATED multiplier
    and the COMMITTED flow direction, so the plastic strain is a pure history quantity
        Δγ̃ = Δγ_n · (Δt_{n+1}/Δt_n)        (uniform step ⇒ Δγ̃ = Δγ_n)
        εᵖ̃ = εᵖ_n + Δγ̃ · N_n               (N_n committed, frozen)
        σ̃  = C : (ε − εᵖ̃)
    Because εᵖ̃ carries NO dependence on the current ε, the reported tangent is the
    ELASTIC operator C — CONSTANT, SYMMETRIC, POSITIVE-DEFINITE. That is the whole
    point: a factor-once SPD global stiffness, no consistent-tangent singularity at
    perfect plasticity / softening, and (paired with a softening law or explicit
    dynamics) no Newton stall. The cost is an O(Δt) consistency error that vanishes
    under step refinement, and a one-step lag at first yield (Δγ_n = 0 ⇒ σ̃ = elastic
    trial until the extrapolation has a nonzero increment to carry forward).

This mirrors the fork's existing damage IMPL-EX (`returnMapDamaged`, dScaleOverride =
D_n + ΔD_n): a frozen extrapolated history quantity ⇒ SPD tangent, implicit value
committed. There it freezes the damage D; here it freezes the plastic multiplier Δγ.

NB for PLAIN (undamaged) hardening J2 the implicit tangent is already SPD, so the
robustness payoff is for softening / explicit use; the value on its own is the
constant factor-once tangent. Verified in tests/test_ladrunoJ2_implex_reference.py.
"""
from __future__ import annotations

import numpy as np

import ladrunoj2_reference as lr

_I = np.eye(3)


def elastic_tangent_6x6(p: lr.Params) -> np.ndarray:
    """Constant isotropic elastic modulus C in 3D Voigt {00,11,22,01,12,20} with
    ENGINEERING shear (the order/convention the kernel + element use). This is the
    full IMPL-EX material tangent — SPD and step-independent."""
    K, G = p.K, p.G
    lam = K - 2.0 * G / 3.0
    C = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            C[i, j] = lam
        C[i, i] += 2.0 * G
    for i in range(3, 6):
        C[i, i] = G               # engineering shear: τ = G·γ
    return C


def _elastic_stress(p: lr.Params, eps: np.ndarray, epsP: np.ndarray) -> np.ndarray:
    """σ = C:(ε − εᵖ) in 3×3 tensor form (εᵖ deviatoric)."""
    K, G = p.K, p.G
    tr = np.trace(eps)
    edev = eps - (tr / 3.0) * _I
    return K * tr * _I + 2.0 * G * (edev - epsP)


class ImplexLadrunoJ2:
    """Stateful material point with the IMPL-EX stress/tangent over the verified
    implicit return map. Prescribed-strain driven (setTrialStrain), mirroring how an
    OpenSees NDMaterial is used at a Gauss point.

    State committed each step: the plastic history (εᵖ, ε̄ᵖ, α) from the IMPLICIT
    return map, plus the IMPL-EX extrapolation history (Δγ_n, N_n)."""

    def __init__(self, p: lr.Params):
        self.p = p
        self.epsP = np.zeros((3, 3))
        self.ebarP = 0.0
        self.alpha = [np.zeros((3, 3)) for _ in range(p.nBack)]
        self.dGamma_n = 0.0                  # committed multiplier increment
        self.N_n = np.zeros((3, 3))          # committed unit flow direction (frozen)
        self._res = None                     # staged implicit result
        self._sigma_implex = None

    def setTrialStrain(self, eps_total: np.ndarray, dt_ratio: float = 1.0):
        """Stages the step. dt_ratio = Δt_{n+1}/Δt_n (1.0 for uniform stepping).
        Returns the IMPL-EX (explicit) Cauchy stress."""
        # 1) implicit pass — committed history + next-step extrapolation source
        self._res = lr.return_map(self.p, eps_total, self.epsP, self.ebarP, self.alpha)

        # 2) IMPL-EX explicit stress with frozen extrapolated multiplier + committed N
        dG_tilde = self.dGamma_n * dt_ratio
        epsP_tilde = self.epsP + dG_tilde * self.N_n
        self._sigma_implex = _elastic_stress(self.p, eps_total, epsP_tilde)
        return self._sigma_implex

    def getStress(self) -> np.ndarray:
        """IMPL-EX (explicit) Cauchy stress — what the global solver consumes."""
        return self._sigma_implex

    def getStressImplicit(self) -> np.ndarray:
        """The fully-implicit Cauchy stress at the same trial strain (for comparison /
        the value a fully-implicit material would report)."""
        return self._res["stress"]

    def getTangent(self) -> np.ndarray:
        """IMPL-EX material tangent = constant elastic C (6×6, SPD)."""
        return elastic_tangent_6x6(self.p)

    def commitState(self):
        """Commit the IMPLICIT history (so the committed state is identical to a pure
        implicit run) and refresh the extrapolation history (Δγ_n, N_n)."""
        dG = self._res["dGamma"]
        if self._res["plastic"] and dG > 0.0:
            self.N_n = (self._res["epsP"] - self.epsP) / dG     # unit deviatoric flow
        else:
            self.N_n = np.zeros((3, 3))
        self.dGamma_n = dG
        self.epsP = self._res["epsP"].copy()
        self.ebarP = self._res["ebarP"]
        self.alpha = [a.copy() for a in self._res["alpha"]]
