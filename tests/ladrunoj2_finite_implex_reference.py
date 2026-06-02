"""IMPL-EX over the FINITE-STRAIN-NATIVE combined-hardening J2 — the numpy oracle
the C++ `nDMaterial LadrunoJ2Finite -implex` must match.

It layers the classic Δγ-extrapolation IMPL-EX (see ladrunoj2_implex_reference) on the
co-rotated finite return map (see ladrunoj2_finite_native_reference): the IMPLICIT
co-rotated return map runs every step and is committed (so the finite plastic history
bᵉ_n, ε̄ᵖ_n, α_n is byte-identical to a fully-implicit NativeFiniteJ2 run), while the
stress reported to the solver freezes the extrapolated multiplier Δγ̃ and the
committed log-strain flow direction N_n:

    εᵉ̃ = εᵉᵗʳ − Δγ̃ · Ñ_n,   Ñ_n = R_Δ N_n R_Δᵀ,   R_Δ = polar(f_Δ)
    τ̃  = elastic-D : εᵉ̃           (Kirchhoff; linear elastic in log strain)
    σ̃  = τ̃ / J

Two things make this clean and consistent with the v2 co-rotation:
  * the flow direction N_n is a spatial deviatoric tensor and is co-rotated by the
    SAME R_Δ = polar(f_Δ) as the backstress — so σ̃ is frame-indifferent;
  * the backstress α drops OUT of the explicit stress entirely (σ̃ depends only on
    εᵉᵗʳ and Δγ̃Ñ_n), so the IMPL-EX material tangent is the constant SPD elastic
    log-strain operator — no co-rotation "channel B", no plastic h term.

Verified (tests/test_ladrunoJ2_finite_implex_reference.py): committed history ==
implicit; σ̃ → implicit under F-step refinement; objective under superposed rotation;
reduces to the small-strain implex at tiny strain.
"""
from __future__ import annotations

import numpy as np

import ladrunoj2_reference as lr
import logstrain_reference as ls
import ladrunoj2_finite_native_reference as nat
import ladrunoj2_implex_reference as ix

_I3 = np.eye(3)


class ImplexNativeFiniteJ2:
    """Finite-strain-native combined-hardening J2 with the IMPL-EX stress/tangent.
    setTrialF(F)-driven, owns the finite state + the IMPL-EX extrapolation history."""

    def __init__(self, params: lr.Params):
        self.p = params
        self.Be_n = _I3.copy()
        self.Fn = _I3.copy()
        self.ebarP = 0.0
        self.alpha = [np.zeros((3, 3)) for _ in range(params.nBack)]
        self.dGamma_n = 0.0
        self.N_n = np.zeros((3, 3))          # committed log-strain flow direction (spatial)
        self._stage = None
        self._sigma_implex = None
        self._sigma_implicit = None

    def setTrialF(self, F: np.ndarray, dt_ratio: float = 1.0):
        """Returns the IMPL-EX (explicit) Cauchy stress. Stages state; commit() to keep."""
        fD = F @ np.linalg.inv(self.Fn)
        Be_tr = fD @ self.Be_n @ fD.T
        eps_tr = ls.hencky_strain(Be_tr)

        R = nat.polar_rotation(fD)
        alpha_pf = [R @ a @ R.T for a in self.alpha]
        N_pf = R @ self.N_n @ R.T            # co-rotate the committed flow direction

        # implicit co-rotated return map (elastic-trial form, epsP_n = 0) — the
        # committed truth + the next extrapolation source
        res = lr.return_map(self.p, eps_tr, np.zeros((3, 3)), self.ebarP, alpha_pf)
        J = float(np.linalg.det(F))
        self._sigma_implicit = res["stress"] / J
        self._plastic = res["plastic"]

        # IMPL-EX explicit stress: elastic log-strain stress at the frozen plastic strain
        dG_tilde = self.dGamma_n * dt_ratio
        epsP_tilde = dG_tilde * N_pf         # epsP_n = 0 in the elastic-trial form
        tau_tilde = ix._elastic_stress(self.p, eps_tr, epsP_tilde)   # Kirchhoff τ̃
        self._sigma_implex = tau_tilde / J

        # stage commit (implicit history) + refresh extrapolation (Δγ_n, N_n)
        eps_e = eps_tr - res["epsP"]
        Be_new = ls.iso_function(2.0 * eps_e, np.exp)
        dG = res["dGamma"]
        N_new = (res["epsP"] / dG) if (res["plastic"] and dG > 0.0) else np.zeros((3, 3))
        self._stage = (Be_new, res["alpha"], res["ebarP"], F.copy(), dG, N_new)
        return self._sigma_implex

    def getStress(self):
        return self._sigma_implex

    def getStressImplicit(self):
        return self._sigma_implicit

    def getTangent6x6(self):
        """IMPL-EX material tangent — the constant SPD elastic operator (the
        constitutive part of the log-strain spatial tangent; no plastic h, no
        co-rotation channel-B term)."""
        return ix.elastic_tangent_6x6(self.p)

    def commit(self):
        self.Be_n, self.alpha, self.ebarP, self.Fn, self.dGamma_n, self.N_n = self._stage
