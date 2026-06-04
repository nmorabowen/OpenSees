"""Reference (oracle) for the FINITE-STRAIN-NATIVE combined-hardening J2 with a
CO-ROTATING backstress — the de Souza Neto §14.11 "v2".

This is the numpy oracle the future C++ FiniteStrainNDMaterial-subclass J2 must
match, written BEFORE the C++ (the same way logstrain_reference.py preceded the
C++ LogStrainNDMaterial). It de-risks the one new ingredient: the spatial
push-forward of the kinematic backstress.

The problem it fixes — why the v1 LogStrain-over-LadrunoJ2 path is non-objective
for combined hardening under large rotation:
  * The log-strain (MATISU) algorithm co-rotates the ELASTIC state automatically,
    bᵉ_tr = f_Δ bᵉ_n f_Δᵀ, so the deviatoric stress s emerges in the CURRENT
    (rotated) frame.
  * Isotropic hardening depends only on ‖s‖ and the scalar ε̄ᵖ — both
    rotation-invariant — so it is exact.
  * Kinematic hardening uses the relative stress ξ = s − α. In v1 the backstress
    α lives inside the unchanged small-strain inner in a FIXED frame; the wrapper
    has no channel to rotate it. So ‖s − α‖ subtracts two tensors in different
    frames ⇒ wrong yield, wrong flow ⇒ a superposed rigid rotation is NOT
    stress-free. (A framework limit, not a wiring bug — the direct chain shows it.)

The fix (this module): own α as a SPATIAL tensor and push it forward by the
incremental material rotation R_Δ = polar(f_Δ) each step, BEFORE the return map:
      α̃_{n,k} = R_Δ · α_{n,k} · R_Δᵀ
Then run the SAME verified small-strain return map (ladrunoj2_reference.return_map)
in the rotated frame. For a pure rigid increment f_Δ = Q this gives α̃ = Q α Qᵀ,
which exactly cancels the rotation of s ⇒ frame indifference restored. As
f_Δ → I, R_Δ → I, so it reduces to the small-strain / non-rotating case exactly.

Verified (tests/test_ladrunoJ2_finite_native_reference.py): objective under a
superposed finite rotation AND along a continuously rotating+stretching plastic
path to ~1e-12; reduces to v1 when there is no material rotation.

The return map itself is reused VERBATIM — the C++ v2 calls LadrunoJ2Kernel.h
directly, exactly as the small-strain LadrunoJ2 does. No OpenSees / openseespy
dependency — numpy only.
"""
from __future__ import annotations

import numpy as np

import ladrunoj2_reference as lr
import logstrain_reference as ls

_I3 = np.eye(3)


def polar_rotation(F: np.ndarray) -> np.ndarray:
    """Rotation part R of the polar decomposition F = R U (R orthogonal, U SPD),
    via the SVD F = W Σ Vᵀ ⇒ R = W Vᵀ. For det F > 0 this is a proper rotation;
    a reflection guard is kept for safety. R is the incremental MATERIAL rotation
    used to objectively transport the backstress."""
    W, _s, Vt = np.linalg.svd(F)
    R = W @ Vt
    if np.linalg.det(R) < 0.0:                 # guard (shouldn't trip for det F>0)
        W = W.copy(); W[:, -1] *= -1.0; R = W @ Vt
    return R


class NativeFiniteJ2:
    """Finite-strain combined-hardening J2 that OWNS its finite-strain state
    (bᵉ_n, ε̄ᵖ_n, and the backstresses α_{n,k} as SPATIAL tensors) and reuses the
    small-strain return map. Driven by setTrialF(F).

    corotate=False reproduces the v1 LogStrain-over-LadrunoJ2 behaviour (α left in
    a fixed frame) — kept so tests can exhibit the non-objectivity it cures.
    """

    def __init__(self, params: lr.Params, corotate: bool = True):
        self.p = params
        self.corotate = corotate
        self.Be_n = _I3.copy()
        self.Fn = _I3.copy()
        self.ebarP = 0.0
        self.alpha = [np.zeros((3, 3)) for _ in range(params.nBack)]
        self._stage = None

    def setTrialF(self, F: np.ndarray):
        """Returns (Cauchy σ, plastic_flag). Stages state; call commit() to keep."""
        fD = F @ np.linalg.inv(self.Fn)               # relative deformation gradient
        Be_tr = fD @ self.Be_n @ fD.T                 # (i)+(ii) elastic predictor
        eps_tr = ls.hencky_strain(Be_tr)              # εᵉᵗʳ = ½ ln bᵉᵗʳ

        R = polar_rotation(fD) if self.corotate else _I3
        alpha_pf = [R @ a @ R.T for a in self.alpha]  # ← the §14.11 push-forward

        # (iii) the UNCHANGED small-strain return map, in the rotated frame.
        # Elastic-trial form: epsP_n = 0 so the input IS the elastic trial.
        res = lr.return_map(self.p, eps_tr, np.zeros((3, 3)), self.ebarP, alpha_pf)

        J = float(np.linalg.det(F))
        sigma = res["stress"] / J                     # (iv) Cauchy σ = τ / J
        eps_e = eps_tr - res["epsP"]                  # εᵉ_{n+1} = εᵉᵗʳ − Δεᵖ
        Be_new = ls.iso_function(2.0 * eps_e, np.exp) # bᵉ = exp[2 εᵉ]
        self._stage = (Be_new, res["alpha"], res["ebarP"], F.copy())
        return sigma, res["plastic"]

    def commit(self):
        self.Be_n, self.alpha, self.ebarP, self.Fn = self._stage

    # ------------------------------------------------------------------ #
    #  Consistent material tangent  ∂σ/∂F  at the current committed state #
    # ------------------------------------------------------------------ #
    # The Cauchy stress depends on the trial F through TWO channels:
    #   (A) the trial elastic strain  εᵉᵗʳ = ½ ln(f_Δ Be_n f_Δᵀ)
    #   (B) the co-rotated backstress  α̃ = R(f_Δ) α_n R(f_Δ)ᵀ,  R = polar(f_Δ)
    # The shipped log-strain spatial tangent (LogStrainKernel/spatial_tangent_a)
    # captures only (A). Channel (B) is a PURE CONSTITUTIVE sensitivity (α̃→σ),
    # small (~2–4e-4 relative) and frame-objective, but it sits at the tangent-FD
    # test tolerance. The C++ material recovers it EXACTLY-to-FD by perturbing only
    # R (≈9 extra scalar-Newton return-map calls) — NO kernel re-derivation. These
    # helpers are the numpy target for that recipe.
    def _cauchy_trial(self, F, R_override=None):
        """Cauchy σ for a trial F at the current committed (n) state; no staging."""
        fD = F @ np.linalg.inv(self.Fn)
        Be_tr = fD @ self.Be_n @ fD.T
        eps_tr = ls.hencky_strain(Be_tr)
        R = polar_rotation(fD) if (R_override is None and self.corotate) else \
            (R_override if R_override is not None else _I3)
        alpha_pf = [R @ a @ R.T for a in self.alpha]
        res = lr.return_map(self.p, eps_tr, np.zeros((3, 3)), self.ebarP, alpha_pf)
        return res["stress"] / float(np.linalg.det(F))

    def tangent_dsigma_dF(self, F, freeze_R=False, h=1e-7):
        """Full consistent ∂σ_ij/∂F_kl by central FD. freeze_R=True holds the
        backstress rotation fixed (= channel A only)."""
        A = np.zeros((3, 3, 3, 3))
        R_fix = polar_rotation(F @ np.linalg.inv(self.Fn)) if freeze_R else None
        for k in range(3):
            for l in range(3):
                Fp = F.copy(); Fp[k, l] += h
                Fm = F.copy(); Fm[k, l] -= h
                A[:, :, k, l] = (self._cauchy_trial(Fp, R_fix)
                                 - self._cauchy_trial(Fm, R_fix)) / (2.0 * h)
        return A

    def channelB_dsigma_dF(self, F, h=1e-7):
        """Channel-B-only ∂σ/∂F: perturb R = polar(f_Δ) while holding εᵉᵗʳ and J
        fixed at the base F (the term the C++ adds via R-perturbation)."""
        fD0 = F @ np.linalg.inv(self.Fn)
        eps0 = ls.hencky_strain(fD0 @ self.Be_n @ fD0.T)
        J0 = float(np.linalg.det(F))

        def sig(R):
            res = lr.return_map(self.p, eps0, np.zeros((3, 3)), self.ebarP,
                                [R @ a @ R.T for a in self.alpha])
            return res["stress"] / J0

        B = np.zeros((3, 3, 3, 3))
        for k in range(3):
            for l in range(3):
                Fp = F.copy(); Fp[k, l] += h
                Fm = F.copy(); Fm[k, l] -= h
                B[:, :, k, l] = (sig(polar_rotation(Fp @ np.linalg.inv(self.Fn)))
                                 - sig(polar_rotation(Fm @ np.linalg.inv(self.Fn)))) \
                    / (2.0 * h)
        return B
