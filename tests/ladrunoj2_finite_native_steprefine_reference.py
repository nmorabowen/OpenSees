"""Step-refinement reference for the FINITE-STRAIN-NATIVE combined-hardening J2
(co-rotating backstress) — closes the ADR-16 P1 caveat.

The open P1 item (ADR Ladruno_implementation/16_finite_native_j2_adr.md, "Residual
open item"): the objectivity tests prove the co-rotated backstress is *frame
indifferent*, but they do not pin the *magnitude accuracy* of the Armstrong–Frederick
evolution when a step carries BOTH a large rotation AND a large stretch. There the
shipped incremental transport

      α̃_{n,k} = R_Δ · α_{n,k} · R_Δᵀ ,   R_Δ = polar(f_Δ) ,   f_Δ = F·F_nⁱⁿᵛ

differs at higher order from the exact §14.11 exponential-map (continuous-spin)
transport. The instruction: *verify by step-refinement against a fine-increment
reference; only if a discrepancy appears, swap to the exp-map transport.*

This module is that verification, numpy-only (the numpy oracle is "enough" per the
ADR). It drives the SAME small-strain return map (`ladrunoj2_reference.return_map`)
and the SAME Hencky kinematics (`logstrain_reference`) the C++ material uses, along a
path with a CONSTANT velocity gradient L = D + W whose symmetric stretch-rate D and
skew spin W do NOT commute (non-coaxial). That is the canonical case where the
rotation transports genuinely differ: F(s) = exp(s·L), so over a coarse step the
relative gradient is f_Δ = exp(Δs·L) and

      polar(f_Δ) = polar(exp(Δs L)) ≠ exp(Δs W)

(they agree only to first order in Δs). A fixed-axis stretch would be degenerate —
all U(s) commute ⇒ every transport coincides — so it is deliberately avoided here.

Two transport laws are provided so the step-refinement can compare them:

  transport="polar"  — R_Δ = polar(f_Δ) over the whole coarse step. This is EXACTLY
                       what the shipped `NativeFiniteJ2` / C++ `LadrunoJ2Finite` does;
                       `run(N, "polar")` is asserted bit-equal to `NativeFiniteJ2`.

  transport="exact"  — the backstress rotation over a coarse step is the time-ordered
                       product of M micro-rotations polar(f_Δ,micro) taken along the
                       REAL path. As M→∞ this converges to exp(Δs W), the exact
                       continuous-spin (vorticity exp-map) transport that polar(f_Δ)
                       approximates to first order. Only the backstress *rotation*
                       differs from "polar"; the elastic predictor bᵉᵗʳ = f_Δ bᵉ_n f_Δᵀ
                       and the return map are identical, so this isolates the
                       transport term under study.

Findings (see tests/test_ladrunoJ2_finite_native_steprefine.py):
  * the shipped polar scheme CONVERGES under refinement at observed order ≈ 1
    (consistent — no O(1) discrepancy in the AF magnitude);
  * the polar and exact-transport schemes converge to the SAME limit, and their
    difference at coarse N is itself O(Δs) and small — so polar(f_Δ) and the exact
    exp-map transport give the same physical answer ⇒ no swap needed.
"""
from __future__ import annotations

import numpy as np

import ladrunoj2_reference as lr
import logstrain_reference as ls
import ladrunoj2_finite_native_reference as nat

_I3 = np.eye(3)


# --------------------------------------------------------------------------- #
#  numpy-only matrix exponential (scaling-and-squaring + Taylor).             #
#  Adequate for the |s·L| ≲ 2 we use; verified against repeated-squaring      #
#  consistency in the tests. No scipy in the test env.                        #
# --------------------------------------------------------------------------- #
def expm(A: np.ndarray, terms: int = 24) -> np.ndarray:
    nrm = float(np.abs(A).sum())
    sq = max(0, int(np.ceil(np.log2(nrm + 1.0e-300))))
    B = A / (2.0 ** sq)
    E = np.eye(A.shape[0])
    term = np.eye(A.shape[0])
    for k in range(1, terms + 1):
        term = term @ B / k
        E = E + term
    for _ in range(sq):
        E = E @ E
    return E


# --------------------------------------------------------------------------- #
#  A severe combined large-rotation + large-stretch path with NON-COAXIAL     #
#  stretch-rate D and spin W (so the transport choice actually matters).      #
#  F(s) = exp(s·L), L = D + W, s ∈ [0, 1].                                     #
#    • W gives ≈2.0 rad (~115°) of rotation about a tilted axis,              #
#    • D is deviatoric (drives yield) with principal axes not aligned to W.   #
# --------------------------------------------------------------------------- #
_D = np.array([[0.25, 0.10, 0.04],
               [0.10, -0.15, 0.06],
               [0.04, 0.06, -0.10]])           # symmetric, traceless (isochoric rate)

_AXIS = np.array([0.2, -0.3, 1.0])
_AXIS = _AXIS / np.linalg.norm(_AXIS)
_ANGLE = 2.0                                   # total spin angle over s ∈ [0,1]
_W = _ANGLE * np.array([[0.0, -_AXIS[2], _AXIS[1]],
                        [_AXIS[2], 0.0, -_AXIS[0]],
                        [-_AXIS[1], _AXIS[0], 0.0]])   # skew

_L = _D + _W


def F_of_s(s: float) -> np.ndarray:
    """Deformation gradient at path parameter s ∈ [0, 1] (constant L = D + W)."""
    return expm(s * _L)


# --------------------------------------------------------------------------- #
#  Drivers                                                                     #
# --------------------------------------------------------------------------- #
# three-term Chaboche AF + Voce+linear iso (same params as the oracle test)
PARAMS = lr.Params(K=1.5e3, G=7.0e2, sig0=10.0, Qinf=6.0, bIso=25.0, Hiso=40.0,
                   C=[600.0, 350.0, 120.0], gam=[120.0, 60.0, 8.0])


def _transport_rotation(s_a, s_b, transport, M):
    """Rotation used to push the backstress forward over the coarse step a→b.

    "polar": single polar(f_Δ) of the whole coarse step (shipped behaviour).
    "exact": time-ordered product of M micro polar-rotations along the real path
             (→ exp(Δs W), the continuous-spin / exp-map transport, as M→∞)."""
    F_a = F_of_s(s_a)
    if transport == "polar":
        return nat.polar_rotation(F_of_s(s_b) @ np.linalg.inv(F_a))
    if transport == "exact":
        R = _I3.copy()
        F_prev = F_a
        for i in range(1, M + 1):
            s_i = s_a + (s_b - s_a) * i / M
            F_cur = F_of_s(s_i)
            R = nat.polar_rotation(F_cur @ np.linalg.inv(F_prev)) @ R
            F_prev = F_cur
        return R
    raise ValueError(transport)


def run(N: int, transport: str = "polar", M: int = 64, params: lr.Params = PARAMS):
    """Integrate F(s) over N uniform steps. Returns the committed final state.

    transport="polar" reproduces the shipped NativeFiniteJ2 exactly. The elastic
    predictor and return map are identical for both transports; only the backstress
    push-forward rotation differs."""
    Be_n = _I3.copy()
    Fn = _I3.copy()
    ebarP = 0.0
    alpha = [np.zeros((3, 3)) for _ in range(params.nBack)]
    sigma = None
    n_plastic = 0

    for n in range(1, N + 1):
        s_a = (n - 1) / N
        s_b = n / N
        F = F_of_s(s_b)

        fD = F @ np.linalg.inv(Fn)
        Be_tr = fD @ Be_n @ fD.T
        eps_tr = ls.hencky_strain(Be_tr)

        R = _transport_rotation(s_a, s_b, transport, M)
        alpha_pf = [R @ a @ R.T for a in alpha]

        res = lr.return_map(params, eps_tr, np.zeros((3, 3)), ebarP, alpha_pf)
        J = float(np.linalg.det(F))
        sigma = res["stress"] / J
        eps_e = eps_tr - res["epsP"]
        Be_n = ls.iso_function(2.0 * eps_e, np.exp)
        alpha = res["alpha"]
        ebarP = res["ebarP"]
        Fn = F.copy()
        if res["plastic"]:
            n_plastic += 1

    alpha_tot = sum(alpha) if alpha else np.zeros((3, 3))
    return dict(sigma=sigma, alpha=alpha, ebarP=ebarP, alpha_total=alpha_tot,
                alpha_norm=float(np.sqrt(np.tensordot(alpha_tot, alpha_tot))),
                n_plastic=n_plastic)


def state_vector(result: dict) -> np.ndarray:
    """Flatten the comparable part of a result (Cauchy σ, all backstresses, ε̄ᵖ)
    into one vector for norm-based convergence measurement."""
    parts = [result["sigma"].ravel()]
    for a in result["alpha"]:
        parts.append(a.ravel())
    parts.append(np.array([result["ebarP"]]))
    return np.concatenate(parts)
