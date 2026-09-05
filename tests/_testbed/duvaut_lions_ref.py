"""ADR 90 / phase P0 leg A0 — Duvaut-Lions viscoplastic regularization ORACLE (numpy only).

Purpose
-------
The fork promised TIMs a *generic* `NDMaterial` Duvaut-Lions wrapper.  A generic wrapper can only
see the inner material's `getStress()`, so it necessarily implements a **two-track (TT)** model:

    inner (inviscid) material evolves on the TOTAL strain path;
    sigma_tr = sigma_vp,n + C_e : deps
    beta     = dt / (tau + dt)
    sigma_vp = (1-beta) sigma_tr + beta * inner.getStress()
    D_vp     = (1-beta) C_e      + beta * inner.getTangent()

**True Duvaut-Lions (TDL)** (Simo & Hughes, *Computational Inelasticity* Sec. 2.7;
Simo, Kennedy & Govindjee 1988) instead projects the *viscoplastic* trial state
(sigma_vp,n + C_e deps, q_vp,n) onto the yield surface and relaxes BOTH the stress and the
internal variables toward that projection:

    (sigma_bar, q_bar) = P( sigma_vp,n + C_e:deps , q_vp,n )        # closest-point projection
    sigma_vp = (sigma_tr + (dt/tau) sigma_bar)/(1 + dt/tau) = (1-beta) sigma_tr + beta sigma_bar
    q_vp     = (q_n      + (dt/tau) q_bar    )/(1 + dt/tau) = (1-beta) q_n      + beta q_bar

For PERFECT plasticity the two coincide.  For hardening / softening / state-dependent laws they
do not, and *nobody had measured whether TT still regularizes localization WIDTH*.  That is what
this module measures (planning brief OQ3, decision D2, risk R1-BLOCKING).

Contents
--------
1. `dl_beta`                    — the fork's blend factor + the `tau<=0 or dt<=0 => beta=1` rule
                                  (`LEDGER_quirks.md` "LadrunoConcrete3D -eta"; kernel
                                  `SRC/material/nD/LadrunoConcrete3DKernel.h:1486-1488`).
2. 1-D elastic-plastic point model with linear hardening/softening H (H<0 allowed) + residual
   floor; exact return map and consistent tangent; fully vectorized over points.
3. 3-D J2 point model, linear isotropic hardening, Voigt-6 in the fork ordering
   {xx, yy, zz, xy, yz, zx} with ENGINEERING shear strain; exact radial return + Simo-Hughes
   algorithmic tangent.
4. `variant="tt" | "tdl" | "inviscid"` steps over both point models.
5. PV1..PV6 falsification battery, cloned from the shipped `-eta` gate
   (`concrete3d_ref.run_p3_eta_gate`), for BOTH variants and BOTH point models.
6. `run_tt_vs_tdl_point` — quantified TT-vs-TDL study on hardening AND softening laws.
7. `a0_bar` — the A0 1-D quasi-static softening bar (Bazant / de Borst mesh-dependence case) with
   three threshold-free-to-threshold-heavy band-width definitions, and `run_a0_sweep` driving it.

Units are MPa / mm / N throughout (E = 20000 MPa, L = 100 mm, A = 1 mm^2), so a "load" is in N.
Time is the analysis pseudo-time; De = tau / T is the Deborah number of a run.

No OpenSees import: this file is pure numpy and runs on a box with no built fork.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "dl_beta",
    "j1d_params", "j1d_yield", "j1d_return", "j1d_state", "j1d_step",
    "j2_params", "j2_return", "j2_state", "j2_step", "j2_elastic_C",
    "dl_1d_ramp", "dl_1d_analytic",
    "run_pv_gate",
    "run_tt_vs_tdl_point",
    "band_widths", "a0_bar", "run_a0_sweep",
]

SQ23 = np.sqrt(2.0 / 3.0)


# ===========================================================================
# 0. The relaxation factor and the fork's dt / tau convention
# ===========================================================================
def dl_beta(tau, dt):
    """beta = dt/(tau+dt) in [0,1].

    FORK CONVENTION (deliberate, `LEDGER_quirks.md` ~line 1612): `tau <= 0` OR `dt <= 0` gives
    `beta = 1`, i.e. the INVISCID limit -- NOT the mathematical `dt -> 0` elastic limit.  A missing
    or zero pseudo-time increment must never silently turn a material elastic.  The elastic limit
    is reached only through a LARGE tau at finite dt (beta -> 0), which is the correct knob.
    """
    tau = float(tau)
    dt = float(dt)
    if tau <= 0.0 or dt <= 0.0:
        return 1.0
    return dt / (tau + dt)


# ===========================================================================
# 1. 1-D elastic-plastic point model, linear hardening/softening + residual floor
# ===========================================================================
def j1d_params(E=20000.0, sigY=20.0, H=0.0, res_frac=0.0, H_res=None,
               law="linear", K_inf=None, alpha_f=None, E_beta=0.0):
    """1-D rate-independent law.  Two hardening laws:

      law="linear"  K(alpha) = max(sigY + H*alpha, res_frac*sigY + H_res*alpha)  -- H<0 allowed
                    (softening).  The residual branch may carry a SMALL POSITIVE H_res: a
                    perfectly-plastic element in a 1-D chain makes the assembled structural
                    tangent exactly SINGULAR (the chain extends at constant force), which kills
                    Newton in the A0 bar.  H_res > 0 removes that without touching the
                    pre-residual branch, and the return map stays EXACT (two linear branches, the
                    active one chosen at the UPDATED alpha).  DEFAULT H_res = 0 (a flat floor);
                    `a0_bar` sets it to E/2e5 explicitly.  Do NOT default it non-zero: with
                    res_frac = 0 the residual line `0 + H_res*alpha` eventually CROSSES sigY and
                    the "perfectly plastic" bar silently starts hardening (it broke PV3 once).
      law="exp"     K(alpha) = K_inf + (sigY - K_inf)*exp(-alpha/alpha_f)
                    (K_inf < sigY => exponential SOFTENING toward a residual K_inf;
                     K_inf > sigY => saturation HARDENING.  This law is NONLINEAR in alpha, which
                     is exactly what breaks the TT == TDL identity proved for the linear law --
                     see `run_tt_vs_tdl_point`.)

    `sigY` may be an ARRAY (one entry per material point) -- that is how the A0 bar carries its
    imperfection.  Well-posedness needs E + K'(alpha) > 0 everywhere; for the linear law that is
    |H| < E, for the exponential law it is E > (sigY - K_inf)/alpha_f.

    `E_beta` makes the ELASTIC OPERATOR STATE-DEPENDENT (leg P0b-a, the constitutive critic's
    BLOCKING finding):  E(sigma) = E * (1 + E_beta * |sigma| / sigY_ref).  This is the 1-D
    caricature of SANISAND's pressure-dependent G(p), K(p) and of any hypoelastic inner.
    `E_beta = 0.0` (the default) returns `p["E"]` by an EARLY RETURN, so every pre-existing
    result is bit-identical.
    """
    sY = np.asarray(sigY, float)
    p = {"E": float(E), "sigY": sY, "H": float(H), "sig_res": res_frac * sY, "law": law,
         "H_res": float(0.0 if H_res is None else H_res),
         "E_beta": float(E_beta), "sigY_ref": float(np.max(sY))}
    if law == "linear":
        if not np.all(np.abs(H) < E):
            raise ValueError("need |H| < E for a well-posed 1-D return map")
    elif law == "exp":
        if K_inf is None or alpha_f is None:
            raise ValueError("law='exp' needs K_inf and alpha_f")
        Ki = np.asarray(K_inf, float) * np.ones_like(sY)
        p["K_inf"] = Ki
        p["alpha_f"] = float(alpha_f)
        if not np.all(np.abs(sY - Ki) / alpha_f < E):
            raise ValueError("need E > |sigY - K_inf|/alpha_f for a well-posed return map")
    else:
        raise ValueError(f"unknown law {law!r}")
    return p


def j1d_Emod(sig, p):
    """The state-dependent elastic modulus E(sigma).  Early-returns the constant `p["E"]` when
    `E_beta == 0`, so the constant-E path is bit-identical to the pre-P0b oracle."""
    if p.get("E_beta", 0.0) == 0.0:
        return p["E"]
    return p["E"] * (1.0 + p["E_beta"] * np.abs(np.asarray(sig, float)) / p["sigY_ref"])


def j1d_yield(alpha, p):
    if p["law"] == "linear":
        return np.maximum(p["sigY"] + p["H"] * alpha, p["sig_res"] + p["H_res"] * alpha)
    return p["K_inf"] + (p["sigY"] - p["K_inf"]) * np.exp(-alpha / p["alpha_f"])


def j1d_hard(alpha, p):
    """K'(alpha).  Linear law: H on the softening branch, H_res on the residual branch."""
    if p["law"] == "linear":
        return np.where(p["sigY"] + p["H"] * alpha > p["sig_res"] + p["H_res"] * alpha,
                        p["H"], p["H_res"])
    return -(p["sigY"] - p["K_inf"]) / p["alpha_f"] * np.exp(-alpha / p["alpha_f"])


def j1d_return(sig_tr, alpha_n, p, newton_tol=1.0e-14, newton_max=50, E=None):
    """EXACT rate-independent return map + consistent tangent.  Vectorized.

    linear law -- two admissible LINEAR branches (K is their upper envelope), the active one
      selected at the UPDATED alpha, so the map is exact including the step that crosses the kink:
      * softening branch  K = sigY + H*alpha            -> dgam = f1/(E+H),     Dep = E*H/(E+H)
      * residual branch   K = sig_res + H_res*alpha     -> dgam = f2/(E+H_res), Dep = E*H_res/(E+H_res)
    exp law    -- a local Newton on g(dgam) = |sig_tr| - E*dgam - K(alpha_n + dgam) = 0 (g is
      strictly decreasing because E + K' > 0), then Dep = E*K'/(E+K') at the UPDATED alpha.

    `E` overrides the elastic modulus for this step (the hypoelastic / state-dependent case, where
    the caller froze E at the committed stress).  `E = None` uses the constant `p["E"]`.

    Returns (sig_bar, alpha_bar, Dep, plastic).
    """
    E = p["E"] if E is None else E
    sig_tr = np.asarray(sig_tr, float)
    alpha_n = np.asarray(alpha_n, float)
    absig = np.abs(sig_tr)
    f = absig - j1d_yield(alpha_n, p)
    plastic = f > 0.0

    if p["law"] == "linear":
        H, sY0, sres, Hr = p["H"], p["sigY"], p["sig_res"], p["H_res"]
        dg_lin = (absig - (sY0 + H * alpha_n)) / (E + H)
        dg_res = (absig - (sres + Hr * alpha_n)) / (E + Hr)
        a1 = alpha_n + dg_lin
        lin_ok = (sY0 + H * a1) >= (sres + Hr * a1)
        dgam = np.where(plastic, np.where(lin_ok, dg_lin, dg_res), 0.0)
        Dep = np.where(plastic, np.where(lin_ok, E * H / (E + H), E * Hr / (E + Hr)), E)
    else:
        dgam = np.where(plastic, f / (E - j1d_hard(alpha_n, p)), 0.0)
        for _ in range(newton_max):
            a = alpha_n + dgam
            g = absig - E * dgam - j1d_yield(a, p)
            gp = -E - j1d_hard(a, p)
            step = np.where(plastic, -g / gp, 0.0)
            dgam = np.maximum(dgam + step, 0.0)
            if float(np.max(np.abs(step))) < newton_tol * max(float(np.max(absig)), 1.0):
                break
        Kp = j1d_hard(alpha_n + dgam, p)
        Dep = np.where(plastic, E * Kp / (E + Kp), E)

    sig = sig_tr - E * dgam * np.sign(sig_tr)
    alpha = alpha_n + dgam
    return sig, alpha, Dep, plastic


def j1d_state(n=1):
    """Committed state.  `sig`/`eps` are the VISCOPLASTIC (wrapper) stress and total strain;
    `sig_in`/`alp_in` are the TT inner material's own track; `alp` is the TDL / inviscid internal
    variable.  One dict carries all three variants so the A0 bar driver is variant-agnostic."""
    z = np.zeros(int(n))
    return {"sig": z.copy(), "eps": z.copy(),
            "sig_in": z.copy(), "alp_in": z.copy(), "alp": z.copy()}


def j1d_step(state, deps, p, tau, dt, variant="tt"):
    """One stress update.  Does NOT mutate `state` (trial calls repeat from the committed state).

    Variants:
      "tt"       the GENERIC wrapper.  The predictor modulus is whatever `inner.getInitialTangent()`
                 returns, i.e. E evaluated at the INNER's committed stress -- the only elastic
                 operator a generic NDMaterial seam can obtain.
      "tt_vp"    the obvious REPAIR candidate: the same blend, but with the predictor modulus
                 evaluated at the WRAPPER's own committed stress.  This is NOT implementable in a
                 generic wrapper (it needs the function E(.), not one sampled number) -- it exists
                 here only to answer "would it even help?"  See `run_pressure_dependent_leg`.
      "tdl"      true Duvaut-Lions: project the viscoplastic trial state, relax sigma AND alpha.
      "inviscid" the rate-independent reference.

    Returns (sig, D, new_state, info).  `info["dvp"]` is the viscoplastic strain increment
    `C_e^-1 (sig_tr - sig_vp) = beta C_e^-1 (sig_tr - sig_bar)`, which equals the continuous
    `C_e^-1 (sig_vp - sig_bar) dt/tau` exactly and stays finite at tau = 0; `info["Dw"]` is the
    wrapper's own dissipation increment `sig_vp . dvp` (leg P0b-b).
    """
    beta = dl_beta(tau, dt)
    deps = np.asarray(deps, float)

    if variant == "tt":
        Ce = j1d_Emod(state["sig_in"], p)                 # == inner.getInitialTangent()
        Ce_in = Ce
    elif variant == "tt_vp":
        Ce = j1d_Emod(state["sig"], p)
        Ce_in = j1d_Emod(state["sig_in"], p)
    else:
        Ce = j1d_Emod(state["sig"], p)
        Ce_in = Ce
    sig_tr = state["sig"] + Ce * deps                     # viscoplastic elastic predictor

    if variant in ("tt", "tt_vp"):
        # TWO-TRACK: the inner material never sees tau; it walks the total strain path inviscidly.
        sig_tr_in = state["sig_in"] + Ce_in * deps
        sig_bar, alp_bar, Dep, plastic = j1d_return(sig_tr_in, state["alp_in"], p, E=Ce_in)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        D = (1.0 - beta) * Ce + beta * Dep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": sig_bar, "alp_in": alp_bar, "alp": state["alp"].copy()}
    elif variant == "tdl":
        # TRUE DL: project the VISCOPLASTIC trial state, relax sigma AND alpha toward it.
        sig_bar, alp_bar, Dep, plastic = j1d_return(sig_tr, state["alp"], p, E=Ce)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        D = (1.0 - beta) * Ce + beta * Dep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": state["sig_in"].copy(), "alp_in": state["alp_in"].copy(),
               "alp": (1.0 - beta) * state["alp"] + beta * alp_bar}
    elif variant == "inviscid":
        sig_bar, alp_bar, Dep, plastic = j1d_return(sig_tr, state["alp"], p, E=Ce)
        sig = sig_bar
        D = Dep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": sig_bar, "alp_in": alp_bar, "alp": alp_bar}
    else:
        raise ValueError(f"unknown variant {variant!r}")

    dvp = beta * (sig_tr - sig_bar) / Ce
    info = {"beta": beta, "sig_tr": sig_tr, "sig_bar": sig_bar, "Ce": Ce,
            "alp_bar": alp_bar, "Dep": Dep, "plastic": plastic,
            "dvp": dvp, "Dw": sig * dvp}
    return sig, D, new, info


# ===========================================================================
# 2. 3-D J2 point model (Voigt-6, fork ordering xx yy zz xy yz zx, engineering shear)
# ===========================================================================
def _eps_v2t(e):
    return np.array([[e[0], 0.5 * e[3], 0.5 * e[5]],
                     [0.5 * e[3], e[1], 0.5 * e[4]],
                     [0.5 * e[5], 0.5 * e[4], e[2]]], float)


def _sig_v2t(s):
    return np.array([[s[0], s[3], s[5]],
                     [s[3], s[1], s[4]],
                     [s[5], s[4], s[2]]], float)


def _sig_t2v(S):
    return np.array([S[0, 0], S[1, 1], S[2, 2], S[0, 1], S[1, 2], S[0, 2]], float)


def j2_params(E=20000.0, nu=0.2, sigY=20.0, H=0.0, res_frac=0.0):
    G = E / (2.0 * (1.0 + nu))
    K = E / (3.0 * (1.0 - 2.0 * nu))
    return {"E": E, "nu": nu, "G": G, "K": K, "sigY": float(sigY), "H": float(H),
            "sig_res": res_frac * float(sigY)}


def _j2_apply(dE, K, G, theta, thetabar, n):
    """Apply the (algorithmic) modulus to a STRAIN TENSOR increment; returns a stress tensor.
    C = K 1(x)1 + 2 G theta [I - 1/3 1(x)1] - 2 G thetabar n(x)n   (Simo & Hughes Box 3.2)."""
    tr = np.trace(dE)
    dev = dE - (tr / 3.0) * np.eye(3)
    out = K * tr * np.eye(3) + 2.0 * G * theta * dev
    if thetabar != 0.0 and n is not None:
        out = out - 2.0 * G * thetabar * float(np.sum(n * dE)) * n
    return out


def _j2_C_matrix(K, G, theta, thetabar, n):
    """Assemble the 6x6 Voigt matrix column by column from the tensor operator -- no hand-rolled
    Voigt bookkeeping, so the engineering-shear factors cannot be got wrong."""
    C = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6)
        e[j] = 1.0
        C[:, j] = _sig_t2v(_j2_apply(_eps_v2t(e), K, G, theta, thetabar, n))
    return C


def j2_elastic_C(p):
    return _j2_C_matrix(p["K"], p["G"], 1.0, 0.0, None)


def j2_return(sig_tr_v, alpha_n, p):
    """Radial return with linear isotropic hardening (+ residual floor) and the ALGORITHMIC
    tangent.  Returns (sig_bar_v, alpha_bar, C_ep 6x6, plastic)."""
    G, K, H, sY0, sres = p["G"], p["K"], p["H"], p["sigY"], p["sig_res"]
    S = _sig_v2t(np.asarray(sig_tr_v, float))
    ph = np.trace(S) / 3.0
    s = S - ph * np.eye(3)
    ns = float(np.sqrt(np.sum(s * s)))

    ky_n = max(sY0 + H * alpha_n, sres)
    f = ns - SQ23 * ky_n
    if f <= 0.0 or ns == 0.0:
        return np.asarray(sig_tr_v, float).copy(), float(alpha_n), j2_elastic_C(p), False

    dg = f / (2.0 * G + (2.0 / 3.0) * H)
    Hloc = H
    if sY0 + H * (alpha_n + SQ23 * dg) < sres:              # floored branch
        dg = (ns - SQ23 * sres) / (2.0 * G)
        Hloc = 0.0

    n = s / ns
    s_new = s - 2.0 * G * dg * n
    sig = s_new + ph * np.eye(3)
    alpha = float(alpha_n) + SQ23 * dg

    theta = 1.0 - 2.0 * G * dg / ns
    thetabar = 1.0 / (1.0 + Hloc / (3.0 * G)) - (1.0 - theta)
    C = _j2_C_matrix(K, G, theta, thetabar, n)
    return _sig_t2v(sig), alpha, C, True


def j2_state():
    z = np.zeros(6)
    return {"sig": z.copy(), "eps": z.copy(),
            "sig_in": z.copy(), "alp_in": 0.0, "alp": 0.0}


def j2_step(state, deps, p, tau, dt, variant="tt"):
    """J2 counterpart of `j1d_step`; identical algebra with C_e in place of E."""
    Ce = j2_elastic_C(p)
    beta = dl_beta(tau, dt)
    deps = np.asarray(deps, float)
    sig_tr = state["sig"] + Ce @ deps

    if variant == "tt":
        sig_tr_in = state["sig_in"] + Ce @ deps
        sig_bar, alp_bar, Cep, plastic = j2_return(sig_tr_in, state["alp_in"], p)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        D = (1.0 - beta) * Ce + beta * Cep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": sig_bar, "alp_in": alp_bar, "alp": state["alp"]}
    elif variant == "tdl":
        sig_bar, alp_bar, Cep, plastic = j2_return(sig_tr, state["alp"], p)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        D = (1.0 - beta) * Ce + beta * Cep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": state["sig_in"].copy(), "alp_in": state["alp_in"],
               "alp": (1.0 - beta) * state["alp"] + beta * alp_bar}
    elif variant == "inviscid":
        sig_bar, alp_bar, Cep, plastic = j2_return(sig_tr, state["alp"], p)
        sig, D = sig_bar, Cep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": sig_bar, "alp_in": alp_bar, "alp": alp_bar}
    else:
        raise ValueError(f"unknown variant {variant!r}")

    dvp = beta * np.linalg.solve(Ce, sig_tr - sig_bar)      # == C_e^-1 (sig_vp - sig_bar) dt/tau
    info = {"beta": beta, "sig_tr": sig_tr, "sig_bar": sig_bar, "Ce": Ce,
            "alp_bar": alp_bar, "Cep": Cep, "plastic": plastic,
            "dvp": dvp, "Dw": float(sig @ dvp)}
    return sig, D, new, info


# ===========================================================================
# 3. The closed-form 1-D overstress oracle (PV2 / PV3)
# ===========================================================================
def dl_1d_ramp(variant, E, sigY, eps_rate, tau, dt, n, H=0.0):
    """Constant-strain-rate 1-D driver from rest.  Returns (t, sigma)."""
    p = j1d_params(E=E, sigY=sigY, H=H)
    st = j1d_state(1)
    t = np.empty(n)
    s = np.empty(n)
    for k in range(n):
        sig, _, st, _ = j1d_step(st, np.array([E * eps_rate * dt / E]), p, tau, dt, variant)
        t[k] = (k + 1) * dt
        s[k] = float(sig[0])
    return t, s


def dl_1d_analytic(E, sigY, eps_rate, tau, t):
    """Continuous Duvaut-Lions solution, 1-D PERFECT plasticity, constant strain rate from rest:
    elastic sigma = E*eps_rate*t until t_y = sigY/(E*eps_rate), then
    sigma(t) = sigY + E*eps_rate*tau*(1 - exp(-(t-t_y)/tau))."""
    t = np.asarray(t, float)
    t_y = sigY / (E * eps_rate)
    return np.where(t <= t_y, E * eps_rate * t,
                    sigY + E * eps_rate * tau * (1.0 - np.exp(-(t - t_y) / tau)))


# ===========================================================================
# 4. PV1..PV6 falsification battery -- BOTH variants, BOTH point models
# ===========================================================================
def _nrm(a):
    return float(np.sqrt(np.sum(np.asarray(a, float) ** 2)))


def _pv1_1d(variant, H):
    """tau = 0 must be BYTE-identical to the inviscid material over a path that yields and softens."""
    p = j1d_params(E=20000.0, sigY=20.0, H=H, res_frac=0.02)
    st_v = j1d_state(1)
    st_i = j1d_state(1)
    worst = 0.0
    deps = np.array([4.0e-5])
    for k in range(200):
        d = deps if k < 160 else -deps            # load into softening, then unload
        sv, _, st_v, _ = j1d_step(st_v, d, p, 0.0, 1.0, variant)
        si, _, st_i, _ = j1d_step(st_i, d, p, 0.0, 1.0, "inviscid")
        worst = max(worst, abs(float(sv[0] - si[0])))
    return worst


def _pv1_j2(variant, H):
    p = j2_params(E=20000.0, nu=0.2, sigY=20.0, H=H, res_frac=0.02)
    st_v, st_i = j2_state(), j2_state()
    worst = 0.0
    d0 = np.array([6.0e-5, -1.0e-5, -1.0e-5, 3.0e-5, 0.0, 0.0])
    for k in range(200):
        d = d0 if k < 160 else -d0
        sv, _, st_v, _ = j2_step(st_v, d, p, 0.0, 1.0, variant)
        si, _, st_i, _ = j2_step(st_i, d, p, 0.0, 1.0, "inviscid")
        worst = max(worst, _nrm(sv - si))
    return worst


def _plastic_state_1d(p, variant, tau, dt, nadv=40, dstep=4.0e-5):
    """Advance a point WELL inside the plastic branch (so PV4/PV5/PV6 are not tautological)."""
    st = j1d_state(1)
    for _ in range(nadv):
        _, _, st, _ = j1d_step(st, np.array([dstep]), p, tau, dt, variant)
    return st


def _plastic_state_j2(p, variant, tau, dt, nadv=40):
    st = j2_state()
    d = np.array([6.0e-5, -1.0e-5, -1.0e-5, 3.0e-5, 0.0, 0.0])
    for _ in range(nadv):
        _, _, st, _ = j2_step(st, d, p, tau, dt, variant)
    return st


def run_pv_gate(verbose=True):
    """PV1..PV6 for BOTH variants over BOTH point models.  Mirrors the shipped `-eta` gate
    (`concrete3d_ref.run_p3_eta_gate`) so the C++ port inherits an already-argued battery."""
    res = {}

    # ---- PV1: tau = 0 byte-exact (== 0.0, not a tolerance) -----------------
    pv1 = {}
    for var in ("tt", "tdl"):
        for H in (0.0, 2000.0, -2000.0):            # perfect / hardening / softening
            pv1[f"{var}_1d_H{H:+.0f}"] = _pv1_1d(var, H)
            pv1[f"{var}_j2_H{H:+.0f}"] = _pv1_j2(var, H)
    res["PV1"] = pv1
    res["PV1_max"] = max(pv1.values())
    res["PV1_byte_exact"] = bool(res["PV1_max"] == 0.0)

    # ---- PV2 / PV3: the closed-form 1-D overstress oracle -------------------
    E, sigY, eps_rate, tau = 20000.0, 20.0, 1.0e-3, 0.5
    over_exact = E * eps_rate * tau
    res["PV3_target"] = over_exact
    pv3 = {}
    for var in ("tt", "tdl"):
        errs = []
        for dt in (1.0, 0.25, 0.05):
            n = int(60.0 / (eps_rate * E) / dt) + 4000
            _, s = dl_1d_ramp(var, E, sigY, eps_rate, tau, dt, n)
            errs.append(abs((s[-1] - sigY) - over_exact) / over_exact)
        pv3[var] = max(errs)
    res["PV3_rel_err"] = pv3
    res["PV3_exact"] = bool(max(pv3.values()) < 1.0e-10)

    t_y = sigY / (E * eps_rate)
    T2 = t_y + 2.0 * tau                        # MID-transient (not the steady state)
    pv2 = {}
    for var in ("tt", "tdl"):
        errs = []
        for dt in (0.5, 0.25, 0.125):
            n = int(round(T2 / dt))
            t, s = dl_1d_ramp(var, E, sigY, eps_rate, tau, dt, n)
            errs.append(abs(s[-1] - float(dl_1d_analytic(E, sigY, eps_rate, tau, t[-1]))))
        pv2[var] = {"errs": errs, "order": float(np.log(errs[0] / errs[-1]) / np.log(4.0))}
    res["PV2"] = pv2
    res["PV2_converges"] = bool(all(v["order"] > 0.8 and v["errs"][-1] < v["errs"][0]
                                    for v in pv2.values()))

    # ---- PV4: the update IS the (1-beta)trial + beta projection blend ------
    pv4 = {}
    tau4, dt4 = 0.3, 1.0
    b4 = dl_beta(tau4, dt4)
    for H in (2000.0, -2000.0):
        p1 = j1d_params(E=20000.0, sigY=20.0, H=H, res_frac=0.02)
        pj = j2_params(E=20000.0, nu=0.2, sigY=20.0, H=H, res_frac=0.02)
        for var in ("tt", "tdl"):
            st = _plastic_state_1d(p1, var, tau4, dt4)
            d = np.array([4.0e-5])
            sig, _, _, info = j1d_step(st, d, p1, tau4, dt4, var)
            blend = (1.0 - b4) * info["sig_tr"] + b4 * info["sig_bar"]
            pv4[f"{var}_1d_H{H:+.0f}"] = abs(float(sig[0] - blend[0])) / max(abs(float(blend[0])), 1.0)
            stj = _plastic_state_j2(pj, var, tau4, dt4)
            dj = np.array([6.0e-5, -1.0e-5, -1.0e-5, 3.0e-5, 0.0, 0.0])
            sj, _, _, ij = j2_step(stj, dj, pj, tau4, dt4, var)
            bj = (1.0 - b4) * ij["sig_tr"] + b4 * ij["sig_bar"]
            pv4[f"{var}_j2_H{H:+.0f}"] = _nrm(sj - bj) / max(_nrm(bj), 1.0)
    res["PV4"] = pv4
    res["PV4_is_blend"] = bool(max(pv4.values()) < 1.0e-12)

    # ---- PV5: the blended tangent matches a CENTRAL finite difference -------
    pv5 = {}
    tau5, dt5 = 0.4, 1.0
    for H in (2000.0, -2000.0):
        p1 = j1d_params(E=20000.0, sigY=20.0, H=H, res_frac=0.02)
        pj = j2_params(E=20000.0, nu=0.2, sigY=20.0, H=H, res_frac=0.02)
        for var in ("tt", "tdl"):
            st = _plastic_state_1d(p1, var, tau5, dt5)
            d0 = np.array([4.0e-5])
            _, D, _, info = j1d_step(st, d0, p1, tau5, dt5, var)
            assert bool(info["plastic"][0]), "PV5 probe must be plastic (non-tautology)"
            hfd = 1.0e-9
            sp, _, _, _ = j1d_step(st, d0 + hfd, p1, tau5, dt5, var)
            sm, _, _, _ = j1d_step(st, d0 - hfd, p1, tau5, dt5, var)
            Dfd = float((sp[0] - sm[0]) / (2.0 * hfd))
            pv5[f"{var}_1d_H{H:+.0f}"] = abs(float(D[0]) - Dfd) / abs(Dfd)

            stj = _plastic_state_j2(pj, var, tau5, dt5)
            dj = np.array([6.0e-5, -1.0e-5, -1.0e-5, 3.0e-5, 0.0, 0.0])
            _, Dj, _, ij = j2_step(stj, dj, pj, tau5, dt5, var)
            assert ij["plastic"], "PV5 J2 probe must be plastic (non-tautology)"
            Dnum = np.zeros((6, 6))
            for j in range(6):
                e = np.zeros(6)
                e[j] = 1.0e-9
                sp, _, _, _ = j2_step(stj, dj + e, pj, tau5, dt5, var)
                sm, _, _, _ = j2_step(stj, dj - e, pj, tau5, dt5, var)
                Dnum[:, j] = (sp - sm) / (2.0e-9)
            pv5[f"{var}_j2_H{H:+.0f}"] = _nrm(Dj - Dnum) / _nrm(Dnum)
    res["PV5"] = pv5
    res["PV5_ok"] = bool(max(pv5.values()) < 1.0e-5)

    # ---- PV6: overstress above the rate-independent backbone is monotone in tau --
    # The committed state is built INVISCIDLY so the backbone is tau-independent for both
    # variants -- otherwise TDL's own relaxed alpha would move the backbone and the comparison
    # would not isolate the viscosity.
    pv6 = {}
    taus = (0.1, 0.5, 2.0, 10.0)
    for H in (0.0, 2000.0, -2000.0):
        p1 = j1d_params(E=20000.0, sigY=20.0, H=H, res_frac=0.02)
        st = _plastic_state_1d(p1, "inviscid", 0.0, 1.0)
        d = np.array([4.0e-5])
        s_inv, _, _, _ = j1d_step(st, d, p1, 0.0, 1.0, "inviscid")
        for var in ("tt", "tdl"):
            over = []
            for tv in taus:
                s, _, _, inf = j1d_step(st, d, p1, tv, 1.0, var)
                assert bool(inf["plastic"][0]), "PV6 probe must be plastic (non-tautology)"
                over.append(abs(float(s[0] - s_inv[0])))
            pv6[f"{var}_H{H:+.0f}"] = over
    res["PV6"] = pv6
    res["PV6_monotone"] = bool(all(o[0] > 0.0 and all(o[i + 1] > o[i] for i in range(len(o) - 1))
                                   for o in pv6.values()))

    res["PASS"] = bool(res["PV1_byte_exact"] and res["PV2_converges"] and res["PV3_exact"]
                       and res["PV4_is_blend"] and res["PV5_ok"] and res["PV6_monotone"])
    if verbose:
        print(f"  PV1 tau=0 BYTE-identical (both variants, 1-D + J2, H<0/=0/>0): "
              f"max|dsig| = {res['PV1_max']:.1e} ({res['PV1_byte_exact']})")
        print(f"  PV2 mid-transient -> continuous closed form: "
              + ", ".join(f"{k}: order {v['order']:.2f}" for k, v in pv2.items())
              + f" ({res['PV2_converges']})")
        print(f"  PV3 [HEADLINE] steady overstress = E*eps_rate*tau = {over_exact:.4f}: "
              f"max rel err {max(pv3.values()):.2e} ({res['PV3_exact']})")
        print(f"  PV4 update == (1-b)trial + b*projection: max rel {max(pv4.values()):.2e} "
              f"({res['PV4_is_blend']})")
        print(f"  PV5 blended tangent == central FD: max rel {max(pv5.values()):.2e} ({res['PV5_ok']})")
        print(f"  PV6 overstress monotone in tau: {res['PV6_monotone']}")
        print(f"  => PV GATE {'PASS' if res['PASS'] else 'FAIL'}")
    return res


# ===========================================================================
# 5. TT vs TDL on a hardening / softening law (the reason this WP exists)
# ===========================================================================
def _ramp_point_1d(variant, p, eps_max, T, nsteps, tau):
    """Monotonic constant-rate 1-D ramp at a single point.  Returns the path + integrals."""
    dt = T / nsteps
    deps = np.array([eps_max / nsteps])
    st = j1d_state(1)
    sig = np.zeros(nsteps + 1)
    epsp = np.zeros(nsteps + 1)
    W = 0.0
    for k in range(nsteps):
        s_old = float(st["sig"][0])
        ep_old = float(st["eps"][0] - st["sig"][0] / p["E"])
        s, _, st, _ = j1d_step(st, deps, p, tau, dt, variant)
        ep_new = float(st["eps"][0] - st["sig"][0] / p["E"])
        W += 0.5 * (s_old + float(s[0])) * (ep_new - ep_old)
        sig[k + 1] = float(s[0])
        epsp[k + 1] = ep_new
    eps = np.linspace(0.0, eps_max, nsteps + 1)
    return {"eps": eps, "sig": sig, "epsp": epsp, "work": W, "sig_end": sig[-1],
            "epsp_end": epsp[-1]}


def _cycle_point_1d(variant, p, eps_max, T, nsteps, tau, unload_frac=0.6, reload_to=1.2):
    """Load / unload / reload path -- breaks the continuous-plastic-flow assumption behind the
    TT == TDL identity for the linear law."""
    dt = T / nsteps
    n1 = int(nsteps * 0.5)
    n2 = int(nsteps * 0.2)
    n3 = nsteps - n1 - n2
    d1 = np.array([eps_max / n1])
    d2 = np.array([-(1.0 - unload_frac) * eps_max / n2])
    d3 = np.array([(reload_to - unload_frac) * eps_max / n3])
    st = j1d_state(1)
    sig = [0.0]
    for d, n in ((d1, n1), (d2, n2), (d3, n3)):
        for _ in range(n):
            s, _, st, _ = j1d_step(st, d, p, tau, dt, variant)
            sig.append(float(s[0]))
    return np.array(sig)


def run_tt_vs_tdl_point(E=20000.0, sigY=20.0, eps_max=6.0e-3, T=1.0, nsteps=2000,
                        H_over_E=(0.10, 0.02, 0.0, -0.02, -0.10),
                        De_list=(0.0, 0.003, 0.01, 0.03, 0.1, 0.3),
                        res_frac=0.02, verbose=True):
    """Quantified TT-vs-TDL study on the 1-D law across H/E and De, plus a J2-hardening check.

    Reports, per (H/E, De): end stress, peak stress, overstress above the inviscid backbone,
    dissipated plastic work, end plastic strain -- for TT and TDL, and their relative difference.

    HEADLINE FINDING -- PROVED, then measured to machine precision.

    THEOREM (1-D, monotonic).  For the 1-D associated model with ANY hardening function K(alpha),
    started from rest and loaded MONOTONICALLY, the two-track wrapper and true Duvaut-Lions
    produce the IDENTICAL stress path.
      Proof.  In TDL the update gives  sig_{n+1} = sig_tr - b*E*(abar - a_n)  and
      a_{n+1} = a_n + b*(abar - a_n), hence  sig_{n+1} + E a_{n+1} = sig_tr + E a_n, so the
      quantity  psi_n = sig_n + E a_n  advances by exactly E*de each step (elastic steps too,
      where abar = a_n).  From rest psi_0 = 0, therefore
                          a_n = eps_n - sig_n / E                                    (*)
      -- TDL's relaxed internal variable is exactly the plastic strain of its OWN stress.  Now
      substitute (*) into the projection equation sig_tr - E(abar - a_n) = K(abar):
                sig_tr - E*abar + E*(eps_{n+1} - sig_tr/E) = K(abar)
            =>  E*eps_{n+1} - E*abar = K(abar)   =>   sigbar = K(abar), abar = eps_{n+1} - sigbar/E
      which is precisely the definition of the INVISCID stress at the total strain eps_{n+1}.
      So TDL's projection target IS `inner.getStress()` on the total strain path -- the very
      quantity the two-track wrapper blends toward.  Both then return
      (1-b) sig_tr + b sig_inviscid(eps_{n+1}).  QED.
    Consequences, all measured below:
      * linear AND exponential (nonlinear) hardening/softening: TT == TDL to ~1e-13 (round-off).
      * J2 under a PROPORTIONAL (radial) path: same, because alpha is then a function of the
        plastic strain tensor alone.
      * The identity BREAKS where (*) stops reproducing the inviscid history: after UNLOADING
        (the inner material unloads elastically at once while TDL keeps relaxing) and on
        NON-PROPORTIONAL multiaxial paths.  Both are measured; the differences are 4.6e-2 ..
        3.3e-1 relative, i.e. large.
    This is the answer to the brief's OQ3 / decision D2, and it is sharper than "two-track is an
    approximation": it says exactly WHERE the approximation is exact and where it is not.
    """
    rows = []
    for hoe in H_over_E:
        p = j1d_params(E=E, sigY=sigY, H=hoe * E, res_frac=res_frac)
        inv = _ramp_point_1d("inviscid", p, eps_max, T, nsteps, 0.0)
        for De in De_list:
            tau = De * T
            tt = _ramp_point_1d("tt", p, eps_max, T, nsteps, tau)
            td = _ramp_point_1d("tdl", p, eps_max, T, nsteps, tau)
            den_s = max(abs(inv["sig_end"]), 1.0e-12)
            den_w = max(abs(inv["work"]), 1.0e-12)
            den_p = max(abs(inv["epsp_end"]), 1.0e-12)
            rows.append({
                "H_over_E": hoe, "De": De, "tau": tau,
                "sig_inv": inv["sig_end"], "sig_tt": tt["sig_end"], "sig_tdl": td["sig_end"],
                "over_tt": tt["sig_end"] - inv["sig_end"],
                "over_tdl": td["sig_end"] - inv["sig_end"],
                "W_inv": inv["work"], "W_tt": tt["work"], "W_tdl": td["work"],
                "epsp_inv": inv["epsp_end"], "epsp_tt": tt["epsp_end"], "epsp_tdl": td["epsp_end"],
                "d_sig_rel": (tt["sig_end"] - td["sig_end"]) / den_s,
                "d_W_rel": (tt["work"] - td["work"]) / den_w,
                "d_epsp_rel": (tt["epsp_end"] - td["epsp_end"]) / den_p,
                "d_path_rel": float(np.max(np.abs(tt["sig"] - td["sig"]))) / den_s,
            })
    # J2 cross-check.  "prop" = a single proportional (radial) strain direction -- the theorem
    # applies, so TT == TDL.  "nonprop" = the SAME total strain reached in two legs with a
    # 90-degree direction change (axial, then shear): alpha is no longer a function of the
    # plastic strain tensor and the theorem's step (*) fails.  This is the multiaxial analogue of
    # what a state-dependent law (SANISAND: alpha-tensor, fabric z, psi-driven M^b) does on EVERY
    # path, so it is the leg that actually bears on the SANISAND question.
    j2rows = []
    d_ax = np.array([6.0e-5, -1.0e-5, -1.0e-5, 3.0e-5, 0.0, 0.0])
    d_sh = np.array([0.0, 0.0, 0.0, 0.0, 5.0e-5, 4.0e-5])
    for hoe in (0.05, -0.05):
        pj = j2_params(E=E, nu=0.2, sigY=sigY, H=hoe * E, res_frac=res_frac)
        for path, legs in (("prop", ((d_ax, 200),)),
                           ("nonprop", ((d_ax, 120), (d_sh, 80)))):
            for De in (0.0, 0.01, 0.1):
                tau = De * T
                dt = T / 200.0
                outs = {}
                for var in ("inviscid", "tt", "tdl"):
                    st = j2_state()
                    for d, n in legs:
                        for _ in range(n):
                            s, _, st, _ = j2_step(st, d, pj, tau, dt, var)
                    outs[var] = s.copy()
                j2rows.append({"H_over_E": hoe, "path": path, "De": De,
                               "sig_inv_norm": _nrm(outs["inviscid"]),
                               "sig_tt_norm": _nrm(outs["tt"]),
                               "sig_tdl_norm": _nrm(outs["tdl"]),
                               "d_rel": _nrm(outs["tt"] - outs["tdl"])
                                        / max(_nrm(outs["tdl"]), 1e-12)})
    # perfect plasticity must give TT == TDL exactly
    p0 = j1d_params(E=E, sigY=sigY, H=0.0, res_frac=0.0)
    tt0 = _ramp_point_1d("tt", p0, eps_max, T, nsteps, 0.05 * T)
    td0 = _ramp_point_1d("tdl", p0, eps_max, T, nsteps, 0.05 * T)
    perfect_identical = float(np.max(np.abs(tt0["sig"] - td0["sig"])))

    # ---- the identity BREAKERS: a NONLINEAR (exponential) law, and a cyclic path ----
    # The exponential legs must actually EXERCISE the nonlinearity: alpha_end/alpha_f >> 1.
    # With eps_max_exp = 2e-2 and alpha_f = 4e-3 the run reaches alpha ~ 4.75 alpha_f, i.e. the
    # law traverses ~99 % of its range.  (A first pass with alpha_f = 0.5 reached
    # alpha/alpha_f ~ 0.01, where exp() IS linear -- and reproduced the linear-law identity to
    # 1e-14.  Recorded because it is exactly the tautological-fixture trap the fork's -eta gate
    # already stepped in once.)
    exp_rows = []
    eps_max_exp = 2.0e-2
    for Kinf_frac, af, tag in ((0.02, 4.0e-3, "exp-soft"), (2.0, 4.0e-3, "exp-hard")):
        pe = j1d_params(E=E, sigY=sigY, law="exp", K_inf=Kinf_frac * sigY, alpha_f=af)
        inve = _ramp_point_1d("inviscid", pe, eps_max_exp, T, nsteps, 0.0)
        for De in De_list:
            tt = _ramp_point_1d("tt", pe, eps_max_exp, T, nsteps, De * T)
            td = _ramp_point_1d("tdl", pe, eps_max_exp, T, nsteps, De * T)
            den = max(abs(inve["sig_end"]), 1.0e-12)
            exp_rows.append({
                "law": tag, "K_inf": Kinf_frac * sigY, "alpha_f": af, "De": De,
                "sig_inv": inve["sig_end"], "sig_tt": tt["sig_end"], "sig_tdl": td["sig_end"],
                "W_inv": inve["work"], "W_tt": tt["work"], "W_tdl": td["work"],
                "d_sig_rel": (tt["sig_end"] - td["sig_end"]) / den,
                "d_W_rel": (tt["work"] - td["work"]) / max(abs(inve["work"]), 1e-12),
                "d_epsp_rel": (tt["epsp_end"] - td["epsp_end"]) / max(abs(inve["epsp_end"]), 1e-12),
                "d_path_rel": float(np.max(np.abs(tt["sig"] - td["sig"]))) / den,
            })
    cyc_rows = []
    for hoe in (0.10, -0.02):
        pc = j1d_params(E=E, sigY=sigY, H=hoe * E, res_frac=res_frac)
        for De in (0.01, 0.1):
            a = _cycle_point_1d("tt", pc, eps_max, T, nsteps, De * T)
            b = _cycle_point_1d("tdl", pc, eps_max, T, nsteps, De * T)
            cyc_rows.append({"H_over_E": hoe, "De": De,
                             "d_path_abs": float(np.max(np.abs(a - b))),
                             "d_path_rel": float(np.max(np.abs(a - b))) / max(abs(a[-1]), 1e-12),
                             "sig_tt_end": float(a[-1]), "sig_tdl_end": float(b[-1])})

    out = {"rows": rows, "j2_rows": j2rows, "exp_rows": exp_rows, "cyclic_rows": cyc_rows,
           "perfect_plastic_max_diff": perfect_identical,
           "perfect_plastic_identical": bool(perfect_identical <= 1.0e-12),
           "linear_max_d_path_rel": max(r["d_path_rel"] for r in rows),
           "exp_max_d_path_rel": max(r["d_path_rel"] for r in exp_rows),
           "cyclic_max_d_path_rel": max(r["d_path_rel"] for r in cyc_rows),
           "j2_prop_max_d_rel": max(r["d_rel"] for r in j2rows if r["path"] == "prop"),
           "j2_nonprop_max_d_rel": max(r["d_rel"] for r in j2rows if r["path"] == "nonprop")}
    if verbose:
        print(f"  TT == TDL for PERFECT plasticity: max|dsig| = {perfect_identical:.1e}")
        print(f"  {'H/E':>6} {'De':>6} {'sig_inv':>9} {'sig_TT':>9} {'sig_TDL':>9} "
              f"{'dW/W':>9} {'d_eps_p':>9} {'d_path':>9}")
        for r in rows:
            print(f"  {r['H_over_E']:6.2f} {r['De']:6.3f} {r['sig_inv']:9.4f} {r['sig_tt']:9.4f} "
                  f"{r['sig_tdl']:9.4f} {r['d_W_rel']:9.2e} {r['d_epsp_rel']:9.2e} "
                  f"{r['d_path_rel']:9.2e}")
        print(f"  LINEAR law  max relative TT-vs-TDL path difference: "
              f"{out['linear_max_d_path_rel']:.2e}  (the identity)")
        print(f"  {'law':>9} {'De':>6} {'sig_inv':>9} {'sig_TT':>9} {'sig_TDL':>9} "
              f"{'dsig':>9} {'dW/W':>9} {'d_path':>9}")
        for r in exp_rows:
            print(f"  {r['law']:>9} {r['De']:6.3f} {r['sig_inv']:9.4f} {r['sig_tt']:9.4f} "
                  f"{r['sig_tdl']:9.4f} {r['d_sig_rel']:9.2e} {r['d_W_rel']:9.2e} "
                  f"{r['d_path_rel']:9.2e}")
        print(f"  EXP law     max relative TT-vs-TDL path difference: "
              f"{out['exp_max_d_path_rel']:.2e}")
        print(f"  CYCLIC (linear law) max relative TT-vs-TDL path difference: "
              f"{out['cyclic_max_d_path_rel']:.2e}")
    return out


# ===========================================================================
# 6. A0 -- the 1-D quasi-static softening bar
# ===========================================================================
def band_widths(prof, h, rel_tol=1.0e-3):
    """Three band-width definitions on a piecewise-constant element profile `prof` (>= 0).

    w1  threshold count   : h * #{e : prof_e > rel_tol * max(prof)}          (de Borst's "how many
                            elements are yielding", made explicit about its threshold)
    w2  second moment     : sqrt(12 * Var), Var = sum p_e[(x_e-xbar)^2 + h^2/12] / sum p_e.
                            THRESHOLD-FREE.  The +h^2/12 is the within-element variance of the
                            piecewise-constant profile: it makes a single-element band return
                            EXACTLY h and a k-element top hat return EXACTLY k*h, so w2 is directly
                            comparable across meshes (planning brief OQ5 / risk R4).
    w3  FWHM              : h * #{e : prof_e >= 0.5 max(prof)} -- the measure of the piecewise-
                            constant profile's half-maximum super-level set.
    """
    p = np.maximum(np.asarray(prof, float), 0.0)
    tot = float(p.sum())
    if tot <= 0.0 or p.max() <= 0.0:
        return 0.0, 0.0, 0.0
    x = (np.arange(len(p)) + 0.5) * h
    xb = float((p * x).sum() / tot)
    var = float((p * ((x - xb) ** 2)).sum() / tot) + h * h / 12.0
    w2 = float(np.sqrt(12.0 * var))
    w1 = h * int(np.count_nonzero(p > rel_tol * p.max()))
    w3 = h * int(np.count_nonzero(p >= 0.5 * p.max()))
    return w1, w2, w3


def _thomas(a, b, c, d):
    """Tridiagonal solve (no pivoting).  a = sub, b = diag, c = super, d = rhs.
    Raises on a vanishing pivot so the caller can fall back to a dense solve."""
    n = len(b)
    cp = np.empty(n)
    dp = np.empty(n)
    piv = b[0]
    if abs(piv) < 1.0e-300:
        raise ZeroDivisionError("tridiagonal pivot")
    cp[0] = c[0] / piv
    dp[0] = d[0] / piv
    for i in range(1, n):
        piv = b[i] - a[i] * cp[i - 1]
        if abs(piv) < 1.0e-300:
            raise ZeroDivisionError("tridiagonal pivot")
        cp[i] = c[i] / piv if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / piv
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def a0_bar(variant="tt", N=40, tau=0.0, T=1.0, nsteps=1000, u_max=0.25,
           imperfection="one_element", L=100.0, A=1.0, E=20000.0, sigY=20.0,
           H=-50.0, res_frac=0.02, imp_frac=0.10, imp_len_frac=0.10, imp_skew=0.013,
           law="linear", kinf_frac=0.02, alpha_f=0.5,
           rel_tol_w1=1.0e-3, tol=1.0e-10, maxit=60):
    """A0: quasi-static (no inertia) 1-D bar of N linear elements, end-displacement ramp
    u(t) = u_max * t/T with uniform dt = T/nsteps, Newton with the CONSISTENT (blended) tangent.

    PARAMETER CHOICE (documented, planning brief leg A0).  The classic Bazant / de Borst
    mesh-dependent bar needs H < 0 and localization into one element at tau = 0.  It must ALSO
    converge under end-displacement control, i.e. no structural SNAP-BACK.  For a band of length
    l_b in a bar of length L the post-peak `du/dsigma` changes sign when

        l_b < L * |H| / E                                    (snap-back)

    (from du = l_b dsig/E_t + (L-l_b) dsig/E with E_t = E H/(E+H) < 0).  With L = 100 mm and the
    finest mesh N = 160 (h = 0.625 mm) this demands |H|/E < 1/160.  The default H = -50 MPa gives
    |H|/E = 1/400, a factor 2.5 of margin, so NO mesh in the A0 sweep snaps back and Newton
    converges on every leg -- including the tau = 0 negative control, which is the leg that must
    run for H1 to mean anything.  A residual floor (`res_frac`) keeps the yield stress positive.

    `imperfection`:
      "one_element"   de Borst's convention (a): ONE central element at (1-imp_frac)*sigY.  Its
                      physical size shrinks with N.
      "fixed_flat"    convention (b) as specified: a FIXED physical zone (imp_len_frac*L, i.e.
                      2/4/8/16 elements) all at (1-imp_frac)*sigY.  Every element in the zone is
                      EXACTLY as weak as its neighbours -- a tie.
      "fixed_graded"  convention (b'): the same fixed zone but with a parabolic notch
                      sigY(x) = sigY[1 - imp_frac(1 - r^2)], r = 2|x-L/2|/(imp_len_frac*L), so the
                      imperfection is a fixed, mesh-CONVERGENT field with a unique weakest point.

    Returns a dict with the load-displacement curve, dissipated work, plastic-strain profile,
    the three widths, and convergence diagnostics.
    """
    h = L / N
    xc = (np.arange(N) + 0.5) * h
    sY = np.full(N, float(sigY))
    if imperfection == "one_element":
        sY[N // 2] *= (1.0 - imp_frac)
    elif imperfection == "fixed_flat":
        sY[np.abs(xc - 0.5 * L) <= 0.5 * imp_len_frac * L] *= (1.0 - imp_frac)
    elif imperfection == "fixed_graded":
        # The notch centre is offset by `imp_skew` of the zone half-width.  Without it the
        # parabola is symmetric about x = L/2, which for an EVEN element count is a NODE, so the
        # two central elements get EXACTLY equal sigY -- an exact tie, and the tau=0 band comes
        # out 2h wide for a numerical reason rather than a physical one.  The offset is
        # mesh-independent, so the imperfection field still converges under refinement.
        xcen = 0.5 * L + imp_skew * 0.5 * imp_len_frac * L
        r = np.abs(xc - xcen) / (0.5 * imp_len_frac * L)
        sY = np.where(r < 1.0, sigY * (1.0 - imp_frac * (1.0 - r * r)), sY)
    else:
        raise ValueError(f"unknown imperfection {imperfection!r}")

    if law == "linear":
        p = j1d_params(E=E, sigY=sY, H=H, res_frac=res_frac, H_res=E / 2.0e5)
    else:
        p = j1d_params(E=E, sigY=sY, law="exp", K_inf=kinf_frac * sY, alpha_f=alpha_f)
    dt = T / nsteps
    u = np.zeros(N + 1)
    st = j1d_state(N)
    eps_c = np.zeros(N)
    sig_c = np.zeros(N)
    epsp_c = np.zeros(N)

    U = np.zeros(nsteps + 1)
    P = np.zeros(nsteps + 1)
    Wcum = np.zeros(nsteps + 1)
    EPSP = np.zeros((nsteps + 1, N))
    iters = []
    converged = True
    ref = tol * A * float(np.max(sY))

    def _solve(D, rhs):
        kd = A / h * (D[:-1] + D[1:])
        ka = np.empty(N - 1)
        kc = np.empty(N - 1)
        ka[0] = 0.0
        ka[1:] = -A / h * D[1:-1]
        kc[:-1] = -A / h * D[1:-1]
        kc[-1] = 0.0
        try:
            return _thomas(ka, kd, kc, rhs)
        except (ZeroDivisionError, FloatingPointError):
            Kd = np.diag(kd) + np.diag(kc[:-1], 1) + np.diag(ka[1:], -1)
            return np.linalg.solve(Kd, rhs)

    def _resid(uu):
        eps = (uu[1:] - uu[:-1]) / h
        sig, D, _, _ = j1d_step(st, eps - eps_c, p, tau, dt, variant)
        R = np.zeros(N + 1)
        R[:-1] += A * sig
        R[1:] -= A * sig                          # R = F_ext - F_int, F_ext = 0
        return R[1:N], sig, D

    D_c = np.full(N, E)                           # tangent at the last converged state
    for step in range(nsteps):
        # PREDICTOR: apply the prescribed end increment through the CONSISTENT tangent of the
        # last converged state.  Without it the first iterate of every step carries an O(E*du/h)
        # artificial residual concentrated at the loaded end; once the band softens, linearizing
        # that residual on a tangent 400x softer than E throws the iterate clean out of the
        # quadratic basin and Newton diverges at the very first post-peak step (measured: N>=40,
        # max|du| 42x the step increment).  This is the same "use the last converged tangent, not
        # the elastic one" rule the fork's own algorithms rely on.
        due = u_max / nsteps
        rhs = np.zeros(N - 1)
        rhs[-1] = A / h * D_c[N - 1] * due
        u[1:N] += _solve(D_c, rhs)
        u[N] = u_max * (step + 1) / nsteps

        ok = False
        for it in range(maxit):
            Ri, sig, D = _resid(u)
            rn = float(np.max(np.abs(Ri)))
            if rn <= max(ref, 1e-14):
                ok = True
                break
            du = _solve(D, Ri)
            # backtracking safeguard: an indefinite tangent may overshoot on the yielding step
            lam, u_try = 1.0, u.copy()
            for _ in range(8):
                u_try = u.copy()
                u_try[1:N] += lam * du
                if float(np.max(np.abs(_resid(u_try)[0]))) < rn:
                    break
                lam *= 0.5
            u = u_try
        iters.append(it + 1)
        if not ok:
            converged = False
            break
        D_c = D
        eps = (u[1:] - u[:-1]) / h
        sig, D, st, _ = j1d_step(st, eps - eps_c, p, tau, dt, variant)
        epsp = eps - sig / E
        Wcum[step + 1] = Wcum[step] + float(np.sum(A * h * 0.5 * (sig_c + sig) * (epsp - epsp_c)))
        eps_c, sig_c, epsp_c = eps, sig.copy(), epsp
        U[step + 1] = u[N]
        P[step + 1] = float(A * np.mean(sig))
        EPSP[step + 1] = epsp

    ns = len(iters) if converged else max(len(iters) - 1, 0)
    U, P, Wcum, EPSP = U[:ns + 1], P[:ns + 1], Wcum[:ns + 1], EPSP[:ns + 1]
    kpk = int(np.argmax(P))
    prof_end = np.maximum(EPSP[-1], 0.0)
    prof_post = np.maximum(EPSP[-1] - EPSP[kpk], 0.0)
    w1e, w2e, w3e = band_widths(prof_end, h, rel_tol_w1)
    w1p, w2p, w3p = band_widths(prof_post, h, rel_tol_w1)

    # A MESH-FAIR energy measure.  Comparing W at a common end displacement is not fair across
    # meshes: at fixed u a coarse band (large h) has softened much less than a fine one, so W(u)
    # mixes "how much energy per unit band volume" with "how far along the softening curve we
    # are".  W_50 is the dissipated work at the instant the load first falls to 50 % of its peak
    # -- a common POINT ON THE SOFTENING CURVE -- which is what "fracture energy scales with h"
    # actually claims.  NaN when the run never drops that far (the viscous, non-localizing runs).
    W50, u50 = float("nan"), float("nan")
    thr = 0.5 * float(np.max(P))
    idx = np.nonzero((np.arange(len(P)) > kpk) & (P <= thr))[0]
    if idx.size:
        k = int(idx[0])
        f = (P[k - 1] - thr) / max(P[k - 1] - P[k], 1e-300)
        W50 = float(Wcum[k - 1] + f * (Wcum[k] - Wcum[k - 1]))
        u50 = float(U[k - 1] + f * (U[k] - U[k - 1]))
    # non-tautology diagnostic: how much of the response is actually inelastic
    W_el_ref = 0.5 * float(np.max(P)) * float(np.max(P)) / (A * E / L)
    return {
        "variant": variant, "law": law, "N": N, "h": h, "tau": tau, "T": T,
        "De": tau / T if T > 0 else 0.0,
        "nsteps": nsteps, "dt": dt, "beta": dl_beta(tau, dt), "imperfection": imperfection,
        "u": U, "P": P, "W": Wcum, "W_end": float(Wcum[-1]), "W_50": W50, "u_50": u50,
        "P_peak": float(np.max(P)), "u_peak": float(U[kpk]), "k_peak": kpk,
        "P_end": float(P[-1]), "epsp_end": prof_end, "epsp_post": prof_post,
        "w1": w1p, "w2": w2p, "w3": w3p,           # HEADLINE widths: post-peak plastic increment
        "w1_end": w1e, "w2_end": w2e, "w3_end": w3e,
        "epsp_max": float(prof_end.max()), "converged": converged,
        "iters_mean": float(np.mean(iters)) if iters else 0.0,
        "iters_max": int(np.max(iters)) if iters else 0,
        "inelastic_ratio": float(Wcum[-1]) / max(W_el_ref, 1e-30),
    }


def min_steps_to_converge(n_lo=250, n_hi=64000, **kw):
    """Smallest uniform step count (doubling search) at which the A0 bar completes.  This is a
    REGULARIZATION SIGNATURE in its own right: the ill-posed tau = 0 problem needs far more steps
    than the same mesh at tau > 0, and the gap widens with refinement.  Returns None if the run
    never completes within n_hi."""
    n = n_lo
    while n <= n_hi:
        if a0_bar(nsteps=n, **kw)["converged"]:
            return n
        n *= 2
    return None


def run_a0_sweep(variants=("tt", "tdl"), Ns=(20, 40, 80, 160),
                 De_list=(0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0),
                 imperfections=("one_element", "fixed_flat", "fixed_graded"),
                 T=1.0, nsteps=1000, u_max=0.25, verbose=True, **kw):
    """The full A0 sweep.  Returns a flat list of `a0_bar` result dicts (curves stripped to keep
    the object small; `u`/`P`/`W` are kept because the H2 curve-convergence check needs them)."""
    out = []
    for imp in imperfections:
        for var in variants:
            for De in De_list:
                for N in Ns:
                    r = a0_bar(variant=var, N=N, tau=De * T, T=T, nsteps=nsteps,
                               u_max=u_max, imperfection=imp, **kw)
                    r.pop("epsp_end", None)
                    out.append(r)
                    if verbose:
                        print(f"  {imp:13s} {var:3s} De={De:<6.3f} N={N:<4d} "
                              f"w2={r['w2']:8.4f} w1={r['w1']:8.4f} w3={r['w3']:8.4f} "
                              f"Ppk={r['P_peak']:7.3f} W={r['W_end']:9.4f} "
                              f"conv={r['converged']} it={r['iters_max']}")
    return out


def curve_distance(rA, rB):
    """L2 distance between two load-displacement curves, normalized by the peak load.
    Both runs use the same u grid only if nsteps matches; otherwise interpolate onto A's grid."""
    uA, PA = rA["u"], rA["P"]
    uB, PB = rB["u"], rB["P"]
    Pi = np.interp(uA, uB, PB)
    den = max(rA["P_peak"], 1.0e-30) * np.sqrt(len(uA))
    return float(np.linalg.norm(PA - Pi) / den)


# ===========================================================================
# 7. P0b leg (a) -- a STATE-DEPENDENT elastic operator breaks the identity on EVERY path
# ===========================================================================
def run_pressure_dependent_leg(E=20000.0, sigY=20.0, H=2000.0, eps_max=6.0e-3, T=1.0,
                               nsteps=2000, E_beta=0.6,
                               De_list=(1.0e-4, 3.0e-4, 1.0e-3, 3.0e-3,
                                        0.01, 0.03, 0.1, 0.3),
                               verbose=True):
    """P0b-(a).  Repeat the TT-vs-TDL identity on a HYPOELASTIC 1-D law,
    `E(sigma) = E0 (1 + E_beta |sigma| / sigY)`, under monotonic proportional ISOTROPIC hardening.

    Why this kills the ADR-90 theorem.  The proof rests on `psi = sigma + E alpha` advancing by
    exactly `E de` every step, which requires E to be the SAME NUMBER on both tracks.  A generic
    wrapper's only elastic operator is `inner.getInitialTangent()` -- E evaluated at the INNER's
    committed stress `sigma_bar_n`.  But the stress it must predict from is its OWN `sigma_vp,n`,
    and under relaxation `sigma_vp != sigma_bar`.  So the wrapper integrates
    `sigma_tr = sigma_vp,n + E(sigma_bar_n) de` while true DL integrates
    `sigma_vp,n + E(sigma_vp,n) de`: the conserved quantity is destroyed at first yield and the
    error accumulates on EVERY path, monotonic and proportional included.

    Three columns are reported:
      tt      the generic wrapper (predictor modulus at the inner's stress)   -- what is buildable
      tt_vp   the same blend with the modulus at the WRAPPER's own stress     -- NOT buildable
      tdl     true Duvaut-Lions
    `tt_vp` is the obvious repair and it is measured here only to establish whether ANY
    wrapper-accessible choice of a single scalar modulus closes the gap.
    """
    rows = []
    for beta_E, tag in ((0.0, "constant E"), (E_beta, f"E(sig), beta_E={E_beta}")):
        p = j1d_params(E=E, sigY=sigY, H=H, E_beta=beta_E)
        inv = _ramp_point_1d("inviscid", p, eps_max, T, nsteps, 0.0)
        for De in De_list:
            tau = De * T
            tt = _ramp_point_1d("tt", p, eps_max, T, nsteps, tau)
            tv = _ramp_point_1d("tt_vp", p, eps_max, T, nsteps, tau)
            td = _ramp_point_1d("tdl", p, eps_max, T, nsteps, tau)
            den = max(abs(td["sig_end"]), 1.0e-12)
            dpath = float(np.max(np.abs(tt["sig"] - td["sig"]))) / den
            dpath_v = float(np.max(np.abs(tv["sig"] - td["sig"]))) / den
            rows.append({
                "case": tag, "E_beta": beta_E, "De": De,
                "sig_inv": inv["sig_end"], "sig_tt": tt["sig_end"],
                "sig_ttvp": tv["sig_end"], "sig_tdl": td["sig_end"],
                "d_end_rel": abs(tt["sig_end"] - td["sig_end"]) / den,
                "d_path_rel": dpath, "d_path_rel_ttvp": dpath_v,
                "W_tt": tt["work"], "W_tdl": td["work"],
                "d_W_rel": abs(tt["work"] - td["work"]) / max(abs(td["work"]), 1e-12),
                "repair_ratio": dpath_v / dpath if dpath > 0 else float("nan"),
            })
    out = {"rows": rows,
           "const_max": max(r["d_path_rel"] for r in rows if r["E_beta"] == 0.0),
           "var_max": max(r["d_path_rel"] for r in rows if r["E_beta"] != 0.0),
           "var_max_ttvp": max(r["d_path_rel_ttvp"] for r in rows if r["E_beta"] != 0.0)}
    # the working window of the A0 bar deck
    win = [r for r in rows if r["E_beta"] != 0.0 and 1.0e-4 <= r["De"] <= 1.0e-2]
    out["window_max"] = max(r["d_path_rel"] for r in win)
    out["window_max_ttvp"] = max(r["d_path_rel_ttvp"] for r in win)
    if verbose:
        print(f"  {'case':<22}{'De':>8}{'sig_inv':>10}{'sig_TT':>10}{'sig_TTvp':>10}"
              f"{'sig_TDL':>10}{'d_path':>11}{'d_path_vp':>11}{'dW':>10}")
        for r in rows:
            print(f"  {r['case']:<22}{r['De']:>8.4g}{r['sig_inv']:>10.4f}{r['sig_tt']:>10.4f}"
                  f"{r['sig_ttvp']:>10.4f}{r['sig_tdl']:>10.4f}{r['d_path_rel']:>11.2e}"
                  f"{r['d_path_rel_ttvp']:>11.2e}{r['d_W_rel']:>10.2e}")
        print(f"  constant E   max rel TT-vs-TDL = {out['const_max']:.2e}   (the ADR-90 theorem)")
        print(f"  E(sigma)     max rel TT-vs-TDL = {out['var_max']:.2e}"
              f"   over the working window 1e-4..1e-2: {out['window_max']:.2e}")
        print(f"  E(sigma)     the NON-buildable repair (modulus at the wrapper's own stress): "
              f"{out['var_max_ttvp']:.2e}")
    return out


# ===========================================================================
# 8. P0b leg (b) -- the DISSIPATION gate
# ===========================================================================
def run_dissipation_gate(E=20000.0, sigY=20.0, nsteps=1200, T=1.0,
                         De_list=(1.0e-3, 1.0e-2, 0.1), verbose=True):
    """P0b-(b).  Accumulate the wrapper's own dissipation increment

        D_w = sigma_vp . dEps^vp_w ,   dEps^vp_w = C_e^-1 (sigma_vp - sigma_bar) dt/tau
                                                 = beta C_e^-1 (sigma_tr - sigma_bar)

    on three paths where the ADR-90 theorem does NOT apply.  In true Duvaut-Lions `sigma_bar` is
    the closest-point projection OF `sigma_vp`, so `sigma_vp - sigma_bar` is an outward normal
    cone direction and D_w >= 0 follows from convexity.  In the two-track wrapper `sigma_bar` is
    the projection of a DIFFERENT trial (the inner's), so nothing forces that geometry and the
    sign is an open question -- which is why it is measured rather than asserted.

    The inner material's own dissipation `sigma_bar . dEps^p_inner` is reported alongside, so the
    sign of the composite is known and a negative wrapper increment cannot hide behind a positive
    total.
    """
    res = {}
    dt = T / nsteps

    # -- path 1: J2, non-proportional (axial leg, then a 90-degree turn into shear) -------
    d_ax = np.array([6.0e-5, -1.0e-5, -1.0e-5, 3.0e-5, 0.0, 0.0])
    d_sh = np.array([0.0, 0.0, 0.0, 0.0, 5.0e-5, 4.0e-5])
    # -- path 3: swipe surrogate -- rotate the strain DIRECTION at fixed |eps| ------------
    for tag, legs in (("j2_nonprop", ((d_ax, 120), (d_sh, 80))),
                      ("j2_swipe", None)):
        for H in (2000.0, -1000.0):
            pj = j2_params(E=E, nu=0.2, sigY=sigY, H=H, res_frac=0.02)
            for De in De_list:
                tau = De * T
                st = j2_state()
                Dw, Dw_min, Din, n_neg, neg_mag, neg_sum = 0.0, float("inf"), 0.0, 0, 0.0, 0.0
                if tag == "j2_swipe":
                    # hold |eps| and rotate the deviatoric direction through 90 degrees
                    amp, nrot = 4.0e-3, 200
                    prev = np.zeros(6)
                    steps = []
                    for k in range(nrot):
                        th = 0.5 * np.pi * k / (nrot - 1)
                        e = amp * np.array([np.cos(th), -np.cos(th), 0.0,
                                            np.sin(th), 0.0, 0.0])
                        steps.append(e - prev)
                        prev = e
                else:
                    steps = [d for d, n in legs for _ in range(n)]
                sig_in_prev = st["sig_in"].copy()
                for d in steps:
                    s, _, stn, inf = j2_step(st, d, pj, tau, dt, "tt")
                    Dw += inf["Dw"]
                    if inf["Dw"] < Dw_min:
                        Dw_min = inf["Dw"]
                    if inf["Dw"] < 0.0:
                        n_neg += 1
                        neg_mag = min(neg_mag, inf["Dw"])
                        neg_sum += inf["Dw"]
                    # inner's own dissipation: sigma_bar . (deps - Ce^-1 dsig_bar)
                    dep_in = d - np.linalg.solve(inf["Ce"], inf["sig_bar"] - sig_in_prev)
                    Din += float(inf["sig_bar"] @ dep_in)
                    sig_in_prev = inf["sig_bar"].copy()
                    st = stn
                res[f"{tag}_H{H:+.0f}_De{De:g}"] = {
                    "path": tag, "H": H, "De": De, "Dw_cum": Dw, "Dw_min": Dw_min,
                    "n_neg": n_neg, "neg_worst": neg_mag, "neg_sum": neg_sum,
                    "D_inner": Din, "neg_frac": n_neg / len(steps)}

    # -- path 2: 1-D cyclic load / unload / reload ----------------------------------------
    for H in (2000.0, -1000.0):
        p = j1d_params(E=E, sigY=sigY, H=H, res_frac=0.02)
        for De in De_list:
            tau = De * T
            st = j1d_state(1)
            n1, n2, n3 = nsteps // 2, nsteps // 5, nsteps - nsteps // 2 - nsteps // 5
            emax = 6.0e-3
            legs = ((np.array([emax / n1]), n1),
                    (np.array([-0.4 * emax / n2]), n2),
                    (np.array([0.6 * emax / n3]), n3))
            Dw, Dw_min, Din, n_neg, neg_mag, ns, neg_sum = 0.0, float("inf"), 0.0, 0, 0.0, 0, 0.0
            sprev = 0.0
            for d, n in legs:
                for _ in range(n):
                    s, _, stn, inf = j1d_step(st, d, p, tau, dt, "tt")
                    w = float(inf["Dw"][0])
                    Dw += w
                    Dw_min = min(Dw_min, w)
                    if w < 0.0:
                        n_neg += 1
                        neg_mag = min(neg_mag, w)
                        neg_sum += w
                    sb = float(inf["sig_bar"][0])
                    Din += sb * (float(d[0]) - (sb - sprev) / float(np.atleast_1d(inf["Ce"])[0]))
                    sprev = sb
                    st = stn
                    ns += 1
            res[f"cyclic1d_H{H:+.0f}_De{De:g}"] = {
                "path": "cyclic1d", "H": H, "De": De, "Dw_cum": Dw, "Dw_min": Dw_min,
                "n_neg": n_neg, "neg_worst": neg_mag, "neg_sum": neg_sum,
                "D_inner": Din, "neg_frac": n_neg / ns}

    # Separate ROUND-OFF negatives from real ones: a negative increment counts only if its
    # magnitude exceeds 1e-9 of the run's largest increment.  Without this the gate would report
    # -1.5e-18 on a path whose cumulative is +0.19 and look like a violation when it is noise.
    for v in res.values():
        scale = max(abs(v["Dw_cum"]), 1e-300)
        v["neg_sum_rel"] = v["neg_sum"] / scale
        v["worst_rel"] = v["neg_worst"] / scale
        v["significant"] = bool(abs(v["neg_worst"]) > 1.0e-9 * scale)
    sig_rows = {k: v for k, v in res.items() if v["significant"]}
    out = {"rows": res,
           "worst_step": min(v["Dw_min"] for v in res.values()),
           "worst_step_rel": min(v["worst_rel"] for v in res.values()),
           "any_negative_significant": bool(sig_rows),
           "neg_sum_rel_worst": min(v["neg_sum_rel"] for v in res.values()),
           "paths_with_significant_negatives": sorted(sig_rows),
           "paths_with_roundoff_only": sorted(k for k, v in res.items()
                                              if v["n_neg"] > 0 and not v["significant"])}
    if verbose:
        print(f"  {'case':<28}{'Dw_cum':>12}{'Dw_min/step':>14}{'n_neg':>7}{'neg_frac':>10}"
              f"{'sum(neg)/cum':>14}{'signif':>8}{'D_inner':>12}")
        for k, v in res.items():
            print(f"  {k:<28}{v['Dw_cum']:>12.4e}{v['Dw_min']:>14.3e}{v['n_neg']:>7d}"
                  f"{v['neg_frac']:>10.3f}{v['neg_sum_rel']:>14.2e}"
                  f"{str(v['significant']):>8}{v['D_inner']:>12.4e}")
        print(f"  worst single-step D_w = {out['worst_step']:.3e} "
              f"({out['worst_step_rel']:.2e} of that run's cumulative)")
        print(f"  significant (non-round-off) violations on: "
              f"{out['paths_with_significant_negatives']}")
    return out


# ===========================================================================
# 9. P0b leg (c) -- is the converged width a TAU property or an IMPERFECTION property?
# ===========================================================================
def run_imperfection_study(De=3.0e-4, T=1.0, nsteps=2000, u_max=2.0,
                           Ns=(80, 160), amps=(0.05, 0.10, 0.20),
                           lens=(0.05, 0.10, 0.20), verbose=True):
    """P0b-(c).  At the ONE De where A0 saw both w2 and w3 converge, vary the imperfection
    amplitude and the imperfection zone length.  If the converged width tracks the imperfection,
    the width is an imperfection property that viscosity merely stops from collapsing -- not a
    material/numerical length set by tau.

    `w2/h` is reported per run: it is the band-resolution floor, and w2 is bounded above by L by
    construction (a uniform profile has w2 = L exactly), so a w2 near L means "no band", not
    "a wide band".
    """
    rows = []
    for amp in amps:
        for ln in lens:
            for N in Ns:
                r = a0_bar(variant="tt", N=N, tau=De * T, T=T, nsteps=nsteps, u_max=u_max,
                           imperfection="fixed_graded", imp_frac=amp, imp_len_frac=ln)
                rows.append({"amp": amp, "len_frac": ln, "N": N, "h": r["h"],
                             "converged": r["converged"], "w2": r["w2"], "w3": r["w3"],
                             "w2_over_h": r["w2"] / r["h"], "w2_over_L": r["w2"] / 100.0,
                             "P_peak": r["P_peak"], "W_end": r["W_end"],
                             "epsp_max": r["epsp_max"]})
    fine = [r for r in rows if r["N"] == max(Ns) and r["converged"]]
    at10 = [r for r in fine if r["len_frac"] == 0.10]
    amp10 = [r for r in fine if r["amp"] == 0.10]
    out = {"rows": rows,
           "w2_span_over_amp": (max(r["w2"] for r in at10) / min(r["w2"] for r in at10)) if at10 else float("nan"),
           "w2_span_over_len": (max(r["w2"] for r in amp10) / min(r["w2"] for r in amp10)) if amp10 else float("nan"),
           "w3_span_over_len": (max(r["w3"] for r in amp10) / min(r["w3"] for r in amp10)) if amp10 else float("nan"),
           "min_w2_over_h": min(r["w2_over_h"] for r in rows if r["converged"])}
    if verbose:
        print(f"  {'amp':>6}{'zone/L':>8}{'N':>5}{'h':>8}{'w2':>10}{'w3':>9}{'w2/h':>8}"
              f"{'w2/L':>8}{'P_peak':>10}{'conv':>6}")
        for r in rows:
            print(f"  {r['amp']:>6.2f}{r['len_frac']:>8.2f}{r['N']:>5}{r['h']:>8.4f}"
                  f"{r['w2']:>10.4f}{r['w3']:>9.4f}{r['w2_over_h']:>8.2f}{r['w2_over_L']:>8.3f}"
                  f"{r['P_peak']:>10.4f}{str(r['converged']):>6}")
        print(f"  w2 spread over IMPERFECTION AMPLITUDE (zone 10 %, finest mesh): "
              f"x{out['w2_span_over_amp']:.2f}")
        print(f"  w2 spread over ZONE LENGTH        (amp 10 %, finest mesh): "
              f"x{out['w2_span_over_len']:.2f}")
    return out


# ===========================================================================
# 10. P0b leg (d) -- plane-strain NON-ASSOCIATED Drucker-Prager + the blended acoustic tensor
# ===========================================================================
def dp_params(E=20000.0, nu=0.25, rho=0.40, rho_bar=0.15, c0=20.0, Hc=-200.0, c_res=2.0):
    """Non-associated Drucker-Prager, tension-positive, Voigt-6 fork ordering.

        f = sqrt(J2) + rho * p - c(alpha)          p = I1/3 (compression negative)
        g = sqrt(J2) + rho_bar * p                 rho_bar < rho  =>  NON-associated
        c(alpha) = max(c0 + Hc*alpha, c_res)       Hc < 0  =>  cohesion softening

    Well-posedness of the return map needs `G + rho*rho_bar*K + Hc > 0`.
    """
    G = E / (2.0 * (1.0 + nu))
    K = E / (3.0 * (1.0 - 2.0 * nu))
    if G + rho * rho_bar * K + Hc <= 0.0:
        raise ValueError("DP return map ill-posed: G + rho*rho_bar*K + Hc <= 0")
    return {"E": E, "nu": nu, "G": G, "K": K, "rho": rho, "rho_bar": rho_bar,
            "c0": c0, "Hc": Hc, "c_res": c_res}


def dp_elastic_C(p):
    return _j2_C_matrix(p["K"], p["G"], 1.0, 0.0, None)


def dp_cohesion(alpha, p):
    return max(p["c0"] + p["Hc"] * alpha, p["c_res"])


def dp_return(sig_tr_v, alpha_n, p, fd=1.0e-7):
    """Closest-point return on the DP cone (smooth part) with an apex flag; the ALGORITHMIC
    tangent is obtained by central differencing the return map itself, which is exact to ~1e-8
    and cannot get the non-associated / unsymmetric structure wrong by hand.

    Returns (sig_bar_v, alpha_bar, C_alg 6x6, plastic, apex)."""
    def _core(str_v, a_n):
        S = _sig_v2t(np.asarray(str_v, float))
        ph = np.trace(S) / 3.0
        s = S - ph * np.eye(3)
        q = float(np.sqrt(0.5 * np.sum(s * s)))                 # sqrt(J2)
        cn = p["c0"] + p["Hc"] * a_n
        floored = cn < p["c_res"]
        c_of = (lambda a: p["c_res"]) if floored else (lambda a: p["c0"] + p["Hc"] * a)
        Hl = 0.0 if floored else p["Hc"]
        f = q + p["rho"] * ph - c_of(a_n)
        if f <= 0.0:
            return np.asarray(str_v, float).copy(), float(a_n), False, False
        dg = f / (p["G"] + p["rho"] * p["rho_bar"] * p["K"] + Hl)
        q_new = q - p["G"] * dg
        if q_new < 0.0:                                          # apex return
            dg = (c_of(a_n) - p["rho"] * ph) / max(p["rho"] * p["rho_bar"] * p["K"] + Hl, 1e-12)
            dg = max(dg, 0.0)
            p_new = ph - p["K"] * p["rho_bar"] * dg
            return _sig_t2v(p_new * np.eye(3)), float(a_n) + dg, True, True
        s_new = s * (q_new / q) if q > 0.0 else s
        p_new = ph - p["K"] * p["rho_bar"] * dg
        return _sig_t2v(s_new + p_new * np.eye(3)), float(a_n) + dg, True, False

    sig, alpha, plastic, apex = _core(sig_tr_v, alpha_n)
    Ce = dp_elastic_C(p)
    if not plastic:
        return sig, alpha, Ce, False, False
    C = np.zeros((6, 6))
    scale = max(float(np.max(np.abs(sig_tr_v))), 1.0)
    for j in range(6):
        e = np.zeros(6)
        e[j] = fd * scale
        sp, _, _, _ = _core(np.asarray(sig_tr_v, float) + Ce @ e, alpha_n)
        sm, _, _, _ = _core(np.asarray(sig_tr_v, float) - Ce @ e, alpha_n)
        C[:, j] = (sp - sm) / (2.0 * fd * scale)
    return sig, alpha, C, True, apex


def dp_state():
    z = np.zeros(6)
    return {"sig": z.copy(), "eps": z.copy(), "sig_in": z.copy(), "alp_in": 0.0, "alp": 0.0}


def dp_step(state, deps, p, tau, dt, variant="tt"):
    """TT / TDL / inviscid over the non-associated DP model.  Same algebra as `j2_step`."""
    Ce = dp_elastic_C(p)
    beta = dl_beta(tau, dt)
    deps = np.asarray(deps, float)
    sig_tr = state["sig"] + Ce @ deps

    if variant == "tt":
        sig_tr_in = state["sig_in"] + Ce @ deps
        sig_bar, alp_bar, Cep, plastic, apex = dp_return(sig_tr_in, state["alp_in"], p)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        D = (1.0 - beta) * Ce + beta * Cep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": sig_bar, "alp_in": alp_bar, "alp": state["alp"]}
    elif variant == "tdl":
        sig_bar, alp_bar, Cep, plastic, apex = dp_return(sig_tr, state["alp"], p)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        D = (1.0 - beta) * Ce + beta * Cep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": state["sig_in"].copy(), "alp_in": state["alp_in"],
               "alp": (1.0 - beta) * state["alp"] + beta * alp_bar}
    elif variant == "inviscid":
        sig_bar, alp_bar, Cep, plastic, apex = dp_return(sig_tr, state["alp"], p)
        sig, D = sig_bar, Cep
        new = {"sig": sig, "eps": state["eps"] + deps,
               "sig_in": sig_bar, "alp_in": alp_bar, "alp": alp_bar}
    else:
        raise ValueError(f"unknown variant {variant!r}")

    dvp = beta * np.linalg.solve(Ce, sig_tr - sig_bar)
    info = {"beta": beta, "sig_tr": sig_tr, "sig_bar": sig_bar, "Cep": Cep, "Ce": Ce,
            "plastic": plastic, "apex": apex, "dvp": dvp, "Dw": float(sig @ dvp)}
    return sig, D, new, info


def acoustic_Q(Cv, n):
    """Q_ik = n_j C_ijkl n_l, built from the 6x6 Voigt matrix WITHOUT hand-rolling the
    engineering-shear factors: for each basis direction g, feed C the strain tensor
    sym(g (x) n) and contract the resulting stress with n."""
    n = np.asarray(n, float)
    Q = np.zeros((3, 3))
    for k in range(3):
        g = np.zeros(3)
        g[k] = 1.0
        Eps = 0.5 * (np.outer(g, n) + np.outer(n, g))
        ev = np.array([Eps[0, 0], Eps[1, 1], Eps[2, 2],
                       2.0 * Eps[0, 1], 2.0 * Eps[1, 2], 2.0 * Eps[0, 2]])
        Q[:, k] = _sig_v2t(Cv @ ev) @ n
    return Q


def min_det_acoustic(Cv, Cv_elastic=None, ntheta=721):
    """Minimum over in-plane band normals of det of the 2x2 IN-PLANE acoustic tensor, normalised
    by the elastic value.  <= 0 means loss of ellipticity (a discontinuous bifurcation is
    possible; Rudnicki & Rice 1975).  Plane strain, so the normal and the jump both lie in
    (x, y) and the relevant block is the 2x2."""
    best, best_th = np.inf, 0.0
    Ce = dp_elastic_C({"K": 1.0, "G": 1.0}) if Cv_elastic is None else Cv_elastic
    for th in np.linspace(0.0, np.pi, ntheta):
        nn = np.array([np.cos(th), np.sin(th), 0.0])
        d = float(np.linalg.det(acoustic_Q(Cv, nn)[:2, :2]))
        de = float(np.linalg.det(acoustic_Q(Ce, nn)[:2, :2])) if Cv_elastic is not None else 1.0
        v = d / de if de != 0.0 else d
        if v < best:
            best, best_th = v, th
    return best, best_th


def _dp_plane_strain_path(p, tau, dt, variant, n_iso=60, n_shear=800,
                          p_target=-100.0, dgam=5.0e-5, rotate=0.0):
    """Isotropic consolidation to `p_target`, then a plane-strain deviatoric push.  `rotate` > 0
    turns the deviatoric strain DIRECTION through `rotate` radians over the shear leg (the
    rotating-principal-direction path); `rotate = 0` is the proportional path.

    Returns the committed history: stresses, cohesion, algorithmic tangents, and the flags."""
    st = dp_state()
    e_iso = p_target / (3.0 * p["K"] * n_iso)      # p = K * eps_vol, eps_vol = 3 * e_iso
    for _ in range(n_iso):
        _, _, st, _ = dp_step(st, np.array([e_iso, e_iso, e_iso, 0, 0, 0]),
                              p, tau, dt, variant)
    hist = []
    for k in range(n_shear):
        th = rotate * k / max(n_shear - 1, 1)
        d = dgam * np.array([-np.cos(2 * th), np.cos(2 * th), 0.0,
                             np.sin(2 * th), 0.0, 0.0])
        s, D, st, inf = dp_step(st, d, p, tau, dt, variant)
        hist.append({"k": k, "sig": s.copy(), "D": D.copy(), "Cep": inf["Cep"].copy(),
                     "beta": inf["beta"], "plastic": inf["plastic"], "apex": inf["apex"],
                     "alp": st["alp_in"] if variant == "tt" else st["alp"],
                     "Dw": inf["Dw"]})
    return st, hist


def run_dp_leg(nsteps=800, T=1.0, De_list=(1.0e-4, 3.0e-4, 1.0e-3), c_res=5.0,
               n_shear=800, ntheta=361, verbose=True):
    """P0b-(d).  The only leg that speaks to WELL-POSEDNESS for the consumer's material class.

    (i)  TT vs TDL on a PROPORTIONAL and on a ROTATING plane-strain path over the non-associated
         Drucker-Prager model.
    (ii) The BLENDED acoustic tensor `det[n . ((1-beta) C_e + beta C^ep) . n]` minimised over the
         in-plane band normal, at committed states along the inviscid path.  The inviscid
         non-associated softening tangent loses ellipticity (Rudnicki & Rice 1975); the question
         is whether the blend restores it, and at which beta -- NOT at which De, because
         `beta = 1/(1 + nsteps*De)` depends on the step count as well.
    """
    p = dp_params(c_res=c_res)
    dt = T / nsteps
    out = {"params": {k: v for k, v in p.items()}, "nsteps": nsteps}

    # ---- (i) TT vs TDL, proportional and rotating -------------------------------------
    tt_vs = []
    for rot, tag in ((0.0, "proportional"), (0.35 * np.pi, "rotating")):
        for De in (0.0,) + tuple(De_list) + (0.01, 0.1):
            tau = De * T
            _, ha = _dp_plane_strain_path(p, tau, dt, "tt", n_shear=n_shear, rotate=rot)
            _, hb = _dp_plane_strain_path(p, tau, dt, "tdl", n_shear=n_shear, rotate=rot)
            A = np.array([h["sig"] for h in ha])
            B = np.array([h["sig"] for h in hb])
            den = max(float(np.max(np.abs(B))), 1e-12)
            tt_vs.append({"path": tag, "De": De,
                          "beta": dl_beta(tau, dt),
                          "n_plastic": int(sum(h["plastic"] for h in ha)),
                          "d_path_rel": float(np.max(np.abs(A - B))) / den,
                          "d_end_rel": float(np.max(np.abs(A[-1] - B[-1]))) / den,
                          "apex": int(sum(h["apex"] for h in ha))})
    out["tt_vs_tdl"] = tt_vs
    out["prop_max"] = max(r["d_path_rel"] for r in tt_vs if r["path"] == "proportional")
    out["rot_max"] = max(r["d_path_rel"] for r in tt_vs if r["path"] == "rotating")

    # ---- (ii) the blended acoustic tensor ----------------------------------------------
    Ce = dp_elastic_C(p)
    _, hist = _dp_plane_strain_path(p, 0.0, dt, "inviscid", n_shear=n_shear)
    ell = []
    k_loss = None
    for h in hist[::5]:
        d0, th0 = min_det_acoustic(h["Cep"], Ce, ntheta=ntheta)
        if d0 <= 0.0 and k_loss is None:
            k_loss = h["k"]
        ell.append({"k": h["k"], "plastic": h["plastic"], "det_inviscid": d0, "theta": th0,
                    "alp": h["alp"]})
    out["ellipticity_path"] = ell
    out["k_loss"] = k_loss
    out["min_det_inviscid"] = min(e["det_inviscid"] for e in ell)

    # at the WORST inviscid state, sweep beta and find where the blend regains ellipticity
    kw = ell[int(np.argmin([e["det_inviscid"] for e in ell]))]["k"]
    Cep_worst = hist[kw]["Cep"]
    out["k_worst"] = kw
    sweep = []
    for b in (1.0, 0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0):
        Cb = (1.0 - b) * Ce + b * Cep_worst
        d, th = min_det_acoustic(Cb, Ce, ntheta=ntheta)
        sweep.append({"beta": b, "min_det": d, "theta": th})
    out["beta_sweep"] = sweep
    lo, hi = 0.0, 1.0                        # bisect for beta_crit (elliptic for beta < b_crit)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        Cb = (1.0 - mid) * Ce + mid * Cep_worst
        if min_det_acoustic(Cb, Ce, ntheta=min(ntheta, 181))[0] > 0.0:
            lo = mid
        else:
            hi = mid
    out["beta_crit"] = lo
    out["De_crit_at_nsteps"] = {n: (1.0 / lo - 1.0) / n for n in (250, 1000, 2000, 8000)}

    # what beta / De the requested ladder actually delivers, and whether it is elliptic
    out["de_ladder"] = []
    for De in De_list:
        b = dl_beta(De * T, dt)
        Cb = (1.0 - b) * Ce + b * Cep_worst
        d, _ = min_det_acoustic(Cb, Ce, ntheta=ntheta)
        out["de_ladder"].append({"De": De, "nsteps": nsteps, "beta": b, "min_det": d,
                                 "elliptic": bool(d > 0.0)})
    if verbose:
        print(f"  {'path':<14}{'De':>8}{'beta':>10}{'d_path_rel':>12}{'d_end_rel':>12}")
        for r in tt_vs:
            print(f"  {r['path']:<14}{r['De']:>8.4g}{r['beta']:>10.4f}"
                  f"{r['d_path_rel']:>12.2e}{r['d_end_rel']:>12.2e}")
        print(f"  inviscid non-associated DP: min normalised det Q = {out['min_det_inviscid']:.4f}"
              f"  (loss of ellipticity first at step {out['k_loss']})")
        print(f"  {'beta':>8}{'min det Q / det Q_e':>22}{'theta [deg]':>14}")
        for r in sweep:
            print(f"  {r['beta']:>8.2f}{r['min_det']:>22.5f}{np.degrees(r['theta']):>14.1f}")
        print(f"  beta_crit (elliptic for beta < this) = {out['beta_crit']:.4f}")
        print(f"  => the De that delivers it depends on the STEP COUNT: "
              + ", ".join(f"n={n}: De>{v:.2e}" for n, v in out["De_crit_at_nsteps"].items()))
        for r in out["de_ladder"]:
            print(f"  De={r['De']:.0e} at nsteps={r['nsteps']}: beta={r['beta']:.4f} "
                  f"min det={r['min_det']:.4f} elliptic={r['elliptic']}")
    return out


if __name__ == "__main__":                                    # pragma: no cover
    run_pv_gate(verbose=True)
