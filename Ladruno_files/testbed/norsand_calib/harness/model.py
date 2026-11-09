"""Where the fitted vector, the fixed sand constants and the energy plug meet: -> oracle Params + initial state.

Fitted vector (task item 4): FIT_NAMES = (chi, h, N, N_bar, rho, rho_bar).
Fixed: the sand's CSL (fork mode, sheet S.22) and M; the elastic part from the energy plug; the Lode shape
(WW default, plan §2.4), the cap (smooth, plan §2.7) and pi_i0 from the initial state (default rule 'unified' =
sheet (S.53), computed by the ORACLE's own helper; see Setup). The p' floor p_min (sheet §9.7) is RESOLVED HERE,
explicitly, and handed to BOTH oracles (round 3c): Setup.p_min None -> P_MIN_FRAC (5e-3) x p_ref, with p_ref = |p0|
(BA06, = the policy's p_rep) or p_a (HAR) as in the sheet's parser default; a number is used as given; 0 = off.
O1's dataclass default is OFF and O2's is ON (5e-3 p_ref), so leaving it to the oracles would make the two differ
silently; common_kwargs() always writes p_min, so they cannot.
"""
from __future__ import annotations

import math
import os
import sys
from dataclasses import dataclass, field

import numpy as np

from . import energy as EN
from .sand import Sand

_HERE = os.path.dirname(os.path.abspath(__file__))
ORACLE_DIR = os.path.normpath(os.path.join(_HERE, "..", "..", "norsand_oracle"))
if ORACLE_DIR not in sys.path:
    sys.path.insert(0, ORACLE_DIR)

import o1_rate  # noqa: E402
import o2_algo  # noqa: E402
import warnings  # noqa: E402

# rho > rho_bar is admissible (condition A holds) and only WARNED by the oracles (sheet §11.2); inside a fit it
# would print once per residual evaluation. Silenced here; the fit reports it per result (fit.soft_flags).
warnings.filterwarnings("ignore", message=r"rho ?= ?\S+ > rho_bar", category=UserWarning)

FIT_NAMES = ("chi", "h", "N", "N_bar", "rho", "rho_bar")
# optional: NorSand's OWN fork CSL (K3c, owner 2026-11-09); absent from theta = the sand's (DM04 Table 1) values
CSL_NAMES = ("e0", "lambda_c", "xi")

P_MIN_FRAC = 5.0e-3     # sheet §9.7 / §1.3: the parser default of the p' floor, p_min = P_MIN_FRAC * p_ref

# common name -> per-oracle field name (tests/conftest.py _RENAME)
_RENAME = {"O1": {}, "O2": {"lambda_tilde": "lam_tilde", "lambda_c": "lam_c"}}


def oracle_module(oracle: str):
    return {"O1": o1_rate, "O2": o2_algo}[oracle]


@dataclass(frozen=True)
class Setup:
    """Everything the fit does not change. p_min: the p' floor (sheet §9.7, kPa). None = the harness default,
    P_MIN_FRAC (5e-3) x p_ref resolved HERE by resolve_p_min() (|p0| under BA06, p_a under HAR) and passed to both
    oracles; a number is used as given; 0 = off. (Round 3c: was "the oracle's default", which differed between O1, off,
    and O2, on.) pi0_rule:
      'unified'    : the sheet's unified rule (S.53, owner decision (d)): the yield surface through (p_init,
                     max(eta_init, c2 M)) -- the oracle's own initial_state(pi_i0=None) (O2 pi0_rule='unified', O1
                     pi_rule='S53'); c2 = the smooth cap's upper bound, c1 for a planar cap, 0 (the apex) for no cap.
                     For the harness's isotropic starts (eta_init = 0) it is the 'ramp_end' surface below, to round-off
                     (checked: harness.check_plugs);
      'on_surface' : pi_i0 puts the yield surface through the isotropic initial stress (hydrostatic: the apex,
                     pi_i0 = p (1-N)^((1-N)/N); N = 0: p e^-1) -- 'pi_i0 from the initial state';
      'ramp_end'   : pi_i0 puts the surface through (p_init, eta_y) with eta_y = c2 M (smooth cap; planar: c1 M; no cap:
                     0 = the apex): the elastic range of the isotropic state reaches just past the cap ramp, so
                     monotonic shearing yields above the ramp. Still no free parameter (p_init and the set only).
                     Why it exists: from the apex, O2 substeps every increment that crosses the smooth-cap ramp
                     (eta < c2 M, nested pi_i fold), and where the substep level changes with theta or with the
                     lateral strain the response jumps (drivers.JUMP_ACCEPT_REL note), which stalls least_squares;
      'ratio'      : pi_i0 = pi0_ratio * p_init (|pi_i0| > |apex| is an elastic start, an overconsolidated surface)."""
    sand: Sand
    energy: str = "BA06"
    policy: EN.ElasticPolicy = field(default_factory=EN.ElasticPolicy)
    zeta: str = "WW"
    cap: str = "smooth"
    c1: float = 0.05
    c2: float = 0.15
    pi0_rule: str = "unified"
    pi0_ratio: float | None = None
    p_min: float | None = None


def resolve_p_min(setup: Setup, kw: dict) -> float:
    """The p' floor value (kPa) one element test hands to BOTH oracles: Setup.p_min if given (0 = off), else
    P_MIN_FRAC x p_ref with p_ref = |p0| under BA06 (kw['p0'] = -p_rep) and p_a under HAR (sheet §2.4: p0 := -p_a).
    The same product the O2 default computes (so every pre-round-3c O2 result is unchanged) and the one O1 only
    reaches via p_min='default'."""
    if setup.p_min is not None:
        if float(setup.p_min) < 0.0:
            raise ValueError(f"Setup.p_min = {setup.p_min} < 0 (0 switches the floor off)")
        return float(setup.p_min)
    p_ref = float(kw["p_a"]) if setup.energy == "HAR" else abs(float(kw["p0"]))
    return P_MIN_FRAC * p_ref


def common_kwargs(setup: Setup, theta: dict, p_init: float, e_init: float) -> dict:
    """The full parameter set in common names for one element test (p_init kPa > 0, e_init void ratio)."""
    s = setup.sand
    kw = dict(M=s.M, csl_mode="fork", e0=s.e0, lambda_c=s.lambda_c, xi=s.xi, p_a=s.p_a,
              zeta=setup.zeta, cap=setup.cap, c1=setup.c1, c2=setup.c2)
    for k in FIT_NAMES:
        kw[k] = float(theta[k])
    for k in CSL_NAMES:
        if k in theta:
            kw[k] = float(theta[k])
    kw.update(EN.get(setup.energy).params(s, p_init, e_init, setup.policy))
    kw["p_min"] = resolve_p_min(setup, kw)          # always explicit: O1 (default off) and O2 (default on) must agree
    return kw


def make_params(oracle: str, setup: Setup, theta: dict, p_init: float, e_init: float):
    kw = common_kwargs(setup, theta, p_init, e_init)
    plug = EN.get(setup.energy)
    ok, why = plug.available(oracle)
    if not ok:
        raise EN.EnergyUnavailable(f"{setup.energy} on {oracle}: {why}")
    kw.update(plug.oracle_extra(oracle))
    kw = {_RENAME[oracle].get(k, k): v for k, v in kw.items()}
    P = oracle_module(oracle).Params(**kw)
    if oracle == "O2":
        P.validate()
    return P


def pi0_for(setup: Setup, theta: dict, p_init: float) -> float:
    p = -abs(p_init)
    N = float(theta["N"])
    if setup.pi0_rule == "on_surface":
        return p * math.exp(-1.0) if N == 0.0 else p * (1.0 - N) ** ((1.0 - N) / N)
    if setup.pi0_rule in ("ramp_end", "unified"):          # (S.53) at an isotropic start: eta_init = 0 < c2 M
        M = setup.sand.M
        eta = (setup.c2 if setup.cap == "smooth" else setup.c1 if setup.cap == "planar" else 0.0) * M
        if N == 0.0:
            return p * math.exp(eta / M - 1.0)
        return p * ((1.0 - N) / (1.0 - eta * N / M)) ** ((1.0 - N) / N)
    if setup.pi0_rule == "ratio":
        return float(setup.pi0_ratio) * p
    raise ValueError(f"unknown pi0_rule {setup.pi0_rule!r}")


def initial(oracle: str, setup: Setup, theta: dict, p_init: float, e_init: float):
    """(Params, State) at the isotropic stress -p_init I, v0 = 1 + e_init."""
    P = make_params(oracle, setup, theta, p_init, e_init)
    sig0 = -abs(p_init) * np.eye(3)
    mod = oracle_module(oracle)
    if setup.pi0_rule == "unified":            # (S.53) by the oracle's own helper, not by the harness's copy
        if oracle == "O2":
            st = mod.initial_state(P, sig0, 1.0 + e_init, None, pi0_rule="unified")
        else:
            st = mod.initial_state(P, sig0, 1.0 + e_init, None, pi_rule="S53")
    else:
        st = mod.initial_state(P, sig0, 1.0 + e_init, pi0_for(setup, theta, p_init))
    return P, st
