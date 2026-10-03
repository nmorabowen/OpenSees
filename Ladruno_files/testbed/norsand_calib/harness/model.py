"""Where the fitted vector, the fixed sand constants and the energy plug meet: -> oracle Params + initial state.

Fitted vector (task item 4): FIT_NAMES = (chi, h, N, N_bar, rho, rho_bar).
Fixed: the sand's CSL (fork mode, sheet S.22) and M; the elastic part from the energy plug; the Lode shape
(WW default, plan §2.4), the cap (smooth, plan §2.7) and pi_i0 from the initial state (pi0_for; default rule
'ramp_end', see Setup).
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

# common name -> per-oracle field name (tests/conftest.py _RENAME)
_RENAME = {"O1": {}, "O2": {"lambda_tilde": "lam_tilde", "lambda_c": "lam_c"}}


def oracle_module(oracle: str):
    return {"O1": o1_rate, "O2": o2_algo}[oracle]


@dataclass(frozen=True)
class Setup:
    """Everything the fit does not change. pi0_rule:
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
    pi0_rule: str = "ramp_end"
    pi0_ratio: float | None = None


def common_kwargs(setup: Setup, theta: dict, p_init: float, e_init: float) -> dict:
    """The full parameter set in common names for one element test (p_init kPa > 0, e_init void ratio)."""
    s = setup.sand
    kw = dict(M=s.M, csl_mode="fork", e0=s.e0, lambda_c=s.lambda_c, xi=s.xi, p_a=s.p_a,
              zeta=setup.zeta, cap=setup.cap, c1=setup.c1, c2=setup.c2)
    for k in FIT_NAMES:
        kw[k] = float(theta[k])
    kw.update(EN.get(setup.energy).params(s, p_init, e_init, setup.policy))
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
    if setup.pi0_rule == "ramp_end":
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
    st = oracle_module(oracle).initial_state(P, sig0, 1.0 + e_init, pi0_for(setup, theta, p_init))
    return P, st
