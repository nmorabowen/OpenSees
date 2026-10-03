"""Parameters of LadrunoNORSAND (equation sheet 144a, table 1.3) and the parser refusals.

Names (ASCII for the sheet's symbols):
  p0 (p_0 < 0), kappa_hat (kappa^), eps_v0 (eps^e_v0), mu0, alpha0      BA06 energy (S.4)
  M, N, N_bar, rho, rho_bar, chi (< 0, BA06's alpha), h                 F, Q, D*, hardening
  lambda_tilde, v_c0                                                    "paper" CSL (AB06 41)
  e0, lambda_c, xi, p_a                                                 "fork" CSL (S.22)
  c1, c2                                                                cap blend bounds (S.35)
  csl_mode in {"paper","fork"}, zeta in {"WW","GA"}, cap in {"none","planar","smooth"},
  cap_blend in {"quintic","cubic"} (S.35), check (run validate() at construction).

v0 (initial specific volume) is a STATE input (initial_state), not a parameter.
Units: kPa.  Compression negative.

Defaults: the AB06 section 6.1 (K2) set, case rho = rho_bar = 1, paper CSL, WW, no cap
(sheet section 14), with the fork-CSL defaults of the TIMs DM04 calibration (section 15).
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass


@dataclass
class Params:
    # BA06 energy (S.4)
    p0: float = -100.0
    kappa_hat: float = 0.01
    eps_v0: float = 0.0
    mu0: float = 5400.0
    alpha0: float = 0.0
    # yield function / potential / hardening
    M: float = 1.2
    N: float = 0.4
    N_bar: float = 0.2
    rho: float = 1.0
    rho_bar: float = 1.0
    chi: float = -3.5
    h: float = 280.0
    # paper CSL
    lambda_tilde: float = 0.0135
    v_c0: float = 1.81
    # fork CSL (DM04 form)
    e0: float = 0.83
    lambda_c: float = 0.027
    xi: float = 0.45
    p_a: float = 101.325
    # Q-cap
    c1: float = 0.05
    c2: float = 0.15
    # modes
    csl_mode: str = "paper"
    zeta: str = "WW"
    cap: str = "none"
    cap_blend: str = "quintic"
    check: bool = True

    def __post_init__(self):
        if self.check:
            self.validate()

    # derived ---------------------------------------------------------------
    @property
    def beta(self) -> float:
        """beta = (1-N)/(1-N_bar)  (sheet 1.3)."""
        return (1.0 - self.N) / (1.0 - self.N_bar)

    @property
    def chi_bar(self) -> float:
        """chi_bar = chi / beta  (S.23)."""
        return self.chi / self.beta

    # refusals --------------------------------------------------------------
    def validate(self) -> "Params":
        """Raise ValueError on a refused set; warnings.warn on rho > rho_bar (sheet 11.2)."""
        if self.csl_mode not in ("paper", "fork"):
            raise ValueError(f"csl_mode must be 'paper' or 'fork', got {self.csl_mode!r}")
        if self.zeta not in ("WW", "GA"):
            raise ValueError(f"zeta must be 'WW' or 'GA', got {self.zeta!r}")
        if self.cap not in ("none", "planar", "smooth"):
            raise ValueError(f"cap must be 'none', 'planar' or 'smooth', got {self.cap!r}")
        if self.cap_blend not in ("quintic", "cubic"):
            raise ValueError(f"cap_blend must be 'quintic' or 'cubic', got {self.cap_blend!r}")
        if not self.p0 < 0.0:
            raise ValueError("p0 must be < 0 (compression negative)")
        if not self.kappa_hat > 0.0:
            raise ValueError("kappa_hat must be > 0")
        if not self.mu0 > 0.0:
            raise ValueError("mu0 must be > 0")
        if not self.M > 0.0:
            raise ValueError("M must be > 0")
        if not (0.0 <= self.N < 1.0 and 0.0 <= self.N_bar < 1.0):
            raise ValueError("need 0 <= N < 1 and 0 <= N_bar < 1")
        if not self.chi < 0.0:
            raise ValueError("chi must be < 0 (D* = chi psi_i, BA06 2.26)")
        if not self.h >= 0.0:
            raise ValueError("h must be >= 0")
        # Lode shape: admissible ranges (sheet 4.1, 4.2, K1.3); applied to rho AND rho_bar.
        # WW: (1/2, 1] -- owner decision 2026-10-01: rho = 1/2 exactly is REFUSED, because there
        # zeta = 2 cos(theta) and the compression meridian is a VERTEX (zeta'(pi/3) = -sqrt3 != 0,
        # zeta_y unbounded).  GA: [7/9, 1] (closed; convexity limit).
        for name in ("rho", "rho_bar"):
            r = getattr(self, name)
            if self.zeta == "WW":
                if not (0.5 < r <= 1.0):
                    why = ("WW at 1/2 has a vertex at the compression meridian, owner decision "
                           "2026-10-01" if r == 0.5 else "non-convex deviatoric section")
                    raise ValueError(f"{name} = {r} outside (0.5, 1] for zeta = WW ({why})")
            elif not (7.0 / 9.0 <= r <= 1.0):
                raise ValueError(f"{name} = {r} outside [{7.0 / 9.0:.6g}, 1] for zeta = GA "
                                 "(non-convex deviatoric section)")
        # dissipation, condition A (S.39), owner-approved refusal
        if not self.N_bar <= self.N:
            raise ValueError(f"refused: N_bar = {self.N_bar} > N = {self.N} (condition A, S.39)")
        if not self.rho / self.rho_bar >= self.beta * (1.0 - 1e-15):
            raise ValueError(
                f"refused: rho/rho_bar = {self.rho / self.rho_bar:.6g} < (1-N)/(1-N_bar) = "
                f"{self.beta:.6g} (condition A, S.39: negative dissipation possible)")
        if self.rho > self.rho_bar:
            warnings.warn(f"rho = {self.rho} > rho_bar = {self.rho_bar}: violates AB06's "
                          "psi_c <= phi_c reading (sheet 11.2) but not dissipation (condition A)",
                          stacklevel=2)
        if self.csl_mode == "paper":
            if not self.lambda_tilde > 0.0:
                raise ValueError("lambda_tilde must be > 0")
        else:
            if not (self.lambda_c > 0.0 and self.xi > 0.0 and self.p_a > 0.0):
                raise ValueError("fork CSL needs lambda_c > 0, xi > 0, p_a > 0")
        if self.cap == "planar":
            if not self.c1 >= 0.0:
                raise ValueError("planar cap needs c1 >= 0 (chi_cap)")
        elif self.cap == "smooth":
            if not (0.0 <= self.c1 < self.c2):
                raise ValueError("smooth cap needs 0 <= c1 < c2")
            if not self.c2 < 1.0:
                warnings.warn("smooth cap with c2 >= 1 (eta_2 >= M): the D*-identity no longer "
                              "holds at drained peaks (sheet 10.2)", stacklevel=2)
        return self
