"""Params for the O2 algorithmic oracle of LadrunoNORSAND (sheet 144a §1.3).

Names follow the equation sheet's §1.3 table. `validate()` implements the
owner-approved refusals (sheet §11.2, G0 decision 2):
  * hard ValueError unless  N_bar <= N  and  rho/rho_bar >= (1-N)/(1-N_bar)
  * warnings.warn if rho > rho_bar (violates the paper's psi_c <= phi_c reading only)
  * GA refused outside [7/9, 1]; WW refused outside (1/2, 1] -- rho = 1/2 exactly is refused
    (owner decision 2026-10-01: the compression corner becomes a vertex); both rho and rho_bar
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass


@dataclass
class Params:
    # --- BA06 energy (sheet §2.2) ---
    p0: float = -100.0          # reference pressure (< 0)
    kappa_hat: float = 0.01     # elastic compressibility
    eps_v0: float = 0.0         # reference elastic volumetric strain (at p = p0)
    mu0: float = 5400.0         # shear modulus
    alpha0: float = 0.0         # pressure/shear coupling
    # --- yield surface / potential (sheet §5) ---
    M: float = 1.2              # critical stress ratio in compression (theta = pi/3)
    N: float = 0.4              # curvature of F on the meridian plane
    N_bar: float = 0.2          # curvature of Q
    rho: float = 0.7            # ellipticity of F
    rho_bar: float = 0.8        # ellipticity of Q
    zeta: str = "WW"            # "WW" (Willam-Warnke, default) | "GA" (Gudehus-Argyris)
    # --- dilatancy and hardening (sheet §7-8) ---
    chi: float = -3.5           # maximum-dilatancy coefficient (BA06's alpha)
    h: float = 280.0            # hardening constant
    # --- CSL (sheet §6) ---
    csl_mode: str = "paper"     # "paper": v_c = v_c0 - lam_tilde ln(-p) | "fork": e_c = e0 - lam_c (-p/p_a)^xi
    lam_tilde: float = 0.0135
    v_c0: float = 1.81
    e0: float = 0.83
    lam_c: float = 0.027
    xi: float = 0.45
    p_a: float = 101.325
    # --- Q-cap (sheet §10) ---
    cap: str = "none"           # "none" | "planar" | "smooth"
    c1: float = 0.05            # eta_1 = c1 M
    c2: float = 0.15            # eta_2 = c2 M   (planar: c1 = c2)

    # derived
    @property
    def beta(self) -> float:
        return (1.0 - self.N) / (1.0 - self.N_bar)

    @property
    def chi_bar(self) -> float:
        return self.chi / self.beta

    def validate(self) -> "Params":
        if self.zeta not in ("WW", "GA"):
            raise ValueError(f"zeta must be 'WW' or 'GA', got {self.zeta!r}")
        if self.csl_mode not in ("paper", "fork"):
            raise ValueError(f"csl_mode must be 'paper' or 'fork', got {self.csl_mode!r}")
        if self.cap not in ("none", "planar", "smooth"):
            raise ValueError(f"cap must be 'none'|'planar'|'smooth', got {self.cap!r}")
        if not (self.p0 < 0.0):
            raise ValueError("p0 must be negative (compression negative)")
        if not (self.kappa_hat > 0.0):
            raise ValueError("kappa_hat must be > 0")
        if not (self.mu0 > 0.0):
            raise ValueError("mu0 must be > 0")
        if not (self.M > 0.0):
            raise ValueError("M must be > 0")
        if not (0.0 <= self.N < 1.0):
            raise ValueError("N must satisfy 0 <= N < 1")
        if not (0.0 <= self.N_bar < 1.0):
            raise ValueError("N_bar must satisfy 0 <= N_bar < 1")
        if self.h < 0.0:
            raise ValueError("h must be >= 0")
        # shape-function admissibility (sheet §4.1-4.2, K1.3). Applied to BOTH rho and rho_bar
        # (choice; the sheet states the range for rho only, but zeta_bar uses the same function).
        # GA: [7/9, 1] (closed). WW: (1/2, 1] -- rho = 1/2 EXACTLY is refused (owner decision
        # 2026-10-01): at rho = 1/2 the Willam-Warnke ellipse degenerates and the compression
        # corner theta = pi/3 becomes a vertex of the deviatoric section.
        for name, r in (("rho", self.rho), ("rho_bar", self.rho_bar)):
            if self.zeta == "GA":
                if not (7.0 / 9.0 - 1e-15 <= r <= 1.0 + 1e-15):
                    raise ValueError(f"{name}={r} outside the convex range [7/9, 1] of zeta='GA'")
            else:
                if not (0.5 < r <= 1.0 + 1e-15):
                    raise ValueError(f"{name}={r} outside the admissible range (1/2, 1] of zeta='WW' "
                                     "(rho = 1/2 exactly is refused: the compression corner is a vertex)")
        # dissipation condition A (sheet S.39), owner-approved hard refusal
        if self.N_bar > self.N + 1e-15:
            raise ValueError(f"dissipation refusal: N_bar={self.N_bar} > N={self.N} (sheet S.39)")
        if self.rho / self.rho_bar < self.beta - 1e-15:
            raise ValueError(
                f"dissipation refusal: rho/rho_bar={self.rho/self.rho_bar:.6g} < (1-N)/(1-N_bar)={self.beta:.6g} (sheet S.39)")
        if self.rho > self.rho_bar + 1e-15:
            warnings.warn(f"rho={self.rho} > rho_bar={self.rho_bar}: violates AB06's psi_c <= phi_c reading "
                          "(dissipation under reading A is still guaranteed)", stacklevel=2)
        if self.chi > 0.0:
            warnings.warn(f"chi={self.chi} > 0: the model expects chi < 0 (D* = chi psi_i)", stacklevel=2)
        if self.csl_mode == "paper":
            if not (self.lam_tilde > 0.0):
                raise ValueError("lam_tilde must be > 0")
        else:
            if not (self.lam_c > 0.0 and self.xi > 0.0 and self.p_a > 0.0):
                raise ValueError("fork CSL needs lam_c > 0, xi > 0, p_a > 0")
        if self.cap != "none":
            if not (0.0 <= self.c1 <= self.c2 < 1.0):
                raise ValueError("cap needs 0 <= c1 <= c2 < 1")
            if self.cap == "planar" and self.c1 != self.c2:
                raise ValueError("planar cap needs c1 == c2 (= chi_cap)")
            if self.cap == "smooth" and not (self.c2 > self.c1):
                raise ValueError("smooth cap needs c2 > c1")
        return self
