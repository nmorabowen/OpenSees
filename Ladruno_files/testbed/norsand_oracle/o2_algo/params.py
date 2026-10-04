"""Params for the O2 algorithmic oracle of LadrunoNORSAND (sheet 144a §1.3, §2.4, §9.7, §10.2).

Names follow the equation sheet's §1.3 table. `validate()` implements the
owner-approved refusals (sheet §11.2, G0 decision 2; §2.4 parser refusals; §9.7 p_min; (S.56)):
  * hard ValueError unless  N_bar <= N  and  rho/rho_bar >= (1-N)/(1-N_bar)
  * warnings.warn if rho > rho_bar (violates the paper's psi_c <= phi_c reading only)
  * GA refused outside [7/9, 1]; WW refused outside (1/2, 1] -- rho = 1/2 exactly is refused
    (owner decision 2026-10-01: the compression corner becomes a vertex); both rho and rho_bar
  * energy option (sheet §2.3-§2.4, owner decision (a) 2026-10-02): energy='BA06' (default, paper mode)
    uses p0, kappa_hat, eps_v0, mu0, alpha0; energy='HAR' uses k, g, n_e, p_a and REFUSES any of the
    five BA06 values given explicitly (never ignored: HAR replaces alpha0, the other four are not read);
    BA06 refuses any of k, g, n_e given. HAR: k > 0, g > 0, 0 <= n_e < 1, p_a > 0.
    p_a is ONE flag shared by the HAR energy and the fork CSL (round 3b, A4).
  * p_min (sheet §9.7, owner decision (b)): the p' floor; None -> default 5e-3 * p_ref with
    p_ref = |p0| (BA06) or p_a (HAR); 0 switches the floor off (pre-round-3 refusals); < 0 refused.
  * (S.56), cap='smooth' ONLY (round 3b, A3): PI_SCAN_REL <= W_ramp/10 with
    W_ramp = 1 - pi_i(eta_1)/pi_i(eta_2) (> 0; 0.0606 on the K2 defaults); planar/none have no ramp.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

# the five BA06 energy fields (sheet §2.2) and the three HAR fields (sheet §2.3; p_a is shared)
BA06_FIELDS = ("p0", "kappa_hat", "eps_v0", "mu0", "alpha0")
BA06_DEFAULTS = dict(p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0)
HAR_FIELDS = ("k", "g", "n_e")
P_MIN_DEFAULT_FRAC = 5.0e-3      # default p_min = 5e-3 p_ref (sheet §1.3, §9.7)


@dataclass
class Params:
    # --- BA06 energy (sheet §2.2); None = "not given" (filled with the paper defaults under energy='BA06') ---
    p0: float | None = None          # reference pressure (< 0)            [default -100.0]
    kappa_hat: float | None = None   # elastic compressibility             [default 0.01]
    eps_v0: float | None = None      # reference elastic volumetric strain [default 0.0]
    mu0: float | None = None         # shear modulus                       [default 5400.0]
    alpha0: float | None = None      # pressure/shear coupling             [default 0.0]
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
    p_a: float = 101.325        # ONE flag: fork-CSL reference pressure AND the HAR reference pressure (§2.4, A4)
    # --- Q-cap (sheet §10) ---
    cap: str = "none"           # "none" | "planar" | "smooth"
    c1: float = 0.05            # eta_1 = c1 M
    c2: float = 0.15            # eta_2 = c2 M   (planar: c1 = c2)
    # --- energy option (sheet §2.3-§2.4) ---
    energy: str = "BA06"        # "BA06" (default, paper mode) | "HAR" (Houlsby-Amorosi-Rojas 2005, n = 1/2 for TIMs)
    k: float | None = None      # HAR bulk stiffness factor (dimensionless, > 0)
    g: float | None = None      # HAR shear stiffness factor (dimensionless, > 0)
    n_e: float | None = None    # HAR pressure exponent, 0 <= n_e < 1 (n = 1 is another closed form, not shipped)
    # --- p' floor (sheet §9.7) ---
    p_min: float | None = None  # None -> 5e-3 p_ref; 0 = off; > 0: every trial and committed state has p <= -p_min

    def __post_init__(self):
        # record what was GIVEN (sheet §2.4: refused, never ignored), then fill the BA06 paper defaults
        self._ba06_given = tuple(n for n in BA06_FIELDS if getattr(self, n) is not None)
        self._har_given = tuple(n for n in HAR_FIELDS if getattr(self, n) is not None)
        if self.energy == "BA06":
            for n, d in BA06_DEFAULTS.items():
                if getattr(self, n) is None:
                    setattr(self, n, d)
        self._p_min_given = self.p_min is not None
        if self.p_min is None:
            try:
                self.p_min = P_MIN_DEFAULT_FRAC * self.p_ref
            except (TypeError, ValueError):
                self.p_min = 0.0          # unvalidated garbage (e.g. energy typo): validate() reports it

    # derived
    @property
    def beta(self) -> float:
        return (1.0 - self.N) / (1.0 - self.N_bar)

    @property
    def chi_bar(self) -> float:
        return self.chi / self.beta

    @property
    def p_ref(self) -> float:
        """Reference pressure of the §9.1 scalings (F_tol, r4/p_ref) and of the p_min default (sheet §2.4):
        |p0| under BA06, p_a under HAR (p0 := -p_a)."""
        if self.energy == "HAR":
            return float(self.p_a)
        return abs(float(self.p0))

    @property
    def W_ramp(self) -> float:
        """(S.56) ramp width of the smooth cap in pi_i at fixed p, W_ramp = 1 - pi_i(eta_1)/pi_i(eta_2) > 0
        (round 3b, A2); 0 for planar (c1 = c2) and no cap."""
        if self.cap != "smooth":
            return 0.0
        if self.N == 0.0:
            return 1.0 - math.exp(-(self.c2 - self.c1))
        return 1.0 - ((1.0 - self.c2 * self.N) / (1.0 - self.c1 * self.N)) ** ((1.0 - self.N) / self.N)

    def validate(self) -> "Params":
        if self.energy not in ("BA06", "HAR"):
            raise ValueError(f"energy must be 'BA06' or 'HAR', got {self.energy!r}")
        if self.zeta not in ("WW", "GA"):
            raise ValueError(f"zeta must be 'WW' or 'GA', got {self.zeta!r}")
        if self.csl_mode not in ("paper", "fork"):
            raise ValueError(f"csl_mode must be 'paper' or 'fork', got {self.csl_mode!r}")
        if self.cap not in ("none", "planar", "smooth"):
            raise ValueError(f"cap must be 'none'|'planar'|'smooth', got {self.cap!r}")
        # --- energy option, sheet §2.4 parser refusals ---
        if self.energy == "HAR":
            if self._ba06_given:
                raise ValueError(f"energy='HAR' refuses the BA06 parameters {self._ba06_given}: HAR replaces alpha0 and "
                                 "p0/kappa_hat/eps_v0/mu0 are not read (sheet §2.4; refused, never ignored)")
            if self.k is None or self.g is None or self.n_e is None:
                raise ValueError("energy='HAR' needs k, g and n_e (sheet §2.3)")
            if not (self.k > 0.0):
                raise ValueError("HAR k must be > 0")
            if not (self.g > 0.0):
                raise ValueError("HAR g must be > 0")
            if not (0.0 <= self.n_e < 1.0):
                raise ValueError(f"HAR n_e={self.n_e} must satisfy 0 <= n_e < 1 (n = 1 is HAR05 eq 47-48, not shipped)")
            if not (self.p_a > 0.0):
                raise ValueError("HAR p_a must be > 0")
        else:
            if self._har_given:
                raise ValueError(f"energy='BA06' refuses the HAR parameters {self._har_given} (sheet §2.4)")
            if not (self.p0 < 0.0):
                raise ValueError("p0 must be negative (compression negative)")
            if not (self.kappa_hat > 0.0):
                raise ValueError("kappa_hat must be > 0")
            if not (self.mu0 > 0.0):
                raise ValueError("mu0 must be > 0")
        # --- p' floor (sheet §9.7) ---
        if self.p_min < 0.0:
            raise ValueError(f"p_min={self.p_min} < 0 refused (sheet §9.7; 0 switches the floor off)")
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
            if self.cap == "smooth":
                if not (self.c2 > self.c1):
                    raise ValueError("smooth cap needs c2 > c1")
                # (S.56) scan-step contract, smooth cap ONLY (round 3b, A3: planar/none have W_ramp = 0)
                from . import kernel as _K      # lazy: kernel imports this module
                if _K.PI_SCAN_REL > self.W_ramp / 10.0 + 1e-15:
                    raise ValueError(f"smooth cap too narrow for the nested scan (S.56): W_ramp={self.W_ramp:.4g} "
                                     f"needs PI_SCAN_REL={_K.PI_SCAN_REL:g} <= W_ramp/10={self.W_ramp / 10.0:.4g} "
                                     "(widen c2 - c1)")
        return self
