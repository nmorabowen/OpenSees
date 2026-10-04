"""Parameters of LadrunoNORSAND (equation sheet 144a, table 1.3) and the parser refusals.

Names (ASCII for the sheet's symbols):
  energy in {"BA06","HAR"}                                              energy option (sheet 2.4)
  p0 (p_0 < 0), kappa_hat (kappa^), eps_v0 (eps^e_v0), mu0, alpha0      BA06 energy (S.4)
  k, g, n_e (HAR's n), p_a                                              HAR energy (S.4h)
  M, N, N_bar, rho, rho_bar, chi (< 0, BA06's alpha), h                 F, Q, D*, hardening
  lambda_tilde, v_c0                                                    "paper" CSL (AB06 41)
  e0, lambda_c, xi, p_a                                                 "fork" CSL (S.22)
  c1, c2                                                                cap blend bounds (S.35)
  p_min                                                                 p' floor (sheet 9.7)
  csl_mode in {"paper","fork"}, zeta in {"WW","GA"}, cap in {"none","planar","smooth"},
  cap_blend in {"quintic","cubic"} (S.35), check (run validate() at construction).

p_a is ONE parameter shared by the HAR energy and the fork CSL (sheet 2.4 / 15, round-3b
amendment A4).  The TIMs value is 101 kPa (campaign 'Patm 101'); the dataclass default stays the
pre-round-3b 101.325 only so that no existing fork-CSL result moves (every G1 test passes p_a
explicitly) -- a HAR or TIMs run must pass p_a = 101 explicitly.

The five BA06 values default to None and are resolved at construction: under energy = "BA06" a
None takes the K2 value (-100, 0.01, 0, 5400, 0); under "HAR" any of them GIVEN is refused (sheet
2.4: refused, never ignored), and k, g, n_e given under BA06 are refused.  Under HAR p0 stays
None and every reference-pressure scaling uses `p_ref` = p_a (sheet 2.4: p0 := -p_a).

p_min (sheet 9.7): 0 = floor off (the O1 default, so that every pre-round-3 result is unchanged);
"default" = 5e-3 p_ref (the parser default of the sheet: 0.5 kPa K2, 0.505 kPa at p_a = 101);
a float > 0 is used as given; < 0 is refused.

v0 (initial specific volume) is a STATE input (initial_state), not a parameter.
Units: kPa.  Compression negative.

Defaults: the AB06 section 6.1 (K2) set, case rho = rho_bar = 1, paper CSL, WW, no cap, BA06
energy, no floor (sheet section 14), with the fork-CSL defaults of the TIMs DM04 calibration
(section 15).
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

# (S.56): the nested pi_i scan step the O2 oracle and the kernel ship; O1 has no scan, but the
# validate() refusal is part of the shared Params contract (smooth cap only, amendment A3).
PI_SCAN_REL = 1.0e-3
_BA06_DEFAULTS = dict(p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0)
_HAR_NAMES = ("k", "g", "n_e")


@dataclass
class Params:
    # BA06 energy (S.4); None = "not given" (resolved in __post_init__)
    p0: float | None = None
    kappa_hat: float | None = None
    eps_v0: float | None = None
    mu0: float | None = None
    alpha0: float | None = None
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
    # energy option (sheet 2.3-2.4); HAR parameters (None = not given)
    energy: str = "BA06"
    k: float | None = None
    g: float | None = None
    n_e: float | None = None
    # p' floor (sheet 9.7): 0 = off, "default" = 5e-3 p_ref
    p_min: float | str = 0.0

    def __post_init__(self):
        given_ba06 = [nm for nm in _BA06_DEFAULTS if getattr(self, nm) is not None]
        given_har = [nm for nm in _HAR_NAMES if getattr(self, nm) is not None]
        self._given = (tuple(given_ba06), tuple(given_har))
        if self.energy != "HAR":
            for nm, val in _BA06_DEFAULTS.items():
                if getattr(self, nm) is None:
                    setattr(self, nm, val)
        if self.check:
            self.validate()

    # derived ---------------------------------------------------------------
    @property
    def p_ref(self) -> float:
        """Reference pressure of the 9.1 scalings and of the p_min default: |p0| (BA06), p_a (HAR)."""
        return self.p_a if self.energy == "HAR" else abs(self.p0)

    @property
    def pmin(self) -> float:
        """The p' floor value in kPa (0 = off)."""
        if isinstance(self.p_min, str):
            if self.p_min != "default":
                raise ValueError(f"p_min must be a number or 'default', got {self.p_min!r}")
            return 5.0e-3 * self.p_ref
        return float(self.p_min)

    @property
    def W_ramp(self) -> float:
        """(S.56) cap-ramp width in pi_i at fixed p: 1 - pi_i(eta_1)/pi_i(eta_2) (amendment A2); 0 unless
        cap = smooth with c1 < c2."""
        if self.cap != "smooth" or not self.c1 < self.c2:
            return 0.0
        if self.N == 0.0:
            return 1.0 - math.exp(-(self.c2 - self.c1))
        return 1.0 - ((1.0 - self.c2 * self.N) / (1.0 - self.c1 * self.N)) ** ((1.0 - self.N) / self.N)

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
        if self.energy not in ("BA06", "HAR"):
            raise ValueError(f"energy must be 'BA06' or 'HAR', got {self.energy!r}")
        given_ba06, given_har = getattr(self, "_given", ((), ()))
        if self.energy == "HAR":
            # sheet 2.4: BA06 values given with -energy HAR are refused, never ignored
            if given_ba06:
                raise ValueError(f"refused: {', '.join(given_ba06)} given with energy = 'HAR' (sheet 2.4: HAR "
                                 "replaces the BA06 parameters, alpha0 included; refused, never ignored)")
            for nm in _HAR_NAMES:
                if getattr(self, nm) is None:
                    raise ValueError(f"energy = 'HAR' needs {nm}")
            if not self.k > 0.0:
                raise ValueError("HAR: k must be > 0")
            if not self.g > 0.0:
                raise ValueError("HAR: g must be > 0")
            if not 0.0 <= self.n_e < 1.0:
                raise ValueError("HAR: need 0 <= n < 1 (n = 1 is HAR05 eq 47-48, not shipped)")
            if not self.p_a > 0.0:
                raise ValueError("HAR: p_a must be > 0")
        else:
            if given_har:
                raise ValueError(f"refused: {', '.join(given_har)} given with energy = 'BA06' (sheet 2.4)")
            if not self.p0 < 0.0:
                raise ValueError("p0 must be < 0 (compression negative)")
            if not self.kappa_hat > 0.0:
                raise ValueError("kappa_hat must be > 0")
            if not self.mu0 > 0.0:
                raise ValueError("mu0 must be > 0")
        if not self.pmin >= 0.0:
            raise ValueError(f"p_min = {self.pmin} < 0 refused (sheet 9.7; 0 = floor off)")
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
            # (S.56), cap = smooth only (amendment A3): planar / none have W_ramp = 0 and are not refused
            if not PI_SCAN_REL <= self.W_ramp / 10.0:
                raise ValueError(f"refused: smooth-cap ramp width W_ramp = {self.W_ramp:.4g} < 10 PI_SCAN_REL = "
                                 f"{10 * PI_SCAN_REL:g} (S.56)")
        return self
