"""Independent 1D reference integrator for Lemaitre coupled ductile damage.

This is an OpenSees-free, from-the-equations re-implementation of the constitutive
model documented in Ladruno_implementation/15_lemaitre_ductile_damage_adr.md (de
Souza Neto, Peric & Owen 2008, Ch.12 sec 12.3/12.4) reduced to uniaxial stress:

  * combined hardening: Voce+linear isotropic  sig_y(p) = s0 + Qinf(1-e^{-b p}) + Hiso p
    and N-term Armstrong-Frederick (Chaboche) kinematic, saturation C_k/g_k ;
  * strain-equivalence effective stress  sigma~ = sigma/(1-D) ;
  * the effective return map is D-INDEPENDENT, so the plastic multiplier is solved
    on the undamaged effective problem (scalar backward-Euler Newton on dp) and the
    damage is recovered from the converged effective state:
      Y = sigma~^2/(2E)   (triaxiality 1/3 => R_v == 1 in 1D, exact, nu-independent),
      D_{n+1} = clamp( D_n + lchScale (Y/r)^s * max(0, p_{n+1} - max(p_n, pD)), 0, Dc),
      sigma = (1-D) sigma~ .

It is written ONLY from the equations (no peeking at OpenSees state), so reproducing
the LadrunoUniaxialJ2 output to ~1e-8 is a genuine oracle: it proves the shipped C++
solves the documented discrete scheme. The discrete AF update used here
(X_k = (X_k,n + C_k dp s)/(1 + g_k dp)) is the same backward-Euler form the material
ships; agreement is therefore tight, not merely asymptotic.

Closed-form sanity (perfectly plastic, N=0, s=1, pD=0): on the yield surface
sigma~ == s0, Y == s0^2/(2E), so D(p) = s0^2/(2E r) * p, linear in p.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field


@dataclass
class Lemaitre1D:
    """Stateful single-material-point reference. Call ``step(eps)`` with a monotone
    or arbitrary engineering-strain history; each call commits the state."""

    E: float
    # isotropic Voce+linear: sig_y(p) = s0 + Qinf(1-e^{-b p}) + Hiso p
    s0: float = 0.0
    Qinf: float = 0.0
    b: float = 0.0
    Hiso: float = 0.0
    # Chaboche AF kinematic: lists of (C_k, g_k); saturation C_k/g_k
    C: list = field(default_factory=list)
    g: list = field(default_factory=list)
    # Lemaitre damage (None => no damage, exact undamaged path)
    r: float | None = None
    s: float = 1.0
    pD: float = 0.0
    Dc: float = 1.0
    lchScale: float = 1.0

    # committed state
    epsp: float = 0.0       # plastic strain
    p: float = 0.0          # accumulated (equivalent) plastic strain pbar
    X: list = field(default_factory=list)   # backstresses
    D: float = 0.0

    def __post_init__(self):
        if not self.X:
            self.X = [0.0] * len(self.C)

    def sig_y(self, p):
        return self.s0 + self.Qinf * (1.0 - math.exp(-self.b * p)) + self.Hiso * p

    def sig_y_slope(self, p):
        return self.Qinf * self.b * math.exp(-self.b * p) + self.Hiso

    @property
    def dmg_on(self):
        return self.r is not None and self.r > 0.0

    def step(self, eps):
        """Advance to total strain ``eps``; commit; return (sigma, D, p, sigma_eff)."""
        E = self.E
        sig_tr = E * (eps - self.epsp)
        Xtot = sum(self.X)
        eta = sig_tr - Xtot
        abs_eta = abs(eta)
        syn = self.sig_y(self.p)
        f = abs_eta - syn
        tol = 1.0e-12 * max(abs_eta, syn, 1.0e-300)

        if f <= tol:
            # elastic: state frozen, prior damage still degrades
            sig_eff = sig_tr
            sigma = (1.0 - self.D) * sig_eff if self.dmg_on else sig_eff
            return sigma, self.D, self.p, sig_eff

        # plastic corrector: scalar Newton on dp (backward-Euler), fixed direction
        s_dir = 1.0 if eta >= 0.0 else -1.0
        h0 = self.sig_y_slope(self.p) + sum(self.C)
        dp = max(0.0, f / (E + h0))
        a = [self.C[k] - self.g[k] * s_dir * self.X[k] for k in range(len(self.C))]
        for _ in range(50):
            R = abs_eta - E * dp - self.sig_y(self.p + dp)
            dR = -E - self.sig_y_slope(self.p + dp)
            for k in range(len(self.C)):
                den = 1.0 + self.g[k] * dp
                R -= dp * a[k] / den
                dR -= a[k] / (den * den)
            if abs(R) <= tol:
                break
            if dR == 0.0:
                break
            dp = max(0.0, dp - R / dR)

        p_new = self.p + dp
        epsp_new = self.epsp + dp * s_dir
        X_new = [(self.X[k] + self.C[k] * dp * s_dir) / (1.0 + self.g[k] * dp)
                 for k in range(len(self.C))]
        sig_eff = sig_tr - E * dp * s_dir   # effective (undamaged) stress

        D_new = self.D
        if self.dmg_on:
            Y = sig_eff * sig_eff / (2.0 * E) if E > 0.0 else 0.0
            p_start = max(self.p, self.pD)
            dp_dmg = max(0.0, p_new - p_start)
            dDraw = self.lchScale * (Y / self.r) ** self.s * dp_dmg if (Y > 0.0 and dp_dmg > 0.0) else 0.0
            D_new = min(self.Dc, max(0.0, self.D + dDraw))

        sigma = (1.0 - D_new) * sig_eff if self.dmg_on else sig_eff

        # commit
        self.epsp, self.p, self.X, self.D = epsp_new, p_new, X_new, D_new
        return sigma, D_new, p_new, sig_eff

    def run(self, eps_target, nsteps):
        """Monotone displacement-controlled push; return list of (eps, sigma, D)."""
        out = []
        for i in range(1, nsteps + 1):
            eps = eps_target * i / nsteps
            sigma, D, _p, _se = self.step(eps)
            out.append((eps, sigma, D))
        return out


def integrate_rate(eps_hist, *, E, s0=0.0, Qinf=0.0, b=0.0, Hiso=0.0,
                   C=(), g=(), r=None, s=1.0, pD=0.0, Dc=1.0, nsub=4000):
    """INDEPENDENT forward-Euler explicit integration of the SAME constitutive *rate*
    equations (dSNPO §12.3) — deliberately a DIFFERENT numerical method from the shipped
    backward-Euler implicit return map (Lemaitre1D above / LadrunoUniaxialJ2):

      * Armstrong-Frederick is integrated **forward-Euler** here
        (dX_k = (C_k·s - g_k·X_k)·dp), NOT the implicit closed form
        X_k=(X_k,n + C_k·dp·s)/(1+g_k·dp) the material ships.
      * the plastic increment is a **single forward step** dp = Phi_trial/(E+h),
        NOT an iterated scalar Newton to consistency.
      * damage is forward-Euler on this independent (p, X, sigma~) trajectory.

    Both methods are O(Δ) consistent, so they converge to the SAME exact solution as the
    step → 0. Asserting the implicit material converges to this fine-step explicit
    integrator therefore validates that the shipped *return map* correctly solves the
    rate equations (a genuine cross-method check), not merely that two copies of one
    discrete scheme agree. (It cannot certify the *equations themselves* — that needs an
    external oracle, e.g. LS-DYNA *MAT_105 — which stays deferred.)

    Drives the piecewise-linear strain path through `eps_hist` (list of total strains)
    using `nsub` equal substeps PER SEGMENT. Returns (sigma, D) at the end of the path.
    """
    Ck, gk = list(C), list(g)
    epsp = 0.0
    p = 0.0
    X = [0.0] * len(Ck)
    D = 0.0
    dmg_on = (r is not None and r > 0.0)
    prev = 0.0

    def sig_y(pp):
        return s0 + Qinf * (1.0 - math.exp(-b * pp)) + Hiso * pp

    def sig_y_slope(pp):
        return Qinf * b * math.exp(-b * pp) + Hiso

    sig_eff = 0.0
    for target in eps_hist:
        for i in range(1, nsub + 1):
            eps = prev + (target - prev) * i / nsub
            sig_eff = E * (eps - epsp)
            xi = sig_eff - sum(X)
            sy = sig_y(p)
            if abs(xi) > sy:                       # plastic — one forward step
                sdir = 1.0 if xi >= 0.0 else -1.0
                h = sig_y_slope(p) + sum(Ck[k] - gk[k] * sdir * X[k] for k in range(len(Ck)))
                dp = (abs(xi) - sy) / (E + h)
                if dp < 0.0:
                    dp = 0.0
                # forward-Euler Armstrong-Frederick (distinct from the implicit form)
                X = [X[k] + (Ck[k] * sdir - gk[k] * X[k]) * dp for k in range(len(Ck))]
                p += dp
                epsp += dp * sdir
                sig_eff = E * (eps - epsp)
                if dmg_on and p > pD:
                    Y = sig_eff * sig_eff / (2.0 * E)
                    D = min(Dc, D + (Y / r) ** s * dp)
        prev = target
    sigma = (1.0 - D) * sig_eff if dmg_on else sig_eff
    return sigma, D
