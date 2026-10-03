"""Round 3b (2026-10-03): the HAR energy (sheet §2.3, (S.4h)-(S.5h'')) as a drop-in replacement for the O2
oracle's BA06 `elastic()` so that the O2 return map, Jacobian, chain data and tangents run on the HAR law.

The plastic part of O2 consumes the energy only through p, q, D11, D12, D22, q/eps_s and a^e of (S.3) (sheet §2.1,
§2.4 table), all read from the `Elastic` record, so `install(P_har)` monkey-patches `o2_algo.kernel.elastic`,
`energy_psi`, `invert_elastic` and nothing else. The §9.1 scalings (F_tol, r4/|p0|) take p0 := -p_a (§2.4).

The O2 oracle is imported READ-ONLY from the WP-144 worktree (ORACLE below).

Sheet signs: compression negative, p < 0; eps* := 1/(k(1-n)) - eps_v > 0 is the domain.
"""
import math
import os
import sys
import warnings
from dataclasses import dataclass

import numpy as np

ORACLE = (r"C:/Users/nmora/Documents/Github/OpenSees/.claude/worktrees/"
          r"ladrunonorsand-implementation-review-7bbd75/Ladruno_files/testbed/norsand_oracle")
sys.path.insert(0, ORACLE)
sys.path.insert(0, os.path.join(ORACLE, "tests"))
warnings.simplefilter("ignore")
from o2_algo import kernel as K      # noqa: E402

I3 = np.eye(3)
ONES = np.ones(3)
SQ23 = math.sqrt(2.0 / 3.0)
SQ32 = math.sqrt(1.5)

# TIMs elastic constants (sheet §2.3 gate table): n = 1/2, G0 = 264.32, nu = 0.3129, e_ref = 0.6944, p_a = 101 kPa
TIMS_HAR = dict(k=1889.48104361, g=807.80387674, n=0.5, p_a=101.0)


@dataclass
class HarParams:
    k: float
    g: float
    n: float
    p_a: float

    @property
    def edge(self):            # domain edge in eps_v: eps* = 0
        return 1.0 / (self.k * (1.0 - self.n))


class _DomainError(K.EvalError):
    pass


def har_pq(H: HarParams, ev: float, es: float):
    """(S.5h) strain form; raises EvalError('elastic_domain') for eps* <= 0."""
    est = H.edge - ev
    if est <= 0.0:
        raise K.EvalError("elastic_domain")
    kn = H.k * (1.0 - H.n)
    u = math.sqrt(est * est + 3.0 * H.g * es * es / kn)
    w = (kn * u) ** (H.n / (1.0 - H.n))
    p = -H.p_a * kn * est * w
    q = 3.0 * H.g * H.p_a * es * w
    return p, q


def har_D(H: HarParams, p: float, q: float):
    """(S.5h') stress form: D11, D12, D22, q/eps_s."""
    kn = H.k * (1.0 - H.n)
    varpi2 = p * p + kn * q * q / (3.0 * H.g)
    varpi = math.sqrt(varpi2)
    fac = H.p_a * (varpi / H.p_a) ** H.n
    Z = varpi2 / (p * p)
    D11 = H.k * fac * (1.0 - H.n + H.n / Z)
    D22 = (3.0 * H.g / (1.0 - H.n)) * fac * (1.0 - H.n / Z)
    D12 = H.n * H.k * p * q * fac / varpi2
    ratio = 3.0 * H.g * fac
    return D11, D12, D22, ratio


def make_elastic(H: HarParams):
    def elastic(P, eps_e):
        ev = float(eps_e.sum())
        e = eps_e - ev / 3.0
        ne = float(np.linalg.norm(e))
        es = SQ23 * ne
        nh = e / ne if es > K.EPS_S_TOL else np.zeros(3)
        p, q = har_pq(H, ev, es)
        if p >= 0.0:
            raise K.EvalError("p_nonneg")
        D11, D12, D22, ratio = har_D(H, p, q)
        sig = p * ONES + SQ23 * q * nh
        ae = (D11 * np.outer(ONES, ONES)
              + SQ23 * D12 * (np.outer(ONES, nh) + np.outer(nh, ONES))
              + (2.0 / 3.0) * D22 * np.outer(nh, nh)
              + (2.0 * ratio / 3.0) * (I3 - np.outer(ONES, ONES) / 3.0 - np.outer(nh, nh)))
        return K.Elastic(ev, es, nh, p, q, D11, D12, D22, sig, ae)
    return elastic


def make_energy(H: HarParams):
    def energy_psi(P, eps_e):
        ev = float(eps_e.sum())
        e = eps_e - ev / 3.0
        es = SQ23 * float(np.linalg.norm(e))
        est = H.edge - ev
        kn = H.k * (1.0 - H.n)
        u = math.sqrt(est * est + 3.0 * H.g * es * es / kn)
        return H.p_a / (H.k * (2.0 - H.n)) * (kn * u) ** ((2.0 - H.n) / (1.0 - H.n))
    return energy_psi


def make_invert(H: HarParams):
    def invert_elastic(P, sig):
        """(S.5h'') closed form."""
        p = float(sig.mean())
        xi = sig - p
        R = float(np.linalg.norm(xi))
        q = SQ32 * R
        nh = xi / R if R > 0.0 else np.zeros(3)
        if p >= 0.0:
            raise ValueError("initial stress must have p < 0")
        kn = H.k * (1.0 - H.n)
        varpi = math.sqrt(p * p + kn * q * q / (3.0 * H.g))
        ev = (1.0 / kn) * (1.0 - (abs(p) / H.p_a) ** (1.0 - H.n) * (abs(p) / varpi) ** H.n)
        es = q / (3.0 * H.g * H.p_a * (varpi / H.p_a) ** H.n)
        return ev * ONES / 3.0 + SQ32 * es * nh
    return invert_elastic


_ORIG = dict(elastic=K.elastic, energy_psi=K.energy_psi, invert_elastic=K.invert_elastic)


def install(H: HarParams):
    K.elastic = make_elastic(H)
    K.energy_psi = make_energy(H)
    K.invert_elastic = make_invert(H)


def uninstall():
    K.elastic = _ORIG["elastic"]
    K.energy_psi = _ORIG["energy_psi"]
    K.invert_elastic = _ORIG["invert_elastic"]


# ----------------------------------------------------------------------------------------------- floor (S.50)
def har_floor_evf(H: HarParams, es: float, pmin: float):
    """(S.50), n = 1/2 closed form (general n: safeguarded Newton/bisection on the exact bracket).
    Returns x = varpi_f/p_a, eps_v,f, q_f, eps'_f, eps*_f."""
    kn = H.k * (1.0 - H.n)
    a = 3.0 * kn * H.g * es * es
    b = (pmin / H.p_a) ** 2
    if abs(H.n - 0.5) < 1e-15:
        x = 0.5 * (a + math.sqrt(a * a + 4.0 * b))
    else:
        f = lambda x_: x_ * x_ - a * x_ ** (2 * H.n) - b          # noqa: E731
        xs = (H.n * a) ** (1.0 / (2.0 - 2.0 * H.n)) if a > 0 else 0.0
        xhi = max((2.0 * a) ** (1.0 / (2.0 - 2.0 * H.n)) if a > 0 else 0.0, 2.0 * math.sqrt(b))
        lo, hi = xs, xhi
        x = 0.5 * (lo + hi)
        for _ in range(100):
            fx = f(x)
            if abs(fx) <= 1e-14 * (b + a * x ** (2 * H.n)):
                break
            fpx = 2 * x - 2 * H.n * a * x ** (2 * H.n - 1)
            xn = x - fx / fpx if fpx != 0 else 0.5 * (lo + hi)
            if not (lo < xn < hi):
                xn = 0.5 * (lo + hi)
            if f(xn) < 0:
                lo = xn
            else:
                hi = xn
            x = xn
    est_f = (pmin / H.p_a) / (kn * x ** H.n)
    ev_f = H.edge - est_f
    q_f = 3.0 * H.g * H.p_a * es * x ** H.n
    epsp = H.n * pmin * q_f / ((1.0 - H.n) * H.p_a ** 2 * x * x + H.n * pmin * pmin)
    return x, ev_f, q_f, epsp, est_f


def floor_op(H: HarParams, eps_e, pmin):
    """Pi_f (S.48) under HAR: returns (eps_f, active, Phi_block, deps_f_v). Out-of-domain trial is a floor event."""
    ev = float(eps_e.sum())
    e = eps_e - ev / 3.0
    ne = float(np.linalg.norm(e))
    es = SQ23 * ne
    nh = e / ne if es > K.EPS_S_TOL else np.zeros(3)
    active = True
    if ev < H.edge:                                   # in the domain: test p
        p, _ = har_pq(H, ev, es)
        if p <= -pmin * (1.0 - 1e-12):
            active = False
    if not active:
        return eps_e.copy(), False, I3.copy(), 0.0
    _, ev_f, _, epsp, _ = har_floor_evf(H, es, pmin)
    eps_f = eps_e - (ev - ev_f) / 3.0
    Phi = I3 - 1.0 / 3.0 + (1.0 / 3.0) * epsp * SQ23 * np.outer(ONES, nh)
    return eps_f, True, Phi, ev - ev_f


def off_corner_sig(P, p_s, eta_s, direction=(-1.0, -0.35, 1.35)):
    """principal stress at (p, eta) on F = 0 with a deviatoric direction away from both WW corners."""
    xi = np.array(direction, dtype=float)
    xi -= xi.mean()
    nh = xi / np.linalg.norm(xi)
    inv = K.invariants(p_s + nh)
    z, _, _ = K.zeta_y(inv.theta, P.rho, P.zeta)
    q_s = eta_s * abs(p_s) / z
    return p_s + SQ23 * q_s * nh, inv.theta, nh
