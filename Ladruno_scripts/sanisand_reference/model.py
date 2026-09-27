"""SANISAND (Dafalias & Manzari 2004) constitutive equations, written from the paper.

Reference: Y.F. Dafalias and M.T. Manzari (2004), "Simple plasticity sand model
accounting for fabric change effects", J. Eng. Mech. 130(6):622-634,
doi:10.1061/(ASCE)0733-9399(2004)130:6(622).  Equation labels below ("DM04 ...")
name the quantity as the paper's multiaxial summary table does; the equation
numbers are collected in Ladruno_implementation/134_sanisand_reference_integrator.md
(section 2), together with every place the paper is ambiguous and how it was
resolved here.

INDEPENDENCE.  Nothing in this file is transcribed from ManzariDafalias.cpp or from
WP-128's md_port.py.  The C++ was read ONLY to list the U. Washington additions to
the paper; each one is a separately switchable option in `Options` below, OFF by
default (the paper-pure model), and `Options.uw()` switches them all on for a
like-for-like comparison with the fork's C++.

CONVENTIONS (the ones the fork's C++ uses, so states can be exchanged unchanged):
  * compression POSITIVE for stress and strain (soil mechanics sign);
  * tensors are stored as 6-vectors of TENSOR components in the order
    xx, yy, zz, xy, yz, zx (sigma, alpha, alpha_in, z are all stress-like);
  * strain INPUT is Voigt with ENGINEERING shear (gamma = 2 eps_ij), as in
    `ladrunoSANISANDReplay -dStrain`;
  * kPa, with P_atm a parameter (101 kPa in the fork's campaign set).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

import numpy as np

SQ23 = math.sqrt(2.0 / 3.0)
SQ6 = math.sqrt(6.0)
I3 = np.eye(3)


# ----------------------------------------------------------------------------
# tensor helpers
# ----------------------------------------------------------------------------
def v2t(v):
    """6-vector of tensor components (xx yy zz xy yz zx) -> symmetric 3x3."""
    return np.array([[v[0], v[3], v[5]],
                     [v[3], v[1], v[4]],
                     [v[5], v[4], v[2]]], dtype=float)


def t2v(t):
    """symmetric 3x3 -> 6-vector of tensor components (xx yy zz xy yz zx)."""
    return np.array([t[0, 0], t[1, 1], t[2, 2], t[0, 1], t[1, 2], t[0, 2]],
                    dtype=float)


def strain_v2t(v):
    """Voigt strain with ENGINEERING shear -> symmetric 3x3 strain tensor."""
    return np.array([[v[0], 0.5 * v[3], 0.5 * v[5]],
                     [0.5 * v[3], v[1], 0.5 * v[4]],
                     [0.5 * v[5], 0.5 * v[4], v[2]]], dtype=float)


def strain_t2v(t):
    """symmetric 3x3 strain tensor -> Voigt with ENGINEERING shear."""
    return np.array([t[0, 0], t[1, 1], t[2, 2],
                     2.0 * t[0, 1], 2.0 * t[1, 2], 2.0 * t[0, 2]], dtype=float)


_R2 = math.sqrt(2.0)


def t2m(t):
    """symmetric 3x3 -> Mandel 6-vector (orthonormal basis, shear * sqrt2)."""
    return np.array([t[0, 0], t[1, 1], t[2, 2],
                     _R2 * t[0, 1], _R2 * t[1, 2], _R2 * t[0, 2]], dtype=float)


def m2t(m):
    s = 1.0 / _R2
    return np.array([[m[0], s * m[3], s * m[5]],
                     [s * m[3], m[1], s * m[4]],
                     [s * m[5], s * m[4], m[2]]], dtype=float)


def ddot(a, b):
    return float(np.tensordot(a, b, axes=2))


def norm(a):
    return math.sqrt(max(ddot(a, a), 0.0))


def dev(a):
    return a - (np.trace(a) / 3.0) * I3


def mac(x):
    """Macaulay bracket <x>."""
    return x if x > 0.0 else 0.0


# ----------------------------------------------------------------------------
# parameters and options
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class Params:
    """The 17 DM04 constants (+ e_init, used only by UW-convention options).

    Order of `from_opensees` = the `nDMaterial ManzariDafalias / LadrunoSANISAND`
    argument order: G0 nu e_init Mc c lambda_c e0 ksi P_atm m h0 ch nb A0 nd
    z_max cz [Den]."""
    G0: float
    nu: float
    e_init: float
    Mc: float
    c: float
    lambda_c: float
    e0: float        # e_c0 of the critical state line (DM04 notation)
    xi: float
    P_atm: float
    m: float
    h0: float
    ch: float
    nb: float
    A0: float
    nd: float
    z_max: float
    cz: float

    @classmethod
    def from_opensees(cls, vals):
        v = [float(x) for x in vals]
        return cls(*v[:17])

    def as_opensees(self, den=2.0):
        return [self.G0, self.nu, self.e_init, self.Mc, self.c, self.lambda_c,
                self.e0, self.xi, self.P_atm, self.m, self.h0, self.ch, self.nb,
                self.A0, self.nd, self.z_max, self.cz, den]


# TIMs 2d-model campaign set (attachments README; nu = 0.312885, the Jaky K0
# substitution).  e_init only matters under the UW-convention options.
CAMPAIGN = Params.from_opensees(
    [264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005,
     1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0])

# Toyoura sand, DM04 Table (calibration against Verdugo & Ishihara 1996).
# P_atm: the paper normalises by atmospheric pressure; 100 kPa used here (the
# OpenSees documentation example for this set uses 100) -- ambiguity A11 in the doc.
TOYOURA = Params.from_opensees(
    [125.0, 0.05, 0.8, 1.25, 0.712, 0.019, 0.934, 0.7, 100.0, 0.01,
     7.05, 0.968, 1.1, 0.704, 3.5, 4.0, 600.0])


@dataclass(frozen=True)
class Options:
    """Switches.  Defaults = the paper (DM04).  Every non-default value is a
    University of Washington addition/convention found in ManzariDafalias.cpp
    (or a fork seam in LadrunoSANISAND), named in the doc, section 3."""
    # --- UW constitutive additions -------------------------------------
    d_factor: bool = False          # U1: low-p dilatancy sigmoid, p < 0.05 P_atm
    p_residual: float = 0.0         # U2: p -> p + p_r in f, n, r, psi, b0, D, Kp
    p_min: float = 0.0              # U3: G, K use max(p, p_min)  (0 = off)
    g_void_ratio: str = "current"   # U4: "current" (paper) | "initial" (UW: e_init)
    void_ratio_law: str = "current"  # U5: de = -(1+e) deps_v | "initial": -(1+e_init)
    alpha_in_rule: str = "paper"    # U6: "paper" | "uw" (once per increment,
    #                                  on (alpha_n - alpha_in_n):(Ce:deps) < 0)
    h_cap: float | None = None      # U7: UW h = 1e10 when |(alpha-alpha_in):n| < 1e-10
    elastic_moduli: str = "continuous"  # U8: "continuous" (paper: G(p, e) at every
    #                                  instant) | "frozen" (UW elastic predictor and
    #                                  intersection: Ce of the committed state on the
    #                                  elastic part) | "frozen_increment" (U9, UW
    #                                  ModifiedEuler: K, G of the committed state for
    #                                  the WHOLE increment, plastic substeps included)
    # --- integrator policy (not constitutive) ---------------------------
    start_outside: str = "stop"     # f0 > tol at t = 0: "stop" | "plastic"
    p_floor: float = 1.0e-6         # stop (status p_floor) when p_eff <= this (kPa)
    ftol_rel: float = 1.0e-8        # on-surface tolerance, relative to the cone radius
    ftol_abs_start: float = 1.0e-7  # accept a given start with f <= this as on-surface
    #                                  (the C++ TolF, absolute, kPa)
    kink_events: bool = True        # restart at <z:n> and <-dε_v^p> kinks

    @classmethod
    def uw(cls, p_min=0.0101, p_residual=0.0, **kw):
        """All UW additions ON: the like-for-like comparator for the fork's C++
        RungeKutta45 (IntScheme 45), which recomputes K, G at every stage."""
        base = dict(d_factor=True, p_residual=p_residual, p_min=p_min,
                    g_void_ratio="initial", void_ratio_law="initial",
                    alpha_in_rule="uw", h_cap=1.0e10, elastic_moduli="frozen")
        base.update(kw)
        return cls(**base)

    @classmethod
    def uw_me(cls, p_min=0.0101, p_residual=0.0, **kw):
        """As `uw`, plus U9: the comparator for the fork's ModifiedEuler
        (IntScheme 1), which never recomputes K, G inside the increment."""
        kw.setdefault("elastic_moduli", "frozen_increment")
        return cls.uw(p_min=p_min, p_residual=p_residual, **kw)

    def with_(self, **kw):
        return replace(self, **kw)


# ----------------------------------------------------------------------------
# the DM04 functions
# ----------------------------------------------------------------------------
def g_lode(cos3t, c):
    """DM04 Lode-angle interpolation g(theta, c) = 2c / ((1+c) - (1-c) cos3theta).
    g = 1 in triaxial compression (cos3theta = 1), g = c in extension (-1)."""
    return 2.0 * c / ((1.0 + c) - (1.0 - c) * cos3t)


def e_critical(p, P):
    """DM04 critical state line  e_c = e_c0 - lambda_c (p_c/p_at)^xi."""
    return P.e0 - P.lambda_c * (max(p, 0.0) / P.P_atm) ** P.xi


def elastic_moduli(p, e, P, O):
    """DM04 hypo-elastic moduli:
         G = G0 p_at (2.97 - e)^2/(1+e) (p/p_at)^(1/2),
         K = 2(1+nu)/(3(1-2nu)) G.
    U3 (p_min) and U4 (e_init in place of e) are the UW conventions."""
    pG = max(p, O.p_min) if O.p_min > 0.0 else p
    pG = max(pG, 0.0)
    eG = P.e_init if O.g_void_ratio == "initial" else e
    G = P.G0 * P.P_atm * (2.97 - eG) ** 2 / (1.0 + eG) * math.sqrt(pG / P.P_atm)
    K = 2.0 * (1.0 + P.nu) / (3.0 * (1.0 - 2.0 * P.nu)) * G
    return G, K


def uw_d_factor(p, P):
    """U1: UW's low-p damping of the dilatancy, NOT in DM04.
    1 / (1 + exp(7.6349 - 7.2713 p)) for p < 0.05 P_atm (p in kPa, constants
    non-dimensionalised by 101/P_atm as in the fork, exact at P_atm = 101)."""
    if p < 0.05 * P.P_atm:
        return 1.0 / (1.0 + math.exp(7.6349 - 7.2713 * 101.0 / P.P_atm * p))
    return 1.0


def yield_f(sig, alpha, P, O):
    """DM04 yield surface  f = ||s - p alpha|| - sqrt(2/3) m p  (stress units)."""
    p = np.trace(sig) / 3.0 + O.p_residual
    return norm(dev(sig) - p * alpha) - SQ23 * P.m * p


@dataclass
class Quantities:
    p: float          # p_eff = tr(sigma)/3 + p_r
    p_true: float     # tr(sigma)/3
    s: np.ndarray
    r: np.ndarray
    x: np.ndarray     # s - p alpha
    f: float
    n: np.ndarray
    cos3t: float
    g: float
    psi: float
    ab: float         # scalar alpha^b_theta = g Mc exp(-nb psi) - m
    ad: float         # scalar alpha^d_theta = g Mc exp(nd psi) - m
    b: np.ndarray     # sqrt(2/3) ab n - alpha
    d: np.ndarray     # sqrt(2/3) ad n - alpha
    b0: float
    a: float          # (alpha - alpha_in):n
    A: float
    D: float
    B: float
    C: float
    Rdev: np.ndarray  # R' = B n - C (n^2 - I/3)
    R: np.ndarray     # R' + D/3 I
    nr: float         # n:r
    G: float
    K: float
    X: float          # Q:E:R = 2G n:R' - K D n:r
    Hs: float         # regularised denominator a*H = (2/3) p b0 b:n + a X
    rho_b: float      # ||alpha|| / (sqrt(2/3) ab)  (ab at the Lode angle of n; WP-128's measure)
    bn: float         # b:n
    rho_alpha: float  # ||alpha|| / (sqrt(2/3) ab(theta_alpha)): ab at alpha's OWN Lode
    #                   angle -- the geometric test "alpha inside the bounding surface"


def quantities(sig, alpha, z, e, alpha_in, P, O, moduli=None):
    """Every DM04 state function at one state.  `moduli` = (G, K) overrides the
    elastic moduli (U8, the frozen-predictor convention)."""
    p_true = float(np.trace(sig)) / 3.0
    p = p_true + O.p_residual
    s = sig - p_true * I3
    pp = p if p > 0.0 else float("nan")
    r = s / pp
    x = s - p * alpha
    nx = norm(x)
    f = nx - SQ23 * P.m * p
    n = x / nx if nx > 0.0 else np.zeros((3, 3))
    c3 = SQ6 * float(np.trace(n @ n @ n))
    c3 = max(-1.0, min(1.0, c3))          # Lode: cos3theta = sqrt6 tr(n^3), clamped (A4)
    g = g_lode(c3, P.c)
    psi = e - e_critical(p, P)
    ab = g * P.Mc * math.exp(-P.nb * psi) - P.m
    ad = g * P.Mc * math.exp(P.nd * psi) - P.m
    b = SQ23 * ab * n - alpha
    d = SQ23 * ad * n - alpha
    b0 = P.G0 * P.h0 * (1.0 - P.ch * e) / math.sqrt(pp / P.P_atm)
    a = ddot(alpha - alpha_in, n)
    zn = ddot(z, n)
    A = P.A0 * (1.0 + mac(zn))
    D = A * ddot(d, n)
    if O.d_factor:
        D *= uw_d_factor(p, P)
    k = (1.0 - P.c) / P.c
    B = 1.0 + 1.5 * k * g * c3
    C = 3.0 * math.sqrt(1.5) * k * g
    Rdev = B * n - C * (n @ n - I3 / 3.0)
    R = Rdev + (D / 3.0) * I3
    nr = ddot(n, r)
    if moduli is None:
        G, K = elastic_moduli(p_true, e, P, O)
    else:
        G, K = moduli
    X = 2.0 * G * ddot(n, Rdev) - K * D * nr
    bn = ddot(b, n)
    Hs = (2.0 / 3.0) * p * b0 * bn + a * X
    rho_b = norm(alpha) / (SQ23 * ab) if ab > 0.0 else float("inf")
    rho_alpha = rho_alpha_of(alpha, psi, P)
    return Quantities(p=p, p_true=p_true, s=s, r=r, x=x, f=f, n=n, cos3t=c3, g=g,
                      psi=psi, ab=ab, ad=ad, b=b, d=d, b0=b0, a=a, A=A, D=D, B=B,
                      C=C, Rdev=Rdev, R=R, nr=nr, G=G, K=K, X=X, Hs=Hs,
                      rho_b=rho_b, bn=bn, rho_alpha=rho_alpha)


def plastic_weights(q, O):
    """Return (w, hw, sgnH): L = w N, h L = hw N, sgnH = sign of the loading
    denominator H = Kp + Q:E:R.

    Paper form (h_cap None): with a = (alpha-alpha_in):n, h = b0/a,
    Kp = (2/3) p h b:n, so H = Hs/a with Hs = (2/3) p b0 b:n + a X, and
        L = N/H = a N/Hs,   h L = b0 N/Hs
    -- finite at a = 0 (h = infinity right after an alpha_in reseat): the exact
    continuous extension, no cap needed (ambiguity A7).
    U7 (h_cap): UW's h = 1e10 when |a| < 1e-10, else b0/a; L = N/(Kp + X)."""
    if O.h_cap is None:
        Hs = q.Hs
        if Hs == 0.0:
            return 0.0, float("inf"), 0
        w = q.a / Hs
        hw = q.b0 / Hs
        if q.a > 0.0:
            sg = 1 if Hs > 0.0 else -1
        elif q.a == 0.0:
            sg = 1 if Hs > 0.0 else -1
        else:
            sg = -1 if Hs > 0.0 else 1
        return w, hw, sg
    h = O.h_cap if abs(q.a) < 1.0e-10 else q.b0 / q.a
    Kp = (2.0 / 3.0) * q.p * h * q.bn
    H = Kp + q.X
    if H == 0.0:
        return 0.0, float("inf"), 0
    return 1.0 / H, h / H, (1 if H > 0.0 else -1)


def elastic_mandel(G, K):
    one = np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])
    return K * np.outer(one, one) + 2.0 * G * (np.eye(6) - np.outer(one, one) / 3.0)


def rho_alpha_of(alpha, psi, P):
    """||alpha|| over the bounding surface's radius IN alpha's OWN direction:
    the surface in alpha-space is {sqrt(2/3) ab(theta(u)) u : u unit deviatoric},
    so alpha is inside iff this is < 1 (ambiguity A9 in the doc)."""
    na = norm(alpha)
    if na == 0.0:
        return 0.0
    u = alpha / na
    c3 = max(-1.0, min(1.0, SQ6 * float(np.trace(u @ u @ u))))
    ab = g_lode(c3, P.c) * P.Mc * math.exp(-P.nb * psi) - P.m
    return na / (SQ23 * ab) if ab > 0.0 else float("inf")


def bounding_report(sig, alpha, z, e, alpha_in, P, O):
    """alpha relative to the bounding surface at the point's own Lode angle."""
    q = quantities(sig, alpha, z, e, alpha_in, P, O)
    return dict(p=q.p, f=q.f, psi=q.psi, Mb=q.ab + P.m, alpha_b=q.ab,
                rho_b=q.rho_b, rho_alpha=q.rho_alpha, b_dot_n=q.bn, alpha_dot_n=ddot(alpha, q.n),
                eta=math.sqrt(1.5) * norm(q.s) / q.p if q.p > 0 else float("nan"),
                cos3t=q.cos3t)


@dataclass
class State:
    """One material-point state (compression positive, tensor components)."""
    sigma: np.ndarray     # 3x3
    alpha: np.ndarray     # 3x3
    z: np.ndarray         # 3x3
    e: float
    alpha_in: np.ndarray  # 3x3

    @classmethod
    def from_voigt(cls, sigma, alpha, z, e, alpha_in, project=True):
        a = v2t(alpha)
        ai = v2t(alpha_in)
        zz = v2t(z)
        if project:   # alpha, alpha_in, z are deviatoric by construction
            a, ai, zz = dev(a), dev(ai), dev(zz)
        return cls(v2t(sigma), a, zz, float(e), ai)

    def copy(self):
        return State(self.sigma.copy(), self.alpha.copy(), self.z.copy(),
                     float(self.e), self.alpha_in.copy())

    def to_dict(self):
        return dict(sigma=t2v(self.sigma).tolist(), alpha=t2v(self.alpha).tolist(),
                    z=t2v(self.z).tolist(), e=self.e,
                    alpha_in=t2v(self.alpha_in).tolist())


def on_yield_state(p, eta, e, P, direction=None, alpha_in=None, z=None):
    """A state ON the yield surface with stress ratio eta (q/p) along a unit
    deviatoric `direction` (default: triaxial compression on axis x), alpha
    placed so that f = 0 exactly and n = direction."""
    if direction is None:
        direction = np.diag([2.0, -1.0, -1.0])
    nd = dev(np.asarray(direction, dtype=float))
    nd = nd / norm(nd)
    s = eta * p * math.sqrt(2.0 / 3.0) * nd        # q = sqrt(3/2)||s|| = eta p
    sig = p * I3 + s
    r = s / p
    alpha = r - SQ23 * P.m * nd
    return State(sig, alpha, np.zeros((3, 3)) if z is None else z, float(e),
                 np.zeros((3, 3)) if alpha_in is None else alpha_in)
