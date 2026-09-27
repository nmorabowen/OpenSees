"""WP-129 helpers for the SAS-ME tests (not collected: no `test_` prefix).

Material prototypes (the TIMs campaign set), the WP-128 state builders and
chain driver re-expressed over `ladrunoSANISANDReplay`, and the two
bounding-surface metrics:

  alpha_over_b_n     WP-128's: sqrt(3/2)|alpha| / alpha^b(theta of n)
                     (n the yield normal at (sigma, alpha)) -- the number its
                     tables report;
  alpha_over_b       SAS-ME's own check: the same with the Lode angle of
                     alpha itself (n-independent), returned by the replay tail.

Sign convention: everything here is the INTERNAL compression-positive one.
"""
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS = os.path.join(os.path.dirname(_HERE), "Ladruno_scripts")
if _SCRIPTS not in sys.path:
    sys.path.insert(0, _SCRIPTS)

import sanisand_replay as sr  # noqa: E402

P = list(sr.CAMPAIGN_PARAMS)
(G0, NU, E_INIT, MC, C_, LAMC, E0, KSI, PATM, M_, H0, CH, NB, A0, ND, ZMAX, CZ, DEN) = P

TAG_ME = 1        # IntScheme 1, the campaign options (today)
TAG_SAS = 2       # IntScheme 129, reseat (default), TolR 1e-4
TAG_SAS_BR = 3    # IntScheme 129, -sasAlphaIn bracket
TAG_SAS_ABL = 4   # IntScheme 129 with BOTH attribution switches (stale + stress-only)
TAG_SAS_NOE = 5   # IntScheme 129, stress-only error (G fix kept)
TAG_SAS_PRJ = 6   # IntScheme 129, -alphaProject 1
KAPPA = 0.1       # SAS-ME's default -alphaBoundTol

_COMMON = ("-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", 20000,
           "-flipAlphaIn", "init")


def sas_opts(tolr=1.0e-4, tantype=0, extra=()):
    return (129, tantype, 1, 1.0e-7, tolr) + _COMMON + tuple(extra)


def define_prototypes(ops, tolr=1.0e-4):
    ops.wipe()
    ops.nDMaterial("LadrunoSANISAND", TAG_ME, *P, *sr.CAMPAIGN_OPTS)
    ops.nDMaterial("LadrunoSANISAND", TAG_SAS, *P, *sas_opts(tolr))
    ops.nDMaterial("LadrunoSANISAND", TAG_SAS_BR, *P, *sas_opts(tolr, extra=("-sasAlphaIn", "bracket")))
    ops.nDMaterial("LadrunoSANISAND", TAG_SAS_ABL, *P,
                   *sas_opts(tolr, extra=("-sasAlphaIn", "stale", "-sasErrorVars", "stress",
                                          "-alphaBoundTol", 1.0e6)))
    # attribution prototypes: the alpha backstop is out of the way (kappa 1e6),
    # so what is measured is the ablated mechanism itself
    ops.nDMaterial("LadrunoSANISAND", TAG_SAS_NOE, *P,
                   *sas_opts(tolr, extra=("-sasErrorVars", "stress", "-alphaBoundTol", 1.0e6)))
    ops.nDMaterial("LadrunoSANISAND", TAG_SAS_PRJ, *P, *sas_opts(tolr, extra=("-alphaProject", 1)))


# ------------------------------------------------------------------ tensors
def tr(v):
    return v[0] + v[1] + v[2]


def dev(v):
    t = tr(v) / 3.0
    return [v[0] - t, v[1] - t, v[2] - t, v[3], v[4], v[5]]


def ddot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2.0 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5])


def norm(a):
    return math.sqrt(max(ddot(a, a), 0.0))


def _cos3(n):
    # sqrt(6) tr(n.n.n) for a symmetric Voigt tensor
    a = [[n[0], n[3], n[5]], [n[3], n[1], n[4]], [n[5], n[4], n[2]]]
    t = 0.0
    for i in range(3):
        for j in range(3):
            for k in range(3):
                t += a[i][j] * a[j][k] * a[k][i]
    return max(-1.0, min(1.0, math.sqrt(6.0) * t))


def _alpha_b(c3, e, p):
    g = 2 * C_ / ((1 + C_) - (1 - C_) * c3)
    psi = e - (E0 - LAMC * (max(p, 1e-10) / PATM) ** KSI)
    return g * MC * math.exp(-NB * psi) - M_


def yield_f(sig, alpha):
    p = tr(sig) / 3.0
    s = dev(sig)
    return norm([s[i] - p * alpha[i] for i in range(6)]) - math.sqrt(2.0 / 3.0) * M_ * p


def alpha_over_b_n(sig, alpha, e):
    p = tr(sig) / 3.0
    if not p > 0:
        return float("nan")
    s = dev(sig)
    x = [s[i] - p * alpha[i] for i in range(6)]
    nx = norm(x)
    n = [xi / nx for xi in x] if nx > 1e-300 else [0.0] * 6
    ab = _alpha_b(_cos3(n), e, p)
    return math.sqrt(1.5) * norm(alpha) / ab if ab > 0 else float("inf")


def alpha_over_b(sig, alpha, e):
    na = norm(alpha)
    if na < 1e-14:
        return 0.0
    ab = _alpha_b(_cos3([a / na for a in alpha]), e, tr(sig) / 3.0)
    return math.sqrt(1.5) * na / ab if ab > 0 else float("inf")


def eta(sig):
    p = tr(sig) / 3.0
    return math.sqrt(1.5) * norm(dev(sig)) / p if p > 0 else float("nan")


# ------------------------------------------------------------------ states
def k0_state(p0, K0=None, e=None):
    """WP-128 drive.k0_state: plane-strain K0 start, alpha = dev/p = alpha_in, z = 0."""
    if K0 is None:
        K0 = NU / (1.0 - NU)
    sv = 3.0 * p0 / (1.0 + 2.0 * K0)
    sig = [K0 * sv, sv, K0 * sv, 0.0, 0.0, 0.0]
    p = tr(sig) / 3.0
    al = [x / p for x in dev(sig)]
    return dict(sigma=sig, alpha=al, alpha_in=list(al), z=[0.0] * 6,
                e=E_INIT if e is None else e)


def ncov(d):
    return math.sqrt(sum(x * x for x in d[:3]) + 0.5 * sum(x * x for x in d[3:]))


def step(ops, tag, st, de, prev=0.0, trace=0):
    o = sr.replay(ops, tag, st["sigma"], st["alpha"], st["alpha_in"], st["z"],
                  st["e"], de, "compressionPositive", trace=trace, prev_incr_norm=prev)
    new = dict(sigma=list(o["sigma"]), alpha=list(o["alpha"]),
               alpha_in=list(o["alpha_in"]), z=list(o["z"]), e=o["e"])
    return new, o


# WP-128 q1_attrib chains
PATHS = {
    "vertUnload": [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    "extShear": [0.3, -1.0, 0.0, 0.5, 0.0, 0.0],
}


def incs_for(d, delta, n):
    out = []
    for _ in range(3):
        out += [[delta * x for x in d]] * n
        out += [[-delta * x for x in d]] * (n // 2)
    return out


def run_chain(ops, tag, st, incs):
    """Chain of committed replays; a refused increment is NOT committed (the
    chain stays put, as a cut step would).  Returns per-increment records."""
    hist, prev = [], 0.0
    for k, de in enumerate(incs):
        new, o = step(ops, tag, st, de, prev)
        rec = dict(k=k, rc=o["rc"], substeps=int(o["stats"]["substeps"]),
                   sas=o["sas"], f=o["f_after"], p=o["p"],
                   ab_n=alpha_over_b_n(new["sigma"], new["alpha"], new["e"]),
                   ab=alpha_over_b(new["sigma"], new["alpha"], new["e"]),
                   forced=int(o["stats"]["forcedAtDTmin"]),
                   clamp=int(o["stats"]["forcedClampMc"]))
        hist.append(rec)
        if o["rc"] == 0:
            st, prev = new, ncov(de)
    return hist


def ring_probes(delta):
    """WP-128's 8-probe set, plane-strain admissible, compression-positive."""
    return {
        "isoComp": [delta, delta, 0, 0, 0, 0],
        "isoExt": [-delta, -delta, 0, 0, 0, 0],
        "shear+": [0, 0, 0, delta, 0, 0],
        "shear-": [0, 0, 0, -delta, 0, 0],
    }
