"""WP-128: single-point strain-path driver = a chain of committed replays.

Each increment is ONE `ladrunoSANISANDReplay` call (the C++ build) -- or one
md_port.Material.update (the validated Python port) -- started from the state
the previous call returned, with prevIncrNorm = the previous increment's
covariant norm (feeds the P2-5 reversal guard exactly as a converged step
would).  WP-127 pinned replay(step k+1 from committed k) == analysis step k+1
to 1e-9, so a chain IS a material-point analysis with a prescribed strain
history (no element, no Newton).
"""
import math

import _boot as B
from _boot import ops, sr


def k0_state(p0, K0=None, e=None):
    """Plane-strain K0 start: sigma_v = sigma_yy; sigma_xx = sigma_zz = K0 sigma_v,
    alpha = dev(sigma)/p (what -flipAlphaIn init leaves), alpha_in = alpha, z = 0."""
    if K0 is None:
        K0 = B.NU / (1.0 - B.NU)
    sv = 3.0 * p0 / (1.0 + 2.0 * K0)
    sig = [K0 * sv, sv, K0 * sv, 0.0, 0.0, 0.0]
    p = B.tr(sig) / 3.0
    al = [x / p for x in B.dev(sig)]
    return dict(sigma=sig, alpha=al, alpha_in=list(al), z=[0.0] * 6,
                e=B.E_INIT if e is None else e)


def ncov(d):
    return math.sqrt(sum(x * x for x in d[:3]) + 0.5 * sum(x * x for x in d[3:]))


def step_cpp(tag, st, de, prev, trace=0):
    o = sr.replay(ops, tag, st["sigma"], st["alpha"], st["alpha_in"], st["z"],
                  st["e"], de, "compressionPositive", trace=trace, prev_incr_norm=prev)
    s = o["stats"]
    new = dict(sigma=list(o["sigma"]), alpha=list(o["alpha"]),
               alpha_in=list(o["alpha_in"]), z=list(o["z"]), e=o["e"])
    info = dict(rc=o["rc"], path=o["path"], substeps=int(s["substeps"]),
                acc=int(s["accepted"]), rej=int(s["rejectedErr"]),
                forced=int(s["forcedAtDTmin"]), clamp=int(s["forcedClampMc"]),
                abandon=int(s["abandonedLowP"]), cap=int(s["capHits"]),
                pnReset=int(s["pnResets"]), entryPmin=int(s["entryPminClamps"]),
                f_after=o["f_after"], trace=o["trace"])
    return new, info


def step_port(mat, st, de, prev):
    o = mat.update(st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], de,
                   prev_incr_norm=prev)
    new = dict(sigma=o["sigma"], alpha=o["alpha"], alpha_in=o["alpha_in"],
               z=o["z"], e=o["e"])
    info = dict(rc=o["rc"], path=o["path"], substeps=o["substeps"], acc=o["acc"],
                rej=o["rej"], forced=o["forced"], clamp=o["clamp"],
                abandon=o["abandon"], cap=o["cap"], pnReset=o["pnReset"],
                entryPmin=o["entryPmin"], f_after=o["f_after"], trace=o["trace"],
                corrGiveUp=o["corrGiveUp"], corrLowP=o["corrLowP"],
                corrMaxIter=o["corrMaxIter"])
    return new, info


def diag(st):
    sig, al = st["sigma"], st["alpha"]
    p = B.tr(sig) / 3.0
    bd = B.bounding(sig, al, st["e"]) if p > 0 else None
    return dict(p=p, eta=B.eta_sigma(sig) if p > 0 else float("nan"),
                eta_alpha=B.eta_alpha(al),
                alpha_over_b=bd["alpha_over_b"] if bd else float("nan"),
                Mb=bd["Mb"] if bd else float("nan"),
                f=B.yield_f(sig, al), tr_alpha=B.tr(al), tr_z=B.tr(st["z"]),
                norm_z=B.norm(st["z"]))


def run(st, incs, backend="cpp", tag=None, mat=None, stop_on_fail=False,
        keep_trace=False):
    hist = []
    prev = 0.0
    for k, de in enumerate(incs):
        if backend == "cpp":
            new, info = step_cpp(tag, st, de, prev, trace=20000 if keep_trace else 0)
        else:
            new, info = step_port(mat, st, de, prev)
        if not keep_trace:
            info.pop("trace", None)
        d = diag(new)
        hist.append(dict(k=k, de=de, before=st, after=new, **info, **d))
        if info["rc"] != 0 and stop_on_fail:
            break
        if info["rc"] == 0:
            st = new
            prev = ncov(de)
        # a refused update is not committed: the chain stays at st (as a cut step would)
    return hist
