"""ADR-93 P0 -- replay the Esmeralda footing's free-surface ring Gauss point.

Input: `data/ring_point_<leg>_<kind>.csv`, ONE Gauss point (element 4047, GP 5) of the
TIMs Esmeralda strip-footing legs `D-L-dl-vt-{dense,gorini}-q10-sp-{implex,implicit}`,
every COMMITTED step: the six `implexDetail` slots plus the committed stress and strain
6-vectors.  Engine c162833ed, `LoadControl -ds` on the prescribed-settlement sp.

The committed strain path is exact, so the committed stress is a pure function of the
committed state at the first dumped step plus the strain increments.  The dump carries
NO internal state -- its own header says `available at S4: detail, refusals, stress,
strain` -- so `alpha`, `z` and `alpha_in` at the first step are RECONSTRUCTED here, by a
constrained fit, and the reproduction gate measures what that reconstruction is worth.

    python3.12 replay_ring.py --gate path        # what the dumped path is
    python3.12 replay_ring.py --gate identity    # the sigma_ref refactor is a no-op
    python3.12 replay_ring.py --gate calibrate   # fit the seed, write data/seed_*.json
    python3.12 replay_ring.py --gate repro       # the reproduction gate
    python3.12 replay_ring.py --gate I1          # substep norm floor
    python3.12 replay_ring.py --gate II1         # decoupled floors + D_factor

No OpenSees binary, no C++; numpy only, `python3.12`.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "adr92_p0_oracle"))

from sanisand_implex_oracle import (  # noqa: E402
    I1, ONE3, ROOT23, SMALL, Abandoned, Sanisand, compliance, dev, dd_contr,
    norm_contr, stiffness, to_contra, trace,
)

DATA = os.path.join(HERE, "data")

# From `run.json`'s `declared.material_common` of the sibling leg
# D-L-dl-vt-dense-q10-grow-implex (the four `-sp-` legs were still running and had not
# written their own).  `e_init` is the per-leg one; everything else is shared.
# `IntScheme 1` = INT_ModifiedEuler (ManzariDafalias.cpp:41); `TanType 2` is the
# tangent only.  `-Pmin 0.0101` and `-Presidual 1.01` are 1e-4*P_atm and 1e-2*P_atm,
# i.e. exactly what `initialize()` (:918-919) sets: these legs are the VANILLA twin.
MATERIAL = dict(G0=264.32, nu=0.3129, Mc=1.3309, c=0.71, lambda_c=0.027, e0=0.83,
                ksi=0.45, P_atm=101.0, m=0.005, h0=1.3, ch=0.968, nb=3.5, A0=0.05,
                nd=5.75, z_max=12.5, cz=1100.0)
E_INIT = dict(dense=0.6271, gorini=0.6944)
P_MIN = 0.0101
P_R = 1.01

SIG_COLS = ["sig11", "sig22", "sig33", "sig12", "sig23", "sig13"]
EPS_COLS = ["eps11", "eps22", "eps33", "eps12", "eps23", "eps13"]
LEGS = (("dense", "implex"), ("dense", "implicit"),
        ("gorini", "implex"), ("gorini", "implicit"))


def read_ring(leg, kind):
    path = os.path.join(DATA, f"ring_point_{leg}_{kind}.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#")]
    out = []
    for r in csv.DictReader(lines):
        def g(k):
            v = r[k]
            return float("nan") if v in ("", "nan", "NaN") else float(v)
        out.append(dict(step=int(float(r["step"])), s_over_B=g("s_over_B"),
                        implex_err=g("implex_err"), implex_f=g("implex_f"),
                        clamp=g("clamp_fired"), ref_total=g("ref_total"),
                        # element convention -> material convention (compression +)
                        sig=-np.array([g(c) for c in SIG_COLS]),
                        eps=-np.array([g(c) for c in EPS_COLS])))
    return out


# ------------------------------------------------------------------ the seed
#
# What the dump does NOT carry, and how each piece is recovered:
#   e        -- `e_init - (1+e_init) tr(eps)`, the material's own rule; eps IS dumped.
#   epsE     -- irrelevant to the stress path (it is carried additively and never fed
#               back), so it is set to `Ce^-1 sigma` and never fitted.
#   alpha    -- the point is loading plastically, so `f = 0` holds to TolF, which pins
#               `alpha = r - sqrt(2/3) m n` with `r = dev(sigma)/p`.  Only the
#               direction of the unit normal `n` is free (5 dof, carried as a 6-vector
#               and normalised).
#   alpha_in -- free (6).  It enters only through `aain = (alpha - alpha_in):n`, which
#               sets `h = b0/aain` and therefore the whole plastic modulus.
#   z        -- free (6); zero at `updateMaterialStage 1` and small this early.

def _seed_path(leg, kind):
    return os.path.join(DATA, f"seed_{leg}_{kind}.json")


def _state_from_theta(th, sig0, m_m, n_par=18):
    v = dev(th[:6])
    n = v / max(norm_contr(v), 1e-30)
    p = ONE3 * trace(sig0) + P_R
    alpha = dev(sig0) / p - ROOT23 * m_m * n
    if n_par == 7:
        alpha_in = alpha - th[6] * n
        fabric = np.zeros(6)
    else:
        alpha_in = th[6:12]
        fabric = th[12:18]
    return alpha, fabric, alpha_in


def _make(leg, **kw):
    return Sanisand(consts=dict(MATERIAL, e_init=E_INIT[leg]), scheme=1,
                    Pmin=P_MIN, Presidual=kw.pop("p_r", P_R), **kw)


def build(rows, leg, alpha, fabric, alpha_in, **kw):
    mat = _make(leg, **kw)
    sig0, eps0 = rows[0]["sig"], rows[0]["eps"]
    K, G = mat.elastic_moduli(sig0, E_INIT[leg])
    epsE0 = np.linalg.solve(stiffness(K, G), sig0)
    mat.set_committed(eps0, sig0, epsE0, alpha, fabric, alpha_in, 0.0)
    return mat


def naive_seed(rows, leg):
    """The seed available without any fitting: n from dev(sigma), z = 0, alpha_in = 0."""
    sig0 = rows[0]["sig"]
    p = ONE3 * trace(sig0) + P_R
    n = dev(sig0) / max(norm_contr(dev(sig0)), SMALL)
    alpha = dev(sig0) / p - ROOT23 * MATERIAL["m"] * n
    return alpha, np.zeros(6), np.zeros(6)


def _n0_from_data(rows, leg):
    """The unit normal implied by the first dumped step's plastic strain increment."""
    mat = _make(leg)
    K, G = mat.elastic_moduli(rows[0]["sig"], E_INIT[leg])
    de = rows[1]["eps"] - rows[0]["eps"]
    ds = rows[1]["sig"] - rows[0]["sig"]
    dep = de - compliance(K, G) @ ds
    n = to_contra(dev(dep))
    return n / norm_contr(n)


def calibrate(leg, kind, m_steps=(1, 3), verbose=True):
    rows = read_ring(leg, kind)
    sig0 = rows[0]["sig"]

    def residual(th, M, n_par):
        a, f, ain = _state_from_theta(th, sig0, MATERIAL["m"], n_par)
        mat = build(rows, leg, a, f, ain)
        res = []
        for k in range(1, M + 1):
            try:
                mat.integrate(rows[k]["eps"])
            except Abandoned:
                pass
            mat.commit()
            res.append(mat.sig_n - rows[k]["sig"])
        return np.concatenate(res)

    def lm(th, M, n_par, iters=25):
        Rv = residual(th, M, n_par)
        f = np.linalg.norm(Rv)
        lam = 1e-3
        for _ in range(iters):
            J = np.empty((len(Rv), n_par))
            for j in range(n_par):
                h = 1e-7 * max(abs(th[j]), 1e-4)
                tp = th.copy()
                tp[j] += h
                J[:, j] = (residual(tp, M, n_par) - Rv) / h
            A, g = J.T @ J, J.T @ Rv
            moved = False
            for _ in range(25):
                try:
                    d = np.linalg.solve(
                        A + lam * np.diag(np.maximum(np.diag(A), 1e-16)), -g)
                except np.linalg.LinAlgError:
                    lam *= 10
                    continue
                R2 = residual(th + d, M, n_par)
                if np.linalg.norm(R2) < f:
                    th, Rv, f = th + d, R2, np.linalg.norm(R2)
                    lam = max(lam * 0.3, 1e-14)
                    moved = True
                    break
                lam *= 10
            if not moved or f < 1e-13:
                break
        return th, f

    # stage 1: 6 free dof for 6 residuals -- the unit normal `n` (5) and the single
    # scalar `aain = (alpha - alpha_in):n` that sets the plastic modulus. Multi-start,
    # because `h = b0/aain` is stiff in `aain` and a single start finds local minima.
    n0 = _n0_from_data(rows, leg)
    best = None
    for t0 in (1e-3, 1e-2, 0.1, 0.43, 1.0):
        th7, r7 = lm(np.concatenate([n0, [t0]]), 1, 7, iters=25)
        if best is None or r7 < best[1]:
            best = (th7, r7)
    th7, r7 = best
    if verbose:
        nrm1 = np.linalg.norm(rows[1]["sig"])
        print(f"  {leg}/{kind}  stage 1 (n + aain, 1 step): rel = {r7 / nrm1:.3e}, "
              f"aain = {th7[6]:.4g}", flush=True)
    a, f0, ain = _state_from_theta(th7, sig0, MATERIAL["m"], 7)
    th = np.concatenate([dev(th7[:6]) / max(norm_contr(dev(th7[:6])), 1e-30),
                         ain, f0])
    hist = []
    for M in m_steps:
        th, r = lm(th, M, 18)
        nrm = np.linalg.norm(np.concatenate([rows[k]["sig"] for k in range(1, M + 1)]))
        hist.append((M, r / nrm))
        if verbose:
            print(f"  {leg}/{kind}  fit over {M} step(s): rel = {r / nrm:.3e}",
                  flush=True)
    a, fab, ain = _state_from_theta(th, sig0, MATERIAL["m"], 18)
    with open(_seed_path(leg, kind), "w", encoding="utf-8") as fh:
        json.dump(dict(theta=list(th), alpha=list(a), fabric=list(fab),
                       alpha_in=list(ain), fit=hist), fh, indent=1)
    return a, fab, ain


def load_seed(rows, leg, kind):
    p = _seed_path(leg, kind)
    if not os.path.exists(p):
        return naive_seed(rows, leg)
    d = json.load(open(p, encoding="utf-8"))
    return (np.array(d["alpha"]), np.array(d["fabric"]), np.array(d["alpha_in"]))


# ------------------------------------------------------------------ the replay

def replay(rows, leg, kind="implex", seed=None, sigma_ref=0.5, p_r=P_R, p_r_e=0.0,
           D_factor=True, n_max=None):
    a, f, ain = seed if seed is not None else load_seed(rows, leg, kind)
    mat = build(rows, leg, a, f, ain, p_r=p_r, Presidual_e=p_r_e,
                substep_stress_ref=sigma_ref, D_factor=D_factor)
    out = []
    lim = len(rows) if n_max is None else min(len(rows), n_max)
    for k in range(1, lim):
        s0, c0, fa0 = (mat.cnt.me_substeps, mat.cnt.me_lowp_clamp,
                       mat.cnt.me_force_accept)
        Kc, Gc = mat.mK, mat.mG
        ab = False
        try:
            mat.integrate(rows[k]["eps"])
        except Abandoned:
            ab = True
        mat.commit()
        sd = rows[k]["sig"]
        out.append(dict(step=rows[k]["step"], s_over_B=rows[k]["s_over_B"],
                        substeps=mat.cnt.me_substeps - s0,
                        clamps=mat.cnt.me_lowp_clamp - c0,
                        force_accept=mat.cnt.me_force_accept - fa0,
                        rel=norm_contr(mat.sig_n - sd) / norm_contr(sd),
                        sig=mat.sig_n.copy(), p=ONE3 * trace(mat.sig_n),
                        lam_min=float(min(Kc - 2 * Gc / 3 + 2 * Gc, 2 * Gc)),
                        Gc=Gc, Kc=Kc, D=_D_of(mat), abandoned=ab))
    return mat, out


def _D_of(mat):
    try:
        return mat.state_dependent(mat.sig_n, mat.alpha_n, mat.fabric_n,
                                   mat.void_ratio, mat.alpha_in_n).D
    except Abandoned:
        return float("nan")


def _lam_min(K, G):
    """Smallest eigenvalue of the isotropic elastic tangent: min(3K, 2G) -> 2G here."""
    return min(3.0 * K, 2.0 * G)


def _tab(hdr, rows):
    print("| " + " | ".join(hdr) + " |")
    print("|" + "|".join("---" for _ in hdr) + "|")
    for r in rows:
        print("| " + " | ".join(str(x) for x in r) + " |")
    print()


# ------------------------------------------------------------------ gates

def gate_path():
    print("## the dumped path\n")
    rows = read_ring("dense", "implex")
    tab = []
    for i in [0, 1, 2, 4, 9, 19, 39, 79, len(rows) - 1]:
        r = rows[i]
        p = ONE3 * trace(r["sig"])
        q = math.sqrt(1.5 * dd_contr(dev(r["sig"]), dev(r["sig"])))
        tab.append([r["step"], f"{r['s_over_B']:.3e}", f"{p:.4g}", f"{q:.4g}",
                    f"{norm_contr(r['sig']):.4g}", f"{r['implex_err']:.3e}",
                    int(r["clamp"])])
    _tab(["step", "s/B", "p [kPa]", "q [kPa]", "|sigma|", "implex_err", "clamp"], tab)
    tab = []
    for leg, kind in LEGS:
        rr = read_ring(leg, kind)
        pm = [ONE3 * trace(x["sig"]) for x in rr]
        cl = [x["clamp"] for x in rr if x["clamp"] == x["clamp"]]
        tab.append([f"{leg}/{kind}", len(rr), f"{rr[-1]['s_over_B']:.3e}",
                    f"{min(pm):.4g}", f"{max(pm):.4g}", int(sum(cl))])
    _tab(["leg", "steps dumped", "s/B reached", "min p [kPa]", "max p [kPa]",
          "clamp firings"], tab)


def gate_identity():
    print("## identity -- sigma_ref = 0.5 reproduces the hardcoded branch\n")
    out = []
    for leg, kind in (("dense", "implex"), ("gorini", "implex")):
        rows = read_ring(leg, kind)
        seed = naive_seed(rows, leg)
        _, a = replay(rows, leg, seed=seed, sigma_ref=0.5)
        _, b = replay(rows, leg, seed=seed, sigma_ref=None)
        same = all(x["sig"].tobytes() == y["sig"].tobytes() for x, y in zip(a, b))
        sn = [norm_contr(x["sig"]) for x in a]
        out.append([f"ring {leg}", len(a), f"{min(sn):.3g}", f"{max(sn):.3g}",
                    "BIT-IDENTICAL" if same else "DIFFERS"])
    for tag, p0 in (("synthetic p0 = 0.05 kPa", 0.05),
                    ("synthetic p0 = 0.15 kPa", 0.15)):
        hs = []
        for sr in (0.5, None):
            mat = Sanisand(consts=dict(MATERIAL, e_init=0.6944), scheme=1, Pmin=P_MIN,
                           Presidual=P_R, substep_stress_ref=sr)
            mat.set_committed(np.zeros(6), p0 * I1, np.zeros(6), np.zeros(6),
                              np.zeros(6), np.zeros(6), 0.0)
            e, h = np.zeros(6), []
            for _ in range(40):
                e = e + np.array([2.5e-5, 2.5e-5, -5e-5, 0.0, 0.0, 0.0])
                try:
                    mat.integrate(e)
                except Abandoned:
                    pass
                mat.commit()
                h.append(mat.sig_n.copy())
            hs.append(h)
        same = all(x.tobytes() == y.tobytes() for x, y in zip(*hs))
        sn = [norm_contr(x) for x in hs[0]]
        out.append([tag, 40, f"{min(sn):.3g}", f"{max(sn):.3g}",
                    "BIT-IDENTICAL" if same else "DIFFERS"])
    _tab(["path", "steps", "min |sigma|", "max |sigma|", "verdict"], out)


def gate_calibrate():
    print("## seed calibration\n")
    for leg, kind in LEGS:
        calibrate(leg, kind)


def gate_repro():
    print("## reproduction gate\n")
    tab = []
    for leg, kind in LEGS:
        rows = read_ring(leg, kind)
        for label, seed in (("naive", naive_seed(rows, leg)),
                            ("fitted", load_seed(rows, leg, kind))):
            _, out = replay(rows, leg, kind, seed=seed)
            rel = np.array([o["rel"] for o in out])
            tab.append([f"{leg}/{kind}", label, len(out), f"{rel[0]:.2e}",
                        f"{np.median(rel):.2e}", f"{rel.max():.2e}",
                        f"{rel[-1]:.2e}"])
    _tab(["leg", "seed", "steps", "step 1", "median", "max", "last"], tab)


def gate_I1():
    print("## I.1 -- the substep norm's stress reference\n")
    tab = []
    for leg, kind in (("dense", "implex"), ("gorini", "implex")):
        rows = read_ring(leg, kind)
        seed = naive_seed(rows, leg)      # common-mode: every arm shares this seed
        base = None
        for sr in (0.5, 5.0, 50.0):
            _, out = replay(rows, leg, kind, seed=seed, sigma_ref=sr)
            sub = np.array([o["substeps"] for o in out], float)
            if base is None:
                base, dv, ratio = out, np.zeros(len(out)), 1.0
                bsub = sub
            else:
                dv = np.array([norm_contr(o["sig"] - b["sig"]) / norm_contr(b["sig"])
                               for o, b in zip(out, base)])
                ratio = bsub.sum() / max(sub.sum(), 1)
            tab.append([f"{leg}", f"{sr:g}", f"{sr / MATERIAL['P_atm']:.3g}",
                        int(sub.sum()), f"{sub.mean():.1f}", int(sub.max()),
                        f"{ratio:.2f}x", f"{np.median(dv):.1e}", f"{dv.max():.1e}"])
    _tab(["leg", "sigma_ref [kPa]", "/P_atm", "substeps total", "mean/step",
          "max/step", "substep cut", "median d(sigma)", "max d(sigma)"], tab)


def gate_II1():
    print("## II.1 -- decoupled floors (and II.2 / D5a)\n")
    arms = [("p_r = 1.01 coupled (control)", P_R, 0.0),
            ("p_r,p = 0, p_r,e = 0", 0.0, 0.0),
            ("p_r,p = 0, p_r,e = 0.1", 0.0, 0.1),
            ("p_r,p = 0, p_r,e = 1.0", 0.0, 1.0)]
    tab, tab2 = [], []
    for leg, kind in (("dense", "implex"), ("gorini", "implex")):
        rows = read_ring(leg, kind)
        seed = naive_seed(rows, leg)      # common-mode: every arm shares this seed
        base = None
        for name, prp, pre in arms:
            _, out = replay(rows, leg, kind, seed=seed, p_r=prp, p_r_e=pre)
            sub = np.array([o["substeps"] for o in out], float)
            dv = (np.zeros(len(out)) if base is None else
                  np.array([norm_contr(o["sig"] - b["sig"]) / norm_contr(b["sig"])
                            for o, b in zip(out, base)]))
            if base is None:
                base = out
            lam = np.array([_lam_min(o["Kc"], o["Gc"]) for o in out])
            pp = np.array([o["p"] for o in out])
            tab.append([leg, name, int(sub.sum()), f"{sub.mean():.1f}",
                        f"{np.median(dv):.1e}", f"{dv.max():.1e}",
                        f"{lam.min():.5g}", f"{pp.min():.4g}", f"{pp.max():.4g}"])
            for df in (True, False):
                _, o2 = replay(rows, leg, kind, seed=seed, p_r=prp, p_r_e=pre,
                               D_factor=df)
                p2 = np.array([o["p"] for o in o2])
                D2 = np.array([o["D"] for o in o2])
                tab2.append([leg, name, "on" if df else "off",
                             int((p2 + prp < 0.05 * MATERIAL["P_atm"]).sum()),
                             f"{np.nanmin(D2):+.4f}", f"{np.nanmax(D2):+.4f}",
                             f"{p2.min():.4g}"])
    _tab(["leg", "arm", "substeps", "mean/step", "median d(sigma)", "max d(sigma)",
          "min lambda(Ce) [kPa]", "min p", "max p"], tab)
    print("### II.2 / D5a -- what the `D_factor` sigmoid does on each arm\n")
    _tab(["leg", "arm", "D_factor", "steps below 0.05 P_atm", "min D", "max D",
          "min p"], tab2)



def gate_implex():
    """Slot 0 split into numerator and denominator -- EXACT, no oracle, no seed.

    The dump's own header gives the identity the material uses
    (`LadrunoSANISAND.cpp:1533-1555`):
        implex_err = ||sigma~ - sigma_impl|| / den,
        den        = ||sigma_impl||_contr + P_atm ||eps||_contr,
    with `eps` the NEW committed strain, both norms the doubled-shear contravariant
    one -- and both `sigma_impl` and `eps` are dumped.  So the ABSOLUTE discrepancy
    `||sigma~ - sigma_impl|| = implex_err * den` is recoverable offline and exactly.
    """
    print("## IMPL-EX slot 0, numerator vs denominator (exact from the dump)" + chr(10))
    tab = []
    for leg in ("dense", "gorini"):
        rows = read_ring(leg, "implex")
        prev = 0.0
        for i, r in enumerate(rows):
            den = norm_contr(r["sig"]) + MATERIAL["P_atm"] * norm_contr(r["eps"])
            num = r["implex_err"] * den
            ds = r["s_over_B"] - prev
            prev = r["s_over_B"]
            if i in (0, 1, 2, 4, 9, 19, 39, 79, len(rows) - 1):
                tab.append([leg, r["step"], f"{ds:.3e}", f"{r['implex_err']:.3e}",
                            f"{den:.4g}", f"{num:.3e}", f"{num / max(ds, 1e-30):.3g}"])
    _tab(["leg", "step", "d(s/B)", "implex_err (slot 0)", "den [kPa]",
          "|sigma~ - sigma_impl| [kPa]", "num / d(s/B)"], tab)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", default="all")
    a = ap.parse_args()
    for name, fn in (("path", gate_path), ("identity", gate_identity),
                     ("calibrate", gate_calibrate), ("repro", gate_repro),
                     ("I1", gate_I1), ("II1", gate_II1),
                     ("implex", gate_implex)):
        if a.gate in ("all", name):
            fn()


if __name__ == "__main__":
    main()
