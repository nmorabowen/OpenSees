"""ADR-93 / ADR-92 -- replay THE SEAT of the O(1) IMPL-EX event at one Gauss point.

Input: `data/census/argmax_events.csv` from the TIMs Esmeralda census leg
`D-L-dl-vt-dense-q10-sp-implex-census` (engine c162833ed, `LoadControl -ds` on the
prescribed-settlement sp).  The seat is **step 331, element 4095 GP 8** (p 57 kPa,
implexDetail[0] = 0.5309, f = 0.5 exact, 240 substeps, clamp never fired).

Unlike `ring_point.csv`, the argmax rows carry the FULL 26-slot `getState`, so the
committed state at the seat is READ, not fitted -- except the deviatoric part of the
TOTAL strain, which is not dumped anywhere per step.  Its trace comes from the void
ratio exactly (`e = e_init - (1+e_init) tr(eps)`), and its deviator is taken from the
zero-increment wall census (`census_wall_ele4095.csv`, s/B 0.017842) -- an INFERENCE
that enters ONLY the error denominator (`P_atm |eps|` is ~1.7 % of `den`).

    python3.12 seat_replay.py            # everything
    python3.12 seat_replay.py --only gate|diag|probe|repro|prevent

numpy only, python3.12, no binary.
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "adr92_p0_oracle"))

from sanisand_implex_oracle import (  # noqa: E402
    I1, ONE3, ROOT23, SMALL, Abandoned, Implex, Sanisand, compliance, dd_contr,
    dev, norm_contr, single_dot, stiffness, to_contra, to_cov, trace,
)

CENSUS = os.path.join(HERE, "data", "census")

# run.json `declared.material_common` + `declared.column.flags` of the census leg.
MATERIAL = dict(G0=264.32, nu=0.3129, Mc=1.3309, c=0.71, lambda_c=0.027, e0=0.83,
                ksi=0.45, P_atm=101.0, m=0.005, h0=1.3, ch=0.968, nb=3.5, A0=0.05,
                nd=5.75, z_max=12.5, cz=1100.0, e_init=0.6271)
P_MIN, P_R = 0.0101, 1.01
SEAT_ELE, SEAT_GP = 4095, 8
SIG_COLS = ["sig11", "sig22", "sig33", "sig12", "sig23", "sig13"]
EPS_COLS = ["eps11", "eps22", "eps33", "eps12", "eps23", "eps13"]


def _rows(path):
    lines = [ln for ln in open(path, encoding="utf-8") if not ln.startswith("#")]
    return list(csv.DictReader(lines))


def seat_rows():
    rs = _rows(os.path.join(CENSUS, "argmax_events.csv"))
    return {int(r["step"]): r for r in rs
            if int(r["ele"]) == SEAT_ELE and int(r["gp"]) == SEAT_GP}


def unpack(r):
    """One argmax row -> the material's committed state (material sign convention).

    `sig*` come through `LadrunoSANISAND3D`, which negates; the 26 `getState` slots are
    the material's own internals and are NOT negated (verified: `F(sig, alpha)` is
    4.6e-11 < TolF with this pairing, 1.2e+02 with any sign flip).
    """
    g = lambda k: float(r[k])                                          # noqa: E731
    v = lambda p: np.array([g(f"{p}_{i+1}") for i in range(6)])        # noqa: E731
    return dict(step=int(g("step")), err=g("err"), f=g("f"), ds=g("ds_m"),
                s_over_B=g("s_over_B"), substeps=int(g("substeps")),
                clamp=int(g("clamp_fired")), p=-g("p_kPa"),
                sig=-np.array([g(c) for c in SIG_COLS]),
                epsE=v("epsE"), alpha=v("alpha"), z=v("z"), alpha_in=v("alpha_in"),
                e=g("e_void"), dgamma=g("dgamma"))


def wall_eps():
    """Deviatoric total strain at the seat, from the zero-increment wall census."""
    for r in _rows(os.path.join(CENSUS, "census_wall_ele4095.csv")):
        if int(r["gp"]) == SEAT_GP:
            return -np.array([float(r[c]) for c in EPS_COLS])
    raise SystemExit("seat not in wall census")


def make(**kw):
    return Sanisand(consts=MATERIAL, scheme=1, Pmin=P_MIN,
                    Presidual=kw.pop("p_r", P_R), **kw)


def eps_of(st, eps_dev_ref):
    """Total strain whose TRACE is exact (from `e`) and whose deviator is inferred."""
    tr = (MATERIAL["e_init"] - st["e"]) / (1 + MATERIAL["e_init"])
    return dev(eps_dev_ref) + ONE3 * tr * I1


def seeded(st, eps, **kw):
    m = make(**kw)
    m.set_committed(eps, st["sig"], st["epsE"], st["alpha"], st["z"], st["alpha_in"],
                    st["dgamma"])
    return m


def run_implicit(st, eps, deps, **kw):
    """ModifiedEuler from the committed state over `deps`. Returns the new state."""
    m = seeded(st, eps, **kw)
    Ce = stiffness(m.mK, m.mG)
    m.integrate(eps + deps)
    m.commit()
    return dict(sig=m.sig_n.copy(), alpha=m.alpha_n.copy(), z=m.fabric_n.copy(),
                alpha_in=m.alpha_in_n.copy(), epsE=m.epsE_n.copy(),
                dgamma=m.dGamma_n, e=m.void_ratio,
                deps_p=(m.eps_n - m.epsE_n) - (np.array(eps) - st["epsE"]),
                substeps=m.cnt.me_substeps, force=m.cnt.me_force_accept,
                clamp=m.cnt.me_lowp_clamp, Ce=Ce, K=m.mK, G=m.mG)


# ------------------------------------------------------------------ 1. the gate

def invert_deps(st, eps, sig_target, deps0=None, itmax=60, tol=1e-11):
    """Solve ModifiedEuler(state_n, deps) = sig_target for `deps` (6 unknowns).

    The strain increment at the seat is dumped nowhere, so it is RECOVERED by
    inversion.  It is over-determined in the sense that matters: the same `deps` must
    also reproduce `alpha`, `z`, `alpha_in`, `e` and `dGamma` at the next commit, none
    of which the inversion is allowed to touch -- that is the reproduction gate.
    """
    m0 = seeded(st, eps)
    Ce = stiffness(m0.mK, m0.mG)
    x = np.linalg.solve(Ce, sig_target - st["sig"]) if deps0 is None else deps0.copy()
    x = to_cov(x) if deps0 is None else x
    best, bestr = x.copy(), np.inf
    for _ in range(itmax):
        r = run_implicit(st, eps, x)["sig"] - sig_target
        rn = norm_contr(r) / max(norm_contr(sig_target), 1.0)
        if rn < bestr:
            best, bestr = x.copy(), rn
        if rn < tol:
            break
        h = max(1e-10, 1e-6 * norm_contr(to_contra(x)))
        J = np.empty((6, 6))
        for j in range(6):
            xp = x.copy()
            xp[j] += h
            J[:, j] = (run_implicit(st, eps, xp)["sig"] - sig_target - r) / h
        try:
            dx = np.linalg.lstsq(J, -r, rcond=None)[0]
        except np.linalg.LinAlgError:
            break
        lam, ok = 1.0, False
        for _ in range(12):                      # damped
            xt = x + lam * dx
            rt = run_implicit(st, eps, xt)["sig"] - sig_target
            if norm_contr(rt) < norm_contr(r):
                x, ok = xt, True
                break
            lam *= 0.5
        if not ok:
            break
    return best, bestr


# ------------------------------------------------------------------ 2. diagnosis

def diagnose(st, eps, alpha_in=None, label=""):
    m = seeded(st, eps)
    ain = m.alpha_in_n if alpha_in is None else alpha_in
    sd = m.state_dependent(st["sig"], st["alpha"], st["z"], m.void_ratio, ain)
    p = ONE3 * trace(st["sig"]) + P_R
    K, G = m.mK, m.mG
    bn = dd_contr(sd.b, sd.n)                       # (alpha_b - alpha):n
    dn = dd_contr(sd.d, sd.n)                       # (alpha_d - alpha):n
    Kp = 2.0 / 3.0 * p * sd.h * bn
    r = dev(st["sig"]) / p
    t4 = (Kp + 2.0 * G * (sd.B - sd.C * trace(single_dot(sd.n, single_dot(sd.n, sd.n))))
          - K * sd.D * dd_contr(sd.n, r))
    aain = dd_contr(st["alpha"] - ain, sd.n)
    eta = math.sqrt(1.5) * norm_contr(dev(st["sig"])) / p
    out = dict(label=label, p=p, K=K, G=G, psi=sd.psi, h=sd.h, b0=sd.b0, aain=aain,
               bn=bn, dn=dn, Kp=Kp, t4=t4, D=sd.D, A=sd.A, aB=sd.aB, aD=sd.aD,
               Mb=sd.aB + MATERIAL["m"], Md=sd.aD + MATERIAL["m"], eta=eta,
               dist_b=norm_contr(sd.b), dist_d=norm_contr(sd.d),
               e=m.void_ratio, zn=dd_contr(st["z"], sd.n), sig_norm=norm_contr(st["sig"]),
               Kp_over_G=Kp / G, t4_over_G=t4 / G)
    return out


# ------------------------------------------------------------------ helpers

def den_of(sig, eps):
    return norm_contr(sig) + MATERIAL["P_atm"] * norm_contr(eps)


def implex_err(sig_tilde, sig_impl, eps_new):
    return norm_contr(sig_tilde - sig_impl) / den_of(sig_impl, eps_new)


def fmt(x, p=4):
    return f"{x:.{p}e}" if (x and abs(x) < 1e-3 or abs(x) >= 1e5) else f"{x:.{p}f}"


# ------------------------------------------------------------------ main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="all")
    a = ap.parse_args()
    only = a.only

    rows = seat_rows()
    st331, st332 = unpack(rows[331]), unpack(rows[332])
    eps_ref = wall_eps()
    eps331, eps332 = eps_of(st331, eps_ref), eps_of(st332, eps_ref)

    print("=" * 78)
    print("SEAT  element 4095 GP 8 (1.39, 1.11, -0.39 m)  --  dense q10 census leg")
    print("=" * 78)
    for st in (st331, st332):
        m = seeded(st, eps_of(st, eps_ref))
        print(f"  row {st['step']}: s/B {st['s_over_B']:.6f}  ds {st['ds']:.4e} m  "
              f"err {st['err']:.4f}  f {st['f']:.6f}  substeps {st['substeps']}  "
              f"clamp {st['clamp']}")
        print(f"           p {ONE3*trace(st['sig']):.3f} kPa  |sig| {norm_contr(st['sig']):.3f}"
              f"  e {st['e']:.6f}  dGamma {st['dgamma']:.4e}"
              f"  F {m.get_F(st['sig'], st['alpha']):.3e}")
    print(f"  |alpha| 331 {norm_contr(st331['alpha']):.4f} -> alpha_in 332 "
          f"{norm_contr(st332['alpha_in']):.4f}   (alpha_in was RESET at step 332: "
          f"|alpha_in_332 - alpha_331| = {norm_contr(st332['alpha_in']-st331['alpha']):.2e})")
    print(f"  inferred eps deviator (wall census): |eps|_contr {norm_contr(eps331):.5e}"
          f"  -> P_atm|eps| {MATERIAL['P_atm']*norm_contr(eps331):.3f} kPa "
          f"= {100*MATERIAL['P_atm']*norm_contr(eps331)/den_of(st331['sig'], eps331):.2f} % of den")

    deps = None
    if only in ("all", "gate", "probe", "repro", "prevent"):
        print("\n" + "-" * 78)
        print("1. REPRODUCTION GATE -- invert d_eps for step 332 from state 331")
        print("-" * 78)
        deps, res = invert_deps(st331, eps331, st332["sig"])
        got = run_implicit(st331, eps331, deps)
        tr_meas = (st331["e"] - st332["e"]) / (1 + MATERIAL["e_init"])
        print(f"  residual |sig - sig_332| / |sig_332|      {res:.3e}")
        print(f"  |d_eps|_contr {norm_contr(deps):.5e}   |d_eps| / ds "
              f"{norm_contr(to_contra(deps))/st332['ds']:.3e}   ds {st332['ds']:.3e} m")
        print(f"  tr(d_eps) recovered {trace(deps):+.6e}  vs from e_void "
              f"{tr_meas:+.6e}   rel {abs(trace(deps)-tr_meas)/max(abs(tr_meas),1e-30):.3e}")
        print("  the inversion was NOT allowed to touch these -- they are the gate:")
        for k in ("alpha", "z", "alpha_in", "epsE"):
            d = norm_contr(got[k] - st332[k])
            print(f"    |{k:9s} oracle - binary| {d:.3e}   (|binary| "
                  f"{norm_contr(st332[k]):.4e},  rel {d/max(norm_contr(st332[k]),1e-30):.3e})")
        print(f"    dGamma oracle {got['dgamma']:.6e}  binary {st332['dgamma']:.6e}"
              f"   rel {abs(got['dgamma']-st332['dgamma'])/max(abs(st332['dgamma']),1e-30):.3e}")
        print(f"    e      oracle {got['e']:.8f}  binary {st332['e']:.8f}")
        print(f"    substeps oracle {got['substeps']}  binary {st332['substeps']}"
              f"   (force-accepts {got['force']}, low-p clamps {got['clamp']})")
        print(f"  |d_eps_p| {norm_contr(to_contra(got['deps_p'])):.5e}   "
              f"|d_eps_p|/|d_eps| {norm_contr(to_contra(got['deps_p']))/norm_contr(to_contra(deps)):.4f}")
        # how anomalous is that strain increment?  Same norm, same step, at the RING.
        rp = {int(float(r["step"])): -np.array([float(r[c]) for c in EPS_COLS])
              for r in _rows(os.path.join(CENSUS, "ring_point_slice_320_345.csv"))}
        dring = norm_contr(rp[332] - rp[331])
        print(f"  SEAT   |d_eps|/ds {norm_contr(deps)/st332['ds']:9.2f} /m")
        print(f"  RING   |d_eps|/ds {dring/st332['ds']:9.2f} /m  (same step 332, "
              f"same norm; its steady rows 327-328 give "
              f"{norm_contr(rp[328]-rp[327])/5.0e-5:.2f} /m)")
        print(f"  nominal ds/element (0.5 m) {st332['ds']/0.5:.3e} -> the seat carries "
              f"{norm_contr(deps)/(st332['ds']/0.5):.0f}x the kinematic strain of the step")

    if only in ("all", "diag"):
        print("\n" + "-" * 78)
        print("2. DIAGNOSIS at the row-331 committed state (no integration)")
        print("-" * 78)
        for ain, lab in ((None, "committed alpha_in_331"),
                         (st331["alpha"], "alpha_in := alpha (the reset step 332 did)")):
            d = diagnose(st331, eps331, ain, lab)
            print(f"  [{lab}]")
            print(f"    psi {d['psi']:+.6f}   e {d['e']:.6f}   p {d['p']:.3f} kPa   "
                  f"eta = q/p {d['eta']:.4f}")
            print(f"    M^b = g*Mc*exp(-nb*psi) {d['Mb']:.4f}   M^d {d['Md']:.4f}   "
                  f"eta/M^b {d['eta']/d['Mb']:.4f}   eta/M^d {d['eta']/d['Md']:.4f}")
            print(f"    |alpha_b - alpha| {d['dist_b']:.4e}   |alpha_d - alpha| "
                  f"{d['dist_d']:.4e}")
            print(f"    (alpha_b - alpha):n {d['bn']:+.6e}   (alpha_d - alpha):n "
                  f"{d['dn']:+.6e}")
            print(f"    aain = (alpha - alpha_in):n {d['aain']:+.6e}   b0 {d['b0']:.4f}"
                  f"   h = b0/aain {d['h']:.6e}")
            print(f"    Kp = 2/3 p h (alpha_b-alpha):n {d['Kp']:+.6e} kPa   Kp/G "
                  f"{d['Kp_over_G']:+.6e}")
            print(f"    denominator t4 = Kp + 2G(B - C tr n^3) - K D (n:r) "
                  f"{d['t4']:+.6e}   t4/G {d['t4_over_G']:+.6e}")
            print(f"    D {d['D']:+.6f}   A {d['A']:.4f}   z:n {d['zn']:+.4f}   "
                  f"G {d['G']:.1f}  K {d['K']:.1f} kPa")

    if only in ("all", "probe") and deps is not None:
        print("\n" + "-" * 78)
        print("3. PROBE -- implicit return on a small d_eps along the step-332 direction")
        print("-" * 78)
        u = deps / norm_contr(to_contra(deps))
        sig_n = st331["sig"]
        _m = seeded(st331, eps331)
        Ce = stiffness(_m.mK, _m.mG)
        print(f"  {'|d_eps|':>10s} {'|d_eps_p|':>12s} {'ratio':>9s} {'|Ce:d_eps_p|':>13s}"
              f" {'/|sig|':>9s} {'|Ce:d_eps|':>12s} {'substeps':>9s} {'F_end':>10s}")
        for t in (1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1.28061e-4, 1e-3):
            g = run_implicit(st331, eps331, t * u)
            dp = g["deps_p"]
            nd = norm_contr(to_contra(t * u))
            ndp = norm_contr(to_contra(dp))
            cp = norm_contr(g["Ce"] @ dp)
            ce = norm_contr(g["Ce"] @ (t * u))
            print(f"  {nd:10.3e} {ndp:12.4e} {ndp/nd:9.4f} {cp:13.4e} "
                  f"{cp/norm_contr(sig_n):9.4f} {ce:12.4e} {g['substeps']:9d} "
                  f"{seeded(st331, eps331).get_F(g['sig'], g['alpha']):10.2e}")
        print("  (Ce is the FROZEN committed operator; |sig_n| = "
              f"{norm_contr(sig_n):.3f} kPa, G = {Ce[3,3]:.1f} kPa)")

    if only in ("all", "repro", "prevent") and deps is not None:
        print("\n" + "-" * 78)
        print("4. REPRODUCE the committed error of step 332 (f = 1 exact)")
        print("-" * 78)
        got = run_implicit(st331, eps331, deps)
        Ce, sig_n = got["Ce"], st331["sig"]
        den = den_of(st332["sig"], eps332)
        num_binary = st332["err"] * den
        dep_332 = got["deps_p"]
        print(f"  den = |sig_impl| + P_atm|eps| = {norm_contr(st332['sig']):.3f} + "
              f"{MATERIAL['P_atm']*norm_contr(eps332):.3f} = {den:.3f} kPa")
        print(f"  binary numerator |sig~ - sig_impl| = err*den = {num_binary:.3f} kPa "
              f"( = {100*num_binary/norm_contr(st332['sig']):.1f} % of |sig| )")
        # The EXACT identity. ManzariDafalias holds G, K at their committed values for
        # the whole substep loop, so d_sigma = Ce:(d_eps - d_eps_p) exactly and
        #   sigma~ - sigma_impl = Ce:( d_eps_p(n+1) - f d_eps_p(n) ).
        # Setting d_eps_p(n) := d_eps_p(n+1) with f = 1 must give EXACTLY zero.
        ident = norm_contr(sig_n + Ce @ (deps - dep_332) - st332["sig"])
        print(f"  identity check  sigma~ - sigma_impl = Ce:(d_eps_p(n+1) - f d_eps_p(n))"
              f"   residual with d_eps_p(n):=d_eps_p(n+1), f=1: {ident:.3e} kPa")
        print(f"  |Ce:d_eps|          {norm_contr(Ce@deps):9.4f} kPa   "
              f"|d_eps| {norm_contr(deps):.4e}")
        print(f"  |Ce:d_eps_p(332)|   {norm_contr(Ce@dep_332):9.4f} kPa   "
              f"= {100*norm_contr(Ce@dep_332)/norm_contr(st332['sig']):.1f} % of |sig|")
        lo = max(0.0, num_binary - norm_contr(Ce @ dep_332))
        hi = num_binary + norm_contr(Ce @ dep_332)
        print(f"  => |Ce:d_eps_p(331)| implied by the binary's own 0.4625: "
              f"[{lo:.3f}, {hi:.3f}] kPa = [{100*lo/norm_contr(st332['sig']):.0f}, "
              f"{100*hi/norm_contr(st332['sig']):.0f}] % of |sig|, i.e. "
              f"{lo/norm_contr(Ce@dep_332):.0f}-{hi/norm_contr(Ce@dep_332):.0f}x step 332's own")
        # calibrate: what |d_eps| at this state delivers that plastic increment?
        u = deps / norm_contr(to_contra(deps))
        t332 = norm_contr(to_contra(deps))
        lo_t, hi_t = t332, 4.0e-3
        for _ in range(45):
            mid = math.sqrt(lo_t * hi_t)
            g = run_implicit(st331, eps331, mid * u)
            e = implex_err(sig_n + Ce @ (deps - g["deps_p"]), st332["sig"], eps332)
            if e < st332["err"]:
                lo_t = mid
            else:
                hi_t = mid
        t_cal = math.sqrt(lo_t * hi_t)
        g_cal = run_implicit(st331, eps331, t_cal * u)
        dep_n = g_cal["deps_p"]
        e_cal = implex_err(sig_n + Ce @ (deps - dep_n), st332["sig"], eps332)
        print(f"  CALIBRATED d_eps_p(n): the same state reaches it at |d_eps| "
              f"{norm_contr(t_cal*u):.4e} ({norm_contr(t_cal*u)/norm_contr(deps):.2f}x "
              f"step 332's), {g_cal['substeps']} substeps (binary's step 331: 240)")
        print(f"    -> reproduced err {e_cal:.4f}  vs BINARY {st332['err']:.4f}   "
              f"|Ce:d_eps_p(n)| {norm_contr(Ce@dep_n):.3f} kPa")
        globals()["_CAL"] = (dep_n, t_cal, u, Ce, sig_n, dep_332)

    if only in ("all", "prevent") and deps is not None:
        print("\n" + "-" * 78)
        print("5. WHAT WOULD HAVE PREVENTED IT -- same state, same d_eps, four arms")
        print("-" * 78)
        dep_n, t_cal, u, Ce, sig_n, dep_332 = globals()["_CAL"]

        def rep(name, sig_t, note=""):
            e = implex_err(sig_t, st332["sig"], eps332)
            print(f"  {name:53s} err {e:.4f}  x{e/st332['err']:.3f}  {note}")
            return e

        rep("(0) as shipped: f = 1, history d_eps_p(n) (form A)",
            sig_n + Ce @ (deps - dep_n), f"binary {st332['err']:.4f}")
        rep("(i) alpha = 0.5 (f halved)", sig_n + Ce @ (deps - 0.5 * dep_n))
        # (ii) LANE E variant B == the oracle's form "T": the flow DIRECTION at the
        # elastic trial stress, magnitude f*|d_eps_p(n)| from the history.
        ix = Implex(seeded(st331, eps331), form="T", alpha=1.0)
        ix.d_eps_p = dep_n.copy()
        ix.g_n = norm_contr(to_contra(dep_n))
        rep("(ii) lane-E variant B (trial-flow direction, form T)",
            ix.extrapolate(eps331 + deps, f=1.0), f"overflow {ix.trial_overflow}")
        rep("(iii) f capped at 1 after a refusal chain",
            sig_n + Ce @ (deps - min(1.0, st332["f"]) * dep_n),
            f"f = {st332['f']:.4f} <= 1 -> INERT, identical to (0)")
        # (iv) history + dt from the last step at the leg's FULL ds that was under the
        # control tolerance: step 329, ds 2.5e-5 m, implex_err_max 4.7e-3.
        full_ds, f_full = 2.5e-5, st332["ds"] / 2.5e-5
        for lab, t_full in (("strain scales with ds", t_cal * full_ds / st332["ds"]),
                            ("strain at the ring's own 0.77/m ratio", 0.77 * full_ds)):
            try:
                gf = run_implicit(st331, eps331, t_full * u)
                rep(f"(iv) history+dt from step 329 ({lab})",
                    sig_n + Ce @ (deps - f_full * gf["deps_p"]),
                    f"f {f_full:.3e}, |d_eps| {norm_contr(t_full*u):.2e}, "
                    f"f|Ce:d_eps_p| {f_full*norm_contr(Ce@gf['deps_p']):.3f} kPa")
            except (Abandoned, OverflowError) as exc:
                print(f"  (iv) history+dt from step 329 ({lab}): ABANDONED -- {exc}")
        print(f"  floor: with NO extrapolation term at all (f = 0) the error is "
              f"{implex_err(sig_n + Ce @ deps, st332['sig'], eps332):.4f}")


if __name__ == "__main__":
    main()
