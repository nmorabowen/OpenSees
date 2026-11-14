"""K4 (plan 144 §5.2), LadrunoSANISAND leg: the fork's C++ DM04, single material point, chained replays.

LadrunoSANISAND is the DM04 the strip-footing campaign runs, so K4 compares NorSand with BOTH the paper's
equations (WP-134's Python oracle, in k4_crossmodel.py) and the production C++ integrator. Runs under the
build's interpreter (CPython 3.12, BUILD_GOTCHAS §4); k4_crossmodel.py reads the JSON written here.

Driver: `ladrunoSANISANDReplay` (WP-127, Ladruno_scripts/sanisand_replay.py) applies ONE strain increment to a
private copy of the material from a GIVEN state (sigma, alpha, alpha_in, z, e) and returns the new state, with
the process-wide stage flag set for the call and restored (so no state leaks between cases, unlike PM4Sand).
A path is a chain of replays, each from the previous returned state.
  TCU / PSU  pure strain: (da, -da/2, -da/2) / (da, 0, -da).
  TC         (da, x, x): the lateral strain x solves sigma_yy(x) = sigma_0 (sigma_zz equal by symmetry).
  PS         (da, 0, x): x solves sigma_zz(x) = sigma_0.
The scalar x is found by solve_lateral (secant, then bracket + Illinois; each trial is a fresh replay from the
committed state) to |sigma - sigma_0| <= 1e-6 sigma_0; steps on a response jump (the integrator's substep count
changing with x) take the closer side, up to 1e-3 sigma_0, and are counted. Compression positive throughout.

Material: Toyoura (DM04 Table 1, P_atm 101 = the harness p_a), e_init = the case's e0, with one of CONFIGS
(--config):
  sasme   SAS-ME (IntScheme 129, WP-129) with the replay wrapper's options (TolF = TolR = 1e-7, -Pmin 0.0101,
          -maxSubsteps 20000, -flipAlphaIn init): the DM04 equations, no regularization.
  me      the same with IntScheme 1 (ModifiedEuler). At 1 kPa with 0.05 % steps it returns wrong states with
          rc 0 (PS sr 32.7, TC lateral tension; quirk WP-144-intscheme-1-at-1-kpa-...); it converges at 0.005 %.
  tims    the TIMs strip baseline, verbatim from tims-workbench Tries/2d-model/model/strip.py
          (kimura_v856_dm04_fp on a build with SAS R1): IntScheme 129, TanType 0, TolF 1e-7, TolR 1e-4,
          -Presidual 0 -honorTolR 0 -flipAlphaIn init -Pmin 0.0101 -maxSubsteps 2000, and the R1 (WP-152) flags
          -sasHFloor 1.0 -sasReseatHyst 1 -sasSoftCap 0.5 -sasTensionCutoff 0.5 1.0. P_atm stays 101 here (the
          preset's TOYOURA_DM04 has 100; owner 2026-11-08: keep 101 so the CSL stays NorSand's).

    py -3.12 -S k4_sanisand.py --dist <worktree>/dist/bin --cases cases.json --out sanisand.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time

TOL_ACCEPT = 1.0e-6      # |sigma_lateral - sigma_0| / sigma_0 accepted (K4 compares at the percent level)
TOL_JUMP = 1.0e-3        # largest residual taken on a response jump (counted in jump_steps, reported)


CONFIGS = {
    "sasme": (129, 0, 1, 1.0e-7, 1.0e-7, "-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", 20000,
              "-flipAlphaIn", "init"),
    "me": (1, 0, 1, 1.0e-7, 1.0e-7, "-Presidual", 0.0, "-Pmin", 0.0101, "-maxSubsteps", 20000,
           "-flipAlphaIn", "init"),
    "tims": (129, 0, 1, 1.0e-7, 1.0e-4, "-Presidual", 0.0, "-honorTolR", 0, "-flipAlphaIn", "init",
             "-Pmin", 0.0101, "-maxSubsteps", 2000,
             "-sasHFloor", 1.0, "-sasReseatHyst", 1, "-sasSoftCap", 0.5, "-sasTensionCutoff", 0.5, 1.0),
}


def toyoura_params(e):
    # G0, nu, e_init, Mc, c, lambda_c, e0, xi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, Den
    return [125.0, 0.05, float(e), 1.25, 0.712, 0.019, 0.934, 0.7, 101.0, 0.01,
            7.05, 0.968, 1.1, 0.704, 3.5, 4.0, 600.0, 2.0]


def run_case(ops, SR, case, n, tag):
    kind, p0, e0, eps_end = case["kind"], case["p0"], case["e"], case["eps_a_end"]
    ops.nDMaterial("LadrunoSANISAND", tag, *toyoura_params(e0), *SR.CAMPAIGN_OPTS)
    st = dict(sigma=[p0] * 3 + [0.0] * 3, alpha=[0.0] * 6, alpha_in=[0.0] * 6, z=[0.0] * 6, e=e0)
    da = eps_end / n
    eps = [0.0, 0.0, 0.0]
    rows = [dict(eps=list(eps), sigma=list(st["sigma"][:3]), e=e0)]
    S = dict(st=st, prev_norm=0.0, ratio=-0.5)  # ratio: lateral/axial strain guess (incompressible start)
    n_it_max, n_jump, n_cut, status = 0, 0, 0, "ok"
    j = 1 if kind == "TC" else 2

    def increment(dax):
        """One increment of axial strain dax from S; commits into S and returns the total strain increment, or
        None (S untouched) when the lateral solve fails or the material refuses."""
        nonlocal n_it_max, n_jump
        st0 = S["st"]

        def rep(d):
            return SR.replay(ops, tag, st0["sigma"], st0["alpha"], st0["alpha_in"], st0["z"], st0["e"], d,
                             "compressionPositive", trace=0, prev_incr_norm=S["prev_norm"])

        if kind == "TCU":
            d = [dax, -0.5 * dax, -0.5 * dax, 0, 0, 0]
            r = rep(d)
        elif kind == "PSU":
            d = [dax, 0.0, -dax, 0, 0, 0]
            r = rep(d)
        else:
            def mk(x):
                return [dax, x, x, 0, 0, 0] if kind == "TC" else [dax, 0.0, x, 0, 0, 0]

            def f(x):
                rr = rep(mk(x))
                return (rr["sigma"][j] - p0 if rr["rc"] == 0 else float("nan")), rr

            x1, _, r, it, jump = solve_lateral(f, S["ratio"] * dax, dax, p0)
            n_it_max = max(n_it_max, it)
            if r is None:
                return None
            n_jump += jump
            d = mk(x1)
            S["ratio"] = x1 / dax
        if r["rc"] != 0:
            return None
        S["st"] = dict(sigma=list(r["sigma"]), alpha=list(r["alpha"]), alpha_in=list(r["alpha_in"]),
                       z=list(r["z"]), e=float(r["e"]))
        S["prev_norm"] = math.sqrt(sum(v * v for v in d))
        return d

    def advance(dax, depth=0):
        """dax, cut into 4 sub-increments (up to twice: 16) when it fails, as a solver would cut the step."""
        nonlocal n_cut
        d = increment(dax)
        if d is not None or depth == 2:
            return d
        n_cut += 1
        tot = [0.0] * 6
        for _ in range(4):
            ds = advance(dax / 4.0, depth + 1)
            if ds is None:
                return None
            tot = [a + b for a, b in zip(tot, ds)]
        return tot

    for k in range(n):
        d = advance(da)
        if d is None:
            status = f"failed at step {k + 1} after cutting it to 1/16"
            break
        eps = [eps[i] + d[i] for i in range(3)]
        rows.append(dict(eps=list(eps), sigma=list(S["st"]["sigma"][:3]), e=S["st"]["e"]))
    ops.wipe()
    return dict(status=status, rows=rows, secant_iterations_max=n_it_max, jump_steps=n_jump, cut_steps=n_cut,
                n_target=n)


def solve_lateral(f, x_guess, da, p0):
    """Root of f(x) = sigma_lateral(x) - sigma_0 for the lateral strain increment x. f increases with x (more
    lateral compression, more lateral stress) but is only PIECEWISE smooth: the C++ integrator's substep count
    can change with x, so the response can jump. Secant first (fast on the smooth pieces), then a bracket and
    Illinois false position. Accepted at |f| <= TOL_ACCEPT p0; if the bracket shrinks to round-off without
    reaching it, the step sits on a jump: the side with the smaller |f| is taken and counted (jump = 1), provided
    |f| <= TOL_JUMP p0. Returns (x, f, replay result or None, evaluations, jump)."""
    tol, n_ev = TOL_ACCEPT * p0, 0
    best = None

    def ev(x):
        nonlocal n_ev, best
        n_ev += 1
        fx, rx = f(x)
        if math.isfinite(fx) and (best is None or abs(fx) < abs(best[1])):
            best = (x, fx, rx)
        return fx, rx

    x0, x1 = x_guess, x_guess - 0.01 * abs(da) - 1e-12
    f0, _ = ev(x0)
    f1, r1 = ev(x1)
    for _ in range(12):                                     # secant
        if math.isfinite(f1) and abs(f1) <= tol:
            return x1, f1, r1, n_ev, 0
        if not (math.isfinite(f0) and math.isfinite(f1)) or f1 == f0:
            break
        x0, f0, x1 = x1, f1, x1 - f1 * (x1 - x0) / (f1 - f0)
        f1, r1 = ev(x1)
    if math.isfinite(f1) and abs(f1) <= tol:
        return x1, f1, r1, n_ev, 0
    # bracket around the best point so far, expanding by |da| steps
    xb = best[0] if best else x_guess
    fb = best[1] if best else ev(xb)[0]
    step = abs(da)
    lo = hi = None
    for _ in range(30):
        xn = xb - step if fb > 0 else xb + step
        fn, _ = ev(xn)
        if not math.isfinite(fn):
            step *= 0.5
            continue
        if (fn > 0) != (fb > 0):
            lo, hi = (xn, xb) if xn < xb else (xb, xn)
            break
        xb, fb, step = xn, fn, step * 2.0
    if lo is None:
        x, fx, rx = best if best else (x1, f1, None)
        return x, fx, None, n_ev, 0
    flo, _ = ev(lo)
    fhi, _ = ev(hi)
    side = 0
    for _ in range(80):                                     # Illinois false position
        xm = hi - fhi * (hi - lo) / (fhi - flo) if fhi != flo else 0.5 * (lo + hi)
        fm, rm = ev(xm)
        if not math.isfinite(fm):
            xm = 0.5 * (lo + hi)
            fm, rm = ev(xm)
        if abs(fm) <= tol:
            return xm, fm, rm, n_ev, 0
        if (fm > 0) == (fhi > 0):
            hi, fhi = xm, fm
            if side == 1:
                flo *= 0.5
            side = 1
        else:
            lo, flo = xm, fm
            if side == -1:
                fhi *= 0.5
            side = -1
        if hi - lo <= 1e-15 * max(1.0, abs(hi)):
            break
    x, fx, rx = best
    if abs(fx) <= TOL_JUMP * p0:
        return x, fx, rx, n_ev, 1
    return x, fx, None, n_ev, 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", required=True)
    ap.add_argument("--scripts", required=True, help="Ladruno_scripts (sanisand_replay.py)")
    ap.add_argument("--cases", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--config", default="sasme", choices=sorted(CONFIGS))
    a = ap.parse_args(argv)
    os.add_dll_directory(a.dist)
    sys.path[:0] = [a.dist, a.scripts]
    import opensees as ops
    import sanisand_replay as SR
    SR.CAMPAIGN_OPTS = CONFIGS[a.config]
    spec = json.load(open(a.cases, encoding="utf-8"))
    out = dict(build=ops.ladrunoBuild(), config=a.config, opts=[str(o) for o in SR.CAMPAIGN_OPTS], results={})
    for i, case in enumerate(spec["cases"]):
        t0 = time.perf_counter()
        res = run_case(ops, SR, case, int(case.get("n", a.n)), 100 + i)
        res["seconds"] = time.perf_counter() - t0
        out["results"][case["id"]] = res
        print(case["id"], res["status"], f"{res['seconds']:.1f} s", flush=True)
    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
