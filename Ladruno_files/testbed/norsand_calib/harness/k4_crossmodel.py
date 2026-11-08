"""K4, the cross-model check (plan 144 §5.2, soft gate): single-point drained and undrained triaxial and plane-strain
compression of Toyoura, NorSand against DM04: the paper's equations (Python) and the fork's C++ (LadrunoSANISAND).

Models
  NorSand  O2 (the C++ kernel's contract) through drivers.simulate, with K3's setup: HAR energy, p' floor
           P_MIN_FRAC x p_a, unified pi_i0 rule, WW, smooth cap; parameters theta from --theta (a K3 result JSON or
           a dict). Until K3 reports, --theta provisional uses smoke_recovery.THETA_TRUE: a DRY RUN of the plumbing,
           not a K4 verdict.
  DM04     WP-134's independent reference integrator (Ladruno_scripts/sanisand_reference), Toyoura Table 1 with
           P_atm = the harness p_a (101 kPa), so its CSL e_c = e0 - lambda_c (p/p_a)^xi and M are EXACTLY NorSand's
           (sand.toyoura_dm04, sheet S.22; checked below). Two option sets: 'paper' (the published model, primary)
           and 'uw_model' (the plan's named oracle: UW additions U1-U5).
  LadrunoSANISAND  the fork's C++ DM04 (the material the strip campaign runs), one material point driven by
           chained ladrunoSANISANDReplay calls in k4_sanisand.py under Python 3.12; same Toyoura set and P_atm,
           the campaign's integration options. Its distance from DM04 'paper' is reported too: a check of the C++
           against the paper's equations on these paths.
  (PM4Sand: k4_pm4sand.py is kept as a standalone driver; owner 2026-11-07: K4 sticks to LadrunoSANISAND.)

Cases: p0 in {100, 500} kPa x e0 in {0.735, 0.833, 0.907} (dense, medium, loose Toyoura; DM04's Verdugo &
Ishihara range) x {TC, TCU, PS, PSU}, axial strain to --eps-end.

Gate (soft, plan: 'peak q and eps_v must agree within the triaxial scatter; differences are explained, not tuned
away'). The scatter is K3's measurement model (objective.Weights: sigma_sr 0.10 on sigma1'/sigma3', sigma_ev
0.10 % on eps_v). Drained cases: |d(sr_peak)| <= 0.10 and |d(eps_v at peak)| <= 0.10 % against DM04-paper.
Undrained cases have no Toyoura lab scatter in the project (plan §5.2 K3): reported, not gated.

    python -m harness.k4_crossmodel --theta provisional --out out/k4 [--no-sanisand] [--dist <worktree>/dist/bin]
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import math
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from . import energy as EN
from .drivers import simulate
from .model import Setup
from .sand import TIMS_P_A, toyoura_dm04

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, *[os.pardir] * 4))
SANISAND_REF = os.path.join(ROOT, "Ladruno_scripts")

P0S = (10.0, 25.0, 100.0, 500.0)    # 10 and 25 kPa: the low-p range of the strip deck's near-surface soil (K5)
E0S = (0.735, 0.833, 0.907)
KINDS = ("TC", "TCU", "PS", "PSU")
SIG_SR, SIG_EV_PCT = 0.10, 0.10          # objective.Weights defaults: K3's scatter model
PEAK_TOL = 1.0e-3


# Critical-state convergence (independent of theta): loose undrained compression from both sides of the CSL point
# of its void ratio, to large strain. DM04 and NorSand share the CSL, so both must head to the same p_cs.
CS_E0 = 0.907
CS_P0S = (100.0, 500.0)
CS_EPS = 0.60


def case_list(eps_end, n):
    out = [dict(id=f"{k}_p{p:g}_e{e:g}", kind=k, p0=p, e=e, eps_a_end=eps_end, n=n)
           for k in KINDS for p in P0S for e in E0S]
    n_cs = int(round(n * CS_EPS / eps_end))
    out += [dict(id=f"CS_{k}_p{p:g}_e{CS_E0:g}", kind=k, p0=p, e=CS_E0, eps_a_end=CS_EPS, n=n_cs, cs=True)
            for k in ("TCU", "PSU") for p in CS_P0S]
    return out


def p_cs(sand, e):
    """p on the shared CSL e_c = e0 - lambda_c (p/p_a)^xi at void ratio e (S.22)."""
    return sand.p_a * ((sand.e0 - e) / sand.lambda_c) ** (1.0 / sand.xi)


# ---------------------------------------------------------------------------------------------- common metrics
def metrics(kind, ea_pct, s1, s2, s3, ev_pct):
    """ea, ev in % (compression positive); s1 >= s2 >= s3 principal effective stresses (kPa, compression positive).
    Drained: peak sigma1/sigma3 and the strains there, end values. Undrained: the minimum p (phase-transformation
    proxy) and q at 5 % and at the end."""
    ea, ev = np.asarray(ea_pct), np.asarray(ev_pct)
    s1, s3 = np.asarray(s1), np.asarray(s3)
    sr = s1 / s3
    out = dict(ea_end=float(ea[-1]), n=int(len(ea)))
    if s2 is not None:
        s2 = np.asarray(s2)
        p = (s1 + s2 + s3) / 3.0
        q = np.sqrt(0.5 * ((s1 - s2) ** 2 + (s2 - s3) ** 2 + (s3 - s1) ** 2))
    else:
        p, q = None, None
    t, s = 0.5 * (s1 - s3), 0.5 * (s1 + s3)
    if kind in ("TC", "PS"):
        # the peak POINT is the first one within PEAK_TOL of the maximum: on a plateau (loose sand, low p) the
        # argmax wanders along it (TC p10 e0.907: 2.9 % in one model, 13.9 % in the other, same sr to 0.002)
        i = int(np.argmax(sr >= (1.0 - PEAK_TOL) * sr.max()))
        out.update(sr_peak=float(sr.max()), ea_peak=float(ea[i]), ev_peak=float(ev[i]), ev_end=float(ev[-1]),
                   sr_end=float(sr[-1]), peak_interior=bool(i < len(sr) - 1))
        if q is not None:
            out.update(q_peak=float(q[int(np.argmax(q))]))
    else:
        pp, dev = p, q
        j = int(np.argmin(pp))
        k5 = int(np.argmin(np.abs(ea - 5.0)))
        out.update(p_min=float(pp[j]), ea_at_p_min=float(ea[j]), dev_at_5=float(dev[k5]), dev_end=float(dev[-1]),
                   p_end=float(pp[-1]), sr_end=float(sr[-1]))
    out["t_end"], out["s_end"] = float(t[-1]), float(s[-1])
    return out


# ---------------------------------------------------------------------------------------------- NorSand (O2)
def k3_setup(p_a=TIMS_P_A):
    sand = toyoura_dm04(p_a=p_a)
    return Setup(sand, energy="HAR", policy=EN.ElasticPolicy(mode="per_test"), pi0_rule="unified", p_min=None)


def run_norsand(args):
    case, theta, n = args
    n = int(case.get("n", n))
    setup = k3_setup()
    t0 = time.perf_counter()
    c = simulate(setup, theta, case["kind"], case["p0"], case["e"], case["eps_a_end"], n)
    s = np.sort(-c.sig, axis=1)[:, ::-1]                      # principal, compression positive, s1 >= s2 >= s3
    # metrics only from a COMPLETE curve: a partial one would report its last point as a 'peak'
    m = metrics(case["kind"], c.eps_a_pct, s[:, 0], s[:, 1], s[:, 2], c.eps_v_pct) if c.complete else {}
    return case["id"], dict(status=c.status, complete=c.complete, seconds=time.perf_counter() - t0, metrics=m)


# ---------------------------------------------------------------------------------------------- DM04 (WP-134)
def dm04_params(e, p_a):
    if SANISAND_REF not in sys.path:
        sys.path.insert(0, SANISAND_REF)
    from sanisand_reference.model import TOYOURA
    return dataclasses.replace(TOYOURA, P_atm=p_a, e_init=e)


def run_dm04(args):
    case, variant, p_a = args
    if SANISAND_REF not in sys.path:
        sys.path.insert(0, SANISAND_REF)
    from sanisand_reference.driver import isotropic_state
    from sanisand_reference.integrator import Control, integrate, path_table
    from sanisand_reference.ring import ring_variants
    P = dm04_params(case["e"], p_a)
    O = ring_variants()[variant]
    ea = case["eps_a_end"]
    ctl = {"TC": Control((True, False, False, True, True, True), (ea, 0, 0, 0, 0, 0)),
           "TCU": Control.strain([ea, -0.5 * ea, -0.5 * ea, 0, 0, 0]),
           "PS": Control((True, True, False, True, True, True), (ea, 0, 0, 0, 0, 0)),
           "PSU": Control.strain([ea, 0, -ea, 0, 0, 0])}[case["kind"]]
    t0 = time.perf_counter()
    res = integrate(isotropic_state(case["p0"], case["e"]), ctl, P, O)
    tab = path_table(res, P, O)
    eps = np.array([r["eps"][:3] for r in tab])               # compression positive
    sig = np.array([r["sigma"][:3] for r in tab])
    s = np.sort(sig, axis=1)[:, ::-1]
    done = res.status == "ok" and len(tab) > 1 and abs(eps[-1, 0] - ea) <= 1e-9 * abs(ea)
    m = metrics(case["kind"], 100.0 * eps[:, 0], s[:, 0], s[:, 1], s[:, 2], 100.0 * eps.sum(axis=1)) if done else {}
    return (case["id"], variant), dict(status=res.status, reached_eps_pct=float(100 * eps[-1, 0]) if len(tab) else 0.0,
                                       seconds=time.perf_counter() - t0, metrics=m)


# ---------------------------------------------------------------------------------------------- LadrunoSANISAND (3.12 child)
def run_sanisand(cases, dist, outdir, n):
    import tempfile
    cpath = os.path.join(outdir, "sanisand_cases.json")
    opath = os.path.join(tempfile.mkdtemp(prefix="k4_sanisand_"), "sanisand_raw.json")   # raw curves: not committed
    json.dump(dict(cases=cases), open(cpath, "w", encoding="utf-8"), indent=1)
    p = subprocess.run(["py", "-3.12", "-S", os.path.join(HERE, "k4_sanisand.py"), "--dist", dist,
                        "--scripts", SANISAND_REF, "--cases", cpath, "--out", opath, "--n", str(n)],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    if p.returncode != 0 or not os.path.exists(opath):
        return dict(error=p.stderr[-800:]), {}
    raw = json.load(open(opath, encoding="utf-8"))
    out = {}
    for cid, r in raw["results"].items():
        rows = r.get("rows", [])
        kind = next(c["kind"] for c in cases if c["id"] == cid)
        if r["status"] != "ok" or len(rows) != r.get("n_target", -1) + 1:
            out[cid] = dict(status=r["status"], metrics={})      # incomplete: no metrics (no false 'peak')
            continue
        eps = np.array([x["eps"] for x in rows])               # compression positive
        sig = np.sort(np.array([x["sigma"] for x in rows]), axis=1)[:, ::-1]
        out[cid] = dict(status=r["status"], seconds=r.get("seconds"), solver_evals_max=r.get("secant_iterations_max"),
                        jump_steps=r.get("jump_steps"),
                        metrics=metrics(kind, 100.0 * eps[:, 0], sig[:, 0], sig[:, 1], sig[:, 2], 100.0 * eps.sum(axis=1)))
    return dict(build=raw.get("build"), opts=raw.get("opts")), out


# ---------------------------------------------------------------------------------------------- report
def compare(ns, ref, kind):
    if not ns or not ref:
        return None
    if kind in ("TC", "PS"):
        d_sr = ns["sr_peak"] - ref["sr_peak"]
        d_ev = ns["ev_peak"] - ref["ev_peak"]
        return dict(d_sr_peak=d_sr, d_ev_peak=d_ev, d_ea_peak=ns["ea_peak"] - ref["ea_peak"],
                    within=bool(abs(d_sr) <= SIG_SR and abs(d_ev) <= SIG_EV_PCT))
    return dict(r_p_min=ns["p_min"] / ref["p_min"] if ref["p_min"] else None,
                r_dev_end=ns["dev_end"] / ref["dev_end"] if ref["dev_end"] else None)


def csl_check(p_a):
    sand = toyoura_dm04(p_a=p_a)
    P = dm04_params(0.8, p_a)
    pairs = dict(M=(sand.M, P.Mc), e0=(sand.e0, P.e0), lambda_c=(sand.lambda_c, P.lambda_c), xi=(sand.xi, P.xi),
                 p_a=(sand.p_a, P.P_atm), c=(sand.c_ext, P.c))
    return dict(pairs=pairs, identical=all(abs(a - b) <= 1e-12 for a, b in pairs.values()))


def cs_summary(rep, cases):
    """Critical-state convergence: DM04 and LadrunoSANISAND must END at p_cs (|p/p_cs - 1| <= 1 %); NorSand
    must approach it from both sides (p_end from p0 below p_cs <= p_cs <= p_end from p0 above) with eta -> M
    in TCU (within 1 %). How fast NorSand gets there depends on theta; WHERE it goes does not."""
    sand = toyoura_dm04(p_a=TIMS_P_A)
    pc = p_cs(sand, CS_E0)
    out = dict(e0=CS_E0, p_cs=pc, q_cs=sand.M * pc, rows=[], verdicts=[])
    src = [("NorSand", lambda cid: rep["norsand"].get(cid, {})), ("DM04 paper", lambda cid: rep["dm04"].get(f"{cid}|paper", {})),
           ("DM04 uw_model", lambda cid: rep["dm04"].get(f"{cid}|uw_model", {})),
           ("LadrunoSANISAND", lambda cid: rep["sanisand"].get(cid, {}))]
    ends = {}
    for c in (c for c in cases if c.get("cs")):
        for name, get in src:
            m = get(c["id"]).get("metrics") or {}
            if "p_end" not in m:
                out["verdicts"].append(f"{name} {c['id']}: no complete curve")
                continue
            row = dict(case=c["id"], model=name, p=m["p_end"], q=m["dev_end"], eta=m["dev_end"] / m["p_end"],
                       ratio=m["p_end"] / pc)
            out["rows"].append(row)
            ends[(name, c["kind"], c["p0"])] = row
    for kind in ("TCU", "PSU"):
        for name in ("DM04 paper", "DM04 uw_model", "LadrunoSANISAND"):
            rr = [ends.get((name, kind, p)) for p in CS_P0S]
            if all(rr):
                worst = max(abs(r["ratio"] - 1.0) for r in rr)
                out["verdicts"].append(f"{name} {kind}: ends at p_cs within {100 * worst:.2f} % "
                                       f"({'PASS' if worst <= 0.01 else 'FAIL'}, 1 %)")
        lo, hi = ends.get(("NorSand", kind, min(CS_P0S))), ends.get(("NorSand", kind, max(CS_P0S)))
        if lo and hi:
            brk = lo["p"] <= pc <= hi["p"]
            msg = (f"NorSand {kind}: p_end {lo['p']:.1f} (from {min(CS_P0S):g} kPa) and {hi['p']:.1f} (from "
                   f"{max(CS_P0S):g} kPa) {'bracket' if brk else 'DO NOT bracket'} p_cs {pc:.1f}")
            if kind == "TCU":
                de = max(abs(r["eta"] / sand.M - 1.0) for r in (lo, hi))
                msg += f"; eta within {100 * de:.2f} % of M"
                brk = brk and de <= 0.01
            out["verdicts"].append(msg + f" ({'PASS' if brk else 'FAIL'})")
    return out


def write_md(path, rep):
    L = [f"# K4 cross-model check ({rep['label']})", "",
         f"theta: `{rep['theta']}` ({rep['theta_source']}). {rep['verdict_note']}", "",
         f"CSL shared with DM04: {rep['csl']['identical']} ({', '.join(f'{k} {a:g}/{b:g}' for k, (a, b) in rep['csl']['pairs'].items())}).",
         "", "## Drained (gate: |d sr_peak| <= 0.10, |d ev_peak| <= 0.10 %, NorSand against DM04-paper)", "",
         "| case | NorSand sr_peak / ea / ev | DM04 paper | DM04 uw_model | LadrunoSANISAND | d sr | d ev % | within |",
         "|---|---|---|---|---|---|---|---|"]

    def fmt_d(m):
        return f"{m['sr_peak']:.3f} / {m['ea_peak']:.2f} / {m['ev_peak']:.3f}" if m and "sr_peak" in m else "—"

    def fmt_u(m):
        return f"{m['p_min']:.1f} @ {m['ea_at_p_min']:.2f} ; {m['dev_end']:.1f}" if m and "p_min" in m else "—"

    for c in rep["cases"]:
        cid, kind = c["id"], c["kind"]
        if c.get("cs"):
            continue
        ns = rep["norsand"].get(cid, {}).get("metrics")
        dp = rep["dm04"].get(f"{cid}|paper", {}).get("metrics")
        du = rep["dm04"].get(f"{cid}|uw_model", {}).get("metrics")
        pm = rep["sanisand"].get(cid, {}).get("metrics")
        cmp_ = rep["compare"].get(cid) or {}
        if kind in ("TC", "PS"):
            L.append(f"| {cid} | {fmt_d(ns)} | {fmt_d(dp)} | {fmt_d(du)} | {fmt_d(pm)} | "
                     f"{cmp_.get('d_sr_peak', float('nan')):+.3f} | {cmp_.get('d_ev_peak', float('nan')):+.3f} | "
                     f"{cmp_.get('within', '—')} |")
    L += ["", "## Undrained (reported: min p @ ea % ; q at the end)", "",
          "| case | NorSand | DM04 paper | DM04 uw_model | LadrunoSANISAND | NorSand/DM04 p_min | NorSand/DM04 dev_end |",
          "|---|---|---|---|---|---|---|"]
    for c in rep["cases"]:
        cid, kind = c["id"], c["kind"]
        if kind in ("TC", "PS") or c.get("cs"):
            continue
        ns = rep["norsand"].get(cid, {}).get("metrics")
        dp = rep["dm04"].get(f"{cid}|paper", {}).get("metrics")
        du = rep["dm04"].get(f"{cid}|uw_model", {}).get("metrics")
        pm = rep["sanisand"].get(cid, {}).get("metrics")
        cmp_ = rep["compare"].get(cid) or {}
        f = lambda x: f"{x:.3f}" if isinstance(x, float) else "—"  # noqa: E731
        L.append(f"| {cid} | {fmt_u(ns)} | {fmt_u(dp)} | {fmt_u(du)} | {fmt_u(pm)} | {f(cmp_.get('r_p_min'))} | "
                 f"{f(cmp_.get('r_dev_end'))} |")
    cs = rep.get("cs")
    if cs:
        L += ["", f"## Critical state (independent of theta): loose undrained, e0 {cs['e0']:g}, to {100 * CS_EPS:g} % axial strain",
              "", f"Shared CSL at e {cs['e0']:g}: p_cs = {cs['p_cs']:.1f} kPa (TCU: q_cs = M p_cs = {cs['q_cs']:.1f} kPa). "
              "Undrained, e stays e0, so every model must end at p_cs; PSU's q at critical state depends on each "
              "model's Lode-angle shape, so only its p is checked.", "",
              "| case | model | p_end | q_end | eta_end | p_end / p_cs |", "|---|---|---|---|---|---|"]
        for row in cs["rows"]:
            L.append(f"| {row['case']} | {row['model']} | {row['p']:.1f} | {row['q']:.1f} | {row['eta']:.3f} | {row['ratio']:.4f} |")
        L += [""] + [f"- {v}" for v in cs["verdicts"]]
    cv = [v for v in rep.get("cpp_vs_paper", {}).values() if v and "d_sr_peak" in v]
    if cv:
        L += ["", f"LadrunoSANISAND (C++) against DM04 paper (Python), drained: max |d sr_peak| "
                  f"{max(abs(v['d_sr_peak']) for v in cv):.4f}, max |d ev_peak| {max(abs(v['d_ev_peak']) for v in cv):.4f} %."]
    bad = [f"{k}: {v['status']}" for k, v in {**rep["norsand"], **rep["dm04"], **rep["sanisand"]}.items()
           if v.get("status") not in ("ok",)]
    L += ["", "## Runs that did not finish cleanly", ""] + ([f"- {b}" for b in bad] or ["- none"])
    open(path, "w", encoding="utf-8").write("\n".join(L) + "\n")


def load_theta(spec):
    if spec == "provisional":
        from .smoke_recovery import THETA_TRUE
        return dict(THETA_TRUE), "provisional: smoke_recovery.THETA_TRUE (synthetic), DRY RUN"
    if os.path.exists(spec):
        rec = json.load(open(spec, encoding="utf-8"))          # fit_tatsuoka's record: best.theta + fixed
        th = dict(rec.get("fixed") or {})
        th.update(rec.get("theta") or rec.get("best", {}).get("theta"))
        if rec.get("energy") not in (None, "HAR"):
            raise SystemExit(f"K4 runs K3's HAR setup; {spec} was fitted under {rec.get('energy')}")
        return th, f"K3 result {os.path.basename(spec)} ({rec.get('date', '?')}, cost {rec.get('best', {}).get('cost')})"
    return json.loads(spec), "given on the command line"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--theta", default="provisional", help="'provisional', a K3 result JSON, or a JSON dict")
    ap.add_argument("--out", default=os.path.join(HERE, os.pardir, "out", "k4"))
    ap.add_argument("--eps-end", type=float, default=0.20)
    ap.add_argument("--n", type=int, default=400, help="NorSand steps (0.05 %% at 20 %%)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--no-sanisand", action="store_true", help="skip the LadrunoSANISAND leg (needs dist/bin)")
    ap.add_argument("--dist", default=os.path.join(ROOT, "dist", "bin"))
    ap.add_argument("--only", default=None, help="comma list of case ids (smoke)")
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    theta, theta_src = load_theta(a.theta)
    cases = case_list(a.eps_end, a.n)
    if a.only:
        keep = set(a.only.split(","))
        cases = [c for c in cases if c["id"] in keep]
    rep = dict(label="DRY RUN, provisional theta" if a.theta == "provisional" else "K4", theta=theta,
               theta_source=theta_src, cases=cases, csl=csl_check(TIMS_P_A),
               verdict_note=("Not a K4 verdict: NorSand runs on synthetic parameters until K3 reports."
                             if a.theta == "provisional" else "Soft gate: differences are explained, not tuned away."),
               norsand={}, dm04={}, sanisand={}, compare={})
    t0 = time.perf_counter()
    with ProcessPoolExecutor(a.workers) as ex:
        fn = ex.map(run_norsand, [(c, theta, a.n) for c in cases])
        fd = ex.map(run_dm04, [(c, v, TIMS_P_A) for c in cases for v in ("paper", "uw_model")])
        for cid, r in fn:
            rep["norsand"][cid] = r
        for (cid, v), r in fd:
            rep["dm04"][f"{cid}|{v}"] = r
    if not a.no_sanisand:
        info, rep["sanisand"] = run_sanisand(cases, a.dist, a.out, a.n)
        rep["sanisand_info"] = info
    for c in cases:
        rep["compare"][c["id"]] = compare(rep["norsand"].get(c["id"], {}).get("metrics"),
                                          rep["dm04"].get(f"{c['id']}|paper", {}).get("metrics"), c["kind"])
        rep.setdefault("cpp_vs_paper", {})[c["id"]] = compare(rep["sanisand"].get(c["id"], {}).get("metrics"),
                                                              rep["dm04"].get(f"{c['id']}|paper", {}).get("metrics"), c["kind"])
    rep["cs"] = cs_summary(rep, cases)
    rep["seconds"] = time.perf_counter() - t0
    tag = "dryrun" if a.theta == "provisional" else "k4"
    json.dump(rep, open(os.path.join(a.out, f"k4_{tag}.json"), "w", encoding="utf-8"), indent=1, default=float)
    write_md(os.path.join(a.out, f"k4_{tag}.md"), rep)
    print(f"K4 {tag}: {len(cases)} cases in {rep['seconds']:.0f} s -> {a.out}")


if __name__ == "__main__":
    main()
