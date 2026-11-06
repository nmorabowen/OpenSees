"""Entry point for the K3 fit against Tatsuoka et al. (1986) (plan §3 P3, §5.2 K3; data/README.md).

Inputs (data pack, 2026-10-02): the four PS curves data/tatsuoka1986/tats86_*.csv and the point tests
data/tatsuoka1986/point_tests.csv (Fig. 9 phi_peak(e), Fig. 22 eps_peak(e)); the sand from
data/dm04/toyoura_table1.csv (toyoura_dm04, p_a 101 kPa = the TIMs p_atm, sheet §2.4 / round 3b A4: ONE flag for the
HAR energy and the fork CSL; was 100 kPa, the pack's assumption, before 2026-10-03). e is e_0.05 for every test (data/README.md §3.1):
the 49 kPa and >= 98 kPa states are slightly denser than modelled; stated, not corrected.
Defaults: energy BA06 ('per_test'; --energy HAR for the TIMs/Toyoura energy, owner decision (a)), WW, smooth cap,
pi_i0 rule 'unified' (S.53, the oracle's helper; model.Setup), the oracle's default p' floor, rho pinned at
DM04's c = 0.712 (sheet §15 'direct': M_e = rho M_c; there is no Toyoura TC/TE pair to pin it from data, data/README.md
§6 item 9), free (chi, h, N, N_bar, rho_bar); eps_min 0.1 % (seating, data/README.md §3.3).

  python -m harness.fit_tatsuoka --eval chi=-3,h=150,N=0.3,N_bar=0.2,rho_bar=0.75   # one evaluation, breakdown
  python -m harness.fit_tatsuoka --starts 16 --workers 16 [--energy HAR] [--free chi,h,N,N_bar,rho_bar,rho]
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys
import time

import numpy as np

from . import data as DATA
from . import energy as EN
from . import fit as FIT
from .model import Setup, FIT_NAMES
from .objective import Objective, Weights
from .sand import toyoura_dm04, DATA_DIR, TIMS_P_A

OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))


LIT_LAM = ("lt88_fig5a_w00_hw1p0", "lt88_fig9_psc_w00_hw1p0", "lt88_fig9_psc_w00_hw1p9")


def load_literature(sub_weight=0.5, low_kpa=12.0):
    """Round-3b K3 additions: Lam & Tatsuoka 1988 (omega = 0 only: TC H/W 1.0, PSC H/W 1.0 and 1.9; the H/W 0.25 TC and
    the omega 30/60/90 curves are left out: end-restraint geometry / cross-anisotropy the model does not have) and the
    12 Fukushima & Tatsuoka 1984 Fig. 5 TC curves. Fukushima sigma3 = sigma3_f_kPa (the measured effective confining
    stress at failure; the header's sigma_c' is the nominal) where measured, else sigma_c'. Their phi is uncorrected for
    membrane force (biased low at low sigma3'): weight `sub_weight` where sigma3 < `low_kpa`. Post-peak stays a gate
    (Weights.post_peak_weight applies to all curves)."""
    out = []
    for stem in LIT_LAM:
        out.append(DATA.load_curve(os.path.join(DATA_DIR, "lam_tatsuoka1988", "curves", stem + ".csv")))
    for p in sorted(glob.glob(os.path.join(DATA_DIR, "fukushima1984", "curves", "ft84_fig5*.csv"))):
        lc = DATA.load_curve(p)
        try:
            lc.meta["sigma3_kPa"] = float(lc.meta["sigma3_f_kPa"])
        except (KeyError, ValueError):
            pass
        if lc.sigma3 < low_kpa:
            lc.meta["weight"] = sub_weight
        out.append(lc)
    return out


def build(args):
    sand = toyoura_dm04(p_a=args.p_a)
    policy = EN.ElasticPolicy(mode=args.policy, p_rep=args.p_rep, e_rep=args.e_rep)
    setup = Setup(sand, energy=args.energy, policy=policy, pi0_rule=args.pi0_rule, pi0_ratio=args.pi0_ratio,
                  p_min=args.p_min)
    curves = [DATA.load_curve(p) for p in sorted(glob.glob(os.path.join(DATA_DIR, "tatsuoka1986", "tats86_*.csv")))]
    if args.lit:
        curves += load_literature()
    points = [] if args.no_points else DATA.load_points(os.path.join(DATA_DIR, "tatsuoka1986", "point_tests.csv"),
                                                         source="Tatsuoka 1986 Figs. 9 / 22")
    W = Weights(post_peak_pct=args.post_peak, post_peak_weight=args.post_weight, eps_min_pct=args.eps_min)
    return sand, setup, Objective(setup, curves, points, W)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--energy", default="BA06")
    ap.add_argument("--policy", default="per_test")
    ap.add_argument("--p-rep", type=float, default=None)
    ap.add_argument("--e-rep", type=float, default=None, help="global policy: the one void ratio of the BVP body")
    ap.add_argument("--p-a", type=float, default=TIMS_P_A, help="the shared -p_a (HAR energy and fork CSL), kPa")
    ap.add_argument("--p-min", type=float, default=None, help="p' floor, kPa (default: the oracle's 5e-3 p_ref; 0 = off)")
    ap.add_argument("--pi0-rule", default="unified")
    ap.add_argument("--pi0-ratio", type=float, default=None)
    ap.add_argument("--rho", type=float, default=None, help="pinned rho (default: the sand's c)")
    ap.add_argument("--free", default="chi,h,N,N_bar,rho_bar")
    ap.add_argument("--starts", type=int, default=8)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--post-peak", type=float, default=2.0)
    ap.add_argument("--post-weight", type=float, default=0.5)
    ap.add_argument("--eps-min", type=float, default=0.1)
    ap.add_argument("--no-points", action="store_true")
    ap.add_argument("--profiles", action="store_true")
    ap.add_argument("--eval", default=None, help="name=value,... : evaluate once and print the breakdown")
    ap.add_argument("--lit", action="store_true", help="add Lam & Tatsuoka 1988 (omega 0) and Fukushima & Tatsuoka 1984 curves")
    ap.add_argument("--out", default=None, help="output directory (default ../out)")
    ap.add_argument("--tag", default="")
    ap.add_argument("--checkpoint", default=None,
                    help="directory: each start is saved as start_<k>.json when it finishes; existing ones are skipped")
    ap.add_argument("--only", default=None,
                    help="start indices to run, e.g. '0,1,2,3' or '4-7' (Slurm array slices of one Sobol set); "
                         "needs --checkpoint, writes no final record")
    ap.add_argument("--collect", action="store_true",
                    help="no new starts: read --checkpoint, then identifiability/profiles and the final record")
    a = ap.parse_args(argv)
    if (a.only or a.collect) and not a.checkpoint:
        ap.error("--only and --collect need --checkpoint")
    sand, setup, obj = build(a)
    rho = a.rho if a.rho is not None else sand.c_ext
    free = tuple(n for n in a.free.split(",") if n)
    fixed = {n: rho for n in ("rho",) if n not in free}
    print(f"sand: {sand.name}; energy {a.energy} ({EN.get(a.energy).available('O2')[1]}); "
          f"{len(obj.curves)} curves, {len(obj.points)} point tests; free {free}; fixed {fixed}")
    if a.eval:
        th = dict(fixed)
        for kv in a.eval.split(","):
            k, v = kv.split("=")
            th[k.strip()] = float(v)
        t0 = time.perf_counter()
        bd = obj.breakdown(th)
        print(json.dumps(bd, indent=1, default=float))
        print(f"one breakdown: {time.perf_counter() - t0:.1f} s")
        return 0
    spec = FIT.FitSpec(free=free, fixed=fixed)
    t0 = time.perf_counter()
    only = None
    if a.only:
        only = []
        for part in a.only.split(","):
            lo, _, hi = part.partition("-")
            only += list(range(int(lo), int(hi) + 1)) if hi else [int(lo)]
    if not a.collect:
        runs = FIT.multistart(obj, spec, n_starts=a.starts, seed=a.seed, workers=a.workers, only=only,
                              checkpoint_dir=a.checkpoint)
        if a.only:
            for r in sorted(runs, key=lambda r: r["start"]):
                print(f"start {r['start']}: cost {r['cost']:.6g}, {r['residual_evals']} evals, {r['seconds']:.0f} s")
            print(f"slice {a.only} done -> {a.checkpoint}")
            return 0
    if a.checkpoint:
        runs = FIT.load_checkpoints(a.checkpoint)
        missing = sorted(set(range(a.starts)) - {int(r["start"]) for r in runs})
        print(f"collected {len(runs)} starts from {a.checkpoint}; missing {missing if missing else 'none'}")
        if not runs:
            return 1
    best = runs[0]
    rec = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), sand=sand.as_dict(), energy=a.energy, policy=vars(setup.policy),
               setup=dict(zeta=setup.zeta, cap=setup.cap, c1=setup.c1, c2=setup.c2, pi0_rule=setup.pi0_rule),
               free=list(spec.free), fixed=fixed, weights=obj.W.as_dict(), runs=runs, best=best,
               breakdown=obj.breakdown(best["theta"]), wall_seconds=time.perf_counter() - t0)
    if a.checkpoint:
        rec["checkpoint"] = dict(dir=a.checkpoint, starts_planned=a.starts, missing=missing)
    rec["identifiability"] = FIT.identifiability(obj, best["theta"], spec.free)
    if a.profiles:
        rec["profiles"] = {nm: FIT.profile(obj, spec, best["theta"], nm, workers=a.workers) for nm in ("chi", "h")}
    out_dir = a.out or OUT
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"fit_tatsuoka_{a.energy}{('_' + a.tag) if a.tag else ''}.json")
    with open(path, "w") as f:
        json.dump(rec, f, indent=1, default=str)
    print(f"best cost {best['cost']:.4g}: {best['theta']}  -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
