"""Entry point for the K3 fit against Tatsuoka et al. (1986) (plan §3 P3, §5.2 K3; data/README.md).

Inputs (data pack, 2026-10-02): the four PS curves data/tatsuoka1986/tats86_*.csv and the point tests
data/tatsuoka1986/point_tests.csv (Fig. 9 phi_peak(e), Fig. 22 eps_peak(e)); the sand from
data/dm04/toyoura_table1.csv (toyoura_dm04, p_a assumed 100 kPa). e is e_0.05 for every test (data/README.md §3.1):
the 49 kPa and >= 98 kPa states are slightly denser than modelled; stated, not corrected.
Defaults: energy BA06 ('per_test'), WW, smooth cap, pi_i0 rule 'ramp_end' (model.Setup), rho pinned at
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
from .sand import toyoura_dm04, DATA_DIR

OUT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "out"))


def build(args):
    sand = toyoura_dm04(p_a=args.p_a)
    policy = EN.ElasticPolicy(mode=args.policy, p_rep=args.p_rep)
    setup = Setup(sand, energy=args.energy, policy=policy, pi0_rule=args.pi0_rule, pi0_ratio=args.pi0_ratio)
    curves = [DATA.load_curve(p) for p in sorted(glob.glob(os.path.join(DATA_DIR, "tatsuoka1986", "tats86_*.csv")))]
    points = [] if args.no_points else DATA.load_points(os.path.join(DATA_DIR, "tatsuoka1986", "point_tests.csv"),
                                                         source="Tatsuoka 1986 Figs. 9 / 22")
    W = Weights(post_peak_pct=args.post_peak, post_peak_weight=args.post_weight, eps_min_pct=args.eps_min)
    return sand, setup, Objective(setup, curves, points, W)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--energy", default="BA06")
    ap.add_argument("--policy", default="per_test")
    ap.add_argument("--p-rep", type=float, default=None)
    ap.add_argument("--p-a", type=float, default=100.0)
    ap.add_argument("--pi0-rule", default="ramp_end")
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
    ap.add_argument("--tag", default="")
    a = ap.parse_args(argv)
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
    runs = FIT.multistart(obj, spec, n_starts=a.starts, seed=a.seed, workers=a.workers)
    best = runs[0]
    rec = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), sand=sand.as_dict(), energy=a.energy, policy=vars(setup.policy),
               setup=dict(zeta=setup.zeta, cap=setup.cap, c1=setup.c1, c2=setup.c2, pi0_rule=setup.pi0_rule),
               free=list(spec.free), fixed=fixed, weights=obj.W.as_dict(), runs=runs, best=best,
               breakdown=obj.breakdown(best["theta"]), wall_seconds=time.perf_counter() - t0)
    rec["identifiability"] = FIT.identifiability(obj, best["theta"], spec.free)
    if a.profiles:
        rec["profiles"] = {nm: FIT.profile(obj, spec, best["theta"], nm, workers=a.workers) for nm in ("chi", "h")}
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, f"fit_tatsuoka_{a.energy}{('_' + a.tag) if a.tag else ''}.json")
    with open(path, "w") as f:
        json.dump(rec, f, indent=1, default=str)
    print(f"best cost {best['cost']:.4g}: {best['theta']}  -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
