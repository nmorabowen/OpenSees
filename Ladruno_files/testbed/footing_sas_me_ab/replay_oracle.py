"""WP-138: integrate the footing's dumped Gauss-point states through WP-134's
independent SANISAND reference (`sanisand_reference`, Radau IIA at rtol 1e-10,
event-segmented), preset `uw_model` (UW constitutive additions U1-U5, paper
alpha_in rule, continuous moduli -- the oracle a corrected C++ integrator
should reproduce).

Needs numpy + scipy (CPython 3.11 here); the package is found via --ref-dir
(the directory CONTAINING `sanisand_reference/`, e.g. a checkout of
origin/wp/134-sanisand-reference-integrator's Ladruno_scripts).

    python replay_oracle.py --csv <replay.csv> --out <json> --ref-dir <dir>
"""
import argparse
import csv
import json
import math
import sys
import time


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref-dir", required=True)
    ap.add_argument("--preset", default="uw_model")
    ap.add_argument("--max-rows", type=int, default=0)
    a = ap.parse_args()
    sys.path.insert(0, a.ref_dir)
    from sanisand_reference import CAMPAIGN, State, integrate
    from sanisand_reference.ring import ring_variants
    O = ring_variants()[a.preset]
    rows = list(csv.DictReader(open(a.csv, newline="")))
    if a.max_rows:
        rows = rows[:a.max_rows]
    out = []
    for r in rows:
        v = lambda k: [float(r[f"{k}_{i}"]) for i in range(6)]
        st = State.from_voigt(v("sigma"), v("alpha"), v("z"), float(r["e"]), v("alpha_in"))
        t0 = time.time()
        try:
            res = integrate(st, v("dStrain"), CAMPAIGN, O)
            s = res.summary()
            s.pop("notes", None)
            d = dict(element=int(r["element"]), gp=int(r["gp"]), step=int(r["step"]),
                     select=r["select"], wall=time.time() - t0, **s,
                     start_rho_alpha=res.start.get("rho_alpha"),
                     min_H_sign_margin=res.min_H_sign_margin)
        except Exception as exc:     # an oracle failure is a finding, not a crash
            d = dict(element=int(r["element"]), gp=int(r["gp"]), step=int(r["step"]),
                     select=r["select"], status=f"exception: {exc!r}",
                     wall=time.time() - t0)
        out.append(d)
        print(f"{d['element']}/{d['gp']} {d['status']} t_end {d.get('t_end')} "
              f"f_end {d.get('f_end')} rho_a_end {d.get('rho_alpha_end')} "
              f"({d['wall']:.1f}s)", flush=True)
    with open(a.out, "w") as f:
        json.dump(dict(preset=a.preset, csv=a.csv, rows=out), f, indent=1, default=float)


if __name__ == "__main__":
    main()
