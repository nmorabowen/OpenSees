"""WP-138: q-s curve match between two legs (the C-vs-B gate) and cost ratio.

    python compare_curves.py REF_LEG LEG [--tol-kpa 0.28]

q of LEG is compared with REF_LEG's at LEG's own s/B points inside the common
range (REF linearly interpolated). The tolerance default 0.28 kPa is the
ladder's loosest acceptance (KrylovNewton at 10 x 1e-5 x the 4152 kN/m applied
load = 0.415 kN/m residual norm) spread over B = 1.5 m -- a crude bound on how
far two converged states may differ in the footing reaction; the relative
difference is reported too. Cost: cumulative push wall, substeps and Newton
iterations at the common end.
"""
import argparse
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def load(leg):
    return list(csv.DictReader(open(os.path.join(HERE, "runs", leg, "steps.csv"), newline="")))


def cum(rows, key, upto):
    return sum(float(r[key] or 0) for r in rows if float(r["s_over_B"]) <= upto + 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("leg")
    ap.add_argument("--tol-kpa", type=float, default=0.28)
    a = ap.parse_args()
    R, L = load(a.ref), load(a.leg)
    rs = [float(r["s_over_B"]) for r in R]
    rq = [float(r["q_kPa"]) for r in R]
    smax = min(rs[-1], float(L[-1]["s_over_B"]))
    worst = (0.0, 0.0, 0.0)
    worst_rel = (0.0, 0.0, 0.0)
    nbad = 0
    for r in L:
        s, q = float(r["s_over_B"]), float(r["q_kPa"])
        if s > smax:
            break
        # interpolate REF
        j = next(i for i, x in enumerate(rs) if x >= s - 1e-15)
        if j == 0:
            qr = rq[0]
        else:
            t = (s - rs[j - 1]) / (rs[j] - rs[j - 1])
            qr = rq[j - 1] + t * (rq[j] - rq[j - 1])
        d = q - qr
        if abs(d) > a.tol_kpa:
            nbad += 1
        if abs(d) > abs(worst[1]):
            worst = (s, d, d / qr)
        if abs(d / qr) > abs(worst_rel[2]):
            worst_rel = (s, d, d / qr)
    print(f"{a.leg} vs {a.ref}: common s/B <= {smax:.5f}; worst |dq| {worst[1]:+.3f} kPa "
          f"({100*worst[2]:+.3f} %) at s/B {worst[0]:.5f}; worst relative "
          f"{100*worst_rel[2]:+.3f} % ({worst_rel[1]:+.3f} kPa) at s/B {worst_rel[0]:.6f}; "
          f"points beyond {a.tol_kpa} kPa: {nbad}")
    # relative band by s/B range
    for lo, hi in ((0, 0.001), (0.001, 0.005), (0.005, 1.0)):
        rel = []
        for r in L:
            s, q = float(r["s_over_B"]), float(r["q_kPa"])
            if not (lo <= s < hi) or s > smax:
                continue
            j = next(i for i, x in enumerate(rs) if x >= s - 1e-15)
            qr = rq[0] if j == 0 else rq[j - 1] + (s - rs[j - 1]) / (rs[j] - rs[j - 1]) * (rq[j] - rq[j - 1])
            rel.append(abs(q - qr) / qr)
        if rel:
            print(f"  s/B in [{lo}, {min(hi, smax):.4f}): max |dq|/q = {100*max(rel):.3f} % over {len(rel)} points")
    for lab, rows in ((a.ref, R), (a.leg, L)):
        print(f"  {lab:22s} to s/B {smax:.5f}: push wall {cum(rows, 'wall_step_s', smax)/3600:.2f} h, "
              f"substeps {cum(rows, 'sub_step_total', smax):.3e}, iterations "
              f"{cum(rows, 'iters', smax):.0f}, steps "
              f"{sum(1 for r in rows if float(r['s_over_B']) <= smax + 1e-12)}, "
              f"failed attempts {cum(rows, 'fails_before', smax):.0f}, refusals/caps "
              f"{cum(rows, 'cap_step', smax):.0f}")


if __name__ == "__main__":
    main()
