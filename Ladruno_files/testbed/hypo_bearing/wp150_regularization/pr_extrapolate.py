"""WP-150: de-bias footing legs run with -Presidual p_r > 0 (needed to keep the free surface from refusing).

    python pr_extrapolate.py LABEL=RUNDIR:PR [LABEL=RUNDIR:PR ...] [--sb 0.005,0.01,...]

p_r shifts the Mohr-Coulomb envelope, tau = (sigma + p_r) tan(phi), i.e. an apparent cohesion c = p_r tan(phi). The
classical capacity then gains c N_c = p_r (N_q - 1). So, at the limit state, q is ~linear in p_r with a slope N_q - 1.
Hence:
  * at every s/B reached by >= 2 legs: least-squares q = q0 + b p_r  ->  q0 = the p_r -> 0 extrapolation;
  * the slope b gives an INDEPENDENT operative friction angle, phi_b: N_q(phi_b) - 1 = b (exact N_q), to compare with
    the phi' implied by q0 through the classical band (t6_capacity_bands.py);
  * the peak (a maximum of q followed by >= 1 % drop, if any) is extrapolated the same way.
Before a mechanism forms, q is not at a limit state and the slope is only an apparent one; quote b and phi_b at or near
the peak.
"""
import csv, math, os, sys
import numpy as np


def curve(d):
    rows = list(csv.DictReader(open(os.path.join(d, "steps.csv"))))
    return np.array([float(r["s_over_B"]) for r in rows]), np.array([float(r["q_kPa"]) for r in rows])


def n_q(phi):
    t = math.tan(math.radians(phi))
    return math.exp(math.pi * t) * math.tan(math.radians(45 + phi / 2)) ** 2


def phi_from_slope(b):
    if not b > 0:
        return float("nan")
    lo, hi = 1.0, 60.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if n_q(mid) - 1 < b else (lo, mid)
    return 0.5 * (lo + hi)


def peak(sb, q):
    i = int(np.argmax(q))
    return (sb[i], q[i]) if i < len(q) - 1 and q[i:].min() < 0.99 * q[i] else (float("nan"), float("nan"))


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    grid = [0.002, 0.005, 0.01, 0.015, 0.02, 0.03, 0.05, 0.075, 0.1]
    if "--sb" in sys.argv:
        grid = [float(v) for v in sys.argv[sys.argv.index("--sb") + 1].split(",")]
        args = [a for a in args if a != sys.argv[sys.argv.index("--sb") + 1]]
    legs = []
    for a in args:
        lab, rest = a.split("=", 1)
        d, pr = rest.rsplit(":", 1)
        sb, q = curve(d)
        legs.append((lab, float(pr), sb, q))
    legs.sort(key=lambda t: t[1])
    print("| leg | p_r kPa | reached s/B | peak s/B / q_peak |")
    print("|---|---|---|---|")
    for lab, pr, sb, q in legs:
        ps, pq = peak(sb, q)
        print(f"| {lab} | {pr:g} | {sb[-1]:.4f} | {ps:.4f} / {pq:.1f} |")
    print("\n| s/B | legs | q at each p_r (kPa) | q0 (p_r → 0) | slope b = dq/dp_r | φ_b from N_q − 1 = b |")
    print("|---|---|---|---|---|---|")
    for s in grid:
        pts = [(pr, float(np.interp(s, sb, q))) for _, pr, sb, q in legs if s <= sb[-1]]
        if len(pts) < 2:
            continue
        P = np.array([p for p, _ in pts]); Q = np.array([v for _, v in pts])
        b, q0 = np.polyfit(P, Q, 1)
        print(f"| {s:g} | {len(pts)} | " + " / ".join(f"{v:.1f}@{p:g}" for p, v in pts)
              + f" | {q0:.1f} | {b:.1f} | {phi_from_slope(b):.1f}° |")
    pk = [(pr, peak(sb, q)[1]) for _, pr, sb, q in legs]
    pk = [(p, v) for p, v in pk if v == v]
    if len(pk) >= 2:
        P = np.array([p for p, _ in pk]); Q = np.array([v for _, v in pk])
        b, q0 = np.polyfit(P, Q, 1)
        print(f"\nPEAK: q_peak(p_r → 0) = {q0:.1f} kPa from {len(pk)} legs; slope {b:.1f} → φ_b {phi_from_slope(b):.1f}°")
    else:
        print("\nPEAK: fewer than two legs show a peak yet")


if __name__ == "__main__":
    main()
