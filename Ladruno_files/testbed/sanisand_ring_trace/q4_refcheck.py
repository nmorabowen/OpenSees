"""WP-128 Q4 caveat check: q4_baseline's reference (C++ ModifiedEuler at TolR
1e-8) is itself ALPHA-BLIND (finding E) -- where both Heun stages are elastic
in stress it accepts in one substep at any tolerance.  How far is it from an
alpha-aware reference?  Reference B = the port with alpha AND fabric in the
substep error at TolE 1e-8 (dT_min 1e-6, forced accepts counted).
Ring rows x 4 probes @ 1e-6 (the 1e-5 set is too slow in Python at 1e-8).
Output: out/q4_refcheck.txt
"""
import math
import statistics as stx

import _boot as B
from _boot import ops, sr
import md_port

TAG_R8 = 11


def cn(v):
    return math.sqrt(sum(x * x for x in v[:3]) + 2 * sum(x * x for x in v[3:]))


def main():
    B.define_prototypes()
    ops.nDMaterial("LadrunoSANISAND", TAG_R8, *B.P, *B._opts(1, 1.0e-8, 1, 1000000))
    mref = md_port.Material(B.P, TolE=1e-8, alpha_err=True, fabric_err=True, maxSubsteps=1000000)
    rows = [r for p in sr.RING_CSVS for r in sr.read_ring_csv(p)]
    d = 1e-6
    probes = [[d, d, 0, 0, 0, 0], [-d, -d, 0, 0, 0, 0], [0, 0, 0, d, 0, 0], [0, 0, 0, -d, 0, 0]]
    gap_ref, err_today_B, err_today_A = [], [], []
    nforced = 0
    for r in rows:
        al, ai, z = B.dev(r["alpha"]), B.dev(r["alpha_in"]), B.dev(r["z"])
        for de in probes:
            A = sr.replay(ops, TAG_R8, r["sigma"], r["alpha"], r["alpha_in"], r["z"], r["e"], de,
                          "compressionPositive", trace=0)
            T = sr.replay(ops, B.TAG_ME, r["sigma"], r["alpha"], r["alpha_in"], r["z"], r["e"], de,
                          "compressionPositive", trace=0)
            o = mref.update(r["sigma"], al, ai, z, r["e"], de)
            if o["forced"] or o["rc"] != 0 or A["stats"]["forcedAtDTmin"] > 0:
                nforced += 1
                continue
            gap_ref.append(cn([a - b for a, b in zip(A["sigma"], o["sigma"])]))
            err_today_B.append(cn([a - b for a, b in zip(T["sigma"], o["sigma"])]))
            err_today_A.append(cn([a - b for a, b in zip(T["sigma"], A["sigma"])]))
    lines = []
    for name, xs in (("|ref A (C++ ME 1e-8) - ref B (alpha+z-aware 1e-8)|", gap_ref),
                     ("|today - ref A|", err_today_A), ("|today - ref B|", err_today_B)):
        xs = sorted(xs)
        lines.append(f"{name:<52} kPa  median {stx.median(xs):.2e}  p95 {xs[int(0.95 * (len(xs) - 1))]:.2e}  max {xs[-1]:.2e}")
    lines.append(f"cases {len(gap_ref)} (excluded {nforced}: either reference forced/refused)")
    txt = "\n".join(lines)
    print(txt)
    with open(f"{B.OUT}/q4_refcheck.txt", "w") as fh:
        fh.write(txt + "\n")


if __name__ == "__main__":
    main()
