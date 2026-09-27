"""WP-128 Q4: is ModifiedEuler's substep count ACCURACY-limited (then a looser
norm / error floor buys substeps) or STABILITY-limited (then it cannot)?

Take committed states from the constant-p chains (q4_baseline's constp, d_eps_yy
1e-5, eta/M^b near its peak) and re-integrate the NEXT increment on the
validated port while varying TolE (1e-6 .. 1e-1) and the proposed floor
s_ref (1 = today, 20, 100, 1e4 kPa).  An accuracy-limited explicit scheme
scales substeps ~ TolE^-1/2; a stability-limited one does not move.
Also: the same on the ring rows (probes @ 1e-5), median over rows.
Output: out/q4_stiffness.txt
"""
import statistics as stx

import _boot as B
from _boot import sr
import drive as D
import md_port

LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def main():
    B.define_prototypes()
    tols = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]
    floors = [1.0, 20.0, 100.0, 1e4]
    for p0 in (2.0, 20.0):
        h = D.run_const_p(D.k0_state(p0), 1e-5, 600, B.TAG_ME)
        x = max(h, key=lambda x: x["eta"] / x["Mb"])
        st, de = x["before"], x["de"]
        say(f"== constant-p p0={p0}, increment k={x['k']} (eta/M^b {x['eta'] / x['Mb']:.3f}), d_eps {['%.3g' % v for v in de[:2]]}")
        say("   substeps by TolE (rows) x s_ref (cols); today = TolE 1e-4, s_ref 1")
        say(f"   {'TolE':>7} " + " ".join(f"{'s_ref=' + format(f, 'g'):>11}" for f in floors))
        for t in tols:
            row = []
            for f in floors:
                o = md_port.Material(B.P, TolE=t, err_floor=f).update(
                    st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], de)
                row.append(f"{o['substeps']:>6}{'*' if o['forced'] else ' '}{'R' if o['rc'] else ' '}")
            say(f"   {t:>7.0e} " + " ".join(f"{r:>11}" for r in row))
    say("   (* = a forced accept at dT_min occurred; R = refused)")
    # ring rows, probes @ 1e-5: median substeps by TolE
    rows = [r for p in sr.RING_CSVS for r in sr.read_ring_csv(p)]
    d = 1e-5
    probes = [[d, d, 0, 0, 0, 0], [-d, -d, 0, 0, 0, 0], [0, 0, 0, d, 0, 0], [0, 0, 0, -d, 0, 0]]
    say("== ring rows x 4 probes @ 1e-5: substeps median / p95 / max by TolE (s_ref = today)")
    for t in tols:
        m = md_port.Material(B.P, TolE=t)
        subs = []
        for r in rows:
            for de in probes:
                o = m.update(r["sigma"], B.dev(r["alpha"]), B.dev(r["alpha_in"]), B.dev(r["z"]), r["e"], de)
                subs.append(o["substeps"])
        subs.sort()
        say(f"   TolE {t:>7.0e}: {stx.median(subs):>6.0f} / {subs[int(0.95 * (len(subs) - 1))]:>6} / {subs[-1]:>6}")
    with open(f"{B.OUT}/q4_stiffness.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
