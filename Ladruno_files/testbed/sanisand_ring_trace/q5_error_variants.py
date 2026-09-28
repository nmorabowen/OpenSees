"""WP-128 Q5 / finding E: does an alpha (and fabric) aware substep error
suffice ALONE, or must it be combined with refusing at dT_min (and a bound on
alpha)?  All on the validated port (md_port), one switch set per variant:

  today        C++-faithful ModifiedEuler
  +a           substep error = max(stress, alpha)       (RK45's form)
  +a+z         ... and the fabric z
  +a+z,refuse  ... and a substep that fails at dT_min REFUSES the update
               (rc != 0) instead of the forced accept + Mc clamp

Sets: (a) the q1 reversal chains (vertUnload / extShear, delta 1e-4, p0 2 kPa);
(b) one increment from every committed state of the constant-p chains at
p0 = 2 and 20 kPa (C++ chain, eta/M^b to its peak) -- the cost on benign
states; (c) the 80 ring rows x 4 probes x {1e-6, 1e-5}.
Output: out/q5_error_variants.txt / .json
"""
import json
import math
import statistics as stx

import _boot as B
from _boot import sr
import drive as D
import md_port

LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def variants():
    return {
        "today": dict(),
        "+a": dict(alpha_err=True),
        "+a+z": dict(alpha_err=True, fabric_err=True),
        "+a+z,refuse": dict(alpha_err=True, fabric_err=True, forced_policy="refuse"),
    }


def cyc(d, delta, n):
    out = []
    for _ in range(3):
        out += [[delta * x for x in d]] * n
        out += [[-delta * x for x in d]] * (n // 2)
    return out


def pct(xs, q):
    xs = sorted(xs)
    return xs[int(q * (len(xs) - 1))]


def main():
    B.define_prototypes()
    out = {}
    say("== (a) reversal chains, p0 = 2 kPa, delta 1e-4, (20 fwd, 10 back) x 3")
    say(f"{'path':<12} {'variant':<12} {'max a/ab':>9} {'#out':>5} {'substeps':>9} {'forced':>6} {'refused':>7} {'f>1e-6':>6}")
    for pname, d in (("vertUnload", [0, -1.0, 0, 0, 0, 0]), ("extShear", [0.3, -1.0, 0, 0.5, 0, 0])):
        incs = cyc(d, 1e-4, 20)
        for vn, kw in variants().items():
            h = D.run(D.k0_state(2.0), incs, "port", mat=md_port.Material(B.P, **kw))
            r = [x["alpha_over_b"] for x in h if math.isfinite(x["alpha_over_b"])]
            row = dict(max=max(r), n_out=sum(1 for x in r if x > 1), substeps=sum(x["substeps"] for x in h),
                       forced=sum(x["forced"] for x in h), refused=sum(1 for x in h if x["rc"] != 0),
                       f_pos=sum(1 for x in h if x["rc"] == 0 and x["f"] > 1e-6))
            out[f"chain|{pname}|{vn}"] = row
            say(f"{pname:<12} {vn:<12} {row['max']:>9.3f} {row['n_out']:>5} {row['substeps']:>9} "
                f"{row['forced']:>6} {row['refused']:>7} {row['f_pos']:>6}")
    say("")
    say("== (b) one increment from each committed state of the constant-p chains (C++), to the eta/M^b peak")
    say(f"{'p0':>5} {'variant':<12} {'sub med':>8} {'sub max':>8} {'max |dsig| vs today kPa':>24} {'forced':>6} {'refused':>7}")
    for p0 in (2.0, 20.0):
        h = D.run_const_p(D.k0_state(p0), 1e-5, 600, B.TAG_ME)
        kpk = max(range(len(h)), key=lambda k: h[k]["eta"] / h[k]["Mb"])
        base = None
        for vn, kw in variants().items():
            m = md_port.Material(B.P, **kw)
            subs, sig, forced, refused = [], [], 0, 0
            for x in h[:kpk + 1]:
                st = x["before"]
                o = m.update(st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], x["de"])
                subs.append(o["substeps"])
                sig.append(o["sigma"])
                forced += o["forced"]
                refused += (o["rc"] != 0)
            if base is None:
                base = sig
            dmax = max(math.sqrt(sum((a - b) ** 2 for a, b in zip(s1, s0))) for s1, s0 in zip(sig, base))
            out[f"constp|{p0}|{vn}"] = dict(sub_med=stx.median(subs), sub_max=max(subs), dmax=dmax,
                                            forced=forced, refused=refused)
            say(f"{p0:>5} {vn:<12} {stx.median(subs):>8.0f} {max(subs):>8} {dmax:>24.2e} {forced:>6} {refused:>7}")
    say("")
    say("== (c) ring rows x 4 probes x {1e-6, 1e-5}")
    say(f"{'variant':<12} {'sub med/p95/max':>18} {'in->out':>7} {'forced':>6} {'refused':>7} {'f>1e-6':>6}")
    rows = [r for p in sr.RING_CSVS for r in sr.read_ring_csv(p)]
    probes = []
    for dl in (1e-6, 1e-5):
        probes += [[dl, dl, 0, 0, 0, 0], [-dl, -dl, 0, 0, 0, 0], [0, 0, 0, dl, 0, 0], [0, 0, 0, -dl, 0, 0]]
    for vn, kw in variants().items():
        m = md_port.Material(B.P, **kw)
        subs, io, forced, refused, fpos = [], 0, 0, 0, 0
        for r in rows:
            al, ai, z = B.dev(r["alpha"]), B.dev(r["alpha_in"]), B.dev(r["z"])
            a0 = B.bounding(r["sigma"], al, r["e"])["alpha_over_b"]
            for de in probes:
                o = m.update(r["sigma"], al, ai, z, r["e"], de)
                subs.append(o["substeps"])
                forced += o["forced"]
                refused += (o["rc"] != 0)
                if o["rc"] == 0:
                    a1 = B.bounding(o["sigma"], o["alpha"], o["e"])["alpha_over_b"] if B.tr(o["sigma"]) > 0 else float("nan")
                    io += (a0 <= 1.0 < a1)
                    fpos += (o["f_after"] > 1e-6)
        out[f"ring|{vn}"] = dict(sub_med=stx.median(subs), sub_p95=pct(subs, 0.95), sub_max=max(subs),
                                 in_to_out=io, forced=forced, refused=refused, f_pos=fpos)
        say(f"{vn:<12} {stx.median(subs):>6.0f}/{pct(subs, 0.95):>5}/{max(subs):>5} {io:>7} {forced:>6} {refused:>7} {fpos:>6}")
    with open(f"{B.OUT}/q5_error_variants.json", "w") as fh:
        json.dump(out, fh, indent=1)
    with open(f"{B.OUT}/q5_error_variants.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
