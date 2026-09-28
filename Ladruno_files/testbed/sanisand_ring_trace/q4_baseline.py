"""WP-128 Q4: the F18(a) measurement baseline.

For ONE strain increment from a given committed state, compare the campaign
integrator (IntScheme 1, TolE hard-coded 1e-4, today's norm) with a reference
and count substeps.  The same increment is also run through the validated
Python port (md_port) with the proposed norm err/max(2||sigma||, s_ref) for
s_ref in {0, 0.1, 1, 5, 20} kPa -- the numbers F18(a) asked for, offline.

REFERENCE.  RK45 (IntScheme 45) is NOT usable: its dT_min is hard-coded 1e-3
and at dT_min it force-accepts with the Mc clamp (q3: it "returns" eta 1.33
from eta 12.87 on a 1e-7 increment).  The reference is ModifiedEuler with
-honorTolR 1 at TolR 1e-8 (tag 11) and 1e-9 (tag 12), -maxSubsteps 1e6; a
case whose reference force-accepted, capped, or disagrees with the 1e-9 run by
more than 10 % of the campaign error is reported as "no reference".

Sets:
  ring   the 80 attached states x {isoComp, isoExt, shear+, shear-} x {1e-6, 1e-5}
  chain  single points from K0 at p0 = 2, 5, 20 kPa driven along the plane-
         strain `active` path (lateral extension; q1_search: it reaches
         alpha/alpha^b 0.88-0.99) at delta 1e-5 until eta/M^b(theta) peaks;
         every committed state of that chain is re-integrated for its NEXT
         increment.  (p RISES along it: 2 -> 12.5 kPa at the peak.)
  constp the same p0 values held CONSTANT (d_eps_yy = 1e-5 compression,
         d_eps_xx solved per increment by secant for p = p0), until
         eta/M^b(theta) peaks -- the low-confinement version of the ask.
Output: out/q4_baseline.json, out/q4_baseline.txt
"""
import json
import math
import statistics as stx

import _boot as B
from _boot import ops, sr
import drive as D
import md_port

TAG_R8, TAG_R9 = 11, 12
FLOORS = [0.0, 0.1, 1.0, 5.0, 20.0]
LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def define():
    B.define_prototypes()
    ops.nDMaterial("LadrunoSANISAND", TAG_R8, *B.P, *B._opts(1, 1.0e-8, 1, 1000000))
    ops.nDMaterial("LadrunoSANISAND", TAG_R9, *B.P, *B._opts(1, 1.0e-9, 1, 1000000))


def cnorm(v):
    return math.sqrt(sum(x * x for x in v[:3]) + 2 * sum(x * x for x in v[3:]))


def one(st, de, ports):
    rep = lambda tag: sr.replay(ops, tag, st["sigma"], st["alpha"], st["alpha_in"],
                                st["z"], st["e"], de, "compressionPositive", trace=0)
    c, r8, r9 = rep(B.TAG_ME), rep(TAG_R8), rep(TAG_R9)
    s8 = r8["stats"]
    ref_ok = (r8["rc"] == 0 and int(s8["forcedAtDTmin"]) == 0 and int(s8["capHits"]) == 0
              and int(s8["abandonedLowP"]) == 0)
    err_c = cnorm([a - b for a, b in zip(c["sigma"], r8["sigma"])])
    err_ref = cnorm([a - b for a, b in zip(r8["sigma"], r9["sigma"])])
    if err_ref > 0.1 * max(err_c, 1e-12):
        ref_ok = False
    rec = dict(ref_ok=ref_ok, err=err_c, err_ref=err_ref, sub=int(c["stats"]["substeps"]),
               forced=int(c["stats"]["forcedAtDTmin"]), rc=c["rc"],
               sub_ref=int(s8["substeps"]), snorm=cnorm(st["sigma"]),
               p=B.tr(st["sigma"]) / 3.0, f_after=c["f_after"], floors={})
    al, ai, z = B.dev(st["alpha"]), B.dev(st["alpha_in"]), B.dev(st["z"])
    for fl, m in ports.items():
        o = m.update(st["sigma"], al, ai, z, st["e"], de)
        rec["floors"][fl] = dict(sub=o["substeps"], forced=o["forced"], rc=o["rc"],
                                 err=cnorm([a - b for a, b in zip(o["sigma"], r8["sigma"])]))
    return rec


def summarize(tag, recs):
    ok = [r for r in recs if r["ref_ok"] and r["rc"] == 0]
    say(f"-- {tag}: {len(recs)} increments, reference available on {len(ok)} "
        f"(excluded: ref forced/capped/unconverged {len(recs) - len(ok)})")
    if not ok:
        return {}
    out = {}

    def stats(xs):
        xs = sorted(xs)
        return dict(median=stx.median(xs), p95=xs[int(0.95 * (len(xs) - 1))], max=xs[-1])

    row = dict(sub=stats([r["sub"] for r in ok]), err=stats([r["err"] for r in ok]),
               rel=stats([r["err"] / max(r["snorm"], 1e-12) for r in ok]),
               forced=sum(r["forced"] for r in ok))
    out["today"] = row
    say(f"   today (C++)  substeps med/p95/max {row['sub']['median']:.0f}/{row['sub']['p95']:.0f}/{row['sub']['max']:.0f}"
        f"   |dsigma| kPa med/p95/max {row['err']['median']:.2e}/{row['err']['p95']:.2e}/{row['err']['max']:.2e}"
        f"   rel {row['rel']['median']:.1e}/{row['rel']['p95']:.1e}/{row['rel']['max']:.1e}  forced {row['forced']}")
    for fl in FLOORS:
        rs = [r["floors"][fl] for r in ok]
        row = dict(sub=stats([x["sub"] for x in rs]), err=stats([x["err"] for x in rs]),
                   forced=sum(x["forced"] for x in rs), refused=sum(1 for x in rs if x["rc"] != 0))
        out[f"floor_{fl}"] = row
        say(f"   s_ref={fl:>5} (port) substeps med/p95/max {row['sub']['median']:.0f}/{row['sub']['p95']:.0f}/{row['sub']['max']:.0f}"
            f"   |dsigma| kPa med/p95/max {row['err']['median']:.2e}/{row['err']['p95']:.2e}/{row['err']['max']:.2e}"
            f"  forced {row['forced']} refused {row['refused']}")
    return out


def main():
    define()
    ports = {fl: md_port.Material(B.P, err_floor=fl) for fl in FLOORS}
    ports_today = md_port.Material(B.P)
    result = {}
    # ---------------- ring states
    ring = []
    same_as_today = 0
    ncmp = 0
    for path in sr.RING_CSVS:
        for r in sr.read_ring_csv(path):
            st = dict(sigma=r["sigma"], alpha=r["alpha"], alpha_in=r["alpha_in"], z=r["z"], e=r["e"])
            for d in (1e-6, 1e-5):
                for pn, de in (("isoComp", [d, d, 0, 0, 0, 0]), ("isoExt", [-d, -d, 0, 0, 0, 0]),
                               ("shear+", [0, 0, 0, d, 0, 0]), ("shear-", [0, 0, 0, -d, 0, 0])):
                    rec = one(st, de, ports)
                    rec.update(row=f"{r['element']}/{r['gp']}", probe=f"{pn}@{d:.0e}",
                               eta_over_Mb=r["eta_over_Mb_compression"])
                    ring.append(rec)
                    # the norm identity: s_ref = 1 kPa IS today's norm
                    o1 = ports[1.0].update(r["sigma"], B.dev(r["alpha"]), B.dev(r["alpha_in"]), B.dev(r["z"]), r["e"], de)
                    o0 = ports_today.update(r["sigma"], B.dev(r["alpha"]), B.dev(r["alpha_in"]), B.dev(r["z"]), r["e"], de)
                    ncmp += 1
                    same_as_today += (o1["sigma"] == o0["sigma"] and o1["substeps"] == o0["substeps"])
    say(f"norm identity: port with s_ref = 1 kPa reproduces today's norm bit-for-bit on "
        f"{same_as_today}/{ncmp} ring increments")
    result["norm_identity"] = [same_as_today, ncmp]
    result["ring"] = summarize("ring states (80 x 8 probes)", ring)
    for d in ("1e-06", "1e-05"):
        result[f"ring_{d}"] = summarize(f"ring states, probes @ {d}", [x for x in ring if x["probe"].endswith(d)])
    # ---------------- chains
    for p0 in (2.0, 5.0, 20.0):
        st = D.k0_state(p0)
        de = [-1e-5, 1e-5, 0.0, 0.0, 0.0, 0.0]
        hist = D.run(st, [de] * 400, "cpp", tag=B.TAG_ME)
        etab = [x["eta"] / x["Mb"] for x in hist]
        kpk = max(range(len(etab)), key=lambda k: etab[k])
        say(f"== chain p0={p0}: active path, delta 1e-5; eta/M^b(theta) peaks at {etab[kpk]:.3f} "
            f"(k={kpk}, p={hist[kpk]['p']:.3g} kPa)")
        recs = []
        for x in hist[:kpk + 1]:
            rec = one(x["before"], de, ports)
            rec.update(k=x["k"], eta_over_Mb=x["eta"] / x["Mb"])
            recs.append(rec)
        result[f"chain_{p0}"] = summarize(f"chain p0={p0} (all increments to the peak)", recs)
        near = [r for r in recs if r["eta_over_Mb"] >= 0.9 * etab[kpk]]
        result[f"chain_{p0}_near_peak"] = summarize(
            f"chain p0={p0}, increments with eta/M^b >= 0.9 x peak", near)
    # ---------------- constant-p chains: eta/M^b -> 1 AT p0 (the active path lets p rise)
    for p0 in (2.0, 5.0, 20.0):
        st = D.k0_state(p0)
        hist = D.run_const_p(st, 1e-5, 600, B.TAG_ME)
        etab = [x["eta"] / x["Mb"] for x in hist]
        kpk = max(range(len(etab)), key=lambda k: etab[k])
        say(f"== constant-p chain p0={p0}: d_eps_yy 1e-5, d_eps_xx solved for p = p0 "
            f"(max |p miss| {max(abs(x['p_miss']) for x in hist):.1e} kPa); eta/M^b(theta) "
            f"peaks at {etab[kpk]:.3f} (k={kpk}, p={hist[kpk]['p']:.4g} kPa)")
        recs = []
        for x in hist[:kpk + 1]:
            rec = one(x["before"], x["de"], ports)
            rec.update(k=x["k"], eta_over_Mb=x["eta"] / x["Mb"])
            recs.append(rec)
        result[f"constp_{p0}"] = summarize(f"constant-p p0={p0} (all increments to the peak)", recs)
        near = [r for r in recs if r["eta_over_Mb"] >= 0.9 * etab[kpk]]
        result[f"constp_{p0}_near_peak"] = summarize(
            f"constant-p p0={p0}, increments with eta/M^b >= 0.9 x peak", near)
    with open(f"{B.OUT}/q4_baseline.json", "w") as fh:
        json.dump(dict(result=result, ring=ring), fh, indent=1, default=str)
    with open(f"{B.OUT}/q4_baseline.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
