"""WP-129: the SAS-ME evidence tables (printed; out/sasme_report.txt + .json).

  T1   WP-128's reproducer grid, ME vs SAS-ME (reseat / bracket / ablated)
  T2   the vertUnload / extShear reversal chains (p0 2 kPa, delta 1e-4)
  T3   the 80 ring rows x 8 probes: refusals by code, f at exit, alpha/alpha^b,
       substeps against today's ModifiedEuler
  CP   constant-direction smooth chains at p0 2 / 5 / 20 kPa (cost)
  PROF the OPS_PROFILE_SCOPE split of SAS-ME at a ring state and a deep state

Run with the fork's test bootstrap (CPython 3.12 -S, this worktree's
dist/bin + tests/ + Ladruno_scripts on sys.path).
"""
import json
import math
import os
import statistics as stat
import sys
import tempfile
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, os.path.join(ROOT, "tests"))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

from _testbed import ops  # noqa: E402
import sanisand_replay as sr  # noqa: E402
import wp129_sasme_tools as W  # noqa: E402

OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
LINES = []
DATA = {}


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def sub(o):
    return int(o["sas"]["substeps"]) if o["sas"] and o["sas"]["updates"] > 0 else int(o["stats"]["substeps"])


def pct(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * len(xs)))] if xs else float("nan")


def t1():
    say("== T1: sigma = p_s I, alpha = alpha_in = z = 0, d_eps_yy = delta (one increment)")
    say(f"{'p_s':>7} {'delta':>7} | " + " | ".join(f"{n:^30}" for n in ("ME", "SAS reseat", "SAS bracket", "SAS ablated")))
    say(f"{'':>17}" + " | ".join(f"{'rc':>3} {'sub':>6} {'a/ab_n':>7} {'f':>9}" for _ in range(4)))
    rows = []
    for ps in (0.0101, 0.1, 1.0, 5.0):
        for delta in (1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4):
            st = dict(sigma=[ps, ps, ps, 0, 0, 0], alpha=[0.0] * 6, alpha_in=[0.0] * 6,
                      z=[0.0] * 6, e=0.697787979641054)
            cells = []
            for tag in (W.TAG_ME, W.TAG_SAS, W.TAG_SAS_BR, W.TAG_SAS_ABL):
                new, o = W.step(ops, tag, st, [0, delta, 0, 0, 0, 0])
                cells.append(dict(rc=o["rc"], sub=sub(o),
                                  ab_n=W.alpha_over_b_n(new["sigma"], new["alpha"], new["e"]),
                                  f=o["f_after"]))
            rows.append(dict(ps=ps, delta=delta, cells=cells))
            say(f"{ps:>7.4g} {delta:>7.0e} | " + " | ".join(
                f"{c['rc']:>3} {c['sub']:>6} {c['ab_n']:>7.3f} {c['f']:>9.1e}" for c in cells))
    DATA["T1"] = rows


def t2():
    say("\n== T2: WP-128 reversal chains, K0 p0 = 2 kPa, delta 1e-4, (20 fwd, 10 back) x 3")
    res = {}
    for pname, d in W.PATHS.items():
        for label, tag in (("ME", W.TAG_ME), ("SAS reseat", W.TAG_SAS), ("SAS bracket", W.TAG_SAS_BR)):
            t0 = time.perf_counter()
            h = W.run_chain(ops, tag, W.k0_state(2.0), W.incs_for(d, 1e-4, 20))
            dt = time.perf_counter() - t0
            ok = [x for x in h if x["rc"] == 0]
            subs = sum((int(x["sas"]["substeps"]) if x["sas"] and x["sas"]["updates"] > 0
                        else x["substeps"]) for x in h)
            codes = {}
            for x in h:
                if x["rc"] != 0:
                    c = sr.SAS_REFUSE_CODES.get(int(x["sas"]["lastRefuseCode"]), "?") if x["sas"] else "ME"
                    codes[c] = codes.get(c, 0) + 1
            r = dict(max_ab_n=max(x["ab_n"] for x in ok), max_ab=max(x["ab"] for x in ok),
                     f_gt=sum(1 for x in ok if x["f"] > 1e-6), refused=len(h) - len(ok),
                     codes=codes, substeps=subs, clamp=sum(x["clamp"] for x in h), wall=dt)
            res[f"{pname}/{label}"] = r
            say(f"{pname:>10} {label:<11}: max a/ab_n {r['max_ab_n']:.3f} (a/ab alpha-Lode {r['max_ab']:.3f}), "
                f"f>1e-6 ok {r['f_gt']}, refused {r['refused']} {codes}, Mc-clamps {r['clamp']}, "
                f"substeps {subs}, wall {dt:.1f}s")
    DATA["T2"] = res


def t3():
    say("\n== T3: 80 ring rows x 8 probes (+-iso, +-shear at 1e-6, 1e-5)")
    tab = []
    for path in sr.RING_CSVS:
        for r in sr.read_ring_csv(path):
            for delta in (1e-6, 1e-5):
                for pname, de in W.ring_probes(delta).items():
                    nc, oc = W.step(ops, W.TAG_SAS, r, de)
                    nm, om = W.step(ops, W.TAG_ME, r, de)
                    tab.append(dict(set=os.path.basename(path)[12:15], el=r["element"], gp=r["gp"],
                                    probe=pname, delta=delta, rc=oc["rc"], sub=sub(oc),
                                    code=int(oc["sas"]["lastRefuseCode"]), f=oc["f_after"],
                                    ab=W.alpha_over_b(nc["sigma"], nc["alpha"], nc["e"]),
                                    ab_n=W.alpha_over_b_n(nc["sigma"], nc["alpha"], nc["e"]),
                                    me_rc=om["rc"], me_sub=sub(om), me_f=om["f_after"],
                                    me_ab_n=W.alpha_over_b_n(nm["sigma"], nm["alpha"], nm["e"]),
                                    me_clamp=int(om["stats"]["forcedClampMc"]),
                                    me_forced=int(om["stats"]["forcedAtDTmin"])))
    ok = [t for t in tab if t["rc"] == 0]
    codes = {}
    for t in tab:
        if t["rc"] != 0:
            k = sr.SAS_REFUSE_CODES.get(t["code"], "?")
            codes[k] = codes.get(k, 0) + 1
    say(f"SAS-ME: {len(ok)}/{len(tab)} integrated, refused by code {codes}")
    say(f"  f at exit (rc 0): max {max(t['f'] for t in ok):.2e}; count f > 1e-6: {sum(t['f'] > 1e-6 for t in ok)}")
    say(f"  alpha/alpha^b at exit (rc 0): max (alpha Lode) {max(t['ab'] for t in ok):.3f}, "
        f"max (n Lode, WP-128 metric) {max(t['ab_n'] for t in ok):.3f}")
    ss = [t["sub"] for t in ok]
    ms = [t["me_sub"] for t in tab]
    say(f"  substeps SAS (rc 0) med/p95/max {pct(ss, .5)}/{pct(ss, .95)}/{max(ss)}; "
        f"ME med/p95/max {pct(ms, .5)}/{pct(ms, .95)}/{max(ms)}; totals SAS {sum(t['sub'] for t in tab)} vs ME {sum(ms)}")
    say(f"  ME on the same 640: f > 1e-6 as success {sum(1 for t in tab if t['me_rc'] == 0 and t['me_f'] > 1e-6)}, "
        f"Mc-clamp teleports {sum(1 for t in tab if t['me_clamp'] > 0)}, forced@dTmin {sum(1 for t in tab if t['me_forced'] > 0)}")
    say("  refused rows (set el/gp probe delta code):")
    for t in tab:
        if t["rc"] != 0:
            say(f"    {t['set']} {t['el']}/{t['gp']} {t['probe']:<7} {t['delta']:.0e} "
                f"{sr.SAS_REFUSE_CODES.get(t['code'], '?'):<26} (ME: rc {t['me_rc']}, sub {t['me_sub']}, "
                f"f {t['me_f']:.1e}, clamp {t['me_clamp']})")
    DATA["T3"] = tab


def cp():
    say("\n== CP: smooth monotonic chains from K0 (30 x 1e-5 'active' = (0.3, 1, 0, 0, 0, 0))")
    res = {}
    for p0 in (2.0, 5.0, 20.0):
        for label, tag in (("ME", W.TAG_ME), ("SAS", W.TAG_SAS)):
            h = W.run_chain(ops, tag, W.k0_state(p0), [[3e-6, 1e-5, 0, 0, 0, 0]] * 30)
            subs = [int(x["sas"]["substeps"]) if x["sas"] and x["sas"]["updates"] > 0 else x["substeps"] for x in h]
            res[f"{p0}/{label}"] = dict(med=pct(subs, .5), max=max(subs), refused=sum(x["rc"] != 0 for x in h),
                                        end_ab_n=h[-1]["ab_n"])
            say(f"  p0 {p0:>4} {label:<3}: substeps/increment med {pct(subs, .5)} max {max(subs)}, "
                f"refused {res[f'{p0}/{label}']['refused']}, end a/ab_n {h[-1]['ab_n']:.3f}")
    DATA["CP"] = res


def prof():
    import h5py
    say("\n== PROF: OPS_PROFILE_SCOPE split (coarse profiler), 200 replays each")
    ring = next(r for r in sr.read_ring_csv(sr.RING_CSVS[1]))
    deep = W.k0_state(200.0)
    st = deep
    for _ in range(5):   # put the deep state on the cone
        st, _o = W.step(ops, W.TAG_SAS, st, [3e-5, 1e-4, 0, 0, 0, 0])
    cases = {"ring (b16 row 1)": (ring, [0, 0, 0, 1e-5, 0, 0]),
             "deep (p ~ 200 kPa, on the cone)": (st, [3e-6, 1e-5, 0, 0, 0, 0])}
    res = {}
    for name, (s0, de) in cases.items():
        ops.profiler("reset")
        ops.profiler("start")
        for _ in range(200):
            W.step(ops, W.TAG_SAS, s0, de)
        ops.profiler("stop")
        fn = os.path.join(tempfile.gettempdir(), "wp129_prof.h5")
        ops.profiler("report", fn, "-run", "sas")
        ops.profiler("reset")
        times = {}
        with h5py.File(fn, "r") as f:
            def visit(nm, obj):
                if "sanisand.sasME" in nm and isinstance(obj, h5py.Group):
                    key = [p for p in nm.split("/") if p.startswith("sanisand.sasME")][-1]
                    for a in ("total_ns", "self_ns", "calls", "total_s"):
                        if a in obj.attrs:
                            times.setdefault(key, {})[a] = float(obj.attrs[a])
                    for dname in obj:
                        if isinstance(obj[dname], h5py.Dataset) and obj[dname].shape in ((), (1,)):
                            try:
                                times.setdefault(key, {})[dname] = float(obj[dname][()])
                            except Exception:
                                pass
            f.visititems(visit)
        res[name] = times
        say(f"  {name}:")
        for k, v in sorted(times.items()):
            say(f"    {k:<36} {v}")
    DATA["PROF"] = res


def main():
    W.define_prototypes(ops)
    which = sys.argv[1:] or ["t1", "t2", "t3", "cp", "prof"]
    for w in which:
        globals()[w]()
    with open(os.path.join(OUT, "sasme_report_" + "_".join(which) + ".txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(LINES) + "\n")
    with open(os.path.join(OUT, "sasme_report_" + "_".join(which) + ".json"), "w") as fh:
        json.dump(DATA, fh, default=str)


if __name__ == "__main__":
    main()
