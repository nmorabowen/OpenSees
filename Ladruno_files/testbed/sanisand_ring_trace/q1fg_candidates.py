"""WP-128: the orchestrator's candidates F and G for finding B, tested on the
validated port (md_port) alongside E.

F  a loading substep classified ELASTIC because the plastic-multiplier
   denominator temp4 = Kp + 2G(B - C tr n^3) - K D n:r is NEGATIVE (dgamma < 0
   although the numerator n:Ce:deps > 0), then the dgamma<0 branch drags alpha
   by the stress-ratio change; err = 0 -> uncapped q.
G  alpha_in re-seated once per global increment, so (alpha - alpha_in):n can
   cross zero / go negative inside the substeps and h = b0/((alpha-alpha_in):n)
   blows up or flips sign.

1. The stage-by-stage signs (temp4, numerator, Kp, (alpha-alpha_in):n) on the
   reproducer and on 1950/3's f>0 probes.
2. Every accepted substep of the two reversal chains that carries alpha/alpha^b
   across 1: which stage kinds and signs it had.
3. Counterfactuals on the chains + reproducer, one switch each:
   G-fix  h = b0/<(alpha-alpha_in):n>  (Macaulay: 1e10 when <= 0)
   F-q    accepted-substep growth capped q <= 2
   F-drag the dgamma<0 stage leaves alpha alone
   E-fix  alpha in the substep error
Output: out/q1fg_candidates.txt
"""
import math
import _boot as B
from _boot import sr
import drive as D
import md_port

LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def sgn(x):
    return "+" if x > 0 else ("-" if x < 0 else "0")


def fmt_sl(sl):
    return " | ".join(f"{k}: temp4 {sgn(t)}{abs(t):.2g} num {sgn(nm)} Kp {sgn(kp)} aain {sgn(a)}{abs(a):.2g}"
                      for (t, nm, kp, a, k) in sl)


def ab(S, A, e):
    return B.bounding(S, A, e)["alpha_over_b"] if B.tr(S) > 0 else float("nan")


def cyc(d, delta, n):
    out = []
    for _ in range(3):
        out += [[delta * x for x in d]] * n
        out += [[-delta * x for x in d]] * (n // 2)
    return out


def open_update(m, st, de, label, prev=0.0, nshow=4):
    o = m.update(st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], de, prev_incr_norm=prev)
    say(f"-- {label}: path {o['path']}, substeps {o['substeps']}, alpha_in reset {o['alpha_in'] != list(st['alpha_in'])}, "
        f"f_after {o['f_after']:.3g}, alpha/alpha_b {ab(o['sigma'], o['alpha'], o['e']):.3f}, corrGiveUp {o['corrGiveUp']}")
    for i, t in enumerate(o["trace"][:nshow]):
        if t[5] is None:
            continue
        say(f"   sub {i} T={t[0]:.3g} dT={t[1]:.2g} code {t[3]} err {t[2]:.2g}: {fmt_sl(t[5]['sl'])}")
    return o


def main():
    B.define_prototypes()
    m = md_port.Material(B.P)
    say("== 1. stage signs")
    floor = dict(sigma=[0.0101] * 3 + [0.0] * 3, alpha=[0.0] * 6, alpha_in=[0.0] * 6, z=[0.0] * 6, e=0.697787979641054)
    open_update(m, floor, [0, 1e-4, 0, 0, 0, 0], "reproducer: floor, deps_yy 1e-4")
    rows = {(r["element"], r["gp"]): r for r in sr.read_ring_csv(sr.RING_CSVS[0])}
    r = rows[(1950, 3)]
    st = dict(sigma=r["sigma"], alpha=B.dev(r["alpha"]), alpha_in=B.dev(r["alpha_in"]), z=B.dev(r["z"]), e=r["e"])
    st0 = m.state_dep(st["sigma"], st["alpha"], st["z"], st["e"], st["alpha_in"])
    p = B.tr(st["sigma"]) / 3
    say(f"   1950/3 COMMITTED state: Kp = {2/3*p*st0['h']*md_port.contr(st0['b'], st0['n']):.3g}, h {st0['h']:.3g}, "
        f"b:n {md_port.contr(st0['b'], st0['n']):.3g}, (alpha-alpha_in):n {st0['aain']:.3g}")
    for lab, de in (("1950/3 shear+ 1e-5", [0, 0, 0, 1e-5, 0, 0]), ("1950/3 isoExt 1e-6", [-1e-6, -1e-6, 0, 0, 0, 0]),
                    ("1950/3 isoComp 1e-6 (loading)", [1e-6, 1e-6, 0, 0, 0, 0])):
        open_update(m, st, de, lab)

    say("")
    say("== 2. accepted substeps that carry alpha/alpha_b across 1 (reversal chains, delta 1e-4, p0 2)")
    for pname, d in (("vertUnload", [0, -1.0, 0, 0, 0, 0]), ("extShear", [0.3, -1.0, 0, 0.5, 0, 0])):
        incs = cyc(d, 1e-4, 20)
        stc, prev = D.k0_state(2.0), 0.0
        tally = {}
        for k, de in enumerate(incs):
            o = m.update(stc["sigma"], stc["alpha"], stc["alpha_in"], stc["z"], stc["e"], de, prev_incr_norm=prev)
            reset = o["alpha_in"] != list(stc["alpha_in"])
            for t in o["trace"]:
                info = t[5]
                if info is None or t[3] != 0:
                    continue
                a0, a1 = ab(info["S"], info["A"], info["e"]), ab(info["nS"], info["nA"], info["e"])
                if a0 <= 1.0 < a1:
                    sl = info["sl"]
                    key = " + ".join(
                        f"{kd}(temp4{sgn(t4)},num{sgn(nm)},aain{'=0' if abs(a) < 1e-10 else sgn(a)})"
                        for (t4, nm, kp, a, kd) in sl) + f"; alpha_in reset this incr: {reset}; dT {info['dT']:.2g}"
                    tally[key] = tally.get(key, 0) + 1
            if o["rc"] == 0:
                stc = dict(sigma=o["sigma"], alpha=o["alpha"], alpha_in=o["alpha_in"], z=o["z"], e=o["e"])
                prev = D.ncov(de)
        say(f"   {pname}:")
        for kk, v in tally.items():
            say(f"      {v:>3} x  {kk}")

    say("")
    say("== 3. counterfactuals: max alpha/alpha_b on the chains, and the reproducer")
    cfs = {"as built": {}, "G-fix h Macaulay": dict(h_mode="macaulay"), "F-q q<=2": dict(q_max=2.0),
           "F-drag frozen": dict(drag="frozen"), "E-fix alpha in error": dict(alpha_err=True),
           "G-fix + F-drag": dict(h_mode="macaulay", drag="frozen")}
    say(f"   {'variant':<22} {'vertUnload':>11} {'extShear':>9} {'reproducer':>11} {'repro substeps':>14}")
    for name, kw in cfs.items():
        mm = md_port.Material(B.P, **kw)
        res = []
        for d in ([0, -1.0, 0, 0, 0, 0], [0.3, -1.0, 0, 0.5, 0, 0]):
            h = D.run(D.k0_state(2.0), cyc(d, 1e-4, 20), "port", mat=mm)
            res.append(max(x["alpha_over_b"] for x in h if math.isfinite(x["alpha_over_b"])))
        o = mm.update(floor["sigma"], floor["alpha"], floor["alpha_in"], floor["z"], floor["e"], [0, 1e-4, 0, 0, 0, 0])
        say(f"   {name:<22} {res[0]:>11.3f} {res[1]:>9.3f} {ab(o['sigma'], o['alpha'], o['e']):>11.3f} {o['substeps']:>14}")
    with open(f"{B.OUT}/q1fg_candidates.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
