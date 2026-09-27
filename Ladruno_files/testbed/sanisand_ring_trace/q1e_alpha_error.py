"""WP-128 finding E (orchestrator): ModifiedEuler's substep error is
STRESS-ONLY (ManzariDafalias.cpp ModifiedEuler: ||dSigma2 - dSigma1||);
RungeKutta45 (IntScheme **45**, not 4 -- the #define is INT_RungeKutta45 45)
uses max(stress error, alpha error).  Does an alpha-aware error keep alpha
inside the bounding surface where ME lets it escape?

Integrators compared, no C++ change:
  ME        IntScheme 1, campaign (TolE 1e-4 hard-coded, honorTolR 0)
  RK45/1e-4 IntScheme 45, TolR 1e-4  (-honorTolR is ME-only; RK45 always
            reads TolR).  NOTE its dT_min is hard-coded 1e-3 and at dT_min it
            FORCE-ACCEPTS with the same Mc clamp + alpha re-derivation as ME.
  RK45/1e-7 IntScheme 45, TolR 1e-7 (the campaign's TolR)
  ME+aErr   the validated port with alpha ALSO in the substep error
            (RK45's form: absolute below ||alpha|| 0.5, /2||alpha|| above),
            dT_min 1e-6, TolE 1e-4 -- ME with only that one change.
RK45 is not instrumented by WP-127 (substepStats covers ME only), so its cost
is reported as wall time and its forced accepts only through the clamp
signature eta == Mc (the radial clamp sets sqrt(3/2)||s||/p = Mc exactly).

Sets: the q1_attrib reversal chains (vertUnload, extShear, delta 1e-4,
(20 fwd, 10 back) x 3, p0 = 2 kPa); the same vertUnload at delta 1e-5; the
q1_threshold single increments from the floor; the 80 ring rows x 8 probes.
Output: out/q1e_alpha_error.txt / .json
"""
import json
import math
import time

import _boot as B
from _boot import ops, sr
import drive as D
import md_port

TAG_RK4, TAG_RK7 = 5, 6
LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def cyc(d, delta, n):
    out = []
    for _ in range(3):
        out += [[delta * x for x in d]] * n
        out += [[-delta * x for x in d]] * (n // 2)
    return out


def clamp_sig(x):
    return math.isfinite(x["eta"]) and abs(x["eta"] - B.MC) < 1e-9


def main():
    B.define_prototypes()
    ops.nDMaterial("LadrunoSANISAND", TAG_RK4, *B.P, *B._opts(45, 1.0e-4, 1))
    ops.nDMaterial("LadrunoSANISAND", TAG_RK7, *B.P, *B._opts(45, 1.0e-7, 1))
    integ = [("ME", "cpp", B.TAG_ME), ("RK45/1e-4", "cpp", TAG_RK4),
             ("RK45/1e-7", "cpp", TAG_RK7), ("ME+aErr", "port", None)]
    res = {}
    paths = {
        "vertUnload 1e-4": cyc([0, -1.0, 0, 0, 0, 0], 1e-4, 20),
        "extShear 1e-4": cyc([0.3, -1.0, 0, 0.5, 0, 0], 1e-4, 20),
        "shear 1e-4": cyc([0, 0, 0, 1.0, 0, 0], 1e-4, 20),
        "vertUnload 1e-5": cyc([0, -1.0, 0, 0, 0, 0], 1e-5, 200),
    }
    say("== reversal chains from K0 at p0 = 2 kPa")
    say(f"{'path':<17} {'integrator':<10} {'max a/ab':>9} {'first>1':>7} {'#>1':>5} {'ME subst':>9} "
        f"{'clamp-sig':>9} {'f>1e-6 commits':>14} {'wall s':>7}")
    for pname, incs in paths.items():
        st0 = D.k0_state(2.0)
        for name, be, tag in integ:
            t0 = time.perf_counter()
            if be == "cpp":
                h = D.run(st0, incs, "cpp", tag=tag)
            else:
                h = D.run(st0, incs, "port", mat=md_port.Material(B.P, alpha_err=True))
            wall = time.perf_counter() - t0
            r = [x["alpha_over_b"] for x in h if math.isfinite(x["alpha_over_b"])]
            first = next((x["k"] for x in h if x["alpha_over_b"] > 1.0), None)
            nout = sum(1 for x in r if x > 1.0)
            ncl = sum(1 for x in h if clamp_sig(x))
            nf = sum(1 for x in h if x["rc"] == 0 and x["f"] > 1e-6)
            sub = sum(x["substeps"] for x in h) if name.startswith("ME") else float("nan")
            say(f"{pname:<17} {name:<10} {max(r):>9.3f} {str(first):>7} {nout:>5} {sub:>9.0f} "
                f"{ncl:>9} {nf:>14} {wall:>7.2f}")
            res[f"{pname}|{name}"] = dict(max=max(r), first=first, n_out=nout, substeps=sub,
                                         clamp_sig=ncl, f_pos=nf, wall=wall)
    # ---- single increments from the floor (q1_threshold's start state)
    say("")
    say("== one increment from sigma = p_s I, alpha = alpha_in = 0 (d_eps_yy = delta): alpha/alpha_b after")
    say(f"{'p_s':>7} {'delta':>7} " + " ".join(f"{n:>11}" for n, _, _ in integ))
    m_ae = md_port.Material(B.P, alpha_err=True)
    thr = {}
    for ps in (0.0101, 0.1, 1.0):
        for delta in (1e-5, 3e-5, 1e-4, 3e-4):
            st = dict(sigma=[ps, ps, ps, 0, 0, 0], alpha=[0.0] * 6, alpha_in=[0.0] * 6,
                      z=[0.0] * 6, e=0.697787979641054)
            de = [0.0, delta, 0.0, 0.0, 0.0, 0.0]
            row = []
            for name, be, tag in integ:
                if be == "cpp":
                    new, info = D.step_cpp(tag, st, de, 0.0)
                    d = D.diag(new)
                    v = (d["alpha_over_b"], clamp_sig(d))
                else:
                    o = m_ae.update(st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], de)
                    v = (B.bounding(o["sigma"], o["alpha"], o["e"])["alpha_over_b"], False)
                row.append(v)
            thr[f"{ps}|{delta}"] = row
            say(f"{ps:>7.4g} {delta:>7.0e} " + " ".join(
                f"{v:>10.3f}{'c' if c else ' '}" for v, c in row))
    say("   (c = the returned eta equals Mc exactly: the dT_min forced accept's radial clamp fired)")
    # ---- ring replays
    say("")
    say("== ring rows x {isoComp, isoExt, shear+, shear-} x {1e-6, 1e-5}: alpha/alpha_b after one increment")
    say(f"{'integrator':<10} {'inc':>5} {'start in->out':>13} {'start out->in':>13} {'still out':>9} "
        f"{'clamp-sig':>9} {'f>1e-6':>7} {'wall s':>7}")
    rows = [r for p in sr.RING_CSVS for r in sr.read_ring_csv(p)]
    probes = []
    for dlt in (1e-6, 1e-5):
        probes += [[dlt, dlt, 0, 0, 0, 0], [-dlt, -dlt, 0, 0, 0, 0], [0, 0, 0, dlt, 0, 0], [0, 0, 0, -dlt, 0, 0]]
    ring = {}
    for name, be, tag in integ:
        t0 = time.perf_counter()
        n = io = oi = oo = ncl = nf = 0
        for r in rows:
            st = dict(sigma=r["sigma"], alpha=B.dev(r["alpha"]), alpha_in=B.dev(r["alpha_in"]),
                      z=B.dev(r["z"]), e=r["e"])
            a0 = B.bounding(st["sigma"], st["alpha"], st["e"])["alpha_over_b"]
            for de in probes:
                if be == "cpp":
                    new, info = D.step_cpp(tag, st, de, 0.0)
                else:
                    new, info = D.step_port(m_ae, st, de, 0.0)
                d = D.diag(new)
                a1 = d["alpha_over_b"]
                n += 1
                io += (a0 <= 1.0 < a1)
                oi += (a0 > 1.0 >= a1)
                oo += (a0 > 1.0 and a1 > 1.0)
                ncl += clamp_sig(d)
                nf += (info["rc"] == 0 and d["f"] > 1e-6)
        wall = time.perf_counter() - t0
        ring[name] = dict(n=n, in_to_out=io, out_to_in=oi, still_out=oo, clamp_sig=ncl, f_pos=nf, wall=wall)
        say(f"{name:<10} {n:>5} {io:>13} {oi:>13} {oo:>9} {ncl:>9} {nf:>7} {wall:>7.2f}")
    with open(f"{B.OUT}/q1e_alpha_error.json", "w") as fh:
        json.dump(dict(chains=res, threshold=thr, ring=ring), fh, indent=1, default=str)
    with open(f"{B.OUT}/q1e_alpha_error.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
