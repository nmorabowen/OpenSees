"""WP-128 Q1: attribute the alpha-outside-the-bounding-surface escape found by
q1_search.py.  Smallest reproducer: plane-strain vertical extension
d_eps = (0, -delta, 0, 0, 0, 0) from a K0 state at p0 = 2 kPa, delta = 1e-4
(the search's `vertUnload`), plus the `extShear` path.

1. C++ chain vs port chain agree increment by increment (the port is what
   lets the substeps be opened up).
2. Per increment up to the first escape: path, census, alpha/alpha^b.
3. Inside the escaping increment: every substep, with the two Heun stages'
   kind (plastic / elasticDrag = the dgamma < 0 branch that moves alpha by
   the stress-ratio change / neutral), h, Kp, b:n and alpha/alpha^b after it.
4. Counterfactuals on the port (one switch at a time): which mechanism,
   removed, removes the escape.
Output: out/q1_attrib.txt, out/q1_attrib.json
"""
import json
import math

import _boot as B
import drive as D
import md_port

LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def ab_ratio(S, A, e):
    if B.tr(S) <= 0:
        return float("nan")
    return B.bounding(S, A, e)["alpha_over_b"]


def chain_summary(h):
    r = [x["alpha_over_b"] for x in h if math.isfinite(x["alpha_over_b"])]
    first = next((x["k"] for x in h if x["alpha_over_b"] > 1.05), None)
    return max(r) if r else float("nan"), first


PATHS = {
    "vertUnload": [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
    "extShear": [0.3, -1.0, 0.0, 0.5, 0.0, 0.0],
}


def incs_for(d, delta, n):
    """q1_search's cycle: n forward, n/2 back, x3 (total strain per leg n*delta)."""
    out = []
    for _ in range(3):
        out += [[delta * x for x in d]] * n
        out += [[-delta * x for x in d]] * (n // 2)
    return out


def main():
    B.define_prototypes()
    out = {}
    for pname, d in PATHS.items():
        say(f"=================== {pname}, p0 = 2 kPa, delta = 1e-4, (20 fwd, 10 back) x 3")
        st0 = D.k0_state(2.0)
        incs = incs_for(d, 1e-4, 20)
        hc = D.run(st0, incs, "cpp", tag=B.TAG_ME)
        mat = md_port.Material(B.P)
        hp = D.run(st0, incs, "port", mat=mat)
        worst_dev = 0.0
        for a, b in zip(hc, hp):
            num = math.sqrt(sum((x - y) ** 2 for x, y in zip(a["after"]["sigma"], b["after"]["sigma"])))
            den = max(math.sqrt(sum(x * x for x in a["after"]["sigma"])), 1e-3)
            worst_dev = max(worst_dev, num / den)
        say(f"C++ vs port chain: worst rel sigma difference over the chain {worst_dev:.2e}")
        say(" k  path  (the first 20 increments drive p to the p_min floor; k=20 is the first REVERSAL)")
        say(" k  path              sub  forced clamp  corrGiveUp  p_after   eta_after  eta_alpha  alpha/alpha_b  psi     Mb      f_after")
        first = None
        for x, y in zip(hc, hp):
            bd = B.bounding(x["after"]["sigma"], x["after"]["alpha"], x["after"]["e"]) if x["p"] > 0 else {"psi": float("nan"), "Mb": float("nan")}
            say(f"{x['k']:>2}  {x['path']:<16} {x['substeps']:>5} {x['forced']:>6} {x['clamp']:>5} "
                f"{y.get('corrGiveUp', 0):>10}  {x['p']:>8.4g} {x['eta']:>9.4g} {x['eta_alpha']:>9.4g} "
                f"{x['alpha_over_b']:>12.4g}  {bd['psi']:>7.4f} {bd['Mb']:>6.3f} {x['f']:>9.3g}")
            if first is None and x["alpha_over_b"] > 1.0:
                first = x["k"]
        out[pname] = dict(first_cross=first, cpp_vs_port=worst_dev)
        if first is None:
            continue
        # ---- open the escaping increment on the port
        st = hp[first]["before"]
        prev = D.ncov(incs[first - 1]) if first > 0 else 0.0
        o = mat.update(st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"],
                       incs[first], prev_incr_norm=prev)
        say(f"--- increment k={first} opened (port): path {o['path']}, substeps {o['substeps']}, "
            f"alpha_in reset: {o['alpha_in'] != st['alpha_in']}")
        say("   sub   T         dT        code  kinds                      h1          h2          Kp1        b:n(1)    a/ab before  a/ab after")
        kinds_count = {}
        crossings = []
        for i, t in enumerate(o["trace"]):
            T, dT, err, code, atmin, info = t
            if info is None:
                continue
            before = ab_ratio(info["S"], info["A"], info["e"])
            after = ab_ratio(info["nS"], info["nA"], info["e"]) if code in (0,) else before
            kinds_count[info["kinds"]] = kinds_count.get(info["kinds"], 0) + (code == 0)
            if code == 0 and before <= 1.0 < after:
                crossings.append(i)
            if i < 12 or (code == 0 and before <= 1.0 < after) or i == len(o["trace"]) - 1:
                say(f"   {i:>4} {T:.6f} {dT:.3e} {code:>4}  {str(info['kinds']):<26} {info['h'][0]:>11.4g} "
                    f"{info['h'][1]:>11.4g} {info['Kp'][0]:>11.4g} {info['bn'][0]:>9.4g} {before:>10.4g} {after:>10.4g}")
        say(f"   accepted substeps by stage kinds: {kinds_count}")
        say(f"   accepted substeps that carried alpha across the bounding surface: {len(crossings)}")
        out[pname]["kinds_in_first_escape"] = {str(k): v for k, v in kinds_count.items()}
        out[pname]["crossing_substeps"] = len(crossings)

        # ---- counterfactuals on the whole 40-increment chain
        say("--- counterfactuals (port, same chain): max alpha/alpha_b, first k > 1.05")
        cf = {
            "as built (C++-faithful)": md_port.Material(B.P),
            "dgamma<0 stage leaves alpha alone (drag frozen)": md_port.Material(B.P, drag="frozen"),
            "forced accept at dT_min REFUSES": md_port.Material(B.P, forced_policy="refuse"),
            "Stress_Correction off": md_port.Material(B.P, correction=False),
            "error test also on alpha (RK45-style)": md_port.Material(B.P, alpha_err=True),
            "TolE 1e-8": md_port.Material(B.P, TolE=1e-8),
        }
        out[pname]["counterfactuals"] = {}
        for name, m in cf.items():
            h = D.run(st0, incs, "port", mat=m)
            mx, fk = chain_summary(h)
            nref = sum(1 for x in h if x["rc"] != 0)
            say(f"   {name:<50} max a/ab = {mx:8.3f}  first k>1.05 = {fk}  refused increments = {nref}")
            out[pname]["counterfactuals"][name] = dict(max=mx, first=fk, refused=nref)
        for delta, n in ((1e-5, 200), (1e-6, 2000)):
            h = D.run(st0, incs_for(d, delta, n), "cpp", tag=B.TAG_ME)
            mx, fk = chain_summary(h)
            say(f"   C++, same total strain at delta={delta:.0e} ({n} increments): max a/ab = {mx:.3f}, first k>1.05 = {fk}, "
                f"min p = {min(x['p'] for x in h):.3g}")
            out[pname][f"cpp_delta_{delta:.0e}"] = dict(max=mx, first=fk)
    with open(f"{B.OUT}/q1_attrib.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")
    with open(f"{B.OUT}/q1_attrib.json", "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    main()
