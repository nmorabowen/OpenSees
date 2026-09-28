"""WP-128 Q1 (ii): does the dT_min forced accept (finding C) put alpha outside
the bounding surface on its own?  q1_search's vertUnload at delta = 1e-5
(p0 = 2 kPa; 200 fwd / 100 back x 3) reached alpha/alpha^b = 1.193 at k=299
(p 79.9 kPa) with 2 forced accepts, both Mc-clamped.  Re-run it on the port
as built and with the forced accept REFUSED, and open the forced increments.
Output: out/q1_forced.txt
"""
import _boot as B
import drive as D
import md_port

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


def main():
    B.define_prototypes()
    incs = cyc([0.0, -1.0, 0.0, 0.0, 0.0, 0.0], 1e-5, 200)
    st0 = D.k0_state(2.0)
    hc = D.run(st0, incs, "cpp", tag=B.TAG_ME)
    for name, m in (("as built", md_port.Material(B.P)),
                    ("forced accept REFUSED", md_port.Material(B.P, forced_policy="refuse")),
                    ("alpha also error-tested", md_port.Material(B.P, alpha_err=True))):
        h = D.run(st0, incs, "port", mat=m)
        mx = max(h, key=lambda x: x["alpha_over_b"])
        say(f"{name:<24} max alpha/alpha_b {mx['alpha_over_b']:.3f} at k={mx['k']} (p {mx['p']:.3g}, eta {mx['eta']:.3g}); "
            f"forced {sum(x['forced'] for x in h)}, refused increments {sum(1 for x in h if x['rc'] != 0)}")
        if name == "as built":
            dev = max(abs(a["alpha_over_b"] - b["alpha_over_b"]) for a, b in zip(hc, h))
            say(f"   (C++ vs port on this chain: max |alpha/alpha_b| difference {dev:.2e})")
            for x in h:
                if x["forced"]:
                    say(f"   forced increment k={x['k']}: path {x['path']}, forced {x['forced']}, clamp {x['clamp']}, "
                        f"alpha/alpha_b {x['before'] and D.diag(x['before'])['alpha_over_b']:.3f} -> {x['alpha_over_b']:.3f}, "
                        f"eta {D.diag(x['before'])['eta']:.3f} -> {x['eta']:.3f}, p {x['p']:.3g}")
    with open(f"{B.OUT}/q1_forced.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
