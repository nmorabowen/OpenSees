"""WP-128 Q1: the smallest reproducer of finding B, and its threshold.

Start state = what Stress_Correction's low-p branch leaves (sigma = p I,
alpha = 0) after a reversal has set alpha_in := alpha = 0, i.e.
(alpha - alpha_in):n = 0 and h = 1e10 (GetStateDependent's sentinel).
ONE plane-strain increment d_eps = (0, delta, 0, 0, 0, 0) (vertical
compression, compression-positive).  Sweep delta and the start pressure p_s.

Columns: C++ (IntScheme 1 campaign) rc / substeps / returned p, eta,
alpha/alpha^b_theta, f; and the port with the substep error ALSO measured on
alpha (the RK45-style counterfactual) for the same increment.
Output: out/q1_threshold.txt
"""
import _boot as B
import drive as D
import md_port

LINES = []


def say(s=""):
    print(s, flush=True)
    LINES.append(s)


def main():
    B.define_prototypes()
    m_ae = md_port.Material(B.P, alpha_err=True)
    say("start: sigma = p_s I, alpha = alpha_in = z = 0, e = 0.6978 (the vertUnload chain's)")
    say(f"{'p_s':>7} {'delta':>7} | {'rc':>3} {'sub':>5} {'p_out':>8} {'eta':>7} {'a/ab':>7} {'f':>9} {'dp/p_s':>7} | "
        f"alpha-err port: {'sub':>5} {'eta':>7} {'a/ab':>7}")
    for ps in (0.0101, 0.1, 1.0, 5.0):
        for delta in (1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4):
            st = dict(sigma=[ps, ps, ps, 0, 0, 0], alpha=[0.0] * 6, alpha_in=[0.0] * 6,
                      z=[0.0] * 6, e=0.697787979641054)
            de = [0.0, delta, 0.0, 0.0, 0.0, 0.0]
            new, info = D.step_cpp(B.TAG_ME, st, de, 0.0)
            d = D.diag(new)
            o = m_ae.update(st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"], de)
            pa = B.tr(o["sigma"]) / 3
            ab = B.bounding(o["sigma"], o["alpha"], o["e"])["alpha_over_b"]
            say(f"{ps:>7.4g} {delta:>7.0e} | {info['rc']:>3} {info['substeps']:>5} {d['p']:>8.4g} "
                f"{d['eta']:>7.3f} {d['alpha_over_b']:>7.3f} {d['f']:>9.2e} {d['p'] / ps:>7.1f} | "
                f"{'':>16}{o['substeps']:>5} {B.eta_sigma(o['sigma']):>7.3f} {ab:>7.3f}")
    with open(f"{B.OUT}/q1_threshold.txt", "w") as fh:
        fh.write("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
