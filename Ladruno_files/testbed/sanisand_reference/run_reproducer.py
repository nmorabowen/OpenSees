"""WP-134: WP-128's smallest reproducer, and the escape-onset sweep, through the
reference and the C++.

Start sigma = p_s I (compression positive), alpha = alpha_in = z = 0,
e = e_init (campaign), one plane-strain d eps_yy = +delta.
    python Ladruno_files/testbed/sanisand_reference/run_reproducer.py
Writes out/reproducer.json and out/reproducer.md."""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

from sanisand_reference import CAMPAIGN, cxx  # noqa: E402
from sanisand_reference.crosscheck import job_of  # noqa: E402
from sanisand_reference.integrator import Control, integrate  # noqa: E402
from sanisand_reference.model import bounding_report, v2t  # noqa: E402
from sanisand_reference.ring import reproducer_state, ring_variants  # noqa: E402

OUT = os.path.join(HERE, "out")
CASES = [(0.0101, 1e-5), (0.0101, 3e-5), (0.0101, 1e-4), (0.0101, 3e-4),
         (0.1, 1e-4), (1.0, 1e-4), (1.0, 1.5e-4), (1.0, 2e-4), (1.0, 2.5e-4),
         (1.0, 3e-4), (5.0, 3e-4), (5.0, 5e-4)]


def main():
    os.makedirs(OUT, exist_ok=True)
    P = CAMPAIGN
    V = ring_variants()
    res = []
    jobs = []
    for p_s, d in CASES:
        st = reproducer_state(p_s)
        de = [0.0, d, 0.0, 0.0, 0.0, 0.0]
        row = dict(p_s=p_s, delta=d)
        for vn, O in V.items():
            r = integrate(st, Control.strain(de), P, O)
            row[vn] = dict(status=r.status, eta=r.end["eta"], rho_end=r.end["rho_b"],
                           max_rho=r.max_rho_b, max_rhoa=r.max_rho_alpha,
                           f_end=r.f_end, p_end=r.end["p"],
                           dp_over_p=(r.end["p"] - p_s) / p_s, segs=len(r.segments),
                           reseats=len(r.reseats), negh=r.uw_negative_h)
        res.append(row)
        jobs += [job_of(st, de, "ME"), job_of(st, de, "ME8")]
    have_cxx = cxx.available()
    if have_cxx:
        out = cxx.run_jobs(jobs, P.as_opensees())
        for k, row in enumerate(res):
            for j, pr in enumerate(["ME", "ME8"]):
                c = out[2 * k + j]
                sig, al = v2t(c["sigma"]), v2t(c["alpha"])
                br = bounding_report(sig, al, v2t(c["z"]), c["e"], v2t(c["alpha_in"]),
                                     P, ring_variants()["uw_model"])
                row[pr] = dict(rc=c["rc"], substeps=c["substeps"], eta=br["eta"],
                               rho_end=br["rho_b"], rhoa=br["rho_alpha"],
                               f_after=c["f_after"], p=c["p"])
    json.dump(res, open(os.path.join(OUT, "reproducer.json"), "w"), indent=1)
    L = ["| p_s | delta | Δp/p (ref) | paper η / ρ_b | uw_model η / ρ_b (max ρ_b, max ρ_α) | uw_rule η / ρ_b | C++ ME rc/sub η / ρ_b / ρ_α / f | C++ ME8 rc/sub η / ρ_b |",
         "|---|---|---|---|---|---|---|---|"]
    for r in res:
        a, b, c = r["paper"], r["uw_model"], r["uw_rule"]
        s = (f"| {r['p_s']} | {r['delta']:.0e} | {b['dp_over_p']:.2f} | {a['eta']:.3f} / {a['rho_end']:.3f} | "
             f"{b['eta']:.3f} / {b['rho_end']:.3f} ({b['max_rho']:.3f}, {b['max_rhoa']:.3f}) | {c['eta']:.3f} / {c['rho_end']:.3f} {c['status']} |")
        if have_cxx:
            m, m8 = r["ME"], r["ME8"]
            s += (f" {m['rc']}/{m['substeps']:.0f} {m['eta']:.2f} / {m['rho_end']:.2f} / {m['rhoa']:.2f} / {m['f_after']:.1e} |"
                  f" {m8['rc']}/{m8['substeps']:.0f} {m8['eta']:.2f} / {m8['rho_end']:.2f} |")
        else:
            s += " n/a | n/a |"
        L.append(s)
    open(os.path.join(OUT, "reproducer.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
