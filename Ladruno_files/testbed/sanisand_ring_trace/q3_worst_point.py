"""WP-128 Q3 / F18(e): can ANY integrator take the b8 worst point
(el 1950 gp 3, p' 0.352 kPa, eta 12.87)?  Also prints the admissibility
diagnostics of the dumped states (Q1 finding B, Q2 f > 0).

Replays the row under each prototype (_boot.PROTOS) for the probe set
{isoComp, shear, isoExt(-), shear(-)} x {1e-7, 1e-6, 1e-5}, plus a zero
increment, and records rc / substeps / census / returned p, eta, f.
Output: out/q3_worst_point.json + stdout table.
"""
import json
import math

import _boot as B
from _boot import ops, sr
import md_port

ROWS = [(1950, 3), (1950, 2), (1859, 2)]


def probes():
    out = {"zero": [0.0] * 6}
    for d in (1e-7, 1e-6, 1e-5):
        out[f"isoComp@{d:.0e}"] = [d, d, 0, 0, 0, 0]
        out[f"isoExt@{d:.0e}"] = [-d, -d, 0, 0, 0, 0]
        out[f"shear+@{d:.0e}"] = [0, 0, 0, d, 0, 0]
        out[f"shear-@{d:.0e}"] = [0, 0, 0, -d, 0, 0]
    return out


def admissibility(r):
    sig, al, ai, z = r["sigma"], r["alpha"], r["alpha_in"], r["z"]
    bd = B.bounding(sig, al, r["e"])
    return dict(
        p=B.tr(sig) / 3, eta=B.eta_sigma(sig), eta_alpha=B.eta_alpha(al),
        f=B.yield_f(sig, al), f_over_p=B.yield_f(sig, al) / (B.tr(sig) / 3),
        cone_radius=math.sqrt(2 / 3) * B.M_ * B.tr(sig) / 3,
        Mb_theta=bd["Mb"], alpha_b_theta=bd["alpha_b"], psi=bd["psi"],
        alpha_over_alpha_b=bd["alpha_over_b"], b_dot_n=bd["b_dot_n"],
        alpha_minus_in_dot_n=B.ddot([al[i] - ai[i] for i in range(6)], bd["n"]),
        tr_alpha=B.tr(al), tr_alpha_in=B.tr(ai), tr_z=B.tr(z), norm_z=B.norm(z))


def main():
    B.define_prototypes()
    rows = {(r["element"], r["gp"]): r
            for r in sr.read_ring_csv(sr.RING_CSVS[0])}
    res = {}
    for key in ROWS:
        r = rows[key]
        adm = admissibility(r)
        print(f"== b8 el {key[0]} gp {key[1]}: " + ", ".join(
            f"{k}={v:.4g}" for k, v in adm.items()))
        res[str(key)] = dict(admissibility=adm, runs={})
        for tag, name in B.PROTOS.items():
            for pname, de in probes().items():
                o = sr.replay_row(ops, tag, r, de, trace=20000)
                s = o["stats"]
                sig, al = o["sigma"], o["alpha"]
                p = B.tr(sig) / 3
                row = dict(rc=o["rc"], path=o["path"], substeps=int(s["substeps"]),
                           accepted=int(s["accepted"]), rejErr=int(s["rejectedErr"]),
                           forced=int(s["forcedAtDTmin"]), clampMc=int(s["forcedClampMc"]),
                           abandonLowP=int(s["abandonedLowP"]), cap=int(s["capHits"]),
                           pnReset=int(s["pnResets"]), entryPmin=int(s["entryPminClamps"]),
                           p=p, eta=B.eta_sigma(sig), eta_alpha=B.eta_alpha(al),
                           f_before=o["f_before"], f_after=o["f_after"],
                           f_after_over_p=o["f_after"] / p if p else float("nan"),
                           maxErr=max([t["err"] for t in o["trace"]
                                       if math.isfinite(t["err"])], default=float("nan")),
                           sigma=sig, alpha=al)
                res[str(key)]["runs"][f"{name}|{pname}"] = row
                print(f"  {name:<15} {pname:<14} rc={row['rc']:>3} {row['path']:<18}"
                      f" sub={row['substeps']:>6} acc={row['accepted']:>6} "
                      f"forced={row['forced']:>5} clampMc={row['clampMc']:>5} "
                      f"abandon={row['abandonLowP']} cap={row['cap']} pnR={row['pnReset']} "
                      f"p={p:.4g} eta={row['eta']:.4g} eta_a={row['eta_alpha']:.4g} "
                      f"f0={row['f_before']:.3g} f1={row['f_after']:.3g} maxErr={row['maxErr']:.2g}")
        # the validated port, ME with alpha ALSO in the substep error (finding E)
        m = md_port.Material(B.P, alpha_err=True)
        for pname, de in probes().items():
            o = m.update(r["sigma"], B.dev(r["alpha"]), B.dev(r["alpha_in"]), B.dev(r["z"]), r["e"], de)
            sig, al = o["sigma"], o["alpha"]
            p = B.tr(sig) / 3
            row = dict(rc=o["rc"], path=o["path"], substeps=o["substeps"], accepted=o["acc"],
                       forced=o["forced"], clampMc=o["clamp"], corrGiveUp=o["corrGiveUp"],
                       p=p, eta=B.eta_sigma(sig), eta_alpha=B.eta_alpha(al),
                       f_after=o["f_after"], sigma=sig, alpha=al)
            res[str(key)]["runs"][f"ME+aErr(port)|{pname}"] = row
            print(f"  {'ME+aErr(port)':<15} {pname:<14} rc={row['rc']:>3} path={row['path']:<3}"
                  f" sub={row['substeps']:>6} acc={row['accepted']:>6} forced={row['forced']:>5} "
                  f"clampMc={row['clampMc']:>5} corrGiveUp={row['corrGiveUp']} "
                  f"p={p:.4g} eta={row['eta']:.4g} eta_a={row['eta_alpha']:.4g} f1={row['f_after']:.3g}")
    with open(f"{B.OUT}/q3_worst_point.json", "w") as fh:
        json.dump(res, fh, indent=1)


if __name__ == "__main__":
    main()
