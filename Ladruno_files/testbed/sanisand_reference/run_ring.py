"""WP-134 validation (4): replay the 80 TIMs ring states through the reference
and the C++ with the documented WP-127 probes (isoComp, shear) at
delta = 1e-5 / 1e-4 / 1e-3.

    python Ladruno_files/testbed/sanisand_reference/run_ring.py [--deltas ...] [--rows N]

Writes out/ring.json, out/ring_rows.md (per row) and out/ring_summary.md."""
import argparse
import json
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))

import numpy as np  # noqa: E402

from sanisand_reference import CAMPAIGN, cxx  # noqa: E402
from sanisand_reference.crosscheck import job_of  # noqa: E402
from sanisand_reference.integrator import Control, integrate  # noqa: E402
from sanisand_reference.model import bounding_report, t2v, v2t  # noqa: E402
from sanisand_reference.ring import load_ring_csv, probes, ring_variants, row_state  # noqa: E402

OUT = os.path.join(HERE, "out")


def _one(args):
    st, de, O = args
    r = integrate(st, Control.strain(de), CAMPAIGN, O, record=False)
    return dict(status=r.status, t_end=r.t_end, f_end=r.f_end, max_rho=r.max_rho_b,
                rho_start=r.start["rho_b"], rho_end=r.end["rho_b"],
                rhoa_start=r.start["rho_alpha"], rhoa_end=r.end["rho_alpha"],
                max_rhoa=r.max_rho_alpha, eta_end=r.end["eta"],
                p_end=r.end["p"], sigma=t2v(r.state.sigma).tolist(),
                alpha=t2v(r.state.alpha).tolist(), z=t2v(r.state.z).tolist(),
                e=r.state.e, segs=[(s["mode"], s["event"]) for s in r.segments],
                reseats=len(r.reseats), negh=r.uw_negative_h, min_a=r.min_a_plastic,
                notes=r.notes[:2] + r.notes[-1:],
                f_start=r.start["f"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--deltas", type=float, nargs="*", default=[1e-5, 1e-4, 1e-3])
    ap.add_argument("--rows", type=int, default=0)
    ap.add_argument("--meshes", nargs="*", default=["b8", "b16"])
    args = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    V = ring_variants()
    cases = []
    for mesh in args.meshes:
        rows = load_ring_csv(mesh)
        if args.rows:
            rows = rows[:args.rows]
        for row in rows:
            st = row_state(row)
            for d in args.deltas:
                for pn, de in probes(d).items():
                    cases.append(dict(mesh=mesh, element=row["element"], gp=row["gp"],
                                      p0=row["p_kPa"], eta0=row["eta"],
                                      screen=row["eta_over_Mb_compression"], probe=pn,
                                      delta=d, deps=de, state=st,
                                      raw=dict(sigma=row["sigma"], alpha=row["alpha"],
                                               alpha_in=row["alpha_in"], z=row["z"],
                                               e=row["e"])))
    tasks = [(c["state"], c["deps"], O) for vn, O in V.items() for c in cases]
    with ProcessPoolExecutor() as ex:
        flat = list(ex.map(_one, tasks, chunksize=4))
    for k, vn in enumerate(V):
        for i, c in enumerate(cases):
            c[vn] = flat[k * len(cases) + i]
    have = cxx.available()
    if have:
        jobs = [dict(proto=pr, deps=c["deps"], **c["raw"]) for pr in ("ME", "ME8")
                for c in cases]
        cres = cxx.run_jobs(jobs, CAMPAIGN.as_opensees())
        for j, pr in enumerate(("ME", "ME8")):
            for i, c in enumerate(cases):
                x = cres[j * len(cases) + i]
                br = bounding_report(v2t(x["sigma"]), v2t(x["alpha"]), v2t(x["z"]),
                                     c["raw"]["e"], v2t(x["alpha_in"]), CAMPAIGN,
                                     V["uw_model"])
                ref = c["uw_model"]
                dsig = float(np.linalg.norm(np.array(x["sigma"]) - np.array(ref["sigma"])))
                dal = float(np.linalg.norm(np.array(x["alpha"]) - np.array(ref["alpha"])))
                c[pr] = dict(rc=x["rc"], substeps=x["substeps"], forced=x["forced"],
                             cap=x["cap"], f_after=x["f_after"], rho_end=br["rho_b"],
                             rhoa_end=br["rho_alpha"],
                             eta_end=br["eta"], p_end=x["p"], dsig_vs_ref=dsig,
                             dalpha_vs_ref=dal, sigma=x["sigma"], alpha=x["alpha"])
    for c in cases:
        c.pop("state")
    json.dump(cases, open(os.path.join(OUT, "ring.json"), "w"), indent=0, default=float)

    # ---- per-row table (uw_model = the oracle; paper and uw_rule for attribution) ----
    L = ["| mesh | el/gp | p0 | η0 | ρ_α0 (ρ_b0) | probe | δ | ref uw_model: status / f_end / max ρ_α / p_end / η_end | paper: status / max ρ_α | uw_rule: status / max ρ_α | C++ ME: rc / sub / ρ_α end / f | ‖Δσ‖ ME−ref (kPa) | C++ ME8: rc / sub / ρ_α end | ‖Δσ‖ ME8−ref |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in cases:
        u, p, w = c["uw_model"], c["paper"], c["uw_rule"]
        s = (f"| {c['mesh']} | {c['element']}/{c['gp']} | {c['p0']:.3g} | {c['eta0']:.2f} | {u['rhoa_start']:.2f} ({u['rho_start']:.2f}) | "
             f"{c['probe']} | {c['delta']:.0e} | {u['status']} / {u['f_end']:.1e} / {u['max_rhoa']:.2f} / {u['p_end']:.3g} / {u['eta_end']:.2f} | "
             f"{p['status']} / {p['max_rhoa']:.2f} | {w['status']} / {w['max_rhoa']:.2f} |")
        if have:
            m, m8 = c["ME"], c["ME8"]
            s += (f" {m['rc']} / {m['substeps']:.0f} / {m['rhoa_end']:.2f} / {m['f_after']:.1e} | {m['dsig_vs_ref']:.2e} |"
                  f" {m8['rc']} / {m8['substeps']:.0f} / {m8['rhoa_end']:.2f} | {m8['dsig_vs_ref']:.2e} |")
        L.append(s)
    open(os.path.join(OUT, "ring_rows.md"), "w", encoding="utf-8").write("\n".join(L) + "\n")

    # ---- summary ----
    S = []
    for d in args.deltas:
        for pn in ("isoComp", "shear"):
            sub = [c for c in cases if c["delta"] == d and c["probe"] == pn]
            adm = [c for c in sub if c["uw_model"]["rhoa_start"] <= 1.0]
            line = dict(delta=d, probe=pn, n=len(sub), n_adm=len(adm))
            for vn in V:
                st = {}
                for c in sub:
                    st[c[vn]["status"]] = st.get(c[vn]["status"], 0) + 1
                line[vn + "_status"] = st
                line[vn + "_max_rho_adm"] = max([c[vn]["max_rhoa"] for c in adm], default=float("nan"))
                line[vn + "_n_escape_adm"] = sum(1 for c in adm if c[vn]["max_rhoa"] > 1.0)
            if have:
                for pr in ("ME", "ME8"):
                    line[pr + "_n_escape_adm"] = sum(1 for c in adm if c[pr]["rhoa_end"] > 1.0)
                    line[pr + "_max_rho_adm"] = max([c[pr]["rhoa_end"] for c in adm], default=float("nan"))
                    ok = [c for c in adm if c["uw_model"]["status"] == "ok"]
                    ds = [c[pr]["dsig_vs_ref"] / max(np.linalg.norm(np.array(c["uw_model"]["sigma"]) - np.array(c["raw"]["sigma"])), 1e-300) for c in ok]
                    line[pr + "_dsig_rel_median"] = float(np.median(ds)) if ds else float("nan")
                    line[pr + "_dsig_rel_p95"] = float(np.percentile(ds, 95)) if ds else float("nan")
                    line[pr + "_n_fpos"] = sum(1 for c in sub if c[pr]["f_after"] > 1e-6)
            S.append(line)
    json.dump(S, open(os.path.join(OUT, "ring_summary.json"), "w"), indent=1, default=float)
    T = ["| δ | probe | rows (admissible) | uw_model status | uw_model max ρ_α (adm) / escapes | uw_rule status | uw_rule escapes | C++ ME escapes / max ρ_α / f>1e-6 | ME ‖Δσ−ref‖/‖Δσ_ref‖ median / p95 | C++ ME8 escapes / max ρ_α | ME8 rel median / p95 |",
         "|---|---|---|---|---|---|---|---|---|---|---|"]
    for l in S:
        s = (f"| {l['delta']:.0e} | {l['probe']} | {l['n']} ({l['n_adm']}) | {l['uw_model_status']} | "
             f"{l['uw_model_max_rho_adm']:.3f} / {l['uw_model_n_escape_adm']} | {l['uw_rule_status']} | {l['uw_rule_n_escape_adm']} |")
        if have:
            s += (f" {l['ME_n_escape_adm']} / {l['ME_max_rho_adm']:.2f} / {l['ME_n_fpos']} | {l['ME_dsig_rel_median']:.1e} / {l['ME_dsig_rel_p95']:.1e} |"
                  f" {l['ME8_n_escape_adm']} / {l['ME8_max_rho_adm']:.2f} | {l['ME8_dsig_rel_median']:.1e} / {l['ME8_dsig_rel_p95']:.1e} |")
        T.append(s)
    open(os.path.join(OUT, "ring_summary.md"), "w", encoding="utf-8").write("\n".join(T) + "\n")
    print("\n".join(T))


if __name__ == "__main__":
    main()
