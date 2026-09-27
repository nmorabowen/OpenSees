"""Digest of out/ring.json -> out/ring_analysis.md: mechanism counts, the
inadmissible rows, and the attribution of the C++-vs-reference differences."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")


def main():
    d = json.load(open(os.path.join(OUT, "ring.json")))
    L = []
    n = len(d)
    for vn in ("paper", "uw_model", "uw_rule"):
        negh = sum(1 for c in d if c[vn]["negh"])
        st = {}
        for c in d:
            st[c[vn]["status"]] = st.get(c[vn]["status"], 0) + 1
        adm = [c for c in d if c[vn]["rhoa_start"] <= 1.0]
        esc = [c for c in adm if c[vn]["max_rhoa"] > 1.0]
        esc_moved = [c for c in esc if abs(np.linalg.norm(c[vn]["alpha"]) /
                                           np.linalg.norm(c["raw"]["alpha"]) - 1) > 1e-6]
        L.append(f"- **{vn}**: {n} runs; statuses {st}; runs with a plastic sample at "
                 f"(α−α_in):n < −1e-12 (h < 0): **{negh}**; admissible starts {len(adm)}, of which "
                 f"max ρ_α > 1: {len(esc)} (α itself moved in {len(esc_moved)}; the rest are elastic "
                 f"paths on which α^b(ψ) contracted past a fixed α), max ρ_α "
                 f"{max([c[vn]['max_rhoa'] for c in adm], default=float('nan')):.3f}.")
    # C++ vs reference, by the C++ substep count
    L += ["", "C++ vs `uw_model` on admissible starts with a reference status ok "
          "(‖Δσ_C++ − Δσ_ref‖/‖Δσ_ref‖):", "",
          "| C++ | substeps ≤ 7 (err = 0 path): n / median / p95 | substeps > 7: n / median / p95 | escapes ρ_α > 1 (≤ 7 / > 7) |",
          "|---|---|---|---|"]
    for pr in ("ME", "ME8"):
        rows = [c for c in d if c["uw_model"]["rhoa_start"] <= 1 and c["uw_model"]["status"] == "ok"]

        def rel(c):
            dref = np.linalg.norm(np.array(c["uw_model"]["sigma"]) - np.array(c["raw"]["sigma"]))
            return c[pr]["dsig_vs_ref"] / max(dref, 1e-300)
        few = [rel(c) for c in rows if c[pr]["substeps"] <= 7]
        many = [rel(c) for c in rows if c[pr]["substeps"] > 7]
        ef = sum(1 for c in rows if c[pr]["substeps"] <= 7 and c[pr]["rhoa_end"] > 1)
        em = sum(1 for c in rows if c[pr]["substeps"] > 7 and c[pr]["rhoa_end"] > 1)
        fmt = lambda v: f"{len(v)} / {np.median(v):.1e} / {np.percentile(v, 95):.1e}" if v else "0 / - / -"
        L.append(f"| {pr} | {fmt(few)} | {fmt(many)} | {ef} / {em} |")
    # inadmissible rows
    L += ["", "Inadmissible start rows (ρ_α > 1):", "",
          "| row | ρ_α0 (ρ_b0) | f0 | probe | δ | ref uw_model: status, max ρ_α, end ρ_α, η_end, p_end, f_end | C++ ME: rc, substeps, ρ_α end, η_end, f_after |",
          "|---|---|---|---|---|---|---|"]
    for c in d:
        u = c["uw_model"]
        if u["rhoa_start"] > 1.0:
            m = c["ME"]
            note = ""
            if u["status"] != "ok" and u["notes"]:
                note = " — " + u["notes"][-1].split(": Required")[0]
            L.append(f"| {c['mesh']} {c['element']}/{c['gp']} | {u['rhoa_start']:.2f} ({u['rho_start']:.2f}) | "
                     f"{u['f_start']:.1e} | {c['probe']} | {c['delta']:.0e} | {u['status']}, {u['max_rhoa']:.2f}, "
                     f"{u['rhoa_end']:.2f}, {u['eta_end']:.2f}, {u['p_end']:.3g}, {u['f_end']:.1e}{note} | "
                     f"{m['rc']}, {m['substeps']:.0f}, {m['rhoa_end']:.2f}, {m['eta_end']:.2f}, {m['f_after']:.1e} |")
    txt = "\n".join(L) + "\n"
    open(os.path.join(OUT, "ring_analysis.md"), "w", encoding="utf-8").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
