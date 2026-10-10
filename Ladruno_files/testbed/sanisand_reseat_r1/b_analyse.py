"""(b) analysis: every variant vs DM04 on the calibrated-behaviour battery.
Output: out/b_summary.md (+ prints)."""
from __future__ import annotations

import glob
import json
import os

import numpy as np

import r1common as C

OUTB = os.path.join(C.OUT, "b")
MONO = {"TCd_p25": 0, "TCd_p100": 0, "TCd_p400": 0, "TEd_p100": 0, "TCu_p100_e0.6944": 0,
        "TCu_p100_e0.80": 0, "PSd_p100": 0, "SSd_p100": 3, "SSu_p100": 3}
CYCS = ["CTXu_e0.6944_csr0.20", "CTXu_e0.80_csr0.10", "CSSu_e0.80_csr0.08", "CSSu_e0.6944_csr0.25"]
CYCD = ["CSSd_g0.001", "CSSd_g0.005", "CTXd_e0.0005"]


def load(test, v):
    p = os.path.join(OUTB, f"{test}__{v}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def variants_present():
    vs = set()
    for f in glob.glob(os.path.join(OUTB, "*__*.json")):
        vs.add(os.path.basename(f)[:-5].split("__")[1])
    order = ["DM04", "B1", "Badd1", "T1B1", "T2B1", "T2B0.25", "T2B1S", "T1B1S", "T1B0.5S", "T1B2S", "T1B0.5", "T1B2", "R150", "R150+T1"]
    return [v for v in order if v in vs]


def mono_metrics(ref, r, si):
    x0 = np.abs(np.array(ref["eps"])[:, si]); q0 = np.array(ref["q"])
    x1 = np.abs(np.array(r["eps"])[:, si]); q1 = np.array(r["q"])
    ev0 = np.array(ref["eps"])[:, :3].sum(1); ev1 = np.array(r["eps"])[:, :3].sum(1)
    xs = np.linspace(0, min(x0.max(), x1.max()), 2001)
    i0 = np.argsort(x0); i1 = np.argsort(x1)
    qa = np.interp(xs, x0[i0], q0[i0]); qb = np.interp(xs, x1[i1], q1[i1])
    va = np.interp(xs, x0[i0], ev0[i0]); vb = np.interp(xs, x1[i1], ev1[i1])
    dq = np.abs(qb - qa)
    # small-strain window (first 0.1 % of the controlling strain)
    w = xs <= 1e-3
    return dict(max_dq_rel=float(dq.max() / max(np.abs(qa).max(), 1e-12)),
                max_dq=float(dq.max()), x_at_max=float(xs[int(dq.argmax())]),
                dq_small_rel=float(dq[w].max() / max(np.abs(qa[w]).max(), 1e-12)) if w.any() else float("nan"),
                dqpeak_rel=float((np.abs(q1).max() - np.abs(q0).max()) / np.abs(q0).max()),
                max_dev=float(np.abs(vb - va).max()), status=r["status"], n_cap=r.get("n_cap"),
                reseats=r.get("reseats"))


def main():
    V = variants_present()
    lines = ["# (b) calibrated behaviour vs DM04 (campaign set)\n"]
    lines.append("## Monotonic: max |q - q_DM04| / max q_DM04 over the path  [small-strain window eps <= 1e-3]  (peak change)  {max |eps_v diff|}  caps(none/floor/soft)\n")
    hdr = "| test | " + " | ".join(V[1:]) + " |"
    lines += [hdr, "|" + "---|" * (len(V))]
    for t, si in MONO.items():
        ref = load(t, "DM04")
        if ref is None:
            continue
        cells = []
        for v in V[1:]:
            r = load(t, v)
            if r is None:
                cells.append("—"); continue
            if r.get("status") != "ok":
                cells.append(r.get("status")); continue
            m = mono_metrics(ref, r, si)
            cells.append(f"{m['max_dq_rel']:.1e} [{m['dq_small_rel']:.1e}] ({m['dqpeak_rel']:+.1e}) {{{m['max_dev']:.0e}}} c{m['n_cap']}")
        lines.append(f"| {t} (DM04 {ref['status']}, q_max {max(np.abs(ref['q'])):.1f}) | " + " | ".join(cells) + " |")
    lines.append("\n## Cyclic stress-controlled (undrained): N at DA / runaway (half-cycles/2), DA, p_end, re-seats, caps\n")
    lines += ["| test | " + " | ".join(V) + " |", "|" + "---|" * (len(V) + 1)]
    for t in CYCS:
        cells = []
        for v in V:
            r = load(t, v)
            if r is None:
                cells.append("—"); continue
            if "hist" not in r:
                cells.append(r.get("status", "?")[:40]); continue
            cells.append(f"{r['status']} N={r['n_liq']} DA={r['da']:.3f} p={r['p_end']:.1f} rs={r['reseats']} c{r['n_cap']}")
        lines.append(f"| {t} | " + " | ".join(cells) + " |")
    lines.append("\n## Cyclic strain-controlled (drained): cycle 1 / cycle 10 secant stiffness (kPa), damping, eps_v end\n")
    lines += ["| test | " + " | ".join(V) + " |", "|" + "---|" * (len(V) + 1)]
    for t in CYCD:
        cells = []
        for v in V:
            r = load(t, v)
            if r is None or "cycles" not in r or not r["cycles"]:
                cells.append("—" if r is None else r.get("status", "?")); continue
            c1, cN = r["cycles"][0], r["cycles"][-1]
            cells.append(f"k {c1['k_sec']:.4g}/{cN['k_sec']:.4g} D {c1['damping']:.4f}/{cN['damping']:.4f} ev {cN['epsv_end']:.3e} rs={r['reseats']}")
        lines.append(f"| {t} | " + " | ".join(cells) + " |")
    txt = "\n".join(lines)
    open(os.path.join(C.OUT, "b_summary.md"), "w", encoding="utf-8").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
