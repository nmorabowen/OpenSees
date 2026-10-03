"""Figures for the R1 memo (WP-151).  Static PNGs for a markdown memo: reference
palette in fixed categorical order, 2 px lines, hairline grid, one y-axis per panel,
legend whenever >= 2 series.  Run after the result JSONs exist:
    py -3.11 r1plots.py [zeno fan mono gate continuity objectivity chains]"""
from __future__ import annotations

import glob
import json
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
FIG = os.path.join(OUT, "fig")
os.makedirs(FIG, exist_ok=True)

# reference palette (dataviz skill, references/palette.md), light mode
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 9,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
    "axes.titlesize": 10, "axes.titleweight": "semibold", "axes.titlelocation": "left",
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False,
    "lines.linewidth": 2.0, "legend.frameon": False, "legend.fontsize": 8,
})


def save(fig, name):
    p = os.path.join(FIG, name)
    fig.savefig(p, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", p)


# ---------------------------------------------------------------------------
def fig_zeno():
    Z = {z["variant"]: z for z in json.load(open(os.path.join(OUT, "zeno_trace.json")))}
    fig, axs = plt.subplots(1, 3, figsize=(10.5, 3.2))
    for v, col, lab in (("DM04", BLUE, "DM04 (h = b0/a)"), ("B1", ORANGE, "floor only (c_A = 1)")):
        rs = Z[v]["reseats"]
        t = np.array([r["t"] for r in rs]); bn = np.array([r["bn"] for r in rs])
        da = np.array([r["dalpha"] for r in rs])
        dt = np.diff(t)
        m = dt > 0
        k = np.arange(1, len(t))[m]
        axs[0].semilogy(k, dt[m], "o-", color=col, ms=4, label=lab)
        kk = np.arange(len(bn))
        uniq = np.r_[True, np.abs(np.diff(bn)) > 1e-15]
        axs[1].semilogy(kk[uniq], np.abs(bn[uniq]), "o-", color=col, ms=4, label=lab)
        axs[2].semilogy(kk[uniq], da[uniq], "o-", color=col, ms=4, label=lab)
    axs[0].set_title("Re-seat intervals shrink geometrically")
    axs[0].set_xlabel("re-seat k"); axs[0].set_ylabel("t(k) - t(k-1)  (pseudo-time)")
    axs[1].set_title("b:n at each re-seat -> 0+")
    axs[1].set_xlabel("re-seat k"); axs[1].set_ylabel("b:n")
    axs[2].set_title("|d alpha / dt| at each re-seat")
    axs[2].set_xlabel("re-seat k"); axs[2].set_ylabel("|d alpha/dt|  (per unit t)")
    for ax in axs:
        ax.legend(loc="best")
    t1 = Z["T1B1S"]
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.27)
    fig.text(0.01, 0.01, f"E_B 1880/1 committed wall state, one trial increment |d eps| = 3e-6 (a chatter "
             f"direction), exact Radau. DM04 stops ('chatter') at t* = {Z['DM04']['t_end']:.9f} after "
             f"{Z['DM04']['n_reseats']} re-seats; floor-only at t = {Z['B1']['t_end']:.6f} after "
             f"{Z['B1']['n_reseats']}; floor + hysteresis (+cap): {t1['n_reseats']} re-seats, t = 1 reached.",
             color=INK2, fontsize=8, wrap=True)
    save(fig, "fig1_zeno.png")


# ---------------------------------------------------------------------------
def fig_fan():
    R = []
    for f in ("a3_fan.json", "a3_fan_r150.json", "a3_fan_cA.json"):
        p = os.path.join(OUT, f)
        if os.path.exists(p):
            R += json.load(open(p))
    order = ["DM04", "R150", "R150+T1", "B1", "B0.25", "Badd1", "B1S", "T2B0.25", "T1B0.5",
             "T0.5B1", "T1B1", "T2B1", "T1B2", "T1B0.5S", "T1B1S", "T2B1S", "T1B2S"]
    label = {"DM04": "DM04", "R150": "WP-150 R1 (floor if b:n<=0)", "R150+T1": "WP-150 R1 + R1b",
             "B1": "floor c_A=1 only", "B0.25": "floor c_A=0.25 only", "Badd1": "additive b0/(a+eps) only",
             "B1S": "floor c_A=1 + cap", "T2B0.25": "hyst 2 + floor 0.25", "T1B0.5": "hyst 1 + floor 0.5",
             "T0.5B1": "hyst 0.5 + floor 1", "T1B1": "hyst 1 + floor 1", "T2B1": "hyst 2 + floor 1",
             "T1B2": "hyst 1 + floor 2", "T1B0.5S": "hyst 1 + floor 0.5 + cap",
             "T1B1S": "hyst 1 + floor 1 + cap", "T2B1S": "hyst 2 + floor 1 + cap",
             "T1B2S": "hyst 1 + floor 2 + cap"}
    kinds = [("chatter", "Zeno chatter", BLUE), ("solver_failed", "0/0 (solver stops)", ORANGE),
             ("H_nonpositive", "H <= 0 (no rate solution)", AQUA)]
    present = [v for v in order if any(r["variant"] == v for r in R)]
    cnt = {v: {k: 0 for k, _, _ in kinds} for v in present}
    tot = {v: 0 for v in present}
    for r in R:
        if r["variant"] in cnt:
            tot[r["variant"]] += 1
            if r["status"] in cnt[r["variant"]]:
                cnt[r["variant"]][r["status"]] += 1
    fig, ax = plt.subplots(figsize=(7.5, 0.34 * len(present) + 1.0))
    y = np.arange(len(present))[::-1]
    left = np.zeros(len(present))
    for k, lab, col in kinds:
        w = np.array([cnt[v][k] for v in present], dtype=float)
        ax.barh(y, w, left=left, color=col, height=0.62, label=lab, edgecolor=SURF, linewidth=1.5)
        left += w
    for yi, v, L in zip(y, present, left):
        ax.text(L + 3, yi, f"{int(L)}/{tot[v]}", va="center", color=INK2, fontsize=8)
    ax.set_yticks(y, [label[v] for v in present])
    ax.set_xlabel("failed trial increments (5 refuser states x 32 directions x 2 magnitudes)")
    ax.set_xlim(0, 130)
    ax.grid(axis="y", visible=False)
    ax.legend(loc="lower right", ncols=1)
    ax.set_title("Only floor-everywhere + hysteresis integrates every trial")
    save(fig, "fig2_fan.png")


# ---------------------------------------------------------------------------
MONO = [("TCd_p25", 0, "drained TC, p0 25"), ("TCd_p100", 0, "drained TC, p0 100"),
        ("TCd_p400", 0, "drained TC, p0 400"), ("TEd_p100", 0, "drained TE, p0 100"),
        ("PSd_p100", 0, "drained plane strain, p0 100"), ("SSd_p100", 3, "drained simple shear, p0 100"),
        ("TCu_p100_e0.6944", 0, "undrained TC, e0 0.694"), ("TCu_p100_e0.80", 0, "undrained TC, e0 0.80"),
        ("SSu_p100", 3, "undrained simple shear, p0 100")]


def _load_b(test, v):
    p = os.path.join(OUT, "b", f"{test}__{v}.json")
    return json.load(open(p)) if os.path.exists(p) else None


def fig_mono(variant="T1B1S"):
    fig, axs = plt.subplots(3, 3, figsize=(10.5, 8.2))
    for ax, (t, si, title) in zip(axs.flat, MONO):
        ref = _load_b(t, "DM04")
        var = _load_b(t, variant) or _load_b(t, "T2B1S")
        add = _load_b(t, "Badd1")
        if ref is None or var is None:
            ax.set_visible(False)
            continue
        x0 = np.abs(np.array(ref["eps"])[:, si]) * 100
        q0 = np.array(ref["q"])
        x1 = np.abs(np.array(var["eps"])[:, si]) * 100
        q1 = np.array(var["q"])
        ax.plot(x0, q0, color=BLUE, label="DM04")
        ax.plot(x1, q1, color=ORANGE, lw=1.4, ls=(0, (4, 3)), label="R1 (hyst 1 + floor 1 + cap)")
        xs = np.linspace(0, min(x0.max(), x1.max()), 1500)
        i0, i1 = np.argsort(x0), np.argsort(x1)
        d1 = np.abs(np.interp(xs, x1[i1], q1[i1]) - np.interp(xs, x0[i0], q0[i0])).max() / np.abs(q0).max()
        txt = f"max |dq|/q_max: R1 {d1:.1e}"
        if add is not None:
            x2 = np.abs(np.array(add["eps"])[:, si]) * 100; q2 = np.array(add["q"]); i2 = np.argsort(x2)
            xs2 = np.linspace(0, min(x0.max(), x2.max()), 1500)
            d2 = np.abs(np.interp(xs2, x2[i2], q2[i2]) - np.interp(xs2, x0[i0], q0[i0])).max() / np.abs(q0).max()
            txt += f"\nadditive b0/(a+eps): {d2:.1e}"
        ax.text(0.97, 0.05, txt, transform=ax.transAxes, ha="right", va="bottom", color=INK2, fontsize=7.5)
        ax.set_title(title)
        ax.set_xlabel("gamma (%)" if si == 3 else "|eps_a| (%)")
        ax.set_ylabel("q (kPa)")
    axs.flat[0].legend(loc="upper left")
    fig.suptitle("Monotonic element tests, campaign set: R1 and DM04 coincide (curves overlap)",
                 x=0.01, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout()
    save(fig, "fig3_monotonic.png")


# ---------------------------------------------------------------------------
def fig_gate():
    G = json.load(open(os.path.join(OUT, "cyc_gate.json")))
    Tt = json.load(open(os.path.join(OUT, "cyc_toyoura.json")))
    fig, axs = plt.subplots(1, 3, figsize=(12.5, 3.8), gridspec_kw=dict(width_ratios=[1.3, 0.8, 1.1]))
    models = ["DM04", "B1", "T1B1", "T1B1S", "T2B1"]
    col = dict(zip(models, SERIES))
    ax = axs[0]
    perts = [-1e-8, -1e-9, 0.0, 1e-9, 1e-8]
    for i, v in enumerate(models):
        for j, pz in enumerate(perts):
            r = [g for g in G if g["tag"] == "c0.71" and g["variant"] == v and g["pert"] == pz]
            if r:
                L = r[0]["status"] == "DA_reached"
                ax.scatter(j + (i - 2) * 0.13, 8.0 if L else 20.5, s=46, color=col[v], edgecolor=SURF,
                           linewidth=1.2, zorder=3, label=v if j == 0 else None)
    ax.set_xticks(range(len(perts)), ["-1e-8", "-1e-9", "0", "+1e-9", "+1e-8"])
    ax.set_yticks([8, 20.5], ["5% DA at N = 8", "none by N = 20"]); ax.set_ylim(5, 23.5)
    ax.set_xlabel("relative perturbation of sigma_zz at the start")
    ax.set_title("campaign c = 0.71: same branch for every model")
    ax.legend(loc="center right", fontsize=7.5); ax.grid(axis="x", visible=False)
    ax = axs[1]
    for i, v in enumerate(["DM04", "T1B1S"]):
        for j, pz in enumerate([-1e-9, 0.0, 1e-9]):
            r = [g for g in G if g["tag"] == "c0.80" and g["variant"] == v and g["pert"] == pz]
            if r:
                ax.scatter(j + (i - 0.5) * 0.18, r[0]["n_liq"] or 20.5, s=46, color=col[v], edgecolor=SURF,
                           linewidth=1.2, zorder=3, label=v if j == 0 else None)
    ax.set_xticks(range(3), ["-1e-9", "0", "+1e-9"]); ax.set_ylim(5, 23.5)
    ax.set_yticks([8, 16, 20.5], ["8", "16", "none by 20"]); ax.set_ylabel("cycles to 5% DA")
    ax.set_xlabel("perturbation"); ax.set_title("control c = 0.80 (> 7/9)")
    ax.legend(loc="upper right", fontsize=7.5); ax.grid(axis="x", visible=False)
    ax = axs[2]
    for k, csr in enumerate((0.15, 0.20)):
        for i, v in enumerate(["DM04", "B1", "T1B1", "T1B1S"]):
            for j, pz in enumerate([-1e-9, 0.0, 1e-9]):
                r = [x for x in Tt if x["stage"] == "gate" and x["csr"] == csr and x["variant"] == v and x["pert"] == pz]
                if r:
                    ax.scatter(3 * k + j + (i - 1.5) * 0.14, r[0]["n_liq"] or 20.5, s=40, color=col[v],
                               edgecolor=SURF, linewidth=1.2, zorder=3, label=v if (k == 0 and j == 0) else None)
    ax.set_xticks(range(6), ["-1e-9", "0", "+1e-9"] * 2)
    ax.text(1, 21.5, "CSR 0.15", ha="center", color=INK2, fontsize=8); ax.text(4, 21.5, "CSR 0.20", ha="center", color=INK2, fontsize=8)
    ax.set_ylim(5, 23.5); ax.set_ylabel("cycles to 5% DA"); ax.set_xlabel("perturbation")
    ax.set_title("Toyoura (DM04 set, c = 0.712): R1 = DM04")
    ax.legend(loc="lower left", fontsize=7.5); ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.26)
    fig.text(0.01, 0.01, "Undrained cyclic triaxial, p0 100 kPa, exact Radau; campaign e0 0.6944 CSR 0.2, Toyoura e0 0.808. "
             "With c < 7/9 the axisymmetric extension path is unstable: an imposed 1e-9 perturbation picks the branch, "
             "identically for DM04 and every R1 variant; unperturbed runs are picked by round-off. At c = 0.80 the path is stable.",
             color=INK2, fontsize=8, wrap=True)
    save(fig, "fig4_ctxu_gate.png")


# ---------------------------------------------------------------------------
def fig_continuity():
    R = json.load(open(os.path.join(OUT, "c3.json")))
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 6.2), sharex=True)
    states = [("E_B", 7516, "E_B 1880/1"), ("E_B16", 31279, "E_B16 7820/4")]
    for col, (leg, k, name) in enumerate(states):
        for row, (comp, lab) in enumerate(((0, "d sigma_xx (kPa)"), (3, "d sigma_xy (kPa)"))):
            ax = axs[row, col]
            rug = []
            for i, (v, vl) in enumerate((("DM04", "DM04"), ("T1B1S", "R1 (hyst 1 + floor 1 + cap)"),
                                         ("R150", "WP-150 R1"))):
                rr = sorted([r for r in R if r["leg"] == leg and r["k"] == k and r["variant"] == v],
                            key=lambda r: r["theta"])
                th = np.array([math.degrees(r["theta"]) for r in rr])
                ok = np.array([r["status"] == "ok" for r in rr])
                d = np.array([r["dsig"][comp] if "dsig" in r else np.nan for r in rr])
                c = [BLUE, ORANGE, AQUA][i]
                dd = d.copy(); dd[~ok] = np.nan
                ax.plot(th, dd, color=c, lw=[2.0, 1.4, 1.2][i], ls=["-", (0, (4, 3)), ":"][i], label=vl)
                if (~ok).any():   # a rug of 'no answer' directions below the curves, one row per model
                    ax.scatter(th[~ok], np.full((~ok).sum(), np.nan), s=0)  # keep autoscale honest
                    rug.append((th[~ok], c, f"{vl}: no answer"))
            y0, y1 = ax.get_ylim()
            span = y1 - y0
            for j, (tt, cc, ll) in enumerate(rug):
                ax.scatter(tt, np.full(len(tt), y0 - (0.06 + 0.07 * j) * span), marker="|", s=40,
                           color=cc, linewidth=1.4, label=ll, clip_on=False)
            ax.set_ylim(y0 - (0.06 + 0.07 * max(len(rug), 1)) * span, y1)
            ax.set_title(f"{name}: {lab.split(' (')[0]}")
            ax.set_ylabel(lab)
            if row == 1:
                ax.set_xlabel("trial direction theta (deg) around the last converged increment")
    axs[0, 0].legend(loc="best", fontsize=7)
    fig.suptitle("Trial response vs direction (|d eps| = 3e-6): DM04 has holes; R1 is continuous",
                 x=0.01, ha="left", fontsize=11, fontweight="semibold")
    fig.tight_layout()
    save(fig, "fig5_continuity.png")


def fig_objectivity():
    R = json.load(open(os.path.join(OUT, "c2.json")))
    fig, axs = plt.subplots(1, 2, figsize=(10.5, 3.6))
    for ax, (leg, k, name) in zip(axs, (("E_B", 7516, "E_B 1880/1"), ("E_B16", 31279, "E_B16 7820/4"))):
        for i, v in enumerate(["DM04", "B1", "R150", "T1B1", "T1B1S"]):
            rr = {r["lam"]: r for r in R if r["leg"] == leg and r["k"] == k and r["variant"] == v}
            if 0.0 not in rr or rr[0.0]["status"] != "ok":
                n0 = rr.get(0.0, {}).get("n_ok", 0)
                ax.text(0.02, 0.95 - 0.07 * i, f"{v}: smooth path fails after {n0}/60 increments",
                        transform=ax.transAxes, color=SERIES[i], fontsize=7.5, va="top")
                continue
            s0 = np.array(rr[0.0]["sigma_end"])
            lams, dev = [], []
            for lam in (0.01, 0.03, 0.1, 0.3, 1.0):
                r = rr.get(lam)
                if r is None or r["status"] != "ok":
                    continue
                lams.append(lam); dev.append(np.linalg.norm(np.array(r["sigma_end"]) - s0))
            if lams:
                ax.loglog(lams, dev, "o-", color=SERIES[i], ms=4, label=v)
        ref = [r for r in R if r["leg"] == leg and r["k"] == k and r["variant"] == "T1B1S" and r["lam"] == 0.01]
        s0 = [r for r in R if r["leg"] == leg and r["k"] == k and r["variant"] == "T1B1S" and r["lam"] == 0.0]
        if ref and s0 and ref[0]["status"] == "ok" and s0[0]["status"] == "ok":
            d1 = np.linalg.norm(np.array(ref[0]["sigma_end"]) - np.array(s0[0]["sigma_end"]))
            lam = np.array([0.01, 1.0]); ax.loglog(lam, d1 / 0.01 * lam * 0.5, color=MUTED, lw=1, ls="--",
                                                   label="slope 1 (guide)")
        ax.text(0.98, 0.30, "T1B1 and T1B1S coincide", transform=ax.transAxes, ha="right", color=INK2, fontsize=7.5)
        ax.set_title(f"{name}: end stress vs jitter amplitude")
        ax.set_xlabel("jitter amplitude lambda (relative to the increment)")
        ax.set_ylabel("|sigma_end(lambda) - sigma_end(0)| (kPa)")
        ax.legend(loc="lower right", fontsize=7.5)
    fig.tight_layout()
    save(fig, "fig6_objectivity.png")


if __name__ == "__main__":
    which = sys.argv[1:] or ["zeno", "fan", "mono", "gate"]
    for w in which:
        {"zeno": fig_zeno, "fan": fig_fan, "mono": fig_mono, "gate": fig_gate,
         "continuity": fig_continuity, "objectivity": fig_objectivity}[w]()
