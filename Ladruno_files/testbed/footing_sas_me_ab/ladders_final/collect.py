"""Pull ladder legs + E_B from esmeralda, tabulate, and plot q-s per ladder.
    python collect.py            (run from anywhere; writes next to this file)"""
import csv, json, os, re, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SCR = os.path.dirname(os.path.abspath(__file__))
R = "ladruno_wp138/deck/runs"
PRES = [("E_B", 0.0), ("L_pres_0p5", 0.5), ("L_pres_1", 1), ("L_pres_2", 2),
        ("L_pres_5", 5), ("L_pres_10", 10), ("L_pres_20", 20)]
EINI = [("L_e_0p65", 0.65), ("E_B", 0.6944), ("L_e_0p75", 0.75), ("L_e_0p80", 0.80),
        ("L_e_0p85", 0.85)]
ABL = [("E_B", "base"), ("S1_nofabric", "zmax0"), ("S2_nopeak", "+nb0"),
       ("S3_critstate", "+nd0"), ("S4_nodilat", "+A0=.001"), ("ctrl_dp38", "DP psi=0")]
AH = [("E_B", "A0=.05 h0=1.3"), ("A0_0p02", "A0=.02"), ("A0_0p10", "A0=.10"),
      ("h0_x3", "h0=3.9")]
LOCAL = {"ctrl_dp38"}          # already in SCR/runs (from the WP-138 analysis folder)
legs = sorted({l for l, _ in PRES + EINI + ABL + AH} - LOCAL)
SSH = ["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", "esmeralda"]

if "--nopull" not in sys.argv:
    for l in legs:
        d = os.path.join(SCR, "runs", l, "logs")
        os.makedirs(d, exist_ok=True)
        for f in ("steps.csv", "summary.json", "logs/log.log"):
            dst = os.path.join(SCR, "runs", l, f)
            subprocess.run(["scp", "-q", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20",
                            f"esmeralda:{R}/{l}/{f}", dst], timeout=300)


def load(l):
    p = os.path.join(SCR, "runs", l, "steps.csv")
    if not os.path.exists(p):
        return None
    rows = list(csv.DictReader(open(p)))
    if not rows:
        return None
    s = np.array([float(r["s_over_B"]) for r in rows])
    q = np.array([float(r["q_kPa"]) for r in rows])
    return s, q


def refusals(l):
    """sum the per-step 'refusals step N: ... by code {...}' dicts in the log."""
    tot = {}
    p = os.path.join(SCR, "runs", l, "logs", "log.log")
    if not os.path.exists(p):
        return tot
    for m in re.finditer(r"refusals step \d+: \d+ by code (\{[^}]*\})", open(p, errors="replace").read()):
        for k, v in eval(m.group(1)).items():
            tot[k] = tot.get(k, 0) + v
    return tot


def first_nonposh(l):
    """s/B of the first step whose refusal census contains loadingNonPosH."""
    p = os.path.join(SCR, "runs", l, "logs", "log.log")
    if not os.path.exists(p) or D.get(l) is None:
        return float("nan")
    for m in re.finditer(r"refusals step (\d+): \d+ by code (\{[^}]*\})",
                         open(p, errors="replace").read()):
        if "loadingNonPosH" in m.group(2):
            k = int(m.group(1))
            return float(D[l][0][k - 1]) if k <= len(D[l][0]) else float("nan")
    return float("nan")


def status(l):
    p = os.path.join(SCR, "runs", l, "summary.json")
    return json.load(open(p)) if os.path.exists(p) else None


def q_at(sq, sB):
    s, q = sq
    return float(np.interp(sB, s, q)) if s[-1] >= sB else float("nan")


def tail_slope(sq, frac=0.2):
    """dq/d(s/B) over the last `frac` of the reached s/B, normalised by q_end/s_end."""
    s, q = sq
    m = s >= s[-1] * (1 - frac)
    if m.sum() < 3:
        return float("nan"), float("nan")
    k = np.polyfit(s[m], q[m], 1)[0]
    return k, k / (q[-1] / s[-1])


D = {l: load(l) for l in legs + sorted(LOCAL)}
out = []
for name, lad in (("presidual", PRES), ("einit", EINI), ("ablation", ABL), ("A0_h0", AH)):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for l, v in lad:
        sq = D[l]
        if sq is None:
            continue
        ax.plot(sq[0], sq[1], lw=2.2 if l == "E_B" else 1.4,
                color="k" if l == "E_B" else None,
                label=f"{l} ({ {'presidual': 'Presidual = ', 'einit': 'e_init = '}.get(name, '')}{v})")
    ax.set_xlabel("s/B"); ax.set_ylabel("q (kPa)")
    ax.set_title(f"WP-138 strip footing, B/8, SAS-ME 129: {name} ladder (E_B = reference)")
    ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(SCR, f"q_s_{name}.png"), dpi=130)
    out.append(f"\n=== {name} ladder")
    out.append("leg         val     s/B_end   q_end    q_max   q@.01   q@.02   q@.03   q@.05  tailslope(norm) mode  steps fails  sB_1stNonPosH  refusals")
    for l, v in lad:
        sq, sm = D[l], status(l)
        if sq is None:
            out.append(f"{l:11s} {v:<7} (no data)"); continue
        ts = tail_slope(sq)[1]
        out.append(f"{l:11s} {v:<7} {sq[0][-1]:.5f} {sq[1][-1]:8.1f} {sq[1].max():8.1f} "
                   + " ".join(f"{q_at(sq, x):7.1f}" for x in (0.01, 0.02, 0.03, 0.05))
                   + f"  {ts:7.3f}  {sm['mode'] if sm else ('LOCAL' if l in LOCAL else 'RUNNING'):8s}"
                   + (f" {sm['steps']:5d} {sm['failed_attempts']:4d}" if sm else "           ")
                   + f"  {first_nonposh(l):8.5f}"
                   + f"  {(sm or {}).get('refusals_by_code') or refusals(l)}")
txt = "\n".join(out)
print(txt)
open(os.path.join(SCR, "ladder_table.txt"), "w").write(txt)
