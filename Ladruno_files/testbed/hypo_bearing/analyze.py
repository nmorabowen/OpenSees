"""ADR-79 bearing campaign — reduce the three backbones to the verdict.

Reads backbone_<leg>.csv, emits a checkpoint table, the cross-method ratios,
a limit-point / hardening diagnosis per leg, and bearing_backbone.png.
Truncation-honest: a leg that did not reach a checkpoint prints '--', and the
limit-point verdict is only issued on the span a leg actually covered.

Run with the apeGmsh env (it has matplotlib):
    C:/Users/nmb/venv/opensees_env/Scripts/python.exe analyze.py
"""
import csv
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
B_FOOT = 2.0
Q_VESIC = 637.5          # kPa — recomputed below from the runner's constants
LEGS = [("corot", "-geom corot"), ("hypo", "-geom hypo"),
        ("hypo_kc", "-geom hypo -kozenyCarman"),
        ("corot_std", "-geom corot -form std (rerun)"),
        ("corot_bbar", "-geom corot -form bbar"),
        ("corot_std_und", "corot std UNDRAINED"),
        ("corot_bbar_und", "corot bbar UNDRAINED")]
# The two levers this campaign can put on the SAME problem. GEOMETRY is the
# ladder rung (hypo vs corot, both `std`). LOCKING is the formulation
# (corot bbar vs corot std, both rotation-only). Under displacement control a
# softer element carries LESS q at the same settlement, so locking < 1 means
# bbar relieved locking and the std backbone was that much artificial stiffness.
# The locking ratio uses corot_std, NOT the archived corot: those two legs are
# the engine-matched pair. corot stays in the table as the §8 baseline.
RATIOS = [("hypo/corot (geometry)", "hypo", "corot"),
          ("bbar/std (locking)", "corot_bbar", "corot_std"),
          ("std rerun/archived", "corot_std", "corot"),
          ("bbar/std UNDRAINED", "corot_bbar_und", "corot_std_und")]
CHECKPOINTS = (0.01, 0.025, 0.05, 0.075, 0.10, 0.125, 0.15)


def load(leg):
    path = os.path.join(HERE, f"backbone_{leg}.csv")
    if not os.path.exists(path):
        return None
    rows = []
    with open(path) as f:
        for r in list(csv.reader(f))[1:]:
            try:                       # tolerate a torn final line mid-run
                rows.append([float(v) for v in r])
            except ValueError:
                continue
    if not rows:
        return None
    # legs run before the undrained work wrote 7 columns; p_probe_kPa is the
    # 8th and only exists on later runs. Tolerate both widths.
    width = min(len(r) for r in rows)
    a = np.array([r[:width] for r in rows])
    return dict(s=a[:, 0], q=a[:, 3], qn=a[:, 4], ds=a[:, 5], wall=a[:, 6],
                p=a[:, 7] if width > 7 else None)


def at(d, target):
    """q/q_Vesic at settlement `target`, or None if the leg never got there."""
    if d is None or d["s"][-1] < target - 1e-9:
        return None
    return float(np.interp(target, d["s"], d["qn"]))


def main():
    data = {leg: load(leg) for leg, _ in LEGS}
    have = {k: v for k, v in data.items() if v is not None}
    if not have:
        raise SystemExit("no backbones yet")

    print(f"q_Vesic = {Q_VESIC:.1f} kPa\n")
    print("| s/B | " + " | ".join(lbl for _, lbl in LEGS) + " | " +
          " | ".join(lbl for lbl, _, _ in RATIOS) + " |")
    print("|---|" + "---|" * (len(LEGS) + len(RATIOS)))
    for frac in CHECKPOINTS:
        tgt = frac * B_FOOT
        vals = {leg: at(data[leg], tgt) for leg, _ in LEGS}
        cells = [f"{vals[leg]:.2f}" if vals[leg] is not None else "--"
                 for leg, _ in LEGS]
        for _, num, den in RATIOS:
            cells.append(f"{vals[num] / vals[den]:.3f}"
                         if vals[num] and vals[den] else "--")
        print(f"| {frac:.3f} | " + " | ".join(cells) + " |")

    print("\nper-leg diagnosis (over the span each leg actually reached):")
    for leg, lbl in LEGS:
        d = data.get(leg)
        if d is None:
            print(f"  {lbl}: NO DATA")
            continue
        smaxr = d["s"][-1]
        complete = smaxr >= 0.15 * B_FOOT - 1e-9
        imax = int(np.argmax(d["q"]))
        qmax = d["q"][imax]
        # tangent stiffness dq/ds over the last quarter, normalized by the
        # initial one — the honest "is it still hardening" measure.
        n = len(d["s"])
        k0 = np.polyfit(d["s"][:max(n // 20, 2)], d["q"][:max(n // 20, 2)], 1)[0]
        kt = np.polyfit(d["s"][3 * n // 4:], d["q"][3 * n // 4:], 1)[0]
        softening = d["q"][-1] < 0.995 * qmax
        verdict = ("SOFTENING past a limit point" if softening else
                   "still hardening, NO limit point")
        print(f"  {lbl}: {'complete' if complete else 'TRUNCATED'} at "
              f"s/B={smaxr / B_FOOT:.4f}, {n} steps, "
              f"{d['wall'][-1] / 60:.0f} min")
        print(f"      q_end={d['q'][-1]:.0f} kPa ({d['qn'][-1]:.2f}x Vesic), "
              f"q_max={qmax:.0f} kPa at s/B={d['s'][imax] / B_FOOT:.4f}")
        print(f"      dq/ds: initial {k0:.3g} -> final {kt:.3g} kPa/m "
              f"({kt / k0 * 100:.1f}% of initial) => {verdict}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("\n(matplotlib unavailable — skipping the plot)")
        return
    # hypo and hypo+kc coincide to ~1e-4, so they MUST be drawn with different
    # line styles or one hides the other entirely.
    STYLE = {"corot": dict(color="C0", ls="-", lw=1.8),
             "hypo": dict(color="C3", ls="-", lw=1.8),
             "hypo_kc": dict(color="C2", ls="--", lw=1.8),
             "corot_std": dict(color="C0", ls=":", lw=1.4),
             "corot_bbar": dict(color="C4", ls="-.", lw=1.8),
             "corot_std_und": dict(color="C1", ls=":", lw=1.6),
             "corot_bbar_und": dict(color="C5", ls="-.", lw=1.6)}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for leg, lbl in LEGS:
        d = data.get(leg)
        if d is None:
            continue
        x = d["s"] / B_FOOT * 100
        ax[0].plot(x, d["q"], label=lbl, **STYLE[leg])
        ax[1].plot(x, d["qn"], label=lbl, **STYLE[leg])
        # tangent stiffness, smoothed over a sliding window (the raw
        # difference quotient is noisy because ds itself varies)
        w = 15
        if len(x) > 2 * w:
            k = np.gradient(d["q"], d["s"])
            ker = np.ones(w) / w
            ks = np.convolve(k, ker, mode="valid")
            ax[2].plot(x[w - 1:], ks / ks[0], label=lbl, **STYLE[leg])
        ax[0].plot(x[-1], d["q"][-1], "o", ms=5, color=STYLE[leg]["color"])
    ax[0].set_ylabel("footing pressure q [kPa]")
    ax[1].set_ylabel("q / q_Vesic")
    ax[2].set_ylabel("tangent dq/ds, normalized to its initial value")
    ax[1].axhline(1.0, color="k", ls=":", lw=0.9)
    ax[2].axhline(0.0, color="k", ls=":", lw=0.9)
    ax[2].set_ylim(bottom=0)
    ax[0].set_title("markers = last converged point (truncation)", fontsize=8)
    ax[2].set_title("a limit point would be dq/ds crossing zero", fontsize=8)
    for a in ax:
        a.set_xlabel("settlement s/B [%]")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
    fig.suptitle("ADR-79 bearing backbone — geometry ladder + the bbar locking "
                 "leg (B=2 m, saturated PDMY, 4.5B clearance)", fontsize=10)
    fig.tight_layout()
    out = os.path.join(HERE, "bearing_backbone.png")
    fig.savefig(out, dpi=140)
    print(f"\nplot -> {out}")


if __name__ == "__main__":
    main()
