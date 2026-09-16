"""WP-106 / ADR-93 II.1 -- the ring replay, measured against the pRe = 0 arm.

`adr93_p0/replay_ring.py --gate II1` measures every arm against the `p_r = 1.01`
VANILLA control, so its `d(sigma)` column mixes the strength floor's effect with
the elastic floor's. WP-106 needs the other difference: what `-pRe` alone does to
a leg that is otherwise the fork's own default (`p_r,p = 0`).

    python3.12 ring_delta.py

Reports, per leg and per `p_r,e in {0, 0.1, 1}` kPa:
  * substeps total / mean / max per step, and the ratio against `p_r,e = 0`
  * the committed-stress change against the `p_r,e = 0` arm (median and max
    relative), which is the "curve change" the certificate tolerance bounds
  * the same two numbers restricted to the 20 steps of LOWEST committed p, which
    is where the floor can act at all.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
P0 = os.path.join(os.path.dirname(HERE), "adr93_p0")
sys.path.insert(0, P0)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "adr92_p0_oracle"))

os.chdir(P0)   # replay_ring resolves `data/` off its own dir, but be explicit

from replay_ring import naive_seed, read_ring, replay  # noqa: E402
from sanisand_implex_oracle import norm_contr  # noqa: E402

ARMS = (0.0, 0.1, 1.0)


def main():
    rows_out = []
    for leg in ("dense", "gorini"):
        rows = read_ring(leg, "implex")
        seed = naive_seed(rows, leg)
        base = None
        bsub = None
        for pre in ARMS:
            _, out = replay(rows, leg, "implex", seed=seed, p_r=0.0, p_r_e=pre)
            sub = np.array([o["substeps"] for o in out], float)
            p = np.array([o["p"] for o in out], float)
            if base is None:
                base, bsub = out, sub
                dv = np.zeros(len(out))
            else:
                dv = np.array([norm_contr(o["sig"] - b["sig"]) / norm_contr(b["sig"])
                               for o, b in zip(out, base)])
            lo = np.argsort(p)[:20]          # the 20 lowest-p steps
            rows_out.append([
                leg, f"{pre:g}", int(sub.sum()), f"{sub.mean():.1f}", int(sub.max()),
                f"{bsub.sum() / max(sub.sum(), 1):.3f}x",
                f"{np.median(dv):.2e}", f"{dv.max():.2e}",
                f"{p.min():.3f}", f"{sub[lo].mean():.1f}",
                f"{bsub[lo].sum() / max(sub[lo].sum(), 1):.3f}x",
                f"{dv[lo].max():.2e}",
            ])
    head = ["leg", "p_r,e [kPa]", "substeps", "mean/step", "max/step",
            "substep cut vs pRe0", "median d(sig)", "max d(sig)",
            "min p [kPa]", "mean/step @20 lowest p", "cut @20 lowest p",
            "max d(sig) @20 lowest p"]
    print("| " + " | ".join(head) + " |")
    print("|" + "|".join("---" for _ in head) + "|")
    for r in rows_out:
        print("| " + " | ".join(str(x) for x in r) + " |")


if __name__ == "__main__":
    main()
