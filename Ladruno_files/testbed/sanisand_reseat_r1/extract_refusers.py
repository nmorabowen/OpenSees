"""One-off: extract the committed wall states of the loadingNonPosH refusers from the
WP-138 Esmeralda field checkpoints (orchestrator's analysis/ck/<leg>/ckpt/*.npz,
not in the repo) into data/refuser_states.csv, so the testbed runs without them.

Columns (compression POSITIVE, tensor components xx yy zz xy yz zx, kPa):
leg, k, element, gp, step, s_over_B, gx, gy, p, psi, e, sigma_0..5, alpha_0..5,
alpha_in_0..5, z_0..5, deps_0..5 (the GP's last converged strain increment from
the previous checkpoint, Voigt ENGINEERING shear, compression positive; blank
where no earlier checkpoint was used).  sigma_zz is recovered from psi and e as
the WP-138 deck does (p_r = 0)."""
from __future__ import annotations

import csv
import os

import numpy as np

import r1common as C

PREV = {"E_B": "field_step00375.npz", "E_D": "field_step00385.npz", "E_B16": "field_step00070.npz"}
EL = {7516: (1880, 1), 7512: (1879, 1), 7844: (1962, 1), 8228: (2058, 1), 31279: (7820, 4)}


def main():
    os.makedirs(os.path.join(C.HERE, "data"), exist_ok=True)
    rows = []
    for (leg, k, el, gp, _prev) in C.REFUSERS:
        st, info = C.ckpt_state_npz(leg, k)
        de = C.last_increment(leg, k, PREV[leg])
        r = dict(leg=leg, k=k, element=el, gp=gp, step=info["step"], s_over_B=info["s_over_B"],
                 gx=info["gx"], gy=info["gy"], p=info["p"], psi=info["psi"], e=st.e)
        # the RAW vectors State.from_voigt was built from (pre-projection), so the
        # CSV path rebuilds a bit-identical State
        for name in ("sigma", "alpha", "alpha_in", "z"):
            for i, v in enumerate(info["raw"][name]):
                r[f"{name}_{i}"] = repr(float(v))
        for i, v in enumerate(de):
            r[f"deps_{i}"] = repr(float(v))
        rows.append(r)
    p = os.path.join(C.HERE, "data", "refuser_states.csv")
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("wrote", p, len(rows), "rows")


if __name__ == "__main__":
    main()
