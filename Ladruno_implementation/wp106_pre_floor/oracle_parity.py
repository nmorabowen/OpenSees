"""WP-106 / ADR-93 II.1 -- does the C++ `-pRe` compute what the oracle computes?

The ADR-93 numbers for candidate II.1 were all produced by the numpy oracle
(`adr92_p0_oracle/sanisand_implex_oracle.py`, whose `Presidual_e` seam predates
this work package). Those numbers only transfer to the fork if the C++ implements
the SAME floor. This is the ADR-92 G0 gate with one argument added: replay the
binary's OWN recorded strain sequence through the oracle at the SAME `pRe`, and
report the worst relative departure in committed stress / e / alpha / z / epsE.

    python3.12 oracle_parity.py data/tx_p20_pre1.csv data/tx_p0.5_pre1.csv ...

G0's bar is 1e-8 raw on a path the oracle reproduces; the same bar is used here.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "adr92_p0_oracle"))

from sanisand_implex_oracle import (  # noqa: E402
    CONSTS, committed_from_row, make_material, norm_contr, read_probe_csv,
    seed_from_csv,
)


def mat_from_meta(meta):
    consts = dict(CONSTS)
    consts["e_init"] = float(meta["e_init"])
    rec = meta.get("consts")
    if rec:
        import ast
        try:
            rec = ast.literal_eval(rec)
        except (ValueError, SyntaxError):
            rec = None
        if rec:
            for k, v in rec.items():
                if k in consts:
                    consts[k] = float(v)
    return make_material(consts,
                         scheme=int(meta["scheme"]),
                         Pmin=float(meta["Pmin"]),
                         Presidual=float(meta["Presidual"]),
                         Presidual_e=float(meta.get("pRe", 0.0)),
                         TolF=float(meta["tol"]), TolR=float(meta["tol"]),
                         honor_tolR=bool(int(meta["honorTolR"])))


def replay(path):
    meta, rows = read_probe_csv(path)
    if len(rows) < 3:
        return meta, None
    mat = mat_from_meta(meta)
    seed_from_csv(mat, rows[0])
    worst = dict(sig=0.0, e=0.0, alpha=0.0, z=0.0, epsE=0.0, step=0)
    for r in rows[1:]:
        eps = committed_from_row(r)[0]
        mat.integrate(eps)
        mat.commit()
        _, sig, epsE, alpha, fabric, _, e_ref, _ = committed_from_row(r)
        ds = norm_contr(mat.sig_n - sig) / norm_contr(sig)
        if ds > worst["sig"]:
            worst.update(sig=ds, step=int(r["step"]))
        worst["e"] = max(worst["e"], abs(mat.void_ratio - e_ref) / abs(e_ref))
        na = max(norm_contr(alpha), 1e-30)
        worst["alpha"] = max(worst["alpha"], norm_contr(mat.alpha_n - alpha) / na)
        nz = norm_contr(fabric)
        worst["z"] = max(worst["z"], (norm_contr(mat.fabric_n - fabric) / nz)
                         if nz > 1e-30 else norm_contr(mat.fabric_n - fabric))
        worst["epsE"] = max(worst["epsE"],
                            norm_contr(mat.epsE_n - epsE) / norm_contr(epsE))
    return meta, worst


def main():
    print(f"{'file':<28}{'pRe':>6}{'n':>5}{'sigma':>11}{'e':>11}{'alpha':>11}"
          f"{'z':>11}{'epsE':>11}  raw<=1e-8")
    for path in sys.argv[1:]:
        meta, w = replay(path)
        if w is None:
            print(f"{os.path.basename(path):<28} too few rows -- SKIPPED")
            continue
        raw = max(w["sig"], w["e"], w["alpha"], w["z"], w["epsE"])
        _, rows = read_probe_csv(path)
        print(f"{os.path.basename(path):<28}{float(meta.get('pRe', 0)):>6g}"
              f"{len(rows) - 1:>5}"
              f"{w['sig']:>11.2e}{w['e']:>11.2e}{w['alpha']:>11.2e}"
              f"{w['z']:>11.2e}{w['epsE']:>11.2e}  "
              f"{'yes' if raw <= 1e-8 else 'NO (' + f'{raw:.1e}' + ')'}")


if __name__ == "__main__":
    main()
