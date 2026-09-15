"""ADR-94 addendum / F8 -- the ASSOCIATED (psi = phi) leg of the ADR-95 R3
Prandtl-Reissner collapse gate, run on BOTH of the fork's Drucker-Pragers.

WHY
---
The TIMs act measured, on its own plane-strain strip, that `ASDPlasticMaterial3D`
Drucker-Prager tracks the repaired UW `DruckerPrager` to within 1 % at psi = 0
and WALLS at s/B = 0.005-0.006 while STILL HARDENING when the flow rule is made
ASSOCIATED, where the UW material runs clean to s/B = 0.15.  This driver runs
the gate's OWN associated control (`tests/test_r3_prandtl_collapse_gate.py`,
`CONTROL_H0` leg) on the ASD material, so the comparison is made on the deck the
gate actually pins rather than on a second deck.

The gate's leg machinery is IMPORTED, not copied: `_run_leg` grew one
`material=` argument and one read-only `gp_probe=` hook for exactly this, and
`material="UW"` is its default, so every gate leg is unchanged.

INSTRUMENT
----------
Once per converged step (every `--census-every`-th), every Gauss point's
committed stress is classified in the (p, sqrt(J2)) half-plane of the ASD cone
    f = sqrt(J2) + eta*p - xi_c,     p_apex = xi_c/eta
into the counts that discriminate the candidate causes:
  * `n_tension`     : p > 0                       -- is the footing-edge tensile
                                                     spot present at all?
  * `n_over_apex`   : p > p_apex                  -- committed BEYOND the apex
                                                     (inadmissible; should be 0)
  * `n_apex_pinned` : |p - p_apex| tiny AND q ~ 0 -- took the APEX PROJECTION
  * `n_eucl_wedge`  : eta*q <= p - p_apex < (K*etabar/G)*q
        DIAGNOSTIC ONLY, AND IT CANNOT SEE THE DEFECT -- kept because reading a
        column of zeros and knowing WHY is worth more than not having asked.
        The misclassification happens on the TRIAL state, which is not visible
        from outside the material; this band is evaluated on COMMITTED stresses,
        and a committed state is either ON the cone (p - p_apex < 0, outside the
        band by construction) or AT the apex (p - p_apex = 0, q = 0, on the
        band's boundary).  So the only thing that can ever land in it is
        round-off around the apex point -- which is exactly what the two runs
        show: UW reads 125-137 and ASD reads 0, and that difference is the two
        materials' apex round-off, NOT a mechanical difference between them.
        What DOES measure the defect's footprint is `n_apex_pinned` at a station
        where the cone return was available.

RUN (each leg is one process; ~10-40 min at h0 = 1.0)
    python3.12 -u r3_assoc_probe.py --material ASD --h0 1.0 --assoc 1 \
        --out <dir> > asd_assoc_h1.0.log 2>&1
"""
import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, os.pardir, os.pardir, os.pardir))
TESTS = os.path.join(ROOT, "tests")


def _bind(dist):
    """Bind the engine BEFORE the gate module imports `opensees`."""
    sys.path.insert(0, TESTS)
    from _engine import bind_worktree_engine          # noqa: E402
    return bind_worktree_engine(dist)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dist", required=True, help="dist/bin holding opensees.pyd")
    ap.add_argument("--material", default="ASD", choices=("UW", "ASD"))
    ap.add_argument("--h0", type=float, default=1.0)
    ap.add_argument("--assoc", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--census-every", type=int, default=5)
    ap.add_argument("--wall", type=float, default=7200.0)
    a = ap.parse_args()

    ops = _bind(os.path.abspath(a.dist))
    print(f"[probe] ladrunoBuild() = {ops.ladrunoBuild()}", flush=True)
    print(f"[probe] opensees.pyd   = {os.path.abspath(ops.__file__)}", flush=True)

    import test_r3_prandtl_collapse_gate as G          # noqa: E402

    os.makedirs(a.out, exist_ok=True)
    assoc = bool(a.assoc)
    alpha = G._alpha_from_phi_txc(G.PHI_TXC)
    eta = 3.0 * alpha
    xi_c = G.SY / math.sqrt(3.0)
    p_apex = xi_c / eta
    g_el = 3.0 * G.K_EL * (1.0 - 2.0 * G.NU) / (2.0 * (1.0 + G.NU))
    etabar = eta if assoc else 0.0
    slope_exact = G.K_EL * etabar / g_el               # the elastic-metric slope
    tagname = (f"{a.material.lower()}_h{a.h0}_"
               f"{'assoc' if assoc else 'nonassoc'}")
    print(f"[probe] material={a.material} h0={a.h0} assoc={assoc}", flush=True)
    print(f"[probe] eta={eta:.6f} xi_c={xi_c:.6f} p_apex={p_apex:.6f} "
          f"K/G={G.K_EL/g_el:.4f} euclidean_slope={eta:.4f} "
          f"elastic_metric_slope={slope_exact:.4f}", flush=True)

    cpath = os.path.join(a.out, f"census_{tagname}.csv")
    ch = open(cpath, "w", newline="")
    cw = csv.writer(ch)
    cw.writerow(["step", "s_over_B", "q_kPa", "n_gp", "n_yield", "n_tension",
                 "n_over_apex", "n_apex_pinned", "n_eucl_wedge",
                 "p_max", "p_min", "q_at_pmax"])

    tol_p = 1.0e-9 * max(1.0, abs(p_apex))

    def probe(step, s_over_b, q_foot, nele):
        if step % a.census_every:
            return
        sig = np.concatenate([np.asarray(ops.eleResponse(e, "stress"),
                                         dtype=float).reshape(-1, 6)
                              for e in range(1, nele + 1)])
        p = sig[:, :3].mean(axis=1)
        dv = sig[:, :3] - p[:, None]
        j2 = 0.5 * (dv ** 2).sum(axis=1) + (sig[:, 3:6] ** 2).sum(axis=1)
        q = np.sqrt(np.maximum(j2, 0.0))
        f = q + eta * p - xi_c
        n_yield = int((f >= -1.0e-6 * xi_c).sum())
        n_tension = int((p > 0.0).sum())
        n_over = int((p > p_apex + tol_p).sum())
        pinned = (np.abs(p - p_apex) <= 1.0e-7 * max(1.0, abs(p_apex))) & \
                 (q <= 1.0e-7 * xi_c)
        wedge = ((p - p_apex) >= eta * q) & ((p - p_apex) < slope_exact * q)
        imax = int(np.argmax(p))
        cw.writerow([step, f"{s_over_b:.9g}", f"{q_foot:.9g}", len(p), n_yield,
                     n_tension, n_over, int(pinned.sum()), int(wedge.sum()),
                     f"{p.max():.9g}", f"{p.min():.9g}", f"{q[imax]:.9g}"])
        ch.flush()
        if step % (20 * a.census_every) == 0:
            print(f"[census] step {step} s/B {s_over_b:.5f} q {q_foot:8.3f} "
                  f"yield {n_yield} tens {n_tension} overapex {n_over} "
                  f"apexpin {int(pinned.sum())} wedge {int(wedge.sum())} "
                  f"pmax {p.max():.4f}", flush=True)

    t0 = time.time()
    r = G._run_leg(a.h0, assoc, a.out, wall_budget=a.wall,
                   material=a.material, gp_probe=probe)
    ch.close()
    r["wall_total_s"] = time.time() - t0
    r["build"] = str(ops.ladrunoBuild())
    r["census_csv"] = cpath
    print("\n[RESULT] " + json.dumps({k: v for k, v in r.items()
                                      if not isinstance(v, np.ndarray)}),
          flush=True)
    print(f"\n[probe] {r['tag']}: ratio {r['ratio']:.4f} q_max {r['qmax']:.3f} "
          f"mode {r['mode']} plateau {r['plateau']} free {r['free']} "
          f"capacity {r['capacity']} tail {r['tail_pct']:.3f} % "
          f"s_end/B {r['s_end_over_B']:.5f} steps {r['steps']} "
          f"nfail {r['nfail']} nsub {r['nsub']} wall {r['wall_s']:.0f} s",
          flush=True)
    print(f"[probe] verdict: {r['verdict']}", flush=True)


if __name__ == "__main__":
    main()
