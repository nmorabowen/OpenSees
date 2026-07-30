"""Field output for the ADR-79 bearing benchmark — deformation and stresses.

The backbone runs recorded only the footing reaction, which answers "how much
load" but shows nothing about WHERE the soil is working. This re-runs the
`-geom hypo` leg and dumps whole-field snapshots at settlement checkpoints so
the mechanism can be looked at:

  * nodal displacements and pore pressure -> the deformed shape and the
    surface heave profile, which distinguishes a forming Prandtl mechanism
    (which heaves beside the footing) from pure punching (which does not);
  * element-averaged effective stress -> mean stress p', and
  * element-averaged MOBILISED FRICTION sin(phi_mob) = (s1-s3)/(-(s1+s3)),
    which is the diagnostic that says whether the soil is at its strength.
    The single-element probe measured the failure value at 0.770, so 0.770 and
    not sin(33 deg) is the line to compare against.
  * per-GP J from `hypoState` -> how much genuine volume change there is.

Stress and hypoState come back per Gauss point (8 x 6 and 8 x 3 on an H8) and
are averaged to one value per element.

Run (base Python 3.12, ~2-3 h to s/B = 0.05):
    python field_run.py [target_sB]
Writes fields_hypo.npz.
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bearing_backbone as bb  # noqa: E402

ops = bb.ops
HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = float(sys.argv[1]) if len(sys.argv) > 1 else 0.05
CHECK = [c for c in (0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10)
         if c <= TARGET + 1e-12]
SIN_FAIL = 0.770          # measured failure value, element_probe.py


def fields(nodes, nel):
    """Element-averaged stress, sin(phi_mob) and J; plus the nodal state."""
    sig = np.zeros((nel, 6))
    sphi = np.zeros(nel)
    Jbar = np.zeros(nel)
    for e in range(1, nel + 1):
        s = np.array(ops.eleResponse(e, "stress")).reshape(-1, 6)   # 8 GP x 6
        sig[e - 1] = s.mean(axis=0)
        # per-GP principal stresses -> per-GP mobilised friction, then average
        vals = []
        for row in s:
            xx, yy, zz, xy, yz, zx = row
            T = np.array([[xx, xy, zx], [xy, yy, yz], [zx, yz, zz]])
            w = np.sort(np.linalg.eigvalsh(T))[::-1]                # s1>=s3
            den = -(w[0] + w[2])
            # sin(phi_mob) = (s1-s3)/-(s1+s3) is only meaningful for a
            # CONFINED, NON-TENSILE state, and near the footing edge neither
            # holds. Two separate pathologies, both met in practice here:
            #   den -> 0 at the free surface, where the ratio explodes; and
            #   s1 > 0 (tension), where the numerator can exceed the
            #     denominator outright and the "sine" runs past 1.
            # A first version guarded only the first and still reported
            # "sphi max = 41.8", which is impossible for a sine. Mask both so
            # the plots show real mobilisation instead of edge artifacts.
            ok = den > 2.0 and w[0] < 0.5
            vals.append((w[0] - w[2]) / den if ok else np.nan)
        sphi[e - 1] = np.nanmean(vals) if np.any(np.isfinite(vals)) else np.nan
        try:
            hs = np.array(ops.eleResponse(e, "hypoState")).reshape(-1, 3)
            Jbar[e - 1] = hs[:, 0].mean()
        except Exception:
            Jbar[e - 1] = np.nan
    u = np.array([[ops.nodeDisp(i + 1, d) for d in (1, 2, 3, 4)]
                  for i in range(len(nodes))])
    return sig, sphi, Jbar, u


def main():
    nodes, hexes, trib, sets = bb.load_mesh()
    bb.build("hypo", False, nodes, hexes, trib, sets)
    bb.log(f"field run: {len(nodes)} nodes, {len(hexes)} hex, target "
           f"s/B={TARGET}, checkpoints {CHECK}")
    bb.gravity("field")
    nel = len(hexes)

    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = bb.Q_SURCH * float(trib[sets["footing"]].sum())
    uz0 = ops.nodeDisp(foot[0], 3)
    ops.loadConst("-time", 0.0)
    smax = TARGET * bb.B_FOOT
    rate = smax / bb.T_PUSH
    ops.timeSeries("Path", 2, "-time", 0.0, bb.T_PUSH,
                   "-values", uz0, uz0 - smax, "-useLast")
    ops.pattern("Plain", 2, 2)
    for tt in foot:
        ops.sp(tt, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    bb.pick_system()
    ops.test("NormDispIncr", 1e-7, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    out = {"nodes": nodes, "hexes": hexes, "sin_fail": SIN_FAIL}
    # gravity state, before any footing load — the datum for every field
    s0, p0, J0, u0 = fields(nodes, nel)
    out["sig_grav"], out["sphi_grav"], out["J_grav"], out["u_grav"] = s0, p0, J0, u0

    ds, good, nfail, done = bb.DS_BASE, 0, 0, []
    t0 = time.time()
    while True:
        s_now = ops.getTime() * rate
        if s_now >= smax - 1e-12:
            break
        nxt = next((c for c in CHECK if c * bb.B_FOOT > s_now + 1e-12), None)
        ds = min(ds, smax - s_now)
        if nxt is not None:                     # never step past a checkpoint
            ds = min(ds, nxt * bb.B_FOOT - s_now)
        if ops.analyze(1, ds / rate) != 0:
            nfail += 1
            good = 0
            ds *= 0.5
            if ds < bb.DS_MIN:
                bb.log(f"field run: NO CONVERGENCE at s/B={s_now / bb.B_FOOT:.4f}")
                break
            continue
        good += 1
        if good >= bb.GROW_AFTER and ds < bb.DS_MAX:
            ds, good = min(2.0 * ds, bb.DS_MAX), 0
        s = ops.getTime() * rate
        for c in CHECK:
            if c in done or s < c * bb.B_FOOT - 1e-9:
                continue
            ops.reactions()
            R = -sum(ops.nodeReaction(tt, 3) for tt in foot) + r_corr
            sig, sphi, Jb, u = fields(nodes, nel)
            tag = f"{c:.4f}"
            out[f"sig_{tag}"], out[f"sphi_{tag}"] = sig, sphi
            out[f"J_{tag}"], out[f"u_{tag}"] = Jb, u
            out[f"q_{tag}"] = R / bb.A_FOOT
            done.append(c)
            bb.log(f"field run: s/B={c:.4f} q={R / bb.A_FOOT:.1f} kPa "
                   f"sphi max={np.nanmax(sphi):.3f} "
                   f"(failure {SIN_FAIL}) frac>=0.95*fail="
                   f"{np.nanmean(sphi >= 0.95 * SIN_FAIL):.3f} "
                   f"[{time.time() - t0:.0f}s, {nfail} retries]")
            np.savez_compressed(os.path.join(HERE, "fields_hypo.npz"),
                                checkpoints=np.array(done), **out)
    bb.log(f"field run done: checkpoints {done}, {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
