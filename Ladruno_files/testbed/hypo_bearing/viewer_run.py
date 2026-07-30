"""Bearing benchmark run that writes a `.ladruno` HDF5 for the apeGmsh viewer.

`field_run.py` dumps fields to npz for matplotlib sections. This instead writes
the fork's own `.ladruno` recorder file, which `apeGmsh.Results.from_ladruno`
reads DIRECTLY — that path is self-sufficient (the file carries its own
MODEL/NODES/ELEMENTS), so no `model.h5` and no node/element tag reconciliation
between the deck and a FEMData snapshot is needed. `from_mpco` would have
required `model_h5=`; `from_ladruno` does not.

Channels: nodal `displacement` (the deformed shape and the surface profile),
and per-Gauss-point `stress` (EFFECTIVE stress on LadrunoUP, response 1) plus
`porePressure`, from which the viewer can build p', deviatoric and
principal-stress fields.

Frames are throttled with `-T dt`: the adaptive stepper takes thousands of
steps and recording all of them would write per-GP tensors for 2816 elements
each time.

Run (base Python 3.12):
    python viewer_run.py [target_sB] [frames]
Writes bearing_view.ladruno .
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bearing_backbone as bb  # noqa: E402

ops = bb.ops
HERE = os.path.dirname(os.path.abspath(__file__))
TARGET = float(sys.argv[1]) if len(sys.argv) > 1 else 0.03
FRAMES = int(sys.argv[2]) if len(sys.argv) > 2 else 24
OUTFILE = os.path.join(HERE, sys.argv[3] if len(sys.argv) > 3
                       else "bearing_view.ladruno")


def main():
    nodes, hexes, trib, sets = bb.load_mesh()
    bb.build("hypo", False, nodes, hexes, trib, sets)
    bb.log(f"viewer run: {len(nodes)} nodes, {len(hexes)} hex, "
           f"target s/B={TARGET}, {FRAMES} frames -> {os.path.basename(OUTFILE)}")
    bb.gravity("viewer")

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

    # the recorder goes on AFTER gravity so the file holds the push only, and
    # frame 0 is the gravity state (the datum the mechanism is measured from).
    # `strain` is recorded alongside `stress` because PDMY exposes NO plastic
    # strain of its own — its setResponse offers only stress / strain / tangent
    # / backbone (PressureDependMultiYield.cpp:1306-1326). Equivalent plastic
    # strain therefore has to be DERIVED, by removing the elastic part from the
    # total strain using PDMY's own pressure-dependent moduli
    # G(p') = G_r (p'/p_ref)^d and K(p') = B_r (p'/p_ref)^d, which needs both
    # tensors at the same Gauss point. See make_gp_figs.py.
    ops.recorder("ladruno", OUTFILE,
                 "-N", "displacement",
                 "-E", "stress", "strain", "porePressure",
                 "-T", "dt", bb.T_PUSH / FRAMES)

    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    bb.pick_system()
    ops.test("NormDispIncr", 1e-7, 100, 0)
    ops.algorithm("KrylovNewton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    ds, good, nfail, n = bb.DS_BASE, 0, 0, 0
    t0 = time.time()
    while True:
        s_now = ops.getTime() * rate
        if s_now >= smax - 1e-12:
            break
        ds = min(ds, smax - s_now)
        if ops.analyze(1, ds / rate) != 0:
            nfail += 1
            good = 0
            ds *= 0.5
            if ds < bb.DS_MIN:
                bb.log(f"viewer run: NO CONVERGENCE at s/B={s_now / bb.B_FOOT:.4f}")
                break
            continue
        n += 1
        good += 1
        if good >= bb.GROW_AFTER and ds < bb.DS_MAX:
            ds, good = min(2.0 * ds, bb.DS_MAX), 0
        if n % 100 == 0:
            ops.reactions()
            R = -sum(ops.nodeReaction(tt, 3) for tt in foot) + r_corr
            bb.log(f"viewer run: s/B={ops.getTime() * rate / bb.B_FOOT:.4f} "
                   f"q={R / bb.A_FOOT:.1f} kPa [{n} steps, {nfail} retries, "
                   f"{time.time() - t0:.0f}s]")
    ops.reactions()
    R = -sum(ops.nodeReaction(tt, 3) for tt in foot) + r_corr
    # wipe closes and flushes the recorder; without it the HDF5 can be left
    # without its final steps.
    ops.remove("recorders")
    ops.wipe()
    bb.log(f"viewer run done: s/B={TARGET}, q={R / bb.A_FOOT:.1f} kPa, "
           f"{n} steps, {nfail} retries, {time.time() - t0:.0f}s")
    if os.path.exists(OUTFILE):
        bb.log(f"wrote {OUTFILE} ({os.path.getsize(OUTFILE) / 1e6:.1f} MB)")
    else:
        bb.log("WARNING: no .ladruno file was produced")


if __name__ == "__main__":
    main()
