"""Quantify the SPATIAL OSCILLATION of the mobilisation field.

Reading the ASCII maps in `quad_path_summary.py`, the linear B-bar field at a
given settlement is smooth and saturated while the H20 `uri` field at the SAME
settlement is visibly mottled -- adjacent elements alternating between mob 0.99
and 0.90.  Eyeballing a checkerboard is exactly the kind of claim this campaign
has had to retract twice, so it gets a number.

THE MEASURE.  For every element with 4 in-plane neighbours, the checkerboard
residual is

    r_e = mob_e - mean(mob over its 4 face neighbours)

which annihilates any linear field exactly (it is a discrete Laplacian up to
sign and scaling) and is maximal for a +/-1 alternating field.  Reported as
RMS(r) over the PLASTIC CORE (mob >= 0.7 and inside the fine region), divided by
the RMS of the field's own spatial variation over the same set,

    OSC = RMS(r_e) / RMS(mob_e - mean(mob))

so a smooth field trends to 0 and a field whose energy is all at the grid
Nyquist frequency trends to O(1).  Scale-free, so it can be compared between
elements and resolutions.

WHAT WOULD HAVE TO BE TRUE FOR A HIGH OSC TO NOT MEAN A CHECKERBOARD?  Two
things, and both are checked/reported:
  (i) a genuinely LOCALISED slip band is also a high-frequency feature, so it
      would raise OSC honestly.  That is why OSC is reported next to the
      localisation metrics -- a band is high-OSC AND thin AND connected; a
      checkerboard is high-OSC AND space-filling.  The measured fields here are
      space-filling (18-23 elements thick), so a band reading is excluded.
  (ii) the graded mesh makes neighbouring elements different SIZES near the
      grading front, which puts real signal into the neighbour difference.  So
      the statistic is restricted to the CORE, where the grid is uniform at h0,
      and the core mask is reported with it.

RESULT: THE CHECKERBOARD HYPOTHESIS IS **NOT SUPPORTED**, AND THIS SCRIPT IS THE
RECORD OF THAT NEGATIVE.  Measured (note 82 §6.3):

  * at h0 = 1.0 the measure looks decisive -- the linear control's plastic core
    is perfectly smooth (RMS residual 0.0000) against the H20 `uri`'s 0.0395;
  * at h0 = 0.5, the FINER mesh where the statistic is better resolved, it
    REVERSES: the linear control reads 0.0398 and the H20 `uri` 0.0203, i.e.
    the element that completes the collapse is the MORE "oscillatory" one.

Both readings have the same cause and it is not a checkerboard: the residual
cannot tell high-frequency NOISE from genuine high-CURVATURE structure, and the
plastic FRONT is strongly curved.  Where the linear core is fully saturated
(every element at mob = 1.0) the residual is identically zero because the field
is CONSTANT, not because it is smooth; where it is still advancing, the front
itself supplies the residual.  The `OSC` ratio additionally degenerates to 0/0
on a saturated field, which is why the raw RMS residual and the field RMS are
both printed and the ratio must not be read alone.

The visible mottling in `quad_path_summary.field_map` that prompted this is a
RENDERING artefact: the map thresholds at 0.99/0.90, so a core sitting at
0.985-0.995 -- a 1 % variation -- alternates between '#' and '+' and reads as a
checkerboard to the eye.  Measured neighbour-to-neighbour variation on the H20
`uri` core is ~2 % of the mean.

CAVEAT that keeps the question open rather than closed.  The classical
checkerboard test for a near-incompressible element is on the PRESSURE field,
not on a deviatoric mobilisation.  Mobilisation sees pressure only through the
I1 in its denominator, so this is an indirect probe and a genuine pressure
checkerboard could hide behind it.  Dumping I1 per Gauss point and repeating the
measure on it is the experiment that would settle it; note 82 §8 names it rather
than assuming its outcome.

Run:  py -3.12 quad_path_osc.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CORE_X, CORE_Z = 4.0, -4.0
PLASTIC = 0.70


def osc(tag):
    p = os.path.join(HERE, f"qpd_{tag}_field.npz")
    if not os.path.exists(p):
        return None
    d = np.load(p)
    mob, cen, ii, kk, x, z = (d["mob_mean"], d["cen"], d["ii"], d["kk"],
                              d["x"], d["z"])
    nx1, nz1 = len(x) - 1, len(z) - 1
    g = np.full((nx1, nz1), np.nan)
    g[ii, kk] = mob
    core = np.full((nx1, nz1), False)
    core[ii, kk] = ((np.abs(cen[:, 0]) <= CORE_X) & (cen[:, 2] >= CORE_Z)
                    & np.isfinite(mob) & (mob >= PLASTIC))
    sel = np.zeros_like(core)
    sel[1:-1, 1:-1] = core[1:-1, 1:-1]
    for di, dk in ((1, 0), (-1, 0), (0, 1), (0, -1)):     # need all 4 neighbours
        sh = np.roll(np.roll(np.isfinite(g), -di, 0), -dk, 1)
        sel &= sh
    if sel.sum() < 6:
        return None
    nb = np.zeros_like(g)
    for di, dk in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        nb += np.roll(np.roll(g, -di, 0), -dk, 1)
    r = g - nb / 4.0
    vals = g[sel]
    rms_r = float(np.sqrt((r[sel] ** 2).mean()))
    rms_f = float(np.sqrt(((vals - vals.mean()) ** 2).mean()))
    return dict(tag=tag, n=int(sel.sum()), rms_resid=rms_r, rms_field=rms_f,
                mob_mean=float(vals.mean()),
                # normalised by the field MEAN, which stays well-defined on a
                # saturated core where `rms_field` -> 0 and the ratio below
                # degenerates to 0/0.
                resid_pct=100.0 * rms_r / max(abs(float(vals.mean())), 1e-12),
                osc=rms_r / max(rms_f, 1e-12))


def main():
    out = []
    for p in sorted(glob.glob(os.path.join(HERE, "qpd_*.json"))):
        r = json.load(open(p))
        if "_smoke" in r["tag"] or "_elastic" in r["tag"]:
            continue
        o = osc(r["tag"])
        if o:
            o.update(h0=r["h0"], elem=r["elem"], mode=r["mode"],
                     s_end=r["s_end_over_B"], ratio=r["ratio"], M=r["M"],
                     thick=r["band_thick_el"])
            out.append(o)
    out.sort(key=lambda r: (r["h0"], r["elem"], r["s_end"]))
    h = (f"{'leg':>26} {'h0':>5} {'s/B end':>9} {'q/q_ex':>7} {'M':>7} "
         f"{'thick_el':>9} {'n_core':>7} {'mob_mean':>9} {'RMS resid':>10} "
         f"{'resid %':>8} {'RMS field':>10} {'OSC':>8}")
    print("SPATIAL-OSCILLATION (CHECKERBOARD) INDEX OF THE MOBILISATION FIELD")
    print("*** NEGATIVE RESULT: this does NOT discriminate. See the docstring. ***")
    print(h)
    print("-" * len(h))
    for r in out:
        print(f"{r['tag']:>26} {r['h0']:5.2f} {r['s_end']:9.5f} "
              f"{r['ratio']:7.4f} {r['M']:7.3f} {r['thick']:9.2f} "
              f"{r['n']:7d} {r['mob_mean']:9.4f} {r['rms_resid']:10.4f} "
              f"{r['resid_pct']:8.2f} {r['rms_field']:10.4f} {r['osc']:8.3f}")
    print("\n  resid % = RMS(mob - mean of 4 neighbours) / mean(mob), over the")
    print("  PLASTIC CORE only (mob >= 0.70, |x| <= 4 m, z >= -4 m, where the grid")
    print("  is uniform at h0).  OSC divides by the field's own RMS variation")
    print("  instead, and DEGENERATES to 0/0 on a saturated core -- read resid %.")
    print("\n  THE MEASURED VERDICT: at h0 = 1.0 the linear control reads 0.00 %")
    print("  against H20 uri's ~4 %, which looks decisive; at h0 = 0.5 it REVERSES")
    print("  (linear 4.3 %, H20 uri 2.1 %).  The residual cannot separate")
    print("  high-frequency NOISE from a high-CURVATURE plastic FRONT, and it is")
    print("  identically zero on a saturated core because that field is CONSTANT.")
    print("  So the checkerboard reading of the ASCII field maps is NOT supported;")
    print("  the mottling there is a thresholding artefact on a ~2 % variation.")


if __name__ == "__main__":
    main()
