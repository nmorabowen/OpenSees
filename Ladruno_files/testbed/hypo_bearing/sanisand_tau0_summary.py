"""Summarise an ADR-90 WP-A2 tau = 0 band campaign into the 3x2 table.

Reads every `a2_*.json` written by `sanisand_tau0_band.py` in a directory and
prints:

  * the provenance block (engine hash, driver, dates) -- refusing to summarise a
    directory whose legs do not all carry the SAME build hash, because a band
    assembled across two binaries is not a band;
  * the controls table (resultant identity, 1-D geostatic patch, eta/M_c at the
    stage flip, `Outside Bounding` and low-p `CLAMPING` counts, the
    base-vs-footing reaction mismatch);
  * the 3x2 q_u table with, per leg, the termination mode, the tail dq/ds, the
    steps and subdivisions against the pinned budget, and the wall time;
  * the width-vs-h numbers;
  * per density, the THREE-MESH BAND: max/min of q_u over the three meshes,
    reported as a percentage, next to the +/- 3 % per-resolution half-width the
    R3 DruckerPrager gate uses.

WHAT THE BAND MAY AND MAY NOT BE COMPARED TO
--------------------------------------------
The R3 gate's +/- 3 % is a per-resolution TOLERANCE ON A KNOWN CENTRE, not a
convergence criterion: R3's own three-mesh spread is 1.0849 -> 0.9514, i.e.
**14.0 %**, and that sequence is nonetheless accepted as converged-from-above
because it is monotone and every leg is a genuine plateau.  So the number this
script prints as "R3 comparison" is context, not a gate.

**The campaign's own tolerance is OQ2 and has not been supplied** (ADR 90 §7.5:
TIMs must set the target band width relative to B, the ramp duration, the De,
and the tolerance bands; as of this campaign none of that exists on either
side).  Any verdict of the form "inside the tolerance" is therefore UNAVAILABLE
by construction, and this script says so rather than substituting a number of
its own.

A leg whose termination mode is `WALL` or `FLOOR` is printed as INADMISSIBLE
and is EXCLUDED from the band, with the exclusion stated -- a seizure is not a
capacity, whatever its q_max says.

    python3.12 sanisand_tau0_summary.py <out_dir>
"""

from __future__ import annotations

import csv
import glob
import json
import math
import os
import sys

import numpy as np

# --------------------------------------------------------------------------
# The width metric lives HERE, not in the driver, for two reasons: this module
# imports no engine (so the metric can be re-run on an old campaign's CSVs on
# any box), and there is then exactly ONE implementation -- the driver imports
# these names.  Widths are ALWAYS recomputed from the field CSVs rather than
# read out of the leg JSON, so a campaign run before a probe-depth fix is
# re-reduced correctly instead of quoting a stale number.
# --------------------------------------------------------------------------

# Odd multiples of 0.125 m: strictly INTERIOR to one element layer on all three
# meshes (whose z-lines in the fine block sit at multiples of h0 = 1.0 / 0.5 /
# 0.25).  A round -0.5 or -1.0 lands ON a mesh line at h0 = 0.5 and 0.25, and a
# nearest-layer rule then ties between two layers -- i.e. the coarse mesh gets
# measured on one layer and the fine meshes on two, which is a mesh-dependent
# measuring rule inside a mesh-dependence study.
Z_PROBES = (-0.625, -1.375)
EPSQ_YIELD = 1.0e-5      # only for the yielding COUNT and VOLUME, never for w2


def w2_metric(x, hx, p):
    """ADR 90 §7.3 threshold-free width `sqrt(12*Var)`.

    `Var = sum p_e[(x_e - xbar)^2 + hx_e^2/12] / sum p_e`.  The `hx^2/12` term
    is the within-element variance of a piecewise-constant profile: it makes a
    one-element band read exactly `hx` and a k-element top hat read exactly
    `k*hx`, so the number is comparable across meshes.  Without it a
    one-element band reads 0.  Returns (nan, nan) on an all-zero profile.
    """
    x, hx, p = np.asarray(x), np.asarray(hx), np.asarray(p)
    tot = float(np.sum(p))
    if tot <= 0.0:
        return float("nan"), float("nan")
    xbar = float(np.sum(p * x) / tot)
    var = float(np.sum(p * ((x - xbar) ** 2 + hx ** 2 / 12.0)) / tot)
    return math.sqrt(12.0 * var), xbar


def widths_from_field(epsq, xc, zc, hx, hz, vol, z_probes=Z_PROBES):
    """w2 at each fixed physical depth, plus the yield extent.

    Row selection is SPAN CONTAINMENT -- the one element layer whose z-extent
    strictly contains the probe depth -- and the profile is restricted to
    x >= 0, because the mesh is symmetric and a two-sided profile would measure
    the footing width rather than the band width.
    """
    epsq, xc, zc = np.asarray(epsq), np.asarray(xc), np.asarray(zc)
    hx, hz, vol = np.asarray(hx), np.asarray(hz), np.asarray(vol)
    out = {}
    for zp in z_probes:
        row = (zc - 0.5 * hz < zp) & (zp < zc + 0.5 * hz) & (xc >= 0.0)
        ww, xb = w2_metric(xc[row], hx[row], epsq[row])
        out[f"w2_z{zp}"] = ww
        out[f"xbar_z{zp}"] = xb
        out[f"nrow_z{zp}"] = int(row.sum())
        out[f"epsqmax_z{zp}"] = float(epsq[row].max()) if row.sum() else 0.0
        out[f"hx_z{zp}"] = float(hx[row].min()) if row.sum() else float("nan")
    y = epsq > EPSQ_YIELD
    out["n_yield_ele"] = int(y.sum())
    out["vol_yield"] = float(vol[y].sum())
    out["epsq_max"] = float(epsq.max())
    return out


def read_field(path):
    """Read a driver field CSV (one `#` provenance line, then a header)."""
    with open(path) as f:
        f.readline()                                   # provenance comment
        rd = csv.DictReader(f)
        cols = {k: [] for k in ("xc", "zc", "hx", "hz", "vol", "eps_q_p")}
        for row in rd:
            for k in cols:
                cols[k].append(float(row[k]))
    return {k: np.array(v) for k, v in cols.items()}


def widths_from_csv(path, z_probes=Z_PROBES):
    c = read_field(path)
    return widths_from_field(c["eps_q_p"], c["xc"], c["zc"], c["hx"], c["hz"],
                             c["vol"], z_probes)


R3_BAND_HALFWIDTH_PCT = 3.0
R3_SEQUENCE = {1.0: 1.0849, 0.5: 0.9977, 0.25: 0.9514}   # measured, DruckerPrager
_ADMISSIBLE = ("TARGET", "PEAK", "BUDGET")


def load(out_dir):
    """Load the leg JSONs and RE-REDUCE every width from its field CSV.

    The driver writes widths into the JSON too, but this module recomputes them
    so that a campaign run before a probe-depth or metric change is reported
    correctly rather than quoting a stale number.  A leg whose field CSV is
    missing keeps its JSON value and is marked.
    """
    legs = []
    for p in sorted(glob.glob(os.path.join(out_dir, "a2_*.json"))):
        with open(p) as f:
            r = json.load(f)
        fld = r.get("field", "")
        if fld and not os.path.isabs(fld):
            fld = os.path.join(out_dir, fld)
        if fld and os.path.exists(fld):
            r.update(widths_from_csv(fld))
            r["_width_source"] = "recomputed from field CSV"
        else:
            r["_width_source"] = "JSON (field CSV missing)"
        for c in r.get("checkpoints", []):
            cp = os.path.join(out_dir,
                              f"a2_{r['tag']}_field_sB{c['cp_target']:g}.csv")
            if os.path.exists(cp):
                c.update(widths_from_csv(cp))
        legs.append(r)
    return legs


def selfcheck():
    """The width metric is CALIBRATED, not fitted: a k-element top hat must
    read exactly k*hx.  Run it before quoting any width."""
    ok = True
    for hx0 in (0.25, 1.0):
        for k in (1, 2, 3, 5, 8):
            x = np.arange(20) * hx0
            p = np.zeros(20)
            p[5:5 + k] = 1.0
            w, _ = w2_metric(x, np.full(20, hx0), p)
            good = abs(w - k * hx0) < 1e-12
            ok &= good
            print(f"  hx={hx0}  k={k}: w2={w:.12g}  expect {k * hx0:.12g}  "
                  f"{'OK' if good else 'FAIL'}")
    print("width metric calibration:", "OK" if ok else "FAILED")
    return 0 if ok else 1


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "--selfcheck":
        return selfcheck()
    out_dir = argv[0]
    legs = load(out_dir)
    if not legs:
        print(f"no a2_*.json in {out_dir}")
        return 1

    builds = sorted({r["build"] for r in legs})
    print("=" * 100)
    print("ADR 90 WP-A2 -- tau = 0 collapse-load band, strip footing on "
          "softening SANISAND")
    print("=" * 100)
    print(f"  engine build : {', '.join(builds)}")
    if len(builds) > 1:
        print("  *** REFUSING: legs span more than one build hash; a band "
              "assembled across two binaries is not a band ***")
        return 1
    print(f"  driver       : {legs[0]['driver']}")
    print(f"  output       : {os.path.abspath(out_dir)}")
    print(f"  dates        : {min(r['date'] for r in legs)} .. "
          f"{max(r['date'] for r in legs)}")
    print(f"  deck         : gamma={legs[0]['gamma']} kN/m3, K0="
          f"{legs[0]['K0']:.4f}, M_c={legs[0]['M_c']}, "
          f"-Presidual {legs[0]['presidual']}, -Pmin {legs[0]['pmin']}, "
          f"solver {legs[0]['solver']}")
    print(f"  controller   : SUBDIV_BUDGET {legs[0]['subdiv_budget']} (pinned, "
          f"R3), DS_MAX {legs[0]['ds_max']} m (pinned here), "
          f"target s/B {legs[0]['sfrac_target']}, wall "
          f"{legs[0]['wall_budget_s']:.0f} s/leg")

    print("\n--- CONTROLS (a red cell invalidates the leg, not just the row) ---")
    hdr = (f"{'leg':>16} {'resultant':>11} {'geostat patch':>14} "
           f"{'eta/M_c':>8} {'OutBnd':>7} {'CLAMP':>6} {'base-foot':>11} "
           f"{'p_min kPa':>10}")
    print(hdr)
    print("-" * len(hdr))
    for r in legs:
        print(f"{r['tag']:>16} {r['resultant_err']:11.2e} "
              f"{r['patch_err']:14.2e} {r['eta_over_Mc']:8.4f} "
              f"{r['n_outside_bounding']:7d} {r['n_clamping']:6d} "
              f"{r['base_foot_mismatch']:11.2e} {r['p_min_grav']:10.4f}")
    print("  (geostatic patch is also the GRAVITY-STATE MESH-SENSITIVITY "
          "measurement: it compares every Gauss point against a\n"
          "   mesh-independent closed form, so its spread across the three "
          "meshes IS the sensitivity.  Requirement: <= 0.1 %.)")

    print("\n--- q_u TABLE (3 meshes x 2 densities) ---")
    hdr = (f"{'leg':>16} {'DOF':>7} {'q_u kPa':>10} {'s_pk/B':>8} "
           f"{'tail %':>10} {'plat':>5} {'peak':>5} {'mode':>7} "
           f"{'ds/floor':>9} {'CAP':>4} {'s_end/B':>8} {'steps':>6} "
           f"{'sub/80':>7} {'wall s':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in legs:
        print(f"{r['tag']:>16} {r['dof']:>7} {r['q_u']:10.2f} "
              f"{r['s_peak_over_B']:8.4f} {r['tail_pct']:10.3f} "
              f"{'yes' if r['plateau'] else 'NO':>5} "
              f"{'yes' if r['peaked'] else 'NO':>5} {r['mode']:>7} "
              f"{r['headroom']:9.1f} {'yes' if r['capacity'] else 'NO':>4} "
              f"{r['s_end_over_B']:8.4f} {r['steps']:6d} "
              f"{r['nsub']:3d}/{r['subdiv_budget']:<3d} {r['wall_s']:8.0f}")

    print("\n--- TERMINATION STORIES ---")
    for r in legs:
        print(f"  [{r['tag']}] {r['mode']}: {r['verdict']}")

    print(f"\n--- WIDTH vs h at the END OF EACH LEG (ADR 90 §7.3 "
          f"threshold-free w2, probe depths {Z_PROBES} m, x >= 0 half) ---")
    print(f"  width source: {legs[0]['_width_source']}")
    hdr = (f"{'leg':>16} {'h0':>6} {'w2(z1)':>11} {'w2/h0':>8} "
           f"{'w2(z2)':>11} {'w2/h0':>8} {'yield ele':>10} "
           f"{'yield vol':>10} {'epsq max':>10}")
    print(hdr)
    print("-" * len(hdr))
    for r in legs:
        w05 = r.get(f"w2_z{Z_PROBES[0]}")
        w10 = r.get(f"w2_z{Z_PROBES[1]}")
        print(f"{r['tag']:>16} {r['h0']:6.2f} "
              f"{w05 if w05 == w05 else float('nan'):11.4f} "
              f"{(w05 / r['h0']) if w05 == w05 else float('nan'):8.2f} "
              f"{w10 if w10 == w10 else float('nan'):11.4f} "
              f"{(w10 / r['h0']) if w10 == w10 else float('nan'):8.2f} "
              f"{r['n_yield_ele']:10d} {r['vol_yield']:10.3f} "
              f"{r['epsq_max']:10.3e}")
    print("  w2/h0 constant across the three meshes == the band is ONE ELEMENT "
          "wide == the classical ill-posedness signature.\n"
          "  The element COUNT is not a mesh-convergent quantity; the yielding "
          "VOLUME is.")

    print("\n--- MATCHED-SETTLEMENT COMPARISON (the deliverable that survives a "
          "leg being stopped by the clock) ---")
    print("  Legs on different meshes do not terminate at the same s/B, so "
          "each leg snapshots q and the width at the\n"
          "  SAME list of s/B values.  A three-mesh spread here is a genuine "
          "h-convergence statement at fixed settlement;\n"
          "  it is NOT a collapse load unless the settlement is past the peak.")
    for ename in sorted({r["e_name"] for r in legs}):
        grp = sorted([r for r in legs if r["e_name"] == ename],
                     key=lambda r: -r["h0"])
        cps = sorted({c["cp_target"] for r in grp
                      for c in r.get("checkpoints", [])})
        if not cps:
            print(f"\n  density {ename}: no checkpoint reached")
            continue
        print(f"\n  density {ename}")
        hdr = "    " + f"{'s/B':>8}" + "".join(
            f"{'q h' + str(r['h0']):>12}" for r in grp) + f"{'band %':>9}" \
            + "".join(f"{'w2 h' + str(r['h0']):>12}" for r in grp)
        print(hdr)
        print("    " + "-" * (len(hdr) - 4))
        for cp in cps:
            qs, ws = [], []
            for r in grp:
                m = [c for c in r.get("checkpoints", []) if c["cp_target"] == cp]
                qs.append(m[0]["q_foot"] if m else None)
                ws.append(m[0].get(f"w2_z{Z_PROBES[0]}") if m else None)
            got = [q for q in qs if q is not None]
            band = (100.0 * (max(got) - min(got)) / (0.5 * (max(got) + min(got)))
                    if len(got) >= 2 else float("nan"))
            print("    " + f"{cp:8.4f}"
                  + "".join(f"{q:12.2f}" if q is not None else f"{'--':>12}"
                            for q in qs)
                  + f"{band:9.2f}"
                  + "".join(f"{w:12.4f}" if w is not None and w == w
                            else f"{'--':>12}" for w in ws))

    print("\n--- THREE-MESH BAND, per density ---")
    for ename in sorted({r["e_name"] for r in legs}):
        grp = [r for r in legs if r["e_name"] == ename]
        adm = [r for r in grp if r["mode"] in _ADMISSIBLE and r["capacity"]]
        excl = [r for r in grp if r not in adm]
        seq = sorted(grp, key=lambda r: -r["h0"])
        pretty = " -> ".join(f"h{r['h0']}: {r['q_u']:.1f}" for r in seq)
        print(f"\n  density {ename}: {pretty}")
        if excl:
            print("    EXCLUDED as not a capacity: "
                  + "; ".join(f"{r['tag']} (mode {r['mode']}, plateau="
                              f"{r['plateau']}, peaked={r['peaked']}, "
                              f"free={r['free']})" for r in excl))
        if len(adm) < 2:
            print("    band UNAVAILABLE: fewer than two admissible legs")
            continue
        q = [r["q_u"] for r in adm]
        band = 100.0 * (max(q) - min(q)) / (0.5 * (max(q) + min(q)))
        print(f"    admissible legs: "
              + ", ".join(f"h{r['h0']}={r['q_u']:.2f}" for r in adm))
        print(f"    three-mesh band = {band:.2f} %  (max-min over the mean)")
        mono = all(seq[i + 1]["q_u"] < seq[i]["q_u"] for i in range(len(seq) - 1))
        print(f"    monotone decreasing under refinement: "
              f"{'YES (converging from above, R3 shape)' if mono else 'NO'}")
        if len(adm) == 3:
            a = sorted(adm, key=lambda r: -r["h0"])
            d1 = abs(a[1]["q_u"] - a[0]["q_u"])
            d2 = abs(a[2]["q_u"] - a[1]["q_u"])
            print(f"    successive differences {d1:.2f} then {d2:.2f} kPa "
                  f"(ratio {d2 / d1:.3f} -- < 1 is contracting)")

    print("\n--- COMPARISON BASIS ---")
    print(f"  R3 (DruckerPrager, perfectly plastic, WEIGHTLESS, SMOOTH footing) "
          f"reads {R3_SEQUENCE[1.0]} / {R3_SEQUENCE[0.5]} / {R3_SEQUENCE[0.25]} "
          f"of exact,")
    print(f"  i.e. a three-mesh spread of "
          f"{100 * (R3_SEQUENCE[1.0] - R3_SEQUENCE[0.25]) / (0.5 * (R3_SEQUENCE[1.0] + R3_SEQUENCE[0.25])):.1f} %, "
          f"accepted as converged because it is MONOTONE and every leg is a")
    print(f"  genuine plateau.  Its +/- {R3_BAND_HALFWIDTH_PCT:.0f} % is a "
          f"per-resolution tolerance on a KNOWN CENTRE, not a convergence "
          f"criterion.")
    print("  THE CAMPAIGN'S OWN TOLERANCE IS OQ2 AND IS NOT YET SUPPLIED "
          "(ADR 90 §7.5). No 'inside tolerance' verdict is available.")
    print("=" * 100)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
