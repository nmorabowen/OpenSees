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

import glob
import json
import os
import sys

R3_BAND_HALFWIDTH_PCT = 3.0
R3_SEQUENCE = {1.0: 1.0849, 0.5: 0.9977, 0.25: 0.9514}   # measured, DruckerPrager
_ADMISSIBLE = ("TARGET", "PEAK", "BUDGET")


def load(out_dir):
    legs = []
    for p in sorted(glob.glob(os.path.join(out_dir, "a2_*.json"))):
        with open(p) as f:
            legs.append(json.load(f))
    return legs


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print(__doc__)
        return 2
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

    print("\n--- WIDTH vs h (ADR 90 §7.3 threshold-free w2, "
          "fixed physical depth, x >= 0 half) ---")
    hdr = (f"{'leg':>16} {'h0':>6} {'w2(z=-0.5)':>11} {'w2/h0':>8} "
           f"{'w2(z=-1.0)':>11} {'w2/h0':>8} {'yield ele':>10} "
           f"{'yield vol':>10} {'epsq max':>10}")
    print(hdr)
    print("-" * len(hdr))
    for r in legs:
        w05, w10 = r.get("w2_z-0.5"), r.get("w2_z-1.0")
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
                ws.append(m[0].get("w2_z-0.5") if m else None)
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
