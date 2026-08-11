"""Assemble the termination-mode table and the mechanism table from the
per-leg JSON that `quad_path_diag.py` writes.

The two tables are deliberately separate because they answer different
questions and only one of them is about a number:

  * the TERMINATION table says what each leg's q_max IS -- a capacity or a
    named controller allowance -- under the merged R3 gate's three clauses;
  * the MECHANISM table says whether a collapse mechanism was forming when the
    leg ended, which is the question note 81 section 5.5 asked for and the one
    that does not depend on any number being a ceiling.

Run:  py -3.12 quad_path_summary.py
"""
import glob
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ORDER = ["h8bbar", "h20uri", "h20std"]


def load():
    out = []
    for p in sorted(glob.glob(os.path.join(HERE, "qpd_*.json"))):
        with open(p) as f:
            r = json.load(f)
        if "_smoke" in r["tag"] or "_elastic" in r["tag"]:
            continue
        out.append(r)
    out.sort(key=lambda r: (r["h0"], ORDER.index(r["elem"])
                            if r["elem"] in ORDER else 9, r["tag"]))
    return out


def field_map(tag, xlim=7.0, zlim=-7.0):
    """An ASCII map of the mobilisation field over the core of the mesh.

    The npz carries the field for real plotting; this exists so the mechanism
    verdict is legible IN THE NOTE, next to the numbers that summarise it -- a
    localisation statistic that nobody can see the field behind is exactly the
    kind of claim this campaign has had to retract twice.

    Rows are depth (z = 0 at the top), columns are x.  The footing occupies
    |x| <= 1 m at the surface and is marked on the axis.
    """
    p = os.path.join(HERE, f"qpd_{tag}_field.npz")
    if not os.path.exists(p):
        return None
    d = np.load(p)
    cen, mob, x, z = d["cen"], d["mob_mean"], d["x"], d["z"]
    ii, kk = d["ii"], d["kk"]
    sel_i = np.where((0.5 * (x[:-1] + x[1:]) >= -xlim)
                     & (0.5 * (x[:-1] + x[1:]) <= xlim))[0]
    sel_k = np.where(0.5 * (z[:-1] + z[1:]) >= zlim)[0]
    grid = {}
    for e in range(len(cen)):
        grid[(int(ii[e]), int(kk[e]))] = mob[e]
    out = []
    for k in sorted(sel_k, reverse=True):                 # surface first
        line = ""
        for i in sel_i:
            m = grid.get((i, k), np.nan)
            if not np.isfinite(m):
                line += " "
            elif m >= 0.99:
                line += "#"
            elif m >= 0.90:
                line += "+"
            elif m >= 0.70:
                line += ":"
            elif m >= 0.40:
                line += "."
            else:
                line += " "
        out.append(f"  z={0.5 * (z[k] + z[k + 1]):+6.2f} |{line}|")
    xs = 0.5 * (x[sel_i] + x[sel_i + 1])
    foot = "".join("F" if abs(v) <= 1.0 else "-" for v in xs)
    out.append(f"          |{foot}|   F = footing, x from {xs[0]:.1f} to "
               f"{xs[-1]:.1f} m")
    out.append("          legend: '#' mob>=0.99  '+' >=0.90  ':' >=0.70  "
               "'.' >=0.40  ' ' elastic/far")
    return "\n".join(out)


def main():
    rows = load()
    if not rows:
        raise SystemExit("no qpd_*.json found")
    b = {r["build"] for r in rows}
    print(f"engine build(s): {', '.join(sorted(b))}")
    print(f"exact q_u = {rows[0]['q_exact']:.3f} kPa\n")

    h = (f"{'leg':>26} {'h0':>5} {'DOF':>6} {'q/q_ex':>7} {'tail%':>8} "
         f"{'MODE':>9} {'ds_end':>9} {'x floor':>8} {'subdiv':>9} {'halv':>5} "
         f"{'plat':>5} {'free':>5} {'CAP':>4} {'wall s':>7}")
    print("TERMINATION-MODE TABLE")
    print(h)
    print("-" * len(h))
    for r in rows:
        print(f"{r['tag']:>26} {r['h0']:5.2f} {r['ndof']:6d} {r['ratio']:7.4f} "
              f"{r['tail_pct']:8.3f} {r['mode']:>9} {r['ds_end_mm']:9.2e} "
              f"{r['headroom']:8.1f} "
              f"{str(r['nsub']) + '/' + str(r['budget']):>9} "
              f"{r['halvings']:5.2f} {'yes' if r['plateau'] else 'NO':>5} "
              f"{'yes' if r['free'] else 'NO':>5} "
              f"{'yes' if r['capacity'] else 'NO':>4} {r['wall_s']:7.0f}")
    print("\n  MODE: TARGET/BUDGET may be a capacity; FLOOR/WALL/STALL/TRUNCATED "
          "are SEIZURE and never are.")
    print("  CAP requires all three clauses: plateau AND admissible mode AND "
          "tail step >= 100x the floor.")
    print("  'halv' = log2(ds_base/ds_floor) = how many halvings the controller "
          "was ALLOWED -- the actual allowance.\n")

    h2 = (f"{'leg':>26} {'chi_el':>7} {'chi':>7} {'M':>7} {'M_tot':>7} "
          f"{'y90 vol%':>9} {'core%':>7} {'comp':>5} {'thick_el':>9} "
          f"{'fill%':>7} {'surf?':>6} {'|x|max':>7}")
    print("MECHANISM-FORMATION TABLE (at the last converged step)")
    print(h2)
    print("-" * len(h2))
    for r in rows:
        print(f"{r['tag']:>26} {r['chi_el']:7.3f} {r['chi']:7.3f} {r['M']:7.3f} "
              f"{r['M_total']:7.3f} {100 * r['vfrac_y90']:9.2f} "
              f"{100 * r['vfrac_core_y90']:7.1f} {r['ncomp']:5d} "
              f"{r['band_thick_el']:9.2f} {100 * r['band_fill']:7.1f} "
              f"{'yes' if r['touches_surface'] else 'NO':>6} "
              f"{r['xmax_any']:7.2f}")
    print("\n  chi = isochoric span ratio (1 = every bit of the footing "
          "increment came back up beside it).")
    print("  M   = (chi - chi_el)/(1 - chi_el): 0 = the increment is purely "
          "ELASTIC, 1 = a fully isochoric mechanism.")
    print("  thick_el / fill% describe the largest 8-connected yielded "
          "component: a slip SURFACE is ~1-3 elements")
    print("  thick and fills little of its bounding box; smeared plasticity is "
          "thick and fills most of it.\n")

    print("MOBILISATION FIELD AT THE LAST CONVERGED STEP (core of the mesh)")
    for r in rows:
        m = field_map(r["tag"])
        if m:
            print(f"\n  [{r['tag']}]  s/B = {r['s_end_over_B']:.5f}, "
                  f"q/q_exact = {r['ratio']:.4f}, mode {r['mode']}, "
                  f"M = {r['M']:+.3f}")
            print(m)
    print()

    cond = [r for r in rows if r.get("cond")]
    if cond:
        print("TANGENT CONDITIONING ON APPROACH TO TERMINATION")
        h3 = (f"{'leg':>26} {'s/B':>9} {'q':>8} {'ds mm':>9} "
              f"{'sig_min/sc':>11} {'cond':>10} {'lam_min_sym/sc':>15} "
              f"{'n_neg':>6}")
        print(h3)
        print("-" * len(h3))
        for r in cond:
            for c in r["cond"]:
                print(f"{r['tag']:>26} {c['s_over_B']:9.5f} {c['q']:8.2f} "
                      f"{c['ds_mm']:9.2e} {c['s_min_rel']:11.3e} "
                      f"{c['cond']:10.3e} {c['lam_min_sym_rel']:15.3e} "
                      f"{c['n_neg_sym']:6d}")
        # NOTE: ASCII only in every print below.  This script gets piped through
        # PowerShell, which decodes as cp1252 and turns a stray em-dash into a
        # UnicodeEncodeError and an exit code of 255 -- the banner capture trap
        # already in LEDGER_quirks.
        print("\n  sig_min/sc: smallest SINGULAR value over the mean diagonal - "
              "the conditioning measure for an")
        print("  UNSYMMETRIC tangent (rho_bar != rho).  lam_min_sym going "
              "NEGATIVE is loss of positive definiteness,")
        print("  which is what a genuine limit point / localisation produces; "
              "sig_min collapsing while lam_min_sym")
        print("  stays positive is ill-conditioning WITHOUT a mechanism.\n")

    # the allowance pairs: legs identical but for ds_floor / budget
    print("ALLOWANCE PAIRS (same element, same mesh, same loads; only the "
          "controller's allowance differs)")
    # A pair is two legs that differ ONLY in the controller's allowance.  The
    # settlement TARGET must match too: the matched-settlement mechanism controls
    # (`_at0103`, `_at0112`) share element, mesh and ladder with the full legs
    # but were deliberately stopped early, and pairing those would report the
    # truncation as an allowance effect (-86 % "reach", measured, on the first
    # cut of this block).  The allowance is `halvings` and `budget`; everything
    # else in the key must be equal.
    key = {}
    for r in rows:
        key.setdefault((r["elem"], r["h0"], r["ladder"], r["s_target"]),
                       []).append(r)
    any_pair = False
    for k, v in key.items():
        if len(v) < 2 or len({(r["halvings"], r["budget"]) for r in v}) < 2:
            continue
        any_pair = True
        v.sort(key=lambda r: r["halvings"])
        print(f"  {k[0]} h0={k[1]}, ladder {k[2]}, s/B target {k[3]}:")
        for r in v:
            print(f"     floor {r['ds_floor_mm']:.2e} mm "
                  f"({r['halvings']:.2f} halvings), budget {r['budget']:>4} "
                  f"-> q/q_ex {r['ratio']:.4f}, reach s/B {r['s_end_over_B']:.5f}, "
                  f"tail {r['tail_pct']:.2f} %, mode {r['mode']}, "
                  f"M {r['M']:+.3f}")
        a, z = v[0], v[-1]
        print(f"     DELTA: capacity {100 * (z['ratio'] / a['ratio'] - 1):+.1f} %, "
              f"reach {100 * (z['s_end_over_B'] / a['s_end_over_B'] - 1):+.1f} %, "
              f"tail {a['tail_pct']:.2f} -> {z['tail_pct']:.2f} %")
    if not any_pair:
        print("  (none yet — an allowance pair needs two legs differing only in "
              "floor/budget)")


if __name__ == "__main__":
    main()
