"""Read the dp_collapse legs and answer the question the campaign left open:
what IS the collapse load of this footing on the cone PDMY actually has, and
how far off was the Chen-Han plane-strain route that estimated it at 52x the
quoted anchor?

Outputs a table to stdout and `dp_collapse.png`.

Run:  py -3.12 dp_analyze.py
"""
import os

import numpy as np

# analysis only — never load the engine (see dp_collapse.py). This lets the
# figure be produced by the one interpreter here that HAS matplotlib.
os.environ["ADR79_NO_ENGINE"] = "1"
import dp_collapse as dc  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PDMY_Q_AT_15 = 3384.0            # kPa, the campaign's -geom hypo backbone end


def tangent(s, q, lo, hi):
    m = (s >= lo) & (s <= hi)
    return np.polyfit(s[m], q[m], 1)[0] if m.sum() > 3 else np.nan


def read(name):
    p = os.path.join(HERE, f"dpcollapse_{name}.csv")
    if not os.path.exists(p):
        return None
    a = np.genfromtxt(p, delimiter=",", names=True)
    if a.size < 20:
        return None
    return dict(s=np.atleast_1d(a["s_m"]), q=np.atleast_1d(a["q_kPa"]))


def summarize(name):
    d = read(name)
    if d is None:
        return None
    s, q = d["s"], d["q"]
    cfg = dc.LEGS[name]
    mc = dc.mc_from_cone(cfg["alpha"])
    qmax = float(q.max())
    # the tangent over the last 10 % of the settlement range, and over the
    # decade before it, so "plateau" is a comparison and not an assertion
    s_end = s[-1]
    t_last = tangent(s, q, 0.90 * s_end, s_end)
    t_prev = tangent(s, q, 0.70 * s_end, 0.80 * s_end)
    n0 = max(5, len(s) // 50)          # a point count, not a settlement window:
    t_init = np.polyfit(s[:n0], q[:n0], 1)[0]   # the early steps are the finest
    return dict(name=name, cfg=cfg, mc=mc, qmax=qmax, s_end=s_end,
                s_at_max=float(s[int(q.argmax())]), n=len(s),
                t_last=t_last, t_prev=t_prev, t_init=t_init,
                vesic_ps=dc.vesic(mc["ps"])["q_u"],
                vesic_txc=dc.vesic(mc["txc"])["q_u"], s=s, q=q)


def field(name):
    p = os.path.join(HERE, f"dpcollapse_{name}_field.npz")
    if not os.path.exists(p):
        return None
    z = np.load(p)
    return z["centroid"], z["mob"]


def collapse_estimate(s, q, s_ref=0.10 * dc.B_FOOT):
    """A number to quote, and an honest account of what it is.

    None of these legs ends in a sharp limit point. The tangent decays
    monotonically — 41 432 -> 605 kPa/m across s/B = 0.5 % -> 10 % on the base
    leg — but as a POWER LAW in settlement rather than to zero at a finite s,
    which is the signature of a punching mechanism rather than general shear.
    So two things are reported:

      q_ref  : q at s/B = 10 %, the settlement-based ultimate-capacity
               criterion conventional for punching (de Beer / Vesic).
      q_inf  : the tail extrapolation. Fit dq/ds = C s^-p over the last decade;
               if p > 1 that integrates to a finite limit
               q_inf = q_ref + C s_ref^(1-p)/(p-1). If p <= 1 there is no
               finite limit and the extrapolation is refused.

    q_inf is an EXTRAPOLATION off a two-decade fit and is quoted as a bracket,
    never as the measurement.
    """
    out = dict(q_ref=float(np.interp(s_ref, s, q)), s_ref=s_ref,
               p=np.nan, q_inf=np.nan)
    if s[-1] < s_ref:
        out["q_ref"] = np.nan
        return out
    # tangents in log-spaced windows over the last decade of settlement
    lo = max(s[0], 0.1 * s[-1])
    edges = np.geomspace(lo, s[-1], 6)
    sm, km = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (s >= a) & (s <= b)
        if m.sum() > 3:
            sm.append(np.sqrt(a * b))
            km.append(np.polyfit(s[m], q[m], 1)[0])
    sm, km = np.array(sm), np.array(km)
    if len(sm) < 3 or (km <= 0).any():
        return out
    p, lnC = np.polyfit(np.log(sm), np.log(km), 1)
    p = -p
    out["p"] = float(p)
    if p > 1.0:
        C = np.exp(lnC)
        out["q_inf"] = float(out["q_ref"]
                             + C * s_ref ** (1.0 - p) / (p - 1.0))
    return out


def anchor_ladder():
    """Every candidate anchor for this footing, in one place.

    The project quotes Vesic at the NOMINAL phi = 33 deg. cone_probe.py argued
    for the plane-strain equivalent of the measured cone, which is 42x higher.
    Neither is right: the plane-strain matching assumes ASSOCIATED flow, and
    this PDMY set is non-dilatant (dil1 = dil2 = 0). Applying Davis's reduction
    for psi = 0 AFTER the plane-strain matching is the anchor that actually
    corresponds to the material, and it is the one the numerics land on.
    """
    import dp_strip as dst
    mc = dc.mc_from_cone(dc.ALPHA_PDMY)
    rows = [("Vesic at the nominal phi = 33 (what is quoted)", 33.0),
            ("cone, triaxial-compression equivalent", mc["txc"]),
            ("cone, plane-strain equivalent (Chen-Han)", mc["ps"]),
            ("cone, plane-strain THEN Davis psi = 0", dst.davis_phi(mc["ps"]))]
    print(f"\nanchor ladder for the 2 x 2 m square footing "
          f"(q0 = {dc.Q_SURCH} kPa, gamma' = {dc.GAMMA_B:.2f} kN/m3)")
    hdr = (f"{'basis':>46} {'phi':>7} {'N_q':>9} {'N_gamma':>10} "
           f"{'q-term':>9} {'gamma-term':>11} {'q_u':>10}")
    print(hdr)
    print("-" * len(hdr))
    for lbl, phi in rows:
        v = dc.vesic(phi)
        print(f"{lbl:>46} {phi:7.2f} {v['Nq']:9.1f} {v['Ngamma']:10.1f} "
              f"{v['q_term']:9.1f} {v['g_term']:11.1f} {v['q_u']:10.1f}")
    return {lbl: dc.vesic(phi)["q_u"] for lbl, phi in rows}


def matching_sensitivity():
    """How stable is the Chen-Han plane-strain route at THIS cone slope?

    The plane-strain matching is alpha = tan(phi)/sqrt(9 + 12 tan^2 phi), which
    inverts to tan(phi) = sqrt(9 alpha^2 / (1 - 12 alpha^2)). That denominator
    vanishes at alpha = 1/sqrt(12) = 0.28868: ABOVE that slope a Drucker-Prager
    cone has NO plane-strain Mohr-Coulomb equivalent at all. The measured cone
    sits at 0.2436, which is 84 % of the way there, and the map is nearly
    vertical — so a small uncertainty in alpha is a large uncertainty in N_q.
    cone_probe.py's own txc-vs-txe spread is 3.9 %, which is not small.
    """
    a0 = dc.ALPHA_PDMY
    print(f"\nChen-Han plane-strain matching, sensitivity at alpha = {a0:.5f}")
    print(f"  the matching has a POLE at alpha = 1/sqrt(12) = "
          f"{1 / np.sqrt(12):.5f}; this cone is at "
          f"{100 * a0 / (1 / np.sqrt(12)):.1f} % of it")
    print(f"{'alpha':>9} {'d(alpha)':>9} {'phi_ps':>8} {'N_q':>10} "
          f"{'N_gamma':>10} {'q_u square':>11}")
    for f in (0.94, 0.98, 1.00, 1.02, 1.06, 1.10):
        a = a0 * f
        mc = dc.mc_from_cone(a)
        if not np.isfinite(mc["ps"]):
            print(f"{a:9.5f} {100 * (f - 1):+8.1f} %  "
                  f"NO plane-strain equivalent exists (alpha > 1/sqrt(12))")
            continue
        v = dc.vesic(mc["ps"])
        print(f"{a:9.5f} {100 * (f - 1):+8.1f} % {mc['ps']:8.2f} "
              f"{v['Nq']:10.1f} {v['Ngamma']:10.1f} {v['q_u']:11.0f}")


def strip_table():
    """The plane-strain legs, including the mesh-refinement series.

    Two things are being read here at once. `nq20_*` has an EXACT oracle
    (Prandtl-Reissner q0*N_q for associated Mohr-Coulomb), so its column says
    how much of any gap is discretisation. The cone legs then inherit that
    reading. A leg that never plateaus has not measured a collapse load at all,
    and the `plateau` column says so rather than letting q_max pass as one.
    """
    import dp_strip as dst
    rows = []
    for name, cfg in dst.SLEGS.items():
        p = os.path.join(HERE, f"dpstrip_{name}.csv")
        if not os.path.exists(p):
            continue
        a = np.genfromtxt(p, delimiter=",", names=True)
        if a.size < 20:
            continue
        s, q = np.atleast_1d(a["s_m"]), np.atleast_1d(a["q_kPa"])
        mc = dc.mc_from_cone(cfg["alpha"])
        vs = dst.vesic_strip(mc["ps"], cfg["q0"], cfg["gamma"])
        m = s >= 0.9 * s[-1]
        n0 = max(5, len(s) // 50)
        t_last = np.polyfit(s[m], q[m], 1)[0]
        t_init = np.polyfit(s[:n0], q[:n0], 1)[0]
        rows.append((name, cfg, vs, float(q.max()), float(s[-1]),
                     t_last, t_init, dst.k0_mobilisation(cfg["alpha"],
                                                         cfg["nu"])))
    if not rows:
        return
    hdr = (f"\n{'strip leg':>22} {'h0':>6} {'el/B':>5} {'flow':>6} {'m0':>6} "
           f"{'s_end/B':>8} {'q_num':>9} {'oracle':>9} {'num/orc':>8} "
           f"{'dq/ds':>9} {'/init':>7} {'plateau':>8}")
    print(hdr)
    print("-" * (len(hdr) - 1))
    for name, cfg, vs, qmax, s_end, t_last, t_init, m0 in rows:
        plateau = abs(t_last) < 0.02 * abs(t_init)
        print(f"{name:>22} {cfg['h0']:6.3f} "
              f"{int(round(2.0 / cfg['h0'])):5d} "
              f"{'assoc' if cfg['assoc'] else 'non':>6} {m0:6.3f} "
              f"{s_end / dc.B_FOOT:8.4f} {qmax:9.1f} {vs['q_u']:9.1f} "
              f"{qmax / vs['q_u']:8.4f} {t_last:9.0f} "
              f"{t_last / t_init:7.3f} {'yes' if plateau else 'NO':>8}")
    print("  m0 = how close the 1-D initial state already is to yield; a leg "
          "with m0 > 0.8 is void.")
    print("  'NO' in plateau means the leg walled out on convergence, so q_num "
          "is a LOWER bound on that mesh, not a collapse load.")


def main():
    anch = anchor_ladder()
    matching_sensitivity()
    try:
        strip_table()
    except Exception as e:                       # strips are optional
        print(f"[strip table skipped: {e}]")
    names = [n for n in dc.LEGS
             if os.path.exists(os.path.join(HERE, f"dpcollapse_{n}.csv"))]
    res = [r for r in (summarize(n) for n in names) if r]
    if not res:
        print("\nno square-footing legs to analyze yet")
        return

    v33 = dc.vesic(33.0)
    print(f"Vesic at the NOMINAL phi = 33 deg (what the project quotes): "
          f"q_u = {v33['q_u']:.1f} kPa  "
          f"(q-term {v33['q_term']:.1f} + gamma-term {v33['g_term']:.1f})")
    print(f"PDMY benchmark backbone, -geom hypo, s/B = 0.15: "
          f"{PDMY_Q_AT_15:.0f} kPa\n")

    hdr = (f"{'leg':>16} {'flow':>6} {'form':>5} {'phi_ps':>7} {'s_end/B':>8} "
           f"{'q_max':>9} {'dq/ds_end':>10} {'/init':>7} {'q/Vesic_ps':>11} "
           f"{'q/PDMY':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(res, key=lambda r: r["name"]):
        print(f"{r['name']:>16} "
              f"{'assoc' if r['cfg']['assoc'] else 'non' :>6} "
              f"{r['cfg']['form']:>5} {r['mc']['ps']:7.2f} "
              f"{r['s_end'] / dc.B_FOOT:8.4f} {r['qmax']:9.1f} "
              f"{r['t_last']:10.0f} "
              f"{r['t_last'] / r['t_init'] if r['t_init'] else np.nan:7.3f} "
              f"{r['qmax'] / r['vesic_ps']:11.3f} "
              f"{r['qmax'] / PDMY_Q_AT_15:8.3f}")

    print("\ncollapse estimate — q at the s/B = 10 % ultimate-capacity "
          "criterion, and the tail extrapolation:")
    for r in sorted(res, key=lambda r: r["name"]):
        ce = collapse_estimate(r["s"], r["q"])
        r["ce"] = ce
        inf = ("no finite limit (p <= 1)" if not np.isfinite(ce["q_inf"])
               else f"{ce['q_inf']:8.0f} kPa")
        ref = "--" if not np.isfinite(ce["q_ref"]) else f"{ce['q_ref']:8.1f}"
        print(f"  {r['name']:>16}: q(s/B=10 %) = {ref} kPa, "
              f"dq/ds ~ s^-{ce['p']:.2f} -> q_inf = {inf}")

    print("\nplateau evidence (tangent over the last 10 % of settlement vs the "
          "70-80 % window):")
    for r in sorted(res, key=lambda r: r["name"]):
        ratio = r["t_last"] / r["t_prev"] if r["t_prev"] else np.nan
        print(f"  {r['name']:>16}: dq/ds = {r['t_prev']:8.0f} -> "
              f"{r['t_last']:8.0f} kPa/m  ({ratio:.3f} of the earlier window), "
              f"{'COLLAPSE' if abs(r['t_last']) < 0.02 * abs(r['t_init']) else 'still hardening'}")

    print("\nyielded-zone extent at the end of each leg (m > 0.99). 'boundary'"
          " is tested by ELEMENT COLUMN, not by a fraction of the domain — the"
          " mesh is graded and its outermost hex is 3.1 m wide:")
    for r in sorted(res, key=lambda r: r["name"]):
        f = field(r["name"])
        if f is None:
            continue
        cen, mob = f
        y = mob > 0.99
        if not y.any():
            print(f"  {r['name']:>16}: nothing at m > 0.99")
            continue
        c = cen[y]
        cols = np.unique(np.round(np.abs(cen[:, 0]), 4))
        rows = np.unique(np.round(cen[:, 2], 4))
        outer = np.isclose(np.abs(np.round(cen[:, 0], 4)), cols.max())
        base = np.isclose(np.round(cen[:, 2], 4), rows.min())
        print(f"  {r['name']:>16}: {y.sum():4d}/{len(mob)} elements "
              f"({100 * y.mean():4.1f} %), |x| <= {np.abs(c[:, 0]).max():5.2f} m, "
              f"z >= {c[:, 2].min():6.2f} m  |  outermost column "
              f"{(y & outer).sum():3d}/{outer.sum():3d} yielded, base row "
              f"{(y & base).sum():3d}/{base.sum():3d}"
              + ("  <-- REACHES THE SIDE" if (y & outer).any() else ""))

    # --- the regularizer's price -----------------------------------------
    # sigma_y is a numerical apex offset, not a soil property, so the collapse
    # load has to be extrapolated back to sigma_y = 0. A cohesion enters a
    # bearing capacity LINEARLY (the c*N_c term), so two points do it — and
    # the slope measured here is worth comparing against what Vesic's N_c at
    # phi_ps = 53.7 would predict, which is a further test of that formula.
    a = next((r for r in res if r["name"] == "nonassoc"), None)
    b = next((r for r in res if r["name"] == "nonassoc_sy2"), None)
    if a and b:
        sa, sb = dc.LEGS["nonassoc"]["sy"], dc.LEGS["nonassoc_sy2"]["sy"]
        slope = (b["qmax"] - a["qmax"]) / (sb - sa)
        q0 = a["qmax"] - slope * sa
        mc = a["mc"]
        ph = np.radians(mc["ps"])
        nq = np.exp(np.pi * np.tan(ph)) * np.tan(np.pi / 4 + ph / 2) ** 2
        nc = (nq - 1.0) / np.tan(ph)
        sc = 1.0 + nq / nc
        # sigma_y -> Mohr-Coulomb c matched in triaxial compression
        pt = np.radians(mc["txc"])
        c_per_sy = (3.0 - np.sin(pt)) / (6.0 * np.cos(pt))
        print(f"\nnumerical-cohesion sweep: q_max = {a['qmax']:.1f} kPa at "
              f"sigma_y = {sa}, {b['qmax']:.1f} at {sb}")
        print(f"  slope {slope:.1f} kPa per kPa of sigma_y -> extrapolated to "
              f"sigma_y = 0: q_u = {q0:.1f} kPa "
              f"({100 * (a['qmax'] / q0 - 1):+.2f} % on the base leg)")
        print(f"  Vesic would predict c*N_c*s_c = "
              f"{c_per_sy * nc * sc:.0f} kPa per kPa of sigma_y "
              f"(N_c = {nc:.0f}, s_c = {sc:.2f}) — "
              f"{c_per_sy * nc * sc / slope:.0f}x the measured slope")

    # --- the answer -------------------------------------------------------
    head = [r for r in res if r["name"] in ("nonassoc", "assoc")]
    if len(head) == 2:
        na = next(r for r in head if not r["cfg"]["assoc"])
        aa = next(r for r in head if r["cfg"]["assoc"])
        print(f"\nnon-associated / associated = {na['qmax'] / aa['qmax']:.3f} "
              f"— the size of the flow-rule assumption Chen & Han's matching "
              f"makes")
        print(f"Chen-Han plane-strain Vesic said {na['vesic_ps']:.0f} kPa; "
              f"the 3D square footing measures {na['qmax']:.0f} kPa "
              f"(non-associated) / {aa['qmax']:.0f} kPa (associated)")
        print(f"benchmark PDMY backbone / measured collapse load = "
              f"{PDMY_Q_AT_15 / na['qmax']:.3f} (non-assoc), "
              f"{PDMY_Q_AT_15 / aa['qmax']:.3f} (assoc)")
    na = next((r for r in res if r["name"] == "nonassoc"), None)
    if na and "ce" in na:
        dv = anch["cone, plane-strain THEN Davis psi = 0"]
        print(f"\nagainst the Davis-reduced anchor ({dv:.0f} kPa):")
        for lbl, v in (("q at s/B = 10 %", na["ce"]["q_ref"]),
                       ("tail extrapolation q_inf", na["ce"]["q_inf"]),
                       ("PDMY backbone at s/B = 15 %", PDMY_Q_AT_15)):
            if np.isfinite(v):
                print(f"  {lbl:>28} = {v:8.1f} kPa -> {v / dv:.3f} x")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    for r in sorted(res, key=lambda r: r["name"]):
        ax[0].plot(r["s"] / dc.B_FOOT * 100, r["q"], label=r["name"], lw=1.4)
    ax[0].axhline(PDMY_Q_AT_15, color="k", ls="--", lw=1,
                  label="PDMY backbone @ s/B=15 %")
    ax[0].axhline(v33["q_u"], color="0.5", ls=":", lw=1,
                  label="Vesic, nominal phi=33")
    ax[0].set_xlabel("s / B  [%]")
    ax[0].set_ylabel("q  [kPa]")
    ax[0].set_title("Drucker-Prager collapse, the cone PDMY has")
    ax[0].legend(fontsize=7)
    ax[0].grid(alpha=0.3)
    for r in sorted(res, key=lambda r: r["name"]):
        s, q = r["s"], r["q"]
        k = max(len(s) // 60, 3)
        sm = s[::k]
        qm = q[::k]
        ax[1].plot(sm[1:] / dc.B_FOOT * 100, np.diff(qm) / np.diff(sm),
                   label=r["name"], lw=1.2)
    ax[1].set_xlabel("s / B  [%]")
    ax[1].set_ylabel("dq/ds  [kPa/m]")
    ax[1].set_yscale("symlog", linthresh=100)
    ax[1].axhline(0, color="k", lw=0.8)
    ax[1].set_title("tangent — a collapse load needs dq/ds -> 0")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(HERE, "dp_collapse.png")
    fig.savefig(out, dpi=130)
    print(f"\n[fig] {out}")


if __name__ == "__main__":
    main()
