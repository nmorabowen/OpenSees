"""ADR-95 visualization harness: deformed-mesh + mean-stress snapshots for the
12-leg Prandtl-Reissner strip-footing survey (pre-fix vs fixed UW-DP, ASD-DP
cross-check, LadrunoSANISAND).

Reuses `h20_prandtl.py` (mesh / consistent surcharge / DS_LADDER / the D16
guards), `tet_path_diag.py` (tet10 mesh + Bernstein/Lagrange surcharge),
`asd_path_diag.py` (ASDPlasticMaterial3D model builder) and
`sanisand_path_diag.py` (LadrunoSANISAND model builder + staging) BY IMPORT --
none of those files are modified.  Each of the 12 legs is one process
invocation (`--leg <NAME>`) so the build (dist_fixed vs dist_p6, via
ADR95_DIST) is fixed for the life of the interpreter, matching how those
harnesses already select a binary.

Per leg this script dumps `snap_<LEG>.npz`:
    nodes (N,3), cells (E,nen), elem_type (str),
    disp_end / disp_001 (N,3)      -- nodal displacement at the end state and
                                       at the first converged step with
                                       s/B >= 0.01 (NaN-filled if the leg
                                       walled before reaching 0.01),
    p_end / p_001 (E,)             -- per-element mean stress (tension +ve)
    mob_end / mob_001 (E,)         -- per-element plastic-mobilisation
                                       fraction of Gauss points (UW: the
                                       ladrunoBranch[0] > 0 census; ASD/
                                       SANISAND: |dev|/(-p) vs the cone/
                                       critical ratio, threshold 0.9)
    mode / verdict / build / q_end / s_end_over_B / q_exact / ...

NO tangent/branch/cond-number sampling -- the push loop here is the fast
allowance-only ladder (the reference 3-rung ladder + the D16 ds guards),
nothing else.

Run:
    py -3.12 deformed_snapshot.py --leg A1
    py -3.12 deformed_snapshot.py --leg B2 --tmax 1200
"""
import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

_SCRATCH = ("C:/Users/nmb/AppData/Local/Temp/claude/"
            "C--Users-nmb-Documents-Github-OpenSees--claude-worktrees-"
            "prandtl-reissner-bezier-954a90/5408bcef-bf62-4967-964e-"
            "a179fa61676a/scratchpad")
DIST = {
    "prefix": os.path.join(_SCRATCH, "dist_fixed", "bin"),  # pre-fix UW-DP
    "fixed":  os.path.join(_SCRATCH, "dist_p6", "bin"),      # repaired build
}
Q_EXACT = 138.907
SNAP_S_OVER_B = 0.01

LEGS = {
    "A1": dict(dist="prefix", mat="uwdp",     mesh="hex", elem="h20uri"),
    "A2": dict(dist="prefix", mat="uwdp",     mesh="tet", elem="tet10"),
    "A3": dict(dist="prefix", mat="uwdp",     mesh="tet", elem="beziertet10"),
    "B1": dict(dist="fixed",  mat="uwdp",     mesh="hex", elem="h8bbar"),
    "B2": dict(dist="fixed",  mat="uwdp",     mesh="hex", elem="h20uri"),
    "B3": dict(dist="fixed",  mat="uwdp",     mesh="tet", elem="tet10"),
    "B4": dict(dist="fixed",  mat="uwdp",     mesh="tet", elem="beziertet10"),
    "B5": dict(dist="fixed",  mat="uwdp",     mesh="tet", elem="beziertet10bbar"),
    "C1": dict(dist="fixed",  mat="asddp",    mesh="hex", elem="h8bbar"),
    "C2": dict(dist="fixed",  mat="asddp",    mesh="hex", elem="h20uri"),
    "D1": dict(dist="fixed",  mat="sanisand", mesh="hex", elem="h20uri", implex=True),
    "D2": dict(dist="fixed",  mat="sanisand", mesh="hex", elem="h8bbar", implex=False),
}

# the y=0 side face's LOCAL node indices, H20_OFF ordering (h20_prandtl.py):
# corners [0,1,5,4] (H8 and the H20 corner block); H20 mid-edges 8 (e01),
# 17 (e15), 12 (e54), 16 (e40) -- Q8 order c0 c1 c2 c3 m01 m12 m23 m30.
Y0_LOC = {8: [0, 1, 5, 4], 20: [0, 1, 5, 4, 8, 17, 12, 16]}


def elem_stats(ops, n_cells, mat, params):
    """Cheap per-element census from the COMMITTED state only: mean stress p
    (tension-positive) and a plastic-mobilisation fraction of Gauss points.
    No tangent sampling, no extra analyze() calls."""
    p_arr = np.full(n_cells, np.nan)
    mob_arr = np.full(n_cells, np.nan)
    for e in range(1, n_cells + 1):
        try:
            st = np.asarray(ops.eleResponse(e, "stress"), dtype=float)
        except Exception:
            st = np.array([])
        if st.size < 6:
            continue
        st = st.reshape(-1, 6)
        p_gp = (st[:, 0] + st[:, 1] + st[:, 2]) / 3.0
        p_arr[e - 1] = float(np.nanmean(p_gp))
        ngp = st.shape[0]
        if mat == "uwdp":
            nb = nseen = 0
            for gp in range(1, ngp + 1):
                try:
                    br = ops.eleResponse(e, "material", gp, "ladrunoBranch")
                except Exception:
                    br = None
                if br is None or len(br) < 1:
                    continue
                nseen += 1
                if br[0] > 0.5:
                    nb += 1
            mob_arr[e - 1] = (nb / nseen) if nseen else float("nan")
        else:
            dev0 = st[:, 0] - p_gp
            dev1 = st[:, 1] - p_gp
            dev2 = st[:, 2] - p_gp
            j2 = (0.5 * (dev0 ** 2 + dev1 ** 2 + dev2 ** 2)
                  + st[:, 3] ** 2 + st[:, 4] ** 2 + st[:, 5] ** 2)
            q = np.sqrt(3.0 * j2)
            pc = -p_gp
            if mat == "asddp":
                denom = params["eta"] * pc + params["xic"]
                f_norm = np.divide(q, denom, out=np.full_like(q, np.nan),
                                    where=denom > 1e-9)
                mob_arr[e - 1] = float(np.nanmean(f_norm > 0.9))
            elif mat == "sanisand":
                eta_ratio = np.divide(q, pc, out=np.full_like(q, np.nan),
                                       where=pc > 1e-9)
                mob_arr[e - 1] = float(np.nanmean((eta_ratio / params["mc"]) > 0.9))
    return p_arr, mob_arr


def run_leg(name, cfg, tmax, budget, sfrac):
    os.environ["ADR95_DIST"] = DIST[cfg["dist"]]
    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

    import h20_prandtl as HP
    ops = HP.ops
    stamp = ops.ladrunoBuild()
    print(f"[snapshot] leg {name}: build={stamp} dist={cfg['dist']} "
          f"mat={cfg['mat']} elem={cfg['elem']} mesh={cfg['mesh']}", flush=True)

    mat = cfg["mat"]
    mesh_kind = cfg["mesh"]
    elem_key = cfg["elem"]
    mat_params = {}
    y0_loc = None

    if mesh_kind == "hex":
        order = 1 if elem_key == "h8bbar" else 2
        form = "bbar" if elem_key == "h8bbar" else "uri"
        nodes, cells, sets, xg, yg, zg = HP.strip_mesh(1.0, order)
        w, nfaces = HP.consistent_surcharge(nodes, cells, order)
        HP.verify_surcharge(nodes, cells, w, nfaces, order)
        n_cells = len(cells)
        elem_type = f"hex{cells.shape[1]}"
        y0_loc = Y0_LOC[cells.shape[1]]
        ladder_kind = "quad" if order == 2 else "linear"

        if mat == "uwdp":
            HP.build_model(nodes, cells, sets, w, form, assoc=False)
            tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
            HP.surcharge_step(nodes, sets, w, tol)
            ladder = HP.attempts(tol)
        elif mat == "asddp":
            import asd_path_diag as ASD
            _, _, eta_asd, xi_c_asd, _ = ASD.build_model_asd(nodes, cells, sets, w, form)
            mat_params = dict(eta=eta_asd, xic=xi_c_asd)
            tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
            HP.surcharge_step(nodes, sets, w, tol)
            ladder = ([("KrylovNewton", tol, 25, 0),
                       ("NewtonLineSearch", tol, 40, 0),
                       ("KrylovNewton", 10.0 * tol, 60, 1)]
                      if ladder_kind == "quad" else
                      [("Newton", tol, 25, 0),
                       ("NewtonLineSearch", tol, 40, 0),
                       ("KrylovNewton", 10.0 * tol, 60, 1)])
        elif mat == "sanisand":
            import sanisand_path_diag as SAN
            max_substeps = 1000 if cfg.get("implex", False) else 0
            SAN.build_model_sanisand(nodes, cells, sets, w, form, max_substeps, cfg.get("implex", False))
            mat_params = dict(mc=SAN.M_C)
            ops.updateMaterialStage("-material", 1, "-stage", 0)
            tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
            HP.surcharge_step(nodes, sets, w, tol)
            ops.updateMaterialStage("-material", 1, "-stage", 1)
            ladder = ([("KrylovNewton", tol, 25, 0),
                       ("NewtonLineSearch", tol, 40, 0),
                       ("KrylovNewton", 10.0 * tol, 60, 1)]
                      if ladder_kind == "quad" else
                      [("Newton", tol, 25, 0),
                       ("NewtonLineSearch", tol, 40, 0),
                       ("KrylovNewton", 10.0 * tol, 60, 1)])
        else:
            raise SystemExit(f"unknown mat {mat!r}")

        base, floor, dmax = HP.DS_LADDER[order]

    else:  # mesh_kind == "tet"
        import tet_path_diag as TPD
        nodes, tets, sets, meta = TPD.load_mesh()
        n_cells = len(tets)
        cells = tets
        basis = TPD.ELEMS[elem_key]["basis"]
        w, nfaces = TPD.consistent_surcharge_tet(nodes, tets, basis)
        area_top = 2.0 * HP.XLIM * HP.THICK
        TPD.verify_surcharge_tet(nodes, w, nfaces, basis, area_top)
        elem_type = f"tet10_{elem_key}"

        assert mat == "uwdp", "the tet legs in this survey are UW-DP only"
        TPD.build_model(nodes, tets, sets, w, elem_key, assoc=False)
        tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
        ops.constraints("Transformation")
        ops.numberer("RCM")
        try:
            ops.system("Pardiso")
        except Exception:
            ops.system("UmfPack")
        ops.test("NormUnbalance", tol, 25, 0)
        ops.algorithm("Newton")
        ops.integrator("LoadControl", 1.0)
        ops.analysis("Static")
        if ops.analyze(1) != 0:
            raise SystemExit("surcharge step failed")

        base, floor, dmax = HP.DS_LADDER[2]
        ladder = [("KrylovNewton", tol, 25, 0),
                  ("NewtonLineSearch", tol, 40, 0),
                  ("KrylovNewton", 10.0 * tol, 60, 1)]

    # ---- the push, common to every leg --------------------------------------
    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = HP.Q0 * float(w[sets["footing"]].sum())
    area = HP.B_FOOT * HP.THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = sfrac * HP.B_FOOT
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso")
    except Exception:
        ops.system("UmfPack")
    ops.analysis("Static")

    def snapshot():
        n = len(nodes)
        disp = np.zeros((n, 3))
        for i in range(1, n + 1):
            d = ops.nodeDisp(i)
            disp[i - 1, :min(3, len(d))] = d[:3]
        p_arr, mob_arr = elem_stats(ops, n_cells, mat, mat_params)
        return disp, p_arr, mob_arr

    rows = []
    ds, good, nfail, nrelax, nsub = base, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"
    snap001 = None
    t0 = time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > tmax:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now / HP.B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", -ds)
        ok, relaxed = False, 0
        for algo, tl, it, rl in ladder:
            ops.test("NormUnbalance", tl, it, 0)
            ops.algorithm(algo)
            rc = ops.analyze(1)
            if rc == 0:
                ok, relaxed = True, rl
                break
            nfail += 1
        if not ok:
            good, nsub = 0, nsub + 1
            ds *= 0.5
            if nsub > budget:
                mode = "BUDGET"
                verdict = f"subdivision budget of {budget} spent at s/B = {s_now / HP.B_FOOT:.5f}"
                break
            if ds < floor:
                mode = "FLOOR"
                verdict = f"step collapsed to the {floor * 1e3:.4g} mm floor at s/B = {s_now / HP.B_FOOT:.5f}"
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= HP.GROW_AFTER and ds < dmax:
            ds, good = min(2 * ds, dmax), 0
        ops.reactions()
        q = (-sum(ops.nodeReaction(t, 3) for t in foot) + r_corr) / area
        s = uz0 - ops.getTime()
        s_ob = s / HP.B_FOOT
        rows.append((s, s_ob, q, ds * 1e3, relaxed, time.time() - t0))
        if snap001 is None and s_ob >= SNAP_S_OVER_B:
            d001, p001, m001 = snapshot()
            snap001 = dict(disp=d001, p=p001, mob=m001, s_over_B=s_ob, q=q)
            print(f"    [snap] leg {name}: s/B={SNAP_S_OVER_B} crossing captured "
                  f"at s/B={s_ob:.5f} q={q:.3f}", flush=True)
        if len(rows) % 20 == 0:
            print(f"    [{name}] step {len(rows)} s/B={s_ob:.5f} q={q:.3f} "
                  f"ds={ds * 1e3:.4g}mm wall={time.time() - t0:.0f}s", flush=True)
        if len(rows) > HP.STALL_WINDOW and \
                rows[-1][0] - rows[-1 - HP.STALL_WINDOW][0] < HP.STALL_ADVANCE * smax:
            mode = "STALL"
            verdict = f"stalled at s/B = {s_ob:.5f}"
            break

    wall = time.time() - t0
    dend, pend, mend = snapshot()
    q_end = rows[-1][2] if rows else float("nan")
    s_end = rows[-1][1] if rows else 0.0
    print(f"[snapshot] leg {name} MODE={mode}  [{verdict}]", flush=True)
    print(f"[snapshot] leg {name}: s_end/B={s_end:.5f} q_end={q_end:.3f} "
          f"wall={wall:.0f}s steps={len(rows)} nsub={nsub} nfail={nfail}", flush=True)

    q_exact = Q_EXACT if mat in ("uwdp", "asddp") else float("nan")

    out = dict(
        leg=name, elem=elem_key, mesh=mesh_kind, material=mat, dist=cfg["dist"],
        build=str(stamp), mode=mode, verdict=verdict,
        nodes=nodes, cells=cells, elem_type=elem_type,
        disp_end=dend, p_end=pend, mob_end=mend,
        disp_001=(snap001["disp"] if snap001 else np.full_like(dend, np.nan)),
        p_001=(snap001["p"] if snap001 else np.full_like(pend, np.nan)),
        mob_001=(snap001["mob"] if snap001 else np.full_like(mend, np.nan)),
        s_end_over_B=s_end, q_end=q_end,
        s_001_over_B=(snap001["s_over_B"] if snap001 else np.nan),
        q_001=(snap001["q"] if snap001 else np.nan),
        q_exact=q_exact, wall_s=wall, nsub=nsub, nfail=nfail, nrelax=nrelax,
        steps=len(rows), budget=budget, sfrac=sfrac,
    )
    if y0_loc is not None:
        out["y0_loc"] = np.array(y0_loc)
    np.savez(os.path.join(HERE, f"snap_{name}.npz"), **out)
    print(f"[snapshot] wrote snap_{name}.npz", flush=True)
    return out


# ---------------------------------------------------------------------------
# plotting -- reads the snap_<LEG>.npz files this same script wrote and builds
# the two survey figures.  No solver imports below this point.
# ---------------------------------------------------------------------------
_TET_EDGES = [(0, 1), (1, 2), (0, 2), (0, 3), (1, 3), (2, 3)]
_DIST_LABEL = {"prefix": "PRE-FIX", "fixed": "FIXED"}
_P_CLIP = 5.0


def _panel(ax, d, scale, xlim, zlim):
    import matplotlib.colors as mcolors
    from matplotlib.collections import LineCollection, PolyCollection

    nodes = d["nodes"]
    cells = d["cells"]
    elem_type = str(d["elem_type"])
    disp = np.nan_to_num(d["disp_end"], nan=0.0)
    p = np.nan_to_num(d["p_end"], nan=0.0)
    mode = str(d["mode"])
    leg = str(d["leg"])
    elem = str(d["elem"])
    material = str(d["material"])
    dist = str(d["dist"])
    s_end = float(d["s_end_over_B"])
    q_end = float(d["q_end"])
    q_exact = float(d["q_exact"])

    xz = nodes[:, [0, 2]] + scale * disp[:, [0, 2]]
    cmap = plt.get_cmap("RdBu_r")
    norm = mcolors.Normalize(vmin=-_P_CLIP, vmax=_P_CLIP)
    pc_clip = np.clip(p, -_P_CLIP, _P_CLIP)

    if elem_type.startswith("hex"):
        y0 = d["y0_loc"][:4]
        quads = xz[cells[:, y0]]
        colors = cmap(norm(pc_clip))
        coll = PolyCollection(quads, facecolors=colors, edgecolors="black",
                               linewidths=0.15)
        ax.add_collection(coll)
    else:
        segs, seg_colors = [], []
        colors = cmap(norm(pc_clip))
        for e in range(cells.shape[0]):
            conn = cells[e]
            col = colors[e]
            for a, b in _TET_EDGES:
                segs.append([xz[conn[a]], xz[conn[b]]])
                seg_colors.append(col)
        ax.add_collection(LineCollection(segs, colors=seg_colors, linewidths=0.4))

    ax.plot([-1.0, 1.0], [0.0, 0.0], color="black", linewidth=3.0,
             solid_capstyle="butt", zorder=5)
    ax.set_xlim(*xlim)
    ax.set_ylim(*zlim)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=6)

    if np.isfinite(q_exact):
        qtxt = f"q={q_end:.1f} (q/qex={q_end / q_exact:.3f})"
    else:
        qtxt = f"q={q_end:.1f} (no q_exact)"
    ax.set_title(
        f"{leg}: {elem}/{material.upper()}/{_DIST_LABEL[dist]}\n"
        f"{mode}  s/B={s_end:.4f}  {qtxt}  x{scale:.0f}",
        fontsize=6.8, linespacing=1.3)
    return norm, cmap


def build_figures(which):
    global plt
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F811 (module-level for _panel too)

    legs = list(LEGS)
    data = {}
    missing = []
    for name in legs:
        path = os.path.join(HERE, f"snap_{name}.npz")
        if os.path.exists(path):
            data[name] = np.load(path, allow_pickle=True)
        else:
            missing.append(name)
    if missing:
        print(f"[plot] WARNING: missing snap_*.npz for legs {missing}; "
              f"those panels will be blank", flush=True)

    def make(fname, xlim, zlim, settle_target):
        fig, axes = plt.subplots(3, 4, figsize=(18, 12), dpi=130)
        fig.subplots_adjust(left=0.035, right=0.90, top=0.90, bottom=0.04,
                             hspace=0.55, wspace=0.25)
        last = None
        for ax, name in zip(axes.flat, legs):
            if name not in data:
                ax.set_visible(False)
                continue
            d = data[name]
            s_m = float(d["s_end_over_B"]) * 2.0   # B_FOOT = 2.0 m
            scale = settle_target / s_m if s_m > 1e-9 else 1.0
            last = _panel(ax, d, scale, xlim, zlim)
        for ax in axes.flat[len(legs):]:
            ax.set_visible(False)
        if last is not None:
            norm, cmap = last
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
            fig.colorbar(sm, cax=cax,
                         label="mean stress p (kPa, tension +), clipped to "
                               f"+/-{_P_CLIP:.0f}")
        fig.suptitle("ADR-95 Prandtl-Reissner strip footing: deformed mesh, "
                      f"{'|'.join(legs)}", fontsize=10, y=0.965)
        out = os.path.join(HERE, fname)
        fig.savefig(out)
        plt.close(fig)
        print(f"[plot] wrote {out}", flush=True)

    if which in ("main", "both"):
        make("adr95_deformed.png", (-8.0, 8.0), (-8.0, 0.8), 0.5)
    if which in ("full", "both"):
        make("adr95_deformed_full.png", (-30.0, 30.0), (-20.0, 0.8), 2.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leg", choices=sorted(LEGS))
    ap.add_argument("--tmax", type=float, default=1200.0)
    ap.add_argument("--budget", type=int, default=200)
    ap.add_argument("--sfrac", type=float, default=0.15)
    ap.add_argument("--plot", choices=["main", "full", "both"], default=None,
                     help="skip the leg run; build the deformed-mesh figures "
                          "from the snap_<LEG>.npz files already on disk")
    args = ap.parse_args()
    if args.plot:
        build_figures(args.plot)
        return
    if not args.leg:
        raise SystemExit("--leg is required unless --plot is given")
    run_leg(args.leg, LEGS[args.leg], args.tmax, args.budget, args.sfrac)


if __name__ == "__main__":
    main()
