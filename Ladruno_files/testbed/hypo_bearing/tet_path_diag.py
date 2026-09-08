"""ADR-95 P3 -- does the H20 corner-branch wall event (P1) TRANSFER to the
tet10 family?  Same Prandtl-Reissner deck as `h20_prandtl.py` / `quad_path_diag.py`
(phi_txc = 20 deg, nu = 0.45, rho_bar = 0, SY = 0.2 kPa apex regulariser, q0
surcharge, rough footing displacement control to s/B, the ADR-63 D16 adaptive
ladder with its subdivision budget/floor), run on `BezierTet10 -bbar`,
`BezierTet10 std` and `TenNodeTetrahedron` instead of the hex family, so the P1
finding (a DruckerPrager corner-branch consistent-tangent event kills the
Newton path, not a loss of ellipticity) can be checked against a genuinely
different element technology rather than just a different hex formulation.

MESH.  `bearing_mesh_tet10.npz` (the file this campaign was first pointed at)
is a DIFFERENT problem -- a 10x10x8 m box under a 2x2 m SQUARE footing with
PDMY sand (ADR-79's bearing_backbone.py), verified by its node coordinate
ranges before writing a line of this file.  It has no q0*Nq oracle and is not
this deck.  `build_mesh_tet10_prandtl.py` (new, this P3) builds the mesh this
script actually needs: the SAME plane-strain strip domain as
`h20_prandtl.strip_mesh` (XLIM=30, ZBOT=-20, B_FOOT=2, THICK=0.5), one element
thick in y, graded on the identical x/z block boundaries at h0=1.0 -- 200
hex-shaped cells, the same cell count as the H20/H8 h0=1.0 legs, each split
into 6 structured tets by gmsh.  2583 nodes / 7749 DOF (vs H20's 4659, H8's
1386 at the same h0) -- tet10 is DOF-expensive per `build_mesh_tet10.py`'s own
note; absolute q is therefore NOT expected to match the hex legs, and the
"matched s/B" table in the P3 results note reads capacity/exact ratios, not
raw q.

THE LOAD IS THE TRAP HERE TOO (BezierTet10.cpp's own docstring, verified
against the code before use): a uniform surface traction is a DIFFERENT nodal
load vector on the two node-conjugate bases that happen to share this mesh.
    BERNSTEIN (BezierTet10, either formulation): every quadratic Bernstein
        face function integrates to A/6 -- q*A/6 on EACH of a face's six nodes
        (3 vertices + 3 mid-edges), uniformly.
    LAGRANGE (TenNodeTetrahedron): the standard consistent T6 rule -- 0 at the
        three vertices, q*A/3 at each of the three mid-edge nodes.
Applying the Lagrange-consistent rule to BezierTet10 (or vice versa) is
exactly the load-basis mismatch note 95's plan flags as a known trap (TIMs
item 1) -- it would put an OSCILLATORY, non-conjugate traction on the element
and can drive surface Gauss points into yield the intended traction never
reaches. So each element gets its OWN consistent surcharge vector, built once
per top face and asserted against the resultant q0*A_top to 1e-9 in an
ELASTIC pre-step (control 3, `h20_prandtl.surcharge_step`'s check, tightened
here per the P3 brief).

FOOTING.  uz is prescribed (Transformation-handler `sp`) on EVERY footing
node, vertices and mid-edges alike -- for BezierTet10 a homogeneous/uniform
control value IS the physical field there too (BezierTet10.cpp: "homogeneous
fixes are equivalent -- zero/uniform control values = uniform field" for a
RIGID footing settling by a single scalar s), so the same `sp(node,3,1.0)`
pattern under one `LoadControl(-ds)` integrator is correct for both bases.

`--branch` reuses `quad_path_diag.py`'s `sample_branch` / `sample_tangent` /
`tangent_health` / `capture_newton_norms` / `collect_top_nodeunbalance` BY
IMPORT -- that file is not modified (another agent's runs depend on its own
CSV/NPZ layout staying put) -- so the branch-census CSV/NPZ column layout is
IDENTICAL to `quad_path_diag.py`'s (`_BRANCH_COLS`, `qpd_<tag>_branch.npz`'s
`stNNN_<field>` station-array convention), just written under a `tpd_` prefix.

Run:
    py -3.12 tet_path_diag.py --elem beziertet10bbar --branch --cond-at 5e-4 --suffix _p3
    py -3.12 tet_path_diag.py --elem tet10 --branch --cond-at 5e-4 --suffix _p3
    py -3.12 tet_path_diag.py --elem beziertet10 --branch --cond-at 5e-4 --suffix _p3
"""
import argparse
import csv
import json
import math
import os
import sys
import time
from collections import deque

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import h20_prandtl as HP                                        # noqa: E402
import quad_path_diag as QPD                                    # noqa: E402 -- NOT modified; reused by import only

ops = HP.ops
SQ23 = math.sqrt(2.0 / 3.0)
log = HP.log

MESH_NPZ = os.path.join(HERE, "bearing_mesh_tet10_prandtl.npz")

_CAPACITY_MODES = QPD._CAPACITY_MODES
PLATEAU_FRAC = QPD.PLATEAU_FRAC
FREE_ADVANCE_FLOOR_FACTOR = QPD.FREE_ADVANCE_FLOOR_FACTOR
DROP_TRUNCATE = QPD.DROP_TRUNCATE
COND_TRIGGER = QPD.COND_TRIGGER

ELEMS = {
    "tet10":           dict(basis="lagrange", label="TenNodeTetrahedron"),
    "beziertet10":      dict(basis="bernstein", bbar=False, label="BezierTet10 std"),
    "beziertet10bbar":  dict(basis="bernstein", bbar=True,  label="BezierTet10 -bbar"),
}
GP_COUNTS = {k: 4 for k in ELEMS}       # verified against SRC (NGAUSS=4 / NumGaussPoints=4)

# tet10 face -> local node indices (3 vertices, 3 mid-edges), TenNodeTetrahedron/
# BezierTet10 node order (v1234, e12 e23 e13 e14 e34 e24) -- IDENTICAL table to
# build_mesh_tet10.py's (same convention, verified there and again here by the
# mesh builder's own edge-midpoint assertion).
TRI_FACES = [(0, 1, 2, 4, 5, 6), (0, 1, 3, 4, 9, 7),
             (1, 2, 3, 5, 8, 9), (0, 2, 3, 6, 8, 7)]


# ---------------------------------------------------------------------------
def load_mesh():
    if not os.path.exists(MESH_NPZ):
        raise SystemExit(
            f"{MESH_NPZ} is missing. Build it with the opensees_env "
            f"interpreter: C:\\Users\\nmb\\venv\\opensees_env\\Scripts\\python.exe "
            f"build_mesh_tet10_prandtl.py")
    d = np.load(MESH_NPZ)
    nodes, tets = d["nodes"], d["tets"]
    sets = {k[4:]: d[k] for k in d.files if k.startswith("set_")}
    meta = dict(h0=float(d["h0"]), b_foot=float(d["b_foot"]),
                thick=float(d["thick"]), xlim=float(d["xlim"]),
                zbot=float(d["zbot"]))
    assert meta["b_foot"] == HP.B_FOOT and meta["thick"] == HP.THICK \
        and meta["xlim"] == HP.XLIM and meta["zbot"] == HP.ZBOT, \
        "the tet mesh was built at different deck parameters than h20_prandtl.py"
    return nodes, tets, sets, meta


def face_area(nodes, face_local, tet):
    a, b, c = (nodes[tet[face_local[i]]] for i in range(3))
    return 0.5 * float(np.linalg.norm(np.cross(b - a, c - a)))


def consistent_surcharge_tet(nodes, tets, basis):
    """int t . B_a dA over every z=0 face, EITHER basis.

    basis='bernstein': q*A/6 to ALL SIX face nodes (BezierTet10.cpp).
    basis='lagrange'  : 0 at the 3 vertices, q*A/3 at the 3 mid-edges (the
                        standard T6 consistent rule for a uniform traction).
    Returns per-node weights w (m^2, the same convention as
    h20_prandtl.consistent_surcharge -- multiply by -Q0 to get the load).
    """
    w = np.zeros(len(nodes))
    nfaces = 0
    for tet in tets:
        for face in TRI_FACES:
            zc = nodes[tet[list(face)], 2]
            if np.any(np.abs(zc) > 1e-9):
                continue
            nfaces += 1
            A = face_area(nodes, face, tet)
            if basis == "bernstein":
                for loc in face:
                    w[tet[loc]] += A / 6.0
            else:
                for loc in face[3:]:            # mid-edges only
                    w[tet[loc]] += A / 3.0
    return w, nfaces


def verify_surcharge_tet(nodes, w, nfaces, basis, area_exact):
    tot = float(w.sum())
    err = abs(tot - area_exact) / area_exact
    log(f"    [control 2] {basis} tet10 surcharge over {nfaces} top faces: "
        f"sum {tot:.9f} m^2 vs {area_exact:.9f} exact (rel {err:.2e})")
    assert err < 1e-9, (tot, area_exact)


# ---------------------------------------------------------------------------
def build_model(nodes, tets, sets, w, elem_key, assoc):
    spec = ELEMS[elem_key]
    alpha = HP.alpha_from_phi_txc(HP.PHI_TXC)
    rho = math.sqrt(2.0) * alpha
    g_el = 3.0 * HP.K_EL * (1.0 - 2.0 * HP.NU) / (2.0 * (1.0 + HP.NU))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))
    ops.nDMaterial("DruckerPrager", 1, HP.K_EL, g_el, HP.SY, rho,
                   rho if assoc else 0.0,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for e, conn in enumerate(tets, start=1):
        tags = [int(c) + 1 for c in conn]
        if elem_key == "tet10":
            ops.element("TenNodeTetrahedron", e, *tags, 1)
        else:
            args = ["BezierTet10", e, *tags, 1]
            if spec["bbar"]:
                args.append("-bbar")
            ops.element(*args)
    bt, xf = set(sets["bottom"].tolist()), set(sets["xface"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if n in xf else 0, 1, 0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if abs(w[n]) > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -HP.Q0 * float(w[n]))
    return alpha, rho, g_el


def build_model_elastic(nodes, tets, sets, w):
    """Control 3: ElasticIsotropic on the same mesh/element/BCs, for the
    reaction-resultant check ONLY (no stress-distribution claim is made for
    the tet basis -- unlike h20_prandtl's control 1, which needs an EXACT
    per-GP 1-D state that a graded tet mesh with two different node-conjugate
    bases is not asserted to reproduce here)."""
    e_el = 9.0 * HP.K_EL * (3.0 * HP.K_EL * (1.0 - 2.0 * HP.NU)
                           / (2.0 * (1.0 + HP.NU))) / (
        3.0 * HP.K_EL + 3.0 * HP.K_EL * (1.0 - 2.0 * HP.NU) / (2.0 * (1.0 + HP.NU)))
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))
    ops.nDMaterial("ElasticIsotropic", 1, e_el, HP.NU, 0.0)
    return e_el


# ---------------------------------------------------------------------------
def run(args):
    if ops is None:
        raise SystemExit("H20_NO_ENGINE is set; a leg needs the engine")
    spec = ELEMS[args.elem]
    basis = spec["basis"]
    assoc = args.assoc

    alpha = HP.alpha_from_phi_txc(HP.PHI_TXC)
    rho = math.sqrt(2.0) * alpha
    mc = HP.mc_from_cone(alpha)
    nq = HP.n_q(mc["ps"])
    q_exact = HP.Q0 * nq
    k0 = HP.NU / (1.0 - HP.NU)
    m0 = (1.0 - k0) / (math.sqrt(3.0) * alpha * (1.0 + 2.0 * k0))

    stamp = ops.ladrunoBuild()
    log(f"[control 0] engine {ops.__file__}")
    log(f"[control 0] ladrunoBuild() = {stamp}")

    nodes, tets, sets, meta = load_mesh()
    ndof = 3 * len(nodes)
    area_top = 2.0 * HP.XLIM * HP.THICK
    w, nfaces = consistent_surcharge_tet(nodes, tets, basis)
    verify_surcharge_tet(nodes, w, nfaces, basis, area_top)

    tag = f"{args.elem}_h{str(meta['h0']).replace('0.', '')}" \
          f"{'_assoc' if assoc else ''}{args.suffix}"
    log(f"=== leg {tag}: {spec['label']} ({basis} surcharge basis), "
        f"tet mesh h0 = {meta['h0']} m, {len(nodes)} nodes / {ndof} DOF, "
        f"{len(tets)} tet10 elements (vs H20 h0=1.0: 4659 DOF; H8 h0=1.0: "
        f"1386 DOF -- tet10 is DOF-expensive per note, not directly comparable)")
    log(f"    ORACLE q_u = q0*N_q = {q_exact:.3f} kPa (phi_ps = {mc['ps']:.3f} deg); "
        f"{'ASSOCIATED (control)' if assoc else 'non-associated (gate)'}")
    log(f"    initial 1-D state m0 = {m0:.4f} of yield"
        + ("  (SAFE)" if m0 < 0.8 else "  *** VOID ***"))
    assert m0 < 0.8

    # --- control 3: elastic reaction resultant, THIS mesh/basis, 1e-9 -------
    build_model_elastic(nodes, tets, sets, w)
    for i in range(1, len(nodes) + 1):
        ops.fix(i, 0, 0, 0)
    bt = set(sets["bottom"].tolist())
    xf = set(sets["xface"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if n in xf else 0, 1, 0)
    for e, conn in enumerate(tets, start=1):
        tags = [int(c) + 1 for c in conn]
        if args.elem == "tet10":
            ops.element("TenNodeTetrahedron", e, *tags, 1)
        else:
            eargs = ["BezierTet10", e, *tags, 1]
            if spec.get("bbar"):
                eargs.append("-bbar")
            ops.element(*eargs)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if abs(w[n]) > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -HP.Q0 * float(w[n]))
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso")
    except Exception:
        ops.system("UmfPack")
    ops.test("NormUnbalance", 1e-6 * HP.Q0 * area_top, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise SystemExit("elastic control-3 pre-step failed to converge")
    ops.reactions()
    rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    want = HP.Q0 * area_top
    err3 = abs(rz / want - 1.0)
    log(f"    [control 3] elastic pre-step reaction {rz:.9f} vs {want:.9f} kN "
        f"({100 * (rz / want - 1):+.3e} %)  -> {'PASS' if err3 < 1e-9 else '*** FAIL ***'}")
    assert err3 < 1e-9, "surcharge resultant does not match q0*A_top to 1e-9"

    # --- the real (plastic) model ---------------------------------------
    build_model(nodes, tets, sets, w, args.elem, assoc)
    tol = 1.0e-5 * max(HP.Q0 * float(w.sum()), 1.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    try:
        ops.system("Pardiso")
        sysname = "Pardiso"
    except Exception:
        ops.system("UmfPack")
        sysname = "UmfPack"
    ops.test("NormUnbalance", tol, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    t0 = time.time()
    if ops.analyze(1) != 0:
        raise SystemExit("surcharge step failed")
    ops.reactions()
    rz2 = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    want2 = HP.Q0 * float(w.sum())
    log(f"    [surcharge] plastic-model equilibrium {rz2:.6f} vs {want2:.6f} kN "
        f"({100 * (rz2 / want2 - 1):+.6f} %)   [{sysname}, {time.time() - t0:.1f}s]")
    assert abs(rz2 / want2 - 1) < 1e-6

    mob0_mean, mob0_max, ngp = QPD.mobilisation(len(tets), rho)
    log(f"    [control M0] elements at mob >= {QPD.MOB_YIELD} before the push: "
        f"{int((mob0_mean >= QPD.MOB_YIELD).sum())} of {len(tets)}  (must be 0), "
        f"{ngp} GP/element (expect 4)")
    assert ngp == 4, f"unexpected GP count {ngp} (expected 4 for tet10)"

    # ---------------- the push ---------------------------------------------
    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = HP.Q0 * float(w[sets["footing"]].sum())
    area = HP.B_FOOT * HP.THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = args.sfrac * HP.B_FOOT
    ops.loadConst("-time", uz0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)

    base, floor, dmax = HP.DS_LADDER[2]      # the quadratic-order D16 ladder

    def setup(sysname):
        ops.wipeAnalysis()
        ops.constraints("Transformation")
        ops.numberer("Plain" if sysname == "FullGeneral" else "RCM")
        if sysname == "FullGeneral":
            ops.system("FullGeneral")
        else:
            try:
                ops.system("Pardiso")
            except Exception:
                ops.system("UmfPack")
        ops.analysis("Static")

    setup("sparse")
    ladder = [("KrylovNewton", tol, 25, 0),
              ("NewtonLineSearch", tol, 40, 0),
              ("KrylovNewton", 10.0 * tol, 60, 1)]

    log(f"    ALLOWANCES: ds base/floor/max = {base * 1e3:.4g}/{floor * 1e3:.4g}/"
        f"{dmax * 1e3:.4g} mm, subdivision budget = {args.budget}, s/B target = "
        f"{args.sfrac}, wall cap = {args.tmax:.0f} s, cond-at = "
        f"{(args.cond_at if args.cond_at is not None else COND_TRIGGER * base) * 1e3:.4g} mm")

    cond_at = args.cond_at if args.cond_at is not None else COND_TRIGGER * base

    path = os.path.join(HERE, f"tpd_{tag}.csv")
    fh = open(path, "w", newline="")
    wr = csv.writer(fh)
    wr.writerow(["s_m", "s_over_B", "q_kPa", "ds_mm", "relaxed", "wall_s"])
    rows, ds, good, nfail, nrelax, nsub = [], base, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"

    cen = nodes[tets[:, :4]].mean(axis=1)          # (Ntet, 3), for GP location tags

    branch_fh = branch_wr = None
    branch_on, next_dense = False, None
    dense_lo = dense_hi = dense_step = None
    n_branch_empty = n_branch_seen = 0
    branch_station_s, branch_station_at_wall, branch_station_arrays = [], [], []
    branch_ngp = GP_COUNTS[args.elem]
    _BRANCH_COLS = (
        ["s_over_B", "q_kPa", "s_min_rel", "cond", "lam_min_sym_rel",
         "n_gp_seen", "n_gp_empty", "n_branch0", "n_branch1", "n_branch2",
         "n_branch3", "n_forced", "detAmin_min", "n_detAmin_nonpos", "I1_min",
         "gamma0_max", "gamma1_max"]
        + [f"lowdet{k}_{f}" for k in range(1, 6)
           for f in ("ele", "gp", "val", "x", "z")]
        + [f"corner{k}_{f}" for k in range(1, 6)
           for f in ("ele", "gp", "x", "z")])
    if args.branch:
        if args.dense is not None:
            dense_lo, dense_hi, dense_step = args.dense
            next_dense = dense_lo
        branch_path = os.path.join(HERE, f"tpd_{tag}_branch.csv")
        branch_fh = open(branch_path, "w", newline="")
        branch_fh.write(
            "# ladrunoBranch census, one row per sampling station (tet10 "
            "family, ADR-95 P3). SAME column layout as quad_path_diag.py's "
            "qpd_<tag>_branch.csv.\n"
            "# station = cadence (every "
            f"{args.cond_every} converged steps once ds < {cond_at * 1e3:.6g} "
            "mm)" + (f" UNIONED with --dense every {dense_step} of s/B in "
                     f"[{dense_lo}, {dense_hi}]" if dense_step else "") + ".\n")
        branch_wr = csv.writer(branch_fh)
        branch_wr.writerow(_BRANCH_COLS)

    def write_branch_station(s_over_B, q_now, th, at_wall=False):
        nonlocal n_branch_empty, n_branch_seen
        bsum, barr, top5, corner5 = QPD.sample_branch(len(tets), branch_ngp, cen)
        n_branch_empty += bsum["n_gp_empty"]
        n_branch_seen += bsum["n_gp_seen"]
        row = [s_over_B, q_now,
               th["s_min_rel"] if th else float("nan"),
               th["cond"] if th else float("nan"),
               th["lam_min_sym_rel"] if th else float("nan"),
               bsum["n_gp_seen"], bsum["n_gp_empty"], bsum["n_branch0"],
               bsum["n_branch1"], bsum["n_branch2"], bsum["n_branch3"],
               bsum["n_forced"], bsum["detAmin_min"], bsum["n_detAmin_nonpos"],
               bsum["I1_min"], bsum["gamma0_max"], bsum["gamma1_max"]]
        for k in range(5):
            if k < len(top5):
                detA, e, gp, x, z = top5[k]
                row += [e, gp, detA, x, z]
            else:
                row += [float("nan")] * 5
        for k in range(5):
            if k < len(corner5):
                e, gp, x, z = corner5[k]
                row += [e, gp, x, z]
            else:
                row += [float("nan")] * 4
        branch_wr.writerow([f"{v:.9g}" for v in row])
        branch_fh.flush()
        branch_station_s.append(s_over_B)
        branch_station_at_wall.append(at_wall)
        branch_station_arrays.append(barr)
        return bsum

    t0 = time.time()
    while True:
        s_now = uz0 - ops.getTime()
        if s_now >= smax - 1e-12:
            break
        if time.time() - t0 > args.tmax:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now / HP.B_FOOT:.5f}"
            break
        ds = min(ds, smax - s_now)
        if args.branch and not branch_on and ds < cond_at:
            branch_on = True
            log(f"    [branch] ds fell to {ds * 1e3:.6g} mm (trigger "
                f"{cond_at * 1e3:.6g} mm) at s/B = {s_now / HP.B_FOOT:.5f}: "
                f"arming the branch-census cadence sampler")
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
            if nsub > args.budget:
                mode = "BUDGET"
                verdict = (f"subdivision budget of {args.budget} spent at "
                           f"s/B = {s_now / HP.B_FOOT:.5f}")
                break
            if ds < floor:
                mode = "FLOOR"
                verdict = (f"step collapsed to the {floor * 1e3:.4g} mm floor at "
                           f"s/B = {s_now / HP.B_FOOT:.5f}")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= HP.GROW_AFTER and ds < dmax:
            ds, good = min(2 * ds, dmax), 0
        ops.reactions()
        q = (-sum(ops.nodeReaction(t, 3) for t in foot) + r_corr) / area
        s = uz0 - ops.getTime()
        rows.append((s, s / HP.B_FOOT, q, ds * 1e3, relaxed, time.time() - t0))
        wr.writerow([f"{v:.9g}" for v in rows[-1]])
        fh.flush()

        th_this_step = None
        if args.branch:
            s_ob = s / HP.B_FOOT
            do_dense = next_dense is not None and s_ob >= next_dense
            do_cadence = branch_on and len(rows) % args.cond_every == 0
            if do_dense or do_cadence:
                th_this_step = QPD.sample_tangent(setup, nsub)
                bsum = write_branch_station(s_ob, q, th_this_step)
                log(f"    [branch] s/B {s_ob:.6f}  q {q:8.2f}  "
                    f"GP seen/empty {bsum['n_gp_seen']}/{bsum['n_gp_empty']}  "
                    f"branch(0,1,2,3)=({bsum['n_branch0']},{bsum['n_branch1']},"
                    f"{bsum['n_branch2']},{bsum['n_branch3']})  forced "
                    f"{bsum['n_forced']}  detAmin_min {bsum['detAmin_min']:.3e}")
                if do_dense:
                    while next_dense is not None and s_ob >= next_dense:
                        next_dense += dense_step
                        if next_dense > dense_hi + 1e-12:
                            next_dense = None
        if len(rows) > HP.STALL_WINDOW and \
                rows[-1][0] - rows[-1 - HP.STALL_WINDOW][0] < HP.STALL_ADVANCE * smax:
            mode = "STALL"
            verdict = (f"stalled at s/B = {s / HP.B_FOOT:.5f}")
            break
    fh.close()
    wall = time.time() - t0
    assert len(rows) > 4, f"only {len(rows)} converged steps ({verdict})"

    th_final = QPD.sample_tangent(setup, nsub)
    if th_final:
        th_final.update(s_over_B=float(rows[-1][1]), q=float(rows[-1][2]), at_wall=True)
        log(f"    [diag 3] AT THE WALL, s/B {th_final['s_over_B']:.5f}: "
            f"sigma_min/scale {th_final['s_min_rel']:.3e}, cond {th_final['cond']:.3e}, "
            f"lam_min(sym)/scale {th_final['lam_min_sym_rel']:.3e}, negative "
            f"eigenvalues of the symmetric part: {th_final['n_neg_sym']} of {th_final['n']}")

    if args.branch:
        bsum = write_branch_station(float(rows[-1][1]), float(rows[-1][2]),
                                     th_final, at_wall=True)
        log(f"    [branch] AT THE WALL, s/B {rows[-1][1]:.6f}: GP seen/empty "
            f"{bsum['n_gp_seen']}/{bsum['n_gp_empty']}  branch(0,1,2,3)="
            f"({bsum['n_branch0']},{bsum['n_branch1']},{bsum['n_branch2']},"
            f"{bsum['n_branch3']})  forced {bsum['n_forced']}  detAmin_min "
            f"{bsum['detAmin_min']:.3e}")
        branch_fh.close()
        np.savez(os.path.join(HERE, f"tpd_{tag}_branch.npz"),
                 station_s_over_B=np.array(branch_station_s),
                 station_at_wall=np.array(branch_station_at_wall),
                 **{f"st{i:03d}_{k}": v
                    for i, arrs in enumerate(branch_station_arrays)
                    for k, v in arrs.items()})
        log(f"    [branch] {len(branch_station_s)} station(s) written to "
            f"tpd_{tag}_branch.csv/.npz; {n_branch_seen} live GP responses, "
            f"{n_branch_empty} empty")

    a = np.array([r[:4] for r in rows])
    s, q, dsm = a[:, 0], a[:, 2], a[:, 3]
    peak = np.maximum.accumulate(q)
    bad = np.where(q < (1.0 - DROP_TRUNCATE) * peak)[0]
    if len(bad) and bad[0] > 4:
        mode = "TRUNCATED"
        verdict += f"; curve TRUNCATED ({len(q) - bad[0]} steps cut)"
        s, q, dsm = s[:bad[0]], q[:bad[0]], dsm[:bad[0]]

    msk = s >= 0.9 * s[-1]
    t_last = float(np.polyfit(s[msk], q[msk], 1)[0]) if msk.sum() > 2 else float("nan")
    n0 = max(4, len(s) // 50)
    t_init = float(np.polyfit(s[:n0], q[:n0], 1)[0])
    plateau = bool(abs(t_last) < PLATEAU_FRAC * abs(t_init))
    qmax = float(q.max())
    floor_mm = floor * 1e3
    ds_tail_min = float(dsm[msk].min()) if msk.sum() else float(dsm[-1])
    headroom = ds_tail_min / floor_mm
    free = bool(headroom >= FREE_ADVANCE_FLOOR_FACTOR)
    capacity = bool(plateau and free and mode in _CAPACITY_MODES)

    log("")
    log(f"    --- termination -------------------------------------------------")
    log(f"    {len(rows)} steps, {nfail} failed attempts, {nsub}/{args.budget} "
        f"subdivisions, {nrelax} relaxed, {wall:.0f} s")
    log(f"    MODE = {mode}   [{verdict}]")
    log(f"    q_max = {qmax:.2f} kPa = {qmax / q_exact:.4f} of exact "
        f"({q_exact:.2f} kPa); end s/B = {s[-1] / HP.B_FOOT:.5f} of {args.sfrac}")
    log(f"    tail dq/ds = {t_last:.1f} kPa/m = {100 * t_last / t_init:.3f} % of "
        f"initial -> {'PLATEAU' if plateau else 'STILL HARDENING'}")
    log(f"    CAPACITY = {'YES' if capacity else 'NO'}"
        + ("" if capacity else f"  -> ALLOWANCE named {mode}, not a capacity"))

    res = dict(tag=tag, elem=args.elem, label=spec["label"], basis=basis,
               bbar=spec.get("bbar"), h0=meta["h0"], assoc=assoc, build=stamp,
               nel=int(len(tets)), ndof=ndof, ngp=int(ngp),
               qmax=qmax, q_exact=q_exact, ratio=qmax / q_exact,
               mode=mode, verdict=verdict, plateau=plateau, free=free,
               capacity=capacity, tail_pct=100 * t_last / t_init,
               s_end_over_B=float(s[-1]) / HP.B_FOOT, s_target=args.sfrac,
               nsub=nsub, budget=args.budget, steps=len(rows), nfail=nfail,
               nrelax=nrelax, wall_s=wall, ladder="quad")
    with open(os.path.join(HERE, f"tpd_{tag}.json"), "w") as f:
        json.dump(res, f, indent=1, default=float)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--elem", default="beziertet10bbar", choices=sorted(ELEMS))
    ap.add_argument("--h0", type=float, default=None,
                    help="informational only -- the tet mesh's own h0 (from "
                         "the npz) is authoritative; a mismatch is asserted")
    ap.add_argument("--sy", type=float, default=None)
    ap.add_argument("--sfrac", type=float, default=0.15)
    ap.add_argument("--assoc", action="store_true")
    ap.add_argument("--budget", type=int, default=200)
    ap.add_argument("--tmax", type=float, default=5400.0)
    ap.add_argument("--cond-every", type=int, default=4)
    ap.add_argument("--cond-at", type=float, default=None, metavar="M")
    ap.add_argument("--branch", action="store_true")
    ap.add_argument("--dense", type=float, nargs=3, default=None,
                    metavar=("LO", "HI", "STEP"))
    ap.add_argument("--suffix", default="")
    args = ap.parse_args()
    if args.sy is not None:
        HP.SY = args.sy
        print(f"[adr95] SY override -> {HP.SY} kPa")
    r = run(args)
    print()
    print(f"{r['tag']:>28} {r['elem']:>16} mode={r['mode']:<9} "
          f"ratio={r['ratio']:.4f} tail={r['tail_pct']:.3f}% "
          f"sub={r['nsub']}/{r['budget']} CAP={'yes' if r['capacity'] else 'NO'}")


if __name__ == "__main__":
    main()
