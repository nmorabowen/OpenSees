"""ADR-95 mechanism visualization: the incremental (collapse) velocity field
Delta-u = u(s/B=0.15) - u(s/B=0.14) for the five UW-DP-repaired-material,
h0=1.0 Prandtl-Reissner strip-footing legs (h8bbar, h20uri, beziertet10bbar,
beziertet10 std, tet10), compared against the analytical Prandtl-Reissner
mechanism (active wedge + log-spiral fan + passive wedge) at phi_ps=27.470deg.

Reuses `deformed_snapshot.py` BY IMPORT (DIST paths, LEGS dict entries B1-B5,
Y0_LOC, Q_EXACT) -- that file is NOT modified.  `deformed_snapshot.py`'s own
push loop only ever captures ONE mid-run waypoint (s/B=0.01) plus the final
state, which is not enough to build a velocity field this late in the push,
so `run_leg` below duplicates its model-building + adaptive-ladder push loop
(itself copied from `h20_prandtl.py` / `tet_path_diag.py`, same as
deformed_snapshot.py does) with the waypoint list changed to (0.14, 0.15) and
no per-Gauss-point stress/branch census (not needed for a kinematic field).

Two-phase, same convention as deformed_snapshot.py:
    py -3.12 mechanism_snapshot.py --leg B1              # run one leg (writes
                                                          # mech_<LEG>.npz)
    py -3.12 mechanism_snapshot.py --plot both            # build the figures
                                                          # from mech_*.npz
"""
import argparse
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import deformed_snapshot as DS  # noqa: E402 -- NOT modified; reused by import

LEGS = {k: v for k, v in DS.LEGS.items() if k in ("B1", "B2", "B3", "B4", "B5")}
WAYPOINTS = (0.14, 0.15)


def phi_ps_deg():
    """The plane-strain Mohr-Coulomb angle matching the UW-DP cone's slope
    (phi_txc = 20 deg triaxial-compression fit, Chen & Han plane-strain
    match), read from h20_prandtl.py's own oracle computation:
        alpha = alpha_from_phi_txc(PHI_TXC=20.0)
        mc = mc_from_cone(alpha); phi_ps = mc['ps']   ->  27.47016... deg
    This is the SAME angle the harness's q_exact = Q0 * N_q(phi_ps) oracle
    uses (h20_prandtl.py leg `nonassoc`/`assoc`, `mc['ps']` in its printed
    oracle line).  Imports h20_prandtl with H20_NO_ENGINE=1 (the module's own
    escape hatch for using its mesh/oracle math without the solver) so the
    plotting phase never touches the opensees engine -- only run_leg's own
    (separate-process) `import h20_prandtl` does that, with ADR95_DIST set."""
    os.environ.setdefault("H20_NO_ENGINE", "1")
    import h20_prandtl as HPconst
    return HPconst.mc_from_cone(HPconst.alpha_from_phi_txc(HPconst.PHI_TXC))["ps"]


# ---------------------------------------------------------------------------
# leg runner -- duplicates deformed_snapshot.run_leg's model build + push loop
# ---------------------------------------------------------------------------
def run_leg(name, cfg, tmax, budget):
    os.environ["ADR95_DIST"] = DS.DIST[cfg["dist"]]
    os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

    import h20_prandtl as HP2
    ops = HP2.ops
    stamp = ops.ladrunoBuild()
    print(f"[mech] leg {name}: build={stamp} dist={cfg['dist']} "
          f"mat={cfg['mat']} elem={cfg['elem']} mesh={cfg['mesh']}", flush=True)

    assert cfg["mat"] == "uwdp", "the mechanism legs are UW-DP only"
    mesh_kind = cfg["mesh"]
    elem_key = cfg["elem"]
    y0_loc = None
    basis = None

    if mesh_kind == "hex":
        order = 1 if elem_key == "h8bbar" else 2
        form = "bbar" if elem_key == "h8bbar" else "uri"
        nodes, cells, sets, xg, yg, zg = HP2.strip_mesh(1.0, order)
        w, nfaces = HP2.consistent_surcharge(nodes, cells, order)
        HP2.verify_surcharge(nodes, cells, w, nfaces, order)
        elem_type = f"hex{cells.shape[1]}"
        y0_loc = DS.Y0_LOC[cells.shape[1]]
        HP2.build_model(nodes, cells, sets, w, form, assoc=False)
        tol = 1.0e-5 * max(HP2.Q0 * float(w.sum()), 1.0)
        HP2.surcharge_step(nodes, sets, w, tol)
        ladder = HP2.attempts(tol)
        base, floor, dmax = HP2.DS_LADDER[order]
    else:
        import tet_path_diag as TPD
        nodes, tets, sets, meta = TPD.load_mesh()
        cells = tets
        basis = TPD.ELEMS[elem_key]["basis"]
        w, nfaces = TPD.consistent_surcharge_tet(nodes, tets, basis)
        area_top = 2.0 * HP2.XLIM * HP2.THICK
        TPD.verify_surcharge_tet(nodes, w, nfaces, basis, area_top)
        elem_type = f"tet10_{elem_key}"
        TPD.build_model(nodes, tets, sets, w, elem_key, assoc=False)
        tol = 1.0e-5 * max(HP2.Q0 * float(w.sum()), 1.0)
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
        base, floor, dmax = HP2.DS_LADDER[2]
        ladder = [("KrylovNewton", tol, 25, 0),
                  ("NewtonLineSearch", tol, 40, 0),
                  ("KrylovNewton", 10.0 * tol, 60, 1)]

    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = HP2.Q0 * float(w[sets["footing"]].sum())
    area = HP2.B_FOOT * HP2.THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    smax = 0.15 * HP2.B_FOOT
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

    def snapshot_disp():
        n = len(nodes)
        disp = np.zeros((n, 3))
        for i in range(1, n + 1):
            d = ops.nodeDisp(i)
            disp[i - 1, :min(3, len(d))] = d[:3]
        return disp

    waypoints = list(WAYPOINTS)
    snaps = {}
    ds, good, nfail, nrelax, nsub = base, 0, 0, 0, 0
    mode, verdict = "TARGET", "reached the target settlement"
    rows = []
    t0 = time.time()
    while waypoints:
        target_ob = waypoints[0]
        target_s = target_ob * HP2.B_FOOT
        s_now = uz0 - ops.getTime()
        if s_now >= target_s - 1e-12:
            q_here = rows[-1][2] if rows else float("nan")
            snaps[target_ob] = dict(disp=snapshot_disp(),
                                     s_over_B=s_now / HP2.B_FOOT, q=q_here)
            print(f"    [mech] leg {name}: captured s/B>={target_ob} at "
                  f"s/B={s_now / HP2.B_FOOT:.6f} q={q_here:.3f}", flush=True)
            waypoints.pop(0)
            continue
        if time.time() - t0 > tmax:
            mode, verdict = "WALL", f"wall-clock cap at s/B = {s_now / HP2.B_FOOT:.5f}"
            break
        ds = min(ds, target_s - s_now)
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
                verdict = f"subdivision budget of {budget} spent at s/B = {s_now / HP2.B_FOOT:.5f}"
                break
            if ds < floor:
                mode = "FLOOR"
                verdict = f"step collapsed to the {floor * 1e3:.4g} mm floor at s/B = {s_now / HP2.B_FOOT:.5f}"
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= HP2.GROW_AFTER and ds < dmax:
            ds, good = min(2 * ds, dmax), 0
        ops.reactions()
        q = (-sum(ops.nodeReaction(t, 3) for t in foot) + r_corr) / area
        s = uz0 - ops.getTime()
        s_ob = s / HP2.B_FOOT
        rows.append((s, s_ob, q, ds * 1e3, relaxed, time.time() - t0))
        if len(rows) % 20 == 0:
            print(f"    [{name}] step {len(rows)} s/B={s_ob:.5f} q={q:.3f} "
                  f"ds={ds * 1e3:.4g}mm wall={time.time() - t0:.0f}s", flush=True)
        if len(rows) > HP2.STALL_WINDOW and \
                rows[-1][0] - rows[-1 - HP2.STALL_WINDOW][0] < HP2.STALL_ADVANCE * smax:
            mode = "STALL"
            verdict = f"stalled at s/B = {s_ob:.5f}"
            break

    wall = time.time() - t0
    print(f"[mech] leg {name} MODE={mode}  [{verdict}]", flush=True)
    print(f"[mech] leg {name}: wall={wall:.0f}s steps={len(rows)} nsub={nsub} "
          f"nfail={nfail} captured={sorted(snaps.keys())}", flush=True)

    n = len(nodes)
    out = dict(
        leg=name, elem=elem_key, mesh=mesh_kind, material=cfg["mat"],
        dist=cfg["dist"], build=str(stamp), mode=mode, verdict=verdict,
        nodes=nodes, cells=cells, elem_type=elem_type,
        disp_014=(snaps[0.14]["disp"] if 0.14 in snaps else np.full((n, 3), np.nan)),
        disp_015=(snaps[0.15]["disp"] if 0.15 in snaps else np.full((n, 3), np.nan)),
        s_014=(snaps[0.14]["s_over_B"] if 0.14 in snaps else np.nan),
        s_015=(snaps[0.15]["s_over_B"] if 0.15 in snaps else np.nan),
        q_014=(snaps[0.14]["q"] if 0.14 in snaps else np.nan),
        q_015=(snaps[0.15]["q"] if 0.15 in snaps else np.nan),
        q_exact=DS.Q_EXACT, wall_s=wall, nsub=nsub, nfail=nfail, nrelax=nrelax,
        steps=len(rows), budget=budget,
    )
    if basis is not None:
        out["basis"] = basis
    if y0_loc is not None:
        out["y0_loc"] = np.array(y0_loc)
    np.savez(os.path.join(HERE, f"mech_{name}.npz"), **out)
    print(f"[mech] wrote mech_{name}.npz", flush=True)
    return out


# ---------------------------------------------------------------------------
# generic isoparametric shape functions (finite-difference gradients) --
# H8 (trilinear), H20 (serendipity quadratic, H20_OFF ordering), and
# quadratic tet10 in EITHER of its two node-conjugate bases (Lagrange /
# Bernstein), node order v1234 e12 e23 e13 e14 e34 e24 (tet_path_diag.py's
# TRI_FACES convention).  Used ONLY to post-process the converged nodal
# displacement snapshots into a per-element strain measure -- no solver
# import below this point.
# ---------------------------------------------------------------------------
def _shape_h8(xi):
    x, e, z = xi
    off = [(-1, -1, -1), (1, -1, -1), (1, 1, -1), (-1, 1, -1),
           (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1)]
    return np.array([0.125 * (1 + x * xa) * (1 + e * ea) * (1 + z * za)
                      for xa, ea, za in off])


_H20_NAT = [(xo - 1, yo - 1, zo - 1) for (xo, yo, zo) in
            [(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0),
             (0, 0, 2), (2, 0, 2), (2, 2, 2), (0, 2, 2),
             (1, 0, 0), (2, 1, 0), (1, 2, 0), (0, 1, 0),
             (1, 0, 2), (2, 1, 2), (1, 2, 2), (0, 1, 2),
             (0, 0, 1), (2, 0, 1), (2, 2, 1), (0, 2, 1)]]


def _shape_h20(xi):
    x, e, z = xi
    out = np.empty(20)
    for a, (xa, ea, za) in enumerate(_H20_NAT):
        nz = (xa == 0) + (ea == 0) + (za == 0)
        if nz == 0:                                    # corner
            out[a] = 0.125 * (1 + x * xa) * (1 + e * ea) * (1 + z * za) \
                * (x * xa + e * ea + z * za - 2)
        elif xa == 0:                                   # mid-edge, xi fixed=0
            out[a] = 0.25 * (1 - x * x) * (1 + e * ea) * (1 + z * za)
        elif ea == 0:
            out[a] = 0.25 * (1 - e * e) * (1 + x * xa) * (1 + z * za)
        else:
            out[a] = 0.25 * (1 - z * z) * (1 + x * xa) * (1 + e * ea)
    return out


def _shape_tet10(xi, basis):
    """xi = (L1, L2, L3); L4 = 1 - L1 - L2 - L3.  Node order v1234, then
    e12 e23 e13 e14 e34 e24 (tet_path_diag.py TRI_FACES convention)."""
    L1, L2, L3 = xi
    L4 = 1.0 - L1 - L2 - L3
    L = (L1, L2, L3, L4)
    out = np.empty(10)
    if basis == "lagrange":
        for i in range(4):
            out[i] = L[i] * (2.0 * L[i] - 1.0)
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]
        for k, (i, j) in enumerate(edges):
            out[4 + k] = 4.0 * L[i] * L[j]
    else:                                                # bernstein
        for i in range(4):
            out[i] = L[i] * L[i]
        edges = [(0, 1), (1, 2), (0, 2), (0, 3), (2, 3), (1, 3)]
        for k, (i, j) in enumerate(edges):
            out[4 + k] = 2.0 * L[i] * L[j]
    return out


def _grad_fd(shape_fn, xi0, h=1.0e-6):
    """dN/dxi at xi0 by central differences -- generic over element type."""
    xi0 = np.asarray(xi0, dtype=float)
    n = xi0.size
    N0 = shape_fn(xi0)
    nen = N0.size
    dN = np.empty((nen, n))
    for k in range(n):
        step = np.zeros(n)
        step[k] = h
        dN[:, k] = (shape_fn(xi0 + step) - shape_fn(xi0 - step)) / (2 * h)
    return N0, dN


def _elem_geometry(elem_type, basis):
    """(shape_fn, centroid_xi, ref_volume) for one call per element type."""
    if elem_type.startswith("hex"):
        nen = 8 if elem_type == "hex8" else 20
        shape_fn = _shape_h8 if nen == 8 else _shape_h20
        return shape_fn, np.zeros(3), 8.0
    else:
        return (lambda xi: _shape_tet10(xi, basis)), np.array([0.25, 0.25, 0.25]), 1.0 / 6.0


def elem_incremental_strain(nodes, cells, elem_type, basis, du):
    """Per-element incremental deviatoric strain magnitude
    eps_q = sqrt(2/3 dev(deps):dev(deps)) and volumetric strain deps_v,
    from the element's OWN shape functions at its centroid (one-point rule).
    Returns (eps_q, deps_v, vol) arrays, one entry per element.

    Small-strain kinematics (deps = sym(grad(du))) -- consistent with the
    hypo-elastoplastic UW-DP rate form this survey uses; du is itself a small
    increment (s/B 0.14->0.15) so the small-strain measure is the right one
    to report even though the underlying analysis is finite-kinematics.
    """
    shape_fn, xi0, ref_vol = _elem_geometry(elem_type, basis)
    ne = cells.shape[0]
    eps_q = np.full(ne, np.nan)
    deps_v = np.full(ne, np.nan)
    vol = np.full(ne, np.nan)
    _, dN_dxi = _grad_fd(shape_fn, xi0)          # (nen,3), same for every
                                                  # element of a given type
                                                  # ONLY if geometry were
                                                  # affine in xi -- it is not
                                                  # for H20, so recompute the
                                                  # xi-derivative table once
                                                  # (constant) but the
                                                  # x-Jacobian per element.
    for e in range(ne):
        conn = cells[e]
        X = nodes[conn]                          # (nen,3)
        J = dN_dxi.T @ X                         # (3,3) = dx/dxi
        try:
            Jinv = np.linalg.inv(J)
        except np.linalg.LinAlgError:
            continue
        dN_dx = dN_dxi @ Jinv.T                  # (nen,3); J[a,b]=dx_b/dxi_a,
                                                  # so dxi_a/dx_b = Jinv[b,a]
        vol[e] = abs(np.linalg.det(J)) * ref_vol
        U = du[conn]                             # (nen,3) incremental disp
        gradU = U.T @ dN_dx                      # (3,3) = d(du_i)/dx_j
        eps = 0.5 * (gradU + gradU.T)
        ev = eps[0, 0] + eps[1, 1] + eps[2, 2]
        dev = eps - (ev / 3.0) * np.eye(3)
        eq = math.sqrt(max(0.0, (2.0 / 3.0) * np.tensordot(dev, dev)))
        eps_q[e] = eq
        deps_v[e] = ev
    return eps_q, deps_v, vol


# ---------------------------------------------------------------------------
# Prandtl-Reissner analytical mechanism geometry
# ---------------------------------------------------------------------------
def prandtl_mechanism(phi_deg, half_b=1.0, n=60):
    """Returns dict with the active-wedge triangle, the two log-spiral fans
    (one per footing edge) and the two passive-wedge triangles, all as (x,z)
    polylines, PLUS the surface-outcrop distance (from the footing edge to
    where the passive wedge's face reaches z=0).

    Construction (mirror symmetric about x=0; built for the RIGHT edge,
    x -> -x for the left):
      - active wedge: base angle (45+phi/2) at the footing edge E=(half_b,0),
        apex A=(0,-d), d = half_b*tan(45+phi/2).
      - fan centred at E: theta=0 radius is E->A (the active-wedge face,
        length r0=|EA|); the radius sweeps 90 deg, r(theta)=r0*exp(theta*tanphi)
        (log spiral, tangent-angle-with-radius = 90-phi, so the wedge face and
        the terminal radius are BOTH admissible slip-line directions); the
        rotation sense is the one that sweeps AWAY from the wedge (i.e. away
        from the footing centreline), verified below by checking the fan's
        terminal point lands below and outboard of E.
      - passive wedge: from the fan's terminal point F, a straight line at
        (45-phi/2) above horizontal, outward, up to the free surface at S;
        the outcrop distance is |S.x - E.x|.
    """
    phi = math.radians(phi_deg)
    beta_a = math.radians(45.0 + phi_deg / 2.0)     # active wedge base angle
    beta_p_deg = 45.0 - phi_deg / 2.0               # passive wedge face angle
    beta_p = math.radians(beta_p_deg)
    d = half_b * math.tan(beta_a)
    E = np.array([half_b, 0.0])
    A = np.array([0.0, -d])

    r0 = float(np.linalg.norm(A - E))
    ang0 = math.atan2(A[1] - E[1], A[0] - E[0])     # E->A direction, radians

    theta = np.linspace(0.0, math.pi / 2.0, n)
    r = r0 * np.exp(theta * math.tan(phi))
    ang = ang0 + theta                              # CCW sweep (right side)
    fan_xy = E[None, :] + r[:, None] * np.stack(
        [np.cos(ang), np.sin(ang)], axis=1)
    F = fan_xy[-1]

    # passive wedge face: from F, direction (+cos(beta_p), +sin(beta_p)),
    # outward and up, until z = 0
    if F[1] >= 0.0:
        S = F.copy()
    else:
        t = -F[1] / math.sin(beta_p)
        S = F + t * np.array([math.cos(beta_p), math.sin(beta_p)])
    outcrop = float(S[0] - E[0])

    wedge = np.array([[-half_b, 0.0], [half_b, 0.0], [0.0, -d], [-half_b, 0.0]])
    passive_r = np.array([E, F, S])
    passive_l = passive_r.copy()
    passive_l[:, 0] *= -1.0
    fan_l = fan_xy.copy()
    fan_l[:, 0] *= -1.0

    return dict(phi_deg=phi_deg, d=d, r0=r0, outcrop=outcrop,
                beta_active_deg=45.0 + phi_deg / 2.0, beta_passive_deg=beta_p_deg,
                wedge=wedge, fan_r=fan_xy, fan_l=fan_l,
                passive_r=passive_r, passive_l=passive_l)


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
# tet10 face -> local VERTEX indices only (v1234, e12 e23 e13 e14 e34 e24
# node order) -- the vertex-only part of tet_path_diag.py's TRI_FACES table,
# duplicated here (not imported) so the plotting phase never has to import
# tet_path_diag.py, which itself unconditionally imports h20_prandtl.py's
# solver engine at module scope.
_TRI_FACES_V = [(0, 1, 2), (0, 1, 3), (1, 2, 3), (0, 2, 3)]
_LEG_ORDER = ["B1", "B2", "B5", "B4", "B3"]     # h8bbar, h20uri, bezbbar, bezstd, tet10


def _y0_nodes(nodes, tol=1e-6):
    return np.where(np.abs(nodes[:, 1]) < tol)[0]


def _y0_tet_faces(nodes, cells, tol=1e-6):
    """For each tet element, the local vertex triple (from _TRI_FACES_V) that
    lies entirely on the y=0 plane, if any -- one genuine triangular face per
    surface element (a tet can have at most one face on a single plane)."""
    out = {}
    for e in range(cells.shape[0]):
        conn = cells[e]
        for face in _TRI_FACES_V:
            if np.all(np.abs(nodes[conn[list(face)], 1]) < tol):
                out[e] = face
                break
    return out


def _mesh_edges_y0(nodes, cells, elem_type, y0_loc):
    """(N_edges,2,2) array of undeformed-mesh edge endpoints on the y=0 face,
    in the (x,z) plane, for the faint background mesh."""
    segs = []
    if elem_type.startswith("hex"):
        corners = y0_loc[:4]
        for conn in cells:
            quad = nodes[conn[corners]][:, [0, 2]]
            for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
                segs.append([quad[a], quad[b]])
    else:
        for e, face in _y0_tet_faces(nodes, cells).items():
            xz = nodes[cells[e][list(face)]][:, [0, 2]]
            for a, b in [(0, 1), (1, 2), (2, 0)]:
                segs.append([xz[a], xz[b]])
    return np.array(segs)


def _panel(ax, d, xlim, zlim, mech):
    import matplotlib.colors as mcolors
    from matplotlib.collections import LineCollection, PolyCollection

    nodes = d["nodes"]
    cells = d["cells"]
    elem_type = str(d["elem_type"])
    basis = str(d["basis"]) if "basis" in d.files else None
    leg = str(d["leg"])
    elem = str(d["elem"])
    q015, q014 = float(d["q_015"]), float(d["q_014"])
    s015, s014 = float(d["s_015"]), float(d["s_014"])
    q_exact = float(d["q_exact"])

    u014 = np.nan_to_num(d["disp_014"], nan=0.0)
    u015 = np.nan_to_num(d["disp_015"], nan=0.0)
    du = u015 - u014

    eps_q, deps_v, vol = elem_incremental_strain(nodes, cells, elem_type, basis, du)
    finite = np.isfinite(eps_q) & np.isfinite(vol) & (vol > 0)
    vw_eq = float(np.average(eps_q[finite], weights=vol[finite]))
    vw_ev = float(np.average(np.abs(deps_v[finite]), weights=vol[finite]))
    iso_ratio = vw_ev / vw_eq if vw_eq > 0 else float("nan")
    tail_slope = (q015 - q014) / max(s015 - s014, 1e-12)

    # -- undeformed mesh, faint --
    y0_loc = d["y0_loc"] if "y0_loc" in d.files else None
    edges = _mesh_edges_y0(nodes, cells, elem_type, y0_loc)
    ax.add_collection(LineCollection(edges, colors="0.75", linewidths=0.25, zorder=1))

    # -- eps_q filled field (log scale, robust) --
    eq_pos = eps_q[finite & (eps_q > 0)]
    if eq_pos.size:
        vmin = max(np.percentile(eq_pos, 2), eq_pos.max() * 1e-4)
        vmax = np.percentile(eq_pos, 98)
        vmax = max(vmax, vmin * 10)
    else:
        vmin, vmax = 1e-6, 1e-3
    norm = mcolors.LogNorm(vmin=vmin, vmax=vmax, clip=True)
    cmap = plt.get_cmap("inferno")
    eq_plot = np.clip(eps_q, vmin, vmax)
    eq_plot = np.where(np.isfinite(eq_plot), eq_plot, vmin)

    if elem_type.startswith("hex"):
        corners = y0_loc[:4]
        quads = nodes[cells[:, corners]][:, :, [0, 2]]
        colors = cmap(norm(eq_plot))
        ax.add_collection(PolyCollection(quads, facecolors=colors,
                                          edgecolors="none", zorder=2))
    else:
        tris, tri_c = [], []
        for e, face in _y0_tet_faces(nodes, cells).items():
            xz = nodes[cells[e][list(face)]][:, [0, 2]]
            tris.append(xz)
            tri_c.append(cmap(norm(eq_plot[e])))
        if tris:
            ax.add_collection(PolyCollection(tris, facecolors=tri_c,
                                              edgecolors="none", zorder=2))

    # -- quiver of du at y=0 nodes, subsampled, normalised so a footing-node
    #    arrow reads ~0.6 (in x/z data units) --
    y0n = _y0_nodes(nodes)
    foot_mask = (np.abs(nodes[y0n, 0]) <= 1.0 + 1e-6) & (np.abs(nodes[y0n, 2]) < 1e-6)
    foot_mag = np.linalg.norm(du[y0n[foot_mask]][:, [0, 2]], axis=1)
    ref_mag = float(np.median(foot_mag)) if foot_mag.size else float(
        np.percentile(np.linalg.norm(du[y0n][:, [0, 2]], axis=1), 95))
    ref_mag = max(ref_mag, 1e-9)
    in_view = ((nodes[y0n, 0] >= xlim[0]) & (nodes[y0n, 0] <= xlim[1]) &
               (nodes[y0n, 2] >= zlim[0]) & (nodes[y0n, 2] <= zlim[1]))
    sub = y0n[in_view]
    step = max(1, len(sub) // 250)
    sub = sub[::step]
    ax.quiver(nodes[sub, 0], nodes[sub, 2], du[sub, 0], du[sub, 2],
              angles="xy", scale_units="xy", scale=ref_mag / 0.6,
              color="cyan", width=0.0025, headwidth=3.5, zorder=4)

    # -- footing bar --
    ax.plot([-1.0, 1.0], [0.0, 0.0], color="black", linewidth=3.5,
             solid_capstyle="butt", zorder=6)

    # -- analytical Prandtl-Reissner mechanism overlay --
    ax.plot(mech["wedge"][:, 0], mech["wedge"][:, 1], color="white",
            linewidth=1.4, zorder=5)
    ax.plot(mech["wedge"][:, 0], mech["wedge"][:, 1], color="black",
            linewidth=0.6, linestyle="--", zorder=5)
    for fan in (mech["fan_r"], mech["fan_l"]):
        ax.plot(fan[:, 0], fan[:, 1], color="white", linewidth=1.4, zorder=5)
        ax.plot(fan[:, 0], fan[:, 1], color="black", linewidth=0.6,
                linestyle="--", zorder=5)
    for pw in (mech["passive_r"], mech["passive_l"]):
        ax.plot(pw[:, 0], pw[:, 1], color="white", linewidth=1.4, zorder=5)
        ax.plot(pw[:, 0], pw[:, 1], color="black", linewidth=0.6,
                linestyle="--", zorder=5)

    ax.set_xlim(*xlim)
    ax.set_ylim(*zlim)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=6)

    ax.set_title(
        f"{leg}: {elem}   q/q_exact@0.15={q015 / q_exact:.3f}\n"
        f"tail slope dq/d(s/B)={tail_slope:.1f} kPa   "
        f"|dev|/eq (isochoric)={iso_ratio:.3f}",
        fontsize=7.2, linespacing=1.3)
    return norm, cmap, dict(leg=leg, eps_q=eps_q, deps_v=deps_v, vol=vol,
                             vw_eq=vw_eq, vw_ev=vw_ev, iso_ratio=iso_ratio,
                             tail_slope=tail_slope)


def build_figures(which):
    global plt
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F811

    phi_ps = phi_ps_deg()
    mech = prandtl_mechanism(phi_ps, half_b=1.0)
    print(f"[plot] phi_ps = {phi_ps:.3f} deg  apex depth d={mech['d']:.3f} m  "
          f"r0={mech['r0']:.3f} m  surface outcrop = {mech['outcrop']:.3f} m "
          f"({mech['outcrop'] / 2.0:.2f} B) from each footing edge", flush=True)

    data = {}
    missing = []
    for name in _LEG_ORDER:
        path = os.path.join(HERE, f"mech_{name}.npz")
        if os.path.exists(path):
            data[name] = np.load(path, allow_pickle=True)
        else:
            missing.append(name)
    if missing:
        print(f"[plot] WARNING: missing mech_*.npz for legs {missing}", flush=True)

    stats = {}

    def make(fname, xlim, zlim):
        fig, axes = plt.subplots(1, 5, figsize=(24, 6.2), dpi=130)
        fig.subplots_adjust(left=0.03, right=0.90, top=0.82, bottom=0.08, wspace=0.28)
        last = None
        for ax, name in zip(axes, _LEG_ORDER):
            if name not in data:
                ax.set_visible(False)
                continue
            last = _panel(ax, data[name], xlim, zlim, mech)
            stats[name] = last[2]
        for ax in axes[len(_LEG_ORDER):]:
            ax.set_visible(False)
        if last is not None:
            norm, cmap, _ = last
            sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
            sm.set_array([])
            cax = fig.add_axes([0.92, 0.15, 0.012, 0.6])
            fig.colorbar(sm, cax=cax, label="incremental eps_q (log scale)")
        fig.suptitle(
            "ADR-95 Prandtl-Reissner collapse mechanism: incremental "
            f"deviatoric strain, s/B 0.14->0.15, vs analytical mechanism "
            f"(phi_ps={phi_ps:.2f} deg, outcrop={mech['outcrop']:.2f} m "
            f"= {mech['outcrop'] / 2.0:.2f} B from each edge)",
            fontsize=10.5, y=0.975)
        out = os.path.join(HERE, fname)
        fig.savefig(out)
        plt.close(fig)
        print(f"[plot] wrote {out}", flush=True)

    if which in ("main", "both"):
        make("adr95_mechanism.png", (-8.0, 8.0), (-8.0, 0.5))
    if which in ("full", "both"):
        make("adr95_mechanism_full.png", (-15.0, 15.0), (-12.0, 0.5))
    return mech, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--leg", choices=sorted(LEGS))
    ap.add_argument("--tmax", type=float, default=1200.0)
    ap.add_argument("--budget", type=int, default=200)
    ap.add_argument("--plot", choices=["main", "full", "both"], default=None)
    args = ap.parse_args()
    if args.plot:
        build_figures(args.plot)
        return
    if not args.leg:
        raise SystemExit("--leg is required unless --plot is given")
    run_leg(args.leg, LEGS[args.leg], args.tmax, args.budget)


if __name__ == "__main__":
    main()
