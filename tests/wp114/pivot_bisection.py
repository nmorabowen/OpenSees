"""WP-114 ask 1 -- bisect the Pardiso perturbed-pivot count on a BezierTri6
`-bbar` field (TIMs Workbench fork_request_bezier_2026-09-18.md).

We do NOT have the reporter's model (26 970-element band-refined strip-footing
mesh). This script builds small, structured, own-authored reproducers and
bisects across the axes named in the ask:

  (a) -bbar on/off
  (b) LadrunoKinematicCoupling tie (Kt) onto a footing skin, vs direct `sp`
      displacement control on the same skin nodes, vs no footing at all
  (c) partial footing (free top surface either side) vs full-width footing
      (no free surface)
  (d) LadrunoQuad -bbar on the equivalent corner-node mesh, as a control
  plus: Kt sweep, Poisson ratio sweep, PARDISO -matrixType sweep, and a
  mesh-size sweep (to see whether the count scales with element count or
  with some fixed local DOF set).

For the smallest case that reproduces a nonzero count, it also pulls the
dense elastic K (system FullGeneral, `ops.printA("-ret")`), computes its
eigen-decomposition with numpy, and reports which nodes/DOFs the smallest-
eigenvalue modes live on (mapped back via `ops.nodeDOFs`).

Run: python3.12 tests/wp114/pivot_bisection.py
Output: tests/wp114/out/pivot_bisection_report.txt (also printed to stdout).
"""
import os
import re
import sys
import tempfile

WORKTREE = r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\fork-request-scoping-dee6db"
sys.path.insert(0, os.path.join(WORKTREE, "dist", "bin"))

import opensees as ops  # noqa: E402
import numpy as np  # noqa: E402

OUT_DIR = os.path.join(WORKTREE, "tests", "wp114", "out")
os.makedirs(OUT_DIR, exist_ok=True)

print("ladrunoBuild:", ops.ladrunoBuild())


# ============================================================ fd capture ===
def capture(fn):
    """Redirect the C-level fd 1/2 (where opserr actually writes) to a temp
    file for the duration of fn(), then return the captured text. This is a
    STANDALONE script (no pytest fd-capture already active), so a plain
    os.dup2 pair works -- the trap recorded in test_ladruno_sanisand_implex.py
    (nested redirects racing pytest's OWN capfd machinery) does not apply
    here; there is nothing else holding fd 1/2."""
    stdout_fd, stderr_fd = 1, 2
    saved_out = os.dup(stdout_fd)
    saved_err = os.dup(stderr_fd)
    tmp = tempfile.TemporaryFile(mode="w+b")
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(tmp.fileno(), stdout_fd)
    os.dup2(tmp.fileno(), stderr_fd)
    try:
        fn()
    finally:
        # CRITICAL: flush Python's OWN buffered stdout/stderr BEFORE restoring
        # the fds. Without this, a Python print() inside fn() can sit in
        # Python's internal buffer and only hit the (by-then-restored)
        # ORIGINAL fd later -- so it silently escapes capture and lands on
        # the real terminal instead of in `tmp`. C++ opserr writes are
        # unbuffered/auto-flushed and are unaffected either way, but any
        # Python-side print() (e.g. an "RC=%d" marker) needs this.
        sys.stdout.flush()
        sys.stderr.flush()
        os.dup2(saved_out, stdout_fd)
        os.dup2(saved_err, stderr_fd)
        os.close(saved_out)
        os.close(saved_err)
    tmp.seek(0)
    data = tmp.read().decode("utf-8", errors="replace")
    tmp.close()
    return data


PIVOT_RE = re.compile(r"perturbed pivots = (\d+)")
BLOCK_RE = re.compile(r"PARDISO stats:")


def parse_pivots(text):
    blocks = len(BLOCK_RE.findall(text))
    counts = [int(m) for m in PIVOT_RE.findall(text)]
    return counts, blocks


# ================================================================ mesh =====
def _graded_coords(n, total, ratio):
    """n+1 coordinates over [0, total], n intervals in geometric progression
    with the given ratio (interval i = d0 * ratio**i). ratio < 1 bunches
    points toward the HIGH end (finer spacing there); ratio == 1 is uniform.
    Used to emulate a "band-refined" mesh: finer triangles near the footing
    (top) than near the far-field base."""
    if abs(ratio - 1.0) < 1e-12:
        step = total / n
        return [i * step for i in range(n + 1)]
    d0 = total * (1.0 - ratio) / (1.0 - ratio ** n)
    xs = [0.0]
    d = d0
    for _ in range(n):
        xs.append(xs[-1] + d)
        d *= ratio
    xs[-1] = total  # kill fp drift
    return xs


def build_tri6_mesh(nx, ny, L, H, ry=1.0, rx=1.0, alt_diag=False):
    """nx*ny cell grid (optionally GRADED: `ry`/`rx` < 1 bunch cells toward
    the top / toward x=L, respectively -- emulating a band-refined field
    finer near the footing than the far field), each cell split into 2
    BezierTri6-ordered triangles. With `alt_diag`, the split direction
    alternates in a checkerboard pattern instead of always running
    bottom-left -> top-right (more varied triangle shapes, closer to an
    unstructured band-refined mesh than a single fixed diagonal). Returns:
      nodes: {id: (x, y)}
      elements: [(n1..n6), ...]   (corner, corner, corner, mid12, mid23, mid31)
      corner: {(i, j): id}        i in 0..nx, j in 0..ny
      row_nodes(j): [(x, id), ...] sorted by x -- corners + horizontal midsides
                     at grid row j (used for the base and the top skin)
    """
    xs = _graded_coords(nx, L, rx)
    ys = _graded_coords(ny, H, ry)
    nodes = {}
    corner = {}
    next_id = [1]

    def new_node(x, y):
        nid = next_id[0]
        next_id[0] += 1
        nodes[nid] = (x, y)
        return nid

    for j in range(ny + 1):
        for i in range(nx + 1):
            corner[(i, j)] = new_node(xs[i], ys[j])

    edge_mid = {}

    def mid(a, b):
        key = (a, b) if a < b else (b, a)
        if key not in edge_mid:
            xa, ya = nodes[a]
            xb, yb = nodes[b]
            edge_mid[key] = new_node(0.5 * (xa + xb), 0.5 * (ya + yb))
        return edge_mid[key]

    elements = []
    for j in range(ny):
        for i in range(nx):
            c00 = corner[(i, j)]
            c10 = corner[(i + 1, j)]
            c11 = corner[(i + 1, j + 1)]
            c01 = corner[(i, j + 1)]
            flip = alt_diag and ((i + j) % 2 == 1)
            if not flip:
                # diagonal c00-c11
                elements.append((c00, c10, c11, mid(c00, c10), mid(c10, c11), mid(c11, c00)))
                elements.append((c00, c11, c01, mid(c00, c11), mid(c11, c01), mid(c01, c00)))
            else:
                # diagonal c10-c01
                elements.append((c00, c10, c01, mid(c00, c10), mid(c10, c01), mid(c01, c00)))
                elements.append((c10, c11, c01, mid(c10, c11), mid(c11, c01), mid(c01, c10)))

    def row_nodes(j):
        row = [(nodes[corner[(i, j)]][0], corner[(i, j)]) for i in range(nx + 1)]
        for i in range(nx):
            a, b = corner[(i, j)], corner[(i + 1, j)]
            key = (a, b) if a < b else (b, a)
            if key in edge_mid:
                mnid = edge_mid[key]
                row.append((nodes[mnid][0], mnid))
        row.sort()
        return row

    return nodes, elements, corner, row_nodes


# ============================================================== builder ====
CASE_DEFAULTS = dict(
    nx=6, ny=4, L=6.0, H=4.0,
    ry=1.0, rx=1.0, alt_diag=False,   # mesh grading / diagonal pattern
    E=2.0e4, nu=0.30,
    bbar=True,
    element="tri6",          # "tri6" | "quad"
    footing_frac=0.34,       # fraction of top width covered by the footing
    tie_mode="kinematic",    # "kinematic" | "sp" | "none"
    Kt=5.0e9,
    disp=-0.01,              # prescribed footing settlement (m), tie_mode != "none"
    matrixType=0,            # PARDISO -matrixType (0 unsym/mtype11, 1 spd/mtype2, 2 symgen/mtype-2)
)


def build_model(**opts):
    p = dict(CASE_DEFAULTS)
    p.update(opts)

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    nodes, elements, corner, row_nodes_all = build_tri6_mesh(
        p["nx"], p["ny"], p["L"], p["H"], ry=p["ry"], rx=p["rx"], alt_diag=p["alt_diag"])
    corner_ids = set(corner.values())

    # CRITICAL for the "quad" control case: LadrunoQuad only ever references
    # the 4 corner nodes per cell. build_tri6_mesh() unconditionally also
    # allocates midside-node IDs/coordinates for the tri6 mesh; if those get
    # registered with ops.node() but never attached to ANY element, they are
    # floating zero-stiffness DOFs -- a guaranteed singular/near-singular K
    # that has nothing to do with the -bbar hypothesis being tested. Only add
    # the nodes each element type actually uses.
    if p["element"] == "quad":
        used_ids = corner_ids
    else:
        used_ids = set(nodes.keys())

    for nid in sorted(used_ids):
        x, y = nodes[nid]
        ops.node(nid, x, y)

    def row_nodes(j):
        entries = row_nodes_all(j)
        if p["element"] == "quad":
            entries = [(x, nid) for x, nid in entries if nid in corner_ids]
        return entries

    ops.nDMaterial("ElasticIsotropic", 1, p["E"], p["nu"])

    if p["element"] == "tri6":
        for k, e in enumerate(elements, start=1):
            args = ["BezierTri6", k, *e, 1.0, "PlaneStrain", 1]
            if p["bbar"]:
                args.append("-bbar")
            ops.element(*args)
    elif p["element"] == "quad":
        # equivalent field: same corner grid, one LadrunoQuad per cell, corners only
        k = 1
        for j in range(p["ny"]):
            for i in range(p["nx"]):
                n1, n2 = corner[(i, j)], corner[(i + 1, j)]
                n3, n4 = corner[(i + 1, j + 1)], corner[(i, j + 1)]
                form = "bbar" if p["bbar"] else "std"
                ops.element("LadrunoQuad", k, n1, n2, n3, n4, 1,
                            "-formulation", form, "-type", "PlaneStrain")
                k += 1
    else:
        raise ValueError(p["element"])

    # base: fully fixed
    for x, nid in row_nodes(0):
        ops.fix(nid, 1, 1)

    # top skin under the footing window
    top = row_nodes(p["ny"])
    x0, x1 = 0.5 * p["L"] * (1 - p["footing_frac"]), 0.5 * p["L"] * (1 + p["footing_frac"])
    footing_nodes = [nid for x, nid in top if x0 - 1e-9 <= x <= x1 + 1e-9]

    ref_tag = None
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)

    if p["tie_mode"] == "none" or p["footing_frac"] <= 0.0:
        pass  # free surface everywhere on top, no footing at all
    elif p["tie_mode"] == "sp":
        for nid in footing_nodes:
            ops.sp(nid, 1, 0.0)
            ops.sp(nid, 2, p["disp"])
    elif p["tie_mode"] == "kinematic":
        ref_tag = 999999
        xc = 0.5 * p["L"]
        ops.node(ref_tag, xc, p["H"], "-ndf", 3)
        ops.fix(ref_tag, 1, 0, 1)   # ux locked, theta locked -- pure vertical settlement
        ops.sp(ref_tag, 2, p["disp"])
        ops.element("LadrunoKinematicCoupling", 888888, ref_tag,
                    len(footing_nodes), *footing_nodes, "-k", p["Kt"])
    else:
        raise ValueError(p["tie_mode"])

    # NOTE: constraints/numberer/system/algorithm/integrator/analysis are
    # deliberately NOT set here. `ops.analysis("Static")` freezes the current
    # (possibly still-default) system at CALL time -- calling it before
    # `ops.system(...)` is what produced the "no LinearSOE specified,
    # ProfileSPDLinSOE default will be used" warning on every case in the
    # first pass, silently running every case through ProfileSPD instead of
    # Pardiso. The caller sets every analysis component, in order, and calls
    # `analysis("Static")` LAST, immediately before analyze().

    meta = dict(n_nodes=len(used_ids), n_elements=len(elements) if p["element"] == "tri6"
                else p["nx"] * p["ny"],
                n_footing=len(footing_nodes), footing_nodes=footing_nodes,
                ref_tag=ref_tag, nodes={k: v for k, v in nodes.items() if k in used_ids},
                corner=corner, row_nodes=row_nodes, params=p)
    return meta


def _setup_analysis(system_args):
    # `Transformation` (not `Plain`): the reference node in the kinematic-tie
    # cases carries a NON-homogeneous SP (a prescribed nonzero settlement) on
    # a mixed-ndf (ndf=3) node. `Plain` choked on that combination --
    # "WARNING PlainHandler::handle() - non-homogeneos constraint ... for
    # node 999999" -- and silently zeroed the prescribed displacement instead
    # of refusing outright, which would have made the dense-K probe solve a
    # DIFFERENT (unloaded) problem without saying so. Transformation handles
    # non-homogeneous SPs and mixed ndf correctly (also matches the pattern
    # used by tests/test_bezierTri6_element.py's prescribed-BC decks).
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


# ======================================================== pivot measurement
def measure_pivots(**opts):
    """Fresh model, `system Pardiso -stats`, one analyze(1), parse the
    perturbed-pivot count out of the captured PARDISO stats block."""
    meta = build_model(**opts)
    mtype = meta["params"]["matrixType"]

    if mtype == 1:
        sys_args = ["Pardiso", "-spd", "-stats"]
    elif mtype == 2:
        sys_args = ["Pardiso", "-symmetric", "-stats"]
    else:
        sys_args = ["Pardiso", "-matrixType", 0, "-stats"]

    def run():
        _setup_analysis(sys_args)
        rc = ops.analyze(1)
        print("RC=%d" % rc)

    text = capture(run)
    counts, blocks = parse_pivots(text)
    rc_ok = "RC=0" in text
    return dict(counts=counts, blocks=blocks, rc_ok=rc_ok, meta=meta, raw=text)


# ============================================== dense-K nullspace analysis
def dense_K_modes(n_show=6, **opts):
    meta = build_model(**opts)

    def run():
        _setup_analysis(["FullGeneral"])
        rc = ops.analyze(1)
        print("RC=%d" % rc)

    text = capture(run)
    if "RC=0" not in text:
        return None, None, meta, text

    n = ops.systemSize()
    flat = ops.printA("-ret")
    K = np.array(flat, dtype=float).reshape(n, n, order="F")
    # symmetrize defensively (elastic tangent should already be symmetric to fp noise)
    Ksym = 0.5 * (K + K.T)
    asym = np.max(np.abs(K - K.T)) if n > 0 else 0.0

    w, v = np.linalg.eigh(Ksym)
    idx = np.argsort(w)
    w = w[idx]
    v = v[:, idx]

    # map eqn index -> (node, local dof 1/2)
    eqn_to_dof = {}
    for nid in meta["nodes"]:
        dofs = ops.nodeDOFs(nid)
        for local, eq in enumerate(dofs, start=1):
            if eq >= 0:
                eqn_to_dof[eq] = (nid, local)
    if meta["ref_tag"] is not None:
        for local, eq in enumerate(ops.nodeDOFs(meta["ref_tag"]), start=1):
            if eq >= 0:
                eqn_to_dof[eq] = (meta["ref_tag"], local)

    modes = []
    for m in range(min(n_show, n)):
        vec = v[:, m]
        order = np.argsort(-np.abs(vec))[:8]
        contrib = [(eqn_to_dof.get(int(e), ("?", "?")), float(vec[e])) for e in order]
        modes.append((float(w[m]), contrib))

    return w, modes, meta, dict(asym=asym, n=n)


# ==================================================================== main
def fmt_counts(res):
    if not res["counts"]:
        return "0 (no PARDISO stats block matched)" if res["rc_ok"] else "n/a (analyze FAILED)"
    return "/".join(str(c) for c in res["counts"]) + (
        f"  [{res['blocks']} block(s)]" if res["blocks"] != len(res["counts"]) else "")


def main():
    report = []

    def log(line=""):
        print(line)
        report.append(line)

    log("=" * 78)
    log("WP-114 ask 1: Pardiso perturbed-pivot bisection on BezierTri6 -bbar")
    log("ladrunoBuild: %s" % ops.ladrunoBuild())
    log("=" * 78)

    # ---- table 1: the main bisection axes on a fixed small mesh -----------
    base = dict(nx=6, ny=4, L=6.0, H=4.0, footing_frac=0.34, Kt=5.0e9, nu=0.30)

    cases = [
        ("A tri6 -bbar  + kinematic tie (Kt=5e9)", dict(base, element="tri6", bbar=True, tie_mode="kinematic")),
        ("B tri6 no-bbar + kinematic tie (Kt=5e9)", dict(base, element="tri6", bbar=False, tie_mode="kinematic")),
        ("C tri6 -bbar  + sp-driven footing (no tie element)", dict(base, element="tri6", bbar=True, tie_mode="sp")),
        ("D tri6 -bbar  + full-width footing, no free surface", dict(base, element="tri6", bbar=True, tie_mode="kinematic", footing_frac=1.0)),
        ("E tri6 -bbar  + NO footing (all-free top, no load)", dict(base, element="tri6", bbar=True, tie_mode="none", footing_frac=0.0)),
        ("F quad -bbar  + kinematic tie (control)", dict(base, element="quad", bbar=True, tie_mode="kinematic")),
        ("G quad std    + kinematic tie (control)", dict(base, element="quad", bbar=False, tie_mode="kinematic")),
    ]

    log("\n---- table 1: bisection axes (nx=%d ny=%d, footing_frac default %.2f) ----"
        % (base["nx"], base["ny"], base["footing_frac"]))
    log(f"{'case':55s} {'nodes':>6s} {'elems':>6s} {'pivots':>18s} {'rc':>4s}")
    results = {}
    for name, opts in cases:
        res = measure_pivots(**opts)
        results[name] = res
        ok = "OK" if res["rc_ok"] else "FAIL"
        log(f"{name:55s} {res['meta']['n_nodes']:6d} {res['meta']['n_elements']:6d} "
            f"{fmt_counts(res):>18s} {ok:>4s}")

    # ---- table 2: Kt sweep on case A ---------------------------------------
    log("\n---- table 2: penalty Kt sweep (tri6 -bbar, kinematic tie) ----")
    log(f"{'Kt':>12s} {'pivots':>18s}")
    for Kt in (1.0e6, 1.0e9, 5.0e9, 1.0e12):
        res = measure_pivots(**dict(base, element="tri6", bbar=True, tie_mode="kinematic", Kt=Kt))
        log(f"{Kt:12.1e} {fmt_counts(res):>18s}")

    # ---- table 3: Poisson ratio sweep --------------------------------------
    log("\n---- table 3: Poisson ratio sweep (tri6 -bbar, kinematic tie, Kt=5e9) ----")
    log(f"{'nu':>8s} {'pivots':>18s}")
    for nu in (0.30, 0.45, 0.49, 0.4999):
        res = measure_pivots(**dict(base, element="tri6", bbar=True, tie_mode="kinematic", nu=nu))
        log(f"{nu:8.4f} {fmt_counts(res):>18s}")

    # ---- table 4: PARDISO -matrixType sweep --------------------------------
    log("\n---- table 4: PARDISO matrixType sweep (tri6 -bbar, kinematic tie) ----")
    log(f"{'matrixType':>12s} {'mtype':>8s} {'pivots':>18s}")
    mtype_label = {0: "11 (unsym)", 1: "2 (spd)", 2: "-2 (symgen)"}
    for mt in (0, 1, 2):
        res = measure_pivots(**dict(base, element="tri6", bbar=True, tie_mode="kinematic", matrixType=mt))
        log(f"{mt:12d} {mtype_label[mt]:>8s} {fmt_counts(res):>18s}")

    # ---- table 5: mesh-size scaling -----------------------------------------
    log("\n---- table 5: mesh-size scaling (tri6 -bbar, kinematic tie, Kt=5e9) ----")
    log(f"{'nx':>4s} {'ny':>4s} {'elems':>7s} {'footing_n':>10s} {'pivots':>18s}")
    scale_rows = []
    for nx, ny in ((4, 3), (6, 4), (8, 6), (10, 8), (14, 10)):
        opts = dict(base, element="tri6", bbar=True, tie_mode="kinematic", nx=nx, ny=ny)
        res = measure_pivots(**opts)
        scale_rows.append((nx, ny, res))
        log(f"{nx:4d} {ny:4d} {res['meta']['n_elements']:7d} {res['meta']['n_footing']:10d} "
            f"{fmt_counts(res):>18s}")

    # ---- table 5b: much larger uniform mesh, toward the user's scale ------
    log("\n---- table 5b: larger uniform meshes (tri6 -bbar, kinematic tie, Kt=5e9) ----")
    log(f"{'nx':>4s} {'ny':>4s} {'elems':>7s} {'pivots':>18s}")
    for nx, ny in ((30, 20), (60, 40)):
        opts = dict(base, element="tri6", bbar=True, tie_mode="kinematic", nx=nx, ny=ny)
        res = measure_pivots(**opts)
        log(f"{nx:4d} {ny:4d} {res['meta']['n_elements']:7d} {fmt_counts(res):>18s}")

    # ---- table 6: graded ("band-refined") + criss-cross diagonal mesh -----
    # The user's field is described as "band-refined": finer near the footing
    # than the far field, and (unlike our clean fixed-diagonal split) almost
    # certainly NOT a single repeated triangle shape. `ry`/`rx` < 1 bunch the
    # cells toward the footing/top; `alt_diag` alternates the split direction
    # checkerboard-style for more varied triangle geometry.
    log("\n---- table 6: graded mesh + alternating diagonal (tri6 -bbar, kinematic tie) ----")
    log(f"{'nx':>4s} {'ny':>4s} {'ry':>6s} {'rx':>6s} {'alt':>4s} {'elems':>7s} {'pivots':>18s}")
    graded_cases = [
        (10, 8, 0.75, 1.0, False),
        (10, 8, 0.55, 1.0, False),
        (10, 8, 0.75, 1.0, True),
        (16, 12, 0.6, 0.85, True),
        (30, 20, 0.5, 0.8, True),
    ]
    for nx, ny, ry, rx, alt in graded_cases:
        opts = dict(base, element="tri6", bbar=True, tie_mode="kinematic",
                     nx=nx, ny=ny, ry=ry, rx=rx, alt_diag=alt)
        res = measure_pivots(**opts)
        log(f"{nx:4d} {ny:4d} {ry:6.2f} {rx:6.2f} {str(alt):>4s} {res['meta']['n_elements']:7d} "
            f"{fmt_counts(res):>18s}")

    # ---- dense-K nullspace / near-null modes on the smallest reproducer ---
    log("\n---- dense K: smallest-eigenvalue modes on case A's mesh ----")
    w, modes, meta, info = dense_K_modes(n_show=8, **dict(base, element="tri6", bbar=True, tie_mode="kinematic"))
    if w is None:
        log("dense-K analyze() FAILED for case A -- skipping mode dump: %s" % info)
    else:
        log(f"n={info['n']}  max|K-K^T|={info['asym']:.3e}")
        log("lowest eigenvalues: " + ", ".join(f"{x:.4g}" for x in w[:8]))
        for i, (lam, contrib) in enumerate(modes):
            log(f"  mode {i} (lambda={lam:.4g}):")
            for (nid, local), val in contrib:
                x, y = meta["nodes"].get(nid, (float("nan"), float("nan"))) if isinstance(nid, int) else (float("nan"), float("nan"))
                is_corner = isinstance(nid, int) and any(nid == v for v in meta["corner"].values())
                kind = "corner" if is_corner else ("ref" if nid == meta["ref_tag"] else "midside")
                log(f"    node {nid!s:>8} dof{local} ({kind:7s} x={x:7.3f} y={y:7.3f})  comp={val:+.4f}")

    # ---- no-footing free mesh: same eigen-probe, to see if -bbar alone
    #      (no tie, no load) already creates a near-null mode --------------
    log("\n---- dense K: smallest-eigenvalue modes, NO footing (case E's mesh) ----")
    w2, modes2, meta2, info2 = dense_K_modes(n_show=6, **dict(base, element="tri6", bbar=True,
                                                                tie_mode="none", footing_frac=0.0))
    if w2 is None:
        log("dense-K analyze() FAILED for the no-footing mesh: %s" % info2)
    else:
        log(f"n={info2['n']}  max|K-K^T|={info2['asym']:.3e}")
        log("lowest eigenvalues: " + ", ".join(f"{x:.4g}" for x in w2[:6]))

    log("\n" + "=" * 78)
    log("done.")

    with open(os.path.join(OUT_DIR, "pivot_bisection_report.txt"), "w") as f:
        f.write("\n".join(report) + "\n")


if __name__ == "__main__":
    main()
