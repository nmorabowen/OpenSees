"""LadrunoUP (ELE 33017) — ADR-71 P4 INIT-RECIPE + RECORDER-CHANNEL battery (WP4.D).

The honest-p contract (ADR-71 §3.2) makes the nodal pressure DOF *is* p a
first-class, directly recorded / constrained / initialized quantity. This
battery gates the three initialization *recipes* (there is no bulk `ic` command
for p), the documented displacement-zeroing sequencing foot-gun, and the fork's
own PressureSource / Monitor recorder channels end-to-end.

Gates (one test each per the WP4.D pin i–iv):

  (i)  THREE init routes on a 2x10 Q4 equal-order column, each asserting the
       post-init p field AND that a transient continuation starts WITHOUT a
       pressure shock (bounded first-step |Δp|):
         (a) static steady-seepage stage (drained top p=0, prescribed base head
             via Penalty sp) -> linear steady p profile -> switch to transient
             (wipeAnalysis, keep domain) -> first-step Δp bounded;
         (b) gravity-stage transient (γ=0.6, big Δt, -dynSeepage off, fluid
             body force) -> hydrostatic p = ρ_f·g·depth within ±1% at every
             carrier -> continuation without shock;
         (c) scripted ops.setNodeDisp(node, pDof, value, '-commit') per §3.2
             route (c) -> element GP porePressure sees the committed p (Np·p) ->
             transient continuation without shock.

  (ii) SEQUENCING-TRAP regression (documentation gate, §3.2): under the honest-p
       contract a displacement-zeroing step (ops.reset() == Domain::revertToStart,
       the InitialStateAnalysis-off / revertToStart family — bare
       `InitialStateAnalysis` is NOT a Python command in this build) ALSO zeroes
       nodal p. Demonstrate the WRONG order (establish p, then zero -> p gone)
       and the CORRECT order (zero first, then establish p -> p survives).

  (iii) PressureSource end-to-end (⟨FW-F4⟩): a `recorder ladruno -N pressure` on
       a MIXED model carrying BOTH a LadrunoUP element and an upstream quadUP
       element (disjoint columns). Read the .ladruno HDF5 back and assert the
       recorded pressure == nodeDisp p-slot for LadrunoUP nodes (honest-p) and
       == nodeVel p-slot for the quadUP nodes (legacy rate-DOF path preserved).

  (iv) Monitor channel smoke (LadrunoMonitorRecorder): tail p on the LadrunoUP
       column via `recorder Monitor -resp disp -dof <ndm+1>` (honest-p: p is the
       disp slot). One write + SWMR-HDF5 read-back smoke.

Run:  py -3.12 -m pytest tests/test_ladruno_up_init_recorders.py -x -q
(the worktree dist/bin opensees.pyd is CPython 3.12 — plain `python` is 3.11
and will not import it).
"""
import glob
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

# h5py under HDF5 file-locking can refuse to open the recorder's own file on
# Windows/network paths (fork memory) — disable locking before any h5py import.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

pytestmark = [pytest.mark.zone_b]

_GENERAL = ("UmfPack",)   # the honest-p tangent is unsymmetric (ADR §3.2 <FW-F1>)

# canonical 2x10 Q4 column parameters (the WP4.D pinned rig)
_COL = dict(nely=10, L=10.0, E=1.0e4, nu=0.2, rho=2.0,
            Kf=4444.0, poro=0.4, rhoF=1.0, kbar=1.0e-4)


# --------------------------------------------------------------------------
# rig builders
# --------------------------------------------------------------------------
def _make_column(extra=(), p=_COL):
    """2-node-wide Q4 LadrunoUP column (ndf=3), no BCs/loads applied yet.
    Node grid nid[(i, j)] with i in {0,1} (width) and j in 0..nely (height)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    k = 1
    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.node(k, float(i), p["L"] * j / p["nely"])
            nid[(i, j)] = k
            k += 1
    ops.nDMaterial("ElasticIsotropic", 1, p["E"], p["nu"], p["rho"])
    for j in range(p["nely"]):
        ops.element("LadrunoUP", j + 1,
                    nid[(0, j)], nid[(1, j)], nid[(1, j + 1)], nid[(0, j + 1)],
                    1, "-Kf", p["Kf"], "-poro", p["poro"], "-rhoF", p["rhoF"],
                    "-perm", p["kbar"], p["kbar"], *extra)
    return nid


def _pin_u(nid, p=_COL):
    """Restrain ux everywhere + uy at the base (standard 1D-column kinematics)."""
    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, (1 if j == 0 else 0), 0)


def _p_profile(nid, p=_COL):
    """Committed nodal pore pressures up the left edge (j = 0..nely)."""
    return np.array([ops.nodeDisp(nid[(0, j)], 3) for j in range(p["nely"] + 1)])


def _cv(p=_COL):
    Eoed = p["E"] * (1 - p["nu"]) / ((1 + p["nu"]) * (1 - 2 * p["nu"]))
    return p["kbar"] / (1.0 / Eoed + p["poro"] / p["Kf"])


# ==========================================================================
# (i-a) static steady-seepage stage -> transient continuation
# ==========================================================================
@pytest.mark.t1
def test_init_route_a_static_steady_seepage_then_transient():
    """Route (a): establish p with a STATIC steady-seepage solve (a genuinely
    new capability — upstream parks H in getDamp and a Static solve is
    structurally singular), then continue transient with NO pressure shock.

    Rig: pure seepage (all u fixed), drained top p=0, base head dP via a
    Penalty sp (§3.2 quirk: staged prescribed-head uses Penalty, NOT
    Transformation). Steady p is exactly linear on the P1 mesh:
        p(y) = dP·(1 - y/L)   (dP at the base, 0 at the drained top).
    Continuation: wipeAnalysis (domain kept), Newmark γ=0.6 -dynSeepage off,
    one step — the committed state IS the steady restart, so |Δp| must be a
    small fraction of the profile scale dP (bound 1e-3·dP, generous vs the
    measured ~machine-level residual)."""
    p = _COL
    dP = 12.0
    nid = _make_column(extra=("-stab", "off"))
    # pure-seepage kinematics: fix ALL u
    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, 1, 0)
    for i in range(2):
        ops.fix(nid[(i, p["nely"])], 0, 0, 1)      # drained top p = 0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.sp(nid[(i, 0)], 3, dP)                  # base head p = dP (Penalty)

    ops.system(*_GENERAL)
    ops.numberer("Plain")
    ops.constraints("Penalty", 1e14, 1e14)
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    assert ops.analyze(1) == 0                      # a STATIC u-p solve succeeds

    p_static = _p_profile(nid)
    p_exact = dP * (1.0 - np.arange(p["nely"] + 1) / p["nely"])
    assert np.max(np.abs(p_static - p_exact)) < 1e-9 * dP, (
        f"static steady-seepage profile not linear: {p_static} vs {p_exact}")

    # --- transient continuation: wipe the analysis, KEEP the committed domain
    p_before = _p_profile(nid)
    ops.wipeAnalysis()
    ops.system(*_GENERAL)
    ops.numberer("Plain")
    ops.constraints("Penalty", 1e14, 1e14)
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(1, 0.5) == 0
    dp = np.max(np.abs(_p_profile(nid) - p_before))
    assert dp < 1e-3 * dP, (
        f"route (a) continuation shock: first-step |Δp| {dp:.3e} "
        f"(bound {1e-3 * dP:.3e}, profile scale dP={dP})")


# ==========================================================================
# (i-b) gravity-stage transient -> hydrostatic p -> continuation
# ==========================================================================
@pytest.mark.t1
def test_init_route_b_gravity_stage_hydrostatic():
    """Route (b): a gravity-stage transient (fluid body force always-on,
    -dynSeepage off, big Δt) drives the pore pressure to the hydrostatic
    no-flow steady state p(depth) = ρ_f·g·(L - y), exact on the P1 mesh
    (H·p = f_seep with drive = fluidBody => ∇p = ρ_f·fluidBody).

    Drained top p=0, sealed sides/base (impervious => no p fix). Gate: p at
    every carrier within ±1% of the profile scale ρ_f·g·L; then a continuation
    step with gravity held is shock-free (the field is the steady restart)."""
    p = _COL
    g = 9.81
    nid = _make_column(extra=("-fluidBody", 0.0, -g, "-stab", "off",
                              "-dynSeepage", "off"))
    _pin_u(nid)
    for i in range(2):
        ops.fix(nid[(i, p["nely"])], 0, 0, 1)      # drained top p = 0

    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.algorithm("Linear")
    ops.analysis("Transient")

    # big Δt so ~12 steps >> consolidation time L^2/cv (settles to steady)
    tconsol = p["L"] ** 2 / _cv()
    dt = tconsol / 2.0
    assert ops.analyze(14, dt) == 0

    depth = p["L"] - p["L"] * np.arange(p["nely"] + 1) / p["nely"]
    p_hydro = p["rhoF"] * g * depth
    p_fem = _p_profile(nid)
    scale = p["rhoF"] * g * p["L"]
    err = np.max(np.abs(p_fem - p_hydro))
    assert err < 0.01 * scale, (
        f"route (b) not hydrostatic within ±1%: max|Δp| {err:.4f} "
        f"(bound {0.01 * scale:.4f}); fem={p_fem}, hydro={p_hydro}")

    # continuation without shock: one more step with gravity held
    p_before = _p_profile(nid)
    assert ops.analyze(1, dt) == 0
    dp = np.max(np.abs(_p_profile(nid) - p_before))
    assert dp < 0.01 * scale, (
        f"route (b) continuation shock: |Δp| {dp:.4f} (bound {0.01 * scale:.4f})")


# ==========================================================================
# (i-c) scripted setNodeDisp p -commit -> element GP sees it -> continuation
# ==========================================================================
@pytest.mark.t1
def test_init_route_c_setnodedisp_commit():
    """Route (c): scripted per-node initialization via
    ops.setNodeDisp(node, pDof=3, value, '-commit') — the §3.2 route with no
    bulk `ic` command. A committed uniform p0 must be seen by the element GPs
    (porePressure = Np·p_nodal = p0 exactly, partition of unity) and by
    nodeDisp; a transient continuation on the sealed field is shock-free."""
    p = _COL
    p0 = 40.0
    nid = _make_column(extra=("-stab", "off", "-dynSeepage", "off"))
    # PURE-SEEPAGE rig: fix ALL u so the scripted uniform p0 is a genuine steady
    # restart (H·p0 = 0 for uniform p on a sealed Laplacian; u cannot move, so
    # there is no -Q·p0 coupling transient). Sealed p (impervious) — no p fix.
    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, 1, 0)

    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.setNodeDisp(nid[(i, j)], 3, p0, "-commit")

    # nodeDisp p-slot carries the committed value
    assert np.allclose(_p_profile(nid), p0, atol=1e-12)
    # element GP porePressure sees the committed field (Np·p = p0 on every GP)
    for e in range(1, p["nely"] + 1):
        gp_p = np.array(ops.eleResponse(e, "porePressure"))
        assert gp_p.size == 4, f"ele {e} porePressure len {gp_p.size} (want 4 GP)"
        assert np.allclose(gp_p, p0, atol=1e-9 * p0), (
            f"ele {e} GP porePressure {gp_p} != committed p0 {p0}")

    # transient continuation on the committed (sealed, zero-gradient) field:
    # H·p0 = 0 and ṗ = 0 => steady restart, first-step Δp bounded (no shock)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    assert ops.analyze(1, 0.1) == 0
    dp = np.max(np.abs(_p_profile(nid) - p0))
    assert dp < 1e-2 * p0, (
        f"route (c) continuation shock: first-step |Δp| {dp:.4f} "
        f"(bound {1e-2 * p0:.4f}, profile scale p0={p0})")


# ==========================================================================
# (ii) sequencing trap — displacement-zeroing ALSO zeroes p
# ==========================================================================
@pytest.mark.t1
def test_sequencing_trap_revert_to_start_zeroes_pressure():
    """Documentation gate (§3.2): under the honest-p contract, the
    displacement-zeroing staged trick (`InitialStateAnalysis off` ->
    Domain::revertToStart()) ALSO zeroes nodal p, because p is a displacement
    DOF here. The gate drives the zeroing via ops.reset() (OPS_resetModel ->
    Domain::revertToStart — the SAME revertToStart the InitialStateAnalysis-off
    branch calls, OpenSeesMiscCommands.cpp OPS_InitialStateAnalysis).

    NOTE (reported upward, element-source, NOT fixed here): the ADR-named form
    `ops.InitialStateAnalysis("off")` itself DOES zero p correctly, but then
    HEAP-CORRUPTS (0xc0000374) on the next model build — the crash is in
    addParameter(new InitialStateParameter(false)) forwarding an
    "InitialStateAnalysis" setParameter into the LadrunoUP element, NOT in
    revertToStart (ops.reset() runs the identical revertToStart cleanly). The
    gate therefore uses ops.reset(); minimal repro in the report.

    WRONG order: establish p, THEN zero -> p gone.
    CORRECT order (the guide rule): zero FIRST, then establish p -> survives."""
    p = _COL
    p0 = 40.0

    # WRONG order: establish p, THEN zero displacements -> p destroyed
    nid = _make_column(extra=("-stab", "off"))
    _pin_u(nid)
    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.setNodeDisp(nid[(i, j)], 3, p0, "-commit")
    assert np.allclose(_p_profile(nid), p0, atol=1e-12), "p not established"
    ops.reset()                                     # revertToStart-family zeroing
    p_after_zero = _p_profile(nid)
    assert np.max(np.abs(p_after_zero)) < 1e-9, (
        f"revertToStart should have zeroed p under the honest-p contract; "
        f"got {p_after_zero} (the documented foot-gun)")

    # CORRECT order: zero FIRST, then establish p -> survives
    nid = _make_column(extra=("-stab", "off"))
    _pin_u(nid)
    ops.reset()                                     # zero up front (no-op here)
    for j in range(p["nely"] + 1):
        for i in range(2):
            ops.setNodeDisp(nid[(i, j)], 3, p0, "-commit")
    assert np.allclose(_p_profile(nid), p0, atol=1e-12), (
        "correct-order p field did not survive")


# ==========================================================================
# (iii) PressureSource end-to-end — mixed LadrunoUP + upstream quadUP model
# ==========================================================================
def _read_ladruno_pressure(h5file):
    """Return {node_tag: last-step pressure} from a .ladruno HDF5 file.
    Layout: MODEL_STAGE[*]/RESULTS/ON_NODES/PRESSURE/{ID[nIds,1], DATA[T,nIds,ncomp]}."""
    import h5py
    ids = {"ID": None, "DATA": None}

    def visit(name, obj):
        up = name.upper().rstrip("/")
        if isinstance(obj, h5py.Dataset) and "PRESSURE" in up:
            if up.endswith("/ID") or up.endswith("PRESSURE/ID"):
                ids["ID"] = np.asarray(obj[...]).reshape(-1)
            elif up.endswith("/DATA") or up.endswith("PRESSURE/DATA"):
                ids["DATA"] = np.asarray(obj[...], dtype=float)

    with h5py.File(h5file, "r") as h:
        h.visititems(visit)
    assert ids["ID"] is not None and ids["DATA"] is not None, (
        "PRESSURE ID/DATA not found in the .ladruno output")
    data = ids["DATA"]                              # [T, nIds, ncomp]
    last = data[-1].reshape(len(ids["ID"]), -1)[:, 0]
    return {int(t): float(v) for t, v in zip(ids["ID"], last)}


@pytest.mark.t1
def test_pressuresource_mixed_ladruno_and_quadup(tmp_path):
    """⟨FW-F4⟩ end-to-end: a `recorder ladruno -N pressure` on a MIXED domain
    with BOTH a LadrunoUP element (honest-p: pressure = disp DOF 3) and an
    upstream quadUP element (rate-DOF: physical pressure = vel DOF 3), as two
    disjoint columns. The contract-aware PressureSource must record the
    disp-slot for LadrunoUP nodes and the vel-slot for quadUP nodes."""
    h5py = pytest.importorskip("h5py")
    out = str(tmp_path / "pressure.ladruno")

    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    # column A (LadrunoUP): nodes 1-4 ; column B (quadUP): nodes 5-8 (disjoint)
    ax = 0.0
    for tag, (x, y) in enumerate(
            [(ax, 0), (ax + 1, 0), (ax + 1, 1), (ax, 1),
             (ax + 3, 0), (ax + 4, 0), (ax + 4, 1), (ax + 3, 1)], 1):
        ops.node(tag, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    lad_nodes = [1, 2, 3, 4]
    qup_nodes = [5, 6, 7, 8]
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1, "-Kf", 4444.0, "-poro", 0.4,
                "-rhoF", 1.0, "-perm", 1e-4, 1e-4, "-stab", "off",
                "-dynSeepage", "off")
    # quadUP arg mapping (module donor _equiv.py): bulk = pre-combined Kf/n
    ops.element("quadUP", 2, 5, 6, 7, 8, 1.0, 1, 4444.0 / 0.4, 1.0,
                1e-4, 1e-4, 0.0, 0.0, 0.0)

    # kinematics + drained tops; load the tops so p develops in both columns
    for n in (1, 2, 5, 6):
        ops.fix(n, 1, 1, 0)                          # base: ux, uy fixed
    for n in (3, 4, 7, 8):
        ops.fix(n, 1, 0, 0)                          # top: ux fixed
    for n in (3, 4, 7, 8):
        ops.fix(n, 0, 0, 1)                          # drained top p = 0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in (3, 4, 7, 8):
        ops.load(n, 0.0, -2.5, 0.0)

    ops.recorder("ladruno", out, "-N", "pressure")
    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    for _ in range(5):
        assert ops.analyze(1, 0.05) == 0

    # snapshot the live node values BEFORE wipe (== the last recorded frame)
    disp_p = {n: ops.nodeDisp(n, 3) for n in lad_nodes}
    vel_p = {n: ops.nodeVel(n, 3) for n in qup_nodes}
    ops.wipe()                                       # destroy recorder -> flush/close

    files = glob.glob(out) or glob.glob(out.replace(".ladruno", "*.ladruno"))
    assert files, "no .ladruno output produced"
    rec = _read_ladruno_pressure(files[0])

    # sanity: a non-trivial pressure actually developed somewhere
    assert max(abs(v) for v in rec.values()) > 1e-6, "no pressure recorded"

    # LadrunoUP nodes: recorded == nodeDisp p-slot (honest-p, disp channel)
    for n in lad_nodes:
        assert n in rec, f"LadrunoUP node {n} missing from PRESSURE output"
        assert abs(rec[n] - disp_p[n]) < 1e-9 * max(abs(disp_p[n]), 1.0), (
            f"LadrunoUP node {n}: recorded {rec[n]} != nodeDisp p {disp_p[n]}")
    # quadUP nodes: recorded == nodeVel p-slot (legacy rate-DOF path preserved)
    for n in qup_nodes:
        assert n in rec, f"quadUP node {n} missing from PRESSURE output"
        assert abs(rec[n] - vel_p[n]) < 1e-9 * max(abs(vel_p[n]), 1.0), (
            f"quadUP node {n}: recorded {rec[n]} != nodeVel p {vel_p[n]}")

    # the two contracts genuinely diverge: at least one LadrunoUP node's
    # disp-p differs from its vel-p (else the test proves nothing)
    assert any(abs(ops_disp) > 1e-6 for ops_disp in disp_p.values()), (
        "LadrunoUP disp-slot pressures are all ~0 — contract dispatch untested")


# ==========================================================================
# (iv) Monitor channel smoke — tail p on the LadrunoUP column
# ==========================================================================
@pytest.mark.t1
def test_monitor_channel_pressure_smoke(tmp_path):
    """LadrunoMonitorRecorder smoke: tail the honest-p pressure via
    `recorder Monitor -resp disp -dof 3` (p is the disp slot). One write +
    SWMR-HDF5 read-back of the streamed frames. Cheap enough to gate at P4."""
    h5py = pytest.importorskip("h5py")
    p = _COL
    sink = str(tmp_path / "monitor_p.h5")

    nid = _make_column(extra=("-stab", "off", "-dynSeepage", "off"))
    _pin_u(nid)
    for i in range(2):
        ops.fix(nid[(i, p["nely"])], 0, 0, 1)      # drained top p = 0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(nid[(i, p["nely"])], 0.0, -5.0, 0.0)

    mid = nid[(0, p["nely"] // 2)]                  # mid-column carrier
    ops.recorder("Monitor", "-node", mid, "-dof", 3, "-resp", "disp",
                 "-sink", sink)
    ops.system(*_GENERAL)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.algorithm("Linear")
    ops.analysis("Transient")
    nstep = 8
    for _ in range(nstep):
        assert ops.analyze(1, 0.05) == 0
    p_live = ops.nodeDisp(mid, 3)
    ops.wipe()                                       # flush + close the SWMR sink

    assert os.path.isfile(sink) and os.path.getsize(sink) > 0, "no monitor sink"
    # read the streamed frames back (post-wipe the file is closed cleanly;
    # fall back to swmr mode if the library still flags it)
    def _open(f):
        try:
            return h5py.File(f, "r")
        except (OSError, IOError):
            return h5py.File(f, "r", libver="latest", swmr=True)
    with _open(sink) as h:
        frames = np.asarray(h["FRAMES"][...], dtype=float)
        times = np.asarray(h["TIME"][...], dtype=float)
    assert frames.shape[0] >= 1, f"monitor wrote no frames (shape {frames.shape})"
    assert times.shape[0] == frames.shape[0], "TIME/FRAMES length mismatch"
    # the last streamed p equals the live mid-column pressure
    assert abs(frames[-1].reshape(-1)[0] - p_live) < 1e-9 * max(abs(p_live), 1.0), (
        f"monitored p {frames[-1].reshape(-1)[0]} != live nodeDisp p {p_live}")
