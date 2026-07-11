"""LadrunoUP (ADR-71 P1) — EQUIVALENCE + STABILIZATION Zone-B battery.

Owns the P1 gate legs that compare the honest-p LadrunoUP element against the
upstream u-p elements at RESPONSE level, plus the B4 stabilization gate and the
bbar behavior gate (ADR 71 §7 P1 row; plan §P1 WP table):

  1. Two-leg FourNodeQuadUP equivalence (both with -dynSeepage off — upstream
     drops the dynamic-seepage force entirely, ADR §3.1):
       leg 1 (tight):      Newmark gamma=1/2, beta=1/4 — the unique pair where
                           the p-as-disp / p-as-vel parameterizations produce
                           the identical one-step rule. u AND p trajectories
                           <= 1e-6 rel (achieved ~1e-14, see test header).
       leg 2 (production): gamma=0.6, beta=0.3025 — Dt-halving MUTUAL
                           convergence (1e-6 at gamma=0.6 is mathematically
                           impossible: the two parameterizations differ by a
                           parasitic-root memory term damped only for
                           gamma > 1/2 — ADR §3.2 / plan §11).
  2. SSPquadUP pressure-block cross-check at matched stabilization alpha
     (different solid treatment — single-point stabilized vs full 2x2 —
     so tolerances are honest response-level, not operator-identity).
  3. B4 McGann checkerboard footing (ADR §7.1): CB index gates with an
     alpha0 sweep, plus the alpha = 6.8e-5 stabilization-value pin.
  4. bbar behavior gate: near-incompressible locking relief + FD consistency
     of the assembled tangent (residual vs printA) for std AND bbar.

ARG MAPPING (pinned here so the guide can cite it; verified against the
upstream parsers in this tree):

  FourNodeQuadUP (SRC/element/UP-ucsd/FourNodeQuadUP.cpp, OPS_FourNodeQuadUP):
      element quadUP  tag i j k l  thk matTag  bulk rhoF  permX permY <b1 b2 p>
    * bulk  = the PRE-COMBINED Qbar ~= Kf/n (code: oneOverKc = 1/kc in the
      storage block) -> pass Kf/poro for our (-Kf Kf -poro n) with alpha=1,
      Ks_grain -> inf (upstream hardwires alpha=1).
    * rhoF  = fluid mass density (drives the seepage body-force term only).
    * perm  = k_hydraulic/gammaW per axis — SAME convention as our -perm
      (their H multiplies grad-p directly; seepage balances rhoF*perm*b).
    * python dispatch name is "quadUP" (OpenSeesElementCommands.cpp:703).
    * upstream physical p = nodal VELOCITY of the p-DOF (rate-DOF trick);
      LadrunoUP p = nodal DISPLACEMENT of the p-DOF (honest-p contract).
      NO sign flip between the two (verified: both compression-positive).

  SSPquadUP (SRC/element/UWelements/SSPquadUP.cpp, OPS_SSPquadUP):
      element SSPquadUP tag i j k l matTag  thk fBulk fDen  k1 k2  e alpha <b1 b2>
    * fBulk = RAW fluid bulk Kf (code: oneOverQ = mPorosity/fBulk).
    * e     = void ratio -> n = e/(1+e); so pass e = n/(1-n).
    * alpha = the RAW stabilization parameter (their Kp = -4*alpha*J0*t*dN*dNp)
      == our '-stab <manual>' == alpha0*h^2/(Ks_skel + 4*Gs_skel/3) for the
      matched-'auto' comparison.
    * p records via nodeVel (same rate-DOF trick as quadUP).

CHECKERBOARD INDEX (test 3): two metrics over the interior nodes of the
structured (Ncell+1)^2 grid at t = 1 s (the ADR's check time):

  CB_lat = |<p, chi>| / (||p||*||chi||),  chi_ij = (-1)^(i+j)
           (the ADR §7.1 lattice-cosine; captures the pure alternating mode)
  CB_lap = ||p_ij - mean(4 neighbors)|| / ||p||
           (normalized interior-Laplacian roughness; captures ALL short-wave
           oscillation incl. the stripe mode (-1)^j that dominates this
           footing's instability)

CB_lap is the PRIMARY gate (captures the full oscillation, monotone across
the sweep); CB_lat is asserted only for the off-vs-auto separation — in a
pulse-load scaffolding variant it went NON-monotone in alpha0 (metric
fragility: heavier smoothing shrinks ||p|| faster than the lattice
projection), so it stays secondary. First-passing-run values (this build,
FullGeneral, Newmark 0.6/0.3025, dt=0.01, t=1s, sustained ramp):
    CB_lap: off=0.603, a0=0.05: 0.222, a0=0.25: 0.126, a0=0.5: 0.0946
    CB_lat: off=0.138, a0=0.05: 0.0531, a0=0.25: 0.0264, a0=0.5: 0.0199

Run mechanics (worktree): the module bootstraps <worktree>/dist/bin onto
PATH + DLL dirs + sys.path and asserts the loaded pyd is THIS worktree's —
the boot .pth preloads a stale pyd otherwise, so invoke with `py -3.12 -S`
and PYTHONPATH carrying pytest's site-packages (see project memory
"OpenSees Python test env").
"""
import math
import os
import sys

# --------------------------------------------------------------------------
# hermetic import: this worktree's build, not a stale preloaded pyd
# --------------------------------------------------------------------------
_WT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_BIN = os.path.join(_WT, "dist", "bin")
if os.path.isfile(os.path.join(_BIN, "opensees.pyd")) or os.path.isfile(
        os.path.join(_BIN, "opensees.so")):
    os.environ["PATH"] = _BIN + os.pathsep + os.environ.get("PATH", "")
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(_BIN)
    sys.path.insert(0, _BIN)
    import opensees as ops
    assert os.path.normcase(os.path.dirname(ops.__file__)) == \
        os.path.normcase(_BIN), (
        f"stale opensees import: {ops.__file__} (want {_BIN}); "
        "run with py -3.12 -S so the boot .pth cannot preload another build")
else:  # CI layout: opensees.so copied next to tests/
    from _testbed import ops  # noqa: F401

import numpy as np
import pytest

pytestmark = [pytest.mark.zone_a]


# ==========================================================================
# shared consolidation column (the equivalence rig)
# ==========================================================================
# Small transient consolidation column: 1 m x 4 m, 4 stacked Q4 u-p elements,
# laterally confined (u_x = 0), base fixed, drained top (p = 0), sealed
# elsewhere; vertical top load ramped over one step then held.
COL = dict(E=1.0e4, nu=0.3, rho=2.0,          # skeleton + saturated density
           poro=0.4, Kf=2.2e5, rhoF=1.0,      # fluid
           perm=1.0e-4,                        # k_hydraulic/gammaW
           H=4.0, b=1.0, nEy=4, load=-10.0)


def _column_nodes():
    ny = COL["nEy"] + 1
    nid = {}
    k = 1
    for j in range(ny):
        for i in range(2):
            ops.node(k, i * COL["b"], j * COL["H"] / COL["nEy"])
            nid[(i, j)] = k
            k += 1
    return nid


def _column_bcs(nid):
    ny = COL["nEy"] + 1
    for j in range(ny):
        for i in range(2):
            n = nid[(i, j)]
            if j == 0:
                ops.fix(n, 1, 1, 0)          # base: u fixed, sealed
            elif j == ny - 1:
                ops.fix(n, 1, 0, 1)          # top: drained (p = 0)
            else:
                ops.fix(n, 1, 0, 0)          # confined


def _column_load_and_analysis(nid, gamma, beta):
    ny = COL["nEy"] + 1
    # identical ramp for every dt: 0 -> 1 over [0, 0.05], then held
    ops.timeSeries("Path", 1, "-time", 0.0, 0.05, 1.0e6,
                   "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(nid[(i, ny - 1)], 0.0, COL["load"] / 2.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")                 # general solver: honest-p tangent
    # unbalance-based test (ADR P1 gate: upstream's int-p-dt DOF grows without
    # bound, so disp-increment tests mis-scale between the two contracts)
    ops.test("NormUnbalance", 1e-10, 30)
    ops.algorithm("Newton")
    ops.integrator("Newmark", gamma, beta)
    ops.analysis("Transient")


def _build_column(which, gamma, beta, stab=("off",), dyn_seepage="off"):
    """which in {'ladruno','quadup','ssp'}; returns the node-id dict."""
    c = COL
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = _column_nodes()
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    e = 1
    for j in range(c["nEy"]):
        n1, n2 = nid[(0, j)], nid[(1, j)]
        n3, n4 = nid[(1, j + 1)], nid[(0, j + 1)]
        if which == "ladruno":
            args = ["LadrunoUP", e, n1, n2, n3, n4, 1, "-thick", 1.0,
                    "-Kf", c["Kf"], "-poro", c["poro"], "-rhoF", c["rhoF"],
                    "-perm", c["perm"], c["perm"],
                    "-dynSeepage", dyn_seepage, "-stab"] + list(stab)
            ops.element(*args)
        elif which == "quadup":
            # ARG MAPPING: bulk = pre-combined Qbar ~= Kf/n (see module header)
            ops.element("quadUP", e, n1, n2, n3, n4, 1.0, 1,
                        c["Kf"] / c["poro"], c["rhoF"],
                        c["perm"], c["perm"], 0.0, 0.0, 0.0)
        elif which == "ssp":
            # ARG MAPPING: raw fBulk + void ratio e = n/(1-n) + raw stab alpha
            evoid = c["poro"] / (1.0 - c["poro"])
            assert stab[0] not in ("off", "auto"), \
                "ssp lane needs the resolved numeric alpha"
            ops.element("SSPquadUP", e, n1, n2, n3, n4, 1, 1.0,
                        c["Kf"], c["rhoF"], c["perm"], c["perm"],
                        evoid, float(stab[0]), 0.0, 0.0)
        else:
            raise ValueError(which)
        e += 1
    _column_bcs(nid)
    _column_load_and_analysis(nid, gamma, beta)
    return nid


def _run_column(which, gamma, beta, nsteps, dt, **kw):
    """Return (u_hist, p_hist): per-step arrays over every node.

    u = (u_x, u_y) of all nodes; p read per the CONTRACT CONVERSION —
    nodeDisp dof-3 for LadrunoUP, nodeVel dof-3 for the upstream rate-DOF
    elements (module header).
    """
    nid = _build_column(which, gamma, beta, **kw)
    keys = sorted(nid)
    u_hist, p_hist = [], []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, f"{which} transient step failed"
        u_row, p_row = [], []
        for kk in keys:
            n = nid[kk]
            d = ops.nodeDisp(n)
            u_row.extend(d[:2])
            p_row.append(d[2] if which == "ladruno" else ops.nodeVel(n, 3))
        u_hist.append(u_row)
        p_hist.append(p_row)
    return np.array(u_hist), np.array(p_hist)


def _rel_inf(a, b):
    """NaN-proof max-abs relative difference against the max-abs of b."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    assert np.all(np.isfinite(a)) and np.all(np.isfinite(b))
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(b)), 1e-30))


# ==========================================================================
# 1a. leg 1 (tight): gamma=1/2, beta=1/4 — trajectories <= 1e-6 rel
# ==========================================================================
def test_quadup_equivalence_leg1_tight():
    """LadrunoUP (-dynSeepage off, -stab off) vs FourNodeQuadUP under Newmark
    gamma=1/2 beta=1/4: full u AND p trajectories (all nodes, 40 steps)
    <= 1e-6 relative. First-passing-run achieved: u 3.6e-15, p 1.4e-14 —
    the two parameterizations produce the machine-identical one-step rule at
    this unique (gamma, beta), pinning the arg mapping AND the contract
    conversion (their p = nodeVel, ours = nodeDisp) in one shot."""
    nsteps, dt = 40, 0.05
    uL, pL = _run_column("ladruno", 0.5, 0.25, nsteps, dt)
    uU, pU = _run_column("quadup", 0.5, 0.25, nsteps, dt)
    rel_u = _rel_inf(uL, uU)
    rel_p = _rel_inf(pL, pU)
    # non-trivial response guard (a dead model would pass any rel gate)
    assert np.max(np.abs(uU)) > 1e-6 and np.max(np.abs(pU)) > 1e-2
    assert rel_u <= 1e-6, f"leg-1 u trajectories diverge: rel {rel_u:.3e}"
    assert rel_p <= 1e-6, f"leg-1 p trajectories diverge: rel {rel_p:.3e}"


# ==========================================================================
# 1b. leg 2 (production): gamma=0.6 — mutual Dt-halving convergence
# ==========================================================================
def test_quadup_equivalence_leg2_production_convergence():
    """gamma=0.6, beta=0.3025 (the ZS84 production set): the p-row memory term
    makes the two parameterizations genuinely different one-step rules, so the
    gate is MUTUAL Dt-convergence, not coincidence: the end-state difference
    |u_L - u_U| must shrink under every Dt halving (~first order) and the
    finest-Dt relative difference must be <= 1%.

    First-passing-run numbers (T=2 s, end-state, all nodes):
        dt=0.1:    diff 1.89e-07 (rel_u 1.37e-04, rel_p 2.08e-03)
        dt=0.05:   diff 5.86e-08 (rel_u 4.22e-05, rel_p 1.37e-04)  order 1.69
        dt=0.025:  diff 3.04e-08 (rel_u 2.19e-05, rel_p 2.27e-05)  order 0.95
        dt=0.0125: diff 1.52e-08 (rel_u 1.09e-05, rel_p 1.14e-05)  order 1.00
    (mean observed order 1.21 — first-order mutual convergence, as the
    parasitic-root memory-term analysis predicts.)
    """
    T = 2.0
    diffs_u, rels_u, rels_p = [], [], []
    for dt in (0.1, 0.05, 0.025, 0.0125):
        ns = int(round(T / dt))
        uL, pL = _run_column("ladruno", 0.6, 0.3025, ns, dt)
        uU, pU = _run_column("quadup", 0.6, 0.3025, ns, dt)
        # end-state comparison (trajectory-wide max is dominated by the ramp
        # step where the memory-term discrepancy peaks; the convergence claim
        # is about the solution the schemes settle on)
        diffs_u.append(float(np.max(np.abs(uL[-1] - uU[-1]))))
        rels_u.append(_rel_inf(uL[-1], uU[-1]))
        rels_p.append(_rel_inf(pL[-1], pU[-1]))
    # mutual convergence: every halving shrinks the disagreement
    for k in range(1, len(diffs_u)):
        assert diffs_u[k] < diffs_u[k - 1], (
            f"no mutual Dt-convergence: diffs {diffs_u}")
    # observed order over the sweep (geometric mean) ~ first order or better
    order = math.log(diffs_u[0] / diffs_u[-1], 2) / (len(diffs_u) - 1)
    assert order >= 0.5, f"mutual convergence order {order:.2f} < 0.5"
    assert rels_u[-1] <= 1e-2, f"finest-Dt u rel {rels_u[-1]:.3e} > 1%"
    assert rels_p[-1] <= 1e-2, f"finest-Dt p rel {rels_p[-1]:.3e} > 1%"


# ==========================================================================
# 2. SSPquadUP cross-check (stabilized single-point element)
# ==========================================================================
def _stab_alpha(h, E, nu, alpha0):
    Ks = E / (3.0 * (1.0 - 2.0 * nu))
    Gs = E / (2.0 * (1.0 + nu))
    return alpha0 * h * h / (Ks + 4.0 * Gs / 3.0)


def test_sspquadup_crosscheck():
    """One more equivalence datapoint against SSPquadUP at MATCHED
    stabilization (our '-stab auto 0.25' == their raw alpha computed from the
    same formula, h = element edge = 1 m here). SSP's solid is single-point
    stabilized (different operator), so this is response-level with honest
    tolerances: final settlement <= 1% rel, final base-pressure <= 1% rel,
    whole-trajectory u <= 10% rel (the early ramp steps differ by a few % —
    single-point vs 2x2 storage/coupling integration).

    First-passing-run: settlement rel 2.6e-5, final p rel 5.4e-5,
    trajectory max u rel 2.3e-2, p rel 5.9e-2.

    Also pins the SSP arg conversion: fBulk = raw Kf, e = n/(1-n)."""
    nsteps, dt = 60, 0.05
    h = COL["H"] / COL["nEy"]
    alpha = _stab_alpha(h, COL["E"], COL["nu"], 0.25)
    uL, pL = _run_column("ladruno", 0.6, 0.3025, nsteps, dt,
                         stab=("auto", 0.25))
    uS, pS = _run_column("ssp", 0.6, 0.3025, nsteps, dt, stab=(alpha,))
    assert np.max(np.abs(uS)) > 1e-6 and np.max(np.abs(pS)) > 1e-2
    rel_settle = _rel_inf(uL[-1], uS[-1])
    rel_p_end = _rel_inf(pL[-1], pS[-1])
    rel_u_traj = _rel_inf(uL, uS)
    rel_p_traj = _rel_inf(pL, pS)
    assert rel_settle <= 1e-2, f"final settlement rel {rel_settle:.3e}"
    assert rel_p_end <= 1e-2, f"final pressure rel {rel_p_end:.3e}"
    assert rel_u_traj <= 0.10, f"trajectory u rel {rel_u_traj:.3e}"
    assert rel_p_traj <= 0.15, f"trajectory p rel {rel_p_traj:.3e}"


# ==========================================================================
# 3. B4 McGann checkerboard footing — CB index + alpha0 sweep + alpha pin
# ==========================================================================
B4 = dict(E=25000.0, nu=0.3, rho=2.67, poro=0.4, Kf=2.2e12, rhoF=1.0,
          perm=1.0e-7, Ncell=10, L=30.0, q=0.1, strip=3.0,
          dt=0.01, T=1.0)  # check time t = 1 s (ADR §7.1)
B4_ALPHA_TARGET = 6.8e-5    # ADR B4: alpha0 ~ 0.25, h = 3 m


def _b4_build(stab):
    """McGann footing, ADR §7.1 B4: 30x30 m half-model, uniform 10x10 mesh
    (h = 3 m = largest EDGE), 3 m strip load on the symmetry edge ramped
    0 -> 0.1 kPa over 0.1 s then held; drainage top surface only; lateral
    u_x rollers, fixed base u_y; Rayleigh 0.05M + 0.02K; Newmark 0.6/0.3025.

    stab: 'off' | ('auto', a0) | ('manual', alpha)."""
    b = B4
    N = b["Ncell"]
    h = b["L"] / N
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    nid = {}
    k = 1
    for j in range(N + 1):
        for i in range(N + 1):
            ops.node(k, i * h, j * h)
            nid[(i, j)] = k
            k += 1
    ops.nDMaterial("ElasticIsotropic", 1, b["E"], b["nu"], b["rho"])
    e = 1
    for j in range(N):
        for i in range(N):
            args = ["LadrunoUP", e, nid[(i, j)], nid[(i + 1, j)],
                    nid[(i + 1, j + 1)], nid[(i, j + 1)], 1,
                    "-thick", 1.0, "-Kf", b["Kf"], "-poro", b["poro"],
                    "-rhoF", b["rhoF"], "-perm", b["perm"], b["perm"],
                    "-dynSeepage", "off"]
            if stab == "off":
                args += ["-stab", "off"]
            elif stab[0] == "auto":
                args += ["-stab", "auto", stab[1]]
            else:
                args += ["-stab", stab[1]]
            ops.element(*args)
            e += 1
    for j in range(N + 1):
        for i in range(N + 1):
            fx = 1 if (i == 0 or i == N) else 0   # rollers (i=0 = symmetry)
            fy = 1 if j == 0 else 0               # fixed base
            fp = 1 if j == N else 0               # drained top ONLY
            if fx or fy or fp:
                ops.fix(nid[(i, j)], fx, fy, fp)
    ops.timeSeries("Path", 1, "-time", 0.0, 0.1, 1.0e6,
                   "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    # 3 m strip = the first top element; consistent nodal loads q*strip/2
    ops.load(nid[(0, N)], 0.0, -b["q"] * b["strip"] / 2.0, 0.0)
    ops.load(nid[(1, N)], 0.0, -b["q"] * b["strip"] / 2.0, 0.0)
    ops.rayleigh(0.05, 0.02, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-8, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    return nid


def _b4_cb_indices(nid):
    """(CB_lat, CB_lap) over the interior nodes at the current state
    (definitions in the module header)."""
    N = B4["Ncell"]
    P = {(i, j): ops.nodeDisp(nid[(i, j)], 3)
         for i in range(N + 1) for j in range(N + 1)}
    inter = [(i, j) for i in range(1, N) for j in range(1, N)]
    p = np.array([P[k] for k in inter])
    assert np.all(np.isfinite(p))
    chi = np.array([(-1.0) ** (i + j) for (i, j) in inter])
    cb_lat = abs(float(p @ chi)) / max(
        float(np.linalg.norm(p) * np.linalg.norm(chi)), 1e-30)
    lap = np.array([P[(i, j)] - 0.25 * (P[(i + 1, j)] + P[(i - 1, j)]
                                        + P[(i, j + 1)] + P[(i, j - 1)])
                    for (i, j) in inter])
    cb_lap = float(np.linalg.norm(lap)) / max(float(np.linalg.norm(p)), 1e-30)
    return cb_lat, cb_lap


def _b4_run(stab):
    nid = _b4_build(stab)
    ns = int(round(B4["T"] / B4["dt"]))
    for s in range(ns):
        assert ops.analyze(1, B4["dt"]) == 0, f"B4 {stab} failed at step {s}"
    return _b4_cb_indices(nid)


@pytest.fixture(scope="module")
def b4_sweep():
    """One pass over {off, 0.05, 0.25, 0.5} (+ the manual-alpha twin of 0.25)
    shared by the B4 tests below."""
    out = {"off": _b4_run("off")}
    for a0 in (0.05, 0.25, 0.5):
        out[a0] = _b4_run(("auto", a0))
    # manual twin: same physics with the alpha recomputed OUTSIDE the element
    alpha = _stab_alpha(B4["L"] / B4["Ncell"], B4["E"], B4["nu"], 0.25)
    out["manual"] = _b4_run(("manual", alpha))
    out["alpha_manual"] = alpha
    return out


def test_b4_checkerboard_off_vs_auto(b4_sweep):
    """-stab off checkerboards (CB HIGH); -stab auto 0.25 suppresses (CB LOW).

    Gates pinned from the first passing run per the ADR §7.1 B4 protocol:
    CB_lap(off) >= 0.5 (measured 0.603), CB_lap(auto 0.25) <= 0.3 (measured
    0.126), ratio >= 3 (measured 4.8); lattice-cosine separation
    CB_lat(off)/CB_lat(auto) >= 2 (measured 0.138 / 0.0264 = 5.2)."""
    lat_off, lap_off = b4_sweep["off"]
    lat_auto, lap_auto = b4_sweep[0.25]
    assert lap_off >= 0.5, f"expected strong checkerboarding: {lap_off:.3f}"
    assert lap_auto <= 0.3, f"stabilization failed to suppress: {lap_auto:.3f}"
    assert lap_off / lap_auto >= 3.0
    assert lat_off / max(lat_auto, 1e-30) >= 2.0, (
        f"lattice-mode separation lost: off {lat_off:.3e} auto {lat_auto:.3e}")


def test_b4_alpha0_sweep_monotone(b4_sweep):
    """alpha0 in {0.05, 0.25, 0.5}: monotone CB_lap suppression
    off > 0.05 > 0.25 > 0.5 (measured 0.603 > 0.222 > 0.126 > 0.0946).

    DOCUMENTED, not asserted: (i) CB_lat is also monotone HERE (0.138 >
    0.0531 > 0.0264 > 0.0199) but went non-monotone in a pulse-load
    scaffolding variant — the normalized lattice cosine is fragile when
    heavy smoothing shrinks ||p|| faster than the residual projection, hence
    CB_lap stays the primary gate; (ii) at alpha0 = 0.5 CB_lap = 0.095 is
    the flattest field of the sweep — on this elastic footing no
    over-smoothing pathology is visible at t = 1 s (the papers' own guidance
    stops at 'limited sensitivity studies', ADR §3.3)."""
    seq = [b4_sweep["off"][1], b4_sweep[0.05][1],
           b4_sweep[0.25][1], b4_sweep[0.5][1]]
    for k in range(1, len(seq)):
        assert seq[k] < seq[k - 1], (
            f"CB_lap not monotone across the alpha0 sweep: {seq}")


def test_b4_stab_alpha_value(b4_sweep):
    """The element's auto stabilization equals the ADR-pinned value.

    (a) recomputed alpha = a0*h^2/(Ks + 4G/3) with h = largest EDGE = 3 m
        exactly on this mesh: 6.686e-5, within 2% of the ADR's 6.8e-5 target
        (their alpha0 ~ 0.254 back-calculated; ours pinned at 0.25);
    (b) the element actually USES it: '-stab auto 0.25' and '-stab <manual
        6.686e-5>' produce the same CB indices to 1e-10 (there is no direct
        alpha-reporting response, so the twin run is the observable)."""
    alpha = b4_sweep["alpha_manual"]
    assert abs(alpha - B4_ALPHA_TARGET) <= 0.05 * B4_ALPHA_TARGET, (
        f"recomputed stab alpha {alpha:.4e} vs ADR target {B4_ALPHA_TARGET}")
    lat_a, lap_a = b4_sweep[0.25]
    lat_m, lap_m = b4_sweep["manual"]
    assert abs(lat_a - lat_m) <= 1e-10 * max(abs(lat_m), 1e-30) + 1e-14
    assert abs(lap_a - lap_m) <= 1e-10 * max(abs(lap_m), 1e-30) + 1e-14


# ==========================================================================
# 4. bbar behavior gate
# ==========================================================================
def _drained_footing_settle(form, nu, N=8, L=8.0, q=-1.0):
    """Static drained footing (p fixed everywhere -> pure solid path, the
    honest-p contract's genuinely-new static capability): N x N mesh, strip
    load on the top-left 2 elements, corner settlement returned."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    h = L / N
    nid = {}
    k = 1
    for j in range(N + 1):
        for i in range(N + 1):
            ops.node(k, i * h, j * h)
            nid[(i, j)] = k
            k += 1
    ops.nDMaterial("ElasticIsotropic", 1, 1e4, nu, 0.0)
    e = 1
    for j in range(N):
        for i in range(N):
            ops.element("LadrunoUP", e, nid[(i, j)], nid[(i + 1, j)],
                        nid[(i + 1, j + 1)], nid[(i, j + 1)], 1,
                        "-Kf", 2.2e5, "-poro", 0.4, "-rhoF", 1.0,
                        "-perm", 1e-4, 1e-4, "-stab", "off",
                        "-formulation", form, "-dynSeepage", "off")
            e += 1
    for j in range(N + 1):
        for i in range(N + 1):
            fx = 1 if (i == 0 or i == N) else 0
            fy = 1 if j == 0 else 0
            ops.fix(nid[(i, j)], fx, fy, 1)       # p fixed EVERYWHERE
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(3):
        w = 1.0 if 0 < i < 2 else 0.5
        ops.load(nid[(i, N)], 0.0, q * w * h, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-10, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return ops.nodeDisp(nid[(0, N)], 2)


def test_bbar_relieves_nearly_incompressible_locking():
    """Near-incompressible skeleton (nu = 0.499): std Q4 locks volumetrically,
    '-formulation bbar' (mean dilatation) relieves it. Measured settlement
    ratios bbar/std: 2.31 at nu = 0.499 (gate >= 1.5), 1.008 at nu = 0.3
    (sanity: bbar must NOT be globally softer away from the incompressible
    limit — gate within [0.95, 1.15])."""
    s_std = _drained_footing_settle("std", 0.499)
    s_bbar = _drained_footing_settle("bbar", 0.499)
    assert s_std < 0 and s_bbar < 0
    ratio_inc = s_bbar / s_std
    assert ratio_inc >= 1.5, (
        f"bbar failed to relieve locking: bbar/std = {ratio_inc:.3f}")
    s_std3 = _drained_footing_settle("std", 0.3)
    s_bbar3 = _drained_footing_settle("bbar", 0.3)
    ratio_c = s_bbar3 / s_std3
    assert 0.95 <= ratio_c <= 1.15, (
        f"bbar deviates at compressible nu=0.3: bbar/std = {ratio_c:.3f}")


# --------------------------------------------------------------------------
# FD consistency of the assembled tangent (std AND bbar) — the bbar-Q
# mean-gradient coupling is the 1.B watch-item this pins.
# --------------------------------------------------------------------------
_FD_NODES = {1: (0.0, 0.0), 2: (1.2, 0.1), 3: (1.1, 1.0), 4: (0.05, 0.9)}
_FD_FIX = {1: (1, 1, 1), 2: (0, 1, 0), 3: (0, 0, 0), 4: (0, 0, 0)}
_FD_FREE = [(n, d) for n in sorted(_FD_NODES) for d in range(3)
            if _FD_FIX[n][d] == 0]
_FD_ORDER = [(n, d) for n in sorted(_FD_NODES) for d in range(3)]


def _fd_elem(form):
    ops.nDMaterial("ElasticIsotropic", 1, 1e4, 0.3, 2.0)
    ops.element("LadrunoUP", 1, 1, 2, 3, 4, 1,
                "-Kf", 2.2e5, "-poro", 0.4, "-rhoF", 1.0,
                "-perm", 1e-4, 1e-4, "-stab", "off",
                "-formulation", form, "-dynSeepage", "off")


def _fd_K_analytic(form):
    """Free-free static tangent block from printA + nodeDOFs.

    printA('-ret') returns the RAW OpenSees Matrix buffer, which is
    COLUMN-major (OpenSeesCommands.cpp:2590 passes &A(0,0) flat) — index it
    A[row + N*col] or the unsymmetric -Q block lands transposed and the FD
    comparison false-fails at ~|Q|/|K| ~ 2.6e-5 (bug we hit in scaffolding)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for t, (x, y) in _FD_NODES.items():
        ops.node(t, x, y)
    _fd_elem(form)
    for t in _FD_NODES:
        ops.fix(t, *_FD_FIX[t])
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormUnbalance", 1e-12, 10)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    A = ops.printA("-ret")
    N = int(round(math.sqrt(len(A))))
    assert N * N == len(A)
    eq = {(n, d): ops.nodeDOFs(n)[d] for (n, d) in _FD_FREE}
    nf = len(_FD_FREE)
    K = np.empty((nf, nf))
    for a in range(nf):
        for b in range(nf):
            K[a, b] = A[eq[_FD_FREE[a]] + N * eq[_FD_FREE[b]]]  # column-major
    return K


def _fd_residual(form, u):
    """Element resisting force at an IMPOSED nodal state: sp on every DOF
    under a Penalty handler + one static solve (triggers element update();
    setNodeDisp alone does NOT re-run update, so eleForce would be stale)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for t, (x, y) in _FD_NODES.items():
        ops.node(t, x, y)
    _fd_elem(form)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for (n, d), v in u.items():
        ops.sp(n, d + 1, v)
    ops.constraints("Penalty", 1e14, 1e14)
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return np.array(ops.eleForce(1))


@pytest.mark.parametrize("form", ["std", "bbar"])
def test_tangent_fd_consistency(form):
    """Central-FD of the static residual (getResistingForce via eleForce at
    imposed states) vs the assembled tangent (printA) on one distorted Q4,
    ALL free DOFs incl. the p columns — so it checks K, -Q (u-row/p-col), H,
    AND the structural ZERO of the p-row/u-col block (honest-p contract,
    ADR §3.2). eps = 1e-4 (exactly linear problem — no truncation error).
    Gate 1e-6 rel; first-passing-run: std 1.5e-10, bbar 1.6e-10."""
    Ka = _fd_K_analytic(form)
    idx = {k: i for i, k in enumerate(_FD_ORDER)}
    nf = len(_FD_FREE)
    Kfd = np.empty((nf, nf))
    eps = 1e-4
    for c, kc in enumerate(_FD_FREE):
        up = {k: 0.0 for k in _FD_ORDER}
        um = dict(up)
        up[kc] = eps
        um[kc] = -eps
        Rp = _fd_residual(form, up)
        Rm = _fd_residual(form, um)
        for r, kr in enumerate(_FD_FREE):
            Kfd[r, c] = (Rp[idx[kr]] - Rm[idx[kr]]) / (2.0 * eps)
    scale = float(np.max(np.abs(Ka)))
    rel = float(np.max(np.abs(Kfd - Ka))) / scale
    assert rel <= 1e-6, f"{form}: tangent-vs-FD rel {rel:.3e} > 1e-6"
    # the unsymmetric signature must actually be there: -Q lives in the
    # u-row/p-col block while p-row/u-col is structurally zero
    asym = float(np.max(np.abs(Ka - Ka.T)))
    assert asym > 1e-3, (
        "tangent unexpectedly symmetric — the -Q block is missing "
        f"(asym {asym:.3e}); wrong-solver silent-drop suspect")
