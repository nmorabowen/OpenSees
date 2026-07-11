r"""LadrunoUP (ELE 33017) — ADR-71 P2 T3-lane Zone-B battery (WP2.A, OPUS-T3).

The linear-triangle (T3, ndm=2 ndf=3, 3 nodes, equal-order p) lane of the
unified Biot u-p element. T3 has parsed+run since P1; this battery GATES it
per ADR-71 §7 P2 row ⟨scope-F12⟩ + plan §WP2.A. Three gates, all on the
CROSSED-diagonal ("union-jack") triangulation of a structured quad footing:
each quad cell -> 4 T3s about an added cell-center node.

    quad cell (i,j)          crossed split (center C = cell midpoint)
     NW ---- NE               NW ------ NE
     |        |               | \  T3  / |
     |        |     --->       |  \   /  |   T3 = {SW,SE,C},{SE,NE,C},
     |        |               |T3  C  T3|        {NE,NW,C},{NW,SW,C}
     SW ---- SE               | /  T3 \ |        (all CCW, center last)
                              SW ------ SE

Element size h = LARGEST TRIANGLE EDGE = the cell side (3 m on B4): every
crossed triangle has edges {cell-side, half-diagonal~2.12, half-diagonal~2.12}
so max edge == cell side for all four (LadrunoUPShapes.h T3::elemSizeH). The
§3.3 auto-alpha therefore transfers to simplices with h = 3 m unchanged — the
twin-identity gate below pins that transfer.

===========================================================================
GATE 1  (test_b4_t3_*): B4-T3 crossed checkerboard + alpha0 sweep + twin id
===========================================================================
The P1 B4 McGann footing (ADR §7.1: 30x30 m half-model, 10x10 cells h=3 m,
3 m strip ramp 0->0.1 kPa/0.1 s then hold, drained top ONLY, Rayleigh
0.05M+0.02K, Newmark 0.6/0.3025, Kf=2.2e12 kPa near-incompressible, to t=1 s),
RE-MESHED with the crossed T3 split.

  ** HEADLINE FINDING (reported to MAIN) ** the crossed-T3 spurious pressure
  mode is CENTER-NODE-LOCALIZED, NOT the corner checkerboard Q4 shows. The
  union-jack center node behaves like a bubble that goes unstable without
  stabilization while the CORNER lattice stays smooth. Consequences:

    * The plan-pinned "CB_lap on the CORNER lattice (comparable to Q4)"
      metric UNDER-REPORTS the instability: measured corner CB_lap
      off=0.160 vs auto0.25=0.119 (ratio only 1.34, and NON-monotone at
      alpha0=0.05) — nothing like Q4's off=0.603 (equiv B4). It is computed
      and recorded here (comparable to P1's Q4 numbers, per the pin) but is
      NOT the suppression gate.
    * The PRIMARY gate is the cell-center roughness CB_cell — the natural
      interior-Laplacian analog for the crossed topology:
          r_cell = p_center - mean(p of the 4 cell corners)   (all cells)
          CB_cell = ||r_cell|| / ||p||   (p over interior corners + centers)
      This is exactly where the crossed mode lives and it is clean.

  First-passing-run values (this build, UmfPack, Newmark 0.6/0.3025, dt=0.01,
  t=1 s, sustained ramp):
      CB_cell:   off=0.7975  a0=0.05: 0.2591  a0=0.25: 0.0969  a0=0.5: 0.0616
                 (monotone; off/0.25 = 8.23)
      CB_corner: off=0.1597  a0=0.05: 0.1599  a0=0.25: 0.1191  a0=0.5: 0.0931
                 (documented, under-reports; corner mode is not the crossed
                  instability)
  auto-alpha (h=3 m) = 0.25*9/(Ks+4Gs/3) = 6.686e-5 (Ks=E/3(1-2nu),
  Gs=E/2(1+nu); E=25000 nu=0.3) — within 2% of the ADR B4 target 6.8e-5.
  '-stab auto 0.25' == '-stab 6.686e-5' fields identical to 0.0 (machine).

===========================================================================
GATE 2  (test_t3_drained_locking_*): T3 honest-baseline locking, MEASURED
===========================================================================
Static DRAINED footing settlement (p fixed to 0 everywhere => pure solid
path), crossed-T3 vs Q4 reference, at nu=0.3 AND nu=0.49; corner settlement
under the 3 m strip. DOCUMENTATION gate: the measured ratios are the
deliverable; the asserts only keep them honest.

  ** MEASURED, and it REFUTES the plan-pinned direction against std-Q4 **
  Settlement (uy at the loaded symmetry corner), first-passing run:
      nu=0.3 : Q4std -3.819e-5  Q4bbar -3.970e-5  T3x -3.811e-5
      nu=0.49: Q4std -2.025e-5  Q4bbar -2.697e-5  T3x -2.557e-5
  Ratios:
      T3x / Q4std : 0.998 (nu .3) -> 1.262 (nu .49)   RATIO RISES, exceeds 1
      T3x / Q4bbar: 0.960 (nu .3) -> 0.948 (nu .49)   below 1, degrades
      Q4std/Q4bbar: 0.962 (nu .3) -> 0.751 (nu .49)   std-Q4 locks hard

  The pin expected "T3 settlement < Q4, ratio degrades further below 1 with
  nu" (crossed-T3 locking worse). Against std (full-2x2) Q4 that is FALSE:
  full-integration Q4 volumetrically LOCKS harder than the union-jack
  crossed-T3 (whose center bubble relieves dilatation), so T3/Q4std RISES
  above 1 at nu->0.5. The pin's spirit ONLY holds against a locking-FREE
  reference (bbar Q4): T3x/Q4bbar sits below 1 at both nu and degrades
  0.960 -> 0.948, i.e. crossed-T3 does carry a mild inherited volumetric
  stiffening, just milder than std-Q4's. The asserted gate is therefore
  stated against the locking-free bbar reference; the std-Q4 refutation is
  pinned in-line for the ADR/guide to reword ("T3 crossed locks, but LESS
  than full-integration Q4 — the union-jack bubble helps").

  Undrained-transient flavor (sealed, fast 0.01 s load, near-incompressible,
  T3 with -stab auto vs std-Q4): T3/Q4std = 1.15 (nu .3) -> 1.43 (nu .49) —
  std-Q4 undrained-locks; documented, lightly asserted (ratio rises, >1).

===========================================================================
GATE 3  (test_t3_terzaghi_*): T3 crossed Terzaghi consolidation smoke
===========================================================================
1x10 stacked column of crossed cells (10 cells -> 40 T3 + 10 centers), ux=0
confined, uy=0 base, drained top; instantaneous load; '-stab auto' DEFAULT
ON (equal-order T3 needs it). Loose gates: no-blowup + right decay CLASS,
not accuracy records. Compared to the numpy Terzaghi series.
  First-passing (stab auto 0.25, dt=0.02): max rel p err Tv0.2=0.05%,
  Tv0.5=0.01%; U_fem 0.352/0.505/0.764 vs series 0.357/0.504/0.764;
  settlement monotone. (Runs equally clean with -stab off: 0.11%/0.02%.)

Run mechanics (worktree, per project memory "OpenSees Python test env"):
the module bootstraps <worktree>/dist/bin onto PATH + DLL dirs + sys.path,
evicts any boot-.pth-preloaded stale pyd, and asserts the loaded pyd is THIS
worktree's. Solver = UmfPack throughout (honest-p tangent is unsymmetric,
ADR §3.2). p is read as nodeDisp slot-3 (LadrunoUP honest-p contract).
Invoke: py -3.12 -m pytest tests/test_ladruno_up_element_t3.py -v
"""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

# --------------------------------------------------------------------------
# bootstrap: THIS worktree's engine (evict any boot-.pth preloaded stale pyd)
# --------------------------------------------------------------------------
_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not (os.path.isfile(os.path.join(_DIST, "opensees.pyd"))
        or os.path.isfile(os.path.join(_DIST, "opensees.so"))):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
try:
    os.add_dll_directory(_DIST)
except (FileNotFoundError, OSError, AttributeError):
    pass
if _DIST not in sys.path:
    sys.path.insert(0, _DIST)
for _m in ("opensees", "openseespy", "openseespy.opensees"):
    sys.modules.pop(_m, None)
import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST}); run with "
    "py -3.12 so the boot .pth cannot preload another build")

pytestmark = [pytest.mark.zone_b]

_UMF = ("UmfPack",)   # honest-p tangent is unsymmetric (ADR §3.2)


# --------------------------------------------------------------------------
# shared numpy oracles + auto-alpha formula (matches the element extraction:
# Ks = E/3(1-2nu) [= lambda + 2G/3], Gs = E/2(1+nu), per WP1.A pins)
# --------------------------------------------------------------------------
def _stab_alpha(h, E, nu, alpha0):
    Ks = E / (3.0 * (1.0 - 2.0 * nu))
    Gs = E / (2.0 * (1.0 + nu))
    return alpha0 * h * h / (Ks + 4.0 * Gs / 3.0)


def terzaghi_series(Z, Tv, nterms=300):
    """Excess-pore-pressure ratio p/p0; Z = z/H from the DRAINED face."""
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return np.sum((2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv))


def terzaghi_U(Tv, nterms=300):
    m = np.arange(nterms)
    M = (2 * m + 1) * np.pi / 2.0
    return 1.0 - np.sum((2.0 / M**2) * np.exp(-M * M * Tv))


# --------------------------------------------------------------------------
# crossed-diagonal T3 mesh builder (Nx x Ny quad cells -> 4*Nx*Ny triangles)
# returns (corner{(i,j):tag}, center{(i,j):tag})
# --------------------------------------------------------------------------
def _crossed_mesh(Nx, Ny, h, E, nu, rho, stab, Kf, poro, rhoF, perm,
                  dyn_seepage="off", thick=1.0, matTag=1):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    corner = {}
    k = 1
    for j in range(Ny + 1):
        for i in range(Nx + 1):
            ops.node(k, i * h, j * h)
            corner[(i, j)] = k
            k += 1
    center = {}
    for j in range(Ny):
        for i in range(Nx):
            ops.node(k, (i + 0.5) * h, (j + 0.5) * h)
            center[(i, j)] = k
            k += 1
    ops.nDMaterial("ElasticIsotropic", matTag, E, nu, rho)
    e = 1

    def _tri(a, b, c):
        nonlocal e
        args = ["LadrunoUP", e, a, b, c, matTag, "-thick", thick,
                "-Kf", Kf, "-poro", poro, "-rhoF", rhoF,
                "-perm", perm, perm, "-dynSeepage", dyn_seepage]
        if stab == "off":
            args += ["-stab", "off"]
        elif stab[0] == "auto":
            args += ["-stab", "auto", stab[1]]
        else:                      # ("manual", value)
            args += ["-stab", stab[1]]
        ops.element(*args)
        e += 1

    for j in range(Ny):
        for i in range(Nx):
            SW = corner[(i, j)]
            SE = corner[(i + 1, j)]
            NE = corner[(i + 1, j + 1)]
            NW = corner[(i, j + 1)]
            C = center[(i, j)]
            # all CCW with the center node last (detJ>0 for T3 provider)
            _tri(SW, SE, C)
            _tri(SE, NE, C)
            _tri(NE, NW, C)
            _tri(NW, SW, C)
    return corner, center


# standard bilinear Q4 reference mesh (for the locking gate)
def _quad_mesh(N, h, E, nu, rho, form, Kf, poro, rhoF, perm, matTag=1):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    corner = {}
    k = 1
    for j in range(N + 1):
        for i in range(N + 1):
            ops.node(k, i * h, j * h)
            corner[(i, j)] = k
            k += 1
    ops.nDMaterial("ElasticIsotropic", matTag, E, nu, rho)
    e = 1
    for j in range(N):
        for i in range(N):
            ops.element("LadrunoUP", e, corner[(i, j)], corner[(i + 1, j)],
                        corner[(i + 1, j + 1)], corner[(i, j + 1)], matTag,
                        "-thick", 1.0, "-Kf", Kf, "-poro", poro, "-rhoF", rhoF,
                        "-perm", perm, perm, "-stab", "off",
                        "-formulation", form, "-dynSeepage", "off")
            e += 1
    return corner


# ==========================================================================
# GATE 1 — B4-T3 crossed checkerboard + alpha0 sweep + twin identity
# ==========================================================================
B4 = dict(E=25000.0, nu=0.3, rho=2.67, poro=0.4, Kf=2.2e12, rhoF=1.0,
          perm=1.0e-7, N=10, L=30.0, q=0.1, strip=3.0, dt=0.01, T=1.0)
B4_ALPHA_TARGET = 6.8e-5      # ADR B4: alpha0 ~ 0.25, h = 3 m


def _b4_build(stab):
    b = B4
    N = b["N"]
    h = b["L"] / N
    corner, center = _crossed_mesh(N, N, h, b["E"], b["nu"], b["rho"], stab,
                                   b["Kf"], b["poro"], b["rhoF"], b["perm"])
    for (i, j), n in corner.items():
        fx = 1 if (i == 0 or i == N) else 0    # rollers (i=0 symmetry, i=N far)
        fy = 1 if j == 0 else 0                # fixed base
        fp = 1 if j == N else 0                # drained top ONLY
        if fx or fy or fp:
            ops.fix(n, fx, fy, fp)
    ops.timeSeries("Path", 1, "-time", 0.0, 0.1, 1.0e6, "-values", 0.0, 1.0, 1.0)
    ops.pattern("Plain", 1, 1)
    # 3 m strip = the first top cell; consistent nodal loads q*strip/2
    ops.load(corner[(0, N)], 0.0, -b["q"] * b["strip"] / 2.0, 0.0)
    ops.load(corner[(1, N)], 0.0, -b["q"] * b["strip"] / 2.0, 0.0)
    ops.rayleigh(0.05, 0.02, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system(*_UMF)
    ops.test("NormUnbalance", 1e-8, 50)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    return corner, center


def _b4_cb(corner, center):
    """(CB_corner, CB_cell): corner-lattice CB_lap (documented, comparable to
    Q4) and the cell-center roughness CB_cell (PRIMARY — where the crossed-T3
    spurious mode lives). Definitions in the module docstring."""
    N = B4["N"]
    Pc = {(i, j): ops.nodeDisp(corner[(i, j)], 3)
          for i in range(N + 1) for j in range(N + 1)}
    Pk = {(i, j): ops.nodeDisp(center[(i, j)], 3)
          for i in range(N) for j in range(N)}
    assert all(np.isfinite(v) for v in Pc.values())
    assert all(np.isfinite(v) for v in Pk.values())
    # corner-lattice interior Laplacian roughness (11x11 interior = 9x9)
    cint = [(i, j) for i in range(1, N) for j in range(1, N)]
    pc = np.array([Pc[k] for k in cint])
    lap = np.array([Pc[(i, j)] - 0.25 * (Pc[(i + 1, j)] + Pc[(i - 1, j)]
                                         + Pc[(i, j + 1)] + Pc[(i, j - 1)])
                    for (i, j) in cint])
    cb_corner = float(np.linalg.norm(lap)) / max(float(np.linalg.norm(pc)), 1e-30)
    # cell-center roughness: center residual vs its 4-corner mean, all cells,
    # normalized by all free-p interior nodes (interior corners + all centers)
    cells = [(i, j) for i in range(N) for j in range(N)]
    res = np.array([Pk[(i, j)] - 0.25 * (Pc[(i, j)] + Pc[(i + 1, j)]
                                         + Pc[(i + 1, j + 1)] + Pc[(i, j + 1)])
                    for (i, j) in cells])
    ref = np.array([Pc[k] for k in cint] + [Pk[c] for c in cells])
    cb_cell = float(np.linalg.norm(res)) / max(float(np.linalg.norm(ref)), 1e-30)
    return cb_corner, cb_cell


def _b4_run(stab):
    corner, center = _b4_build(stab)
    ns = int(round(B4["T"] / B4["dt"]))
    for s in range(ns):
        assert ops.analyze(1, B4["dt"]) == 0, f"B4-T3 {stab} failed at step {s}"
    return _b4_cb(corner, center)


@pytest.fixture(scope="module")
def b4_sweep():
    """One pass over {off, 0.05, 0.25, 0.5} + the manual-alpha twin of 0.25."""
    out = {"off": _b4_run("off")}
    for a0 in (0.05, 0.25, 0.5):
        out[a0] = _b4_run(("auto", str(a0)))
    alpha = _stab_alpha(B4["L"] / B4["N"], B4["E"], B4["nu"], 0.25)   # h = 3 m
    out["manual"] = _b4_run(("manual", repr(alpha)))
    out["alpha_manual"] = alpha
    return out


@pytest.mark.t1
def test_b4_t3_cell_checkerboard_monotone_suppression(b4_sweep):
    """The crossed-T3 spurious pressure mode is center-node-localized; the
    cell-center roughness CB_cell is monotonically suppressed by alpha0 and
    off/auto0.25 >> 2.

    First-passing run (module docstring): CB_cell off=0.798 > 0.05:0.259 >
    0.25:0.097 > 0.5:0.062 (fully monotone); off/0.25 = 8.23.
    Gate: monotone off>0.05>0.25>0.5 AND off/0.25 >= 2 (plan §WP2.A pin).
    This pins that the §3.3 alpha transfer to simplices ACTS (h = largest
    triangle edge = 3 m)."""
    cbc = {k: b4_sweep[k][1] for k in ("off", 0.05, 0.25, 0.5)}
    seq = [cbc["off"], cbc[0.05], cbc[0.25], cbc[0.5]]
    for k in range(1, len(seq)):
        assert seq[k] < seq[k - 1], (
            f"CB_cell not monotone across the alpha0 sweep: {seq}")
    ratio = cbc["off"] / max(cbc[0.25], 1e-30)
    assert ratio >= 2.0, f"off/auto0.25 suppression ratio {ratio:.2f} < 2"
    # non-trivial checkerboard actually present with -stab off (a dead field
    # would pass any ratio gate)
    assert cbc["off"] >= 0.3, f"expected a real crossed-T3 checkerboard: {cbc['off']:.3f}"


@pytest.mark.t1
def test_b4_t3_corner_lattice_documented(b4_sweep):
    """DOCUMENTATION leg (plan pin: compute CB on the CORNER lattice, the
    original 11x11 grid, comparable to P1's Q4 numbers). It UNDER-REPORTS
    the crossed instability because that instability is center-localized
    (see the primary gate). Recorded values, first-passing run:
        CB_corner off=0.160  0.05:0.160  0.25:0.119  0.5:0.093
    (contrast Q4 P1: off=0.603). Non-monotone at alpha0=0.05; only the
    endpoint reduction off > alpha0=0.5 is asserted, as an honest floor."""
    cc = {k: b4_sweep[k][0] for k in ("off", 0.05, 0.25, 0.5)}
    assert cc["off"] > cc[0.5], (
        f"corner CB should at least drop off->0.5: {cc['off']:.3f} !> {cc[0.5]:.3f}")
    # documented, not gated: the corner mode is not the crossed instability
    assert cc["off"] < 0.35, (
        "corner CB unexpectedly large — the crossed mode is supposed to be "
        f"center-localized, not a corner checkerboard (got {cc['off']:.3f})")


@pytest.mark.t1
def test_b4_t3_auto_manual_twin_identity(b4_sweep):
    """The §3.3 alpha transfer to simplices, PINNED: '-stab auto 0.25' and
    '-stab <manual 0.25*9/(Ks+4Gs/3)>' produce identical pressure fields
    (there is no direct alpha-reporting response, so the CB indices computed
    from the two runs are the observable). h = largest triangle edge = 3 m.

    (a) recomputed manual alpha = 6.686e-5, within 2% of the ADR B4 target
        6.8e-5 (their alpha0 ~ 0.254 back-calculated; ours pinned at 0.25);
    (b) both CB_corner and CB_cell identical to 1e-10 (measured 0.0)."""
    alpha = b4_sweep["alpha_manual"]
    assert abs(alpha - B4_ALPHA_TARGET) <= 0.05 * B4_ALPHA_TARGET, (
        f"recomputed stab alpha {alpha:.4e} vs ADR target {B4_ALPHA_TARGET}")
    cc_a, cell_a = b4_sweep[0.25]
    cc_m, cell_m = b4_sweep["manual"]
    assert abs(cc_a - cc_m) <= 1e-10 * max(abs(cc_m), 1e-30) + 1e-14, (
        f"auto/manual CB_corner differ: {cc_a:.12e} vs {cc_m:.12e}")
    assert abs(cell_a - cell_m) <= 1e-10 * max(abs(cell_m), 1e-30) + 1e-14, (
        f"auto/manual CB_cell differ: {cell_a:.12e} vs {cell_m:.12e}")


# ==========================================================================
# GATE 2 — T3 honest-baseline locking (documentation gate, MEASURED)
# ==========================================================================
_G2 = dict(E=25000.0, rho=2.67, q=0.1, strip=3.0, N=10, h=3.0,
           Kf=2.2e12, poro=0.4, rhoF=1.0, perm=1.0e-7)


def _static_drained_settle(corner, center=None):
    """Static drained footing (p fixed to 0 everywhere -> pure solid path),
    3 m strip, corner settlement (uy at the loaded symmetry corner)."""
    g = _G2
    N = g["N"]
    for (i, j), n in corner.items():
        fx = 1 if (i == 0 or i == N) else 0
        fy = 1 if j == 0 else 0
        ops.fix(n, fx, fy, 1)               # p fixed EVERYWHERE (drained solid)
    if center:
        for n in center.values():
            ops.fix(n, 0, 0, 1)             # centers: u free, p fixed
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(corner[(0, N)], 0.0, -g["q"] * g["strip"] / 2.0, 0.0)
    ops.load(corner[(1, N)], 0.0, -g["q"] * g["strip"] / 2.0, 0.0)
    ops.system(*_UMF)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormUnbalance", 1e-10, 20)
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0
    return ops.nodeDisp(corner[(0, N)], 2)


@pytest.mark.t1
def test_t3_drained_locking_vs_lockingfree_reference():
    """Crossed-T3 drained settlement vs a LOCKING-FREE reference (bbar Q4) at
    nu=0.3 and nu=0.49 — the honest volumetric-locking measurement.

    Against the locking-free bbar-Q4, crossed-T3 IS stiffer (ratio < 1) and
    the ratio DEGRADES further below 1 as nu->0.49 (the plan-pinned direction,
    against a fair reference). Measured first-passing run:
        T3x/Q4bbar: 0.960 (nu 0.3) -> 0.948 (nu 0.49).
    Generous documentation bands (the numbers are the deliverable)."""
    g = _G2
    out = {}
    for nu in (0.3, 0.49):
        cb = _quad_mesh(g["N"], g["h"], g["E"], nu, g["rho"], "bbar",
                        g["Kf"], g["poro"], g["rhoF"], g["perm"])
        s_bbar = _static_drained_settle(cb)
        co, ce = _crossed_mesh(g["N"], g["N"], g["h"], g["E"], nu, g["rho"],
                               "off", g["Kf"], g["poro"], g["rhoF"], g["perm"])
        s_t3 = _static_drained_settle(co, ce)
        out[nu] = s_t3 / s_bbar
        assert s_bbar < 0 and s_t3 < 0
    # crossed-T3 stiffer than the locking-free reference at both nu
    assert out[0.3] < 1.0, f"T3/bbar at nu=0.3 = {out[0.3]:.4f} (want < 1, stiffer)"
    assert out[0.49] < 1.0, f"T3/bbar at nu=0.49 = {out[0.49]:.4f} (want < 1)"
    # and it degrades (stiffens relatively) toward the incompressible limit
    assert out[0.49] < out[0.3] + 1e-3, (
        f"T3/bbar should not improve toward nu=0.49: {out[0.3]:.4f} -> {out[0.49]:.4f}")
    # generous pins on the measured values
    assert 0.90 <= out[0.3] <= 1.00, f"T3/bbar(0.3) drift: {out[0.3]:.4f}"
    assert 0.88 <= out[0.49] <= 0.99, f"T3/bbar(0.49) drift: {out[0.49]:.4f}"


@pytest.mark.t1
def test_t3_drained_vs_stdq4_refutes_pinned_direction():
    """REFUTATION leg, MEASURED (reported to MAIN). Against the LITERAL pin's
    std (full-2x2) Q4 reference, crossed-T3 is NOT stiffer: full-integration
    Q4 volumetrically LOCKS harder than the union-jack crossed-T3, so the
    T3/Q4std ratio RISES and crosses 1 as nu->0.49. Measured first-passing:
        T3x/Q4std : 0.998 (nu 0.3) -> 1.262 (nu 0.49)
        Q4std/Q4bbar (std-Q4's own locking): 0.962 -> 0.751
    This test PINS the refutation so the ADR/guide reword to 'T3 crossed locks,
    but LESS than full-integration Q4 — the center bubble helps'."""
    g = _G2
    tq = {}
    ql = {}
    for nu in (0.3, 0.49):
        qs = _quad_mesh(g["N"], g["h"], g["E"], nu, g["rho"], "std",
                        g["Kf"], g["poro"], g["rhoF"], g["perm"])
        s_std = _static_drained_settle(qs)
        qb = _quad_mesh(g["N"], g["h"], g["E"], nu, g["rho"], "bbar",
                        g["Kf"], g["poro"], g["rhoF"], g["perm"])
        s_bbar = _static_drained_settle(qb)
        co, ce = _crossed_mesh(g["N"], g["N"], g["h"], g["E"], nu, g["rho"],
                               "off", g["Kf"], g["poro"], g["rhoF"], g["perm"])
        s_t3 = _static_drained_settle(co, ce)
        tq[nu] = s_t3 / s_std
        ql[nu] = s_std / s_bbar
    # the naive T3/std-Q4 ratio RISES with nu and exceeds 1 at nu=0.49
    assert tq[0.49] > tq[0.3], f"T3/Q4std should rise with nu: {tq[0.3]:.3f}->{tq[0.49]:.3f}"
    assert tq[0.49] > 1.0, f"T3 softer than std-Q4 at nu=0.49 expected: {tq[0.49]:.3f}"
    # because std-Q4 itself locks hard at nu->0.5 (its own ratio drops well below 1)
    assert ql[0.49] < 0.85, f"std-Q4 expected to lock at nu=0.49: Q4std/Q4bbar={ql[0.49]:.3f}"
    assert ql[0.49] < ql[0.3], "std-Q4 locking should worsen with nu"


@pytest.mark.t1
def test_t3_undrained_softer_than_stdq4():
    """Undrained-transient flavor (documented): sealed (no p fixed), fast
    0.01 s ramp, near-incompressible; crossed-T3 (-stab auto) vs std-Q4.
    std-Q4 undrained-locks, so T3 settles MORE and the ratio rises with nu.
    Measured: T3/Q4std = 1.15 (nu 0.3) -> 1.43 (nu 0.49). No-blowup smoke +
    honest direction (not a matched-operator identity)."""
    g = _G2
    N = g["N"]

    def _undrained_settle(corner, center=None, matched_stab=False):
        for (i, j), n in corner.items():
            fx = 1 if (i == 0 or i == N) else 0
            fy = 1 if j == 0 else 0
            if fx or fy:
                ops.fix(n, fx, fy, 0)        # sealed: NO p fixed anywhere
        ops.timeSeries("Path", 1, "-time", 0.0, 0.01, 1.0e6,
                       "-values", 0.0, 1.0, 1.0)
        ops.pattern("Plain", 1, 1)
        ops.load(corner[(0, N)], 0.0, -g["q"] * g["strip"] / 2.0, 0.0)
        ops.load(corner[(1, N)], 0.0, -g["q"] * g["strip"] / 2.0, 0.0)
        ops.rayleigh(0.05, 0.02, 0.0, 0.0)
        ops.system(*_UMF)
        ops.numberer("RCM")
        ops.constraints("Transformation")
        ops.test("NormUnbalance", 1e-8, 50)
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
        for _ in range(20):
            assert ops.analyze(1, 0.005) == 0
        return ops.nodeDisp(corner[(0, N)], 2)

    ratio = {}
    for nu in (0.3, 0.49):
        qs = _quad_mesh(N, g["h"], g["E"], nu, g["rho"], "std",
                        g["Kf"], g["poro"], g["rhoF"], g["perm"])
        s_std = _undrained_settle(qs)
        co, ce = _crossed_mesh(N, N, g["h"], g["E"], nu, g["rho"],
                               ("auto", "0.25"), g["Kf"], g["poro"],
                               g["rhoF"], g["perm"])
        s_t3 = _undrained_settle(co, ce)
        assert s_std < 0 and s_t3 < 0 and np.isfinite(s_t3)
        ratio[nu] = s_t3 / s_std
    assert ratio[0.3] > 1.0 and ratio[0.49] > 1.0, (
        f"T3 expected softer than undrained-locking std-Q4: {ratio}")
    assert ratio[0.49] > ratio[0.3], (
        f"ratio should rise with nu (std-Q4 locks harder): {ratio}")


# ==========================================================================
# GATE 3 — T3 crossed Terzaghi consolidation smoke (-stab auto default ON)
# ==========================================================================
@pytest.mark.t1
def test_t3_terzaghi_consolidation_smoke():
    """1x10 crossed-T3 column (40 T3 + 10 centers), ux=0 confined, uy=0 base,
    drained top, instantaneous load; '-stab auto' DEFAULT ON (equal-order T3
    needs it). Loose gates: no-blowup + right decay CLASS, not accuracy
    records.

    First-passing (stab auto 0.25, dt=0.02): max rel p err Tv0.2=0.05%,
    Tv0.5=0.01%; U_fem 0.352/0.505/0.764 vs series 0.357/0.504/0.764;
    settlement monotone increasing. (Also clean under -stab off.)
    Gate: p-profile err <= 5% at Tv in {0.2, 0.5}; U within 0.03 of series;
    settlement monotone."""
    E, nu, rho = 1.0e4, 0.2, 2.0
    Kf, poro, rhoF, kbar = 4444.0, 0.4, 1.0, 1.0e-4
    L, q, nely, h = 10.0, 10.0, 10, 1.0
    Eoed = E * (1 - nu) / ((1 + nu) * (1 - 2 * nu))
    ooQ = poro / Kf
    Qbar = 1.0 / ooQ
    cv = kbar / (1.0 / Eoed + ooQ)
    s0 = q * L / (Eoed + Qbar)        # undrained settlement
    sinf = q * L / Eoed               # drained settlement

    # -stab auto (default ON for equal-order T3) -- the point of this gate
    corner, center = _crossed_mesh(1, nely, h, E, nu, rho, ("auto", "0.25"),
                                   Kf, poro, rhoF, kbar)
    for (i, j), n in corner.items():
        ops.fix(n, 1, (1 if j == 0 else 0), (1 if j == nely else 0))
    for (i, j), n in center.items():
        ops.fix(n, 1, 0, 0)           # ux=0 confined; center uy,p free
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(corner[(i, nely)], 0.0, -q / 2.0, 0.0)
    ops.system(*_UMF)
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")

    assert ops.analyze(1, 1.0e-3) == 0
    t = 1.0e-3
    dt = 0.02
    settle = {}
    perr = {}
    for Tv in (0.1, 0.2, 0.5):
        tt = Tv * L * L / cv
        nst = int(round((tt - t) / dt))
        assert ops.analyze(nst, dt) == 0, f"Terzaghi T3 blew up before Tv={Tv}"
        t += nst * dt
        TvA = t * cv / (L * L)
        errs = [abs(ops.nodeDisp(corner[(0, j)], 3)
                    - (q * Qbar / (Qbar + Eoed))
                    * terzaghi_series((L - L * j / nely) / L, TvA))
                / (q * Qbar / (Qbar + Eoed))
                for j in range(1, nely)]
        perr[Tv] = max(errs)
        settle[Tv] = -ops.nodeDisp(corner[(0, nely)], 2)

    # loose accuracy: right decay CLASS, not records
    assert perr[0.2] <= 0.05, f"Tv=0.2 p err {perr[0.2]:.3%} > 5%"
    assert perr[0.5] <= 0.05, f"Tv=0.5 p err {perr[0.5]:.3%} > 5%"
    # settlement grows monotonically as consolidation proceeds
    assert settle[0.1] < settle[0.2] < settle[0.5], (
        f"settlement not monotone: {settle}")
    # degree of consolidation tracks the series (loose)
    U_fem = {Tv: (settle[Tv] - s0) / (sinf - s0) for Tv in settle}
    for Tv in (0.2, 0.5):
        assert abs(U_fem[Tv] - terzaghi_U(Tv)) < 0.03, (
            f"U(Tv={Tv}) fem {U_fem[Tv]:.4f} vs series {terzaghi_U(Tv):.4f}")
