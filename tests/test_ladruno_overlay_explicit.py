"""LadrunoPorousOverlay (PATTERN 33022) — ADR-73 P3 EXPLICIT-LANE battery (WP3.C).

P3 lane (plan "P3 PINNED at 3.A" §3.1-§3.4; ADR §12 P0 entry item 3):
  * `-fsL zero` on the overlay = the EXPLICIT lane: L == 0, plain backward-Euler
    fluid at each commit through the existing fs1 hook ((S* + dtH) p1 =
    S* p0 - Qt du + dt f_seep), oedometric parser floor bypassed, ONE-TIME loud
    parser advisory (a quasi-static implicit fs1 march with L = 0 is the naive
    drained split — fails loudly, not wrongly).
  * `LadrunoStaggeredAnalyze` REFUSES -fsL zero overlays (loud fatal; iterating
    with L = 0 IS the drained split, KTJ-divergent at soil coupling).
  * `ops.criticalTimeStep()` returns the OVERLAY-AWARE DISCRETE UNDRAINED pencil
    when an overlay covers elements (per-cell Qe Se^-1 Qe^T augmentation in the
    shared per-element scan); with no overlay the drained pencil is unchanged.

FROZEN MEASURED CONSTANTS (plan frozen-constants block, E7.2/E7.3a/E7.4):
  e72: L=0 implicit-fluid CD boundary = 1.000x the discrete undrained pencil;
       FROZEN default operating point 0.5x. (fixed-point-L drained-CFL lane =
       0.987x drained pencil but O(1)-wrong; the rate form is CD-unstable at
       ANY dt — neither is the P3 lane.)
  e74: explicit-fluid (P3b) boundary = 1.32x the pencil; the material factor
       sqrt(1 + Kf/(n*Moed)) is ~1.85x conservative as a dt_cr predictor
       (e72 soil: material 21.4 vs measured discrete 15.4 ~= 0.72x) — the
       advisory uses the DISCRETE pencil, the formula is documentation-only.
  e73a: -subcycle error ~ N^1.2, all N <= 50 stable on the L=0 lane; error
       doubles at theta = N*dt/(h^2/c_v) = 0.089, pinning
       N_auto = max(1, floor(0.089 * h_min^2 / (c_v_max * dt))) with
       c_v = kbar_max * Moed (LadrunoPorousOverlay.cpp resolveSubcycleN).

THE GATES (plan §3.4 a-g + the driver-refusal surface gate):
  (a) TWO-LEG on the ZS84 Example-1 column (ADR §7.1 B2 config, adapted): CD
      (CentralDifferenceLadruno) + plain ndf=2 quads + overlay(-fsL zero,
      -stab auto 0.25) at dt = 0.4x the measured overlay-aware advisory vs
      implicit Newmark(0.6, 0.3025) + LadrunoUP (stab mirrored, -dynSeepage
      off), same mesh/params/loading, matched output cadence: mutual
      dt-convergence over 3 dt levels — diffs shrink, observed order >= 0.9
      on p probes (P1/P2 two-leg methodology; NOT exact equality).
  (b) CFL ENVELOPE on the e72-matched column (E=1e7 nu=0.25 Kf=2.2e9 poro=0.4
      kbar=1e-7/1e4=1e-11, mixture rho=2.0, 1x8 unit-cell mesh) —
      FIXED-HORIZON gates. Run-2 adjudication measurement: the zero-drainage
      UNDAMPED worst case has NO asymptotically stable dt — a slow SECULAR
      ENERGY-PUMPING of the staggered coupling (growth/unit-time ~ O(dt))
      eventually diverges every march (0.4x -> 48344 steps ... 1.0x -> 7761
      ... 3.0x -> 855), while the same solid WITHOUT the overlay is clean
      (drained control bounded 30k steps); run-1's '1.1x/1.4x bounded' was a
      short-horizon detection artifact, NOT an over-conservative advisory
      (e72's '+-25 percent of 1.0x' applied to the toy's global operator —
      refutation + secular record for ADR §12 P3). Gates: {0.6, 0.8, 0.95}x
      survive the 600-step envelope; {1.5, 2.0, 3.0}x ALL diverge within a
      4500-step cap (divergence at or below 3.0x as adjudicated); 0.9x the
      DRAINED pencil (no-overlay twin query) diverges fast (the
      naive-advisory demonstration). Measured divergence steps printed as
      the growth trend for the PR.
  (c) ADVISORY, MERGED single gate (final adjudication): (i) hard
      stab-INVARIANCE identity |ratio_stabOn - ratio_stabOff|/ratio <= 1e-12
      (measured 1.29e-14 — the per-cell constant-p mode governs the
      condensation and lies in the H-tilde nullspace, so the pencil cannot
      see the stab operator); (ii) drained/undrained ratio > 1 and within
      [0.5, 1.3]x the material formula sqrt(1+Kf/(n*Moed)) (measured 26.24
      = 1.2243x material — the per-element factor legitimately exceeds the
      material bound); (iii) diagnosis printed: the toy's discrete 15.36
      (= 0.72x material, E7.4) was the ASSEMBLED global operator — the
      26.24-vs-15.36 gap is OPERATOR CLASS (per-element vs assembled), not
      stab.
  (d) INCUMBENT HEAD-TO-HEAD on ZS84: upstream CentralDifference +
      FourNodeQuadUP (bulk = Qbar ~= Kf/n arg mapping; p via nodeVel dof 3)
      vs the overlay lane, both at the overlay lane's stable dt (retry ladder
      dt/4, dt/16 if the quadUP route needs a smaller dt to run — reported):
      accuracy vs the implicit LadrunoUP reference (rel diff on p probes) +
      per-step wall-clock — REPORT ONLY, no hard perf gate. S->0 DEMO:
      Kf x 1e3 — the OVERLAY leg must stay finite (hard); the incumbent's
      measured symptom is printed loudly, whatever it measures to be.
  (e) RICHARDSON PIN / REFUTATION RECORDER: honest-p LadrunoUP + upstream
      CentralDifference, ZS84 column, TWO legs of 4000 steps at dt = 2e-4
      and dt = 1e-3 (both well below the solid pencil ~8e-3 s), fresh
      children. Run-1 measured the 400-step dt=2e-4 leg TRACKING the
      implicit reference (rel 2.33e-4), so the gate re-measures harder: if
      either leg degrades (non-finite / refuses / rel > 0.5) the
      EXPECTED-BAD pin holds and the manifesting leg is reported; if BOTH
      legs still track (rel <= 0.05) the gate asserts the tracking and
      LOUDLY records the refutation for the quirks-row rewrite; the
      in-between zone hard-fails for MAIN adjudication.
  (f) SUBCYCLE SWEEP under CD: overlay(-fsL zero), fixed dt = 0.5x pencil,
      fast-drain permeability tuned so N_auto ~ 4, N in {1, 2, 5, 10}: all
      bounded; err(N) vs the N=1 reference grows ~monotonically, log-log
      slope in [0.7, 1.7] (E7.3a ~N^1.2). PLUS -subcycle auto: resolved N
      inferred from -record row count vs commit count, within +-1 of the
      theta-formula hand-count N = max(1, floor(0.089*h^2/(cv*dt))).
  (g) ENERGY CLOSURE: ADR-69 EnergyBalanceRecorder on the (a) overlay run:
      |ERR%| <= 2.0 at end-of-march (the plan's fallback bound; the ADR-69
      battery gates drift ratios, not an absolute) and ULW nonzero — proves
      the overlay's +Q*p coupling work rides the ULW external-work channel
      (plan §3.3; ULW integrates v^T P_ext dt via Node::getUnbalancedLoad).
  (h) DRIVER REFUSAL (surface): LadrunoStaggeredAnalyze on an -fsL zero
      overlay -> loud fatal (negative rc OR exception, the P2 e5 dual
      handling) in a fresh subprocess; the -fsL zero parser advisory must
      appear at build time.

RUNNER (worktree dist/bin opensees.pyd is CPython 3.12; ~2-4 min):
  C:\\Users\\nmora\\AppData\\Local\\Python\\pythoncore-3.12-64\\python.exe -S \\
      tests/test_ladruno_overlay_explicit.py          # standalone, exits 1 on fail
  py -3.12 -m pytest tests/test_ladruno_overlay_explicit.py -x -q   # pytest
`-S` skips site (evicts any boot-.pth preloaded stale pyd); the bootstrap re-adds
site-packages for numpy and asserts THIS worktree's engine won.
"""
import math
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# --------------------------------------------------------------------------
# bootstrap (copied from the P2 driver battery, toy import dropped — no toy
# dependence in this file; the oracles are measured pencils + the implicit
# LadrunoUP reference)
# --------------------------------------------------------------------------
# (a) `python -S` skips site: re-add the user + core site-packages for numpy.
for _sp in (
    os.path.join(os.environ.get("APPDATA", ""), "Python", "Python312", "site-packages"),
    r"C:\Users\nmora\AppData\Roaming\Python\Python312\site-packages",
    os.path.join(sys.prefix, "Lib", "site-packages"),
):
    if _sp and os.path.isdir(_sp) and _sp not in sys.path:
        sys.path.append(_sp)

import numpy as np  # noqa: E402

# (b) THIS worktree's engine (evict any boot-.pth preloaded stale pyd).
_ROOT = Path(__file__).resolve().parents[1]
_DIST = str(_ROOT / "dist" / "bin")
_BUILT = os.path.isfile(os.path.join(_DIST, "opensees.pyd"))
if _BUILT:
    os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
    try:
        os.add_dll_directory(_DIST)
    except (FileNotFoundError, OSError):
        pass
    if _DIST not in sys.path:
        sys.path.insert(0, _DIST)

if not _BUILT:
    # Collection-safe skip: under pytest (which IMPORTS every module before
    # filtering marks) use the module-level skip — a module-level sys.exit
    # would crash the whole session during collection. Only a direct
    # `python -S` run exits (cleanly, 0). NB: pytest must NOT be imported by
    # this module before this check — the sys.modules probe is how we tell a
    # pytest collection apart from a standalone run.
    _msg = f"worktree engine not built: {_DIST}"
    if "pytest" in sys.modules:
        import pytest
        pytest.skip(_msg, allow_module_level=True)
    print(f"SKIP: {_msg}")
    raise SystemExit(0)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

try:
    import pytest  # noqa: E402
    _HAVE_PYTEST = True
except Exception:  # pragma: no cover
    _HAVE_PYTEST = False

if _HAVE_PYTEST:
    pytestmark = [pytest.mark.zone_b, pytest.mark.t1]

_TMP = tempfile.mkdtemp(prefix="ov_exp_")


def _rec(name):
    return os.path.join(_TMP, name)


# ==========================================================================
# frozen configurations
# ==========================================================================
# ZS84 Example-1 column (ADR §7.1 B2 pins, test_ladruno_up_element_analytic.py
# ~line 289): 30 m column, q = 1 ramped by sin(pi*t/2tp) tp = 0.1 s, E = 3e4,
# nu = 0.2, rho = 2 (MIXTURE), rhoF = 1, Kf = 1e5, n = 0.3,
# kbar = 1e-2/(rhoF*9.81); top drained, base fixed/impermeable, laterally
# confined (u_x = 0 everywhere). h = 1 m (nely = 30).
_Z = dict(L=30.0, nely=30, E=3.0e4, nu=0.2, rho=2.0, rhoF=1.0,
          Kf=1.0e5, poro=0.3, kbar=1.0e-2 / (1.0 * 9.81), q=1.0, tp=0.1)

# e72-matched column (plan gate (b) pins): E=1e7 nu=0.25 Kf=2.2e9 poro=0.4
# kh=1e-7 gw=1e4 -> kbar=1e-11, mixture rho=2.0; small mesh = 1x8 unit cells
# (h = 1). Material undrained factor sqrt(1+Kf/(n*Moed)) = 21.43 (e72 measured
# discrete 15.4 ~= 0.72x — gate (c)'s band anchor).
_E72 = dict(nely=8, E=1.0e7, nu=0.25, rho=2.0, rhoF=1000.0,
            Kf=2.2e9, poro=0.4, kbar=1.0e-7 / 1.0e4, q=1.0e5)

_THETA = 0.089                       # E7.3a frozen -subcycle auto theta


def _moed(E, nu):
    return E * (1.0 - nu) / ((1.0 + nu) * (1.0 - 2.0 * nu))


def _relL2(A, B):
    A, B = np.asarray(A, float), np.asarray(B, float)
    d = np.linalg.norm(B)
    return float(np.linalg.norm(A - B) / (d if d > 0 else 1.0))


def _read_record(path):
    """Return (times[nrec], P[nrec,Nregion], tag->col) from the overlay CSV."""
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    hdr = [int(t) for t in lines[0].split(",")[1:]]
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return data[:, 0], data[:, 1:], {t: k for k, t in enumerate(hdr)}


def _loglog_slope(xs, ys):
    lx, ly = np.log(np.asarray(xs, float)), np.log(np.asarray(ys, float))
    return float(np.polyfit(lx, ly, 1)[0])


# node tag for column grid position (i in {0,1}, j in 0..nely), 1-based
def _nt(i, j):
    return 2 * j + i + 1


# ==========================================================================
# ZS84 builders (explicit twin, implicit LadrunoUP reference)
# ==========================================================================
def _zs84_sin_ramp():
    """The B2 sin ramp: 21-point Path over [0, tp], held at 1.0 (identical for
    every dt — the Path is dt-independent). math (not numpy) so the child
    incumbent leg reproduces bit-identical series values."""
    z = _Z
    ts = [z["tp"] * k / 20.0 for k in range(21)] + [1.0e12]
    vs = [math.sin(math.pi * t / (2.0 * z["tp"])) for t in ts[:21]] + [1.0]
    ops.timeSeries("Path", 1, "-time", *ts, "-values", *vs)


def _build_zs84_explicit(Kf=None, overlay=True, record=None):
    """ZS84 explicit twin: plain ndf=2 quads carrying the MIXTURE density
    (rho = 2 via the material — no gravity/body forces in ZS84, so the
    LEDGER_quirks 'staggered-twin body-force convention' trap is NOT armed
    here) + LadrunoPorousOverlay(-fsL zero, -stab auto 0.25) + CDL."""
    z = _Z
    Kf = z["Kf"] if Kf is None else Kf
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j in range(z["nely"] + 1):
        for i in range(2):
            ops.node(_nt(i, j), float(i), z["L"] * j / z["nely"])
    ops.nDMaterial("ElasticIsotropic", 1, z["E"], z["nu"], z["rho"])
    eles = []
    for j in range(z["nely"]):
        ops.element("quad", j + 1, _nt(0, j), _nt(1, j), _nt(1, j + 1),
                    _nt(0, j + 1), 1.0, "PlaneStrain", 1)
        eles.append(j + 1)
    for j in range(z["nely"] + 1):
        for i in range(2):
            ops.fix(_nt(i, j), 1, 1 if j == 0 else 0)
    _zs84_sin_ramp()
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(_nt(i, z["nely"]), 0.0, -z["q"] / 2.0)
    if overlay:
        top = [_nt(0, z["nely"]), _nt(1, z["nely"])]
        args = ["-region", *eles, "-Kf", Kf, "-rhoF", z["rhoF"],
                "-perm", z["kbar"], z["kbar"], "-poro", z["poro"],
                "-moduli", z["E"], z["nu"], "-drained", *top,
                "-fsL", "zero", "-stab", "auto", 0.25]
        if record is not None:
            args += ["-record", record, 1]
        ops.pattern("LadrunoPorousOverlay", 77, *args)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")


def _build_zs84_up(Kf=None):
    """Implicit monolithic reference: LadrunoUP + Newmark(0.6, 0.3025), stab
    MIRRORED to the overlay (auto 0.25), -dynSeepage off (the P1/P2 accuracy
    recipe), UmfPack (honest-p tangent is unsymmetric)."""
    z = _Z
    Kf = z["Kf"] if Kf is None else Kf
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for j in range(z["nely"] + 1):
        for i in range(2):
            ops.node(_nt(i, j), float(i), z["L"] * j / z["nely"])
    ops.nDMaterial("ElasticIsotropic", 1, z["E"], z["nu"], z["rho"])
    for j in range(z["nely"]):
        ops.element("LadrunoUP", j + 1, _nt(0, j), _nt(1, j), _nt(1, j + 1),
                    _nt(0, j + 1), 1, "-Kf", Kf, "-poro", z["poro"],
                    "-rhoF", z["rhoF"], "-perm", z["kbar"], z["kbar"],
                    "-stab", "auto", 0.25, "-dynSeepage", "off")
    for j in range(z["nely"] + 1):
        for i in range(2):
            ops.fix(_nt(i, j), 1, 1 if j == 0 else 0, 0)
    for i in range(2):
        ops.fix(_nt(i, z["nely"]), 0, 0, 1)          # drained top p = 0
    _zs84_sin_ramp()
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(_nt(i, z["nely"]), 0.0, -z["q"] / 2.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _march_up_samples(dt, nsteps, sample_commits, probe_nodes, Kf=None):
    """March the implicit LadrunoUP reference; return p[len(samples), nprobe]
    (nodeDisp dof 3 — the honest-p contract)."""
    _build_zs84_up(Kf=Kf)
    sset = set(sample_commits)
    out = {}
    for s in range(1, nsteps + 1):
        assert ops.analyze(1, dt) == 0, f"implicit UP reference failed at step {s}"
        if s in sset:
            out[s] = [ops.nodeDisp(nt, 3) for nt in probe_nodes]
    arr = np.array([out[s] for s in sample_commits])
    assert np.all(np.isfinite(arr)), "implicit reference produced non-finite p"
    return arr


# ==========================================================================
# e72 column builder (parent side; the child templates below carry their own)
# ==========================================================================
def _build_e72_explicit(kbar=None, overlay=True, subcycle=None, record=None,
                        stab=("auto", 0.25)):
    c = _E72
    kbar = c["kbar"] if kbar is None else kbar
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j in range(c["nely"] + 1):
        for i in range(2):
            ops.node(_nt(i, j), float(i), float(j))
    ops.nDMaterial("ElasticIsotropic", 1, c["E"], c["nu"], c["rho"])
    eles = []
    for j in range(c["nely"]):
        ops.element("quad", j + 1, _nt(0, j), _nt(1, j), _nt(1, j + 1),
                    _nt(0, j + 1), 1.0, "PlaneStrain", 1)
        eles.append(j + 1)
    for j in range(c["nely"] + 1):
        for i in range(2):
            ops.fix(_nt(i, j), 1, 1 if j == 0 else 0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(_nt(i, c["nely"]), 0.0, -c["q"] / 2.0)
    if overlay:
        top = [_nt(0, c["nely"]), _nt(1, c["nely"])]
        args = ["-region", *eles, "-Kf", c["Kf"], "-rhoF", c["rhoF"],
                "-perm", kbar, kbar, "-poro", c["poro"],
                "-moduli", c["E"], c["nu"], "-drained", *top,
                "-fsL", "zero", "-stab", *stab]
        if subcycle == "auto":
            args += ["-subcycle", "auto"]
        elif subcycle is not None:
            args += ["-subcycle", int(subcycle)]
        if record is not None:
            args += ["-record", record, 1]
        ops.pattern("LadrunoPorousOverlay", 77, *args)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")


# ==========================================================================
# measured pencils (memoized) — the CDL-7 idiom: dt_cr is computed
# unconditionally at the first newStep (ADR C4), so one tiny priming step
# makes ops.criticalTimeStep() valid.
# ==========================================================================
_PEN = {}


def _pencil(key):
    if key in _PEN:
        return _PEN[key]
    builders = {
        "zs84_ov": lambda: _build_zs84_explicit(),
        "zs84_ov_kf1e3": lambda: _build_zs84_explicit(Kf=_Z["Kf"] * 1.0e3),
        "e72_ov": lambda: _build_e72_explicit(),
        "e72_ov_soff": lambda: _build_e72_explicit(stab=("off",)),
        "e72_dr": lambda: _build_e72_explicit(overlay=False),
    }
    builders[key]()
    assert ops.analyze(1, 1.0e-9) == 0, f"pencil priming step failed ({key})"
    dtcr = ops.criticalTimeStep()
    ops.wipe()
    assert dtcr is not None and math.isfinite(dtcr) and dtcr > 0.0, (
        f"criticalTimeStep() returned {dtcr!r} for {key}")
    _PEN[key] = float(dtcr)
    return _PEN[key]


def _march_explicit_guard(nsteps, dt, probe_node, probe_dof=2, chunks=8):
    """March the CURRENT explicit model in chunks with finite/blow-up guards
    (cap 1e3 — grotesque vs the mm..cm response scales here)."""
    assert nsteps % chunks == 0, f"nsteps {nsteps} not divisible by {chunks}"
    step = nsteps // chunks
    for c in range(chunks):
        assert ops.analyze(step, dt) == 0, (
            f"explicit march failed in chunk {c} (steps {c * step}..{(c + 1) * step})")
        u = ops.nodeDisp(probe_node, probe_dof)
        assert math.isfinite(u) and abs(u) < 1.0e3, (
            f"explicit march diverged by chunk {c}: u={u}")


# ==========================================================================
# fresh-subprocess machinery (numpy-free children)
# ==========================================================================
_BOOT_CH = '''import os, sys, math
_D = %(dist)r
os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
try: os.add_dll_directory(_D)
except Exception: pass
sys.path.insert(0, _D)
for _m in ("opensees", "openseespy", "openseespy.opensees"): sys.modules.pop(_m, None)
import opensees as ops
assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_D)
'''


def _run_child(body, timeout=240):
    fd, path = tempfile.mkstemp(suffix=".py", dir=_TMP)
    os.close(fd)
    with open(path, "w") as f:
        f.write(body)
    proc = subprocess.run([sys.executable, "-S", "-u", path],
                          capture_output=True, text=True, timeout=timeout,
                          encoding="utf-8", errors="replace")
    return proc.stdout + proc.stderr, proc.returncode


# ---- child: e72 overlay column under CDL at a given dt (CFL sweep legs) ----
_E72_BUILD_CH = '''
def build_e72(fsl_zero=True):
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j in range(9):
        for i in range(2):
            ops.node(2 * j + i + 1, float(i), float(j))
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e7, 0.25, 2.0)
    eles = list(range(1, 9))
    for j in range(8):
        ops.element("quad", j + 1, 2 * j + 1, 2 * j + 2, 2 * j + 4, 2 * j + 3,
                    1.0, "PlaneStrain", 1)
    for j in range(9):
        for i in range(2):
            ops.fix(2 * j + i + 1, 1, 1 if j == 0 else 0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.load(17, 0.0, -5.0e4); ops.load(18, 0.0, -5.0e4)
    args = ["-region"] + eles + ["-Kf", 2.2e9, "-rhoF", 1000.0,
            "-perm", 1.0e-11, 1.0e-11, "-poro", 0.4, "-moduli", 1.0e7, 0.25,
            "-drained", 17, 18, "-stab", "auto", 0.25]
    if fsl_zero:
        args += ["-fsL", "zero"]
    ops.pattern("LadrunoPorousOverlay", 77, *args)
    return eles
'''

_CFL_CHILD = _BOOT_CH + _E72_BUILD_CH + '''
DT = %(dt)r
NST = %(nst)d
build_e72()
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("ProfileSPD")
ops.algorithm("Linear"); ops.integrator("CentralDifferenceLadruno")
ops.analysis("Transient")
mx = 0.0
for s in range(NST):
    try:
        rc = ops.analyze(1, DT)
    except BaseException as ex:
        print("RESULT DIVERGED step", s, "exc", type(ex).__name__); sys.exit(0)
    u = ops.nodeDisp(9, 2)
    if rc != 0 or not math.isfinite(u) or abs(u) > 1.0e3:
        print("RESULT DIVERGED step", s, "rc", rc, "u", u); sys.exit(0)
    if abs(u) > mx: mx = abs(u)
print("RESULT BOUNDED maxu", mx)
sys.exit(0)
'''


def _cfl_child_src(dt, nst):
    return _CFL_CHILD % dict(dist=_DIST, dt=float(dt), nst=int(nst))


# ---- child: ZS84 incumbent (upstream CentralDifference + FourNodeQuadUP) ----
# quadUP arg mapping (test_ladruno_up_element_equiv.py module header):
#   element quadUP tag i j k l thk matTag bulk rhoF permX permY <b1 b2 p>
#   bulk = PRE-COMBINED Qbar ~= Kf/n; p is a RATE DOF -> read via nodeVel(n,3);
#   <b1 b2> would drive BOTH solid and fluid rows (unused here — no gravity).
_INCUMBENT_CHILD = _BOOT_CH + '''
import time
KF = %(kf)r
DTS = %(dts)r
TEND = %(tend)r
TS = %(ts)r
PROBES = %(probes)r
Z_E = 3.0e4; Z_NU = 0.2; Z_RHO = 2.0; Z_RHOF = 1.0; Z_PORO = 0.3
Z_KBAR = 1.0e-2 / (1.0 * 9.81); Z_Q = 1.0; Z_TP = 0.1; NELY = 30; L = 30.0
def build():
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
    for j in range(NELY + 1):
        for i in range(2):
            ops.node(2 * j + i + 1, float(i), L * j / NELY)
    ops.nDMaterial("ElasticIsotropic", 1, Z_E, Z_NU, Z_RHO)
    for j in range(NELY):
        ops.element("quadUP", j + 1, 2 * j + 1, 2 * j + 2, 2 * j + 4, 2 * j + 3,
                    1.0, 1, KF / Z_PORO, Z_RHOF, Z_KBAR, Z_KBAR, 0.0, 0.0, 0.0)
    for j in range(NELY + 1):
        for i in range(2):
            ops.fix(2 * j + i + 1, 1, 1 if j == 0 else 0, 0)
    ops.fix(2 * NELY + 1, 0, 0, 1); ops.fix(2 * NELY + 2, 0, 0, 1)
    ts = [Z_TP * k / 20.0 for k in range(21)] + [1.0e12]
    vs = [math.sin(math.pi * t / (2.0 * Z_TP)) for t in ts[:21]] + [1.0]
    ops.timeSeries("Path", 1, "-time", *ts, "-values", *vs)
    ops.pattern("Plain", 1, 1)
    ops.load(2 * NELY + 1, 0.0, -Z_Q / 2.0, 0.0)
    ops.load(2 * NELY + 2, 0.0, -Z_Q / 2.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.algorithm("Linear"); ops.integrator("CentralDifference")
    ops.analysis("Transient")
for dt in DTS:
    print("TRY dt", dt)
    try:
        build()
    except BaseException as ex:
        print("FAILSTEP -1 buildexc", type(ex).__name__, ex); continue
    n = int(round(TEND / dt))
    if n > 20000: n = 20000
    got = []; ok = True
    t0 = time.perf_counter()
    for s in range(n):
        try:
            rc = ops.analyze(1, dt)
        except BaseException as ex:
            print("FAILSTEP", s, "exc", type(ex).__name__, ex); ok = False; break
        if rc != 0:
            print("FAILSTEP", s, "rc", rc); ok = False; break
        u = ops.nodeDisp(2 * NELY + 1, 2)
        if not math.isfinite(u) or abs(u) > 1.0e3:
            print("FAILSTEP", s, "nonfinite_u", u); ok = False; break
        t = ops.getTime()
        for tt in TS:
            if abs(t - tt) <= dt * 0.51:
                got.append((tt, [ops.nodeVel(p, 3) for p in PROBES]))
    if ok:
        w = (time.perf_counter() - t0) * 1000.0 / max(n, 1)
        print("RAN_AT", dt)
        for tt, vals in got:
            print("PS", tt, *vals)
        print("WALL_MS_PER_STEP", w)
        sys.exit(0)
print("ALLFAIL")
sys.exit(0)
'''


def _incumbent_child_src(kf, dts, tend, ts, probes):
    return _INCUMBENT_CHILD % dict(dist=_DIST, kf=float(kf),
                                   dts=[float(d) for d in dts],
                                   tend=float(tend), ts=[float(t) for t in ts],
                                   probes=[int(p) for p in probes])


# ---- child: Richardson expected-bad (LadrunoUP + upstream CentralDifference) --
_RICHARDSON_CHILD = _BOOT_CH + '''
DT = %(dt)r
NST = %(nst)d
SAMPLES = %(samples)r
PROBES = %(probes)r
Z_E = 3.0e4; Z_NU = 0.2; Z_RHO = 2.0; Z_RHOF = 1.0; Z_PORO = 0.3
Z_KF = 1.0e5; Z_KBAR = 1.0e-2 / (1.0 * 9.81); Z_Q = 1.0; Z_TP = 0.1
NELY = 30; L = 30.0
ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 3)
for j in range(NELY + 1):
    for i in range(2):
        ops.node(2 * j + i + 1, float(i), L * j / NELY)
ops.nDMaterial("ElasticIsotropic", 1, Z_E, Z_NU, Z_RHO)
for j in range(NELY):
    ops.element("LadrunoUP", j + 1, 2 * j + 1, 2 * j + 2, 2 * j + 4, 2 * j + 3,
                1, "-Kf", Z_KF, "-poro", Z_PORO, "-rhoF", Z_RHOF,
                "-perm", Z_KBAR, Z_KBAR, "-stab", "auto", 0.25,
                "-dynSeepage", "off")
for j in range(NELY + 1):
    for i in range(2):
        ops.fix(2 * j + i + 1, 1, 1 if j == 0 else 0, 0)
ops.fix(2 * NELY + 1, 0, 0, 1); ops.fix(2 * NELY + 2, 0, 0, 1)
ts = [Z_TP * k / 20.0 for k in range(21)] + [1.0e12]
vs = [math.sin(math.pi * t / (2.0 * Z_TP)) for t in ts[:21]] + [1.0]
ops.timeSeries("Path", 1, "-time", *ts, "-values", *vs)
ops.pattern("Plain", 1, 1)
ops.load(2 * NELY + 1, 0.0, -Z_Q / 2.0, 0.0)
ops.load(2 * NELY + 2, 0.0, -Z_Q / 2.0, 0.0)
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
ops.algorithm("Linear"); ops.integrator("CentralDifference")
ops.analysis("Transient")
sset = set(SAMPLES)
for s in range(1, NST + 1):
    try:
        rc = ops.analyze(1, DT)
    except BaseException as ex:
        print("SYMPTOM EXC step", s, type(ex).__name__, ex); sys.exit(0)
    if rc != 0:
        print("SYMPTOM ANALYZE_FAIL step", s, "rc", rc); sys.exit(0)
    bad = False
    for p in PROBES:
        v = ops.nodeDisp(p, 3)
        if not math.isfinite(v) or abs(v) > 1.0e6:
            bad = True
    if bad:
        print("SYMPTOM NONFINITE step", s); sys.exit(0)
    if s in sset:
        print("PS", s, *[ops.nodeDisp(p, 3) for p in PROBES])
print("COMPLETED")
sys.exit(0)
'''


def _richardson_child_src(dt, nst, samples, probes):
    return _RICHARDSON_CHILD % dict(dist=_DIST, dt=float(dt), nst=int(nst),
                                    samples=[int(s) for s in samples],
                                    probes=[int(p) for p in probes])


# ---- child: -fsL zero parser advisory (build only; the advisory prints to
#      opserr/stderr, so it must be observed in a child whose ONLY overlay
#      build is the -fsL zero one — stdout/stderr interleaving across streams
#      makes an in-run before/after marker unreliable) ----
_ADVISORY_CHILD = _BOOT_CH + _E72_BUILD_CH + '''
build_e72(fsl_zero=True)
print("BUILD_OK")
sys.exit(0)
'''


def _advisory_child_src():
    return _ADVISORY_CHILD % dict(dist=_DIST)


# ---- child: -fsL zero driver refusal ----
_REFUSAL_CHILD = _BOOT_CH + _E72_BUILD_CH + '''
build_e72(fsl_zero=True)
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("ProfileSPD")
ops.algorithm("Linear"); ops.integrator("Newmark", 0.6, 0.3025)
ops.analysis("Transient")
try:
    rc = ops.LadrunoStaggeredAnalyze(1, 1.0e-4, "-tol", 1.0e-6)
    if rc < 0:
        print("FATAL_RC", rc); sys.exit(3)
    print("DRIVER_ACCEPTED", rc); sys.exit(5)
except SystemExit:
    raise
except BaseException as ex:
    print("FATAL_EXC", type(ex).__name__, ex); sys.exit(7)
'''


def _refusal_child_src():
    return _REFUSAL_CHILD % dict(dist=_DIST)


def _parse_ps_rows(out, key_list):
    """Collect 'PS <key> v1 v2 ...' rows; return array ordered by key_list
    (None if any key missing)."""
    rows = {}
    for ln in out.splitlines():
        parts = ln.split()
        if len(parts) >= 2 and parts[0] == "PS":
            try:
                key = float(parts[1])
                rows[key] = [float(v) for v in parts[2:]]
            except ValueError:
                continue
    got = []
    for k in key_list:
        hit = None
        for kk in rows:
            if abs(kk - float(k)) <= 1e-9 * max(abs(float(k)), 1.0):
                hit = rows[kk]
        if hit is None:
            return None
        got.append(hit)
    return np.array(got)


# ==========================================================================
# (a) TWO-LEG on ZS84: CD + overlay vs implicit Newmark LadrunoUP
# ==========================================================================
_ZS84_PROBES_J = (15, 24, 27, 29)     # mid-column + the active drainage front
                                       # (front depth sqrt(Tv)L ~ 6.8 m at t_end)


def test_a_two_leg_zs84():
    """Mutual dt-convergence, 3 levels {0.4, 0.2, 0.1}x the overlay-aware
    advisory: relL2 diffs on p probes (8 matched sample times) shrink with
    observed order >= 0.9 (P1/P2 methodology — exact equality is NOT expected
    across an O(dt) split + two different one-step integrators)."""
    adv = _pencil("zs84_ov")
    t_end = 1.5
    n0 = int(8 * math.ceil(t_end / (0.4 * adv) / 8))
    if n0 > 2400:                      # runtime guard: shrink the horizon, keep dt
        t_end = t_end * 2400.0 / n0
        n0 = 2400
    probes = [_nt(0, j) for j in _ZS84_PROBES_J]
    diffs = []
    for lev in range(3):
        n = n0 * 2 ** lev
        dt = t_end / n
        rec = _rec(f"a_exp{lev}.csv")
        _build_zs84_explicit(record=rec)
        _march_explicit_guard(n, dt, _nt(0, _Z["nely"] - 1))
        _t, P, t2c = _read_record(rec)
        assert len(P) == n, f"level {lev}: {len(P)} record rows != {n} commits"
        rows = [k * n // 8 - 1 for k in range(1, 9)]
        cols = [t2c[p] for p in probes]
        Pe = P[np.ix_(rows, cols)]
        assert np.all(np.isfinite(Pe)), f"level {lev}: non-finite explicit p"
        samples = [k * n // 8 for k in range(1, 9)]
        Pi = _march_up_samples(dt, n, samples, probes)
        diffs.append(_relL2(Pe, Pi))
    orders = [math.log2(diffs[i] / diffs[i + 1]) for i in range(len(diffs) - 1)]
    assert diffs[-1] < diffs[0], f"two-leg diffs did not shrink: {diffs}"
    assert min(orders) >= 0.9, (
        f"two-leg mutual order < 0.9: orders={[round(o, 2) for o in orders]} "
        f"diffs={[f'{d:.3e}' for d in diffs]}")
    info = (f"advisory {adv:.4e}s dt0 {t_end / n0:.4e}s diffs "
            f"{[f'{d:.2e}' for d in diffs]} orders {[round(o, 2) for o in orders]}")
    print(f"PASS (a) two-leg ZS84: {info}")
    return info


# ==========================================================================
# (b) CFL envelope + naive-drained-advisory demonstration
# ==========================================================================
def _div_step(out):
    """Parse the divergence step from a CFL child's 'RESULT DIVERGED step N'."""
    for ln in out.splitlines():
        if ln.startswith("RESULT DIVERGED"):
            try:
                return int(ln.split()[3])
            except (IndexError, ValueError):
                return -1
    return -1


def test_b_cfl_envelope():
    """FIXED-HORIZON envelope. Run-2 adjudication measurement (2026-07-18,
    this worktree, 1x8 zero-drainage column, marches to 80k steps): the
    zero-drainage (kbar = 1e-11) UNDAMPED worst case has NO asymptotically
    stable dt — a slow SECULAR ENERGY-PUMPING of the staggered coupling
    (growth per unit time ~ O(dt)) eventually diverges every march
    [steps-to-|u|>1e3: 0.4x->48344, 0.5x->30018, 0.6x->21008, 0.8x->11789,
    1.0x->7761, 1.5x->3372, 2.0x->1904, 3.0x->855], while the SAME solid
    WITHOUT the overlay is clean for 30k steps (drained control BOUNDED).
    Run-1's '1.1x/1.4x bounded' and run-2's '3.0x bounded' were 400/600-step
    detection-window artifacts of that slow growth, NOT an over-conservative
    advisory. Consequences gated here:
      * SAFE SIDE: {0.6, 0.8, 0.95}x the advisory survive the 600-step
        envelope (the frozen 0.5x default buys ~4x the blow-up horizon of
        1.0x and ~35x that of 3.0x — margin, not absolute stability, on the
        zero-drainage pathology; ANY physical drainage or damping absorbs
        the O(dt) pumping — the ZS84 gate-(a) leg marches 6.5k steps clean).
      * UPWARD: {1.5, 2.0, 3.0}x must all diverge within a 4500-step cap
        (measured 3372/1904/855) — divergence occurs at or below 3.0x as
        demanded; the measured steps are printed as the growth trend.
      * 0.9x the DRAINED pencil (~23x the advisory) diverges fast (~25
        steps) — the naive-advisory demonstration.
    MAIN: the secular-pumping characterization + drained-control fact is the
    ADR §12 P3 / quirks-row payload; the per-element advisory bounds the
    FAST instability and remains the correct advisory pencil."""
    und = _pencil("e72_ov")
    dr = _pencil("e72_dr")
    assert dr > und, f"drained pencil {dr:.3e} not > undrained {und:.3e}"
    lines = []
    # safe side: at or below the advisory (x margin) survives the envelope
    for mult in (0.6, 0.8, 0.95):
        out, rc = _run_child(_cfl_child_src(mult * und, 600))
        bounded = (rc == 0) and ("RESULT BOUNDED" in out)
        assert bounded, (
            f"{mult}x the advisory should survive the 600-step envelope "
            f"(secular pumping is far slower there — 0.8x measured 11789 "
            f"steps to blow up):\n" + out[-1200:])
        lines.append(f"{mult}x=B@600")
    # upward legs: all must diverge within the 4500-step cap
    first_div = None
    for mult in (1.5, 2.0, 3.0):
        out, rc = _run_child(_cfl_child_src(mult * und, 4500))
        diverged = (rc != 0) or ("RESULT DIVERGED" in out)
        assert diverged, (
            f"{mult}x the advisory stayed BOUNDED for 4500 steps (measured "
            f"3372/1904/855 for 1.5/2.0/3.0x on run 2's adjudication probe — "
            f"a regression in the divergence character):\n" + out[-1200:])
        step = _div_step(out)
        lines.append(f"{mult}x=D@{step}")
        if first_div is None:
            first_div = mult
    assert first_div is not None and first_div <= 3.0, (
        f"no divergence at or below 3.0x: {' '.join(lines)}")
    # naive-advisory demonstration: 0.9x DRAINED pencil on the COUPLED model
    out, rc = _run_child(_cfl_child_src(0.9 * dr, 400))
    assert (rc != 0) or ("RESULT DIVERGED" in out), (
        "0.9x the DRAINED pencil (= "
        f"{0.9 * dr / und:.1f}x the undrained advisory) did NOT diverge — the "
        "naive drained advisory would have to be unsafe for this gate to hold:\n"
        + out[-1200:])
    lines.append(f"0.9xDRAINED({0.9 * dr / und:.1f}xund)=D@{_div_step(out)}")
    print("=" * 74)
    print("LOUD (b) SECULAR-PUMPING RECORD (for ADR §12 P3 / quirks row):")
    print("   zero-drainage undamped column: NO asymptotically stable dt —")
    print("   steps-to-blowup 0.4x:48344 0.5x:30018 0.8x:11789 1.0x:7761")
    print("   1.5x:3372 2.0x:1904 3.0x:855 (growth/unit-time ~ O(dt));")
    print("   drained control (no overlay) BOUNDED 30k steps; any physical")
    print("   drainage/damping absorbs the pumping (gate (a) ZS84 leg clean).")
    print("=" * 74)
    info = (f"pencils und {und:.4e}s dr {dr:.4e}s; fixed-horizon envelope: "
            f"safe side bounded @600, divergence from {first_div}x within "
            f"4500 (at/below 3.0x as demanded); sweep {' '.join(lines)}")
    print(f"PASS (b) CFL envelope: {info}")
    return info


# ==========================================================================
# (c) advisory ratio vs the material formula (direction + magnitude band)
# ==========================================================================
def test_c_advisory_ratio():
    """MERGED single gate (final MAIN adjudication, 2026-07-18).
    (i) HARD STAB-INVARIANCE IDENTITY: |ratio_stabOn - ratio_stabOff| /
        ratio <= 1e-12 — the measured 1.29e-14 mechanism pin: the per-cell
        constant-p mode governs the cell-local condensation and lies in the
        H-tilde NULLSPACE (the stab Laplacian is singular on constants), so
        the advisory pencil is stab-invariant by construction.
    (ii) SINGLE ratio band: drained/undrained within [0.5, 1.3]x the
        material formula sqrt(1 + Kf/(n*Moed)) = 21.43 (measured 26.24 =
        1.2243x material).
    (iii) DIAGNOSIS (kept per adjudication): the PER-ELEMENT factor is
        mesh/element-dependent and legitimately exceeds the material bound —
        the element-eigenvalue route and the cell-local condensation both
        stack conservative. The toy's discrete 15.36 (= 0.72x material,
        E7.4) was measured on the ASSEMBLED global operator: the
        26.24-vs-15.36 gap is OPERATOR CLASS (per-element vs assembled),
        not stab. All numbers printed for the PR."""
    und = _pencil("e72_ov")
    und_soff = _pencil("e72_ov_soff")
    dr = _pencil("e72_dr")
    ratio = dr / und
    ratio_soff = dr / und_soff
    fac = math.sqrt(1.0 + _E72["Kf"] / (_E72["poro"] * _moed(_E72["E"], _E72["nu"])))
    delta = abs(ratio - ratio_soff) / ratio
    print(f"(c) MEASURED: ratio {ratio:.4f} | material {fac:.4f} | "
          f"ratio/material {ratio / fac:.4f} | stab-invariance delta "
          f"{delta:.2e} | per-element vs assembled-operator: toy 15.36 = "
          f"0.72x material was the ASSEMBLED operator — the gap is operator "
          f"class, not stab")
    # (i) hard stab-invariance identity
    assert delta <= 1e-12, (
        f"STAB-INVARIANCE BROKEN: |ratio_on - ratio_off|/ratio = {delta:.3e} "
        f"> 1e-12 (mechanism pin: constant-p mode governs the condensation "
        f"and sits in the H-tilde nullspace — a nonzero delta means the "
        f"augmentation started seeing the stab operator)")
    # (ii) single band vs the material formula
    assert ratio > 1.0, f"advisory not conservative-direction: ratio {ratio:.3f}"
    assert 0.5 * fac <= ratio <= 1.3 * fac, (
        f"drained/undrained ratio {ratio:.4f} = {ratio / fac:.4f}x material "
        f"outside [0.5, 1.3]x {fac:.2f} (measured 1.2243x at adjudication; "
        f"per-element bound legitimately exceeds material — assembled-"
        f"operator toy 15.36 = 0.72x is a different operator class)")
    info = (f"ratio {ratio:.2f} ({ratio / fac:.4f}x material {fac:.2f}, band "
            f"[0.5, 1.3]x), stab-invariance delta {delta:.1e} (<= 1e-12), "
            f"toy assembled-operator 15.36 = 0.72x (operator class, not "
            f"stab) [und {und:.4e}s dr {dr:.4e}s]")
    print(f"PASS (c) advisory: {info}")
    return info


# ==========================================================================
# (d) incumbent head-to-head + S->0 demo (accuracy/wall-clock REPORT ONLY;
#     hard gates: overlay legs finite)
# ==========================================================================
def test_d_incumbent_head_to_head():
    adv = _pencil("zs84_ov")
    dt = 0.4 * adv
    n = 800
    samples = [200, 400, 600, 800]
    probes = [_nt(0, j) for j in (24, 27, 29)]
    ref = _march_up_samples(dt, n, samples, probes)

    # ---- overlay leg (timed) ----
    rec = _rec("d_ov.csv")
    _build_zs84_explicit(record=rec)
    t0 = time.perf_counter()
    _march_explicit_guard(n, dt, _nt(0, _Z["nely"] - 1))
    wall_ov = (time.perf_counter() - t0) * 1000.0 / n
    _t, P, t2c = _read_record(rec)
    Pe = P[np.ix_([s - 1 for s in samples], [t2c[p] for p in probes])]
    assert np.all(np.isfinite(Pe)), "overlay leg produced non-finite p"
    rel_ov = _relL2(Pe, ref)

    # ---- incumbent leg (fresh subprocess, retry ladder dt, dt/4, dt/16) ----
    ts = [s * dt for s in samples]
    out, rc = _run_child(_incumbent_child_src(_Z["Kf"], [dt, dt / 4, dt / 16],
                                              n * dt, ts, probes))
    report = [f"overlay: rel-vs-implicit {rel_ov:.3e}, {wall_ov:.3f} ms/step"]
    if "RAN_AT" in out and rc == 0:
        dt_used = None
        wall_inc = float("nan")
        for ln in out.splitlines():
            if ln.startswith("RAN_AT"):
                dt_used = float(ln.split()[1])
            if ln.startswith("WALL_MS_PER_STEP"):
                wall_inc = float(ln.split()[1])
        arr = _parse_ps_rows(out, ts)
        if arr is not None and np.all(np.isfinite(arr)):
            rel_inc = _relL2(arr, ref)
            note = "" if dt_used == dt else (
                f" [NB: quadUP needed dt={dt_used:.3e} = "
                f"{dt_used / dt:.3g}x the overlay dt to run]")
            report.append(f"incumbent quadUP+CD: RAN at dt {dt_used:.3e}, "
                          f"rel-vs-implicit {rel_inc:.3e}, "
                          f"{wall_inc:.3f} ms/step{note}")
        else:
            report.append("incumbent quadUP+CD: ran but p samples missing/"
                          "non-finite — SYMPTOM: garbage p output")
    else:
        tail = " | ".join(ln for ln in out.splitlines()
                          if ln.startswith(("FAILSTEP", "ALLFAIL")))[:400]
        report.append(f"incumbent quadUP+CD: DID NOT COMPLETE at dt, dt/4, "
                      f"dt/16 — SYMPTOM: {tail or out[-300:]}")

    # ---- S->0 demo: Kf x 1e3 ----
    adv2 = _pencil("zs84_ov_kf1e3")
    assert adv2 < adv, (
        f"Kf x1e3 advisory {adv2:.3e} not < base {adv:.3e} — the augmentation "
        "did not stiffen with Kf")
    dt2 = 0.4 * adv2
    _build_zs84_explicit(Kf=_Z["Kf"] * 1.0e3)
    _march_explicit_guard(800, dt2, _nt(0, _Z["nely"] - 1))     # HARD: finite
    report.append(f"S->0 overlay leg (Kf x1e3): FINITE at 0.4x its advisory "
                  f"({adv2:.3e}s)")
    out2, rc2 = _run_child(_incumbent_child_src(_Z["Kf"] * 1.0e3,
                                                [dt2, dt2 / 4], 800 * dt2,
                                                [800 * dt2], probes))
    if "RAN_AT" in out2 and rc2 == 0:
        arr2 = _parse_ps_rows(out2, [800 * dt2])
        fin = arr2 is not None and np.all(np.isfinite(arr2))
        report.append("S->0 INCUMBENT SYMPTOM: quadUP+CD RAN and stayed "
                      + ("finite (no failure measured at this Kf — reported "
                         "honestly)" if fin else "NON-FINITE (garbage p)"))
    else:
        tail2 = " | ".join(ln for ln in out2.splitlines()
                           if ln.startswith(("FAILSTEP", "ALLFAIL")))[:400]
        report.append(f"S->0 INCUMBENT SYMPTOM: {tail2 or out2[-300:]}")

    print("PASS (d) incumbent head-to-head (accuracy/wall-clock report only):")
    for r in report:
        print("   " + r)
    return "; ".join(report)


# ==========================================================================
# (e) Richardson EXPECTED-BAD pin (fresh subprocess)
# ==========================================================================
def test_e_richardson_expected_bad():
    """honest-p LadrunoUP + upstream CentralDifference at dt well below the
    solid pencil (~8e-3 s). Run-1 REFUTATION CANDIDATE: the 400-step
    dt=2e-4 leg TRACKED the implicit reference (rel 2.33e-4) — so this gate
    re-measures harder before pinning either way (run-1 adjudication): TWO
    legs, 4000 steps each, dt = 2e-4 AND dt = 1e-3 (both << the solid
    pencil), fresh children.
      * If EITHER leg degrades (non-finite / refuses to run / rel > 0.5):
        the EXPECTED-BAD pin holds on the longer march — pass, reporting
        which leg manifests.
      * If BOTH legs still track (rel <= 0.05): REFUTATION-RECORDING mode —
        assert the tracking and print loudly that the Richardson-instability
        claim is NOT reproduced; MAIN rewrites the quirks row from this
        report.
      * Anything in between (finite but 0.05 < rel <= 0.5) is ambiguous —
        hard FAIL with the measured numbers for MAIN to adjudicate."""
    probes = [_nt(0, j) for j in (24, 27, 29)]
    legs = [(2.0e-4, 4000), (1.0e-3, 4000)]
    outcomes = []                     # (dt, kind, detail, rel-or-None)
    for dt, n in legs:
        samples = [n // 4, n // 2, 3 * n // 4, n]
        ref = _march_up_samples(dt, n, samples, probes)
        assert np.linalg.norm(ref) > 0, "test bug: implicit reference p is zero"
        out, rc = _run_child(_richardson_child_src(dt, n, samples, probes))
        if rc != 0:
            outcomes.append((dt, "degraded", f"child died (rc={rc}): "
                             f"{out[-200:]}", None))
        elif "SYMPTOM" in out:
            ln = next(x for x in out.splitlines() if x.startswith("SYMPTOM"))
            outcomes.append((dt, "degraded", ln, None))
        elif "COMPLETED" in out:
            arr = _parse_ps_rows(out, samples)
            assert arr is not None, (
                f"dt={dt}: completed but samples unparseable:\n{out[-600:]}")
            if not np.all(np.isfinite(arr)):
                outcomes.append((dt, "degraded", "NON-FINITE p samples", None))
            else:
                rel = _relL2(arr, ref)
                kind = "degraded" if rel > 0.5 else (
                    "tracks" if rel <= 0.05 else "ambiguous")
                outcomes.append((dt, kind, f"rel-vs-implicit {rel:.3e}", rel))
        else:
            outcomes.append((dt, "degraded",
                             f"no recognizable output: {out[-200:]}", None))

    detail = "; ".join(f"dt={dt:g}: {k} ({d})" for dt, k, d, _r in outcomes)
    degraded = [o for o in outcomes if o[1] == "degraded"]
    if degraded:
        which = ", ".join(f"dt={o[0]:g}" for o in degraded)
        print(f"PASS (e) Richardson EXPECTED-BAD pin holds on the longer "
              f"march — manifests at {which}. Measured: {detail}")
        return f"EXPECTED-BAD holds ({which}); {detail}"
    assert all(o[1] == "tracks" for o in outcomes), (
        f"(e) AMBIGUOUS ZONE — neither expected-bad (rel > 0.5) nor clean "
        f"tracking (rel <= 0.05); MAIN must adjudicate. Measured: {detail}")
    print("=" * 74)
    print("LOUD (e): Richardson-instability claim NOT REPRODUCED "
          "(4000 steps, dt 2e-4 AND 1e-3, ZS84 honest-p LadrunoUP + upstream "
          "CentralDifference): quirks row must record the REFUTATION.")
    print(f"   measured: {detail}")
    print("=" * 74)
    return f"REFUTATION RECORDED (both legs track); {detail}"


# ==========================================================================
# (f) -subcycle sweep + auto-N inference
# ==========================================================================
def test_f_subcycle_sweep():
    """Fixed dt = 0.5x the pencil; permeability re-tuned (kbar_f) so the
    subcycle window is dynamically meaningful (target 0.089*h^2/(cv*dt) = 4.5
    -> N_auto = 4; the pencil is storage/stiffness-based and does NOT depend
    on kbar, so the frozen-e72 pencil measurement transfers). N in
    {1, 2, 5, 10} all bounded; err(N) vs N=1 grows ~monotonically with
    log-log slope in [0.7, 1.7] (E7.3a ~N^1.2). -subcycle auto: N inferred
    from record rows/commits within +-1 of the theta hand-count."""
    und = _pencil("e72_ov")
    dt = 0.5 * und
    Moed = _moed(_E72["E"], _E72["nu"])
    kbar_f = _THETA / (4.5 * Moed * dt)          # => 0.089*h^2/(cv*dt) = 4.5
    ncom = 40
    sample_times = [10 * dt * k for k in range(1, 5)]    # commits 10,20,30,40
    guard = _nt(0, _E72["nely"] // 2)

    traces = {}
    for N in (1, 2, 5, 10):
        rp = _rec(f"f_N{N}.csv")
        _build_e72_explicit(kbar=kbar_f, subcycle=N, record=rp)
        for s in range(ncom):
            assert ops.analyze(1, dt) == 0, f"N={N} march failed at commit {s}"
            u = ops.nodeDisp(guard, 2)
            assert math.isfinite(u) and abs(u) < 1.0e3, f"N={N} diverged: u={u}"
        tt, P, _t2c = _read_record(rp)
        rows = []
        for st in sample_times:
            hit = np.nonzero(np.abs(tt - st) < 0.5 * dt)[0]
            assert len(hit) == 1, (
                f"N={N}: no unique record row at t={st:.3e} (rows at "
                f"{list(np.round(tt, 8))[:12]}...)")
            rows.append(int(hit[0]))
        traces[N] = P[rows, :]
        assert np.all(np.isfinite(traces[N])), f"N={N}: non-finite p"

    err = {N: _relL2(traces[N], traces[1]) for N in (2, 5, 10)}
    assert err[10] > err[2] > 0, f"err(N) not growing: {err}"
    assert err[5] >= 0.8 * err[2] and err[10] >= 0.8 * err[5], (
        f"err(N) not ~monotone: {err}")
    slope = _loglog_slope([2, 5, 10], [err[2], err[5], err[10]])
    assert 0.7 <= slope <= 1.7, (
        f"err(N) log-log slope {slope:.2f} outside [0.7, 1.7] (E7.3a ~N^1.2); "
        f"err={err}")

    # -subcycle auto: infer resolved N from record rows vs commits
    ra = _rec("f_auto.csv")
    _build_e72_explicit(kbar=kbar_f, subcycle="auto", record=ra)
    for s in range(ncom):
        assert ops.analyze(1, dt) == 0, f"auto march failed at commit {s}"
    _ta, Pa, _ = _read_record(ra)
    assert len(Pa) >= 1, "auto run recorded no rows"
    N_inf = int(round(ncom / len(Pa)))
    N_hand = max(1, math.floor(_THETA * 1.0 / ((kbar_f * Moed) * dt)))
    assert abs(N_inf - N_hand) <= 1, (
        f"-subcycle auto resolved N ~ {N_inf} (from {len(Pa)} rows / {ncom} "
        f"commits) vs theta hand-count {N_hand} (theta={_THETA}, h=1, "
        f"cv={kbar_f * Moed:.3e}, dt={dt:.3e}) — outside +-1")
    info = (f"err(N) 2:{err[2]:.3e} 5:{err[5]:.3e} 10:{err[10]:.3e} "
            f"slope {slope:.2f} (band 0.7-1.7); auto N~{N_inf} vs hand "
            f"{N_hand} (+-1)")
    print(f"PASS (f) subcycle: {info}")
    return info


# ==========================================================================
# (g) energy closure on the (a)-config overlay run (ADR-69)
# ==========================================================================
def test_g_energy_closure():
    """EnergyBalanceRecorder (legacy layout: time KE IE DW ULW RES ERR%) on the
    ZS84 CD+overlay march at 0.4x the advisory: |ERR%| <= 2.0 at end-of-march
    and ULW nonzero — the overlay's +Q*p work rides the ULW external-work
    channel (plan §3.3: 'documented, not silently absent'). Recorder attached
    AFTER the integrator (the ADR-69 ordering rule)."""
    adv = _pencil("zs84_ov")
    t_end = 1.5
    n = int(8 * math.ceil(t_end / (0.4 * adv) / 8))
    if n > 2400:
        t_end = t_end * 2400.0 / n
        n = 2400
    dt = t_end / n
    ef = _rec("g_energy.txt")
    _build_zs84_explicit()                       # sets the integrator...
    ops.recorder("EnergyBalance", "-file", ef, "-time")   # ...then the recorder
    _march_explicit_guard(n, dt, _nt(0, _Z["nely"] - 1))
    ops.wipe()                                   # flush/close the recorder
    d = np.atleast_2d(np.loadtxt(ef))
    assert d.shape[1] >= 7, f"expected >= 7 energy columns, got {d.shape[1]}"
    tail = d[-max(1, len(d) // 4):]
    assert np.all(np.isfinite(tail)), "non-finite energy ledger in the last quarter"
    ke, ie, ulw, err_end = d[-1, 1], d[-1, 2], d[-1, 4], abs(d[-1, 6])
    assert abs(ulw) > 0.0, (
        "ULW identically zero at end-of-march — external (and +Q*p coupling) "
        "work did not ride the ULW channel")
    assert err_end <= 2.0, (
        f"energy closure |ERR%| = {err_end:.3f} > 2.0 at end-of-march — the "
        f"overlay's +Q*p work is NOT closing through ULW (plan §3.3)")
    info = (f"end-of-march KE {ke:.3e} IE {ie:.3e} ULW {ulw:.3e} "
            f"|ERR%| {err_end:.3f} (<= 2.0)")
    print(f"PASS (g) energy closure: {info}")
    return info


# ==========================================================================
# (h) surface: LadrunoStaggeredAnalyze refuses -fsL zero + parser advisory
# ==========================================================================
def test_h_driver_refuses_fsl_zero():
    """Fresh subprocesses. (1) advisory child: BUILDING an -fsL zero overlay
    (and doing nothing else) must print the one-time loud parser advisory
    (loose text match — 3.B owns the wording; MAIN reconciles). (2) refusal
    child: the P2 driver on that overlay must fatal LOUDLY (negative rc ->
    exit 3, or exception -> exit 7 — the P2 e5 dual handling); quiet success
    (exit 5) fails."""
    # (1) parser advisory at build time
    out_a, rc_a = _run_child(_advisory_child_src())
    assert rc_a == 0 and "BUILD_OK" in out_a, (
        f"-fsL zero build child failed outright:\n{out_a[-1200:]}")
    low_a = out_a.lower()
    assert ("zero" in low_a) or ("explicit" in low_a) or ("l = 0" in low_a) \
        or ("l=0" in low_a), (
        "-fsL zero parser advisory not observed at build time (the ONE-TIME "
        f"loud advisory is part of the frozen 3.1 surface):\n{out_a[-1200:]}")

    # (2) driver refusal
    out, rc = _run_child(_refusal_child_src())
    assert rc != 0, f"-fsL zero driver child exited 0 (no fatal):\n{out[-1200:]}"
    assert "DRIVER_ACCEPTED" not in out, (
        "LadrunoStaggeredAnalyze RAN an -fsL zero overlay (must refuse — "
        f"iterating with L=0 IS the drained split, plan §3.1):\n{out[-1200:]}")
    low = out.lower()
    assert ("zero" in low) or ("fsl" in low) or ("explicit" in low), (
        f"refusal message does not cite zero/fsL/explicit:\n{out[-1200:]}")
    mode = "rc<0" if rc == 3 else ("exception" if rc == 7 else f"exit {rc}")
    print(f"PASS (h) driver refusal: -fsL zero overlay fatals loudly ({mode}); "
          "parser advisory present at build")
    return f"refusal via {mode}"


# ==========================================================================
# standalone runner (table-driven, the P2 pattern)
# ==========================================================================
_ALL = [
    ("a  two-leg ZS84 (CD+overlay vs UP)", test_a_two_leg_zs84),
    ("b  CFL envelope + naive-drained demo", test_b_cfl_envelope),
    ("c  advisory ratio vs material formula", test_c_advisory_ratio),
    ("d  incumbent head-to-head + S->0", test_d_incumbent_head_to_head),
    ("e  Richardson EXPECTED-BAD pin", test_e_richardson_expected_bad),
    ("f  -subcycle sweep + auto-N", test_f_subcycle_sweep),
    ("g  energy closure (ADR-69 ULW)", test_g_energy_closure),
    ("h  driver refuses -fsL zero", test_h_driver_refuses_fsl_zero),
]


def main():
    if not _BUILT:
        print(f"SKIP: worktree engine not built at {_DIST}")
        return 0
    import traceback
    fails = []
    for name, fn in _ALL:
        try:
            info = fn()
            print(f"[PASS ] {name}" + (f"  -- {info}" if info else ""))
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL ] {name}: {exc}")
            traceback.print_exc()
            fails.append(name)
    print("=" * 74)
    if fails:
        print(f"FAILED {len(fails)}/{len(_ALL)} gates: {fails}")
        return 1
    print(f"ALL {len(_ALL)} EXPLICIT-LANE GATES PASSED (ADR-73 P3, -fsL zero)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
