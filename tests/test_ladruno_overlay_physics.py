"""LadrunoPorousOverlay (PATTERN 33022) — ADR-73 P1 PHYSICS battery (WP1.E-i).

The overlay is a persistent pore-fluid pressure field with its own life-cycle,
coupled to plain ndf=ndm solid elements through the effective-stress seam
(force = +Q*p on the skeleton, compression-positive p). It advances the fluid
with the FIXED-STRESS RATE FORM (ADR §12 P0 pin 1):

    (S* + L + dt*H) p_trial = (S* + L) p0 + L*dp_ref - Q^T du + dt*f_seep

with dp_ref = the previous COMMITTED increment on the first advance of a step
(fs1, O(dt), single pass). This file is the PHYSICS half of the P1 battery.

WHY THE OVERLAY CAN BE CHECKED BIT-FOR-BIT AGAINST THE TOY:
  Driven under a Transient Newmark analysis with a MASSLESS solid (rho = 0),
  each step is: solid solve K u = f + Q*p_committed  ->  commit -> overlay
  advances the fluid with du. That is EXACTLY the frozen numpy toy's
  qs_staggered("fs1_extrap") flow (same LadrunoUPKernel operators, same McGann
  stab alpha0=0.25, same rate form). Measured: overlay == toy fs1_extrap to
  ~3e-7 (CG rel-tol 1e-10), and overlay-vs-monolithic reproduces the frozen
  e7_summary.json e71 splitting error to 1.00x. The toy (imported as a FROZEN
  library, never modified) supplies both the reference solutions and the
  measured error envelopes the gates are pinned to.

THE SIX GATES (plan WP1.E-i contract):
  1. Terzaghi vs analytic + toy anchor: overlay == toy fs1_extrap; isochrones
     vs the Terzaghi series within the mesh floor; late-window splitting error
     <= 2x the frozen-toy fs1_extrap|classic value at matched dt; dt-halving
     order >= 0.9.
  2. Two-leg vs monolithic LadrunoUP (33017): mutual dt-convergence, order >= 1.
  3. Sequencing contract (kills the fluid-first p==0 trap): first committed p
     at an undrained interior node ~ the undrained response, p(0+)/q in
     [0.30,1.0], hard-fail < 0.05.
  4. Rate-form-vs-fixed-point contract: overlay error is in the fs1_extrap
     (rate-form) class -- <= 3x the toy extrap error and <= 0.2x the toy
     fixed-point error vs monolithic-UP -- and matches toy fs1_extrap while
     diverging from toy fs1_plain (fixed-point) -- proving the C++ advance is
     the rate form.
  5. E4-in-OpenSees removal: real `remove element` crack, -onRemoval keep vs
     drain: keep drains THROUGH the crack (no trapped plateau), drain kFactor
     100 strictly faster; fixed-window jump finite + march stable; removal
     forces a fluid sync under -subcycle N>1.
  6. -fsL manual floor loud warning: -fsL 0.1 clamps to the oedometric floor
     with a printed WARNING; the run still completes.
  + drained region nodes stay identically 0 in every overlay run.

RUNNER (worktree dist/bin opensees.pyd is CPython 3.12):
  python -S tests/test_ladruno_overlay_physics.py        # standalone, exits 1 on fail
  py -3.12 -m pytest tests/test_ladruno_overlay_physics.py -x -q   # pytest
The `-S` invocation re-adds site-packages below (numpy) and evicts any boot-.pth
preloaded stale opensees.pyd so THIS worktree's engine wins.
"""
import os
import sys
import types
import tempfile
from pathlib import Path

# --------------------------------------------------------------------------
# bootstrap
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

try:
    import pytest  # noqa: E402
    _HAVE_PYTEST = True
except Exception:  # pragma: no cover
    _HAVE_PYTEST = False

# (b) the frozen ADR-73 P0 toy imports scipy.linalg at top; the functions we
#     call (column_setup / qs_staggered / qs_mono / derived) never touch it, so
#     stub scipy.linalg if scipy is absent (do NOT install / modify the toy).
if "scipy" not in sys.modules:
    try:
        import scipy.linalg  # noqa: F401
    except Exception:
        _sc = types.ModuleType("scipy")
        _scl = types.ModuleType("scipy.linalg")
        _sc.linalg = _scl
        sys.modules["scipy"] = _sc
        sys.modules["scipy.linalg"] = _scl

# (c) THIS worktree's engine (evict any boot-.pth preloaded stale pyd).
_ROOT = Path(__file__).resolve().parents[1]
_DIST = str(_ROOT / "dist" / "bin")
_BUILT = os.path.isfile(os.path.join(_DIST, "opensees.pyd"))
if _BUILT:
    os.environ["PATH"] = _DIST + os.pathsep + os.environ.get("PATH", "")
    try:
        os.add_dll_directory(_DIST)
    except (FileNotFoundError, OSError):
        pass
    for _m in ("opensees", "openseespy", "openseespy.opensees"):
        sys.modules.pop(_m, None)
    if _DIST not in sys.path:
        sys.path.insert(0, _DIST)

if _HAVE_PYTEST and not _BUILT:
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST})"
)

# (d) the frozen toy (ADR-71 spike, extended by ADR-73 P0) as a library.
_TOYDIR = str(_ROOT / "Ladruno_implementation" / "adr71_meshless_p_spike")
if _TOYDIR not in sys.path:
    sys.path.insert(0, _TOYDIR)
import staggered_pins_e7 as e7  # noqa: E402

if _HAVE_PYTEST:
    pytestmark = [pytest.mark.zone_b, pytest.mark.t1]

_TMP = tempfile.mkdtemp(prefix="ov_phys_")


def _rec(name):
    return os.path.join(_TMP, name)


# ==========================================================================
# Toy E7.1 consolidation column (the reference geometry / material)
#   E=1e7, nu=0.25, Kf=2.2e9, poro=0.4, kbar=1e-11, q=1e5, H=10, 2x20 Q4.
# ==========================================================================
_E, _NU, _KF, _PORO, _KH, _GW = 1.0e7, 0.25, 2.2e9, 0.4, 1.0e-7, 1.0e4
_KBAR = _KH / _GW
_NX, _NY, _LX, _LY, _Q = 2, 20, 1.0, 10.0, 1.0e5


def _toy_column():
    """Frozen-toy column_setup (numpy operators + Terzaghi geometry)."""
    return e7.column_setup(_E, _NU, _KF, _PORO, _KH, _GW, stab=True,
                           nx=_NX, ny=_NY, Lx=_LX, Ly=_LY, q=_Q)


def _cv():
    return _toy_column()["mat"]["cv"]


def _relL2(A, B):
    return float(np.linalg.norm(A - B) / np.linalg.norm(B))


def _read_record(path):
    """Return (times, P[nrec,Nregion], tag->col map)."""
    with open(path) as fh:
        lines = fh.read().strip().splitlines()
    hdr = [int(t) for t in lines[0].split(",")[1:]]
    data = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    return data[:, 0], data[:, 1:], {t: k for k, t in enumerate(hdr)}


def _build_overlay_column(rho, nodes, elems, top, rec, **ov):
    """Plain ndf=2 Q4 column + LadrunoPorousOverlay over every cell.

    rho=0 makes the solid quasi-static per step (matches the toy exactly)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for i, (x, y) in enumerate(nodes):
        ops.node(i + 1, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, rho)
    for e, c in enumerate(elems):
        ops.element("quad", e + 1, int(c[0] + 1), int(c[1] + 1),
                    int(c[2] + 1), int(c[3] + 1), 1.0, "PlaneStrain", 1)
    for i, (x, y) in enumerate(nodes):
        fx = 1 if (abs(x) < 1e-9 or abs(x - _LX) < 1e-9 or abs(y) < 1e-9) else 0
        fy = 1 if abs(y) < 1e-9 else 0
        ops.fix(i + 1, fx, fy)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    dx = _LX / _NX
    for i in top:
        x = nodes[i, 0]
        w = dx if (1e-9 < x < _LX - 1e-9) else dx / 2
        ops.load(int(i + 1), 0.0, -_Q * w)
    et = [e + 1 for e in range(len(elems))]
    tt = [int(i + 1) for i in top]
    args = ["-region", *et, "-Kf", _KF, "-rhoF", 1000.0,
            "-perm", _KBAR, _KBAR, "-poro", _PORO, "-moduli", _E, _NU,
            "-drained", *tt]
    for k, v in ov.items():
        args += ["-" + k] + (list(v) if isinstance(v, (list, tuple)) else [v])
    args += ["-record", rec, 1]
    ops.pattern("LadrunoPorousOverlay", 77, *args)
    ops.system("ProfileSPD")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _run_overlay_column(nsteps, dt, rec, rho=0.0):
    cs = _toy_column()
    _build_overlay_column(rho, cs["nodes"], cs["elems"], cs["top"], rec)
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "overlay column analyze failed"
    times, P, t2c = _read_record(rec)
    # BC discipline: drained region nodes are identically 0 through the run.
    drcols = [t2c[int(i) + 1] for i in cs["top"]]
    assert np.max(np.abs(P[:, drcols])) == 0.0, "drained node p != 0"
    pfree_cols = [t2c[int(i) + 1] for i in cs["ops"]["pfree"]]
    return times, P, t2c, P[:, pfree_cols], cs


# ==========================================================================
# 1. Terzaghi vs analytic + frozen-toy anchor + dt-halving
# ==========================================================================
def _terzaghi(Z, Tv, n=400):
    m = np.arange(n)
    M = (2 * m + 1) * np.pi / 2.0
    return np.sum((2.0 / M) * np.sin(M * Z) * np.exp(-M * M * Tv))


def test_1_terzaghi_vs_analytic_and_toy_anchor():
    cs = _toy_column()
    ops_toy, ff, mat = cs["ops"], cs["ff"], cs["mat"]
    cv = mat["cv"]
    Moed = _E * (1 - _NU) / ((1 + _NU) * (1 - 2 * _NU))
    p0 = _Q * (_KF / _PORO) / ((_KF / _PORO) + Moed)   # undrained excess pore p

    # (a) overlay == frozen-toy fs1_extrap  (the rate form) at Tv=0.1, 50 steps
    t_end = 0.1 * _LY ** 2 / cv
    nsteps, dt = 50, (0.1 * _LY ** 2 / cv) / 50
    _, _, _, Pf, _ = _run_overlay_column(nsteps, dt, _rec("t1_ext.csv"))
    _, P_ext, _, _ = e7.qs_staggered(ops_toy, ff, dt, nsteps, "fs1_extrap",
                                     mat["lfac_classic"], _Q)
    rel_toy = _relL2(Pf, P_ext)
    assert rel_toy <= 1e-5, f"overlay != toy fs1_extrap: relL2 {rel_toy:.2e}"

    # (a) late-window splitting error <= 2x the frozen-toy fs1_extrap|classic
    #     value at matched dt (e7_summary.json e71 compressible/late).
    skip = nsteps // 5
    _, P_mono = e7.qs_mono(ops_toy, ff, dt, nsteps)
    split = _relL2(Pf[skip:], P_mono[skip:])
    toy_split = _relL2(P_ext[skip:], P_mono[skip:])
    assert split <= 2.0 * toy_split, (
        f"late splitting {split:.3e} > 2x toy {toy_split:.3e}")

    # (a) analytic isochrones vs the Terzaghi series (mesh h=H/20 => ~1% floor)
    nst_a = 200
    dt_a = (0.5 * _LY ** 2 / cv) / nst_a
    times, P, t2c, _, _ = _run_overlay_column(nst_a, dt_a, _rec("t1_ana.csv"))
    line, Z = cs["line"], cs["Z"]
    lcols = [t2c[int(i) + 1] for i in line]
    for Tv in (0.1, 0.2, 0.5):
        k = int(round((Tv * _LY ** 2 / cv) / dt_a)) - 1
        TvA = times[k] * cv / _LY ** 2
        pnum = P[k, lcols]
        pana = np.array([p0 * _terzaghi(z, TvA) for z in Z])
        interior = Z > 1e-6
        rel = _relL2(pnum[interior], pana[interior])
        assert rel < 0.02, f"Tv={Tv}: overlay vs Terzaghi relL2 {rel:.3e}"

    # (b) dt-halving: late-window splitting error convergence order >= 0.9
    errsL = []
    te = 0.1 * _LY ** 2 / cv
    for lev in range(4):
        ns = 25 * 2 ** lev
        d = te / ns
        sk = ns // 5
        _, _, _, Pfl, _ = _run_overlay_column(ns, d, _rec(f"t1_c{lev}.csv"))
        _, Pm = e7.qs_mono(ops_toy, ff, d, ns)
        errsL.append(_relL2(Pfl[sk:], Pm[sk:]))
    orders = [np.log2(errsL[i] / errsL[i + 1]) for i in range(len(errsL) - 1)]
    assert min(orders) >= 0.9, f"dt-halving order < 0.9: {orders}"

    print(f"PASS test_1 Terzaghi/anchor: overlay==toy fs1_extrap relL2={rel_toy:.2e}; "
          f"late split {split:.3e} (toy {toy_split:.3e}, <=2x); "
          f"analytic isochrone relL2 ~1%; dt-halving orders "
          f"{[round(o, 2) for o in orders]} (>=0.9)")


# ==========================================================================
# 2. Two-leg vs monolithic LadrunoUP (33017) — mutual dt-convergence
# ==========================================================================
def _run_up_column(nsteps, dt, nodes, elems, top, pfree):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    for i, (x, y) in enumerate(nodes):
        ops.node(i + 1, float(x), float(y))
    ops.nDMaterial("ElasticIsotropic", 1, _E, _NU, 0.0)   # massless: quasi-static
    for e, c in enumerate(elems):
        ops.element("LadrunoUP", e + 1, int(c[0] + 1), int(c[1] + 1),
                    int(c[2] + 1), int(c[3] + 1), 1, "-Kf", _KF, "-poro", _PORO,
                    "-rhoF", 1000.0, "-perm", _KBAR, _KBAR, "-dynSeepage", "off")
    for i, (x, y) in enumerate(nodes):
        fx = 1 if (abs(x) < 1e-9 or abs(x - _LX) < 1e-9 or abs(y) < 1e-9) else 0
        fy = 1 if abs(y) < 1e-9 else 0
        ops.fix(i + 1, fx, fy, 0)
    for i in top:
        ops.fix(int(i + 1), 0, 0, 1)          # drained top p=0
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    dx = _LX / _NX
    for i in top:
        x = nodes[i, 0]
        w = dx if (1e-9 < x < _LX - 1e-9) else dx / 2
        ops.load(int(i + 1), 0.0, -_Q * w, 0.0)
    ops.system("UmfPack")                     # honest-p tangent is unsymmetric
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
    hist = []
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "UP column analyze failed"
        hist.append([ops.nodeDisp(int(i) + 1, 3) for i in pfree])
    return np.array(hist)


def test_2_two_leg_vs_monolithic_up():
    cs = _toy_column()
    cv = cs["mat"]["cv"]
    pfree = cs["ops"]["pfree"]
    te = 0.1 * _LY ** 2 / cv
    diffs = []
    for lev in range(4):
        ns = 25 * 2 ** lev
        d = te / ns
        sk = ns // 5
        _, _, _, Pov, _ = _run_overlay_column(ns, d, _rec(f"t2_{lev}.csv"))
        Pup = _run_up_column(ns, d, cs["nodes"], cs["elems"], cs["top"], pfree)
        diffs.append(_relL2(Pov[sk:], Pup[sk:]))
    orders = [np.log2(diffs[i] / diffs[i + 1]) for i in range(len(diffs) - 1)]
    assert min(orders) >= 1.0, f"overlay-vs-UP mutual order < 1: {orders}"
    assert diffs[-1] < diffs[0], "overlay-vs-UP difference did not shrink"
    print(f"PASS test_2 two-leg vs LadrunoUP: late relL2 "
          f"{[f'{d:.2e}' for d in diffs]}, mutual orders "
          f"{[round(o, 2) for o in orders]} (>=1)")


# ==========================================================================
# 3. Sequencing contract — first committed p ~ undrained (kills p==0 trap)
# ==========================================================================
def test_3_sequencing_first_commit_is_undrained():
    cs = _toy_column()
    cv = cs["mat"]["cv"]
    dt = (0.1 * _LY ** 2 / cv) / 500        # tiny step: near-undrained first commit
    times, P, t2c, Pf, _ = _run_overlay_column(6, dt, _rec("t3.csv"))
    # frozen-toy fs1_extrap first-commit reference (same rate form)
    _, P_ext, _, _ = e7.qs_staggered(cs["ops"], cs["ff"], dt, 6, "fs1_extrap",
                                     cs["mat"]["lfac_classic"], _Q)
    first = Pf[0]
    first_ref = P_ext[0]
    max_first = float(np.max(first) / _Q)
    # (i) hard-fail on the fluid-first p==0 trap
    assert max_first >= 0.05, (
        f"first-commit p collapsed to ~0 (p_max/q={max_first:.3e}) -- the "
        "fluid-first p==0 sequencing trap")
    # (ii) undrained-response magnitude gate
    assert 0.30 <= max_first <= 1.0, (
        f"first-commit p_max/q={max_first:.3f} outside undrained band [0.30,1.0]")
    # (iii) matches the toy rate-form reference within 30%
    rel = _relL2(first, first_ref)
    assert rel < 0.30, f"first-commit vs toy rate-form reference relL2 {rel:.3e}"
    print(f"PASS test_3 sequencing: first-commit p_max/q={max_first:.3f} in "
          f"[0.30,1.0] (>> 0.05), matches toy rate form (relL2 {rel:.2e})")


# ==========================================================================
# 4. Rate-form vs fixed-point contract (the C++ advance is the rate form)
# ==========================================================================
def test_4_rate_form_not_fixed_point():
    cs = _toy_column()
    ops_toy, ff, mat = cs["ops"], cs["ff"], cs["mat"]
    cv = mat["cv"]
    pfree = cs["ops"]["pfree"]
    ns = 50
    d = (0.1 * _LY ** 2 / cv) / ns
    sk = ns // 5
    _, _, _, Pov, _ = _run_overlay_column(ns, d, _rec("t4.csv"))

    # frozen-toy references at the SAME dt
    _, P_ext, _, _ = e7.qs_staggered(ops_toy, ff, d, ns, "fs1_extrap",
                                     mat["lfac_classic"], _Q)   # rate form
    _, P_pln, _, _ = e7.qs_staggered(ops_toy, ff, d, ns, "fs1_plain",
                                     mat["lfac_classic"], _Q)   # fixed-point
    _, P_mono = e7.qs_mono(ops_toy, ff, d, ns)

    # (i) overlay MATCHES the rate form and DIVERGES from the fixed-point form
    to_ext = _relL2(Pov, P_ext)
    to_pln = _relL2(Pov, P_pln)
    assert to_ext <= 1e-5, f"overlay != toy rate form (relL2 {to_ext:.2e})"
    assert to_pln > 0.5, (
        f"overlay too close to the fixed-point form (relL2 {to_pln:.2e}); "
        "the single-pass fixed-POINT L would be O(1) wrong -- rate form NOT wired")

    # (ii) overlay error vs monolithic-UP is in the fs1_extrap CLASS
    Pup = _run_up_column(ns, d, cs["nodes"], cs["elems"], cs["top"], pfree)
    ov_vs_up = _relL2(Pov[sk:], Pup[sk:])
    e_ext = _relL2(P_ext[sk:], P_mono[sk:])      # toy rate-form splitting error
    e_pln = _relL2(P_pln[sk:], P_mono[sk:])      # toy fixed-point O(1) error
    assert ov_vs_up <= 3.0 * e_ext, (
        f"overlay-vs-UP {ov_vs_up:.3e} > 3x toy extrap {e_ext:.3e}")
    assert ov_vs_up <= 0.2 * e_pln, (
        f"overlay-vs-UP {ov_vs_up:.3e} > 0.2x toy fixed-point {e_pln:.3e}")
    print(f"PASS test_4 rate-form contract: overlay==toy_extrap {to_ext:.2e}, "
          f"far from toy_plain {to_pln:.2f}; ov-vs-UP {ov_vs_up:.3e} in "
          f"[<=3x extrap {e_ext:.3e}, <=0.2x fixed-point {e_pln:.3e}]")


# ==========================================================================
# 5. E4-in-OpenSees removal (fluid life-cycle + -onRemoval + subcycle sync)
# ==========================================================================
_RE, _RNU, _RPORO = 2.5e7, 0.3, 0.4
_RN = 10
_RL = 10.0
_RQ = 1.0e4
_RMOED = _RE * (1 - _RNU) / ((1 + _RNU) * (1 - 2 * _RNU))
_RKBAR = 1.0 / _RMOED           # cv = kbar*Moed = 1.0 m^2/s
_RDT = 0.5


def _rntag(i, j):
    return j * (_RN + 1) + i + 1


def _retag(i, j):
    return j * _RN + i + 1


_CRACK = [_retag(i, 4) for i in range(8)]      # row j=4, cols i=0..7 (E4 geometry)
_PROBE = _rntag(2, 2)                            # node (2,2), below the crack


def _build_removal_grid(on_removal, rec, subcycle=1):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j in range(_RN + 1):
        for i in range(_RN + 1):
            ops.node(_rntag(i, j), i * _RL / _RN, j * _RL / _RN)
    ops.nDMaterial("ElasticIsotropic", 1, _RE, _RNU, 0.0)
    for j in range(_RN):
        for i in range(_RN):
            ops.element("quad", _retag(i, j), _rntag(i, j), _rntag(i + 1, j),
                        _rntag(i + 1, j + 1), _rntag(i, j + 1), 1.0, "PlaneStrain", 1)
    for j in range(_RN + 1):
        for i in range(_RN + 1):
            ops.fix(_rntag(i, j), 1 if (i == 0 or i == _RN or j == 0) else 0,
                    1 if j == 0 else 0)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    dx = _RL / _RN
    for i in range(_RN + 1):
        w = dx if (0 < i < _RN) else dx / 2
        ops.load(_rntag(i, _RN), 0.0, -_RQ * w)
    et = [_retag(i, j) for j in range(_RN) for i in range(_RN)]
    tt = [_rntag(i, _RN) for i in range(_RN + 1)]
    args = ["-region", *et, "-Kf", _KF, "-rhoF", 1000.0, "-perm", _RKBAR, _RKBAR,
            "-poro", _RPORO, "-moduli", _RE, _RNU, "-drained", *tt]
    if on_removal == "keep":
        args += ["-onRemoval", "keep"]
    else:
        args += ["-onRemoval", "drain", float(on_removal)]
    if subcycle > 1:
        args += ["-subcycle", subcycle]
    args += ["-record", rec, 1]
    ops.pattern("LadrunoPorousOverlay", 77, *args)
    ops.system("ProfileSPD")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")


def _march_removal(on_removal, rec, nst, rm_step, subcycle=1):
    _build_removal_grid(on_removal, rec, subcycle=subcycle)
    for s in range(nst):
        if s == rm_step:
            for e in _CRACK:
                ops.remove("element", e)
        assert ops.analyze(1, _RDT) == 0, f"removal march failed at step {s}"
    return _read_record(rec)


def test_5_removal_fluid_lifecycle():
    nst, rm_step = 280, 8            # remove early (Tv~0.04), drainage-dominated
    rm_time = rm_step * _RDT

    def probe_curve(t, P, t2c):
        return P[:, t2c[_PROBE]] / _RQ

    def t_to_thresh(t, pc, thr=0.1):
        post = np.nonzero((t > rm_time + 2 * _RDT) & (pc <= thr))[0]
        return t[post[0]] if len(post) else np.inf

    curves = {}
    for pol, name in (("keep", "keep"), (10, "d10"), (100, "d100")):
        t, P, t2c = _march_removal(pol, _rec(f"t5_{name}.csv"), nst, rm_step)
        assert np.all(np.isfinite(P)), f"{name}: non-finite p (removal blew up)"
        # BC discipline: drained top nodes identically 0
        drcols = [t2c[_rntag(i, _RN)] for i in range(_RN + 1)]
        assert np.max(np.abs(P[:, drcols])) == 0.0, f"{name}: drained p != 0"
        curves[name] = (t, probe_curve(t, P, t2c), t2c)

    tk, pk, _ = curves["keep"]
    td, pd, _ = curves["d100"]
    _, p10, _ = curves["d10"]

    # (a1) KEEP: fluid persists -> p drains THROUGH the crack, no trapped plateau
    p_before = pk[rm_step - 1]
    assert p_before > 0.3, f"probe not loaded before removal (p={p_before:.3f}q)"
    assert pk[-1] < 0.05, f"keep: probe p trapped at {pk[-1]:.3f}q (should drain)"
    t10_keep = t_to_thresh(tk, pk)
    assert np.isfinite(t10_keep), "keep: probe never drained below 0.1q"

    # (a2) drain kFactor 100 strictly FASTER than keep (late-time metric)
    t10_d100 = t_to_thresh(td, pd)
    assert t10_d100 < t10_keep - _RDT, (
        f"drain100 not faster: t(<=0.1q) drain100={t10_d100:.1f} keep={t10_keep:.1f}")
    # monotone at a fixed late time (after the redistribution transient settles)
    k60 = np.searchsorted(tk, 60.0) - 1
    p_keep60, p_d10_60, p_d100_60 = pk[k60], p10[k60], pd[k60]
    assert p_d100_60 < p_keep60, (
        f"drain100 p@t60 {p_d100_60:.3f} not < keep {p_keep60:.3f}")
    assert p_d10_60 <= p_keep60 + 1e-6, "drain10 p@t60 exceeds keep"

    # (b) fixed-window jump (E7.5 metric) finite + march stable
    i_before = np.searchsorted(tk, rm_time - 1e-9) - 1
    i_after = min(np.searchsorted(tk, rm_time + 20 * _RDT - 1e-9), len(tk) - 1)
    jwin = (pk[i_after] - pk[i_before])
    assert np.isfinite(jwin) and abs(jwin) < 2.0, f"window jump {jwin:.3f}q not bounded"

    # (c) removal forces a fluid sync under -subcycle N>1
    rm2 = 17                                 # NOT a multiple of the 5-cadence
    t_sc, _, _ = _march_removal("keep", _rec("t5_sc.csv"), 60, rm2, subcycle=5)
    t_ref, _, _ = _march_removal("keep", _rec("t5_scref.csv"), 60, 10 ** 9, subcycle=5)
    forced_time = (rm2 + 1) * _RDT            # removal detected at this commit's time
    cadence = 5 * _RDT
    off_cadence = abs((forced_time / cadence) - round(forced_time / cadence)) > 1e-6
    assert off_cadence, "test bug: pick a removal step off the subcycle cadence"
    has_forced = np.any(np.abs(t_sc - forced_time) < 1e-6)
    ref_has = np.any(np.abs(t_ref - forced_time) < 1e-6)
    assert has_forced and not ref_has, (
        f"removal did not force an off-cadence fluid sync at t={forced_time} "
        f"(sync rows {list(np.round(t_sc, 1))[:14]})")

    print(f"PASS test_5 removal: keep drains through crack (p {p_before:.2f}q->"
          f"{pk[-1]:.3f}q); drain100 faster (t<=0.1q {t10_d100:.1f}<{t10_keep:.1f}); "
          f"p@t60 keep {p_keep60:.3f}>d10 {p_d10_60:.3f}>=d100 {p_d100_60:.3f}; "
          f"window jump {jwin:.3f}q finite; subcycle sync forced at t={forced_time}")


# ==========================================================================
# 6. -fsL manual floor loud warning
# ==========================================================================
def test_6_fsL_floor_loud_warning():
    logf = _rec("opserr_fsL.log")
    if os.path.exists(logf):
        os.remove(logf)
    ops.logFile(logf)                          # capture opserr (C++ stream)
    try:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        for j in range(3):
            for i in range(3):
                ops.node(3 * j + i + 1, float(i), float(j))
        ops.nDMaterial("ElasticIsotropic", 1, _RE, _RNU, 0.0)
        conn = [(1, 2, 5, 4), (2, 3, 6, 5), (4, 5, 8, 7), (5, 6, 9, 8)]
        for e, c in enumerate(conn, 1):
            ops.element("quad", e, *c, 1.0, "PlaneStrain", 1)
        for n in (1, 2, 3):
            ops.fix(n, 1, 1)
        ops.timeSeries("Constant", 1)
        ops.pattern("Plain", 1, 1)
        ops.load(7, 0.0, -1.0e3)
        ops.load(9, 0.0, -1.0e3)
        # -fsL 0.1 drives L below the oedometric floor -> clamp + loud warning
        ops.pattern("LadrunoPorousOverlay", 88, "-region", 1, 2, 3, 4,
                    "-Kf", _KF, "-rhoF", 1000.0, "-perm", _RKBAR, _RKBAR,
                    "-poro", _RPORO, "-moduli", _RE, _RNU, "-fsL", 0.1,
                    "-drained", 7, 8, 9)
        ops.system("ProfileSPD")
        ops.numberer("Plain")
        ops.constraints("Transformation")
        ops.algorithm("Linear")
        ops.integrator("Newmark", 0.6, 0.3025)
        ops.analysis("Transient")
        rc = ops.analyze(3, 0.5)              # run still COMPLETES
    finally:
        ops.logFile("nul")                     # detach (best-effort)
    assert rc == 0, f"run did not complete with -fsL floor clamp (rc={rc})"
    assert os.path.exists(logf), "opserr logfile not produced"
    txt = open(logf, encoding="utf-8", errors="replace").read()
    warned = any(("floor" in ln.lower() and "warning" in ln.lower())
                 for ln in txt.splitlines())
    assert warned, f"no loud floor WARNING captured; log:\n{txt[:500]}"
    print("PASS test_6 -fsL floor: manual scale 0.1 clamped to oedometric floor "
          "with a loud WARNING; run completed (rc=0)")


# ==========================================================================
# standalone runner (python -S tests/test_ladruno_overlay_physics.py)
# ==========================================================================
_ALL = [
    ("test_1_terzaghi_vs_analytic_and_toy_anchor", test_1_terzaghi_vs_analytic_and_toy_anchor),
    ("test_2_two_leg_vs_monolithic_up", test_2_two_leg_vs_monolithic_up),
    ("test_3_sequencing_first_commit_is_undrained", test_3_sequencing_first_commit_is_undrained),
    ("test_4_rate_form_not_fixed_point", test_4_rate_form_not_fixed_point),
    ("test_5_removal_fluid_lifecycle", test_5_removal_fluid_lifecycle),
    ("test_6_fsL_floor_loud_warning", test_6_fsL_floor_loud_warning),
]


def main():
    if not _BUILT:
        print(f"SKIP: worktree engine not built at {_DIST}")
        return 0
    fails = []
    for name, fn in _ALL:
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            import traceback
            print(f"FAIL {name}: {exc}")
            traceback.print_exc()
            fails.append(name)
    print("=" * 74)
    if fails:
        print(f"FAILED {len(fails)}/{len(_ALL)}: {fails}")
        return 1
    print(f"ALL {len(_ALL)} PHYSICS GATES PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
