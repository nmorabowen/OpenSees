"""LadrunoPorousOverlay (PATTERN 33022) — ADR-73 P3b EXPLICIT-FLUID battery (WP3b.C).

P3b lane (plan "P3b PINNED at 3b.A" §3b.1-§3b.5; E7.4 frozen constants):
  * `-fluidUpdate explicit` -> FU_EXPLICIT: lumped-S* forward p-step at each
    sync commit, NO CG/factorization in the march
    (pTrial = pCommitted + (-dt*(H p0) - Qt du + dt f_seep)/sLump on FREE rows).
  * `-fsL` is INERT under FU_EXPLICIT (no L anywhere; a non-default -fsL gets a
    one-time notice, not a fatal). `-staticMode hold|steady` stay legal (STEADY
    keeps its CG solve at commits — not a march cost). `-stab` is march-inert
    (rowsum(S*) = rowsum(S_phys), H-tilde annihilates constants).
  * LadrunoStaggeredAnalyze REFUSES FU_EXPLICIT overlays (loud fatal, empty
    driven set — LadrunoStaggeredDriver.cpp:118 re-pinned here).
  * criticalTimeStep() under an FU_EXPLICIT overlay min-folds the explicit-fluid
    diffusion bound dt_diff = minEdge^2/(2*ndm*maxCv) into damped+undamped dt
    and the CDL/ExplicitBathe report appends an explicit-fluid diffusion-limit
    line with the slack factor.
  * SMS (BOTH builders) sizes against the UNDRAINED per-element pencil via the
    elementCriticalDt Kadd seam; the blanket "SMS + porous overlay is
    UNSUPPORTED" warning is RETIRED, replaced by the one-time INFO line
    "SMS sizing priced the UNDRAINED pencil for N overlay-owned elements
    (ADR-73 P3b)".

FROZEN MEASURED CONSTANTS (plan frozen-constants block, E7.4):
  boundary = 1.32x the DISCRETE undrained pencil on both benchmark soils
  (gate band = the e72/e74 +-25 percent spread convention -> [0.99, 1.65]);
  diffusion CFL slack 7.2e3x at realistic kbar (gate floor 1e2x — the POINT is
  orders of slack, direction pinned); material factor ~1.85x conservative =
  documentation-only; L not needed (the reference rule is inert).

THE GATE MAP (plan §3b.5; letters match the pins):
  (a) diffusion-slack sweep, realistic kbar (kh 1e-7..1e-3 m/s class): the
      advisory's dt_diff matches the h^2/(2*ndm*c_v) hand formula (c_v =
      kbar_max * M_oed — the overlay's own maxCv_ reduction,
      LadrunoPorousOverlay.cpp assembleStorageAndL), slack >= 1e2x at every
      realistic point; PLUS the min-fold honesty probe: at an ABSURD kbar the
      returned advisory == the hand dt_diff (diffusion governs, fold honest);
      PLUS the report line: the CDL "-cfl" report prints a diffusion line whose
      dt_diff matches the hand value.
  (b) coupled CFL envelope vs the 1.32x pin: fixed-horizon (4500-step) Dt
      bisection on the e74-matched DRAINING column; boundary/undrained-pencil
      in [0.99, 1.65]; 0.8x pencil marches a 3000-step horizon stable.
  (c) implicit-lane equivalence: FU_EXPLICIT vs FU_IMPLICIT + `-fsL zero` at a
      matched Dt sequence (both O(Dt) splits — mutual-convergence): inter-lane
      diff shrinks with observed order >= 0.9 (ZS84 column, P3 methodology).
  (d) SMS composability, BOTH builders (THE gate — the inverted P3 quirk): the
      pre-fix pathology (unscaled march at the drained-priced-certified
      dtTarget) provably blows up; with the fix the SMS-certified dtTarget
      march is stable over a 4000-step horizon, the report prices the
      undrained pencil (scaled-count direction vs a drained/no-overlay
      control), the old UNSUPPORTED warning is absent and the exact INFO line
      present. Lumped CentralDifferenceSMS + CentralDifferenceSMSConsistent.
  (e) MP — TRANSPOSED at 3b.A: NO test. The overlay is per-process serial v1;
      the gate became (i) a partition-locality code audit (panel duty, not
      run under MPI) and (ii) the halo DESIGN documented in the ADR §8
      amendment. Recorded here for the gate map only — do not overclaim.
  (f) recorder CSV-twin under FU_EXPLICIT: `recorder ladruno -overlay` HDF5
      DATA == `-record` CSV exact (<= 1e-12 rel), ID == header == REGION_NODES,
      drained topology rows intact, drained columns identically 0.
  (g) zero-drainage secular pumping MEASURED under FU_EXPLICIT: the P3 quirks
      scenario (kbar ~ 1e-11, undamped) at 0.4x/0.5x/1.0x pencil,
      steps-to-blowup (or full-horizon survival) printed in a loud RECORD
      block; asserts only what is honest (coherence/monotonicity of the
      measured horizons); the numbers feed the quirks-row amendment.
  Smokes: driver refusal fatal; -fsL inert under explicit (notice + marches
      bit-identical across default/classic/oedometric, 1e-12 band); SM_STEADY
      gravity-init then transient march runs; DB save/restore round-trip
      mid-march resumes bit-consistently (frozen-solid FileDatastore child).

CONTRACT AMBIGUITIES RESOLVED HERE (flag for MAIN):
  * "-fsL classic" IS the parser default — an explicitly-passed classic may
    not count as "non-default", so the one-time notice is asserted on the
    unambiguously non-default `-fsL oedometric`; the classic leg is asserted
    march-identical only.
  * The diffusion report-line wording is 3b.B's; gate (a) parses any line
    containing "diffusion" (case-insensitive) and gates on the NUMBER.
  * SM_STEADY is a whole-overlay mode (onDomainCommit has no static-vs-
    transient detection — the flag governs every commit), so the smoke runs a
    static stage then a transient march on the SAME steady overlay.
  * Gate (d)'s stable march uses a DRAINING kbar (1e-8) + small alphaM
    Rayleigh: the zero-drainage secular pumping is a lane-inherited pathology
    owned by gate (g), not an SMS-sizing defect (quirks row: any drainage or
    damping absorbs it).
  * C++ stderr/stdout observation is done via fresh child subprocesses (the
    house pattern in this family's batteries), not capfd — works under the
    standalone runner too. (Memory pin "capfd not capsys" honored by not
    using capsys anywhere.)

RUNNER (worktree dist/bin opensees.pyd is CPython 3.12; ~5-10 min):
  C:\\Users\\nmora\\AppData\\Local\\Python\\pythoncore-3.12-64\\python.exe -S \\
      tests/test_ladruno_overlay_explicit_fluid.py      # standalone, exit 1 on fail
  py -3.12 -m pytest tests/test_ladruno_overlay_explicit_fluid.py -x -q  # pytest
`-S` skips site (evicts any boot-.pth preloaded stale pyd); the bootstrap re-adds
site-packages for numpy/h5py and asserts THIS worktree's engine won.
"""
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

# --------------------------------------------------------------------------
# bootstrap (copied VERBATIM from tests/test_ladruno_overlay_explicit.py —
# the pinned runner discipline; h5py added for gate (f), LEDGER_quirks
# precedent: disable libhdf5 file locking BEFORE importing h5py)
# --------------------------------------------------------------------------
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

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
    for _m in ("opensees", "openseespy", "openseespy.opensees"):
        sys.modules.pop(_m, None)
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

import opensees as ops  # noqa: E402

assert os.path.normcase(os.path.dirname(ops.__file__)) == os.path.normcase(_DIST), (
    f"wrong opensees.pyd imported: {ops.__file__} (want {_DIST})"
)

try:
    import pytest  # noqa: E402
    _HAVE_PYTEST = True
except Exception:  # pragma: no cover
    _HAVE_PYTEST = False

if _HAVE_PYTEST:
    pytestmark = [pytest.mark.zone_b, pytest.mark.t1]

# h5py is needed ONLY by gate (f); missing h5py skips that gate, not the module.
_HAVE_H5PY = True
_H5_ERR = ""
try:
    import h5py  # noqa: E402
except Exception as _exc:  # pragma: no cover
    _HAVE_H5PY = False
    _H5_ERR = repr(_exc)

_TMP = tempfile.mkdtemp(prefix="ov_expfl_")


def _rec(name):
    return os.path.join(_TMP, name)


class _Skip(Exception):
    """Standalone-runner analog of pytest.skip (a gate that cannot run, not a fail)."""


def _skip(msg):
    if _HAVE_PYTEST:
        pytest.skip(msg)
    raise _Skip(msg)


# ==========================================================================
# frozen configurations (copied from the P3 battery — e72/ZS84 pins)
# ==========================================================================
# e72/e74-matched column: E=1e7 nu=0.25 Kf=2.2e9 poro=0.4, kh=1e-7 gw=1e4 ->
# kbar=1e-11, mixture rho=2.0; 1x8 unit-cell mesh (h = 1). The E7.4 explicit-
# fluid boundary was measured on THIS soil class: 1.32x the discrete undrained
# pencil, diffusion slack 7.2e3x.
_E72 = dict(nely=8, E=1.0e7, nu=0.25, rho=2.0, rhoF=1000.0,
            Kf=2.2e9, poro=0.4, kbar=1.0e-7 / 1.0e4, q=1.0e5)

# ZS84 Example-1 column (ADR §7.1 B2 pins): 30 m column, q = 1 ramped by
# sin(pi*t/2tp) tp = 0.1 s, E = 3e4, nu = 0.2, rho = 2 (MIXTURE), rhoF = 1,
# Kf = 1e5, n = 0.3, kbar = 1e-2/(rhoF*9.81); top drained, base fixed,
# laterally confined. h = 1 m (nely = 30). DRAINING — the gate-(c) rig.
_Z = dict(L=30.0, nely=30, E=3.0e4, nu=0.2, rho=2.0, rhoF=1.0,
          Kf=1.0e5, poro=0.3, kbar=1.0e-2 / (1.0 * 9.81), q=1.0, tp=0.1)


def _moed_cpp(E, nu):
    """M_oed EXACTLY as the overlay computes it (moduliOf,
    LadrunoPorousOverlay.cpp:75): Kdr + 4G/3 — same floating-point path as the
    C++ so the gate-(a) fold identity can be tight."""
    Kdr = E / (3.0 * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    return Kdr + 4.0 * G / 3.0


def _dt_diff_hand(kbar, E, nu, h=1.0, ndm=2):
    """The FROZEN 3b.3 advisory formula: dt_diff = minEdge^2/(2*ndm*maxCv),
    c_v = kbar_max * M_oed (the overlay's assembleStorageAndL reduction)."""
    return h * h / (2.0 * ndm * (kbar * _moed_cpp(E, nu)))


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


# node tag for column grid position (i in {0,1}, j in 0..nely), 1-based
def _nt(i, j):
    return 2 * j + i + 1


# ==========================================================================
# parent-side builders (copied from the P3 battery, -fluidUpdate added)
# ==========================================================================
def _build_e72(fluid="explicit", kbar=None, fsl=None, record=None,
               overlay=True, stab=("auto", 0.25), static_steady=False,
               series="Constant"):
    """e72 column, plain ndf=2 quads (MIXTURE rho via the material — no body
    forces, the staggered-twin trap is NOT armed) + overlay + CDL-ready mesh.
    fluid: 'explicit' -> -fluidUpdate explicit (default -fsL, classic);
           'implicit_fsl0' -> the P3 lane (-fsL zero, default FU_IMPLICIT)."""
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
    ops.timeSeries(series, 1)
    ops.pattern("Plain", 1, 1)
    for i in range(2):
        ops.load(_nt(i, c["nely"]), 0.0, -c["q"] / 2.0)
    if overlay:
        top = [_nt(0, c["nely"]), _nt(1, c["nely"])]
        args = ["-region", *eles, "-Kf", c["Kf"], "-rhoF", c["rhoF"],
                "-perm", kbar, kbar, "-poro", c["poro"],
                "-moduli", c["E"], c["nu"], "-drained", *top,
                "-stab", *stab]
        if fluid == "explicit":
            args += ["-fluidUpdate", "explicit"]
        elif fluid == "implicit_fsl0":
            args += ["-fsL", "zero"]
        else:
            raise AssertionError(f"test bug: unknown fluid lane {fluid!r}")
        if fsl is not None:
            args += ["-fsL", fsl]
        if static_steady:
            args += ["-staticMode", "steady"]
        if record is not None:
            args += ["-record", record, 1]
        ops.pattern("LadrunoPorousOverlay", 77, *args)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")
    return eles


def _zs84_sin_ramp():
    """The B2 sin ramp: 21-point Path over [0, tp], held at 1.0."""
    z = _Z
    ts = [z["tp"] * k / 20.0 for k in range(21)] + [1.0e12]
    vs = [math.sin(math.pi * t / (2.0 * z["tp"])) for t in ts[:21]] + [1.0]
    ops.timeSeries("Path", 1, "-time", *ts, "-values", *vs)


def _build_zs84_lane(lane, record=None):
    """ZS84 column, plain ndf=2 quads + overlay(-stab auto 0.25) + CDL.
    lane: 'explicit' (FU_EXPLICIT, default -fsL) or 'implicit_fsl0'
    (FU_IMPLICIT + -fsL zero — the P3 explicit-lane twin, gate (c))."""
    z = _Z
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
    top = [_nt(0, z["nely"]), _nt(1, z["nely"])]
    args = ["-region", *eles, "-Kf", z["Kf"], "-rhoF", z["rhoF"],
            "-perm", z["kbar"], z["kbar"], "-poro", z["poro"],
            "-moduli", z["E"], z["nu"], "-drained", *top,
            "-stab", "auto", 0.25]
    if lane == "explicit":
        args += ["-fluidUpdate", "explicit"]
    elif lane == "implicit_fsl0":
        args += ["-fsL", "zero"]
    else:
        raise AssertionError(f"test bug: unknown lane {lane!r}")
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
# measured pencils (memoized) — the CDL-7 idiom: one tiny priming step makes
# ops.criticalTimeStep() valid. Under P3b the returned dt is
# min(undrained pencil, dt_diff); at these kbar the diffusion slack is orders,
# so the return IS the discrete undrained pencil (gate (a) pins the fold).
# ==========================================================================
_PEN = {}


def _pencil(key):
    if key in _PEN:
        return _PEN[key]
    builders = {
        "e72_exp":    lambda: _build_e72(fluid="explicit"),
        "e72_exp_k8": lambda: _build_e72(fluid="explicit", kbar=1.0e-8),
        "e72_imp0":   lambda: _build_e72(fluid="implicit_fsl0"),
        "zs84_exp":   lambda: _build_zs84_lane("explicit"),
    }
    builders[key]()
    assert ops.analyze(1, 1.0e-9) == 0, f"pencil priming step failed ({key})"
    dtcr = ops.criticalTimeStep()
    ops.wipe()
    assert dtcr is not None and math.isfinite(dtcr) and dtcr > 0.0, (
        f"criticalTimeStep() returned {dtcr!r} for {key}")
    _PEN[key] = float(dtcr)
    return _PEN[key]


def _march_guard(nsteps, dt, probe_node, probe_dof=2, chunks=8):
    """March the CURRENT explicit model in chunks with finite/blow-up guards
    (cap 1e3 — grotesque vs the mm..cm response scales here)."""
    if nsteps % chunks:
        chunks = 1
    step = nsteps // chunks
    for c in range(chunks):
        assert ops.analyze(step, dt) == 0, (
            f"march failed in chunk {c} (steps {c * step}..{(c + 1) * step})")
        u = ops.nodeDisp(probe_node, probe_dof)
        assert math.isfinite(u) and abs(u) < 1.0e3, (
            f"march diverged by chunk {c}: u={u}")


# ==========================================================================
# fresh-subprocess machinery (numpy-free children; copied from the P3 battery)
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


def _run_child(body, timeout=300):
    fd, path = tempfile.mkstemp(suffix=".py", dir=_TMP)
    os.close(fd)
    with open(path, "w") as f:
        f.write(body)
    proc = subprocess.run([sys.executable, "-S", "-u", path],
                          capture_output=True, text=True, timeout=timeout)
    return proc.stdout + proc.stderr, proc.returncode


# ---- e72 column child builder (parametric kbar + overlay tail args) --------
_E72_CH = '''
KBAR = %(kbar)r
OVARGS = %(ovargs)r
def build_e72(overlay=True):
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
    if overlay:
        args = ["-region"] + eles + ["-Kf", 2.2e9, "-rhoF", 1000.0,
                "-perm", KBAR, KBAR, "-poro", 0.4, "-moduli", 1.0e7, 0.25,
                "-drained", 17, 18, "-stab", "auto", 0.25] + list(OVARGS)
        ops.pattern("LadrunoPorousOverlay", 77, *args)
    return eles
'''

# ---- generic march child (gates b, d, g): RESULT DIVERGED step N / BOUNDED --
_MARCH_CH = _BOOT_CH + _E72_CH + '''
DT = %(dt)r
NST = %(nst)d
INTG = %(intg)r
ALPHAM = %(alpham)r
OV = %(ov)r
build_e72(overlay=OV)
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("ProfileSPD")
ops.algorithm("Linear")
if ALPHAM:
    ops.rayleigh(ALPHAM, 0.0, 0.0, 0.0)
ops.integrator(*INTG)
ops.analysis("Transient")
mx = 0.0
for s in range(NST):
    try:
        rc = ops.analyze(1, DT)
    except BaseException as ex:
        print("RESULT DIVERGED step", s, "exc", type(ex).__name__); sys.exit(0)
    u = ops.nodeDisp(17, 2)
    if rc != 0 or not math.isfinite(u) or abs(u) > 1.0e3:
        print("RESULT DIVERGED step", s, "rc", rc, "u", u); sys.exit(0)
    if abs(u) > mx: mx = abs(u)
print("RESULT BOUNDED maxu", mx)
sys.exit(0)
'''


def _march_child_src(dt, nst, ovargs, intg=("CentralDifferenceLadruno",),
                     kbar=_E72["kbar"], alpham=0.0, ov=True):
    return _MARCH_CH % dict(dist=_DIST, kbar=float(kbar),
                            ovargs=list(ovargs), dt=float(dt), nst=int(nst),
                            intg=tuple(intg), alpham=float(alpham),
                            ov=bool(ov))


def _div_step(out):
    for ln in out.splitlines():
        if ln.startswith("RESULT DIVERGED"):
            try:
                return int(ln.split()[3])
            except (IndexError, ValueError):
                return -1
    return -1


def _bounded(out, rc):
    return (rc == 0) and ("RESULT BOUNDED" in out)


# ---- gate (a) report child: "-cfl" makes CDL print the CTS report ----------
_REPORT_CH = _BOOT_CH + _E72_CH + '''
build_e72()
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("ProfileSPD")
ops.algorithm("Linear"); ops.integrator("CentralDifferenceLadruno", "-cfl")
ops.analysis("Transient")
assert ops.analyze(1, 1.0e-9) == 0
print("DTCR", repr(ops.criticalTimeStep()))
sys.exit(0)
'''

# ---- smoke: driver refusal (FU_EXPLICIT overlay, P2 e5 dual handling) ------
_REFUSAL_CH = _BOOT_CH + _E72_CH + '''
build_e72()
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

# ---- smoke: non-default -fsL notice under FU_EXPLICIT (build + short march) -
_NOTICE_CH = _BOOT_CH + _E72_CH + '''
build_e72()
ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("ProfileSPD")
ops.algorithm("Linear"); ops.integrator("CentralDifferenceLadruno")
ops.analysis("Transient")
for s in range(3):
    assert ops.analyze(1, %(dt)r) == 0
print("BUILD_OK")
sys.exit(0)
'''

# ---- smoke: FileDatastore DB round-trip (framework battery's frozen-solid ---
#      pattern, -fluidUpdate explicit added; continuous 6 vs restart 3+3) ----
_DB_CH = _BOOT_CH + '''
import tempfile, os as _os
FR = tempfile.mkdtemp(prefix="ovdbx_")
NELY = 6; L = 6.0; GW = 9810.0
def col2():
    ops.wipe(); ops.model("basic", "-ndm", 2, "-ndf", 2)
    nid = {}; k = 1
    for j in range(NELY + 1):
        for i in range(2):
            ops.node(k, float(i), L * j / NELY); nid[(i, j)] = k; k += 1
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    eles = []
    for j in range(NELY):
        et = j + 1
        ops.element("quad", et, nid[(0, j)], nid[(1, j)], nid[(1, j + 1)],
                    nid[(0, j + 1)], 1.0, "PlaneStrain", 1)
        eles.append(et)
    for j in range(NELY + 1):
        for i in range(2):
            ops.fix(nid[(i, j)], 1, 1)          # FROZEN solid: u == 0 exactly
    return eles, [nid[(0, NELY)], nid[(1, NELY)]]
def run():
    ops.constraints("Penalty", 1e14, 1e14); ops.numberer("RCM")
    ops.system("ProfileSPD"); ops.test("NormDispIncr", 1e-10, 25)
    ops.algorithm("Newton"); ops.integrator("Newmark", 0.6, 0.3025)
    ops.analysis("Transient")
def mk(rec):
    eles, top = col2()
    ops.pattern("LadrunoPorousOverlay", 77, "-region", *eles,
                "-Kf", 2.2e7, "-rhoF", 1000.0, "-perm", 1e-6, 1e-6,
                "-poro", 0.4, "-moduli", 1.0e4, 0.2, "-drained", *top,
                "-fluidUpdate", "explicit",
                "-pInit", "hydrostatic", GW, L, "-record", rec, 1)
def rows(p):
    with open(p) as f: LNS = f.read().strip().splitlines()
    return [[float(x) for x in r.split(",")] for r in LNS[1:]]
recC = _os.path.join(FR, "c.csv"); mk(recC); run()
for _ in range(6): assert ops.analyze(1, 0.1) == 0
recR = _os.path.join(FR, "r.csv"); mk(recR); run()
for _ in range(3): assert ops.analyze(1, 0.1) == 0
dbp = _os.path.join(FR, "db")
ops.database("File", dbp); ops.save(1); ops.wipe()
ops.database("File", dbp); ops.restore(1); run()
for _ in range(3): assert ops.analyze(1, 0.1) == 0
rC = rows(recC); rR = rows(recR); tC = {round(r[0], 6): r for r in rC}
dmax = max(max(abs(a - b) for a, b in zip(r, tC[round(r[0], 6)])) for r in rR)
scale = max(abs(x) for r in rC for x in r[1:]) or 1.0
print("DMAX", dmax, "REL", dmax / scale)
sys.exit(0)
'''

_FLOAT_RE = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?")


def _floats_in(text):
    out = []
    for m in _FLOAT_RE.finditer(text):
        try:
            out.append(float(m.group(0)))
        except ValueError:
            pass
    return out


# ==========================================================================
# (a) diffusion-slack sweep + min-fold honesty + report line
# ==========================================================================
def test_a_diffusion_slack():
    """PREDICTION (before assert): at realistic kbar the FU_EXPLICIT advisory
    is untouched by the diffusion fold (slack >= 1e2x, E7.4 measured 7.2e3x at
    kh=1e-7-class on the assembled toy operator; the per-element C++ pencil is
    ~1.7x smaller so the slack here is LARGER), the hand formula
    dt_diff = h^2/(2*ndm*c_v) with c_v = kbar*M_oed and h = minEdge = 1 (unit
    cells, elemSizeH = largest edge) matches the report's printed dt_diff, and
    at an ABSURD kbar (1e-1) the returned advisory == the hand dt_diff exactly
    (the min-fold is honest when diffusion governs). The FU_EXPLICIT and
    FU_IMPLICIT(-fsL zero) advisories agree at realistic kbar (the augmentation
    is lane-independent; the fold does not bind)."""
    c = _E72
    und = _pencil("e72_exp")
    und_imp = _pencil("e72_imp0")
    assert abs(und - und_imp) <= 1e-9 * und, (
        f"FU_EXPLICIT advisory {und:.6e} != FU_IMPLICIT(-fsL zero) advisory "
        f"{und_imp:.6e} at realistic kbar — the diffusion fold bound at a "
        f"point where its slack is orders, or the augmentation drifted by lane")

    # realistic sweep: kh in {1e-7, 1e-5, 1e-3} m/s (gw = 1e4 -> kbar 1e-11..1e-7)
    lines = []
    for kh in (1.0e-7, 1.0e-5, 1.0e-3):
        kbar = kh / 1.0e4
        _build_e72(fluid="explicit", kbar=kbar)
        assert ops.analyze(1, 1.0e-9) == 0
        dtcr = float(ops.criticalTimeStep())
        ops.wipe()
        dtd = _dt_diff_hand(kbar, c["E"], c["nu"])
        slack = dtd / dtcr
        assert slack >= 1.0e2, (
            f"kbar={kbar:.1e}: diffusion slack {slack:.2e} < 1e2x (dt_diff "
            f"hand {dtd:.3e} vs advisory {dtcr:.3e}) — realistic-kbar slack "
            f"must be ORDERS (E7.4 class 7.2e3x)")
        # pencil must be kbar-independent (storage/stiffness-based) and the
        # fold must not bind here
        assert abs(dtcr - und) <= 1e-6 * und, (
            f"kbar={kbar:.1e}: advisory {dtcr:.6e} moved off the undrained "
            f"pencil {und:.6e} — the diffusion fold bound at a realistic kbar")
        lines.append(f"kbar={kbar:.0e}: slack={slack:.2e}")

    # min-fold honesty: ABSURD kbar — diffusion governs, return == hand dt_diff
    kbar_abs = 1.0e-1
    dtd_abs = _dt_diff_hand(kbar_abs, c["E"], c["nu"])
    assert dtd_abs < 0.5 * und, "test bug: absurd kbar did not push dt_diff below the pencil"
    _build_e72(fluid="explicit", kbar=kbar_abs)
    assert ops.analyze(1, 1.0e-10) == 0
    dtcr_abs = float(ops.criticalTimeStep())
    ops.wipe()
    assert abs(dtcr_abs - dtd_abs) <= 1e-9 * dtd_abs, (
        f"min-fold DISHONEST at absurd kbar: returned {dtcr_abs:.9e} != hand "
        f"dt_diff {dtd_abs:.9e} (= minEdge^2/(2*ndm*c_v), the FROZEN 3b.3 "
        f"formula) while diffusion governs")

    # report line: the CDL "-cfl" report must print a diffusion line whose
    # dt_diff matches the hand value (wording is 3b.B's; the NUMBER is gated)
    out, rc = _run_child(_REPORT_CH % dict(
        dist=_DIST, kbar=float(c["kbar"]),
        ovargs=["-fluidUpdate", "explicit"]))
    assert rc == 0 and "DTCR" in out, f"report child failed:\n{out[-1500:]}"
    dtd0 = _dt_diff_hand(c["kbar"], c["E"], c["nu"])
    diff_lines = [ln for ln in out.splitlines() if "diffusion" in ln.lower()]
    assert diff_lines, (
        "CDL '-cfl' report printed NO line containing 'diffusion' under an "
        "FU_EXPLICIT overlay (3b.3: the report appends an explicit-fluid "
        f"diffusion-limit line):\n{out[-1500:]}")
    cand = [v for ln in diff_lines for v in _floats_in(ln)]
    hit = [v for v in cand if v > 0 and abs(v - dtd0) <= 5e-3 * dtd0]
    assert hit, (
        f"no float on the report's diffusion line(s) matches the hand dt_diff "
        f"{dtd0:.4e} within 0.5% (floats seen: {cand[:12]}); lines:\n"
        + "\n".join(diff_lines))

    info = ("; ".join(lines)
            + f"; fold@kbar=1e-1: return {dtcr_abs:.4e} == hand {dtd_abs:.4e}"
            + f"; report dt_diff match {hit[0]:.4e} vs hand {dtd0:.4e}")
    print(f"PASS (a) diffusion slack: {info}")
    return info


# ==========================================================================
# (b) coupled CFL envelope vs the 1.32x pin (fixed-horizon bisection)
# ==========================================================================
_B_HORIZON = 4500       # divergence-detection horizon per probe (P3's upward cap)


def _probe_b(mult, und, nst=_B_HORIZON):
    out, rc = _run_child(_march_child_src(
        mult * und, nst, ["-fluidUpdate", "explicit"]))
    if _bounded(out, rc):
        return True, -1
    return False, _div_step(out)


def test_b_cfl_envelope():
    """PREDICTION (before assert): the explicit-fluid coupled boundary on the
    e74-matched column sits at 1.32x the DISCRETE undrained pencil (E7.4,
    both benchmark soils; band = the e72/e74 +-25% spread convention ->
    [0.99, 1.65]). Bisection is FIXED-HORIZON (4500 steps — the P3 upward-leg
    cap; slow-onset divergence just above the true boundary is absorbed by
    the band's upper end, honestly documented per §12 P3 item 3), endpoints:
    0.99x must survive the horizon (the P3 implicit-lane secular pumping at
    1.0x needed 7761 steps — margin ~1.7x; a faster explicit-lane pumping
    would fail here LOUDLY for MAIN adjudication), 2.0x must diverge inside
    it (P3 implicit lane: 1904 steps at 2.0x). A 0.8x march survives a
    3000-step horizon (implicit lane blew at 11789 there)."""
    und = _pencil("e72_exp")

    lo, hi = 0.99, 2.0
    ok_lo, _ = _probe_b(lo, und)
    assert ok_lo, (
        f"0.99x the undrained pencil DIVERGED within {_B_HORIZON} steps — the "
        f"boundary sits below the [0.99, 1.65] band (or the explicit-fluid "
        f"secular pumping is far faster than the implicit lane's 7761-step "
        f"1.0x horizon): MAIN must adjudicate")
    ok_hi, step_hi = _probe_b(hi, und)
    assert not ok_hi, (
        f"2.0x the undrained pencil BOUNDED for {_B_HORIZON} steps — no "
        f"upper divergence: the boundary claim (1.32x) is unmeasurable here")

    probes = [f"0.99x=B", f"2.0x=D@{step_hi}"]
    for _ in range(6):
        mid = 0.5 * (lo + hi)
        ok, st = _probe_b(mid, und)
        probes.append(f"{mid:.3f}x=" + ("B" if ok else f"D@{st}"))
        if ok:
            lo = mid
        else:
            hi = mid
    boundary = 0.5 * (lo + hi)
    assert 0.99 <= boundary <= 1.65, (
        f"measured explicit-fluid boundary {boundary:.3f}x the undrained "
        f"pencil outside [0.99, 1.65] (E7.4 pin 1.32x +-25%); probes: "
        f"{' '.join(probes)}")

    # 0.8x pencil full-horizon leg
    out, rc = _run_child(_march_child_src(
        0.8 * und, 3000, ["-fluidUpdate", "explicit"]))
    assert _bounded(out, rc), (
        f"0.8x the pencil did NOT march the 3000-step horizon stable:\n"
        + out[-1000:])

    print("=" * 74)
    print("LOUD (b) EXPLICIT-FLUID CFL RECORD (fixed-horizon, 4500 steps):")
    print(f"   boundary = {boundary:.3f}x the discrete undrained pencil "
          f"(pencil {und:.4e} s; E7.4 pin 1.32x, band [0.99, 1.65])")
    print(f"   probes: {' '.join(probes)}")
    print("   0.8x marched 3000 steps bounded")
    print("=" * 74)
    info = f"boundary {boundary:.3f}x in [0.99, 1.65]; {' '.join(probes)}"
    print(f"PASS (b) CFL envelope: {info}")
    return info


# ==========================================================================
# (c) implicit-lane equivalence (mutual convergence, order >= 0.9)
# ==========================================================================
_ZS84_PROBES_J = (15, 24, 27, 29)      # mid-column + the active drainage front


def test_c_lane_equivalence():
    """PREDICTION (before assert): FU_EXPLICIT (lumped-S* forward step) and
    FU_IMPLICIT + `-fsL zero` (backward-Euler fluid, L = 0) are BOTH O(Dt)
    splits of the same monolith, so their inter-lane p difference on the
    DRAINING ZS84 column shrinks ~linearly with Dt: matched-Dt runs at
    {0.4, 0.2, 0.1}x the advisory give diffs with observed order >= 0.9
    (log2 of successive ratios — the family's mutual-convergence
    methodology; exact equality is NOT expected at finite Dt)."""
    adv = _pencil("zs84_exp")
    t_end = 1.5
    n0 = int(8 * math.ceil(t_end / (0.4 * adv) / 8))
    if n0 > 1600:                       # runtime guard: shrink horizon, keep dt
        t_end = t_end * 1600.0 / n0
        n0 = 1600
    probes = [_nt(0, j) for j in _ZS84_PROBES_J]
    diffs = []
    for lev in range(3):
        n = n0 * 2 ** lev
        dt = t_end / n
        P_lane = {}
        for lane in ("explicit", "implicit_fsl0"):
            rp = _rec(f"c_{lane}{lev}.csv")
            _build_zs84_lane(lane, record=rp)
            _march_guard(n, dt, _nt(0, _Z["nely"] - 1))
            _t, P, t2c = _read_record(rp)
            assert len(P) == n, (
                f"{lane} level {lev}: {len(P)} record rows != {n} commits")
            rows = [k * n // 8 - 1 for k in range(1, 9)]
            cols = [t2c[p] for p in probes]
            P_lane[lane] = P[np.ix_(rows, cols)]
            assert np.all(np.isfinite(P_lane[lane])), (
                f"{lane} level {lev}: non-finite p")
        diffs.append(_relL2(P_lane["explicit"], P_lane["implicit_fsl0"]))
    orders = [math.log2(diffs[i] / diffs[i + 1]) for i in range(len(diffs) - 1)]
    assert diffs[-1] < diffs[0], f"inter-lane diffs did not shrink: {diffs}"
    assert min(orders) >= 0.9, (
        f"inter-lane mutual order < 0.9: orders={[round(o, 2) for o in orders]} "
        f"diffs={[f'{d:.3e}' for d in diffs]}")
    info = (f"advisory {adv:.4e}s dt0 {t_end / n0:.4e}s diffs "
            f"{[f'{d:.2e}' for d in diffs]} orders {[round(o, 2) for o in orders]}")
    print(f"PASS (c) lane equivalence: {info}")
    return info


# ==========================================================================
# (d) SMS composability — BOTH builders (THE gate)
# ==========================================================================
_SMS_INFO = "priced the UNDRAINED pencil"     # frozen 3b.4 INFO fragment
_SMS_INFO_TAG = "(ADR-73 P3b)"
_SMS_OLD = "UNSUPPORTED"                       # retired blanket warning marker
_D_KBAR = 1.0e-8       # DRAINING (pumping absorbed; diffusion slack still ~1e5)
_D_ALPHAM = 400.0      # ~2% xi at the column fundamental — quirks-row absorber
_D_HORIZON = 4000      # the P3 Richardson-battery step class


def _scaled_count(out):
    m = re.search(r"scaled\s+(\d+)\s*/\s*(\d+)", out)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)


def test_d_sms_composability():
    """PREDICTION (before assert): with dtTarget = 3x the undrained pencil on
    the e72-soil overlay column (drained/undrained per-element factor ~26x —
    the quirks-row scenario), DRAINED pricing sees every drained dt_e >>
    dtTarget and scales NOTHING (the no-overlay control measures exactly
    that), so a pre-fix march at the 'certified' dtTarget blows up fast
    (> 1.32x the true limit — the pathology child). With the P3b Kadd seam,
    BOTH builders price the UNDRAINED per-element pencil: all 8 overlay cells
    scale (loose gate >= 4 — direction pinned), the injected mass is large
    (s ~ T^2 = 9 per cell), the 4000-step march at dtTarget is STABLE
    (draining kbar=1e-8 + alphaM=400 absorb the lane-inherited secular
    pumping, which is gate (g)'s subject, not an SMS defect), the retired
    'UNSUPPORTED' text never prints, and the exact one-time INFO line
    'SMS sizing priced the UNDRAINED pencil for N overlay-owned elements
    (ADR-73 P3b)' does."""
    und = _pencil("e72_exp_k8")
    dtT = 3.0 * und
    ov = ["-fluidUpdate", "explicit"]

    # pathology leg: plain CDL at the would-be-certified dtTarget must diverge
    out_p, rc_p = _run_child(_march_child_src(
        dtT, 2500, ov, kbar=_D_KBAR, alpham=_D_ALPHAM))
    assert (rc_p != 0) or ("RESULT DIVERGED" in out_p), (
        f"pre-fix pathology NOT reproduced: an UNSCALED march at dtTarget = "
        f"3x the undrained pencil stayed bounded 2500 steps — gate (d) cannot "
        f"refute the old drained-priced certification:\n{out_p[-800:]}")
    step_p = _div_step(out_p)

    reports = [f"unscaled@dtT=D@{step_p}"]
    for name in ("CentralDifferenceSMS", "CentralDifferenceSMSConsistent"):
        out, rc = _run_child(_march_child_src(
            dtT, _D_HORIZON, ov, intg=(name, dtT, "-verbose"),
            kbar=_D_KBAR, alpham=_D_ALPHAM), timeout=420)
        assert _bounded(out, rc), (
            f"{name}: the SMS-certified dtTarget march DIVERGED "
            f"(step {_div_step(out)}) — undrained pricing did not make the "
            f"certification honest:\n{out[-1200:]}")
        assert _SMS_OLD not in out, (
            f"{name}: the retired blanket '...UNSUPPORTED...' warning still "
            f"prints (3b.4: warnIfOverlayPresentSMS must be REMOVED):\n"
            + out[-1200:])
        assert _SMS_INFO in out and _SMS_INFO_TAG in out, (
            f"{name}: the one-time INFO line 'SMS sizing priced the UNDRAINED "
            f"pencil for N overlay-owned elements (ADR-73 P3b)' is MISSING "
            f"(battery-greppable honesty, 3b.4):\n{out[-1200:]}")
        k, n = _scaled_count(out)
        assert k is not None, (
            f"{name}: no 'scaled k/n' report line found (-verbose):\n"
            + out[-1200:])
        assert k >= 4, (
            f"{name}: only {k}/{n} elements scaled at dtTarget = 3x the "
            f"undrained pencil — undrained pricing should catch (all) 8 "
            f"overlay cells (~26x factor class)")
        reports.append(f"{name}: scaled {k}/{n}, bounded@{_D_HORIZON}")

    # drained/no-overlay control: same dtTarget scales NOTHING, no INFO line
    out_c, rc_c = _run_child(_march_child_src(
        dtT, 500, [], intg=("CentralDifferenceSMS", dtT, "-verbose"),
        kbar=_D_KBAR, alpham=_D_ALPHAM, ov=False))
    assert _bounded(out_c, rc_c), (
        f"no-overlay control diverged (dtTarget is far below the drained "
        f"pencil — a solid-side regression):\n{out_c[-800:]}")
    kc, nc = _scaled_count(out_c)
    assert kc == 0, (
        f"no-overlay control scaled {kc}/{nc} elements at the same dtTarget — "
        f"the overlay-vs-control scaled-count DIRECTION (undrained pricing) "
        f"is not demonstrated")
    assert _SMS_INFO not in out_c, (
        "no-overlay control printed the overlay INFO line — the one-time line "
        "must fire only when sizing augmented >= 1 element")
    reports.append(f"control(no overlay): scaled {kc}/{nc}")

    info = "; ".join(reports) + f" [dtTarget {dtT:.4e} = 3x und {und:.4e}]"
    print(f"PASS (d) SMS composability: {info}")
    return info


# ==========================================================================
# (f) recorder CSV-twin under FU_EXPLICIT
# ==========================================================================
_TW_RTOL = 1e-12    # the recorder battery's pin (a) twin tolerance


def _attr(obj, name):
    v = obj.attrs[name]
    if isinstance(v, bytes):
        return v.decode("utf-8")
    if isinstance(v, np.ndarray):
        if v.dtype.kind == "S":
            d = [x.decode("utf-8") for x in v.reshape(-1)]
            return d[0] if v.size == 1 else d
        if v.size == 1:
            return v.reshape(-1)[0].item()
    return v


def _stages(f):
    return sorted(k for k in f if k.startswith("MODEL_STAGE"))


def _read_overlay_h5(f, tag):
    name = f"overlayPressure_{tag}"
    for st in _stages(f):
        g = f.get(f"{st}/RESULTS/ON_NODES/{name}")
        if g is not None and "DATA" in g:
            ids = np.asarray(g["ID"][...]).reshape(-1)
            data = np.asarray(g["DATA"][...])
            d2 = data[:, :, 0] if data.ndim == 3 else np.atleast_2d(data)
            comps = _attr(g, "COMPONENTS") if "COMPONENTS" in g.attrs else ""
            return ids, d2, comps
    return None


def _read_topology_h5(f, tag):
    for st in _stages(f):
        g = f.get(f"{st}/MODEL/OVERLAYS/OVERLAY_{tag}")
        if g is not None and "REGION_NODES" in g:
            rn = np.asarray(g["REGION_NODES"][...]).reshape(-1)
            dn = np.asarray(g["DRAINED_NODES"][...]).reshape(-1)
            return rn, dn
    return None


def test_f_recorder_csv_twin():
    """PREDICTION (before assert): the P4 recorder channels are lane-blind
    (they read getCommittedP() through the const seam), so under FU_EXPLICIT
    the `recorder ladruno -overlay` HDF5 DATA equals the overlay's own
    `-record` CSV row-for-row (<= 1e-12 rel — both read the same committed
    field at the same commits), ID order == CSV header == REGION_NODES, the
    topology rows (REGION_NODES/DRAINED_NODES) match the model as built, and
    the drained columns are identically 0."""
    if not _HAVE_H5PY:
        _skip(f"h5py unavailable after site re-add: {_H5_ERR}")
    ns = 25
    csv, h5p = _rec("f_twin.csv"), _rec("f_twin.ladruno")
    for p in (csv, h5p):
        if os.path.exists(p):
            os.remove(p)
    _build_e72(fluid="explicit", record=csv)
    ops.recorder("ladruno", h5p, "-overlay", 77)
    dt = 0.5 * _pencil("e72_exp")
    # NB: _pencil wiped; rebuild (pencil is memoized so this is the only rebuild)
    _build_e72(fluid="explicit", record=csv)
    ops.recorder("ladruno", h5p, "-overlay", 77)
    for _ in range(ns):
        assert ops.analyze(1, dt) == 0, "gate (f) march failed"
    ops.wipe()                                   # flush + close recorder files

    _t, cP, ct2c = _read_record(csv)
    hdr = sorted(ct2c, key=ct2c.get)
    with h5py.File(h5p, "r") as f:
        res = _read_overlay_h5(f, 77)
        topo = _read_topology_h5(f, 77)
    assert res is not None, "no RESULTS/ON_NODES/overlayPressure_77 group written"
    assert topo is not None, "no MODEL/OVERLAYS topology row written"
    ids, d2, comps = res
    rn, dn = topo

    assert cP.shape[0] == ns and d2.shape[0] == ns, (
        f"row counts CSV {cP.shape[0]} / DATA {d2.shape[0]} != commits {ns}")
    assert np.max(np.abs(cP)) > 0.0, "CSV pressures all zero (nothing to compare)"
    assert list(ids) == hdr, f"HDF5 ID order != CSV header ({list(ids)[:6]} vs {hdr[:6]})"
    assert list(rn) == hdr, f"REGION_NODES != CSV header ({list(rn)[:6]} vs {hdr[:6]})"
    drained_expect = {_nt(0, _E72["nely"]), _nt(1, _E72["nely"])}
    assert set(int(x) for x in dn) == drained_expect, (
        f"DRAINED_NODES {sorted(int(x) for x in dn)} != model drained "
        f"{sorted(drained_expect)}")
    rel = _relL2(d2, cP)
    assert rel <= _TW_RTOL, (
        f"HDF5 DATA vs -record CSV relL2 {rel:.2e} > {_TW_RTOL:.0e} under "
        f"FU_EXPLICIT — the recorder twin broke on this lane")
    col = {int(t): k for k, t in enumerate(ids)}
    worst = max(float(np.max(np.abs(d2[:, col[int(t)]]))) for t in dn)
    assert worst == 0.0, f"drained node p != 0 under FU_EXPLICIT (max {worst:.3e})"
    info = (f"DATA[{d2.shape[0]}x{d2.shape[1]}] == CSV relL2 {rel:.1e}; "
            f"ID==header==REGION_NODES ({len(ids)} nodes); drained rows 0")
    print(f"PASS (f) recorder twin: {info}")
    return info


# ==========================================================================
# (g) zero-drainage secular pumping under FU_EXPLICIT — MEASURED
# ==========================================================================
_G_LEGS = ((0.4, 60000), (0.5, 45000), (1.0, 15000))    # (mult, step cap)


def test_g_secular_pumping_measured():
    """PREDICTION (before assert; MEASUREMENT gate — the numbers, not a pass
    band, are the deliverable): the P3 quirks scenario (kbar = 1e-11 e72
    column, UNDAMPED) re-run with the fully-explicit fluid at 0.4x/0.5x/1.0x
    the pencil. Same pathology CLASS expected (frozen +Q.p forces + zero
    dissipation -> O(dt) secular pumping; implicit-fluid lane measured
    48344/30018/7761 steps-to-blowup) but the fluid-side dynamics differ —
    measure, never assume. ASSERTED (the honest core only): every leg either
    survives its cap or diverges at a parsed step; among diverged legs the
    horizon is strictly monotone DECREASING in dt; a smaller-dt leg must not
    die earlier than a larger-dt leg's cap-bounded survival allows. The
    RECORD block below is the quirks-row amendment payload."""
    und = _pencil("e72_exp")
    results = []                        # (mult, cap, diverged, step)
    for mult, cap in _G_LEGS:
        out, rc = _run_child(_march_child_src(
            mult * und, cap, ["-fluidUpdate", "explicit"]), timeout=570)
        if _bounded(out, rc):
            results.append((mult, cap, False, -1))
        else:
            st = _div_step(out)
            assert st >= 1, (
                f"{mult}x leg neither BOUNDED nor parseably DIVERGED:\n"
                + out[-800:])
            results.append((mult, cap, True, st))

    # honest coherence gates
    div = [(m, s) for m, c, d, s in results if d]
    for (m1, s1), (m2, s2) in zip(div, div[1:]):
        assert s1 > s2, (
            f"steps-to-blowup NOT monotone in dt: {m1}x died at {s1} vs "
            f"{m2}x at {s2} (secular growth ~ O(dt) predicts strictly "
            f"earlier blowup at larger dt)")
    for i, (m1, c1, d1, s1) in enumerate(results):
        for (m2, c2, d2, s2) in results[i + 1:]:
            if d1 and not d2 and c2 >= s1:
                raise AssertionError(
                    f"incoherent: {m1}x diverged at step {s1} but the larger "
                    f"{m2}x survived its {c2}-step cap (>= {s1}) — divergence "
                    f"must not be slower at larger dt")

    print("=" * 74)
    print("LOUD (g) SECULAR-PUMPING RECORD under FU_EXPLICIT "
          "(quirks-row amendment payload):")
    print(f"   zero-drainage undamped e72 column, pencil {und:.4e} s:")
    for mult, cap, d, st in results:
        print(f"   {mult}x pencil: "
              + (f"DIVERGED at step {st} (cap {cap})" if d
                 else f"SURVIVED the {cap}-step cap"))
    print("   implicit-fluid lane reference (P3 run-2): 0.4x:48344 "
          "0.5x:30018 1.0x:7761")
    print("   drained control fact (P3): same solid without the overlay "
          "bounded 30k steps")
    print("=" * 74)
    info = "; ".join(
        f"{m}x=" + (f"D@{s}" for _ in (0,)).__next__() if d else f"{m}x=B@{c}"
        for m, c, d, s in results)
    # (readable rebuild of the same string)
    info = "; ".join((f"{m}x=D@{s}" if d else f"{m}x=B@{c}")
                     for m, c, d, s in results)
    print(f"PASS (g) secular pumping measured: {info}")
    return info


# ==========================================================================
# smokes
# ==========================================================================
def test_s1_driver_refusal():
    """PREDICTION (before assert): LadrunoStaggeredAnalyze on a domain whose
    ONLY overlay is FU_EXPLICIT refuses loudly — the driven set comes up
    empty and the driver fatals (negative rc -> exit 3, or exception ->
    exit 7; the P2 e5 dual handling). Quiet success (exit 5) fails. The
    refusal message wording is 3b.B's; it must cite the lane
    (explicit/driven/empty/fluidUpdate — loose match)."""
    out, rc = _run_child(_REFUSAL_CH % dict(
        dist=_DIST, kbar=float(_E72["kbar"]),
        ovargs=["-fluidUpdate", "explicit"]))
    assert rc != 0, (
        f"FU_EXPLICIT driver child exited 0 (no fatal):\n{out[-1200:]}")
    assert "DRIVER_ACCEPTED" not in out, (
        "LadrunoStaggeredAnalyze RAN an FU_EXPLICIT overlay (must refuse — "
        f"3b.1: driven set empty, loud fatal):\n{out[-1200:]}")
    low = out.lower()
    assert any(t in low for t in ("explicit", "driven", "empty", "fluidupdate")), (
        f"refusal message does not cite the lane:\n{out[-1200:]}")
    mode = "rc<0" if rc == 3 else ("exception" if rc == 7 else f"exit {rc}")
    print(f"PASS (s1) driver refusal: FU_EXPLICIT overlay fatals loudly ({mode})")
    return f"refusal via {mode}"


def test_s2_fsl_inert_under_explicit():
    """PREDICTION (before assert): -fsL is INERT under FU_EXPLICIT (no L in
    the update), so marches with default -fsL, explicit `-fsL classic`, and
    `-fsL oedometric` are BIT-IDENTICAL (1e-12 band on the recorded p — the
    CSV is full-precision, the recorder battery's twin pin). The one-time
    inert NOTICE is asserted on the unambiguously NON-default oedometric
    child (3b.1: 'a non-default -fsL ... gets a one-time notice'); an
    explicitly-passed classic equals the parser default and may not notice —
    flagged as a contract ambiguity, march-identity still asserted."""
    dt = 0.5 * _pencil("e72_exp")
    ncom = 60
    traces = {}
    for label, fsl in (("default", None), ("classic", "classic"),
                       ("oedometric", "oedometric")):
        rp = _rec(f"s2_{label}.csv")
        _build_e72(fluid="explicit", fsl=fsl, record=rp)
        for s in range(ncom):
            assert ops.analyze(1, dt) == 0, f"{label} march failed at {s}"
        ops.wipe()
        _t, P, _c = _read_record(rp)
        assert P.shape[0] == ncom and np.all(np.isfinite(P))
        traces[label] = P
    scale = max(1.0, float(np.max(np.abs(traces["default"]))))
    for label in ("classic", "oedometric"):
        dmax = float(np.max(np.abs(traces[label] - traces["default"])))
        assert dmax <= 1e-12 * scale, (
            f"-fsL {label} march NOT bit-identical to the default under "
            f"FU_EXPLICIT (max |dp| = {dmax:.3e}, scale {scale:.3e}) — "
            f"L leaked into the explicit update")

    # one-time notice (non-default oedometric, fresh child)
    out, rc = _run_child(_NOTICE_CH % dict(
        dist=_DIST, kbar=float(_E72["kbar"]),
        ovargs=["-fluidUpdate", "explicit", "-fsL", "oedometric"],
        dt=float(dt)))
    assert rc == 0 and "BUILD_OK" in out, (
        f"-fsL oedometric + explicit child failed outright:\n{out[-1200:]}")
    assert "fsl" in out.lower(), (
        "no one-time notice mentioning -fsL observed for a NON-default -fsL "
        f"under -fluidUpdate explicit (3b.1 frozen surface):\n{out[-1200:]}")
    print("PASS (s2) -fsL inert: default/classic/oedometric marches "
          "bit-identical; oedometric notice observed")
    return "bit-identical (1e-12 band); notice on oedometric"


def test_s3_steady_then_transient():
    """PREDICTION (before assert): `-staticMode steady` stays legal under
    FU_EXPLICIT (3b.1: SM_STEADY keeps its CG solve of H p = f_seep at
    commits — a per-commit CG is not a march cost and the gravity-init
    recipe must keep working). SM_STEADY is a whole-overlay mode (the hook
    has no static-vs-transient detection), so this smoke runs a static
    LoadControl stage THEN a transient CDL march on the SAME overlay: both
    stages run (rc 0), the solid response stays finite, and the overlay
    keeps recording (steady branch records at every commit)."""
    dt = 0.5 * _pencil("e72_exp")
    rp = _rec("s3_steady.csv")
    _build_e72(fluid="explicit", static_steady=True, record=rp,
               series="Linear")
    # static gravity-init stage
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.algorithm("Linear")
    ops.integrator("LoadControl", 1.0 / 3.0)
    ops.analysis("Static")
    for s in range(3):
        assert ops.analyze(1) == 0, f"SM_STEADY static stage failed at step {s}"
    # transient march on the same steady overlay
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno")
    ops.analysis("Transient")
    for s in range(400):
        assert ops.analyze(1, dt) == 0, f"post-steady transient failed at {s}"
    u = ops.nodeDisp(_nt(0, _E72["nely"]), 2)
    ops.wipe()
    assert math.isfinite(u) and abs(u) < 1.0e3, f"transient response blew up: u={u}"
    _t, P, _c = _read_record(rp)
    assert P.shape[0] >= 400, (
        f"overlay stopped recording under SM_STEADY ({P.shape[0]} rows for "
        f"403 commits)")
    assert np.all(np.isfinite(P)), "non-finite recorded p under SM_STEADY"
    print(f"PASS (s3) SM_STEADY: static 3 + transient 400 commits ran; "
          f"{P.shape[0]} record rows, all finite")
    return f"{P.shape[0]} rows, finite"


def test_s4_db_roundtrip_midmarch():
    """PREDICTION (before assert): an FU_EXPLICIT overlay's DB save/restore
    mid-march resumes BIT-consistently — sLump_ is derived state rebuilt at
    the recvSelf -> snapshot path, and the explicit update reads only
    pCommitted_/dpCommitted_/uSnapshot_, all serialized (3b.2: zero format
    changes). Frozen-solid child (the framework battery's 1.D isolation:
    u == 0 exactly, so the round-trip depends ONLY on the overlay's fluid
    state; hydrostatic pInit drives an explicit diffusion march, dt = 0.1 <<
    dt_diff = 22.5 s): continuous 6 steps vs restart 3+3, DMAX == 0.0."""
    out, rc = _run_child(_DB_CH % dict(dist=_DIST))
    assert rc == 0, f"DB round-trip child failed (rc={rc}):\n{out[-1500:]}"
    m = re.search(r"DMAX\s+(\S+)\s+REL\s+(\S+)", out)
    assert m, f"child printed no DMAX line:\n{out[-800:]}"
    dmax, rel = float(m.group(1)), float(m.group(2))
    assert dmax == 0.0, (
        f"FU_EXPLICIT DB restart differs from continuous by {dmax:.3e} "
        f"(rel {rel:.2e}) — the explicit lane's fluid state did not "
        f"round-trip bit-exactly (sLump_ rebuild or p/dp serialization)")
    print("PASS (s4) DB round-trip: restart steps 4-6 bit-identical (DMAX 0.0)")
    return "restart bit-identical"


# ==========================================================================
# standalone runner (table-driven, the family pattern)
# ==========================================================================
_ALL = [
    ("a  diffusion slack + fold + report line", test_a_diffusion_slack),
    ("b  coupled CFL envelope (1.32x pin)", test_b_cfl_envelope),
    ("c  implicit-lane equivalence (order>=0.9)", test_c_lane_equivalence),
    ("d  SMS composability BOTH builders", test_d_sms_composability),
    ("f  recorder CSV-twin under FU_EXPLICIT", test_f_recorder_csv_twin),
    ("g  secular pumping MEASURED", test_g_secular_pumping_measured),
    ("s1 driver refusal", test_s1_driver_refusal),
    ("s2 -fsL inert + notice", test_s2_fsl_inert_under_explicit),
    ("s3 SM_STEADY then transient", test_s3_steady_then_transient),
    ("s4 DB round-trip mid-march", test_s4_db_roundtrip_midmarch),
]


def main():
    if not _BUILT:
        print(f"SKIP: worktree engine not built at {_DIST}")
        return 0
    import traceback
    fails, skips = [], []
    for name, fn in _ALL:
        try:
            info = fn()
            print(f"[PASS ] {name}" + (f"  -- {info}" if info else ""))
        except _Skip as exc:
            print(f"[SKIP ] {name}: {exc}")
            skips.append(name)
        except Exception as exc:  # noqa: BLE001
            print(f"[FAIL ] {name}: {exc}")
            traceback.print_exc()
            fails.append(name)
    print("=" * 74)
    if fails:
        print(f"FAILED {len(fails)}/{len(_ALL)} gates: {fails}"
              + (f"  (skipped {skips})" if skips else ""))
        return 1
    print(f"ALL {len(_ALL) - len(skips)}/{len(_ALL)} P3b EXPLICIT-FLUID GATES "
          f"PASSED (gate (e) MP transposed to documentation — see module "
          f"docstring)" + (f"  (skipped: {skips})" if skips else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
