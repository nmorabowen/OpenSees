"""WP-221: the half-storage asymmetry guard, readable and optionally fatal.

`system Pardiso -matrixType 1|2` and `system ProfileSPD` store only the upper
triangle. On a genuinely unsymmetric tangent the lower half is DISCARDED and the
run solves the reflected system. The guard in addA (ADR-75 P1d, re-armed by TIMs
item 2) used to say so only as a one-time stderr WARNING, so a caller that
wanted to fall back to full storage had to scrape stderr for "UNSYMMETRIC"; a
rewording would have turned its fallback off silently. WP-221 adds:

  ladrunoSymmetryGuard [field]   the guard's state, from openseespy or Tcl
  -strictSymmetry                every assembly checked; once latched, every
                                 solve is refused, so analyze/eigen FAIL

Deck: the 2x2x2 stdBrick cube of test_pardiso_asym_rearm.py. Elastic it is
exactly symmetric; with non-associated DruckerPrager it turns unsymmetric only
once a Gauss point yields (late onset, after the 3-assembly startup budget).

Item 3 of the WP is checked by test_eigen_assembles_through_the_guard: `eigen`
assembles K through the active linear SOE's addA (ArpackSOE forwards to it), so
the guard's `checks` count moves on an eigen-only run. The eigen-strict tests
need ArpackSolver to stop on a failed linear solve, which it ignored before
WP-221.

Never assert on the warning TEXT here: the whole point is that callers no longer
depend on it.
"""
import os
import sys

import pytest

os.environ.setdefault("MKL_NUM_THREADS", "1")

from _testbed import ops  # noqa: E402

# ci-coverage: partial -- the Pardiso half needs MKL (the Windows/oneAPI build), like test_pardiso_asym_rearm.py; the ProfileSPD half runs anywhere but rides along in this file
pytestmark = [pytest.mark.zone_a]

needs_pardiso = pytest.mark.skipif(
    sys.platform != "win32", reason="PARDISO requires MKL (Windows/oneAPI build)")

FIELDS = ["latched", "worstDeviation", "scale", "assembly", "checkedAssemblies",
          "checks", "strict", "matrixType", "maxRelative", "refused"]

# ---- deck (same constants as test_pardiso_asym_rearm.py) --------------------
_K, _G = 27777.78, 9259.26
_RHO, _RHO_BAR = 0.398, 0.1
_SIG_Y, _HARD_H = 5.0, 1000.0
_NX = 2
_H_TOTAL, _V_TOTAL = 2.025, 3.0
_NSTEPS = 48
_E = 9.0 * _K * _G / (3.0 * _K + _G)
_NU = (3.0 * _K - 2.0 * _G) / (2.0 * (3.0 * _K + _G))
_MASS = 0.1


def _nid(i, j, k):
    return 1 + i + (_NX + 1) * (j + (_NX + 1) * k)


def _model(elastic):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    if elastic:
        ops.nDMaterial("ElasticIsotropic", 1, _E, _NU)
    else:
        ops.nDMaterial("DruckerPrager", 1, _K, _G, _SIG_Y, _RHO, _RHO_BAR,
                       0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)
    h = 1.0 / _NX
    for k in range(_NX + 1):
        for j in range(_NX + 1):
            for i in range(_NX + 1):
                n = _nid(i, j, k)
                ops.node(n, i * h, j * h, k * h)
                ops.mass(n, _MASS, _MASS, _MASS)
    for j in range(_NX + 1):
        for i in range(_NX + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    tag = 1
    for k in range(_NX):
        for j in range(_NX):
            for i in range(_NX):
                ops.element("stdBrick", tag,
                            _nid(i, j, k), _nid(i + 1, j, k),
                            _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                            _nid(i, j, k + 1), _nid(i + 1, j, k + 1),
                            _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1), 1)
                tag += 1
    ntop = (_NX + 1) * (_NX + 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(_NX + 1):
        for i in range(_NX + 1):
            ops.load(_nid(i, j, _NX), _H_TOTAL / ntop, 0.0, -_V_TOTAL / ntop)


def _static(system_args):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-8, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / _NSTEPS)
    ops.analysis("Static")


def _ramp(nsteps=_NSTEPS):
    done = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        done += 1
    return done


def _guard():
    rec = ops.ladrunoSymmetryGuard()
    assert len(rec) == len(FIELDS), rec
    g = dict(zip(FIELDS, rec))
    # every field is also its own subcommand, and agrees with the record
    for f in FIELDS:
        assert ops.ladrunoSymmetryGuard(f) == pytest.approx(g[f]), f
    return g


def _yielded_then(system_args):
    """Ramp the non-associated cube past first yield on the FULL CSR (nothing
    discarded, so it converges), then install `system_args` for what follows."""
    _model(elastic=False)
    _static(["FullGeneral"])
    done = _ramp()
    assert done >= _NSTEPS // 2, f"full-storage yielding ramp stalled at {done}"
    ops.system(*system_args)


PARDISO_STRICT = pytest.param(["Pardiso", "-matrixType", 2, "-strictSymmetry"],
                              marks=needs_pardiso, id="pardiso")
PROFILE_STRICT = pytest.param(["ProfileSPD", "-strictSymmetry"], id="profilespd")
PARDISO_SAMPLE = pytest.param(["Pardiso", "-matrixType", 2],
                              marks=needs_pardiso, id="pardiso")
PROFILE_SAMPLE = pytest.param(["ProfileSPD"], id="profilespd")


# ---- the query itself -------------------------------------------------------

def test_query_refuses_without_a_guarded_system():
    """No system, or one that does not half-store: an error, never a silent 0."""
    _model(elastic=True)
    with pytest.raises(Exception):
        ops.ladrunoSymmetryGuard()
    _static(["UmfPack"])
    with pytest.raises(Exception):
        ops.ladrunoSymmetryGuard()
    _static(["ProfileSPD"])
    with pytest.raises(Exception):
        ops.ladrunoSymmetryGuard("noSuchField")
    assert ops.ladrunoSymmetryGuard("tolerance") == pytest.approx(1e-8)


@needs_pardiso
def test_full_csr_reports_guard_off():
    """-matrixType 0 discards nothing: the guard reports matrixType 0 and never looks."""
    _model(elastic=False)
    _static(["Pardiso", "-strictSymmetry"])
    assert _ramp(6) == 6
    g = _guard()
    assert g["matrixType"] == 0 and g["strict"] == 1
    assert g["checks"] == 0 and g["latched"] == 0


# ---- symmetric model: no flag needed, no failure, query clean ----------------

@pytest.mark.parametrize("system_args", [PARDISO_SAMPLE, PROFILE_SAMPLE,
                                         PARDISO_STRICT, PROFILE_STRICT])
def test_symmetric_static_is_clean(system_args):
    _model(elastic=True)
    _static(system_args)
    assert _ramp() == _NSTEPS
    g = _guard()
    assert g["latched"] == 0 and g["assembly"] == -1 and g["refused"] == 0
    assert g["checks"] > 0 and g["checkedAssemblies"] > 0, g
    assert g["maxRelative"] < 1e-8, g
    if "-strictSymmetry" in system_args:
        # strict checks EVERY assembly; the sampler only ~3 of every 64
        assert g["strict"] == 1
        assert g["checkedAssemblies"] >= _NSTEPS, g


@pytest.mark.parametrize("system_args", [PARDISO_SAMPLE, PROFILE_SAMPLE,
                                         PARDISO_STRICT, PROFILE_STRICT])
@pytest.mark.parametrize("with_analysis", [True, False], ids=["analysis", "bare"])
def test_symmetric_eigen_is_clean(system_args, with_analysis):
    """Item 3: an eigen-only run assembles K through this SOE's addA.

    `bare` is the no-`analysis` path, where eigen builds a throwaway transient
    analysis around the user's `system` (OpenSeesCommands::eigen)."""
    _model(elastic=True)
    if with_analysis:
        _static(system_args)
    else:
        ops.system(*system_args)
    lam = ops.eigen(3)
    assert len(lam) == 3 and all(x > 0 for x in lam), lam
    g = _guard()
    assert g["checks"] > 0, (
        "eigen's stiffness assembly never reached the guard: the eigen does not "
        f"assemble through this system's addA ({g})")
    assert g["latched"] == 0 and g["refused"] == 0, g


# ---- unsymmetric model, static ----------------------------------------------

@pytest.mark.parametrize("system_args", [PARDISO_SAMPLE, PROFILE_SAMPLE])
def test_unsymmetric_static_latches_and_continues(system_args):
    """Default behaviour unchanged: the run goes on; the query now says why it is wrong."""
    _model(elastic=False)
    _static(system_args)
    _ramp()
    g = _guard()
    assert g["latched"] == 1, g
    assert g["strict"] == 0 and g["refused"] == 0
    assert g["worstDeviation"] > 1e-8 * g["scale"] > 0, g
    assert g["maxRelative"] >= g["worstDeviation"] / g["scale"] * (1 - 1e-12)
    # late onset: elastic through the startup budget, caught at a re-arm
    assert g["assembly"] >= 64, g


@pytest.mark.parametrize("system_args", [PARDISO_STRICT, PROFILE_STRICT])
def test_unsymmetric_static_strict_fails(system_args):
    _model(elastic=False)
    _static(system_args)
    done = _ramp()
    assert done < _NSTEPS, "strict mode let the half-stored unsymmetric run finish"
    g = _guard()
    assert g["latched"] == 1 and g["strict"] == 1 and g["refused"] >= 1, g
    # every assembly is checked, so the latch is at the FIRST unsymmetric
    # assembly: after the elastic startup, before the sampler's first re-arm
    assert 3 < g["assembly"] < 64, g
    # and it stays refused: no later step can sneak through on the bad system
    assert ops.analyze(1) != 0
    assert _guard()["refused"] > g["refused"]


# ---- unsymmetric model, eigen -----------------------------------------------

@pytest.mark.parametrize("system_args", [PARDISO_SAMPLE, PROFILE_SAMPLE])
def test_unsymmetric_eigen_latches(system_args):
    """The yielded tangent is checked by the eigen's own assembly (a fresh SOE,
    so the startup budget is armed). No assertion on the modes: they are of the
    reflected system, which is exactly what the latch reports."""
    _yielded_then(system_args)
    try:
        ops.eigen(3)
    except Exception:
        pass
    g = _guard()
    assert g["latched"] == 1 and g["refused"] == 0, g
    assert g["checkedAssemblies"] >= 1


@pytest.mark.parametrize("system_args", [PARDISO_STRICT, PROFILE_STRICT])
def test_unsymmetric_eigen_strict_fails(system_args):
    _yielded_then(system_args)
    with pytest.raises(Exception):
        ops.eigen(3)
    g = _guard()
    assert g["latched"] == 1 and g["strict"] == 1 and g["refused"] >= 1, g
