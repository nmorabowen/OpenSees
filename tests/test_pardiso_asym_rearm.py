"""PARDISO half-storage asymmetry guard must RE-ARM, not just sample startup.

`system Pardiso -matrixType 1|2` stores only the upper triangle, so a genuinely
unsymmetric element tangent has its lower half DISCARDED and the run solves the
reflected system — converged, plausible, wrong. `PARDISOGenLinSOE::addA` carries
a guard that detects exactly that, but it was BUDGETED to the first few tangent
assemblies of the sparsity pattern (armed in setSize, spent in zeroA). A model
that is symmetric while elastic and only turns unsymmetric LATER — a contact
pair closing at step 50, a non-associated flow rule after first yield
(DruckerPrager with rho_bar != rho here; ManzariDafalias measured at 4-6%) —
burned the whole budget elastic and was never checked again. The in-source
comment admitted it (TIMs defect report, item 2).

The guard now re-arms every ASYM_RESAMPLE_PERIOD (64) assemblies while
`asymWarned` is unlatched, so it samples startup AND periodically.

  test_late_onset_asymmetry_is_caught   non-associated DP, -matrixType 2:
                                        the warning fires, and the assembly
                                        number it reports is PAST the startup
                                        window — the discriminator, since the
                                        pre-change guard could only ever report
                                        from the first three assemblies.
  test_elastic_run_never_warns          same mesh, ElasticIsotropic,
                                        -matrixType 1: no false positive, over a
                                        run long enough that several re-arms
                                        certainly happened.
  test_full_storage_never_warns         the yielding deck on the default
                                        (-matrixType 0) full CSR: nothing is
                                        discarded, so the guard must stay quiet.

The yielding run is NOT asserted to converge: half-storing an unsymmetric
tangent is wrong by construction, so Newton is iterating against a reflected
matrix. The warning is the deliverable; the ramp stops at the first failed step.
"""
import os
import re
import sys

import pytest

# Same rationale as test_pardiso_solver.py: MKL reads this at its first init.
os.environ.setdefault("MKL_NUM_THREADS", "1")

from _testbed import ops  # noqa: E402

# PARDISO IS MKL: linked on the Windows/oneAPI build only.
pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="PARDISO requires MKL (Windows/oneAPI build)"),
]

# ---- deck ------------------------------------------------------------------
# Unit cube of stdBrick, base fixed, shear (+x) plus compression (-z) at the top
# face. Material constants are the PEER DruckerPrager example's, reused from
# test_beziertet10_unsym_tangent.py where this same kPa-scale cube is shown to
# pass through first yield in a MIXED elastic/plastic state (some Gauss points
# yielded, some not) rather than collapsing.
#
# rho_bar != rho is the whole point: non-associated flow makes the consistent
# tangent Cep unsymmetric, and only AFTER a Gauss point yields. While elastic
# the tangent is exactly symmetric, which is what starves the un-re-armed guard.
_K, _G = 27777.78, 9259.26          # -> E = 25000, nu = 0.35
_RHO = 0.398
_RHO_BAR = 0.1
_SIG_Y = 5.0                        # yields on this ramp
_HARD_H = 1000.0                    # linear isotropic hardening (theta = 1)

_NX = 2                             # 2x2x2 stdBrick, 27 nodes, 81 DOF
_H_TOTAL = 2.025                    # shear, spread over the 9 top nodes
_V_TOTAL = 3.0                      # compression (raises the DP margin)

# Enough load steps that the assembly counter crosses the 64-pass re-arm period
# well AFTER yield onset. Each Newton step costs 2-4 tangent assemblies, so 48
# steps is ~135 assemblies here. Calibrated on this box: the deck is elastic
# through the startup budget, yields partway up the ramp, and the guard reports
# at ASSEMBLY 64 — the first re-arm — with ~70 assemblies of run still to come.
# (A lighter load, e.g. 0.8x these, yields later still and is caught at 128; the
# assertions below accept any re-arm, not the specific one.)
_NSTEPS = 48

_ELASTIC_E = 9.0 * _K * _G / (3.0 * _K + _G)
_ELASTIC_NU = (3.0 * _K - 2.0 * _G) / (2.0 * (3.0 * _K + _G))


def _nid(i, j, k):
    return 1 + i + (_NX + 1) * (j + (_NX + 1) * k)


def _build(system_args, elastic=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)

    if elastic:
        ops.nDMaterial("ElasticIsotropic", 1, _ELASTIC_E, _ELASTIC_NU)
    else:
        ops.nDMaterial("DruckerPrager", 1, _K, _G, _SIG_Y, _RHO, _RHO_BAR,
                       0.0, 0.0, 0.0, 0.0, _HARD_H, 1.0, 0.0)

    h = 1.0 / _NX
    for k in range(_NX + 1):
        for j in range(_NX + 1):
            for i in range(_NX + 1):
                ops.node(_nid(i, j, k), i * h, j * h, k * h)
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
                            _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                            1)
                tag += 1

    ntop = (_NX + 1) * (_NX + 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(_NX + 1):
        for i in range(_NX + 1):
            ops.load(_nid(i, j, _NX), _H_TOTAL / ntop, 0.0, -_V_TOTAL / ntop)

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-8, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / _NSTEPS)
    ops.analysis("Static")


def _ramp(system_args, elastic=False, nsteps=_NSTEPS):
    """Run the ramp; return the number of steps that converged.

    Stops at the first failure instead of asserting: on the half-stored
    unsymmetric tangent Newton is solving the reflected system, so divergence
    after yield is an expected consequence of the very defect under test.
    """
    _build(system_args, elastic=elastic)
    done = 0
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        done += 1
    return done


def _drain(capfd):
    """Console output of this test; skip if the build was made without MKL."""
    cap = capfd.readouterr()
    text = cap.out + cap.err
    if "has no MKL" in text:
        pytest.skip("this build was configured without MKL/PARDISO")
    # A build that does not know the command at all falls back to ProfileSPD,
    # which would make every assertion below pass or fail for the wrong reason.
    assert "ProfileSPD" not in text, (
        "`system Pardiso` was not understood — OpenSees fell back to "
        "ProfileSPD, so nothing here tested PARDISO:\n" + text)
    return text


_ASYM = "UNSYMMETRIC"
_AT_ASSEMBLY = re.compile(r"TANGENT ASSEMBLY (\d+)")

# Mirrors ASYM_CHECK_BUDGET / ASYM_RESAMPLE_PERIOD in PARDISOGenLinSOE.cpp.
_STARTUP_BUDGET = 3
_RESAMPLE_PERIOD = 64


# ---- the regression --------------------------------------------------------

def test_late_onset_asymmetry_is_caught(capfd):
    """The point of the change: caught after the startup budget is long gone."""
    done = _ramp(["Pardiso", "-matrixType", 2])
    text = _drain(capfd)

    assert _ASYM in text, (
        f"non-associated DruckerPrager ran {done}/{_NSTEPS} steps under "
        "`-matrixType 2` (upper triangle only) and the half-storage guard never "
        "fired. Either the deck stayed elastic — check the load/sigma_y "
        "calibration — or the guard is still budgeted to startup only.\n" + text)

    m = _AT_ASSEMBLY.search(text)
    assert m is not None, (
        "the warning must report the assembly pass it fired on, so a user can "
        "correlate it with the load step where their material yielded\n" + text)
    at = int(m.group(1))
    assert at > _STARTUP_BUDGET, (
        f"detected at assembly {at}, inside the startup window — this deck is "
        "no longer testing LATE-onset asymmetry (it would have been caught "
        "before this change too)")
    assert at >= _RESAMPLE_PERIOD, (
        f"detected at assembly {at}: past the startup budget but before the "
        f"first re-arm at {_RESAMPLE_PERIOD}, which is unreachable by design — "
        "the guard's arming logic has drifted from what this test asserts")

    # The warning must stay a once-per-SOE latch, not a per-assembly flood.
    assert text.count(_ASYM) == 1, (
        f"asymmetry reported {text.count(_ASYM)} times; the `asymWarned` latch "
        "must survive the re-arm")


def test_elastic_run_never_warns(capfd):
    """No false positive over a run several re-arm periods long."""
    done = _ramp(["Pardiso", "-matrixType", 1], elastic=True)
    text = _drain(capfd)
    assert done == _NSTEPS, (
        f"the elastic control must converge all {_NSTEPS} steps (got {done}); "
        "a short run would make the no-warning assertion vacuous — fewer than "
        f"{_RESAMPLE_PERIOD} assemblies means no re-arm ever happened and the "
        "test would pass without exercising anything\n" + text)
    assert _ASYM not in text, (
        "an ElasticIsotropic tangent is symmetric to round-off; the re-armed "
        "guard must not manufacture a warning\n" + text)


def test_full_storage_never_warns(capfd):
    """matrixType 0 keeps the full CSR — nothing discarded, nothing to warn."""
    done = _ramp(["Pardiso"])
    text = _drain(capfd)
    assert done > 0, "the full-storage yielding deck did not run a single step"
    assert _ASYM not in text, (
        "the guard is gated on half-storage; warning on the full unsymmetric "
        "CSR would be a false alarm on the very configuration users are told "
        "to switch to\n" + text)
