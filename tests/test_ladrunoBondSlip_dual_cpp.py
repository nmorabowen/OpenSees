"""WP-267: the extracted LadrunoBondSlip kernel and its dual-number DDM.

Compiles tests/ladruno_bondslip_dual_check.cpp (pure: the kernel header
SRC/material/uniaxial/LadrunoBondSlipKernel.h + LadrunoDual.h) with g++ and
asserts, on four MC2010 backbones (default and user s0, a Gf override of s3, a
short plateau, a steep power law) over a 14-step cyclic slip path that visits
every branch (linear segment, power law, plateau, softening, residual, elastic
unload/reload, friction sliding at the cap):

  * CALIB / BITS: the kernel in double equals a verbatim transcription of the
    pre-WP class (setDerived, envelope, envelopeTangent, setTrialStrain) bit for
    bit: derived constants, stress, tangent, trial peak slip;
  * TAN: on the envelope and on elastic unload/reload the reported tangent is
    the Dual<1> slip derivative. On the plateau, the residual branch and while
    sliding the class reports k0 (or -1e-12 k0) instead of the consistent 0 to
    keep the pivot finite: the test pins that, so a change to it is seen;
  * SENS: the DDM bond-stress sensitivity, propagated through the committed
    state (Cslip, Cstress, CslipMaxAbs) the way the material does it, equals a
    central FD of the whole path (with updateParameter's re-derivation) for every
    live backbone constant, the power-law exponent alpha included;
  * NOHIST (break-the-gate): dropping the history sensitivities misses the FD.
Skipped if no C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import pytest

from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material"))

CASES = {"mc2010", "gf_override", "short_plateau", "steep_alpha"}
BR_ENVELOPE, BR_FLOORED, BR_UNLOAD, BR_SLIDE = 0, 1, 2, 3
TAN_TOL = 1e-12      # |d tau/d s (dual) - reported tangent| / k0
SENS_TOL = 1e-6      # relative to the path's max |d tau/d theta|; measured <= 1.2e-9
NOHIST_FLOOR = 1e-4  # measured >= 4.1e-4 (s2 on steep_alpha: softening-only parameter)


@pytest.fixture(scope="module")
def run():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_bondslip_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lbsd.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-Wall",
                             "-I", os.path.join(_SRC, "nD"), "-I", os.path.join(_SRC, "uniaxial"),
                             src, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    calib, bits, sens, nohist = {}, [], {}, {}
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "CALIB":
            calib[p[1]] = int(p[2])
        elif p[0] == "BITS":
            bits.append({"case": p[1], "step": int(p[2]), "ok": int(p[3]),
                         "branch": int(p[5]), "tan": float(p[7])})
        elif p[0] == "SENS":
            sens[(p[1], p[2])] = float(p[3])
        elif p[0] == "NOHIST":
            nohist[(p[1], p[2])] = float(p[3])
    return calib, bits, sens, nohist


def test_kernel_is_bit_identical_to_the_pre_wp_code(run):
    calib, bits, _, _ = run
    assert set(calib) == CASES and all(calib.values()), f"derived constants differ: {calib}"
    assert {b["case"] for b in bits} == CASES
    for br in (BR_ENVELOPE, BR_FLOORED, BR_UNLOAD, BR_SLIDE):
        assert sum(b["branch"] == br for b in bits) >= 4, f"branch {br} barely visited"
    bad = [(b["case"], b["step"]) for b in bits if not b["ok"]]
    assert not bad, f"kernel differs from the pre-WP setTrialStrain at {bad}"


def test_tangent_is_consistent_off_the_plateau_and_the_slider(run):
    _, bits, _, _ = run
    for b in bits:
        if b["branch"] in (BR_ENVELOPE, BR_UNLOAD):
            assert b["tan"] <= TAN_TOL, f"{b['case']} step {b['step']}: tangent off by {b['tan']:.3e}"
        else:
            # the deliberate pivot guard: k0 (relative error 1) or -1e-12 k0
            assert b["tan"] >= 0.5e-12, f"{b['case']} step {b['step']}: the pivot guard changed"


def test_ddm_with_history_matches_fd_of_the_whole_path(run):
    _, _, sens, _ = run
    assert len(sens) == 4 * 7          # 8 constants per case, minus s3 (Gf > 0) or Gf (Gf = 0)
    for (case, par), err in sens.items():
        assert err <= SENS_TOL, f"{case}/{par}: DDM vs FD of the path off by {err:.3e}"


def test_dropping_history_sensitivities_fails(run):
    """Break-the-gate: without d(Cslip, Cstress, CslipMaxAbs)/d(theta) the DDM
    is wrong on unload/reload and sliding. If this passes, the history is untested."""
    _, _, sens, nohist = run
    assert set(nohist) == set(sens)
    for key, err in nohist.items():
        assert err >= NOHIST_FLOOR, f"{key}: DDM without history matched FD ({err:.3e})"
    assert max(nohist.values()) > 0.5
