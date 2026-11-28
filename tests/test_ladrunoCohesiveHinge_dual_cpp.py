"""WP-267: the extracted LadrunoCohesiveHinge kernel and its dual-number DDM.

Compiles tests/ladruno_cohesivehinge_dual_check.cpp (pure: the kernel header
SRC/material/uniaxial/LadrunoCohesiveHingeKernel.h + LadrunoDual.h) with g++ and
asserts, on eight hinges (exponential and linear softening; default, user and
floor-clamped penalty; parabolic P-M and the -pmBal balance law in two of its
branches; IMPL-EX with the STEP and a PSEUDO step ratio) over a 12-step cyclic
jump path:

  * CALIB / BITS: the kernel in double equals a verbatim transcription of the
    pre-WP class (computeDerived, getPMMultiplier, envelope, setTrialStrain) bit
    for bit: calibration, P-M multiplier, moment, tangent, trial kappaMax;
  * TAN: the reported tangent equals the Dual<1> jump derivative;
  * SENS: the DDM moment sensitivity, propagated through committed history the
    way the material does it, equals a central FD of the whole path (with
    updateParameter's re-derivation) for Mc, Gf, Kpen and every live P-M constant;
  * NOHIST (break-the-gate): the same DDM with the history sensitivities
    (kappaMax, kappaMax_old) dropped misses the FD.
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

CASES = {"exp_default", "lin_user", "lin_clamped", "exp_pm", "lin_pmbal_mid",
         "exp_pmbal_hi", "exp_implex", "lin_pm_implex_pseudo"}
TAN_TOL = 1e-12     # |dM/dk (dual) - reported tangent| / Kpen
SENS_TOL = 1e-6     # relative to the path's max |dM/dtheta|; measured <= 2.1e-8
NOHIST_FLOOR = 1e-3  # measured >= 1.47e-3


@pytest.fixture(scope="module")
def run():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_cohesivehinge_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lchd.exe")
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
                         "soft": int(p[5]), "damage": float(p[7]), "tan": float(p[9])})
        elif p[0] == "SENS":
            sens[(p[1], p[2])] = float(p[3])
        elif p[0] == "NOHIST":
            nohist[(p[1], p[2])] = float(p[3])
    return calib, bits, sens, nohist


def test_kernel_is_bit_identical_to_the_pre_wp_code(run):
    calib, bits, _, _ = run
    assert set(calib) == CASES and all(calib.values()), f"calibration differs: {calib}"
    assert {b["case"] for b in bits} == CASES
    assert sum(b["soft"] for b in bits) >= 70, "too few softened steps"
    bad = [(b["case"], b["step"]) for b in bits if not b["ok"]]
    assert not bad, f"kernel differs from the pre-WP setTrialStrain at {bad}"


def test_reported_tangent_equals_dual_tangent(run):
    _, bits, _, _ = run
    for b in bits:
        assert b["tan"] <= TAN_TOL, f"{b['case']} step {b['step']}: tangent off by {b['tan']:.3e}"


def test_ddm_with_history_matches_fd_of_the_whole_path(run):
    _, _, sens, _ = run
    # 8 cases x (Mc, Gf, Kpen) + exp_pm/pm_implex (P0, alphaPM) + 2 pmBal x (P0, P0t, Pbal, muBal)
    assert len(sens) == 24 + 4 + 8
    for (case, par), err in sens.items():
        assert err <= SENS_TOL, f"{case}/{par}: DDM vs FD of the path off by {err:.3e}"


def test_dropping_history_sensitivities_fails(run):
    """Break-the-gate: without d(kappaMax)/d(theta) the DDM is wrong on the
    unload/reload and IMPL-EX steps. If this passes, the history leg is untested."""
    _, _, sens, nohist = run
    assert set(nohist) == set(sens)
    for key, err in nohist.items():
        assert err >= NOHIST_FLOOR, f"{key}: DDM without history matched FD ({err:.3e})"
