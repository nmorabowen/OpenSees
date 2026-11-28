"""WP-262: the extracted LadrunoUniaxialJ2 kernel and its dual-number DDM.

Compiles tests/ladruno_uniaxialj2_dual_check.cpp (pure: the kernel header
SRC/material/uniaxial/LadrunoUniaxialJ2Kernel.h + LadrunoDual.h) with g++ and
asserts, on six materials (Voce iso, two-term Chaboche, kinematic-only with
sig0 = 0, Lemaitre damage, Lemaitre under IMPL-EX, regularized s = 2 damage)
over a 12-step cyclic strain path:

  * BITS: the kernel in double equals a verbatim transcription of the pre-WP
    LadrunoUniaxialJ2::setTrialStrain bit for bit (stress, tangent, every trial
    state field). This is what makes the extraction safe for production;
  * TAN: the analytic consistent tangent equals the Dual<1> strain derivative;
  * SENS: the DDM stress sensitivity, propagated through committed history the
    way the material does it (conditional stress sensitivity, then
    commitSensitivity with the total strain gradient), equals a central FD of
    the whole path for E, sig0, Qinf, b and Hiso;
  * NOHIST (break-the-gate): the same DDM with the history sensitivities dropped
    misses the FD by orders of magnitude on the plastic paths.
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

CASES = {"iso_voce", "chaboche2", "kin_only", "lemaitre", "lemaitre_implex", "lemaitre_reg_s2"}
TAN_TOL = 1e-12     # |d sigma/d eps (dual) - E_alg| / E
SENS_TOL = 1e-5     # relative to the path's max |d sigma/d theta|; FD itself resolves ~1e-9
NOHIST_FLOOR = 1e-3


@pytest.fixture(scope="module")
def run():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_uniaxialj2_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "luj2d.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-Wall",
                             "-I", os.path.join(_SRC, "nD"), "-I", os.path.join(_SRC, "uniaxial"),
                             src, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    bits, sens, nohist = [], {}, {}
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "BITS":
            bits.append({"case": p[1], "step": int(p[2]), "ok": int(p[3]),
                         "plastic": int(p[5]), "damage": float(p[7]), "tan": float(p[9])})
        elif p[0] == "SENS":
            sens[(p[1], p[2])] = float(p[3])
        elif p[0] == "NOHIST":
            nohist[(p[1], p[2])] = float(p[3])
    return bits, sens, nohist


def test_kernel_is_bit_identical_to_the_pre_wp_code(run):
    bits, _, _ = run
    assert {b["case"] for b in bits} == CASES
    assert sum(b["plastic"] for b in bits) >= 50
    assert max(b["damage"] for b in bits) > 0.2, "damage paths too weak to exercise dD/deps"
    bad = [(b["case"], b["step"]) for b in bits if not b["ok"]]
    assert not bad, f"kernel differs from the pre-WP setTrialStrain at {bad}"


def test_analytic_tangent_equals_dual_tangent(run):
    bits, _, _ = run
    for b in bits:
        assert b["tan"] <= TAN_TOL, f"{b['case']} step {b['step']}: tangent off by {b['tan']:.3e}"


def test_ddm_with_history_matches_fd_of_the_whole_path(run):
    _, sens, _ = run
    assert len(sens) == 29          # 6 cases x 5 parameters, minus b on kin_only
    for (case, par), err in sens.items():
        assert err <= SENS_TOL, f"{case}/{par}: DDM vs FD of the path off by {err:.3e}"


def test_dropping_history_sensitivities_fails(run):
    """Break-the-gate: without d(history)/d(theta) the DDM is wrong on plastic
    paths. If this passes the comparison, the history leg is not being tested."""
    _, _, nohist = run
    for case in ("iso_voce", "chaboche2", "lemaitre"):
        for par in ("sig0", "E", "Hiso"):
            assert nohist[(case, par)] >= NOHIST_FLOOR, (
                f"{case}/{par}: DDM without history matched FD ({nohist[(case, par)]:.3e})")
