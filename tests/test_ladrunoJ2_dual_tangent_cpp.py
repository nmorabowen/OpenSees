"""WP-262: dual-number tangent oracle for the LadrunoJ2 return-map kernel.

Compiles tests/ladrunoj2_dual_check.cpp (pure: includes only
SRC/material/nD/LadrunoJ2Kernel.h, which pulls in LadrunoDual.h) with g++, and
runs the templated return map on seven strain paths twice: once in double (the
production path, with the hand-derived analytic tangent Dtan) and once in
ladruno_dual::Dual<6> with the six engineering strain directions seeded, which
gives the exact algorithmic tangent with no step size.

Asserts, at every step of every path:

  * the analytic Dtan equals the dual tangent to round-off (<= 1e-12 relative
    to max|Dtan|; measured ~1e-16, worst 2.8e-14 on the strong-Voce path). The
    finite-difference gate this complements (test_ladrunoJ2_kernel_cpp.py)
    can only resolve ~1e-4;
  * the VALUE part of the Dual run equals the double run bit for bit (stress,
    equivalent plastic strain, dGamma, plastic strain, every backstress), which
    is what makes templating the kernel safe for production;

and, as the break-the-gate leg, that the oracle FLAGS a symmetrized Dtan on the
non-proportional Armstrong-Frederick path, where the algorithmic tangent is
genuinely unsymmetric. Skipped if no C++ compiler is available.
"""
import os
import shutil
import subprocess
import tempfile

import pytest

from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_ND = os.path.normpath(os.path.join(_HERE, "..", "SRC", "material", "nD"))

SCENARIOS = {"elastic_iso", "single_af_cyclic", "chaboche3_3d", "kin_linear_cyclic",
             "linkin_cf", "chaboche_nonprop", "voce_strong_3d"}
TAN_TOL = 1e-12   # relative to max|Dtan|
SYM_FLOOR = 1e-6  # the symmetrized mutant must miss by at least this much


@pytest.fixture(scope="module")
def dual_run():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladrunoj2_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lj2d.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-Wall", "-I", _SRC_ND, src, "-o", exe,
                             *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True)
        assert rr.returncode == 0, f"run failed:\n{rr.stderr}"
    out, cur = {}, None
    for line in rr.stdout.splitlines():
        p = line.split()
        if not p:
            continue
        if p[0] == "SCENARIO":
            cur = out.setdefault(p[1], [])
        elif p[0] == "STEP":
            kv = dict(zip(p[2::2], p[3::2]))
            cur.append({"step": int(p[1]), "plastic": int(kv["PLASTIC"]),
                        "tan": float(kv["TAN"]), "fd": float(kv["FD"]),
                        "sym": float(kv["SYM"]), "bits": int(kv["BITS"])})
    return out


def test_every_path_ran_and_yielded(dual_run):
    assert set(dual_run) == SCENARIOS
    assert all(dual_run[s] for s in SCENARIOS)
    assert sum(r["plastic"] for s in SCENARIOS for r in dual_run[s]) >= 30, \
        "too few plastic steps: the oracle would mostly check the elastic tangent"


def test_analytic_tangent_equals_dual_tangent(dual_run):
    for name, rows in dual_run.items():
        for r in rows:
            assert r["tan"] <= TAN_TOL, (
                f"{name} step {r['step']}: analytic Dtan differs from the dual tangent "
                f"by {r['tan']:.3e} (relative), above {TAN_TOL:.0e}")


def test_dual_value_part_is_bit_identical_to_double(dual_run):
    for name, rows in dual_run.items():
        for r in rows:
            assert r["bits"] == 1, (
                f"{name} step {r['step']}: the Dual run's value part differs from the "
                f"double run (templating changed the production arithmetic)")


def test_oracle_flags_a_symmetrized_tangent(dual_run):
    """Break-the-gate leg: on the non-proportional AF path the tangent is
    unsymmetric, so its symmetrized version must fail the oracle by far."""
    worst = max(r["sym"] for r in dual_run["chaboche_nonprop"] if r["plastic"])
    assert worst >= SYM_FLOOR, (
        f"symmetrized Dtan passes the oracle (worst {worst:.3e}): the gate cannot "
        f"see the Armstrong-Frederick cross-term")
    assert worst > 1e6 * TAN_TOL
