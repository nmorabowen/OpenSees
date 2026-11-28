"""WP-262: dual-number tangent oracle for the LadrunoJ2 return-map kernel.

Compiles tests/ladrunoj2_dual_check.cpp (pure: includes only
SRC/material/nD/LadrunoJ2Kernel.h, which pulls in LadrunoDual.h) with g++, and
runs the templated return map on twelve strain paths twice: once in double (the
production path, with the hand-derived analytic tangent Dtan) and once in
ladruno_dual::Dual<6> with the six engineering strain directions seeded, which
gives the exact algorithmic tangent with no step size. Seven paths drive the
undamaged map (returnMap); five drive the Lemaitre layer (returnMapDamaged),
whose tangent carries the hand-derived -sigma~ (x) dD/deps term, one of them
under IMPL-EX (D~ frozen), where the tangent must be (1 - D~) Dtan_eff.

Asserts, at every step of every path:

  * the analytic Dtan equals the dual tangent to round-off (<= 1e-13 relative
    to max|Dtan|; measured ~1e-16, worst 2.8e-14 on the strong-Voce path,
    including the damaged tangent up to D = 0.85). The
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
             "linkin_cf", "chaboche_nonprop", "voce_strong_3d",
             # Lemaitre damage layer (returnMapDamaged, hand-derived dD/deps term)
             "lemaitre_uniaxial", "lemaitre_3d_mixed", "lemaitre_s2", "lemaitre_threshold",
             # IMPL-EX: D~ frozen at the committed damage => (1-D~) Dtan_eff identity
             "lemaitre_implex_3d"}
DAMAGED = {s for s in SCENARIOS if s.startswith("lemaitre_")}
TAN_TOL = 1e-13   # relative to max|Dtan|. Measured worst 2.8e-14 WITH the IFT step
                  # (voce_strong_3d); WITHOUT it 7.5e-13 -- so this tolerance is what
                  # gates the IFT step (WP-262 review finding 2; 1e-12 let it through).
SYM_FLOOR = 1e-6  # the symmetrized mutant must miss by at least this much


@pytest.fixture(scope="module")
def dual_run():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladrunoj2_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "lj2d.exe")
        # -ffp-contract=off: the bit-identity leg compares the Dual value part with
        # the double run; FMA contraction (e.g. -march=haswell) may fuse the two
        # paths differently (WP-262 review finding 5).
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-Wall", "-ffp-contract=off",
                             "-I", _SRC_ND, src, "-o", exe,
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
                        "sym": float(kv["SYM"]), "bits": int(kv["BITS"]),
                        "damage": float(kv.get("DAMAGE", "0"))})
    return out


def test_damage_paths_actually_damage(dual_run):
    """The damaged scenarios must reach real damage (s=2 to D > 0.5), or the
    dD/deps term of the consistent tangent is never exercised."""
    for name in DAMAGED:
        assert max(r["damage"] for r in dual_run[name]) > 1e-4, f"{name}: no damage grew"
    assert max(r["damage"] for r in dual_run["lemaitre_s2"]) > 0.5
    # the threshold path is undamaged before p crosses pD, damaged after
    thr = dual_run["lemaitre_threshold"]
    assert thr[0]["damage"] == 0.0 and thr[-1]["damage"] > 0.1


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
