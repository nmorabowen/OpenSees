"""WP-265 (D-IMPLEX): the IMPL-EX tangent-identity gate for the templated kernels.

An IMPL-EX material evaluates, on the open step, an EXTRAPOLATED stress map whose
extrapolated quantities depend on the committed history only (the clock is
latched), and ships a tangent for it. That tangent must be the exact derivative of
the map it actually evaluates; otherwise Newton loses the constant-tangent
property IMPL-EX exists to buy, and nothing else in the suite would notice (a
wrong tangent still converges to the right stress, just slower or not at all).

Compiles tests/ladruno_implex_tangent_gate_check.cpp (OpenSees-free: it includes
only LadrunoJ2Kernel.h and LadrunoUniaxialJ2Kernel.h) with g++ and runs, on four
damage paths (two per material, each with an unload/reload or a reversal), the
production double map with -implex history threading exactly as the classes
commit it (D~ = D_n + dD_n; LadrunoJ2.cpp:366/463, LadrunoUniaxialJ2.cpp:304)
against the same map in ladruno_dual::Dual seeded on the strain.

Asserts at every step:
  * TAN  shipped IMPL-EX tangent == dual derivative of the IMPL-EX map
         (<= 1e-12 relative; measured <= 1.7e-13);
  * BITS the Dual run's value part equals the double run bit for bit;
  * DOUT the damage IMPL-EX returns for commit equals the fully implicit one bit
         for bit (IMPL-EX commits a true backward-Euler history);
and the two planted mutants (break-the-gate legs), each a "live d/d(eps)" in the
extrapolated map whose stress VALUE is unchanged:
  * LIVED the map degraded by the implicit damage D(eps): flagged on every step
          where the damage grows (measured min 1.3e-5 relative);
  * CLOCK an unlatched strain clock D~ = D_n + f(eps) dD_n: flagged on every step
          where the extrapolation is active (measured min 2.3e-6 relative).
The FD column is reported only: at the reload step that returns exactly to the
previous maximum it straddles the yield kink and misses by 1.3e-1, which is why
an FD cannot be this gate. Skipped if no C++ compiler is available.
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

SCENARIOS = {"j2_lemaitre_mixed", "j2_lemaitre_fast_uniaxial",
             "uj2_lemaitre_tension", "uj2_lemaitre_cyclic"}
TAN_TOL = 1e-12     # relative to max|K|
MUT_FLOOR = 1e-7    # a planted mutant must miss by at least this (1e5 x TAN_TOL)


@pytest.fixture(scope="module")
def gate_run():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_implex_tangent_gate_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "limx.exe")
        cc = subprocess.run([gpp, "-O2", "-std=c++17", "-Wall",
                             "-I", os.path.join(_SRC, "nD"), "-I", os.path.join(_SRC, "uniaxial"),
                             src, "-o", exe, *cxx_runtime_flags()],
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
            cur.append({"step": int(p[1]), "tan": float(kv["TAN"]), "fd": float(kv["FD"]),
                        "lived": float(kv["LIVED"]), "clock": float(kv["CLOCK"]),
                        "bits": int(kv["BITS"]), "dout": int(kv["DOUT"]),
                        "grow": int(kv["GROW"]), "extrap": float(kv["EXTRAP"]),
                        "damage": float(kv["DAMAGE"])})
    return out


def test_paths_exercise_the_extrapolation(gate_run):
    """Every path must damage, and the extrapolation D~ - D_n must be active on
    most steps, or the identity is only checked on the trivial D~ = D_n case."""
    assert set(gate_run) == SCENARIOS
    for name, rows in gate_run.items():
        assert sum(r["grow"] for r in rows) >= 5, f"{name}: too few damaging steps"
        assert sum(1 for r in rows if r["extrap"] > 0.0) >= 5, f"{name}: extrapolation idle"
    assert max(r["damage"] for r in gate_run["j2_lemaitre_fast_uniaxial"]) > 0.3
    assert max(r["damage"] for r in gate_run["uj2_lemaitre_cyclic"]) > 0.4


def test_implex_tangent_is_the_derivative_of_the_extrapolated_map(gate_run):
    for name, rows in gate_run.items():
        for r in rows:
            assert r["tan"] <= TAN_TOL, (
                f"{name} step {r['step']}: shipped IMPL-EX tangent differs from the dual "
                f"derivative of the extrapolated map by {r['tan']:.3e} (relative)")


def test_dual_value_part_is_bit_identical(gate_run):
    for name, rows in gate_run.items():
        for r in rows:
            assert r["bits"] == 1, f"{name} step {r['step']}: Dual value part != double run"


def test_implex_commits_the_implicit_damage(gate_run):
    for name, rows in gate_run.items():
        for r in rows:
            assert r["dout"] == 1, (
                f"{name} step {r['step']}: the damage IMPL-EX returns for commit is not the "
                f"backward-Euler one")


def test_gate_flags_a_live_damage(gate_run):
    """Break-the-gate leg 1: degrading by the implicit D(eps) instead of the frozen D~."""
    for name, rows in gate_run.items():
        grow = [r for r in rows if r["grow"]]
        worst = min(r["lived"] for r in grow)
        assert worst >= MUT_FLOOR, (
            f"{name}: a live-damage map passes the gate on a damaging step "
            f"(min miss {worst:.3e})")


def test_gate_flags_an_unlatched_clock(gate_run):
    """Break-the-gate leg 2: a strain clock read live at the trial, same stress value."""
    for name, rows in gate_run.items():
        active = [r for r in rows if r["extrap"] > 0.0]
        worst = min(r["clock"] for r in active)
        assert worst >= MUT_FLOOR, (
            f"{name}: an unlatched-clock map passes the gate where the extrapolation is "
            f"active (min miss {worst:.3e})")
