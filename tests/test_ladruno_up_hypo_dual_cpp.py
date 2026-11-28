"""WP-266 (D-UP1): the LadrunoUP ``-geom hypo`` element Jacobian by Dual<32>,
checked WITHOUT building OpenSees.

Compiles tests/ladruno_up_hypo_dual_check.cpp (pure: it includes only
SRC/element/ladrunoUP/LadrunoUPHypoKernel.h, LadrunoUPShapes.h,
SRC/element/solidTransformation/LadrunoHypoKernel.h and
SRC/material/nD/LadrunoJ2Kernel.h + LadrunoDual.h), runs it and asserts every
gate it prints:

  A  bit-identity of the double path: the templated LadrunoHypoKernel.h
     functions (shared with LadrunoBrick and the MPM engine) and the new
     LadrunoUP hypo call sequence against a verbatim copy of the pre-WP code;
  B  polar3<Dual<9>> against the closed-form polar derivative;
  C  the Dual<32> element Jacobian = a central FD of the residual on a step
     plateau, on 7 states (large strain, compacted near / at the Kozeny-Carman
     floor, rotated, small strain, J2-plastic, fast-draining);
  D  the dropped-block table (what the shipped dominant-term tangent leaves out);
  E  break-the-gate: the shipped tangent and two dependency-dropping mutants
     must fail the oracle; the incremental coupling is exactly zero on a rigid
     increment where the velocity form is not.

Skipped if no C++ compiler is present.
"""
import os
import re
import shutil
import subprocess
import tempfile

import pytest

from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, ".."))
_INC = [os.path.join(_ROOT, "SRC", "element", "ladrunoUP"),
        os.path.join(_ROOT, "SRC", "element", "solidTransformation"),
        os.path.join(_ROOT, "SRC", "material", "nD")]
_STATES = ("large_strain", "compacted_kc", "kc_floor", "rotated", "small_strain",
           "j2_plastic", "fast_drain")


@pytest.fixture(scope="module")
def run_check():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_up_hypo_dual_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "uphd.exe")
        cmd = [gpp, "-O2", "-std=c++17", "-Wall", "-Wextra"]
        for d in _INC:
            cmd += ["-I", d]
        cc = subprocess.run(cmd + [src, "-o", exe, *cxx_runtime_flags()],
                            stdin=subprocess.DEVNULL, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=900)
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe], stdin=subprocess.DEVNULL, capture_output=True, text=True,
                            encoding="utf-8", errors="replace", timeout=600)
    return rr.returncode, rr.stdout + rr.stderr


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_up_hypo_dual_all_gates(run_check):
    rc, out = run_check
    fails = [ln for ln in out.splitlines() if ln.startswith("FAIL")]
    assert rc == 0 and "RESULT: PASS" in out and not fails, "\n".join(fails) + "\n\n" + out
    m = re.search(r"\[checks\] (\d+) gates", out)
    assert m and int(m.group(1)) >= 60, out


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_up_hypo_dual_bit_identity_legs(run_check):
    _, out = run_check
    a = [ln for ln in out.splitlines() if ln.startswith("ok  : A")]
    # A1: 6 kernel legs; A2: 5 per class state x 3 states; A3: 2 per class state x 3
    assert len(a) >= 6 + 15 + 6, out


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_up_hypo_dual_oracle_every_state(run_check):
    _, out = run_check
    for st in _STATES:
        m = re.search(r"ok  : C \[%s\] dual Jacobian == central FD \(plateau\): uu (\S+)\s+up (\S+)"
                      r"\s+pu (\S+)\s+pp (\S+)" % st, out)
        assert m, f"no oracle line for {st}:\n{out}"
        assert max(float(x) for x in m.groups()) <= 1e-8, m.group(0)
        assert re.search(r"ok  : C \[%s\] Dual<32> residual value part == double residual bit" % st,
                         out), st


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_up_hypo_dual_dropped_block_table(run_check):
    _, out = run_check
    rows = [ln for ln in out.splitlines() if ln.startswith("TABLE ") and " | " in ln]
    body = rows[1:]                                    # rows[0] is the header
    assert len(body) == 16, "\n".join(rows)
    vals = {}
    for ln in body:
        cells = [c.strip() for c in ln[len("TABLE "):].split(" | ")]
        vals[cells[0]] = [float(c) for c in cells[1:]]
        assert len(vals[cells[0]]) == len(_STATES), ln
    exact = [k for k in vals if "exact by design" in k]
    assert len(exact) == 2 and all(max(vals[k]) <= 1e-13 for k in exact), exact
    # the shipped tangent is genuinely inexact at large strain (the reason WP-266 exists)
    total_uu = next(v for k, v in vals.items() if k.startswith("uu  total"))
    assert total_uu[_STATES.index("large_strain")] > 1e-2, total_uu
