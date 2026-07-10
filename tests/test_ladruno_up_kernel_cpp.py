"""Cross-check the C++ shared Biot u-p (LadrunoUP) kernel WITHOUT building
OpenSees. Compiles tests/ladruno_up_kernel_check.cpp (pure — includes only the
two ADR-71 P0 headers, SRC/element/ladrunoUP/LadrunoUPKernel.h +
LadrunoUPShapes.h), runs it, and asserts it exits 0.

The C++ harness does the comparison itself: it (a) always runs standalone
self-tests — provider invariants (partition of unity, Sum(dN)=0, Sum(w*detJ)=
volume, detJ>0, Taylor-Hood pressure partition) plus the scalar kernel helpers
stabAlphaAuto / oneOverQbar / elemSizeH / dofMap against hand-computed contract
values — and (b) if the oracle-emitted cases file (tests/ladruno_up_cases.txt)
is present, drives every case through the pinned GP rule and diffs Q/H/S/H~/
f_seep + the DOF map against the stored values to <=1e-9.

This catches C++ transcription bugs in the coupled u-p block math — the Q/H/S
assembly, the direct-alpha stabilization block, the Taylor-Hood barycentric
pressure basis, and the heterogeneous-ndf DOF map — the ADR-71 P0 novel-math
surface. Skipped if no C++ compiler is present.

Independence (ADR-71 P0, load-bearing): this test never imports the numpy
oracle (ladruno_up_reference.py). The authoritative oracle-emit -> C++-check
cross-run is MAIN's WP0.D; here we only assert the harness's own gate passes.
"""
import os
import re
import shutil
import subprocess
import tempfile

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_UP = os.path.normpath(os.path.join(_HERE, "..", "SRC", "element", "ladrunoUP"))
_CASES = os.path.join(_HERE, "ladruno_up_cases.txt")
_MIN_CASES = 12  # plan contract: >= 12 canonical cases — guards a formatting
                 # glitch silently shrinking coverage (0.E-F4)


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_ladruno_up_kernel_selfcheck():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    src = os.path.join(_HERE, "ladruno_up_kernel_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "up_check.exe")
        # -std=c++17 is LOAD-BEARING: LadrunoUPShapes.h uses inline constexpr
        # static data members; c++11/14 compiles but fails to link.
        cc = subprocess.run(
            [gpp, "-O2", "-std=c++17", "-Wall", "-Wextra", "-I", _SRC_UP, src, "-o", exe],
            capture_output=True, text=True,
        )
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        # Pass the cases-file path explicitly; the harness self-tests unconditionally
        # and additionally validates the cases file when it exists.
        rr = subprocess.run([exe, _CASES], capture_output=True, text=True)
        assert rr.returncode == 0, f"kernel check failed:\n{rr.stdout}\n{rr.stderr}"
        assert "RESULT: PASS" in rr.stdout, rr.stdout
        m = re.search(r"\[cases\] parsed (\d+) case", rr.stdout)
        assert m, f"no parsed-case count in output:\n{rr.stdout}"
        assert int(m.group(1)) >= _MIN_CASES, (
            f"only {m.group(1)} cases parsed (>= {_MIN_CASES} required — "
            f"cases-file formatting glitch?):\n{rr.stdout}"
        )
