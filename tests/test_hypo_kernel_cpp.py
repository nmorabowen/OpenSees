"""Cross-check the C++ `-geom hypo` kernel WITHOUT building OpenSees (ADR 79
P0). Compiles tests/hypo_kernel_check.cpp (pure — includes only
SRC/element/solidTransformation/LadrunoHypoKernel.h), runs it, and asserts it
exits 0.

The C++ harness does the comparison itself: it (a) always runs standalone
self-tests — the kernel contract invariants (rigid-rotation increment gives
deps EXACTLY 0 at any per-step angle, the Hughes-Winget midpoint property;
polar orthogonality; the pure-stretch midpoint closed form 2(l-1)/(l+1);
bond6(I) identity; isotropic-tangent rotation invariance, the sharpest check
of the engineering-shear Voigt convention) — and (b) diffs every case in the
oracle-emitted tests/hypo_cases.txt (from tests/hypo_reference.py, SVD-based
polar — an independent algorithm from the kernel's Higham iteration) to
<= 1e-9: deps / R1 / J1 / stress push / tangent push.

This catches C++ transcription bugs in the objective-increment math — the
midpoint gradients, the de-rotation into the unrotated (Green-Naghdi) frame,
and the Voigt bond rotation — the ADR-79 P0 novel-math surface. Skipped if no
C++ compiler is present.
"""
import os
import re
import shutil
import subprocess
import tempfile

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_SRC_GEOM = os.path.normpath(
    os.path.join(_HERE, "..", "SRC", "element", "solidTransformation"))
_CASES = os.path.join(_HERE, "hypo_cases.txt")
_MIN_CASES = 16  # tests/hypo_reference.py N_CASES — guards a formatting glitch
                 # silently shrinking coverage


@pytest.mark.zone_a
@pytest.mark.t1
def test_cpp_hypo_kernel_selfcheck():
    gpp = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
    if gpp is None:
        pytest.skip("no C++ compiler available")
    assert os.path.exists(_CASES), (
        "tests/hypo_cases.txt missing — regenerate with "
        "`py -3.12 tests/hypo_reference.py`")
    src = os.path.join(_HERE, "hypo_kernel_check.cpp")
    with tempfile.TemporaryDirectory() as td:
        exe = os.path.join(td, "hypo_check.exe")
        cc = subprocess.run(
            [gpp, "-O2", "-std=c++17", "-Wall", "-Wextra", "-I", _SRC_GEOM,
             src, "-o", exe],
            capture_output=True, text=True, encoding="utf-8", errors="replace",
        )
        assert cc.returncode == 0, f"g++ failed:\n{cc.stderr}"
        rr = subprocess.run([exe, _CASES], capture_output=True, text=True,
                            encoding="utf-8", errors="replace")
        assert rr.returncode == 0, f"kernel check failed:\n{rr.stdout}\n{rr.stderr}"
        assert "RESULT: PASS" in rr.stdout, rr.stdout
        m = re.search(r"\[cases\] parsed (\d+) case", rr.stdout)
        assert m, f"no parsed-case count in output:\n{rr.stdout}"
        assert int(m.group(1)) >= _MIN_CASES, (
            f"only {m.group(1)} cases parsed (>= {_MIN_CASES} required):\n{rr.stdout}"
        )
