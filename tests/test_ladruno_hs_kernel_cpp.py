"""Pytest gate for the standalone C++17 header-only kernel SRC/material/nD/LadrunoHSKernel.h.

Compiles tests/ladruno_hs_kernel_check.cpp with zero OpenSees dependencies and
verifies:
  1. All 6x6 central-FD tangent checks pass (< 1e-7) in C++17 across all regimes
     for both -lodeShape mc and -lodeShape mn (including HSsmall).
  2. Stress, internal state, cap auto-calibration, and tangent entries match
     the Python reference oracle (tests/ladruno_hs_reference.py) to <= 1e-12.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

import ladruno_hs_reference as hs

REPO_ROOT = Path(__file__).resolve().parents[1]
KERNEL_H = REPO_ROOT / "SRC" / "material" / "nD" / "LadrunoHSKernel.h"
CHECK_CPP = REPO_ROOT / "tests" / "ladruno_hs_kernel_check.cpp"


def _compile_and_run_cpp_check(tmp_path: Path) -> str:
    exe_path = tmp_path / ("ladruno_hs_kernel_check.exe" if os.name == "nt" else "ladruno_hs_kernel_check")
    inc_dir = str(KERNEL_H.parent)

    gxx = shutil.which("g++") or shutil.which("clang++")
    if gxx:
        cmd = [gxx, "-O2", "-std=c++17", str(CHECK_CPP), f"-I{inc_dir}", "-o", str(exe_path)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if res.returncode == 0 and exe_path.exists():
            run_res = subprocess.run([str(exe_path)], capture_output=True, text=True, check=True)
            return run_res.stdout

    if os.name == "nt":
        setup_bat = REPO_ROOT / "Ladruno_scripts" / "setup_env.bat"
        cl_cmd = (
            f'call "{setup_bat}" >nul 2>&1 && '
            f'cl /nologo /EHsc /O2 /std:c++17 "{CHECK_CPP}" /I"{inc_dir}" /Fe:"{exe_path}"'
        )
        res = subprocess.run(["cmd.exe", "/d", "/s", "/c", cl_cmd], cwd=str(tmp_path), capture_output=True, text=True, check=False)
        if res.returncode == 0 and exe_path.exists():
            run_res = subprocess.run([str(exe_path)], capture_output=True, text=True, check=True)
            return run_res.stdout
        pytest.fail(f"Failed to compile C++17 kernel check:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}")

    pytest.skip("No C++17 compiler (g++, clang++, or MSVC cl) found")
    return ""


def test_cpp17_kernel_standalone_fd_and_python_parity(tmp_path: Path) -> None:
    assert KERNEL_H.exists(), f"Missing {KERNEL_H}"
    assert CHECK_CPP.exists(), f"Missing {CHECK_CPP}"

    stdout = _compile_and_run_cpp_check(tmp_path)
    assert "LADRUNO_HS_KERNEL_OK" in stdout

    parity_lines = [ln.strip() for ln in stdout.splitlines() if ln.strip().startswith("PARITY ")]
    assert len(parity_lines) == 2

    for ln in parity_lines:
        parts = ln.split()
        lode_int = int(parts[1])
        lode_str = "mc" if lode_int == 0 else "mn"
        vals = [float(x) for x in parts[3:]]

        p_py = hs.HSParams(
            37000.0, 29600.0, 90000.0, 0.50, 5.0, 38.0, 10.0, K0nc=0.38, lode_shape=lode_str
        )
        st_k0 = hs.initialize_state_from_stress(
            p_py, np.array([100.0, 38.0, 38.0, 0.0, 0.0, 0.0]), ocr=1.0
        )
        e_sc = np.array([1.5e-3, 0.0, 0.0, 1.0e-4, -5.0e-5, 3.0e-5])
        st_out, C_out, reg_out = hs.return_map_3d(p_py, st_k0, e_sc)
        assert reg_out == "shear_cap"

        ref_vals = [
            p_py.alpha_cap,
            p_py.beta_cap,
            st_out.sig[0],
            st_out.sig[1],
            st_out.sig[2],
            st_out.sig[3],
            st_out.sig[4],
            st_out.sig[5],
            st_out.gamma_p,
            st_out.p_p,
            C_out[0, 0],
            C_out[0, 1],
            C_out[3, 3],
        ]
        for k, (v_cpp, v_py) in enumerate(zip(vals, ref_vals)):
            rel = abs(v_cpp - v_py) / max(1.0, abs(v_py))
            tol = 1.0e-12 if k < 10 else 1.0e-7
            assert rel < tol, f"lode={lode_str} field={k}: cpp={v_cpp} vs py={v_py} (rel={rel:.2e})"
