"""WP-264: LadrunoConcrete3D's damage-driver gradients by dual numbers.

Compiles tests/_testbed/concrete3d_dual_check.cpp (header-only: the kernel
SRC/material/nD/LadrunoConcrete3DKernel.h + LadrunoDual.h; the replay testbed is
included for its fixture loader) with g++ and requires every gate to pass:

  * B: the templated double entry points (invariants, the four drivers, tcPhi,
    scalarDriver, dscalarDsig) are bit-identical to a verbatim transcription of
    the pre-WP-264 code, and every Dual pass's value equals the double result;
  * G: dscalarDsigDual equals a step-plateau central FD to 1e-6 off the drivers'
    kinks, in MPa, kPa and Pa, on random states, exactly repeated principal pairs
    (uniaxial, equal-biaxial, confined), near-repeated pairs and the replay
    record states; it is finite everywhere off the axis, independent of the
    eigenbasis inside a repeated pair, and equal to an eigen-free Dual<6> pass;
  * M: four break-the-gate mutants are flagged (shear factor 2, Lode derivative,
    eigenpair scrambling, the naive dual through acos at the meridians);
  * T: returnMap's damaged tangent with dual drivers vs FD drivers (same build,
    driverGradFD()): stress and returned State bit-identical, the dual tangent
    == the FD of the damaged map to 1e-6 on three damage-advancing paths in three
    unit systems, and nowhere worse than the FD-driver tangent on the records.
The timing section is skipped here (LADRUNO_SKIP_TIMING); run the binary by hand
for the per-call table. Skipped if g++ is unavailable.
"""
import os
import shutil
import subprocess

import pytest

from _cxx import cxx_runtime_flags  # WP-199: static GCC runtime on Windows

pytestmark = [pytest.mark.zone_a]

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, os.pardir))
TESTBED = os.path.join(HERE, "_testbed")


@pytest.mark.skipif(shutil.which("g++") is None, reason="g++ not available")
def test_concrete3d_dual_driver_gradients(tmp_path):
    gpp = shutil.which("g++")
    exe = os.path.join(tmp_path, "c3ddual.exe")
    src = os.path.join(TESTBED, "concrete3d_dual_check.cpp")
    subprocess.run([gpp, "-std=c++17", "-O2", "-I", REPO, src, "-o", exe, *cxx_runtime_flags()],
                   stdin=subprocess.DEVNULL, check=True, cwd=REPO, env=os.environ.copy())
    env = os.environ.copy()
    env["LADRUNO_SKIP_TIMING"] = "1"
    out = subprocess.run([exe, os.path.join(TESTBED, "concrete3d_replay_records.txt")], stdin=subprocess.DEVNULL,
                         cwd=REPO, capture_output=True, text=True, env=env)
    assert out.returncode == 0, f"dual check failed:\n{out.stdout}\n{out.stderr}"
    assert "DUAL CHECK: ALL PASS" in out.stdout
    assert out.stdout.count("  ok  : ") >= 25, out.stdout
