"""HRZ mass-conserving lumped mass (Ladruno ADR 35) — Zone-A integration tests.

The HRZ kernel math (conservation, reduce-to-rowsum on regular elements,
direction-dependent alpha + rotational mean, positivity guard) is verified
OpenSees-free in tests/_hrz_verify/hrz_standalone.cpp. These tests exercise the
INTEGRATION path: that `integrator CentralDifferenceLadruno ... -lump hrz`
parses, that the CTSLumping::HRZ branch builds dofDir from the element node
layout and calls Ladruno::hrzLump, and that the resulting per-element pencil
yields a sane critical time step on real elements.

Theory: Ladruno_implementation/35_ladruno_hrz_lumped_mass_adr.md.
"""
import math
import os
import shutil
import subprocess

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"


def _truss_bar_dtcr(lump):
    """1-element 2D truss (DEFAULT lumped mass => diagonal). All three lumps must
    agree here (the LEDGER quirk: -lump diagonal == rowsum == hrz on a diagonal
    getMass). Returns criticalTimeStep()."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0)
    ops.fix(2, 0, 1)  # single free DOF (x)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)  # cMass=0 default -> lumped
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", "-lump", lump)
    ops.analysis("Transient")
    ops.analyze(1, 1e-3)  # priming step triggers the dt_cr compute
    return ops.criticalTimeStep()


def _cantilever_beam_dtcr(lump):
    """1-element 3D elastic beam with CONSISTENT mass (-cMass): has rotational
    DOFs, so this is where HRZ actually matters (row-sum would zero the rotational
    mass). Returns criticalTimeStep()."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.node(2, 1.0, 0.0, 0.0)
    E, G, A, Jx, Iy, Iz, rho = 200.0, 80.0, 1.0, 1.0, 1.0, 1.0, 1.0
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.element("elasticBeamColumn", 1, 1, 2, A, E, G, Jx, Iy, Iz, 1,
                "-mass", rho, "-cMass")
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator(CDL, "-cfl", "-lump", lump)
    ops.analysis("Transient")
    ops.analyze(1, 1e-6)
    return ops.criticalTimeStep()


# --------------------------------------------------------------------------
# HRZ-1: -lump hrz parses and runs; dt_cr is finite & positive on a
#        consistent-mass beam (exercises dofDir + Ladruno::hrzLump end to end).
# --------------------------------------------------------------------------
def test_hrz_runs_on_consistent_beam():
    dt = _cantilever_beam_dtcr("hrz")
    assert math.isfinite(dt) and dt > 0.0, "HRZ dt_cr not finite/positive: %r" % dt


# --------------------------------------------------------------------------
# HRZ-2: on a DEFAULT (lumped/diagonal) truss, -lump hrz reduces EXACTLY to
#        -lump rowsum and -lump diagonal (non-breakage + the reduction claim,
#        restricted to the regular case in ADR 35 H3).
# --------------------------------------------------------------------------
def test_hrz_reduces_to_rowsum_on_lumped_truss():
    d = _truss_bar_dtcr("diagonal")
    r = _truss_bar_dtcr("rowsum")
    h = _truss_bar_dtcr("hrz")
    assert all(math.isfinite(x) and x > 0 for x in (d, r, h)), (d, r, h)
    assert abs(h - r) <= 1e-9 * (1.0 + abs(r)), "hrz=%.12g != rowsum=%.12g" % (h, r)
    assert abs(h - d) <= 1e-9 * (1.0 + abs(d)), "hrz=%.12g != diagonal=%.12g" % (h, d)


# --------------------------------------------------------------------------
# HRZ-3: on the consistent-mass beam, HRZ stays finite/positive where row-sum's
#        zero rotational mass makes the pencil degenerate. We only assert HRZ is
#        sane (the kernel math is proven in the standalone); row-sum may still
#        return a translational estimate, so we don't over-constrain equality.
# --------------------------------------------------------------------------
def test_hrz_well_posed_vs_rowsum_on_beam():
    h = _cantilever_beam_dtcr("hrz")
    assert math.isfinite(h) and h > 0.0, "HRZ degenerate on consistent beam: %r" % h


# --------------------------------------------------------------------------
# HRZ-4 (regression guard): on the CONSISTENT-mass beam, alpha != 1, so HRZ
#        MUST produce a different dt_cr than diagonal-of-consistent. This test
#        FAILS if hrzLump silently falls back to Diagonal (the bug the milestone
#        review flagged: the earlier "finite & positive" assertions could not).
# --------------------------------------------------------------------------
def test_hrz_differs_from_diagonal_on_consistent_beam():
    h = _cantilever_beam_dtcr("hrz")
    d = _cantilever_beam_dtcr("diagonal")
    assert math.isfinite(h) and math.isfinite(d) and h > 0 and d > 0, (h, d)
    rel = abs(h - d) / d
    assert rel > 1e-3, (
        "HRZ dt_cr=%.12g indistinguishable from diagonal=%.12g (rel=%.2e) "
        "-> HRZ likely silently fell back to Diagonal" % (h, d, rel)
    )


# --------------------------------------------------------------------------
# HRZ-5 (coverage gap): ExplicitBathe must ACCEPT -lump hrz (the milestone
#        review found its parser rejected hrz and its sendSelf dropped lumping).
#        Smoke: the integrator builds, runs a step, and reports a sane dt_cr.
# --------------------------------------------------------------------------
def test_explicitbathe_accepts_hrz():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 1.0, 0.0)
    ops.fix(2, 0, 1)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("Truss", 1, 1, 2, 1.0, 1, "-rho", 1.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    ops.integrator("ExplicitBathe", 0.54, "-cfl", "-lump", "hrz")
    ops.analysis("Transient")
    assert ops.analyze(1, 1e-3) == 0
    dt = ops.criticalTimeStep()
    assert math.isfinite(dt) and dt > 0.0, "ExplicitBathe -lump hrz dt_cr: %r" % dt


# --------------------------------------------------------------------------
# HRZ-6 (CI hook for the kernel oracle): compile + run the OpenSees-free
#        standalone self-test (tests/_hrz_verify/hrz_standalone.cpp). This is the
#        regression protection for the HRZ arithmetic (conservation, alpha!=1,
#        guards) flagged as untracked/un-CI'd by the milestone review. Skips
#        gracefully where g++ is unavailable.
# --------------------------------------------------------------------------
def test_hrz_standalone_kernel():
    gpp = shutil.which("g++")
    if gpp is None:
        pytest.skip("g++ not available")
    here = os.path.dirname(os.path.abspath(__file__))
    src = os.path.join(here, "_hrz_verify", "hrz_standalone.cpp")
    if not os.path.isfile(src):
        pytest.skip("standalone source not present")
    inc = os.path.join(here, "..", "SRC", "analysis", "integrator")
    exe = os.path.join(here, "_hrz_verify",
                       "hrz_standalone_ci" + (".exe" if os.name == "nt" else ""))
    comp = subprocess.run(
        [gpp, "-std=c++14", "-I", inc, "-DLADRUNO_MASSLUMPING_STANDALONE", src, "-o", exe],
        capture_output=True, text=True,
    )
    assert comp.returncode == 0, "standalone compile failed:\n" + comp.stderr
    run = subprocess.run([exe], capture_output=True, text=True)
    assert run.returncode == 0, "standalone HRZ self-test failed:\n" + run.stdout
