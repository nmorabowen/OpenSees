"""ADR-30 Phase 5 — ExplicitBathe / LNVD adoption of LadrunoProjection.

The projection consumer (P0-P3 on CentralDifferenceLadruno) is now implemented by the
Noh-Bathe explicit family: ExplicitBathe and ExplicitBatheLNVD each became a
LadrunoProjectionConsumer, and the 4 SMS / Consistent subclasses inherit it without
further code (all chain their domainChanged() to the base). Each base projects the solved
acceleration onto the constraint manifold at BOTH Noh-Bathe sub-steps (and the committed a0
at the starter), reading diag(M) from the Diagonal SOE once per stage — the same contract
the CentralDifference family already had.

  TP5-1  ExplicitBathe + equalDOF: tie u1==u2, v1==v2, a1==a2 to machine eps over the run.
  TP5-2  ExplicitBathe: LadrunoProjection == Transformation on the SAME integrator (core).
  TP5-3  ExplicitBatheLNVD (alpha=0): LadrunoProjection == Transformation (LNVD base path).
  TP5-4  ExplicitBatheLNVD (alpha=0.8): tie held to machine eps even WITH FLAC damping.
  TP5-5  ExplicitBatheSMS (inherited projection): tie held (mass-scaling subclass).
  TP5-6  ExplicitBatheLNVDSMSConsistent (inherited + refineAccel): tie held (project AFTER
         the consistent-mass refine).
  TP5-7  no-projection regression: plain ExplicitBathe with no LadrunoProjection oscillates
         normally (the projection code is gated on theProjector != 0, which stays 0).
  TP5-8  constraintTieForce recorder readback under ExplicitBathe (h5py).
"""
import glob
import math
import os

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

PROJ = "LadrunoProjection"


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


# --------------------------------------------------------------------------
# shared tied SHM model, parametrized over the integrator + constraint handler
# --------------------------------------------------------------------------
def _run_tied_shm(integrator_args, handler_args, m1, m2, k, d0, dt, nsteps,
                  record=None):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.equalDOF(1, 2, 1)                       # node1 retained, node2 constrained
    # compliant ICs: both tied DOFs start at d0 (equalDOF built at u=0 -> delta=0)
    ops.setNodeDisp(1, 1, d0, "-commit")
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints(*handler_args)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator_args)
    if record is not None:
        ops.recorder("ladruno", record, "-N", "constraintTieForce")
    ops.analysis("Transient")
    u1 = [ops.nodeDisp(1, 1)]
    tie = 0.0
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        u1.append(ops.nodeDisp(1, 1))
        tie = max(tie, abs(ops.nodeDisp(1, 1) - ops.nodeDisp(2, 1)),
                  abs(ops.nodeVel(1, 1) - ops.nodeVel(2, 1)),
                  abs(ops.nodeAccel(1, 1) - ops.nodeAccel(2, 1)))
    return np.array(u1), tie


_M1, _M2, _K, _D0 = 2.0, 3.0, 100.0, 1.0
_OMEGA = math.sqrt(_K / (_M1 + _M2))
_DT = 0.2 * (2.0 / _OMEGA)                       # well inside the Noh-Bathe stable range


# --------------------------------------------------------------------------
# TP5-1 — ExplicitBathe holds the equalDOF tie to machine eps
# --------------------------------------------------------------------------
def test_TP5_1_explicitBathe_tie_exact():
    u, tie = _run_tied_shm(("ExplicitBathe",), (PROJ,),
                           _M1, _M2, _K, _D0, _DT, 300)
    assert len(u) == 301, "run terminated early (instability?)"
    assert tie < 1e-12, f"tie not held to machine eps (err {tie:.3e})"
    assert np.max(np.abs(u)) > 0.5 * _D0, "model is not oscillating"


# --------------------------------------------------------------------------
# TP5-2 — ExplicitBathe: projection reproduces the eliminated (Transformation) reference
# --------------------------------------------------------------------------
def test_TP5_2_explicitBathe_matches_transformation():
    u_proj, _ = _run_tied_shm(("ExplicitBathe",), (PROJ,),
                              _M1, _M2, _K, _D0, _DT, 300)
    u_ref, _ = _run_tied_shm(("ExplicitBathe",), ("Transformation",),
                             _M1, _M2, _K, _D0, _DT, 300)
    n = min(len(u_proj), len(u_ref))
    assert n == 301, "a run terminated early (instability?)"
    diff = float(np.max(np.abs(u_proj[:n] - u_ref[:n])))
    scale = float(np.max(np.abs(u_ref[:n])))
    assert diff < 1e-9 * (1.0 + scale), (
        f"projection vs Transformation mismatch {diff:.3e} (scale {scale:.3e})")
    assert scale > 0.5 * _D0


# --------------------------------------------------------------------------
# TP5-3 — ExplicitBatheLNVD (undamped) reproduces the Transformation reference
# --------------------------------------------------------------------------
def test_TP5_3_explicitBatheLNVD_matches_transformation():
    integ = ("ExplicitBatheLNVD", 0.54, 0.0)     # alpha=0 -> identical to ExplicitBathe
    u_proj, _ = _run_tied_shm(integ, (PROJ,), _M1, _M2, _K, _D0, _DT, 300)
    u_ref, _ = _run_tied_shm(integ, ("Transformation",), _M1, _M2, _K, _D0, _DT, 300)
    n = min(len(u_proj), len(u_ref))
    assert n == 301, "a run terminated early (instability?)"
    diff = float(np.max(np.abs(u_proj[:n] - u_ref[:n])))
    scale = float(np.max(np.abs(u_ref[:n])))
    assert diff < 1e-9 * (1.0 + scale), (
        f"LNVD projection vs Transformation mismatch {diff:.3e} (scale {scale:.3e})")
    assert scale > 0.5 * _D0


# --------------------------------------------------------------------------
# TP5-4 — the tie holds to machine eps even WITH FLAC local damping active
# --------------------------------------------------------------------------
def test_TP5_4_explicitBatheLNVD_tie_exact_with_damping():
    u, tie = _run_tied_shm(("ExplicitBatheLNVD", 0.54, 0.8), (PROJ,),
                           _M1, _M2, _K, _D0, _DT, 300)
    assert len(u) == 301, "run terminated early (instability?)"
    assert tie < 1e-12, f"tie not held with FLAC damping (err {tie:.3e})"
    # FLAC damping drains the SHM toward rest — amplitude should decay below the start
    assert abs(u[-1]) < abs(u[0]) + 1e-9


# --------------------------------------------------------------------------
# TP5-5 — ExplicitBatheSMS inherits projection (mass-scaling subclass, no own newStep)
# --------------------------------------------------------------------------
def test_TP5_5_explicitBatheSMS_inherits_projection():
    integ = ("ExplicitBatheSMS", 0.54, _DT)      # $p $dtTarget
    u, tie = _run_tied_shm(integ, (PROJ,), _M1, _M2, _K, _D0, _DT, 300)
    assert len(u) == 301, "run terminated early (instability?)"
    assert tie < 1e-12, f"SMS tie not held to machine eps (err {tie:.3e})"
    assert np.max(np.abs(u)) > 0.5 * _D0


# --------------------------------------------------------------------------
# TP5-6 — ExplicitBatheLNVDSMSConsistent inherits projection AND refineAccel
#         (the project-AFTER-consistent-refine path)
# --------------------------------------------------------------------------
def test_TP5_6_explicitBatheLNVDSMSConsistent_inherits_projection():
    integ = ("ExplicitBatheLNVDSMSConsistent", 0.54, 0.0, _DT)   # $p $alpha $dtTarget
    u, tie = _run_tied_shm(integ, (PROJ,), _M1, _M2, _K, _D0, _DT, 300)
    assert len(u) == 301, "run terminated early (instability?)"
    assert tie < 1e-12, f"consistent-SMS tie not held to machine eps (err {tie:.3e})"
    assert np.max(np.abs(u)) > 0.5 * _D0


# --------------------------------------------------------------------------
# TP5-7 — no-projection regression: plain ExplicitBathe is untouched (theProjector==0)
# --------------------------------------------------------------------------
def test_TP5_7_no_projection_regression():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, 1.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)    # omega = 10
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(1, 1, 1.0, "-commit")
    ops.constraints("Plain")                      # NO LadrunoProjection -> projector stays 0
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator("ExplicitBathe")
    ops.analysis("Transient")
    dt = 0.2 * (2.0 / 10.0)
    u = [1.0]
    for _ in range(400):
        assert ops.analyze(1, dt) == 0
        u.append(ops.nodeDisp(1, 1))
    u = np.array(u)
    # free SHM: stays bounded, oscillates (sign changes), amplitude ~ the IC
    assert np.all(np.isfinite(u))
    assert np.max(np.abs(u)) < 2.0
    assert np.min(u) < -0.5 and np.max(u) > 0.5


# --------------------------------------------------------------------------
# TP5-8 — constraintTieForce recorder works under ExplicitBathe (h5py readback)
# --------------------------------------------------------------------------
def test_TP5_8_tie_force_recorder_under_explicitBathe(tmp_path):
    h5py = pytest.importorskip("h5py")

    m1, m2, F = 2.0, 3.0, 10.0
    f_exact = F * m2 / (m1 + m2)                  # 6.0 — internal tie force (±)
    out = str(tmp_path / "tieforce_eb.ladruno")

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equalDOF(1, 2, 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, F)
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator("ExplicitBathe")
    ops.recorder("ladruno", out, "-N", "constraintTieForce")
    ops.analysis("Transient")
    for _ in range(50):
        assert ops.analyze(1, 1e-3) == 0
    ops.wipe()                                    # flush/close the HDF5

    files = glob.glob(out) or glob.glob(out.replace(".ladruno", "*.ladruno"))
    assert files, "no .ladruno output produced"
    assert os.path.getsize(files[0]) > 0

    found = [None]
    with h5py.File(files[0], "r") as h:
        def visit(name, obj):
            up = name.upper()
            if (isinstance(obj, h5py.Dataset) and "CONSTRAINT_TIE_FORCE" in up
                    and up.rstrip("/").endswith("DATA")):
                found[0] = np.asarray(obj[...], dtype=float)
        h.visititems(visit)

    assert found[0] is not None, "CONSTRAINT_TIE_FORCE/DATA not found in the .ladruno output"
    data = found[0]
    assert data.size, "tie-force DATA empty"
    # the recorded tie force settles to the analytic ±F*m2/M once the EB startup transient
    # (committed a0 = 0) has passed; check the LATER steps reach it.
    assert abs(np.max(np.abs(data)) - f_exact) < 1e-6, (
        f"recorded tie force max {np.max(np.abs(data))} != expected {f_exact}")
    # per step the two nodes are equal-and-opposite (sum over nodes ~ 0)
    assert np.max(np.abs(data.sum(axis=tuple(range(1, data.ndim))))) < 1e-9
