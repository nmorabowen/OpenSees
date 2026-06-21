"""ADR-30 Phase 4b — prescribed motion (non-homogeneous SP / imposedMotion).

Before P4b the projection handler WARNED and silently dropped every non-homogeneous SP,
marking the DOF eqn=-1 with no prescribed-displacement enforcement — so a base-excitation
model run under `constraints LadrunoProjection` got the WRONG answer (supports acted
fixed). P4b adds a handler applyLoad() override that imposes the prescribed displacement on
the node each step (mirroring TransformationDOF_Group::enforceSPs); ImposedMotionSP supplies
vel/accel itself. The DOF stays eqn-excluded, so the integrator is unchanged — the free DOFs
feel the support through F_int and damping (DOF_Group::setNode* skips eqn<0, preserving the
imposed value).

  TP1   uniform base excitation via imposedMotion (ground accel) — LadrunoProjection+Diagonal
        MATCHES Transformation+FullGeneral on the SAME model, and the response is nontrivial
        (proving the prescribed motion is actually applied, not dropped).
  TP2   constant non-homogeneous SP (prescribed support settlement, disp-only) — projection
        trajectory matches the Transformation reference.
  TP3   a prescribed DOF that is ALSO an MP/EQ slave -> named refusal (overconstraint).
  TP4   a prescribed DOF used as an MP/EQ MASTER -> in P4b this was refused; P4c now SUPPORTS
        it via kinematic imposition (slave driven u_c = C u_master), so TP4 asserts it runs and
        the slave tracks the prescribed master (full Tier-2 coverage in test_adr30_projection_p4c).

TP1/TP2 are the correctness core; TP3 pins the prescribed-slave refusal boundary.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


def _rel(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    d = _rms(a) + _rms(b)
    return 0.0 if d == 0.0 else _rms(a - b) / (0.5 * d)


# --------------------------------------------------------------------------
# TP1 — base excitation via imposedMotion (ground motion on the support DOF)
# --------------------------------------------------------------------------
_M, _K = 1.0, 100.0          # omega = 10 rad/s, dt_cr = 2/omega = 0.2


def _run_base_excitation(handler, system, nsteps=1000, dt=0.005):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)         # support (prescribed by imposedMotion)
    ops.node(2, 0.0)         # mass, free
    ops.mass(2, _M)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    # ground motion: sinusoidal acceleration record on the support
    ops.timeSeries("Trig", 1, 0.0, 5.0, 0.4, "-factor", 2.0)
    ops.pattern("MultipleSupport", 1)
    ops.groundMotion(1, "Plain", "-accel", 1)
    ops.imposedMotion(1, 1, 1)
    ops.constraints(*handler)
    ops.numberer("Plain")
    ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    u2 = [ops.nodeDisp(2, 1)]
    u1 = [ops.nodeDisp(1, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None, None
        u2.append(ops.nodeDisp(2, 1))
        u1.append(ops.nodeDisp(1, 1))
    return np.array(u2), np.array(u1)


def test_TP1_base_excitation_matches_reference_and_is_applied():
    """The headline P4b result: ground motion imposed on a support under
    LadrunoProjection+Diagonal reproduces the Transformation+FullGeneral answer — proving the
    prescribed motion is actually enforced (not dropped as it was before P4b)."""
    proj_u2, proj_u1 = _run_base_excitation((PROJ,), ("Diagonal",))
    ref_u2, ref_u1 = _run_base_excitation(("Transformation",), ("FullGeneral",))
    assert proj_u2 is not None, "projection base-excitation run failed"
    assert ref_u2 is not None, "reference base-excitation run failed"
    # the prescribed motion must actually drive the structure (would be ~0 if dropped)
    assert _rms(proj_u2) > 1e-3, "structure response too small — prescribed motion not applied?"
    assert _rms(proj_u1) > 1e-3, "support not moving — prescribed displacement not imposed?"
    # the free-DOF response matches the dense Transformation reference
    assert _rel(proj_u2, ref_u2) < 1e-6, (
        f"projection base-excitation response != reference (rel {_rel(proj_u2, ref_u2):.3e})")
    # the support tracks the SAME prescribed displacement under both handlers
    assert _rel(proj_u1, ref_u1) < 1e-6, (
        f"support displacement != reference (rel {_rel(proj_u1, ref_u1):.3e})")


# --------------------------------------------------------------------------
# TP2 — constant non-homogeneous SP (prescribed support settlement, disp-only)
# --------------------------------------------------------------------------
def _run_settlement(handler, system, nsteps=800, dt=0.005, settle=0.5):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)         # support (prescribed constant displacement)
    ops.node(2, 0.0)         # mass, free
    ops.mass(2, _M)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(1, 1, settle, "-const")     # constant prescribed displacement on the support
    ops.constraints(*handler)
    ops.numberer("Plain")
    ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    u2 = [ops.nodeDisp(2, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        u2.append(ops.nodeDisp(2, 1))
    return np.array(u2)


def test_TP2_constant_prescribed_settlement_matches_reference():
    """A constant non-homogeneous SP (disp-only, no GroundMotion) imposes a support
    displacement; the undamped mass oscillates about the static-offset equilibrium. The
    projection trajectory matches the Transformation reference."""
    proj = _run_settlement((PROJ,), ("Diagonal",))
    ref = _run_settlement(("Transformation",), ("FullGeneral",))
    assert proj is not None, "projection settlement run failed"
    assert ref is not None, "reference settlement run failed"
    assert _rms(proj) > 1e-3, "no settlement response — prescribed disp not applied?"
    assert _rel(proj, ref) < 1e-6, (
        f"settlement response != reference (rel {_rel(proj, ref):.3e})")


# --------------------------------------------------------------------------
# TP3 / TP4 — scope boundary: prescribed DOF as a constraint slave / master is refused
# --------------------------------------------------------------------------
def _runs_ok(build):
    try:
        build()
        return ops.analyze(1, 1e-3) == 0
    except Exception:
        return False


def test_TP3_prescribed_slave_refused():
    """A DOF carrying a prescribed motion (non-homog SP) that is ALSO an equalDOF slave is an
    overconstraint -> refused with a named error (cannot be both prescribed and tie-driven)."""
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)
        ops.node(2, 1.0); ops.mass(2, _M)
        ops.equalDOF(1, 2, 1)                 # node 2 dof 1 is the slave
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, 0.3, "-const")           # ... and ALSO prescribed -> overconstraint
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert not _runs_ok(build), "prescribed-on-slave should be refused"


def test_TP4_prescribed_master_now_supported_by_p4c():
    """A prescribed DOF used as an equalDOF MASTER was refused in P4b (would silently force the
    slaves' accel to 0). P4c supports it via KINEMATIC imposition: the slave is excluded from the
    equation set and driven u_c = C u_master each step, tracking the prescribed master exactly.
    (Full Tier-2 coverage lives in tests/test_adr30_projection_p4c.py.)"""
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)     # prescribed master
        ops.node(2, 0.0); ops.mass(2, _M)     # equalDOF slave -> derived (kinematic)
        ops.node(3, 0.0); ops.mass(3, _M)     # free downstream DOF (so the analysis has an equation)
        ops.equalDOF(1, 2, 1)                 # node 1 dof 1 is the master (retained)
        ops.uniaxialMaterial("Elastic", 1, _K)
        ops.element("zeroLength", 1, 2, 3, "-mat", 1, "-dir", 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(1, 1, 0.3, "-const")           # ... and prescribed -> prescribed master
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert _runs_ok(build), "prescribed master should now run (P4c kinematic imposition)"
    # the slave (node 2) tracks the prescribed master (node 1 = 0.3) exactly
    assert abs(ops.nodeDisp(2, 1) - 0.3) < 1e-10, "slave did not track the prescribed master"
