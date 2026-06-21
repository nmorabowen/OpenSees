"""ADR-30 Phase 4c — prescribed MASTER (Tier 2, KINEMATIC imposition).

P4b (shipped) handles prescribed motion on DOFs NOT in a constraint group and REFUSES a
prescribed DOF used as a constraint master. P4c supports a prescribed master driving its
slaves. An acceleration-only known-RHS projection was tried first but gives a non-converging
O(1e-3) displacement drift (the master's disp is externally imposed, the slave is leap-frog
integrated — they don't co-integrate). So a slave driven PURELY by prescribed master(s) is
instead KINEMATICALLY IMPOSED: excluded from the equation set and its u/v/a set from the
masters each step (after the masters' own disp is imposed), exactly like Transformation:

    u_c = sum_k C_k u_{m_k} + delta ;   v_c = sum_k C_k v_{m_k} ;   a_c = sum_k C_k a_{m_k}

  TC1   imposedMotion on an equalDOF MASTER -> the slave is kinematically driven and tracks
        the prescribed master to MACHINE PRECISION; a downstream free DOF matches
        Transformation+FullGeneral to machine precision (same support motion under both).
  TC2   constant non-homog SP on an equalDOF master -> slave holds the prescribed offset
        EXACTLY (a_p=0 would have left an accel-projection slave at 0 — kinematic gets it right).
  TC3   rigidLink -bar (2D), master ux prescribed (uy free + loaded): slave ux is a derived
        DOF (tracks the prescribed master ux exactly) while slave uy is a normal PROJECTED DOF
        (tied to the free master uy). Both ties hold to machine precision -- and projection
        gets this RIGHT where Transformation+CDL drops the uy tie (a Transformation quirk when
        a rigidLink master also carries an imposedMotion SP), so TC3 self-verifies the exact
        ties rather than comparing to an unreliable reference.
  TC4   a single slave tied to BOTH a free and a prescribed master (equationConstraint) ->
        MIXED -> refused with a named error.
  TC5   prescribed SLAVE still refused (overconstraint) — the P4b boundary that stays.
"""
import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"
_M, _K = 1.0, 100.0


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


def _rel(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    d = _rms(a) + _rms(b)
    return 0.0 if d == 0.0 else _rms(a - b) / (0.5 * d)


# --------------------------------------------------------------------------
# TC1 — imposedMotion on an equalDOF MASTER (slave kinematically driven)
# --------------------------------------------------------------------------
def _run_equaldof_master_gm(handler, system, nsteps=1000, dt=0.005):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)             # master: prescribed by imposedMotion (eqn-excluded)
    ops.node(2, 0.0); ops.mass(2, _M)   # slave: equalDOF-tied to the prescribed master
    ops.node(3, 0.0); ops.mass(3, _M)   # free downstream DOF, spring to the slave
    ops.equalDOF(1, 2, 1)               # node 1 (master) retained, node 2 (slave) constrained
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 2, 3, "-mat", 1, "-dir", 1)
    ops.timeSeries("Trig", 1, 0.0, 5.0, 0.4, "-factor", 2.0)
    ops.pattern("MultipleSupport", 1)
    ops.groundMotion(1, "Plain", "-accel", 1)
    ops.imposedMotion(1, 1, 1)
    ops.constraints(*handler); ops.numberer("Plain"); ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear"); ops.integrator(CDL)
    ops.analysis("Transient")
    u1 = [ops.nodeDisp(1, 1)]; u2 = [ops.nodeDisp(2, 1)]; u3 = [ops.nodeDisp(3, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        u1.append(ops.nodeDisp(1, 1)); u2.append(ops.nodeDisp(2, 1)); u3.append(ops.nodeDisp(3, 1))
    return np.array(u1), np.array(u2), np.array(u3)


def test_TC1_imposedMotion_on_equalDOF_master_drives_slave():
    proj = _run_equaldof_master_gm((PROJ,), ("Diagonal",))
    ref = _run_equaldof_master_gm(("Transformation",), ("FullGeneral",))
    assert proj is not None, "projection prescribed-master run failed"
    assert ref is not None, "reference run failed"
    p1, p2, p3 = proj
    r1, r2, r3 = ref
    assert _rms(p1) > 1e-3, "master not moving — prescribed motion not applied?"
    # the slave is kinematically driven by the prescribed master -> EXACT tracking (no drift)
    assert np.max(np.abs(p2 - p1)) < 1e-10, (
        f"slave does not track the prescribed master (max {np.max(np.abs(p2-p1)):.2e})")
    # the downstream free DOF sees the SAME support motion under both handlers -> machine match
    assert _rms(p3) > 1e-3, "downstream free DOF response too small to discriminate"
    assert _rel(p3, r3) < 1e-9, f"free-DOF response != reference (rel {_rel(p3,r3):.3e})"


# --------------------------------------------------------------------------
# TC2 — constant non-homogeneous SP on an equalDOF master
# --------------------------------------------------------------------------
def _run_equaldof_master_const(handler, system, nsteps=800, dt=0.005, settle=0.4):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0); ops.mass(2, _M)
    ops.node(3, 0.0); ops.mass(3, _M)
    ops.equalDOF(1, 2, 1)
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.element("zeroLength", 1, 2, 3, "-mat", 1, "-dir", 1)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
    ops.sp(1, 1, settle, "-const")          # prescribed constant displacement on the master
    ops.constraints(*handler); ops.numberer("Plain"); ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear"); ops.integrator(CDL)
    ops.analysis("Transient")
    u1 = [ops.nodeDisp(1, 1)]; u2 = [ops.nodeDisp(2, 1)]; u3 = [ops.nodeDisp(3, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        u1.append(ops.nodeDisp(1, 1)); u2.append(ops.nodeDisp(2, 1)); u3.append(ops.nodeDisp(3, 1))
    return np.array(u1), np.array(u2), np.array(u3)


def test_TC2_constant_SP_on_equalDOF_master():
    proj = _run_equaldof_master_const((PROJ,), ("Diagonal",))
    ref = _run_equaldof_master_const(("Transformation",), ("FullGeneral",))
    assert proj is not None and ref is not None
    p1, p2, p3 = proj
    r1, r2, r3 = ref
    # the slave tracks the prescribed constant master EXACTLY (kinematic — a_p=0 would have
    # frozen an acceleration-projection slave at its initial 0). Both the master and slave move
    # from the recorded initial 0 to the prescribed 0.4 at the first step, so compare slave->master.
    assert abs(p1[-1] - 0.4) < 1e-10, "master did not reach the prescribed constant value"
    assert np.max(np.abs(p2 - p1)) < 1e-10, "slave did not track the constant prescribed master"
    assert _rel(p3, r3) < 1e-9, f"free-DOF response != reference (rel {_rel(p3,r3):.3e})"


# --------------------------------------------------------------------------
# TC3 — rigidLink -bar (2D), master ux prescribed (uy free): derived ux + normal uy
# --------------------------------------------------------------------------
def _run_rigidlink_master(handler, system, nsteps=800, dt=0.005):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0); ops.mass(1, _M, _M)   # master: ux prescribed, uy free
    ops.node(2, 1.0, 0.0); ops.mass(2, _M, _M)   # slave via rigidLink -bar (ties ux,uy)
    ops.node(3, 1.0, 0.0); ops.mass(3, _M, _M)   # free DOF, spring to slave in y
    ops.uniaxialMaterial("Elastic", 1, _K)
    ops.node(4, 0.0, 0.0); ops.fix(4, 1, 1)
    ops.element("zeroLength", 1, 4, 1, "-mat", 1, "-dir", 2)   # master uy spring (keeps uy dynamic)
    ops.element("zeroLength", 2, 2, 3, "-mat", 1, "-dir", 2)   # slave uy -> node 3 uy
    ops.node(5, 1.0, 0.0); ops.fix(5, 1, 1)
    ops.element("zeroLength", 3, 5, 3, "-mat", 1, "-dir", 2)
    ops.rigidLink("bar", 1, 2)
    ops.timeSeries("Trig", 1, 0.0, 5.0, 0.4, "-factor", 2.0)
    ops.pattern("MultipleSupport", 1)
    ops.groundMotion(1, "Plain", "-accel", 1)
    ops.imposedMotion(1, 1, 1)                                 # prescribe master ux only
    # excite uy (decoupled from ux): a constant transverse load on the free master uy
    ops.timeSeries("Constant", 2); ops.pattern("Plain", 2, 2)
    ops.load(1, 0.0, 5.0)
    ops.constraints(*handler); ops.numberer("Plain"); ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear"); ops.integrator(CDL)
    ops.analysis("Transient")
    sx = [ops.nodeDisp(2, 1)]; sy = [ops.nodeDisp(2, 2)]
    m1x = [ops.nodeDisp(1, 1)]; m1y = [ops.nodeDisp(1, 2)]; u3y = [ops.nodeDisp(3, 2)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return None
        sx.append(ops.nodeDisp(2, 1)); sy.append(ops.nodeDisp(2, 2))
        m1x.append(ops.nodeDisp(1, 1)); m1y.append(ops.nodeDisp(1, 2)); u3y.append(ops.nodeDisp(3, 2))
    return (np.array(sx), np.array(sy), np.array(m1x), np.array(m1y), np.array(u3y))


def test_TC3_rigidlink_master_ux_prescribed_uy_free():
    proj = _run_rigidlink_master((PROJ,), ("Diagonal",))
    assert proj is not None, "projection rigidLink-master run failed"
    sx, sy, m1x, m1y, u3y = proj
    assert _rms(m1x) > 1e-3, "master ux (prescribed) not moving"
    assert _rms(m1y) > 1e-3, "master uy (free, loaded) not moving"
    assert _rms(u3y) > 1e-4, "downstream uy not driven"
    # slave ux is a DERIVED-prescribed DOF -> tracks the prescribed master ux exactly
    assert np.max(np.abs(sx - m1x)) < 1e-10, (
        f"slave ux does not track prescribed master ux (max {np.max(np.abs(sx-m1x)):.2e})")
    # slave uy is a NORMAL projected DOF tied to the free master uy -> exact equalDOF-style tie
    # (both are co-integrated free DOFs; projection holds it to machine precision)
    assert np.max(np.abs(sy - m1y)) < 1e-10, (
        f"projected uy tie not held alongside the derived ux (max {np.max(np.abs(sy-m1y)):.2e})")


# --------------------------------------------------------------------------
# TC4 / TC5 — refusals
# --------------------------------------------------------------------------
def _runs_ok(build):
    try:
        build()
        return ops.analyze(1, 1e-3) == 0
    except Exception:
        return False


def test_TC4_mixed_free_and_prescribed_master_refused():
    """A single slave DOF tied to BOTH a free master and a prescribed master (equationConstraint)
    is refused — the prescribed displacement tie cannot be held by the free-DOF projection."""
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)     # prescribed master
        ops.node(2, 0.0); ops.mass(2, _M)     # free master
        ops.node(3, 0.0); ops.mass(3, _M)     # slave tied to BOTH
        # u3 = 0.5 u1(prescribed) + 0.5 u2(free)  ==  1*u3 - 0.5*u1 - 0.5*u2 = 0
        ops.equationConstraint(3, 1, 1.0, 1, 1, -0.5, 2, 1, -0.5)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(1, 1, 0.3, "-const")            # node 1 prescribed -> mixed master on slave 3
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert not _runs_ok(build), "mixed free+prescribed master should be refused"


def test_TC5_prescribed_slave_still_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)
        ops.node(2, 1.0); ops.mass(2, _M)
        ops.equalDOF(1, 2, 1)                 # node 2 = slave
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(2, 1, 0.3, "-const")           # prescribe the SLAVE -> overconstraint
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert not _runs_ok(build), "prescribed-on-slave should still be refused"


def test_TC6_zero_free_dof_refused_not_segfault():
    """A model whose only DOFs are a prescribed master + its derived slave has NO free equation
    to solve. Excluding the derived slave would leave a 0-size SOE -> the handler refuses with a
    named error (ADR §5.2: zero free DOFs -> clean error, not segfault)."""
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, _M)     # prescribed master
        ops.node(2, 0.0); ops.mass(2, _M)     # equalDOF slave -> derived (excluded) ; no free DOF left
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.sp(1, 1, 0.3, "-const")
        ops.constraints(PROJ); ops.numberer("Plain"); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    assert not _runs_ok(build), "zero-free-DOF model should be refused cleanly (no segfault)"
