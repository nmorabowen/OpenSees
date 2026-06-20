"""ADR-30 Phase 3 — tie-force query, EQ_Constraint, energy closure.

P3 adds the "free byproduct" tie-force query f = M(a_raw - a_proj) (the constraint
force the projection applies), exposed via `ladrunoProjectionTieForce nodeTag dof`;
extends the handler to EQ_Constraint groups (a single constrained DOF tied to a
coefficient vector of retained DOFs = the multi-master general-C row); and validates
energy closure (the projection introduces no spurious dissipation).

  T8     tie-force query returns the exact internal constraint force, equal-and-
         opposite across the tie (momentum-conservation corollary): for a two-mass
         equalDOF under constant force F, f_master = +F*m2/(m1+m2), f_slave = -that.
  T9     energy closure: a tied undamped spring-mass free-vibration conserves total
         mechanical energy (no secular drift / dissipation from the projection).
  EQ     an `equationConstraint` (u_c = sum c_k u_{r_k}) is enforced by the projection
         exactly like an MP tie — proving the EQ groups flow through.
  GUARD  the query errors cleanly when the active handler is not LadrunoProjection.
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


# --------------------------------------------------------------------------
# T8 — tie-force query
# --------------------------------------------------------------------------
def test_T8_tie_force_query_is_the_internal_constraint_force():
    m1, m2, F = 2.0, 3.0, 10.0
    Mtot = m1 + m2
    f_exact = F * m2 / Mtot          # internal tie force on the master (= +6.0 here)

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equalDOF(1, 2, 1)            # node1 retained (master), node2 constrained
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, F)
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")
    assert ops.analyze(1, 1e-3) == 0

    f_master = ops.ladrunoProjectionTieForce(1, 1)
    f_slave = ops.ladrunoProjectionTieForce(2, 1)
    assert abs(f_master - f_exact) < 1e-9, f"tie f_master {f_master} != {f_exact}"
    assert abs(f_slave + f_exact) < 1e-9, f"tie f_slave {f_slave} != {-f_exact}"
    # equal-and-opposite (momentum-conservation corollary: sum of tie forces = 0)
    assert abs(f_master + f_slave) < 1e-12


# --------------------------------------------------------------------------
# T9 — energy closure (no spurious dissipation)
# --------------------------------------------------------------------------
def test_T9_energy_closure_no_spurious_dissipation():
    m1, m2, k, d0 = 2.0, 3.0, 100.0, 1.0
    Mtot = m1 + m2
    omega = math.sqrt(k / Mtot)
    dt = 0.02 * (2.0 / omega)        # small dt -> tight CD energy ripple
    nsteps = 2000

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.equalDOF(1, 2, 1)
    ops.setNodeDisp(1, 1, d0, "-commit")     # compliant ICs (both tied DOFs at d0)
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")

    E0 = 0.5 * k * d0 * d0           # all PE at rest
    energies = []
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        u = ops.nodeDisp(1, 1)
        v1 = ops.nodeVel(1, 1); v2 = ops.nodeVel(2, 1)
        ke = 0.5 * m1 * v1 * v1 + 0.5 * m2 * v2 * v2
        pe = 0.5 * k * u * u
        energies.append(ke + pe)

    energies = np.array(energies)
    assert len(energies) == nsteps, "run terminated early (instability?)"
    drift = np.max(np.abs(energies - E0)) / E0
    assert drift < 1e-3, f"energy drift {drift:.3e} exceeds 1e-3 of peak — spurious dissipation/growth"


# --------------------------------------------------------------------------
# EQ — EQ_Constraint enforced by the projection
# --------------------------------------------------------------------------
def test_EQ_equationConstraint_enforced_by_projection():
    m1, m2, F = 2.0, 3.0, 10.0
    a_exact = F / (m1 + m2)

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)        # retained (master)
    ops.node(2, 1.0); ops.mass(2, m2)        # constrained (slave)
    # equationConstraint: cNode cDOF cCoef  rNode rDOF rCoef  -> u_c = -(rCoef/cCoef) u_r
    # tie u_2 = u_1 :  1*u_c + (-1)*u_r = 0
    ops.equationConstraint(2, 1, 1.0, 1, 1, -1.0)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, F)
    ops.constraints(PROJ)
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")

    tie_err = 0.0
    for _ in range(200):
        assert ops.analyze(1, 1e-3) == 0, "EQ-constrained analyze failed"
        tie_err = max(tie_err, abs(ops.nodeDisp(1, 1) - ops.nodeDisp(2, 1)))
    assert tie_err < 1e-12, f"EQ tie u_c=u_r not held ({tie_err:.3e})"
    assert abs(ops.nodeAccel(1, 1) - a_exact) < 1e-9, "combined-mass accel wrong under EQ"
    # the tie-force query also works on an EQ group
    f = ops.ladrunoProjectionTieForce(2, 1)
    assert abs(f + F * m2 / (m1 + m2)) < 1e-9


# --------------------------------------------------------------------------
# GUARD — query errors when the active handler is not LadrunoProjection
# --------------------------------------------------------------------------
def test_tie_force_query_requires_projection_handler():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, 1.0)
    ops.node(2, 1.0); ops.mass(2, 1.0)
    ops.equalDOF(1, 2, 1)
    ops.constraints("Transformation")        # NOT LadrunoProjection
    ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-10, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")
    ops.analyze(1, 1e-3)
    # the query must refuse (raise or signal) rather than silently return a value
    raised = False
    try:
        ops.ladrunoProjectionTieForce(1, 1)
    except Exception:
        raised = True
    assert raised, "tie-force query should error when the handler is not LadrunoProjection"


def test_constraint_to_missing_node_refused():
    """Gate-D: a constraint that references a node not in this domain (the
    partition-crossing case under OpenSeesSP/MP, v1 partition-interior only) must be
    REFUSED, not silently mis-assembled (treated as SP-fixed)."""
    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0); ops.mass(2, 1.0)
        ops.equalDOF(1, 2, 1)
        # a second tie to a node that does NOT exist (tag 99) -> partition-crossing proxy
        ops.equalDOF(1, 99, 1)
        ops.numberer("Plain"); ops.constraints(PROJ); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")
    # either equalDOF rejects the missing node, or the handler guard refuses at analyze —
    # what matters is it does NOT silently run.
    ran = False
    try:
        build()
        ran = (ops.analyze(1, 1e-3) == 0)
    except Exception:
        ran = False
    assert not ran, "a constraint to a missing/cross-partition node must be refused"
