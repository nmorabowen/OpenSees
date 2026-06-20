"""ADR-30 Phase 1 — LadrunoProjectionHandler core (equalDOF) Zone-A battery.

The handler enforces MP constraints in explicit dynamics by M-orthogonal projection
of the solved accelerations (no penalty, no elimination): every slave DOF keeps its
own equation + diagonal mass, and CentralDifferenceLadruno projects a0 (starter) and
a_{n+1} (update) onto the constraint manifold.

  T1  two-mass equalDOF tie under a constant force: combined-mass trajectory, the tie
      holds to machine eps, and linear momentum tracks the impulse F*t to < 1e-12
      (the projection conserves the retained-conjugate momentum exactly).
  T2  CD+LadrunoProjection vs CD+Transformation (the eliminated-system reference, exact
      for equalDOF since C=I gives a diagonal T^T M T): SHM trajectories match < 1e-9.
  T4  dt sweep at 0.99 / 1.02 * dt_cr WITH and WITHOUT the tie: the tie does not move
      the stability boundary (the projection is Dt-neutral, Gate-A C4).
  T5  conflict battery: double-constraint, chain, SP-on-slave, massless slave,
      non-Diagonal SOE, and a non-projection-aware integrator each FAIL loudly
      (analyze != 0) rather than silently producing a wrong answer.

Built against ADR 30 + Ladruno_implementation/_adr30_p1_design.md (Gate-A R1-R10).
"""
import math

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"
PROJ = "LadrunoProjection"


# --------------------------------------------------------------------------
# T1 — momentum-conserving combined-mass tie
# --------------------------------------------------------------------------
def test_T1_equaldof_tie_combined_mass_and_momentum():
    m1, m2, F = 2.0, 3.0, 10.0
    Mtot = m1 + m2
    a_exact = F / Mtot
    dt, nsteps = 1.0e-3, 400

    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.equalDOF(1, 2, 1)                 # node1 retained, node2 constrained
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(1, F)
    ops.constraints(PROJ)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")

    t = 0.0
    tie_err = 0.0
    mom_err = 0.0
    a_err = 0.0
    for _ in range(nsteps):
        assert ops.analyze(1, dt) == 0, "analyze failed"
        t += dt
        u1 = ops.nodeDisp(1, 1); u2 = ops.nodeDisp(2, 1)
        v1 = ops.nodeVel(1, 1);  v2 = ops.nodeVel(2, 1)
        a1 = ops.nodeAccel(1, 1)
        tie_err = max(tie_err, abs(u1 - u2), abs(v1 - v2))
        # linear momentum of the pair tracks the impulse F*t
        p = m1 * v1 + m2 * v2
        mom_err = max(mom_err, abs(p - F * t))
        a_err = max(a_err, abs(a1 - a_exact))

    assert tie_err < 1e-12, f"tie not held to machine eps (err {tie_err:.3e})"
    assert a_err < 1e-10, f"combined accel wrong (err {a_err:.3e}, expected {a_exact})"
    # momentum error scaled by the final impulse
    assert mom_err < 1e-12 * (1.0 + F * t), f"momentum not conserved (err {mom_err:.3e})"


# --------------------------------------------------------------------------
# T2 — projection reproduces the eliminated (Transformation) reference
# --------------------------------------------------------------------------
def _run_tied_shm(handler_args, m1, m2, k, d0, dt, nsteps):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, m1)
    ops.node(2, 1.0); ops.mass(2, m2)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.equalDOF(1, 2, 1)
    # compliant ICs: both tied DOFs start at d0 (equalDOF built at u=0 -> delta=0)
    ops.setNodeDisp(1, 1, d0, "-commit")
    ops.setNodeDisp(2, 1, d0, "-commit")
    ops.constraints(*handler_args)
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    u = [ops.nodeDisp(1, 1)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        u.append(ops.nodeDisp(1, 1))
    return np.array(u)


def test_T2_projection_matches_transformation_reference():
    m1, m2, k, d0 = 2.0, 3.0, 100.0, 1.0
    omega = math.sqrt(k / (m1 + m2))
    dt = 0.2 * (2.0 / omega)
    nsteps = 300

    u_proj = _run_tied_shm((PROJ,), m1, m2, k, d0, dt, nsteps)
    u_ref  = _run_tied_shm(("Transformation",), m1, m2, k, d0, dt, nsteps)

    n = min(len(u_proj), len(u_ref))
    assert n == nsteps + 1, "a run terminated early (instability?)"
    diff = np.max(np.abs(u_proj[:n] - u_ref[:n]))
    scale = np.max(np.abs(u_ref[:n]))
    assert diff < 1e-9 * (1.0 + scale), (
        f"projection vs Transformation mismatch {diff:.3e} (scale {scale:.3e})")
    # sanity: it is actually oscillating (period ~ 2pi/omega)
    assert scale > 0.5 * d0


# --------------------------------------------------------------------------
# T4 — Dt-neutrality: the tie does not move the stability boundary
# --------------------------------------------------------------------------
def _peak_abs_single(M, k, d0, dt, nsteps):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0); ops.fix(0, 1)
    ops.node(1, 0.0); ops.mass(1, M)
    ops.uniaxialMaterial("Elastic", 1, k)
    ops.element("zeroLength", 1, 0, 1, "-mat", 1, "-dir", 1)
    ops.setNodeDisp(1, 1, d0, "-commit")
    ops.constraints("Plain")
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")
    peak = abs(d0)
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            return float("inf")
        peak = max(peak, abs(ops.nodeDisp(1, 1)))
        if not math.isfinite(peak) or peak > 1e6 * abs(d0):
            return float("inf")
    return peak


def _peak_abs_tied(m1, m2, k, d0, dt, nsteps):
    u = _run_tied_shm((PROJ,), m1, m2, k, d0, dt, nsteps)
    if len(u) < nsteps + 1:
        return float("inf")
    peak = float(np.max(np.abs(u)))
    if not math.isfinite(peak) or peak > 1e6 * abs(d0):
        return float("inf")
    return peak


def test_T4_tie_is_dt_neutral():
    m1, m2, k, d0 = 2.0, 3.0, 100.0, 1.0
    M = m1 + m2
    omega = math.sqrt(k / M)
    dt_cr = 2.0 / omega
    nsteps = 600

    # untied single combined mass (reference stability boundary)
    a_stable = _peak_abs_single(M, k, d0, 0.99 * dt_cr, nsteps)
    a_unstab = _peak_abs_single(M, k, d0, 1.02 * dt_cr, nsteps)
    # tied pair via the projection handler
    b_stable = _peak_abs_tied(m1, m2, k, d0, 0.99 * dt_cr, nsteps)
    b_unstab = _peak_abs_tied(m1, m2, k, d0, 1.02 * dt_cr, nsteps)

    assert math.isfinite(a_stable) and a_stable < 10.0 * d0, "ref should be stable @0.99"
    assert not math.isfinite(a_unstab), "ref should diverge @1.02"
    # the TIE must reproduce the SAME boundary (Dt-neutral)
    assert math.isfinite(b_stable) and b_stable < 10.0 * d0, (
        f"tied run unstable at 0.99 dt_cr (peak {b_stable}) — tie shifted the boundary")
    assert not math.isfinite(b_unstab), (
        "tied run did NOT diverge at 1.02 dt_cr — tie changed the stability limit")


# --------------------------------------------------------------------------
# T5 — conflict battery: every ill-posed setup FAILS loudly (analyze != 0)
# --------------------------------------------------------------------------
def _runs_ok(build):
    """build() configures a model+analysis; return True iff one step succeeds."""
    try:
        build()
        return ops.analyze(1, 1e-4) == 0
    except Exception:
        return False


def _common(handler=(PROJ,), system=("Diagonal",), integrator=(CDL,)):
    ops.numberer("Plain")
    ops.constraints(*handler)
    ops.system(*system)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(*integrator)
    ops.analysis("Transient")


def test_T5_double_constraint_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        for n in (1, 2, 3):
            ops.node(n, float(n)); ops.mass(n, 1.0)
        ops.equalDOF(1, 3, 1)
        ops.equalDOF(2, 3, 1)        # node3 constrained TWICE
        _common()
    assert not _runs_ok(build), "double-constrained slave should be refused"


def test_T5_chain_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        for n in (1, 2, 3):
            ops.node(n, float(n)); ops.mass(n, 1.0)
        ops.equalDOF(1, 2, 1)        # node2 slave
        ops.equalDOF(2, 3, 1)        # node2 retained -> chain
        _common()
    assert not _runs_ok(build), "MP chain should be refused"


def test_T5_sp_on_slave_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0); ops.mass(2, 1.0)
        ops.equalDOF(1, 2, 1)
        ops.fix(2, 1)                # slave also SP-fixed
        _common()
    assert not _runs_ok(build), "SP-on-slave should be refused"


def test_T5_massless_slave_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0); ops.mass(2, 0.0)    # massless slave -> singular own eqn
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 1.0)
        _common()
    assert not _runs_ok(build), "massless slave should be refused at buildMass"


def test_T5_non_diagonal_soe_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0); ops.mass(2, 1.0)
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 1.0)
        _common(system=("FullGeneral",))
    assert not _runs_ok(build), "non-Diagonal SOE should be refused"


def test_T3_ic_violation_aborts_and_projectICs_recovers():
    """Gate-B blocker: an off-manifold committed IC must ABORT (analyze != 0), not run
    on silently; -projectICs must snap it compliant and run. Exercises the full path
    CDL.domainChanged -> enforceIC -> -1 -> DirectIntegrationAnalysis honors the return."""
    def build(handler):
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 2.0)
        ops.node(2, 1.0); ops.mass(2, 3.0)
        ops.equalDOF(1, 2, 1)                       # built at u=0 -> delta=0
        ops.setNodeDisp(1, 1, 1.0, "-commit")
        ops.setNodeDisp(2, 1, 2.0, "-commit")        # incompatible: tie requires u1==u2
        ops.numberer("Plain"); ops.constraints(*handler); ops.system("Diagonal")
        ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
        ops.integrator(CDL); ops.analysis("Transient")

    build((PROJ,))
    assert ops.analyze(1, 1e-3) != 0, "off-manifold ICs must abort (not run silently)"

    build((PROJ, "-projectICs"))
    assert ops.analyze(1, 1e-3) == 0, "-projectICs should snap compliant and run"
    assert abs(ops.nodeDisp(1, 1) - ops.nodeDisp(2, 1)) < 1e-9, "slave not snapped to master"


def test_zero_group_rebuild_no_use_after_free():
    """Gate-B blocker: after a re-domainChanged that drops to ZERO constraint groups
    (the MP removed mid-run), the integrator must not dereference the freed old
    projector. doneNumberingDOF now re-binds the consumer unconditionally (0 push)."""
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.mass(1, 2.0)
    ops.node(2, 1.0); ops.mass(2, 3.0)
    ops.equalDOF(1, 2, 1)
    ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 10.0)
    ops.numberer("Plain"); ops.constraints(PROJ); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator(CDL); ops.analysis("Transient")
    assert ops.analyze(5, 1e-3) == 0, "tied run failed"      # builds + pushes 1-group projector
    ops.remove("mp", 2)                                       # drop the tie -> 0 groups next rebuild
    rc = ops.analyze(5, 1e-3)                                 # re-domainChanged: 0-group push, no UAF
    assert rc == 0, "zero-group rebuild should run cleanly (no use-after-free)"


def test_T5_non_projection_aware_integrator_refused():
    def build():
        ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.mass(1, 1.0)
        ops.node(2, 1.0); ops.mass(2, 1.0)
        ops.equalDOF(1, 2, 1)
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1); ops.load(1, 1.0)
        _common(integrator=("Newmark", 0.5, 0.25))
    assert not _runs_ok(build), "non-projection-aware integrator should be refused"
