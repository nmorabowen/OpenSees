"""ADR-30 Phase 0 — falsify & baseline for LadrunoProjectionHandler.

NO SRC change. These tests exercise only already-shipped pieces
(CentralDifferenceLadruno 33003 + the upstream Transformation handler +
Diagonal/Full/Band SOEs) to nail down the premises the ADR-30 design rests on,
BEFORE a line of handler code is written. Hardened after the Gate-0 adversarial
panel (run wf_853e5055): the core T6 falsification was confirmed SOUND (matched
to <1% by an OpenSees-free analytic eigen-solution, now baked in as `_analytic_uy`);
the massless test was rewritten because its original "non-finite (NaN/Inf)" premise
was empirically false.

  T6  (the load-bearing falsification, ADR §2/D2): a constraint whose condensed
      mass `T^T M T` has OFF-DIAGONAL coupling (an offset slave's translational
      mass producing a translation<->rotation coupling on the master) is integrated
      WRONG by `system Diagonal`, which silently keeps only diag(A) and drops those
      off-diagonals. We triangulate against TWO independent references: (a) an
      OpenSees-free closed-form modal solution of the condensed 3x3 system
      (`_analytic_uy`, handler-independent), and (b) implicit Newmark+Full. The
      suspect CD+Transformation+Diagonal run is shown to disagree grossly with both
      (its coupling-induced uy is identically zero), proving Diagonal is *wrong*, not
      merely different — exactly why ADR-30 needs a projection handler.
      NB: ADR §6 names rigidDiaphragm; a 2D `rigidLink -beam` with an offset point
      mass is the minimal faithful stand-in (a diaphragm slave carries the same
      transport coupling; a diaphragm is N such slaves). See ADR §6 note.

  MASSLESS (ADR §2.5): a free DOF that keeps its own equation but has zero lumped
      mass (the situation the projection handler's Plain-style assembly creates for a
      massless *slave*) makes the assembled M singular. The SOE layer cannot be
      trusted to police this: `system Diagonal` REFUSES with an opaque solve-time
      pivot error (analyze -> nonzero), while `FullGeneral`/`BandGeneral` swallow the
      singular LAPACK factorization and return SUCCESS (rc=0) with a non-physical
      result. Neither is acceptable -> the handler must scan for massless DOFs at
      handle() time and fail early with an actionable, recipe-bearing message.

  BASELINE (ADR §7): the dense Full SOE is the correct-but-costly path; the Diagonal
      recipe is O(n). Reported informationally — the dispositive argument is
      asymptotic memory (dense O(n^2), infeasible at the ~50k DOF the ADR targets),
      not a microbenchmark ratio.

Theory: Ladruno_implementation/30_ladruno_explicit_constraint_projection_adr.md.
"""
import math
import time

import numpy as np
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CDL = "CentralDifferenceLadruno"


# --------------------------------------------------------------------------
# T6 model: master node + offset slave tied by a rigid beam link.
#
# 2D, ndf=3. Master M(0,0) grounded by zeroLength springs in ux,uy,rz.
# Slave S(d,0) rigidly linked to M (rigidLink -beam) carries a translational
# point mass m. Under the rigid offset r=(d,0), S's y-displacement is
# uy_S = uy_M + d * rz_M (verified against RigidBeam.cpp transform), so the
# condensed master mass picks up a (uy,rz) coupling term m*d and a rotational
# inertia m*d^2 -- both ABSENT from the bare nodal mass. Master also carries its
# own diag mass (m0,m0,Irz) so the assembled 3x3 is non-singular (the pure-offset
# block alone is rank-1).
#
#   M_master = [[m0+m,   0,     0    ],
#               [ 0,   m0+m,   m*d   ],
#               [ 0,   m*d,  Irz+m*d^2]]      K = diag(kx,ky,kr)
#
# We kick rz with an initial angular velocity and watch the INDUCED uy. The
# grounding springs add NO uy-rz stiffness coupling, so uy is driven ONLY by the
# mass off-diagonal m*d. With the full mass the rz motion drags uy through it;
# with a diagonal mass that coupling vanishes and uy stays IDENTICALLY 0. uy is
# therefore a clean, isolated discriminator.
# --------------------------------------------------------------------------
_M0 = 1.0
_IRZ = 1.0
_MS = 1.0
_D = 2.0
_KX = 100.0
_KY = 100.0
_KR = 100.0
_RZ_V0 = 1.0     # initial angular velocity on the master rz DOF


def _build_offset_link_model():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)

    # ground anchor (fully fixed) for the springs
    ops.node(10, 0.0, 0.0)
    ops.fix(10, 1, 1, 1)

    # master
    ops.node(1, 0.0, 0.0)
    ops.mass(1, _M0, _M0, _IRZ)

    # slave, offset by d in x, carrying a translational point mass only
    ops.node(2, _D, 0.0)
    ops.mass(2, _MS, _MS, 0.0)

    # grounding springs on the master (ux,uy,rz) — independent, no cross-coupling
    ops.uniaxialMaterial("Elastic", 1, _KX)
    ops.uniaxialMaterial("Elastic", 2, _KY)
    ops.uniaxialMaterial("Elastic", 3, _KR)
    ops.element("zeroLength", 1, 10, 1, "-mat", 1, 2, 3, "-dir", 1, 2, 3)

    # rigid beam link: master retained (1), slave constrained (2) -> S's offset
    # translational mass condenses onto M with transport coupling.
    ops.rigidLink("beam", 1, 2)


def _run_offset_link(integrator_args, system_args, dt, nsteps,
                     algo="Linear", test_args=("NormDispIncr", 1e-12, 10)):
    _build_offset_link_model()
    ops.setNodeVel(1, 3, _RZ_V0, "-commit")   # kick the rotational DOF
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(*system_args)
    ops.test(*test_args)
    ops.algorithm(algo)
    ops.integrator(*integrator_args)
    ops.analysis("Transient")
    uy = [ops.nodeDisp(1, 2)]
    rz = [ops.nodeDisp(1, 3)]
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        uy.append(ops.nodeDisp(1, 2))
        rz.append(ops.nodeDisp(1, 3))
    return np.array(uy), np.array(rz)


def _analytic_uy(t):
    """Handler-independent reference: closed-form modal response of the condensed
    (uy,rz) block. M and K are assembled BY HAND from the known physics (not via
    the Transformation handler), so this anchors T6 against a construction with no
    OpenSees constraint code in the loop. IC = master rz velocity kick, uy at rest
    (the reduced-coordinate state; the slave is fully determined by the master)."""
    t = np.asarray(t, dtype=float)
    M22 = np.array([[_M0 + _MS, _MS * _D],
                    [_MS * _D, _IRZ + _MS * _D * _D]])
    K22 = np.array([[_KY, 0.0], [0.0, _KR]])
    w2, phi = np.linalg.eig(np.linalg.solve(M22, K22))   # inv(M)K
    w = np.sqrt(w2.real)
    phi = phi.real
    for i in range(phi.shape[1]):                         # M-orthonormalize modes
        phi[:, i] /= math.sqrt(phi[:, i] @ M22 @ phi[:, i])
    qdot0 = phi.T @ M22 @ np.array([0.0, _RZ_V0])          # modal velocities (u0=0)
    uy = np.zeros_like(t)
    for i in range(len(w)):
        uy += phi[0, i] * (qdot0[i] / w[i]) * np.sin(w[i] * t)   # uy = comp 0
    return uy


def _rms(a):
    return float(np.sqrt(np.mean(a * a)))


def _rel_rms_diff(a, b):
    n = min(len(a), len(b))
    a, b = a[:n], b[:n]
    denom = _rms(a) + _rms(b)
    if denom == 0.0:
        return 0.0
    return _rms(a - b) / (0.5 * denom)


def _all_finite(a):
    return bool(np.all(np.isfinite(a)))


# --------------------------------------------------------------------------
# T6 — the falsification
# --------------------------------------------------------------------------
def test_T6_diagonal_soe_drops_condensed_offdiagonal_mass():
    """ADR-30 D2/T6: Transformation + `system Diagonal` is WRONG when the
    condensed mass has off-diagonal coupling; Full SOE is right. Triangulated
    against an OpenSees-free analytic modal solution AND an implicit reference so
    the verdict is 'Diagonal is wrong', not merely 'two runs differ'."""
    dt = 0.01
    nsteps = 400
    t = np.arange(nsteps + 1) * dt

    uy_diag, rz_diag = _run_offset_link((CDL,), ("Diagonal",), dt, nsteps)
    uy_full, rz_full = _run_offset_link((CDL,), ("FullGeneral",), dt, nsteps)
    uy_imp, rz_imp = _run_offset_link(
        ("Newmark", 0.5, 0.25), ("FullGeneral",), dt, nsteps,
        algo="Newton", test_args=("NormDispIncr", 1e-10, 50))
    uy_anchor = _analytic_uy(t)

    assert _all_finite(uy_diag) and _all_finite(uy_full) and _all_finite(uy_imp)

    # (1) HANDLER-INDEPENDENT anchor: the correct explicit-Full run reproduces the
    #     hand-assembled analytic coupled-mode uy. This is the leg that proves the
    #     reference is physically right, not just self-consistent. (panel: 0.37%)
    anchor_err = _rel_rms_diff(uy_full, uy_anchor)
    assert anchor_err < 2e-2, (
        f"explicit-Full does not match the analytic condensed-mass reference "
        f"(rel-rms {anchor_err:.3e}) — the off-diagonal physics is not what we think")

    # (2) implicit-Full corroborates explicit-Full (scheme independence).
    assert _rel_rms_diff(uy_full, uy_imp) < 5e-2, (
        f"explicit-Full and implicit-Full disagree on uy "
        f"(rel-rms {_rel_rms_diff(uy_full, uy_imp):.3e})")

    # (3) the correct response has REAL coupled uy motion (the physics exists).
    coupled_uy = _rms(uy_full)
    assert coupled_uy > 1e-4, (
        f"reference uy coupling too weak to discriminate ({coupled_uy:.3e})")

    # (4) THE FALSIFICATION: Diagonal drops the coupling -> its coupling-induced uy
    #     is identically dead, so it disagrees grossly with the correct answer.
    assert _rms(uy_diag) < 1e-9, (
        "Diagonal uy not suppressed to ~0 as expected — re-examine the model "
        f"(diag uy-rms {_rms(uy_diag):.3e} vs correct {coupled_uy:.3e})")
    assert _rel_rms_diff(uy_diag, uy_full) > 0.5, (
        "Diagonal SOE unexpectedly tracked the correct coupled response — T6 "
        f"premise not reproduced (rel-rms {_rel_rms_diff(uy_diag, uy_full):.3e})")

    # (5) Diagonal does not merely lose uy — it DETUNES the coupled rz mode too
    #     (dropping m*d corrupts the whole 2x2 block; rz period 1.54s -> 1.40s).
    #     Documented so a future reader does not think 'only uy is affected'.
    assert _rel_rms_diff(rz_diag, rz_full) > 0.1, (
        "Diagonal rz unexpectedly matched the correct rz — the coupling drop "
        f"should detune rz (rel-rms {_rel_rms_diff(rz_diag, rz_full):.3e})")


# --------------------------------------------------------------------------
# MASSLESS — ADR §2.5: the SOE layer cannot police a zero-mass free DOF
# --------------------------------------------------------------------------
def _run_massless(system_args, mass):
    """A free DOF (keeps its own equation, like a projection-handler slave would)
    with `mass` lumped on it; one CD starter step under `system_args`.
    Returns (analyze_rc, u, a)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    ops.node(2, 0.0, 0.0)
    ops.fix(2, 0, 1)           # free in x — keeps its own equation
    ops.mass(2, mass, 0.0)
    ops.uniaxialMaterial("Elastic", 1, 100.0)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    rc = ops.analyze(1, 0.001)
    return rc, ops.nodeDisp(2, 1), ops.nodeAccel(2, 1)


def test_massless_dof_is_not_policeable_by_the_soe_layer():
    """ADR §2.5: zero lumped mass on a free DOF makes the assembled M singular.
    The SOE layer gives no actionable diagnosis — which is why the handler must
    detect massless DOFs at handle() time, not lean on the SOE.

    HISTORY (premise drift, deliberate): when this test was written, Full/Band
    SWALLOWED a first-pivot-singular LAPACK factor (`return -info+1` defect) and
    reported rc == 0 with a non-physical result — the silently-wrong outcome
    that first motivated the handle()-time scan. The ADR-76 audit fixed that
    (all three LAPACK solvers now refuse; see the LEDGER_quirks first-pivot-
    singularity entry), so the SOE layer is no longer silently WRONG — but it
    is still opaque (a solve-time abort that names no DOF), so the handler
    contract stands unchanged.

      * Diagonal  -> opaque solve-time pivot abort (analyze rc != 0).
      * Full/Band -> singular factor now correctly REFUSED (rc != 0) — post
                     ADR-76; pre-fix these returned rc == 0, silently wrong.

    A positive control (valid mass) isolates the zero-mass DOF as the cause."""
    # positive control: a real mass solves cleanly with a finite acceleration.
    rc_ctrl, u_ctrl, a_ctrl = _run_massless(("Diagonal",), 1.0)
    assert rc_ctrl == 0 and math.isfinite(a_ctrl), (
        f"control (mass=1) failed unexpectedly (rc={rc_ctrl}, a={a_ctrl})")

    # Diagonal: clean refusal — the step cannot be solved (rc != 0), state untouched.
    rc_diag, u_diag, a_diag = _run_massless(("Diagonal",), 0.0)
    assert rc_diag != 0, (
        "Diagonal SOE unexpectedly 'solved' a zero-mass DOF — §2.5 not reproduced")
    assert u_diag == 0.0, "rejected starter should leave displacement at its IC"

    # Full / Band: post-ADR-76 the singular factorization is REFUSED. A rc == 0
    # here means the singular-as-success LAPACK defect has REGRESSED.
    rc_full, _, _ = _run_massless(("FullGeneral",), 0.0)
    rc_band, _, _ = _run_massless(("BandGeneral",), 0.0)
    assert rc_full != 0 and rc_band != 0, (
        "Full/Band reported SUCCESS on a singular mass matrix — the ADR-76 "
        "singular-as-success LAPACK defect has regressed "
        f"(rc_full={rc_full}, rc_band={rc_band})")

    # the contract the handler must provide: every SOE refuses, but none can say
    # WHICH DOF is massless — an actionable, named diagnosis still belongs at
    # handle() time.


# --------------------------------------------------------------------------
# BASELINE — per-step cost, Diagonal (O(n)) vs Full (dense) explicit
# --------------------------------------------------------------------------
def _build_mass_chain(n, k=100.0, m=1.0):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(0, 0.0)
    ops.fix(0, 1)
    ops.uniaxialMaterial("Elastic", 1, k)
    for i in range(1, n + 1):
        ops.node(i, float(i))
        ops.mass(i, m)
        ops.element("zeroLength", i, i - 1, i, "-mat", 1, "-dir", 1)


def _time_chain(system_args, n=400, nsteps=200, dt=0.005):
    _build_mass_chain(n)
    ops.setNodeVel(1, 1, 1.0, "-commit")
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", 1e-12, 10)
    ops.algorithm("Linear")
    ops.integrator(CDL)
    ops.analysis("Transient")
    t0 = time.perf_counter()
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
    return time.perf_counter() - t0


def test_baseline_diagonal_vs_full_is_informational():
    """ADR §7 baseline. The dispositive cost argument is asymptotic, not a single
    timing ratio: a dense Full SOE is O(n^2) MEMORY and O(n^3) factor, infeasible at
    the ~50k-DOF models the ADR targets (50k^2 doubles ≈ 20 GB), whereas the
    Diagonal recipe is O(n) in both. We MEASURE and record (Full is already several-
    fold slower even at this small n), and assert only that both complete with a
    finite time — the projection handler's job is to keep the O(n) Diagonal path
    while enforcing constraints correctly."""
    t_diag = _time_chain(("Diagonal",))
    t_full = _time_chain(("FullGeneral",))
    assert math.isfinite(t_diag) and math.isfinite(t_full)
    print(f"\n[ADR30-P0 baseline] n=400 explicit per-200-steps: "
          f"Diagonal={t_diag*1e3:.1f}ms Full={t_full*1e3:.1f}ms "
          f"(ratio {t_full/max(t_diag,1e-9):.1f}x; environment-noisy at this size). "
          f"Real argument: O(n) vs O(n^2) memory at 50k DOF "
          f"(dense ~20 GB, infeasible).")
