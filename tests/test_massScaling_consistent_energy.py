"""V4 — CONSISTENT (Olovsson) mass-scaling energy closure in the EnergyBalanceRecorder.
ADR 38 follow-up. classTags 33008 / 33010 / 33012.

Consistent mass scaling injects a CENTROIDAL element scaling mass
    M_bar_e = beta_e [diag(m_a) - m m^T / M_e]
that lives ONLY inside the integrator (a matrix-free PCG operand); it is NEVER written
to Node/Element mass. The EnergyBalanceRecorder builds its kinetic energy from
Node::getMass()/Element::getMass(), so before V4 it saw only the lumped diagonal M_lump
and the recorded KE was 1/2 v^T M_lump v — short of the true run kinetic energy
1/2 v^T M_tilde v = 1/2 v^T (M_lump + sum_e M_bar_e) v. The consistent-path energy
balance therefore did NOT close (the lumped sibling's does — it injects into Node mass).

V4 routes the per-element node-major M_bar blocks through a process-global registry
(Ladruno::MassScalingEnergyRegistry); the shared energy kernel adds the missing
1/2 v^T M_bar_e v per element. Two genuine change-detectors (each FAILS on the pre-V4
binary, where the recorder KE equals exactly the lumped KE):

  TV4-ORACLE  (all 3 consistent integrators): the recorder's KE matches an INDEPENDENT
              analytic 1/2 v^T M_tilde v to ~machine precision, every step. M_bar is
              built from the closed-form sizing (dt_e = L/c, beta = (dtTarget/dt_e)^2 - 1
              with betaK = 0), assembled per element as 1/2 sum_e beta_e (m_a/2)(v_a-v_b)^2.
              Instantaneous (no time integration), so dissipation/stability play no role
              and it holds identically for the Noh-Bathe families. Non-vacuity: M_bar is a
              large MAJORITY of the KE here, so a pre-V4 (M_bar-blind) recorder is off by
              >50%.
  TV4-CLOSE   (CentralDifference): undamped free vibration conserves mechanical energy
              KE+IE. With a velocity IC the conserved level is E0 = 1/2 v^T M_tilde v; the
              recorded KE+IE stays flat (small drift) ONLY because KE now carries the
              cross-node M_bar energy. A pre-V4 KE = 1/2 v^T M_lump v gives
              KE_lump+IE = E0 - 1/2 v^T M_bar v, which OSCILLATES with the (large) M_bar
              term — reconstructed here and asserted to drift far more than the V4 record.

Model: the oracle's Case-A axial bar (5 bulk L=1, one tiny L=0.05, 5 bulk L=1; Truss
-rho), seeded with an ALTERNATING (deformation-rich) velocity field so M_bar — which by
construction loads only the relative/deformation modes — carries a large share of the KE.
No load, no damping => DW = ULW = 0. dtTarget=0.012 scales EVERY element.

Theory: Ladruno_implementation/38_ladruno_consistent_mass_scaling_adr.md and the
mass-scaling developer handoff (V4 item). Oracle: mass_scaling_consistent/.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

CONS = "CentralDifferenceSMSConsistent"
EBC  = "ExplicitBatheSMSConsistent"
LNVD = "ExplicitBatheLNVDSMSConsistent"
CDL  = "CentralDifferenceLadruno"      # non-publishing base integrator (teardown control)

# Case-A bar: a tiny interior element (~20x shorter -> ~400x smaller dt_e) between bulk.
LENGTHS = [1.0] * 5 + [0.05] + [1.0] * 5
E, RHO, A = 1.0e4, 1.0, 1.0
C_WAVE = math.sqrt(E / RHO)          # axial wave speed -> truss element dt_e = L / c
DT_TARGET = 0.012                    # above even the bulk dt_e (~0.01) -> all 11 scaled


def _build_bar():
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.node(1, 0.0, 0.0)
    ops.fix(1, 1, 1)
    coords = {1: 0.0}
    x = 0.0
    for i, L in enumerate(LENGTHS):
        x += L
        nd = i + 2
        ops.node(nd, x, 0.0)
        ops.fix(nd, 0, 1)            # axial only (y fixed)
        coords[nd] = x
    ops.uniaxialMaterial("Elastic", 1, E)
    for i in range(len(LENGTHS)):
        ops.element("Truss", i + 1, i + 1, i + 2, A, 1, "-rho", RHO)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 1)
    ops.algorithm("Linear")
    return coords


def _ke_lumped(vmap):
    """1/2 v^T M_lump v exactly as the recorder's getMass() leg forms it: Truss::getMass
    lumps m = 0.5*rho*L to EACH node's translational DOF (A-INDEPENDENT — `-rho` is mass per
    unit length, not density*area; see Truss.cpp), summed over elements. This is the KE a
    recorder blind to the consistent M_bar would report."""
    ke = 0.0
    for i, L in enumerate(LENGTHS):
        m_half = RHO * L / 2.0          # Truss lumped node mass = 0.5*rho*L (no A)
        va = vmap.get(i + 1, 0.0)
        vb = vmap.get(i + 2, 0.0)
        ke += 0.5 * m_half * va * va + 0.5 * m_half * vb * vb
    return ke


def _ke_mbar_analytic(vmap):
    """1/2 v^T M_bar v assembled independently from the closed-form Olovsson sizing.
    For a 2-node axial truss M_bar_e = beta_e (m_a/2)[[1,-1],[-1,1]], m_a = rho*A*L/2, so
    1/2 v^T M_bar_e v = 1/2 beta_e (m_a/2) (v_a - v_b)^2. dt_e = L/c, betaK = 0 ->
    beta_e = (dtTarget/dt_e)^2 - 1 for every sub-target element (all are, here)."""
    ke = 0.0
    for i, L in enumerate(LENGTHS):
        dt_e = L / C_WAVE
        if dt_e >= DT_TARGET:
            continue                 # not scaled -> no M_bar
        beta = (DT_TARGET / dt_e) ** 2 - 1.0
        m_a = RHO * L / 2.0             # element lumped node mass = 0.5*rho*L (A-independent)
        dv = vmap.get(i + 1, 0.0) - vmap.get(i + 2, 0.0)
        ke += 0.5 * beta * (m_a / 2.0) * dv * dv
    return ke


def _run(integrator_args, dt, nsteps, efile, alt_vel=1.0):
    """Seed an ALTERNATING velocity field (deformation-rich -> loads M_bar) and free-
    vibrate. Returns (rows, ke_lump, ke_mbar) aligned per committed step:
      rows     : EnergyBalance file rows  [time KE IE DW ULW RES ERR%]
      ke_lump  : analytic 1/2 v^T M_lump v at each step
      ke_mbar  : analytic 1/2 v^T M_bar  v at each step (independent oracle)."""
    coords = _build_bar()
    for nd in coords:
        if coords[nd] > 0.0:
            ops.setNodeVel(nd, 1, alt_vel * ((-1) ** nd), "-commit")
    ops.recorder("EnergyBalance", "-file", efile, "-time")
    ops.integrator(*integrator_args)
    ops.analysis("Transient")

    ke_lump, ke_mbar = [], []
    for _ in range(nsteps):
        if ops.analyze(1, dt) != 0:
            break
        vmap = {nd: ops.nodeVel(nd, 1) for nd in coords}
        if not all(math.isfinite(v) for v in vmap.values()):
            break
        ke_lump.append(_ke_lumped(vmap))
        ke_mbar.append(_ke_mbar_analytic(vmap))
    ops.remove("recorders")

    rows = []
    with open(efile) as fh:
        for line in fh:
            v = line.split()
            if v:
                rows.append([float(x) for x in v])
    return rows, ke_lump, ke_mbar


# --------------------------------------------------------------------------
# TV4-ORACLE: recorder KE == 1/2 v^T M_tilde v = lumped KE + analytic M_bar KE, every
# step, to ~machine precision — for ALL THREE consistent integrators (each must publish
# its blocks; M_bar is integrator-independent so one analytic oracle serves all three).
# Pre-V4 the recorder reports only the lumped KE (M_bar term == 0) => fails hard.
# --------------------------------------------------------------------------
@pytest.mark.parametrize("integ", [
    pytest.param((CONS, DT_TARGET), id="CentralDifference"),
    pytest.param((EBC, 0.5, DT_TARGET), id="ExplicitBathe"),
    pytest.param((LNVD, 0.5, 0.0, DT_TARGET), id="ExplicitBatheLNVD"),
])
def test_consistent_kernel_matches_analytic_mbar(integ, tmp_path):
    rows, ke_lump, ke_mbar = _run(integ, 0.002, 8, str(tmp_path / "e.txt"))
    n = min(len(rows), len(ke_lump))
    assert n >= 5, "too few aligned steps (%d)" % n

    for k in range(n):
        ke_rec = rows[k][1]                 # recorder KE (with -time => col 1)
        ke_pred = ke_lump[k] + ke_mbar[k]   # independent analytic 1/2 v^T M_tilde v
        assert ke_pred > 0.0
        rel = abs(ke_rec - ke_pred) / ke_pred
        assert rel < 1.0e-2, (
            "step %d: recorder KE %.6e != analytic lumped+M_bar %.6e (rel %.2e). The kernel "
            "is not adding exactly 1/2 v^T M_bar v." % (k, ke_rec, ke_pred, rel)
        )

    # non-vacuity: M_bar is a large MAJORITY of the KE on this deformation-rich IC, so a
    # pre-V4 (M_bar-blind) recorder would be off by >50% — this is not a marginal effect.
    min_frac = min(ke_mbar[k] / rows[k][1] for k in range(n))
    assert min_frac > 0.30, (
        "M_bar KE is only %.0f%% of the recorded KE (expected a large majority); the IC is "
        "not loading the scaling mass enough to be a strong change-detector" % (100 * min_frac)
    )


# --------------------------------------------------------------------------
# TV4-CLOSE: with the M_bar KE included, the central-difference (non-dissipative) run
# CONSERVES mechanical energy KE+IE under consistent scaling — the balance closes. A
# velocity IC sets the conserved level to E0 = 1/2 v^T M_tilde v; KE+IE stays flat (small
# drift) only because KE carries the cross-node M_bar energy. The pre-V4 reconstruction
# KE_lump+IE = E0 - 1/2 v^T M_bar v oscillates with the (large) M_bar term and drifts far
# more — proving the V4 term is what closes the balance.
# --------------------------------------------------------------------------
def test_consistent_energy_conserves_central_difference(tmp_path):
    rows, ke_lump, ke_mbar = _run((CONS, DT_TARGET), 0.002, 800, str(tmp_path / "cd.txt"))
    n = min(len(rows), len(ke_lump))
    assert n > 200, "too few energy rows (%d)" % n

    ke  = [rows[k][1] for k in range(n)]
    ie  = [rows[k][2] for k in range(n)]
    dw  = [rows[k][3] for k in range(n)]
    ulw = [rows[k][4] for k in range(n)]
    ke_peak = max(ke)

    # sanity: undamped, unloaded free vibration -> no damping work, no external work
    assert max(abs(v) for v in dw) < 1.0e-6 * ke_peak, "unexpected damping work (DW != 0)"
    assert max(abs(v) for v in ulw) < 1.0e-6 * ke_peak, "unexpected external work (ULW != 0)"

    # post-V4 conservation: mechanical energy KE+IE flat (drift small)
    mech = [ke[k] + ie[k] for k in range(n)]
    mean_mech = sum(mech) / n
    assert mean_mech > 0.0
    drift_v4 = (max(mech) - min(mech)) / mean_mech
    assert drift_v4 < 0.05, (
        "consistent-SMS mechanical energy not conserved: drift %.1f%% (expect < 5%%). "
        "With the M_bar KE included KE+IE should be the flat level E0." % (100 * drift_v4)
    )

    # non-vacuity / pre-V4 reconstruction: KE_lump + IE (what a M_bar-blind recorder would
    # report) drifts far more, because it equals E0 - 1/2 v^T M_bar v and M_bar is large.
    mech_pre = [ke_lump[k] + ie[k] for k in range(n)]
    mean_pre = sum(mech_pre) / n
    drift_pre = (max(mech_pre) - min(mech_pre)) / abs(mean_pre) if mean_pre else float("inf")
    assert drift_pre > 4.0 * drift_v4 and drift_pre > 0.20, (
        "test weak: the pre-V4 (lumped-KE) reconstruction drifts only %.1f%% vs the V4 record "
        "%.1f%% — the M_bar term is not exercising the closure" % (100 * drift_pre, 100 * drift_v4)
    )

    # and the M_bar KE really is a large share (ties the two reconstructions together)
    peak_mbar_frac = max(ke_mbar[k] / ke[k] for k in range(n))
    assert peak_mbar_frac > 0.30, (
        "peak M_bar KE only %.0f%% of KE — closure gate not exercising the scaling mass"
        % (100 * peak_mbar_frac)
    )


# --------------------------------------------------------------------------
# TV4-TEARDOWN: the registry is process-global mutable state, so the owner-guarded
# clear() in the integrator destructor must actually retire the blocks — otherwise a
# later, NON-publishing run on the same element tags would have its KE inflated by stale
# M_bar (double counting). Leg 1: consistent SMS records KE above the lumped-only KE
# (M_bar present). Leg 2: after _build_bar's wipe() destroys that integrator, a plain
# (non-publishing) CentralDifferenceLadruno on the SAME tags must record EXACTLY the
# analytic lumped KE — proving the registry was cleared and no stale M_bar leaks in.
# --------------------------------------------------------------------------
def test_registry_cleared_after_teardown(tmp_path):
    # leg 1 — consistent publishes: recorded KE exceeds the lumped-only KE
    rows_c, ke_lump_c, _ = _run((CONS, DT_TARGET), 0.002, 8, str(tmp_path / "cons.txt"))
    nc = min(len(rows_c), len(ke_lump_c))
    assert nc >= 5
    peak_c = max(rows_c[k][1] for k in range(nc))
    assert max(rows_c[k][1] - ke_lump_c[k] for k in range(nc)) > 0.10 * peak_c, (
        "leg 1: consistent SMS should record KE above the lumped-only KE (M_bar present)"
    )

    # leg 2 — non-publishing integrator on the same tags after teardown: KE == lumped.
    # Plain CentralDifferenceLadruno is stable at a step below the tiny element's dt_cr
    # (~5e-4). _build_bar() inside _run calls wipe(), which destroys the leg-1 integrator
    # (its dtor clear()s the registry).
    rows_p, ke_lump_p, _ = _run((CDL,), 2.0e-4, 8, str(tmp_path / "plain.txt"))
    npl = min(len(rows_p), len(ke_lump_p))
    assert npl >= 5
    peak_p = max(rows_p[k][1] for k in range(npl)) or 1.0
    # A stale M_bar would inflate KE by ~its share of the deformation energy (~77% on this
    # IC); the gate below (1e-3) is >2 orders above the EnergyBalance file's 6-significant-
    # figure rounding (~1e-6) yet ~3 orders BELOW that stale signal, so it cleanly catches a
    # registry that was not cleared without tripping on output precision.
    for k in range(npl):
        rel = abs(rows_p[k][1] - ke_lump_p[k]) / peak_p
        assert rel < 1.0e-3, (
            "step %d: a NON-publishing integrator recorded KE %.6e != analytic lumped %.6e "
            "(rel %.2e) — stale consistent M_bar survived teardown (registry not cleared / "
            "owner guard wrong)" % (k, rows_p[k][1], ke_lump_p[k], rel)
        )
