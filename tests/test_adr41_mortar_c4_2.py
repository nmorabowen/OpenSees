"""ADR-41 C4.2 — mortar MESH-TYING: commit-cycle ALM (the epsTie-independent bond).

Penalty mesh-tying (C4.1) leaves a residual relative displacement ‖r‖ = O(1/epsTie): the weld is
only as tight as epsTie is large, and a large epsTie ill-conditions the system. Augmented Lagrange
drives ‖r‖ → an epsTie-INDEPENDENT tolerance at a FINITE epsTie by carrying a per-GLOBAL-slave-node
tie multiplier λ_tie (a 3-vector), updated once per Domain::commit() with NO clamp (equality):

    λ_tie,I  <-  λ_tie,I + epsTie · (r_I / a_I)          (r_I = rGlobal_I / aGlobal_I)

λ_tie is COMMITTED-ONLY (the λ_N pattern); r_I is the ORDER-INDEPENDENT global accumulator (a linear
sum, unlike the friction slip — the C4 shared-node resolution). No new DOFs: the SOE equation count
is constant across augmentations. The held-load `analyze_augmented(query=ops.ladrunoMortarTieResidual)`
proc re-commits at zero load increment to drive ‖r‖_∞ < augTol within a step.

Gates (capstone mesh-tying row, ALM phase):
  - penalty ‖r‖ = O(1/epsTie); ALM drives ‖r‖ → augTol at a FINITE epsTie, epsTie-INDEPENDENT;
  - the patch test passes to ~machine precision under augmentation (better than any penalty);
  - SOE equation count constant across augmentations (no new DOFs).

Shares the split-column fixture with C4.1 (a NON-MATCHING tied interface with a shared edge node).
"""
import os
import sys

import pytest

from _testbed import ops
from test_adr41_mortar_c4_1 import _split_column, _static_setup, E, PTOT  # reuse the fixture

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "Ladruno_scripts"))
from analyze_augmented import analyze_augmented  # noqa: E402

pytestmark = [pytest.mark.zone_a]


def _split_column_alm(epsTie, augTol=1.0e-12, maxAug=60):
    """Build + load the tied split column (C4.1 fixture), then drive the tie residual ‖r‖ → augTol
    by held-load Uzawa augmentation. Returns (penalty_res, alm_res, n_aug, u_tip_alm, sysN)."""
    r = _split_column(epsTie=epsTie, tension=True)
    assert r is not None, "penalty tie step did not converge"
    penalty_res = r["tie_res"]
    sysN_before = r["sysN"]
    status, n_aug, alm_res = analyze_augmented(
        ops, maxAug=maxAug, augTol=augTol, query=ops.ladrunoMortarTieResidual)
    sysN_after = ops.systemSize()
    assert status == 0, f"tie augmentation did not converge (status {status}, ‖r‖ {alm_res:.2e})"
    assert sysN_before == sysN_after, "SOE size changed across tie augmentations (λ_tie is not a DOF)"
    return penalty_res, alm_res, n_aug, ops.nodeDisp(15, 3), sysN_before


def test_c4_2_residual_epsTie_independent():
    """THE HEADLINE: the penalty bond ‖r‖ is O(1/epsTie) (a 100× stiffer epsTie ⇒ ~100× tighter weld
    — accuracy bought only with conditioning), but the no-clamp Uzawa drives ‖r‖ below a FIXED augTol
    at ANY epsTie, and the converged tip is epsTie-INDEPENDENT (penalty alone cannot reach)."""
    augTol = 1.0e-12
    lo = _split_column_alm(epsTie=1.0e6, augTol=augTol)
    hi = _split_column_alm(epsTie=1.0e8, augTol=augTol)

    # penalty alone: ‖r‖ ~ t/epsTie ⇒ 100× stiffer epsTie ⇒ ~100× smaller bond (epsTie-DEPENDENT).
    assert lo[0] > 50.0 * hi[0], (
        f"penalty bond not O(1/epsTie): epsTie=1e6 -> {lo[0]:.2e}, 1e8 -> {hi[0]:.2e}")
    # ALM: the bond is driven below augTol at BOTH penalties (the win penalty can't reach).
    assert lo[1] < augTol, f"epsTie=1e6 ALM bond {lo[1]:.2e} !< {augTol}"
    assert hi[1] < augTol, f"epsTie=1e8 ALM bond {hi[1]:.2e} !< {augTol}"
    # and the converged tip is epsTie-INDEPENDENT (the headline ALM win — IDENTICAL, not just small).
    assert abs(lo[3] - hi[3]) < 1.0e-9, (
        f"converged tip epsTie-dependent: 1e6 -> {lo[3]:.6e}, 1e8 -> {hi[3]:.6e}")


def test_c4_2_patch_test_machine_precision():
    """Under augmentation the split == monolithic to ~machine precision: u(z) = (P/EA)·z with the
    residual bond driven below augTol, the transmitted field matches the exact uniaxial bar far
    tighter than the C4.1 penalty tolerance (the λ_tie carries the load — the penalty term → 0)."""
    _, alm_res, _, u_tip, _ = _split_column_alm(epsTie=1.0e6, augTol=1.0e-13)
    u_tip_exact = 2.0 * PTOT / E
    assert alm_res < 1.0e-12, f"tie bond not driven to machine precision (‖r‖ {alm_res:.2e})"
    # the transmitted field now matches the monolithic bar to ~1e-9 (vs the 1e-5 penalty tol in C4.1).
    assert abs(u_tip - u_tip_exact) < 1.0e-9 * u_tip_exact + 1.0e-11, (
        f"augmented tip {u_tip:.8e} != monolithic {u_tip_exact:.8e}")
