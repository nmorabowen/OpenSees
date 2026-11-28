"""WP-266 (D-UP1): smoke test of the large-strain consolidation Newton deck.

tests/_testbed/up_hypo_consolidation_newton.py is the deck the orchestrator runs
to measure Newton iterations under LadrunoUP -geom hypo (and, in wave 2, with a
full-Jacobian option). This test runs a short version of it on THIS worktree's
engine and pins two facts the WP-266 oracle predicts:

  * the deck converges every step with the shipped tangent, under both the
    hypo battery's Newmark(0.6, 0.3025) and Newmark(0.5, 0.5);
  * Newmark(0.5, 0.5) needs FEWER Newton iterations: it scales the p-row rate
    block by c2 = gamma/(beta dt) = 1/dt, which is exactly what the incremental
    p-row residual differentiates to, while Newmark(0.6, 0.3025) scales it by
    1.98/dt (quirk WP-266-newmark-c2-vs-incremental-p-row).

Run:  py -3.12 -S Ladruno_scripts/ci_run_pytest.py -- tests/test_ladruno_up_hypo_newton_deck.py
"""
import importlib.util
import os
from pathlib import Path

import pytest

_DIST = str(Path(__file__).resolve().parents[1] / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

pytestmark = [pytest.mark.zone_b]

_DECK = Path(__file__).resolve().parent / "_testbed" / "up_hypo_consolidation_newton.py"


def _deck():
    # loaded by path: importing the _testbed package would run its __init__
    spec = importlib.util.spec_from_file_location("up_hypo_consolidation_newton", _DECK)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_consolidation_deck_converges_and_matched_c2_needs_fewer_iterations():
    deck = _deck()
    args = dict(nz=8, lz=4.0, nsteps=6, dt=0.2, qfrac=0.2, tol=1.0e-8, maxit=100)
    its_b, traj_b, _, failed_b, _ = deck.run(ops, 0.6, 0.3025, [], **args)
    its_m, traj_m, _, failed_m, _ = deck.run(ops, 0.5, 0.5, [], **args)
    assert failed_b == 0 and failed_m == 0, (its_b, its_m)
    assert len(its_b) == 6 and len(its_m) == 6
    settle_b, settle_m = traj_b[-2], traj_m[-2]
    assert settle_b < -0.05 and settle_m < -0.05, (settle_b, settle_m)   # a real large-strain push
    assert abs(settle_b - settle_m) <= 1e-2 * abs(settle_m), (settle_b, settle_m)
    assert sum(its_m) < sum(its_b), (
        f"Newmark(0.5,0.5) should need fewer iterations than Newmark(0.6,0.3025): "
        f"{its_m} vs {its_b}")
