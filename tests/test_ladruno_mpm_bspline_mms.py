"""WP-232i: the manufactured-solution orders through the shipped binaries — the Zone-A legs.

The gated rows (G2, G3, G4, G5) are defined on the ladder 1/h = 8, 16, 32, 64 (and 128 when a run
costs under 10 min); the full ladder is the warrant run (tests/mpm_bspline_warrant.py, launched by
the orchestrator, results under tests/data/mpm/bspline/, checked by
tests/test_mpm_bspline_warrant_results.py). Here each row runs its three coarsest rungs
(1/h = 8, 16, 32) at the full ladder's dt rule — a regression leg, not the gate — and asserts the
row's own thresholds on them (order over the three rungs, e non-increasing).

Cost (Python, nmb box, 232i build, estimates from the g++ core at the same sizes): each mode-F
row ~10-20 s (the 1/h = 32 rung: 4096 particles x ~370 steps); each mode-A quasi-static row a few
seconds (4 load steps). Zone-A total about 1 min.

This file is in the MPM mutation suite (Ladruno_scripts/mutation_gate.py): ZERO (no internal
force) and SCALE (G2P position increment x1.5) must turn rows red here (G6, with the warrant run's
G2 / G3 against the mutant binaries).
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_bspline_harness as B  # noqa: E402

pytestmark = pytest.mark.zone_a

ROWS = B.parse_rows()
SMOKE = [8, 16, 32]


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


@pytest.mark.parametrize("name", ["G2", "G3", "G4", "G5-G2", "G5-G3", "G5-G4"])
def test_mms_order_on_the_coarse_rungs(name):
    row = ROWS[name]
    full = B.as_list(row["ladder"])
    dt = None
    if row.get("mode", "explicit") == "explicit":
        dt = float(row["dt_rule"]) * B.critical_dt(row, full[-1])
    res = B.ladder(_ops(), row, rungs=SMOKE, dt=dt, validity=False, log=lambda *_: None)
    print(f"{name}: e = {res['e']}, p = {res['p']:.3f}, crossings = {res['crossings']}")
    ok, why = B.verdict(row, res)
    assert ok, f"{name} on 1/h = {SMOKE}: {why}; e = {res['e']}"
