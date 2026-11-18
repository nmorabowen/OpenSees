"""WP-232k S2 and S4 (tests/mpm_soil_interface.md; gates in wp-briefs.md "232k" and its freeze).

S2: the K0 column (mode A QS, fixed base) with a surcharge rho g l0 on the top row: syy against
the closed form within 0.03 rho g l0 and sxx RMS against 232k-V's FE column within 0.02 rho g l0,
at the end of stage 0, after the flip and after the surcharge (tests/data/mpm/soil/s2_column.json).
S4: the stage flip is the same per-particle map in mode A and mode F, and the FE's: on the
low-g twin (g 0.5, l0 2; tests/data/mpm/soil/s4_lowg.json) the row with p < Pmin commits
sigma = -Pmin I and alpha = 0 exactly, every other particle keeps its sigma bit for bit, and
alpha_in = alpha (-flipAlphaIn init); >= 1 reset particle in each mode.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("mat", ["hs", "sanisand"])
def test_s2_column_with_surcharge(mat, basis):
    out = S.s2_run(_ops(), mat, basis)
    assert out["passed"], out["checkpoints"]


@pytest.mark.parametrize("mode", ["A", "F"])
def test_s4_flip_map(mode):
    out = S.s4_flip_map(_ops(), mode)
    assert out["reset"] >= 1, out
    assert not out["bad"], out
