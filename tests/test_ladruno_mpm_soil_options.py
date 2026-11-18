"""WP-232k B3 (the AMEND round): S1 one-point legs for the admitted OPTIONS -- SANISAND IntScheme
2 and 129 and `-implex`; HS `-implex`, `-smallStrain` (HSSmall cyclic shear) and `-sigmat`.
References: 232k-V's one-point FE record tests/data/mpm/soil/s1_onepoint_options.json, whose
legs carry their own deck (`cmd`) beside the S1 fields (F, dF, steps, fe_self_consistency).
Same harness and tolerance rule as S1. An option with no passing leg is refused at admission
instead (the round's rule), so every leg here must pass.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

pytestmark = pytest.mark.zone_a
REC = Path(__file__).resolve().parent / "data" / "mpm" / "soil" / "s1_onepoint_options.json"


def _legs():
    if not REC.exists():
        return []
    return sorted(json.loads(REC.read_text(encoding="utf-8"))["legs"])


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


LEGS = _legs()


@pytest.mark.skipif(not LEGS, reason="s1_onepoint_options.json not committed yet (232k-V)")
@pytest.mark.parametrize("mode", ["A", "F"])
@pytest.mark.parametrize("leg", LEGS or ["none"])
def test_option_leg(leg, mode):
    rec = json.loads(REC.read_text(encoding="utf-8"))
    out = S.s1_run(_ops(), leg, mode, rec=rec)
    bad = [r for r in out["steps"] if not (r["err_sig"] <= r["tol"] and r["err_state"] <= r["tol"])]
    assert not bad, (leg, mode, bad[0], out["max_err"])
