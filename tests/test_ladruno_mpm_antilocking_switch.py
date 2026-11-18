"""WP-232j, gate row L7: the Zone-A switch oracle (freeze amendment 3).

footing_j2_r1 (L3's r1 rung: rough strip footing on weightless von Mises soil, mode A
quasi-static, gimp, ppc 2, 1152 particles, 60 steps, K6's ladder) run on and off.
On must lie within +-10 % of the FE B-bar q_lim at r1 (570.10 kPa: [513.09, 627.11]); off must
lie outside it. The oracle goes red if the averaging is the identity (both legs would agree).
Each leg reports freeze F1's counts and its wall time (freeze Q4: both legs together are meant
to stay under about 3 minutes on Zone-A; the measured total is printed).
"""
from __future__ import annotations

import os
import sys
import time

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_antilocking_harness as A  # noqa: E402

pytestmark = pytest.mark.zone_a

_T0 = time.time()
_WALL: dict = {}


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


@pytest.mark.parametrize("name", ["switch_j2_r1_on", "switch_j2_r1_off"])
def test_switch_oracle(name):
    c = A.case_by_name(name)
    r = A.run_ladder(c, _ops(), log=A.flushprint)
    _WALL[name] = r["wall_s"]
    rows = A.evaluate(c, r)
    print(f"\n{name}: build {r['ladrunoBuild']}, wall {r['wall_s']:.1f} s, failed attempts "
          f"{r['failed_attempts']}, relaxed {r['relaxed_steps']}, halvings {r['halvings']}, "
          f"averaged {r['al_averaged']}, fallback {r['al_fallback']}, plastic-zone fallback share "
          f"{r['zone_fallback_share']:.4f}; both legs so far {sum(_WALL.values()):.1f} s")
    for lab, v, det in rows:
        print(f"{v.upper():6s} {lab}: {det}")
    assert r["stopped"] is None, f"{name}: stopped at step {r['stopped']}"
    assert all(v in ("pass", "report") for _, v, _ in rows), rows
