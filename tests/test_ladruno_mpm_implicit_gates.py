"""WP-232f: the fine-resolution P1 gates of MPM mode A (tests/mpm_implicit_cases.txt).

These are the pre-registered gate rows at their finest resolution (cantilever QS and transient
r2 / r3, the mild-cone footing r2 / r3, the plate with a hole): 10^4 - 5 * 10^4 particles over
up to 400 implicit steps, too long for Zone-A's 40-minute job. They run where the warrant is
assembled: set LADRUNO_MPM_GATES=1 (and dist/bin/opensees.pyd on PYTHONPATH). Without it every
test here is skipped with that reason — they are never silently passed.

An expected value still `TBD` (232f-V filling it) is skipped loudly, never passed.
Written from tests/mpm_implicit_interface.md and the cases file only.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_implicit_harness as H  # noqa: E402

pytestmark = [pytest.mark.zone_a,
              pytest.mark.skipif(os.environ.get("LADRUNO_MPM_GATES") != "1",
                                 reason="fine-resolution MPM gates: set LADRUNO_MPM_GATES=1 "
                                        "(too long for Zone-A)")]

GATES = ["cantilever_qs_r2", "cantilever_qs_r3", "cantilever_dyn_r2", "cantilever_dyn_r3",
         "footing_mildcone_r2", "footing_mildcone_r3", "plate_hole_dw050"]

_CACHE: dict = {}


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    return ops


def run(name):
    if name not in _CACHE:
        _CACHE[name] = H.run_python(H.case_by_name(name), _ops())
    return _CACHE[name]


@pytest.mark.parametrize("name", GATES)
def test_mode_a_gate(name):
    c = H.case_by_name(name)
    r = run(name)
    res = H.evaluate(c, r, runs=run)
    for lab, v, det in res:
        print(f"{v.upper():6s} {lab}: {det}")
    failed = [f"{lab}: {det}" for lab, v, det in res if v == "fail"]
    assert not failed, "\n".join(failed)
    tbd = [lab for lab, v, _ in res if v == "tbd"]
    if tbd and all(v in ("tbd", "report") for _, v, _ in res):
        pytest.skip("expected values TBD (232f-V): " + "; ".join(tbd))
