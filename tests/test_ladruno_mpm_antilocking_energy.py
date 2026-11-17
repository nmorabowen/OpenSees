"""WP-232j, gate row L9 (critic round C1): energy with anti-locking on.

The critic's probe, pre-registered in wp-briefs.md "232j" (critic block), as the `kind energy`
cases of tests/mpm_antilocking_cases.txt: a free elastic block on a fixed base under gravity,
run off and on; KE_max(on) <= 1.5 KE_max(off) and every step taken (no particle leaves the
grid). Mode F at 0.15 dt_crit (4000 steps) and with v0 = 1 m/s at 0.05 dt_crit (9000 steps);
mode A with Newmark average acceleration, 400 steps. Plus the S1 energy ledger (bar_flip_2d,
impact_discs_2d) with anti-locking on, at their frozen tolerances.

Zone-A runs every mode-F row (96 particles each, about 1-2 s a run) and the ledger; the
mode-A rows run where the warrant is assembled (LADRUNO_MPM_GATES=1), except gimp at
nu = 0.45, which stays in Zone-A. Each row prints its wall time.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_antilocking_harness as A  # noqa: E402

pytestmark = pytest.mark.zone_a

_GATES = os.environ.get("LADRUNO_MPM_GATES") == "1"
ENERGY = [c["name"] for c in A.parse_cases() if c["kind"] == "energy"]
LEDGER = [c["name"] for c in A.parse_cases() if c["kind"] == "ledger"]


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


@pytest.mark.parametrize("name", ENERGY)
def test_energy_l9(name):
    c = A.case_by_name(name)
    if c["energy_mode"] == "A" and name != "energy_A_gimp_nu045" and not _GATES:
        pytest.skip("mode-A L9 rows other than gimp nu 0.45 run with LADRUNO_MPM_GATES=1 (warrant)")
    r = A.run_energy(c, _ops())
    vals = {"ke_ratio": r["ke_ratio"], "completed": r["completed"]}
    rows = A.evaluate_values(c, vals)
    print(f"\n{name}: KE_max off {r['off']['ke_max']:.4g} ({r['off']['done']} steps), on "
          f"{r['on']['ke_max']:.4g} ({r['on']['done']} steps), ratio {r['ke_ratio']:.4f}, "
          f"dt {r['on']['dt']:.4g}, wall {r['wall_s']:.1f} s")
    for lab, v, det in rows:
        print(f"{v.upper():6s} {lab}: {det}")
    assert all(v in ("pass", "report") for _, v, _ in rows), rows


@pytest.mark.parametrize("name", LEDGER)
def test_energy_ledger_l9(name):
    c = A.case_by_name(name)
    r = A.run_ledger(c, _ops())
    vals = {k[0]: v[0] for k, v in r["obs"].items() if len(k) == 1 and v}
    rows = A.evaluate_values(c, vals)
    print(f"\n{name}: averaged / fallback particle-steps {r['counts']}")
    for lab, v, det in rows:
        print(f"{v.upper():6s} {lab}: {det}")
    assert r["stopped"] is None, r["stopped"]
    assert r["counts"][0] > 0, r["counts"]               # anti-locking was on
    assert all(v in ("pass", "report") for _, v, _ in rows), rows
