"""WP-232j, gate row L2: the mode-A particle tangent with anti-locking on (interface J10).

`mpm debug tangentCheck` compares each particle element's getTangentStiff with central
differences of its resisting force (232f's diagnosis hook), on `tangent_column_on` of
tests/mpm_antilocking_cases.txt: undeformed, and after the 10 quasi-static steps that bring the
column to its 0.2 % strain. The strained state is checked on and with `-antiLocking off` on the
same elements. FD convergence: the first eps of 1e-8, 1e-9, 1e-10 that agrees with the next to
10 % (as 232f measured). Gate (wp-briefs.md "232j" L2): undeformed <= 1e-5; strained on <= 2x off.
Cost on Zone-A: about 15 s for the seven rows (80-160 particles each).
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_antilocking_harness as A  # noqa: E402

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


TANGENT = [c["name"] for c in A.parse_cases() if c["kind"] == "tangent"]


@pytest.mark.parametrize("name", TANGENT)
def test_tangent_l2(name):
    """Every L2 row (round 3 C5 adds nu 0.3 / 0.45, 3D, the B8 image path, 1 % and 5 %)."""
    c = A.case_by_name(name)
    r = A.run_tangent(c, _ops())
    und_on = A._converged(r["undeformed_on"])
    und_off = A._converged(r["undeformed_off"])
    str_on = A._converged(r["strained_on"])
    str_off = A._converged(r["strained_off"])
    print(f"\neps {r['eps']}; max |epsHat| at the strained state {r['strain_max']:.3e}")
    print(f"undeformed: on {und_on:.3e} (off {und_off:.3e}); strained: on {str_on:.3e}, "
          f"off {str_off:.3e}, ratio {str_on / max(str_off, 1e-300):.3f}")
    for key in ("undeformed_on", "strained_on", "strained_off"):
        assert all(row[4] == r["N"] for row in r[key]), (key, r[key])  # every element checked
    assert r["strain_max"] >= 0.75 * c["tangent_strain"], (r["strain_max"], c["tangent_strain"])
    gates = {e["obs"]: e["vals"][0] for e in c["expects"]}
    assert und_on <= gates["tangent_rel_undeformed"], und_on
    assert str_on <= gates["tangent_ratio_strained"] * str_off, (str_on, str_off)
