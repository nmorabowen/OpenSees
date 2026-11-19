"""WP-232l route S (mode F): T7 as a recorded miss, and the run-time limit's cap (interface §7c and
its amendment; owner's decision, wp-briefs "232l", "Route S after the consistency fix").

- T7 (the gate text is unchanged: stable at 0.9 x the drained pencil for 2000 steps, unstable or
  refused at 1.5 x) is a strict xfail. Diagnosis (build 5d64922a4, tests/mpm_twophase_routeS_t7scan.py,
  records route_s_results.json "t7_5d64922a4"): with the USL skeleton strain, T6's column is stable
  for 2000 steps up to 0.7 x the drained pencil on gimp and 0.8 x on bspline2, and blows up at
  0.8 x / 0.9 x; K_f = 1 gives the same thresholds (not the u-p coupling).
- The route-S limit is dt <= s * dt_d: s defaults to 0.5, 0.5 < s <= 0.7 runs with a warning, and
  s > 0.7 is refused.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

pytestmark = pytest.mark.zone_a

HERE = Path(__file__).resolve().parent


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _harness():
    spec = importlib.util.spec_from_file_location("ladruno_mpm_twophase_cases_harness",
                                                  HERE / "ladruno_mpm_twophase_cases_harness.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _as_list(v):
    return [float(x) for x in v] if isinstance(v, (list, tuple)) else [float(v)]


def _deck(ops):
    """A small route-S deck: gimp, two-phase, base fixed, drained top (as the refusals file's)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e4, 0.2, 2.0)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 1.0, "-cells", 2, 2, "-basis", "gimp")
    ops.mpm("mode", "explicit")
    ops.mpm("body", 1, "box", 0.0, 0.0, 2.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("twoPhase", 1, "-Kf", 4444.0, "-poro", 0.4, "-rhoF", 1.0, "-perm", 0.01, 0.01)
    ops.mpm("fix", "ymin", "y")
    ops.mpm("pfix", "ymax", 0.0)


@pytest.mark.xfail(strict=True, reason="T7 recorded miss (owner, 2026-11-19): the USL skeleton is stable "
                   "only to 0.7 x (gimp) / 0.8 x (bspline2) the drained pencil; see the module docstring")
@pytest.mark.parametrize("basis", ["bspline2", "gimp"])
def test_t7_as_written(basis):
    ops = _ops()
    r = _harness().run_t7(ops, basis, nsteps=2000)
    assert r["pass"], r


def test_route_s_safety_above_cap_refused():
    ops = _ops()
    _deck(ops)
    ops.mpm("safety", 0.8)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-4)
    assert all(v == 0.0 for v in _as_list(ops.mpm("particleData", "v")))


def test_route_s_safety_up_to_cap_runs():
    ops = _ops()
    _deck(ops)
    ops.mpm("safety", 0.7)
    dtd = _as_list(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
    ops.mpm("analyze", 2, 0.69 * dtd)


def test_route_s_default_limit_is_half_the_drained_pencil():
    ops = _ops()
    _deck(ops)
    dtd = _as_list(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 0.51 * dtd)
    ops.mpm("analyze", 2, 0.49 * dtd)
