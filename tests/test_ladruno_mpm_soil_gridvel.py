"""WP-232k D17 / D18 (tests/mpm_soil_interface.md): `mpm gridvel ... profile ... translate`, in
both modes, and the mode-F defect #1109 -- the explicit engine accepted `gridvel` and `pforce`
and ignored them. Now mode F applies a gridvel in both projections and refuses `pforce`.

Closed-form checks only: a block whose every grid node is prescribed one velocity moves by
exactly v Phi(t) with zero stress (partition of unity: L = sum v (x) grad S = 0); an elastic column
on a translating floor in mode A drops by D with its self-weight state kept.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _l(v):
    return [float(x) for x in (v if isinstance(v, (list, tuple)) else [v])]


def _phi_integral(t0, tr, t1, t):
    full = (t1 - t0) - tr
    if t <= t0:
        return 0.0
    if t >= t1:
        return full
    if tr > 0 and t < t0 + tr:
        return 0.5 * (t - t0) ** 2 / tr
    if tr > 0 and t > t1 - tr:
        return full - 0.5 * (t1 - t) ** 2 / tr
    return 0.5 * tr + (t - t0 - tr)


def _block_f(ops, gridvel=True):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 100.0, 0.3)
    ops.mpm("grid", 232, "-origin", 0.0, -2.0, "-h", 0.5, "-cells", 4, 8)
    ops.mpm("body", 1, "box", 0.5, 0.5, 1.5, 1.5, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    if gridvel:
        ops.mpm("gridvel", "all", "box", -1.0, -3.0, 3.0, 3.0, "v", 0.2, -0.1,
                "profile", 0.1, 0.2, 1.0)


def test_modef_gridvel_moves_block_exactly():
    """#1109: every node prescribed v * phi(t) -> every particle moves v * Phi(t), sigma = 0."""
    ops = _ops()
    _block_f(ops)
    x0 = _l(ops.mpm("particleData", "x"))
    ops.mpm("analyze", 120, 0.01)
    x1 = _l(ops.mpm("particleData", "x"))
    t = float(ops.mpm("info", "time"))
    P = _phi_integral(0.1, 0.2, 1.0, t)
    assert abs(P - 0.7) < 1e-9
    for k in range(0, len(x0), 2):
        assert abs(x1[k] - x0[k] - 0.2 * P) < 1e-10
        assert abs(x1[k + 1] - x0[k + 1] + 0.1 * P) < 1e-10
    assert max(abs(s) for s in _l(ops.mpm("particleData", "sig"))) < 1e-6   # E = 100: strain 1e-8


def test_modef_without_gridvel_unchanged():
    """No gridvel: the block stays put (no gravity), as before #1109's fix."""
    ops = _ops()
    _block_f(ops, gridvel=False)
    x0 = _l(ops.mpm("particleData", "x"))
    ops.mpm("analyze", 20, 0.01)
    assert _l(ops.mpm("particleData", "x")) == x0


def test_modef_refuses_pforce():
    """#1109: a mode-F pforce was accepted and ignored; now `mpm analyze` refuses it."""
    ops = _ops()
    _block_f(ops, gridvel=False)
    ops.mpm("pforce", 1, "box", 0.0, 0.0, 2.0, 2.0, 1.0, 0.0)
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 0.01)


def test_modef_condforce_refused():
    ops = _ops()
    _block_f(ops)
    ops.mpm("analyze", 1, 0.01)
    with pytest.raises(Exception):
        ops.mpm("condForce", "all")


def test_profile_argument_checks():
    ops = _ops()
    _block_f(ops, gridvel=False)
    with pytest.raises(Exception):     # t1 < t0 + 2 Tr
        ops.mpm("gridvel", "g1", "box", 0.0, 0.0, 1.0, 1.0, "v", 0.0, 1.0, "profile", 0.0, 1.0, 1.5)
    with pytest.raises(Exception):     # profile replaces start
        ops.mpm("gridvel", "g2", "box", 0.0, 0.0, 1.0, 1.0, "v", 0.0, 1.0, "start", 0.0,
                "profile", 0.0, 0.0, 1.0)
    with pytest.raises(Exception):     # track and translate exclude each other
        ops.mpm("gridvel", "g3", "box", 0.0, 0.0, 1.0, 1.0, "v", 0.0, 1.0, "track", "y",
                "translate", "y")


def test_modea_translating_floor_carries_column():
    """D17, mode A (with the D19a floor block): an elastic column on a 2-cell elastic block whose
    lower rows a `profile ... translate y` gridvel carries; the prescribed rows stay 2 cells below
    the column, so their entries and exits never touch a column particle's support. Between the
    end of gravity and the end of the translation the mean column drop is D (to the column's
    compression change, rho g l / E ~ 2e-4), and the stress keeps its self-weight state to the
    cell-crossing error -- a smoke test of the mechanism."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 2.0e5, 0.3)
    ops.mpm("grid", 232, "-origin", 0.0, -2.5, "-h", 0.5, "-cells", 2, 12)
    ops.mpm("mode", "quasiStatic")
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("body", 2, "box", 0.0, -1.0, 1.0, 0.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("gridvel", "floor", "box", 0.0, -1.5, 1.0, -0.5, "v", "free", -1.0,
            "profile", 0.5, 0.0, 1.0, "translate", "y")
    ops.mpm("gravity", 0.0, -9.81)
    ops.timeSeries("Path", 9001, "-time", 0.0, 0.4, 1.0e9, "-values", 0.0, 1.0, 1.0)
    ops.mpm("timeSeries", 9001)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-12, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    ops.mpm("analyze", 5, 0.1)
    n = 16                                   # body 1 is seeded first: its 16 particles lead
    xa = _l(ops.mpm("particleData", "x"))[:2 * n]
    sa = _l(ops.mpm("particleData", "sig"))[:3 * n]
    ops.mpm("analyze", 8, 0.0625)
    ops.mpm("analyze", 2, 0.1)
    xb = _l(ops.mpm("particleData", "x"))[:2 * n]
    sb = _l(ops.mpm("particleData", "sig"))[:3 * n]
    drop = sum(xb[2 * k + 1] - xa[2 * k + 1] for k in range(n)) / n
    assert abs(drop + 0.5) < 2e-4
    scale = max(abs(s) for s in sa)
    assert max(abs(a - b) for a, b in zip(sa, sb)) < 0.05 * scale
