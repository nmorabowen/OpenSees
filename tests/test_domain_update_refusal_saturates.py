"""A summed refusal stays a refusal (apeGmsh #2619) -- Zone-A openseespy gate.

Before this slice every summing host (``Tri31``, ``FourNodeQuad``, ``LadrunoQuad``:
``ret += setTrialStrain(...)``) and ``Domain::update`` (``ok += theEle->update()``)
added the per-point material return codes in an ``int``. ``LADRUNO_MATERIAL_REFUSED``
is -33086, so 64,907 refusing points in ONE update pushed the model-wide sum past
INT_MIN; it wrapped to a POSITIVE value, the integrators' ``< 0`` checks accepted it,
and the step was committed as if nothing had happened (Opus review of fork #1038).
A bare -1 cannot reach that count. The fix keeps the most negative code instead of
summing, in the three hosts and in both ``Domain::update`` loops.

The gate: a square mesh of unit plane-strain quads, every Gauss point driven to the
same strain the material refuses (``_PE_FAIL`` of ``test_ladrunoConcrete3D_element``,
imposed as the linear field ``u = eps . x`` by ``sp`` on every node, Penalty
constraints, ``algorithm Linear``), with the point count chosen INSIDE the wrap
window ``64,907 <= nGP <= 129,812``. The pre-fix sum of the codes, computed in int32
below, is positive, so a summing engine returns ``analyze(1) == 0``; a saturating
one refuses (``< 0``). The single-element case first, as the premise: it refuses on
both engines, so a wrong strain is told apart from the wrap.

Revert evidence: the probe run named in the fork PR body runs this file on
``origin/ladruno`` before the fix; the two mesh cases fail there with
``analyze returned 0`` and the premise cases pass.

Wall time: the two mesh cases run ~70k failing return maps each (the material warns
once per refused point, so the captured stderr is a few MB); measured on the CI
probe run, see the PR body.
"""
import time

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# The sentinel (SRC/material/LadrunoMaterialStatus.h) and the int32 wrap window of its sum.
_SENTINEL = -33086
_INT32_MIN = -(2 ** 31)
_N_WRAP = (-_INT32_MIN) // (-_SENTINEL) + 1          # 64,907: the first count whose sum wraps

# oracle-gate concrete (E, nu, fc, ft, Gf, Gc) and the plane-strain strain its return map cannot
# integrate from the virgin state (test_ladrunoConcrete3D_element._PE_FAIL, fork #1038; robust to
# 1e-4 relative perturbations). lch is Element::getCharacteristicLength = the shortest node pair,
# 1.0 on a unit quad exactly as on #1038's unit tri, so the classification carries over.
_CONCRETE = (30000.0, 0.2, 30.0, 3.0, 0.1, 5.0)
_PE_FAIL = (-2.9e-4, 2.6e-4, -4.4e-4)

# 132 x 132 unit quads, 4 points each: 69,696 refusing points, inside the wrap window with a
# margin of ~4,800 points on the low side.
_NX = 132


def _int32(v):
    """Two's-complement int32 wrap of a Python int (what the summed C++ `int` holds)."""
    return (v + 2 ** 31) % 2 ** 32 - 2 ** 31


def _quad(host, tag, n1, n2, n3, n4):
    if host == "quad":
        ops.element("quad", tag, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
    elif host == "LadrunoQuad":
        ops.element("LadrunoQuad", tag, n1, n2, n3, n4, 1, "-type", "PlaneStrain", "-thick", 1.0)
    else:
        raise ValueError(host)


def _imposed_mesh(host, nx):
    """nx x nx unit quads of `host` under the homogeneous plane-strain field u = eps . x on every
    node (Penalty, Constant series, LoadControl 1.0, Linear): one update, every point refuses."""
    exx, eyy, gxy = _PE_FAIL
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("LadrunoConcrete3D", 1, *_CONCRETE)
    nn = nx + 1
    for j in range(nn):
        for i in range(nn):
            ops.node(j * nn + i + 1, float(i), float(j))
    for j in range(nx):
        for i in range(nx):
            n1 = j * nn + i + 1
            _quad(host, j * nx + i + 1, n1, n1 + 1, n1 + nn + 1, n1 + nn)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(nn):
        for i in range(nn):
            x, y = float(i), float(j)
            t = j * nn + i + 1
            ops.sp(t, 1, exx * x + 0.5 * gxy * y)
            ops.sp(t, 2, 0.5 * gxy * x + eyy * y)
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")


@pytest.mark.parametrize("host", ["quad", "LadrunoQuad"])
def test_premise_one_refusing_quad_is_refused(host):
    """The premise, independent of the wrap: one unit quad at _PE_FAIL refuses on any engine
    (4 x -33086 cannot overflow). If this fails the strain is wrong, not the reduction."""
    _imposed_mesh(host, 1)
    rc = ops.analyze(1)
    assert rc < 0, f"{host}: one quad at {_PE_FAIL} was accepted (analyze returned {rc})"


@pytest.mark.parametrize("host", ["quad", "LadrunoQuad"])
def test_wrap_window_refusal_is_still_refused(host):
    """69,696 refusing points in one update: the summed int32 code is POSITIVE (the pre-fix
    engine accepts the step and analyze returns 0); the saturating engine refuses (< 0)."""
    n_gp = 4 * _NX * _NX
    summed = _int32(n_gp * _SENTINEL)
    assert n_gp >= _N_WRAP and summed > 0, (
        f"test premise: {n_gp} points x {_SENTINEL} must wrap positive in int32, got {summed}")

    t0 = time.perf_counter()
    _imposed_mesh(host, _NX)
    rc = ops.analyze(1)
    dt = time.perf_counter() - t0
    print(f"{host}: {_NX}x{_NX} quads, {n_gp} refusing points, analyze={rc}, {dt:.1f} s")
    assert rc < 0, (
        f"{host}: {n_gp} refusing points, every one LADRUNO_MATERIAL_REFUSED, but analyze "
        f"returned {rc}: the codes were SUMMED and wrapped to {summed} (> 0), so the refusal "
        f"was accepted (wall {dt:.1f} s)")
