"""LadrunoUniaxialFiber -dowel (Ladruno WP-257, program #1131, improvement I3): dowel action of a bar.

``nDMaterial LadrunoUniaxialFiber $tag $uniTag -dowel $fc $fy $db $ls [-dowelK $kd] [-dowelSlip $su] [-dowelDecay $Cd]``
gives the bar fibre a shear stress per shear component from the crack slip ``s = |gamma| ls``:
Dulacska / Vintzeleou-Tassios capacity ``kd (4/pi) sqrt(fc fy) sqrt(1 - (sigma/fy)^2)`` per unit bar
area, an exponential envelope to 95 % at ``su`` (default 0.05 db), an origin-oriented secant inside the
largest committed slip. Off (the default) the bar carries no shear, as before (WP-245's tests).

The law is gated in g++ (test_ladrunoDowel_cpp.py). This file is the wiring: one bar fibre on a
zeroLengthSection against the numpy transcription, the axial coupling through the uniaxial law, the
response, the parser, and the database round trip.
"""
import math

import numpy as np
import pytest

from _testbed import ops
from _testbed.roundtrip import database_roundtrip
from test_ladrunoDowel_cpp import dowel_ref

pytestmark = [pytest.mark.zone_a]

FC, FY, DB, LS = 33.0, 331.0, 25.4, 150.0
ES = 200000.0
AB = math.pi * DB ** 2 / 4.0
DOF_OF = {0: 1, 1: 6, 2: 5, 3: 2, 4: 3, 5: 4}


def _model():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)


def _section(*dowel, uni="Elastic"):
    if uni == "Elastic":
        ops.uniaxialMaterial("Elastic", 5, ES)
    else:
        ops.uniaxialMaterial("Steel01", 5, FY, ES, 0.0)
    ops.nDMaterial("LadrunoUniaxialFiber", 6, 5, *dowel)
    ops.section("LadrunoNDFiber", 1, "-shearProfile", "uniform")
    ops.fiber(0.0, 0.0, AB, 6)


def _drive(rows):
    t = [0.0] + [float(k + 1) for k in range(len(rows))]
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.element("zeroLengthSection", 1, 1, 2, 1)
    for c in range(6):
        ops.timeSeries("Path", c + 1, "-time", *t, "-values", *([0.0] + [r[c] for r in rows]))
        ops.pattern("Plain", c + 1, c + 1)
        ops.sp(2, DOF_OF[c], 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 30, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


def _path():
    gy = np.concatenate([np.linspace(0, 4e-3, 9)[1:], np.linspace(4e-3, -2e-3, 13)[1:], np.linspace(-2e-3, 6e-3, 17)[1:]])
    gz = 0.3 * gy[::-1]
    e0 = 1.2e-3 * np.sin(np.arange(gy.size) / 4.0)          # |sigma| <= 240 MPa < fy
    return [(float(a), 0.0, 0.0, float(b), float(c), 0.0) for a, b, c in zip(e0, gy, gz)]


@pytest.mark.parametrize("extra, kd, su, cd",
                         [((), 1.3, -1.0, 0.0), (("-dowelK", 0.3, "-dowelSlip", 2.0), 0.3, 2.0, 0.0),
                          (("-dowelK", 0.3, "-dowelDecay", 0.8), 0.3, -1.0, 0.8)],
                         ids=["defaults", "kd_su", "decay"])
def test_section_shear_is_the_dowel_law(extra, kd, su, cd):
    """V_y, V_z = A_b tau(gamma, sigma) on a cyclic path (each step committed), both components."""
    rows = _path()
    _model()
    _section("-dowel", FC, FY, DB, LS, *extra)
    _drive(rows)
    got = []
    for _ in rows:
        assert ops.analyze(1) == 0
        f = list(ops.eleResponse(1, "section", "force"))
        got.append((f[0], f[3], f[4]))
    got = np.array(got)
    sig = ES * np.array([r[0] for r in rows])
    for k, comp in ((1, 3), (2, 4)):
        exp = AB * dowel_ref(FC, FY, DB, LS, kd, su, [(r[comp], s, 1) for r, s in zip(rows, sig)], cd=cd)
        assert np.abs(got[:, k] - exp).max() < 1e-9 * np.abs(exp).max()
    assert got[:, 0] == pytest.approx(AB * sig, rel=1e-12)


def test_off_carries_no_shear_and_the_same_axial_force():
    rows = _path()[:10]
    res = {}
    for on in (False, True):
        _model()
        _section(*(("-dowel", FC, FY, DB, LS) if on else ()))
        _drive(rows)
        for _ in rows:
            assert ops.analyze(1) == 0
        res[on] = list(ops.eleResponse(1, "section", "force"))
    assert res[False][3] == 0.0 and res[False][4] == 0.0
    assert res[True][3] != 0.0
    assert res[True][0] == res[False][0]


def test_yielded_bar_has_no_dowel_and_the_response():
    """Steel01 past yield: sigma = fy, so the capacity and the dowel stress are zero; the response
    `dowel` = (tau_y, tau_z, smax_y, smax_z, cap)."""
    ey = FY / ES
    rows = [(0.5 * ey, 0, 0, 2e-3, 0, 0), (3.0 * ey, 0, 0, 3e-3, 0, 0)]
    _model()
    _section("-dowel", FC, FY, DB, LS, uni="Steel01")
    _drive(rows)
    assert ops.analyze(1) == 0
    d = list(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 6, "dowel"))
    a = 0.5
    cap = 1.3 * 4 / math.pi * math.sqrt(FC * FY) * math.sqrt(1 - a * a)
    assert d[4] == pytest.approx(cap, rel=1e-12)
    assert d[0] == pytest.approx(cap * (1 - math.exp(-2e-3 * LS / (0.05 * DB / 3))), rel=1e-12)
    assert d[2] == pytest.approx(2e-3 * LS, rel=1e-12)    # committed slip
    assert ops.analyze(1) == 0
    d = list(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 6, "dowel"))
    assert d[0] == 0.0 and d[4] == 0.0
    assert list(ops.eleResponse(1, "section", "force"))[3] == 0.0


def test_parser():
    _model()
    ops.uniaxialMaterial("Elastic", 5, ES)
    ops.nDMaterial("LadrunoUniaxialFiber", 6, 5, "-dowel", FC, FY, DB, LS, "-dowelK", 0.5, "-dowelSlip", 1.0)
    ops.nDMaterial("LadrunoUniaxialFiber", 7, 5, "-rho", 7.85e-9, "-dowel", FC, FY, DB, LS)
    bad = [("-dowel", FC, FY, DB),                       # too few values
           ("-dowel", FC, FY, DB, 0.0),                  # ls must be > 0
           ("-dowel", -FC, FY, DB, LS),
           ("-dowelK", 0.5),                             # without -dowel
           ("-dowel", FC, FY, DB, LS, "-dowelK", 0.0),
           ("-dowel", FC, FY, DB, LS, "-dowelSlip", -1.0),
           ("-dowelDecay", 0.5),                         # without -dowel
           ("-dowel", FC, FY, DB, LS, "-dowelDecay", -0.1)]
    for i, b in enumerate(bad):
        with pytest.raises(Exception):
            ops.nDMaterial("LadrunoUniaxialFiber", 20 + i, 5, *b)


def test_serialization_round_trip():
    """The dowel parameters and the committed slips survive a database round trip (two rounds): after
    the restore the fibre is on the secant of its largest slip, not on the virgin envelope."""
    rows = _path()[:20]

    def build():
        _model()
        _section("-dowel", FC, FY, DB, LS, "-dowelK", 0.4, "-dowelSlip", 1.5, "-dowelDecay", 0.6)
        _drive(rows)
        for _ in rows:
            assert ops.analyze(1) == 0

    database_roundtrip(build, probe_nodes=[2], ndf=6, dbname="dowel_wp257", rounds=2,
                       probe_fn=lambda: list(ops.eleResponse(1, "section", "force"))
                       + list(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 6, "dowel")))
