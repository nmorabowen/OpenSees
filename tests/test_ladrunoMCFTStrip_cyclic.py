"""LadrunoMCFTStrip -cyclicDegradation (Ladruno WP-256, program #1131): opt-in cyclic strut degradation.

With the option on, the strut stress of the Model C strip carries 1/(1 + Ccyc P), P the committed cumulative
crack re-opening (the increases of (e1 - ecr)+ that stay below the previous maximum opening; zero on a monotonic
path). Off (the default) is Model C's law: every WP-254 gate is unchanged (test_ladrunoMCFTStrip{,_cpp}.py).

The g++ gate (option on vs the reference, self-checks C1-C3) is test_ladrunoMCFTStrip_cyclic_cpp.py. This file:
- OpenSees: the Model C column section with the option on against the reference Section on a cyclic path;
  the `cyclicDamage` response; the parser; a two-round database round trip of the new memory.
"""
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from _testbed import ops
from _testbed.roundtrip import database_roundtrip

sys.path.insert(0, str(Path(__file__).resolve().parent / "_testbed"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import mcft_strip_ref as ref
import test_ladrunoMCFTStrip as W  # the WP-254 OpenSees helpers (model, section, drive, path)

pytestmark = [pytest.mark.zone_a]

CCYC = 300.0


# ------------------------------------------------------------------------------------------ OpenSees
SEC_C, DATA_C = ref.column_section(c_cyc=CCYC)


def _materials_on():
    ops.nDMaterial("LadrunoMCFTStrip", 1, *W._strip_args(W.CORE), "-cyclicDegradation", CCYC)
    ops.nDMaterial("LadrunoMCFTStrip", 2, *W._strip_args(W.COVER), "-cyclicDegradation", CCYC)
    ops.uniaxialMaterial("LadrunoMCFTSteel", 10, W.SEC.bar_steel.fy, W.SEC.bar_steel.fu)
    ops.nDMaterial("LadrunoUniaxialFiber", 3, 10)


def test_section_with_the_option_matches_the_reference():
    W._model()
    _materials_on()
    W._section()
    rows = W._path()
    W._drive(rows)
    st = SEC_C.new_state()
    sc = np.array([SEC_C.conc.fc * SEC_C.b * SEC_C.h, SEC_C.conc.fc * SEC_C.b * SEC_C.h ** 2,
                   SEC_C.conc.fc * SEC_C.b * SEC_C.h])
    worst = np.zeros(3)
    for _ in rows:
        assert ops.analyze(1) == 0
        d = ops.eleResponse(1, "section", "deformation")
        f = ops.eleResponse(1, "section", "force")
        e = (d[0], d[1], d[3])
        s = SEC_C.forces(e, st)
        SEC_C.commit(e, st)
        worst = np.maximum(worst, np.abs(np.array([f[0], f[1], f[3]]) - s) / sc)
    assert worst.max() < W.TOL, f"N, M, V max rel. diff {worst}"
    assert st.reopen.max() > 1e-3                              # the option acted on this path
    # the response reports the committed memory
    y = float(SEC_C.y[20])
    P, bc = ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "cyclicDamage")
    assert P == pytest.approx(float(st.reopen[20]), rel=1e-12, abs=1e-15)
    assert bc == pytest.approx(1.0 / (1.0 + CCYC * P), rel=1e-15)


def test_option_off_and_zero_are_model_c():
    """-cyclicDegradation 0 is the default: the same section forces bit for bit along a cyclic path."""
    out = []
    for flag in ([], ["-cyclicDegradation", 0.0]):
        W._model()
        ops.nDMaterial("LadrunoMCFTStrip", 1, *W._strip_args(W.CORE), *flag)
        ops.nDMaterial("LadrunoMCFTStrip", 2, *W._strip_args(W.COVER), *flag)
        ops.uniaxialMaterial("LadrunoMCFTSteel", 10, W.SEC.bar_steel.fy, W.SEC.bar_steel.fu)
        ops.nDMaterial("LadrunoUniaxialFiber", 3, 10)
        W._section()
        rows = W._path()
        W._drive(rows)
        f = []
        for _ in rows:
            assert ops.analyze(1) == 0
            f.append(list(ops.eleResponse(1, "section", "force")))
        out.append(f)
    assert out[0] == out[1]


@pytest.mark.parametrize("vals", [(-1.0,), (float("inf"),), (3e4, -1.0), (3e4, 5.0)])
def test_parser_refuses(vals):
    W._model()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoMCFTStrip", 1, *W.GOOD, "-cyclicDegradation", *vals)


def test_parser_accepts():
    W._model()
    ops.nDMaterial("LadrunoMCFTStrip", 1, *W.GOOD, "-cyclicDegradation", 0.0)
    ops.nDMaterial("LadrunoMCFTStrip", 2, *W.GOOD, "-cyclicDegradation", 250.0)
    ops.nDMaterial("LadrunoMCFTStrip", 3, *W.GOOD, "-cyclicDegradation", 3e4, 1.0, "-betaC", 170.0)   # with m
    ops.nDMaterial("LadrunoMCFTStrip", 4, *W.GOOD, "-cyclicDegradation", 250.0, "-betaC", 170.0)      # without


def test_serialization_round_trip_with_the_option():
    """The parameter and the new committed memory (crack opening, cumulative re-opening) survive a database
    round trip, two rounds: the restored strips re-solve to the same forces and damage bit for bit."""
    rows = W._path()[:90]
    y = float(SEC_C.y[20])

    def build():
        W._model()
        _materials_on()
        W._section()
        W._drive(rows)
        for _ in rows:
            assert ops.analyze(1) == 0

    def probe():
        P = list(ops.eleResponse(1, "section", "fiber", y, 0.0, 1, "cyclicDamage"))
        assert P[0] > 0.0 and math.isfinite(P[1])
        return list(ops.eleResponse(1, "section", "force")) + P

    database_roundtrip(build, probe_nodes=[2], ndf=6, dbname="mcft_wp256", rounds=2, probe_fn=probe)
