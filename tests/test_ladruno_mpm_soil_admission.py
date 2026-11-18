"""WP-232k K1-K4, K7 (tests/mpm_soil_interface.md D1-D16): soil particle admission, the stage
walk and its guard, the symmetric-SOE refusal, the soil queries. Small, closed-form or
behavioural checks; the gated numbers (K0, S1-S5) live in their own files.

Decks (tests/mpm_soil_cases.txt, kPa / t / m / s): SANISAND tests/test_ladruno_sanisand.py:245-263,
HS the cases file's k0_ref_column material 4.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(__file__))

pytestmark = pytest.mark.zone_a

SANISAND = [264.32, 0.3129, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005, 1.3, 0.968,
            3.5, 0.05, 5.75, 12.5, 1100.0, 2.0]
HS = [35000.0, 35000.0, 105000.0, 0.5, 10.0, 38.0, 8.0, "-nu", 0.2, "-pref", 100.0, "-Rf", 0.9,
      "-OCR", 1.5, "-stage", 0]
ND_SANISAND_PS, ND_HS = 33021, 33029


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


def _grid(ops, mode="explicit"):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.mpm("grid", 232, "-origin", 0.0, -1.0, "-h", 0.5, "-cells", 2, 6)
    ops.mpm("mode", mode)


def _column(ops, mat_args, tag=3, mode="explicit", kind="LadrunoSANISAND"):
    _grid(ops, mode)
    ops.nDMaterial(kind, tag, *mat_args)
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", tag, "-rho", 2.0)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y")
    ops.mpm("gravity", 0.0, -9.81)


def test_k1_admits_sanisand_and_hs():
    ops = _ops()
    _column(ops, SANISAND, 3)
    assert int(ops.mpm("info", "numParticles")) == 16
    info = _l(ops.mpm("info", "materials"))
    assert info[:5] == [3.0, ND_SANISAND_PS, 1.0, 0.0, 0.0] and info[5] == 16.0
    _column(ops, HS, 4, kind="LadrunoHS")
    info = _l(ops.mpm("info", "materials"))
    assert info[0] == 4.0 and info[1] == ND_HS


def test_d2_refuses_unadmitted_scheme():
    """IntScheme 3 (no error control) is refused at particles; 1, 2 and 129 are admitted."""
    ops = _ops()
    _grid(ops)
    ops.nDMaterial("LadrunoSANISAND", 3, *SANISAND, 3)
    with pytest.raises(Exception):
        ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 3, "-rho", 2.0)
    for scheme in (1, 2, 129):
        _grid(ops)
        ops.nDMaterial("LadrunoSANISAND", 3, *SANISAND, scheme)
        ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 3, "-rho", 2.0)


def test_d6_refuses_hs_sig0():
    ops = _ops()
    _grid(ops)
    ops.nDMaterial("LadrunoHS", 4, *HS, "-sig0", 50.0)
    with pytest.raises(Exception):
        ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 4, "-rho", 2.0)


def test_d7_stage_walk_counts_and_refusals():
    ops = _ops()
    _column(ops, SANISAND, 3)
    assert int(ops.mpm("stage", 3, 1)) == 16
    assert _l(ops.mpm("info", "materials"))[4] == 1.0          # SANISAND's class-wide flag
    assert int(ops.mpm("stage", 3, 0)) == 16
    with pytest.raises(Exception):                             # no particle material of the tag
        ops.mpm("stage", 99, 1)
    _grid(ops)
    ops.nDMaterial("ElasticIsotropic", 1, 1.0e5, 0.3)
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 1, "-rho", 2.0)
    with pytest.raises(Exception):                             # no stage: never a silent no-op
        ops.mpm("stage", 1, 1)


def test_d9_parameter_route_walks_host():
    """`parameter ... loadPattern <mpm> updateMaterialStage <tag>` + `updateParameter` (a DOUBLE)
    reaches every particle through the host's walk, not SANISAND's theInt read."""
    ops = _ops()
    _column(ops, SANISAND, 3)
    ops.parameter(1, "loadPattern", 232, "updateMaterialStage", 3)
    ops.updateParameter(1, 1.0)
    assert _l(ops.mpm("info", "materials"))[4] == 1.0
    ops.updateParameter(1, 0.0)
    assert _l(ops.mpm("info", "materials"))[4] == 0.0


def test_d10_guard_refuses_a_flip_that_missed_the_particles():
    """An FE SANISAND element (all DOFs fixed) answers `updateMaterialStage` first; the class-wide
    flag flips without the particles: `mpm analyze` is refused until `mpm stage` walks them."""
    ops = _ops()
    _column(ops, SANISAND, 3)
    ops.nDMaterial("LadrunoSANISAND", 7, *SANISAND)
    for i, (x, y) in enumerate([(5.0, 5.0), (6.0, 5.0), (6.0, 6.0), (5.0, 6.0)]):
        ops.node(9001 + i, x, y)
        ops.fix(9001 + i, 1, 1)
    ops.element("quad", 9001, 9001, 9002, 9003, 9004, 1.0, "PlaneStrain", 7)
    ops.updateMaterialStage("-material", 7, "-stage", 1)
    with pytest.raises(Exception):          # the guard runs before any step (0 steps here)
        ops.mpm("analyze", 0, 1.0e-4)
    ops.mpm("stage", 3, 1)
    assert int(ops.mpm("analyze", 0, 1.0e-4)) == 0


def test_d16_sigzz_and_response_sizes():
    ops = _ops()
    _column(ops, SANISAND, 3)
    ops.mpm("analyze", 5, 1.0e-4)
    szz = _l(ops.mpm("particleData", "sigZZ"))
    assert len(szz) == 16 and all(v == v for v in szz)          # no NaN: D4's getStressZZ
    sizes = _l(ops.mpm("particleResponse", "-size", "state"))
    assert sizes == [26.0] * 16
    assert len(_l(ops.mpm("particleResponse", "state"))) == 16 * 26


def test_d15_symmetric_soe_refused_with_soil():
    ops = _ops()
    _column(ops, SANISAND, 3, mode="quasiStatic")
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("ProfileSPD")
    ops.test("NormDispIncr", 1.0e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 0.1)


def test_s7_q2_offset_refused_on_soil():
    """Q2: a soil particle takes no initial-stress offset (K3)."""
    ops = _ops()
    _column(ops, HS, 4, kind="LadrunoHS")
    with pytest.raises(Exception):
        ops.mpm("stress", 1, "uniform", -10.0, -10.0, 0.0)
    with pytest.raises(Exception):
        ops.mpm("stress", 1, "geostatic", 1.0, 2.0, 9.81, 0.5)


def _column_a(ops, threads):
    # the column stands on the grid's bottom face (fix ymin y at y = 0)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 2, 5)
    ops.mpm("mode", "quasiStatic")
    ops.nDMaterial("LadrunoHS", 4, *HS)
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 2.0, "-ppc", 2, "-mat", 4, "-rho", 2.0)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y")
    ops.mpm("gravity", 0.0, -9.81)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    ops.ladrunoThreads(threads)
    ops.mpm("analyze", 3, 0.1)
    ops.mpm("stage", 4, 1)
    ops.mpm("analyze", 2, 0.1)
    return _l(ops.mpm("particleData", "sig"))


def test_s7_k2_threads_honoured_as_one(capfd):
    """K2 (D5): soil particles run serial; ladrunoThreads 4 prints the existing 'not on the
    thread-safe allowlist ... SERIAL' warning and gives the bit-identical answer."""
    ops = _ops()
    one = _column_a(ops, 1)
    capfd.readouterr()
    four = _column_a(ops, 4)
    cap = capfd.readouterr()
    ops.ladrunoThreads(1)
    assert four == one
    assert "SERIAL" in (cap.err + cap.out) or "serial" in (cap.err + cap.out).lower()
