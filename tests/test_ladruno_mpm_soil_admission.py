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
    # bspline2: SANISAND particles in mode F need a B-spline basis (the owner's K0 decision)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.mpm("grid", 232, "-origin", 0.0, -1.0, "-h", 0.5, "-cells", 2, 6, "-basis", "bspline2")
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


_D2_CHILD = '''
import sys
sys.path.insert(0, {here!r})
from _testbed import ops
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.mpm("grid", 232, "-origin", 0.0, -1.0, "-h", 0.5, "-cells", 2, 6, "-basis", "bspline2")
ops.mpm("mode", "explicit")
ops.nDMaterial("LadrunoSANISAND", 3, *{mat!r}, 3)
try:
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 3, "-rho", 2.0)
    print("D2-ADMITTED")
except Exception:
    print("D2-REFUSED")
'''


def test_d2_refuses_unadmitted_scheme():
    """IntScheme 4 and 3 (no error control) are refused at particles; 1, 2 and 129 are admitted.

    Scheme 3 runs in a child process: building a scheme-3/5 ManzariDafalias spends the base
    class's once-per-process "no error estimate" warning latch, which
    tests/test_manzari_safety_pack.py::test_scheme3_warns_no_error_control asserts later in the
    same pytest process (Zone-A red at ffec57a09; the same trap as LadrunoSANISAND.cpp's ADR-92
    note)."""
    ops = _ops()
    _grid(ops)
    ops.nDMaterial("LadrunoSANISAND", 3, *SANISAND, 4)
    with pytest.raises(Exception):
        ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 3, "-rho", 2.0)
    from _testbed.subprocess_run import run_python_script
    here = os.path.dirname(os.path.abspath(__file__))
    rc, out = run_python_script(_D2_CHILD.format(here=here, mat=SANISAND), cwd=here, timeout=300)
    assert "D2-REFUSED" in out and "not admitted at particles" in out, (rc, out[-2000:])
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


def _san_block(ops, basis, mode):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 2, 4, "-basis", basis)
    ops.mpm("mode", mode)
    ops.nDMaterial("LadrunoSANISAND", 3, *SANISAND)
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 3, "-rho", 2.0)
    ops.mpm("fix", "xmin", "x")
    ops.mpm("fix", "xmax", "x")
    ops.mpm("fix", "ymin", "y")
    ops.mpm("gravity", 0.0, -9.81)


@pytest.mark.parametrize("basis", ["gimp", "linear"])
def test_s7_k0_sanisand_mode_f_needs_bspline(basis, capfd):
    """The owner's K0 decision: SANISAND particles in explicit mode F are refused on gimp and
    linear at `mpm analyze` (they fail the cell-crossing meter, A = 7.95 on gimp), with K0 named."""
    ops = _ops()
    _san_block(ops, basis, "explicit")
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-4)
    cap = capfd.readouterr()
    assert "K0" in (cap.err + cap.out) and "B-spline" in (cap.err + cap.out)


def test_s7_k0_sanisand_admitted_where_decided():
    """...and admitted in mode F on bspline2, in mode A on gimp; HS stays admitted on gimp in F."""
    ops = _ops()
    _san_block(ops, "bspline2", "explicit")
    assert int(ops.mpm("analyze", 1, 1.0e-4)) == 0
    _san_block(ops, "gimp", "quasiStatic")
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    assert int(ops.mpm("analyze", 1, 0.1)) == 0
    _column(ops, HS, 4, kind="LadrunoHS")
    assert int(ops.mpm("analyze", 1, 1.0e-4)) == 0


# --------------------------------------------------------------------------- #
#  The AMEND round (critic at 926d27153): A1-A5                                 #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mat", ["san", "hs"])
def test_a1_no_particles_after_a_nonzero_walk(mat):
    """A1: a particle added after its tag was walked to a non-zero stage would miss the walk;
    it is refused ('add all bodies before `mpm stage`'). Stage 0 does not lock it."""
    ops = _ops()
    if mat == "san":
        _column(ops, SANISAND, 3)
        tag = 3
    else:
        _column(ops, HS, 4, kind="LadrunoHS")
        tag = 4
    ops.mpm("stage", tag, 0)
    ops.mpm("body", 2, "box", 0.0, 1.0, 1.0, 1.5, "-ppc", 2, "-mat", tag, "-rho", 2.0)
    ops.mpm("stage", tag, 1)
    with pytest.raises(Exception):
        ops.mpm("body", 3, "box", 0.0, 1.5, 1.0, 2.0, "-ppc", 2, "-mat", tag, "-rho", 2.0)


def test_a2_stage_values():
    """A2: SANISAND and HS know stages 0 and 1 only; DruckerPrager keeps 0/1/2."""
    ops = _ops()
    _column(ops, SANISAND, 3)
    with pytest.raises(Exception):
        ops.mpm("stage", 3, 2)
    _column(ops, HS, 4, kind="LadrunoHS")
    with pytest.raises(Exception):
        ops.mpm("stage", 4, 2)
    _grid(ops)
    ops.nDMaterial("DruckerPrager", 2, 1.5e5, 15517.24, 0.2, 0.21, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0)
    ops.mpm("body", 1, "box", 0.0, 0.0, 1.0, 1.0, "-ppc", 2, "-mat", 2, "-rho", 2.0)
    assert int(ops.mpm("stage", 2, 2)) == 16


def test_a3_guard_after_reset_says_stage_0(capfd):
    """A3: after `mpm reset` the particles are at stage 0 under a flipped flag; the guard asks for
    `mpm stage <tag> 0` and never for a flip the particles did not take."""
    ops = _ops()
    _column(ops, SANISAND, 3)
    ops.mpm("stage", 3, 1)
    ops.mpm("analyze", 2, 1.0e-4)
    ops.mpm("reset")
    capfd.readouterr()
    with pytest.raises(Exception):
        ops.mpm("analyze", 1, 1.0e-4)
    cap = capfd.readouterr()
    msg = cap.err + cap.out
    assert "mpm stage <tag> 0" in msg and "mpm stage <tag> 1" not in msg, msg[-1500:]
    ops.mpm("stage", 3, 0)
    assert int(ops.mpm("analyze", 1, 1.0e-4)) == 0


def test_a4_old_pack_refused_with_soil():
    """A4: restart packs 2 and 3 carry no stage record: refused with soil particles."""
    ops = _ops()
    _column(ops, HS, 4, mode="quasiStatic", kind="LadrunoHS")
    with pytest.raises(Exception):
        ops.mpm("debug", "packVersion", 3)
    with pytest.raises(Exception):
        ops.mpm("debug", "packVersion", 2)


def test_a5_moving_box_into_fix_ghost_refused():
    """A5: on a B-spline grid, a `translate` box that is legal at the start but moves into the
    ghost layer behind `fix ymin` is refused at the step it gets there (B10 where the box IS)."""
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 100.0, 0.3)
    ops.mpm("grid", 232, "-origin", 0.0, 0.0, "-h", 0.5, "-cells", 4, 6, "-basis", "bspline2")
    ops.mpm("body", 1, "box", 0.5, 1.0, 1.5, 2.0, "-ppc", 2, "-mat", 1, "-rho", 1.0)
    ops.mpm("fix", "ymin", "y")
    # the box's bottom starts at y = 0.1 (inside the grid) and comes down at 1 m/s
    ops.mpm("gridvel", "g", "box", 0.0, 0.1, 2.0, 0.6, "v", 0.0, -1.0, "profile", 0.0, 0.0, 10.0,
            "translate", "y")
    ops.mpm("analyze", 5, 0.01)            # bottom at 0.05: still inside
    with pytest.raises(Exception):
        ops.mpm("analyze", 10, 0.01)       # crosses y = 0 at t = 0.1
    assert float(ops.mpm("info", "time")) <= 0.115     # refused at the step it gets there
