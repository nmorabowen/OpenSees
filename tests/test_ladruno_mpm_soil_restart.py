"""WP-232k S6 (mode A restart): the database round trip with soil particles mid-stage-1. The S2
column (fixed base) is staged, flipped, and loaded by part of the surcharge; then `save`, `wipe`,
`restore`, the analysis rebuilt, and the rest of the surcharge applied. The continued run equals
the uninterrupted one bit for bit (positions, stress, the material state). For SANISAND the
restored host re-asserts the class-wide stage flag that `wipe` zeroed (restart pack v4); the
materials restore their own flipped state (LadrunoSANISAND.cpp:2659-2667).
"""
from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

pytestmark = pytest.mark.zone_a
L0 = 4.0


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _analysis(ops):
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 25)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")


def _surcharge(ops, i, n=6, q=50.0):
    y = S._l(ops.mpm("particleData", "x"))[1::2]
    m = S._l(ops.mpm("particleData", "m"))
    ytop = max(y)
    bf = []
    for k in range(len(m)):
        f = -(q * 0.25) * i / n if ytop - y[k] < 0.125 else 0.0
        bf += [0.0, f / m[k]]
    ops.mpm("particleBodyForce", *bf)
    ops.mpm("analyze", 1, 0.1)


def _fields(ops):
    out = {f: S._l(ops.mpm("particleData", f)) for f in ("x", "sig", "F", "epsHat")}
    out["state"] = S._l(ops.mpm("particleResponse", "state"))
    return out


@pytest.mark.parametrize("mat", ["hs", "sanisand"])
def test_s6_restart_mid_stage1(mat):
    ops = _ops()
    cmd = S.load("s2_column.json")["decks"][mat]["cmd"]

    def staged():
        S.column_build(ops, cmd, "gimp", "A", l0=L0)
        ops.mpm("analyze", 10, 0.1)
        ops.mpm("stage", 1, 1)
        ops.mpm("analyze", 2, 0.1)
        for i in (1, 2, 3):
            _surcharge(ops, i)

    staged()
    for i in (4, 5, 6):
        _surcharge(ops, i)
    ref = _fields(ops)

    staged()
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        ops.database("File", str(Path(td) / "soildb"))
        ops.save(3)
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.database("File", str(Path(td) / "soildb"))
        ops.restore(3)
        _analysis(ops)
        for i in (4, 5, 6):
            _surcharge(ops, i)
        got = _fields(ops)
    for f in ref:
        assert got[f] == ref[f], f"{mat}: field {f} differs after the restart"


# --------------------------------------------------------------------------- #
#  Review items 1-3 (/code-review high on #1108): the stage record in the pack  #
# --------------------------------------------------------------------------- #
def _save_restore(ops, td):
    ops.database("File", str(Path(td) / "soildb"))
    ops.save(3)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.database("File", str(Path(td) / "soildb"))
    ops.restore(3)
    _analysis(ops)


def test_restart_keeps_hs_add_after_walk_refusal():
    """Item 1: AMEND A1's per-tag HS stage is in the v4 pack: after `mpm stage <hs> 1`, save,
    wipe and restore, a body with that tag is still refused; another material's body is not."""
    ops = _ops()
    cmd = S.load("s2_column.json")["decks"]["hs"]["cmd"]
    S.column_build(ops, cmd, "gimp", "A", l0=L0)
    ops.mpm("analyze", 2, 0.1)
    ops.mpm("stage", 1, 1)
    ops.mpm("analyze", 1, 0.1)
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        _save_restore(ops, td)
        try:                                   # the restore does not re-register nDMaterials
            ops.nDMaterial(*S.deck_args(cmd))
        except Exception:
            pass
        n0 = len(S._l(ops.mpm("particleData", "m")))
        with pytest.raises(Exception):
            ops.mpm("body", 5, "box", 0.0, 0.0, 1.0, 0.5, "-ppc", 2, "-mat", 1, "-rho", 2.0)
        assert len(S._l(ops.mpm("particleData", "m"))) == n0
        ops.nDMaterial("ElasticIsotropic", 77, 1.0e4, 0.3)
        ops.mpm("body", 6, "box", 0.0, 0.0, 1.0, 0.5, "-ppc", 2, "-mat", 77, "-rho", 2.0)
        assert len(S._l(ops.mpm("particleData", "m"))) > n0


def test_pack_version_3_then_soil_writes_v4():
    """Item 3: `debug packVersion 3` set BEFORE the soil body is added used to write a v3 pack
    (no stage record). packState re-checks: it writes v4 with a WARNING, so the SANISAND restart
    continues bit for bit (a v3 pack would lose the class-wide stage and be refused)."""
    ops = _ops()
    cmd = S.load("s2_column.json")["decks"]["sanisand"]["cmd"]

    class _PV3:                                # issue packVersion 3 before the first `mpm body`
        def __init__(self, o):
            self._o, self._done = o, False

        def __getattr__(self, n):
            return getattr(self._o, n)

        def mpm(self, *a):
            if a[:1] == ("body",) and not self._done:
                self._done = True
                self._o.mpm("debug", "packVersion", 3)
            return self._o.mpm(*a)

    def staged(o):
        S.column_build(o, cmd, "gimp", "A", l0=L0)
        ops.mpm("analyze", 10, 0.1)
        ops.mpm("stage", 1, 1)
        ops.mpm("analyze", 2, 0.1)
        for i in (1, 2, 3):
            _surcharge(ops, i)

    staged(ops)
    for i in (4, 5, 6):
        _surcharge(ops, i)
    ref = _fields(ops)
    staged(_PV3(ops))
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        _save_restore(ops, td)
        for i in (4, 5, 6):
            _surcharge(ops, i)
        got = _fields(ops)
    for f in ref:
        assert got[f] == ref[f], f"field {f} differs after the restart"
