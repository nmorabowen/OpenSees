"""WP-232l: the two-phase gates asserted by pytest (review of record, wp-briefs "232l", 2026-11-19).

The bands are read from tests/mpm_twophase_cases.txt (`band` lines) and, for T3 / T6, from the FE
records the cases file names. The harness is tests/ladruno_mpm_twophase_cases_harness.py
(and tests/ladruno_mpm_twophase_harness.py for T0 (iii)).
- zone_a: short forms. T0 (iii); T13 translation (8 legs); T5 route S; the fix-2 patch check;
  T10 bit for bit against cfe478254; T2p, T3, T4p and T6 on one leg to their first gated row.
- zone_a + slow (opt in with --runslow or LADRUNO_RUN_SLOW=1): full forms. T2p, T3, T4p, T6 on every
  leg; T5 route E. (Not zone_b: this repo's zone_b means "needs gmsh" and is skipped without it.)
"""
from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "mpm" / "twophase"
AUTO = ["auto", 0.25]


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


H = _load("ladruno_mpm_twophase_cases_harness", "ladruno_mpm_twophase_cases_harness.py")


def _band(case, key):
    return H.band_value(H.case_bands(case)[key])


# ---------------------------------------------------------------------------- zone_a
@pytest.mark.zone_a
def test_t0_iii_route_s_oracle():
    ops = _ops()
    T = _load("ladruno_mpm_twophase_harness", "ladruno_mpm_twophase_harness.py")
    rec = json.loads((DATA / "t0_oracle.json").read_text())
    tol = _band("t0_routeS", "all")
    rows = T.run_t0_iii(ops, rec)
    assert all(r["max"] <= tol for r in rows), rows


@pytest.mark.zone_a
@pytest.mark.parametrize("mode", ["implicit", "explicit"])
@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("plate", [False, True])
def test_t13_translation_patch(mode, basis, plate):
    ops = _ops()
    r = H.run_t13(ops, mode, basis, plate)
    for k in ("v_rel", "sig_rel", "p_rel", "x_rel"):
        assert r[k] <= _band("t13_translation", k), r


@pytest.mark.zone_a
def test_t5_route_s_checkerboard():
    ops = _ops()
    r = H.run_t5_routeS(ops, AUTO)
    band = H.band_value(H.case_band_line("t5_ref_b4", "route S").split("<=")[1])
    assert r["CB_lap"] <= band, r


@pytest.mark.zone_a
def test_fix2_patch_check():
    ops = _ops()
    with tempfile.TemporaryDirectory() as tmp:
        r = H.patch_fix2(ops, tmp, selfcheck=False)
    assert r["steady"], r
    assert r["grid_coef_err"] <= _band("fix2_patch", "grid_coef_err"), r
    assert r["field_at_particles_err"] <= _band("fix2_patch", "field_at_particles_err"), r


@pytest.mark.zone_a
def test_t10_single_phase_bit_for_bit(tmp_path):
    """T10: the single-phase MPM answers (24 cases) against the cfe478254 baseline. Bit for bit when
    the platform matches the baseline's (Windows / MSVC); otherwise rel 1e-8, a compiler difference,
    reported (as 232k recorded for S8). T10's gated claim is the bit-for-bit result established on
    the Windows warrant run. The child imports opensees through _testbed, as this process does."""
    _ops()
    from _testbed.subprocess_run import run_python_script
    out = tmp_path / "t10.json"
    script = ("import runpy, sys; sys.argv = ['mpm_twophase_t10_single.py', '--testbed', %r]; "
              "runpy.run_path(%r, run_name='__main__')" % (str(out), str(HERE / "mpm_twophase_t10_single.py")))
    rc, text = run_python_script(script, cwd=str(HERE), timeout=600)
    assert rc == 0 and out.exists(), "T10 child failed (rc %d):\n%s" % (rc, text[-4000:])
    new = json.loads(out.read_text())
    base = json.loads((DATA / "t10_baseline_cfe478254.json").read_text())
    cases = [k for k in base if k not in ("build", "platform", "compiler")]
    if sys.platform == base.get("platform"):
        assert [k for k in cases if new.get(k) != base[k]] == []
        return
    worst = 0.0
    for k in cases:
        for f, b in base[k].items():
            a = new[k][f]
            den = max(max((abs(x) for x in b), default=0.0), 1e-300)
            worst = max(worst, max((abs(x - y) for x, y in zip(a, b)), default=0.0) / den)
    assert worst <= 1e-8, "compiler difference, reported: worst rel %.3e" % worst


@pytest.mark.zone_a
def test_t2p_short():
    ops = _ops()
    r = H.run_t2(ops, "gimp", "off", upto=0.1)
    row = r["rows"]["0.1"]
    assert row["p_err_max"] <= _band("t2_ref_terzaghi", "p_err_max@0.1"), row
    assert row["U_err"] <= _band("t2_ref_terzaghi", "U_err@0.1"), row


@pytest.mark.zone_a
def test_t3_short():
    ops = _ops()
    rec = json.loads((DATA / "t3_mandel_fe.json").read_text())
    r = H.run_t3(ops, "gimp", "off", rec, upto=0.04)
    assert all(v["err"] <= v["band"] for v in r["rows"].values()), r


@pytest.mark.zone_a
def test_t4p_short():
    ops = _ops()
    r = H.run_t4(ops, "gimp", "off")
    assert r["p_err"] <= _band("t4_ref_hydro", "hydro_p_err"), r


@pytest.mark.zone_a
def test_t6_short():
    ops = _ops()
    fe = json.loads((DATA / "t6_fe_dt0.005.json").read_text())
    leg = [f for f in fe["terzaghi"] if f["stab"] == "off"][0]
    r = H.run_t6_terzaghi(ops, "gimp", "off", leg, upto=0.1)
    row = r["rows"]["0.1"]
    assert row["p_err_max"] <= row["band_p"] and row["U_err"] <= row["band_U"], row


# ---------------------------------------------------------------------------- zone_a, slow tier (full forms)
@pytest.mark.zone_a
@pytest.mark.slow
@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("stab", ["off", "auto"])
def test_t2p_full(basis, stab):
    ops = _ops()
    r = H.run_t2(ops, basis, "off" if stab == "off" else AUTO)
    for tv in ("0.1", "0.2", "0.5"):
        assert r["rows"][tv]["p_err_max"] <= _band("t2_ref_terzaghi", "p_err_max@" + tv), r
        assert r["rows"][tv]["U_err"] <= _band("t2_ref_terzaghi", "U_err@" + tv), r


@pytest.mark.zone_a
@pytest.mark.slow
@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("stab", ["off", "auto"])
def test_t3_full(basis, stab):
    ops = _ops()
    rec = json.loads((DATA / "t3_mandel_fe.json").read_text())
    r = H.run_t3(ops, basis, "off" if stab == "off" else AUTO, rec)
    assert all(v["err"] <= v["band"] for v in r["rows"].values()), r
    assert r["peak"] >= _band("t3_ref_mandel", "p_center_peak"), r


@pytest.mark.zone_a
@pytest.mark.slow
@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("stab", ["off", "auto"])
def test_t4p_full(basis, stab):
    ops = _ops()
    r = H.run_t4(ops, basis, "off" if stab == "off" else AUTO)
    assert r["p_err"] <= _band("t4_ref_hydro", "hydro_p_err"), r


@pytest.mark.zone_a
@pytest.mark.slow
@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
@pytest.mark.parametrize("stab", ["off", "auto"])
@pytest.mark.parametrize("case", ["terzaghi", "mandel"])
def test_t6_full(basis, stab, case):
    ops = _ops()
    fe = json.loads((DATA / "t6_fe_dt0.005.json").read_text())
    leg = [f for f in fe[case] if (f["stab"] == "off") == (stab == "off")][0]
    st = "off" if stab == "off" else AUTO
    if case == "terzaghi":
        r = H.run_t6_terzaghi(ops, basis, st, leg)
    else:
        r = H.run_t6_mandel(ops, basis, st, leg, _band("t3_ref_mandel", "p_center_peak"))
    assert r["pass"], r


@pytest.mark.zone_a
@pytest.mark.slow
def test_t5_route_e_checkerboard():
    ops = _ops()
    lo = H.band_value(H.case_band_line("t5_ref_b4", "route E, stab auto").split("<=")[1])
    hi = H.band_value(H.case_band_line("t5_ref_b4", "route E, stab off").split(">=")[1])
    assert H.run_t5(ops, AUTO)["CB_lap"] <= lo
    assert H.run_t5(ops, "off")["CB_lap"] >= hi
