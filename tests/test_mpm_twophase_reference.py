"""WP-232l (232l-V): checks the REFERENCES of the two-phase gates -- never MPM.

Round 1: the closed forms (Terzaghi copy against its source, the ramp Duhamel, Mandel's series
against the independent Laplace solution) and the provenance and shape of the stored LadrunoUP
records (tests/data/mpm/twophase/*.json). numpy + scipy only; no fork binary needed. The FE legs
are re-run by `python3.12 tests/mpm_twophase_reference.py <leg>`.
"""
from __future__ import annotations

import ast
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpm_twophase_reference as R  # noqa: E402

pytestmark = pytest.mark.zone_a


def _json(name):
    p = R.DATA / f"{name}.json"
    if not p.exists():
        pytest.skip(f"{p.name} not recorded yet")
    return json.loads(p.read_text(encoding="utf-8"))


def _provenance_ok(d):
    assert len(d["build"]) == 40 and all(c in "0123456789abcdef" for c in d["build"])
    assert d["deck"].startswith("tests/mpm_twophase_reference.py::")
    assert d["units"].startswith("kN, t, m, s, kPa")


def _source_functions(path, names):
    """Exec the named top-level functions of a test module without importing it (its import
    binds an engine, test_ladruno_up_element_analytic.py:74-79)."""
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    ns = {"np": np}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            exec(compile(ast.Module([node], []), str(path), "exec"), ns)
    return ns


def test_terzaghi_copy_matches_its_source():
    ns = _source_functions(HERE / "test_ladruno_up_element_analytic.py",
                           {"terzaghi_series", "terzaghi_U"})
    for Z in (0.0, 0.3, 0.9):
        for Tv in (0.01, 0.1, 0.5):
            assert ns["terzaghi_series"](Z, Tv) == R.terzaghi_series(Z, Tv)
    for Tv in (0.01, 0.1, 0.5):
        assert ns["terzaghi_U"](Tv) == R.terzaghi_U(Tv)


def test_terzaghi_case_is_the_fe_tests():
    c = R.TERZ
    assert (c["E"], c["nu"], c["rho"], c["Kf"], c["poro"], c["rhoF"], c["kbar"], c["L"],
            c["q"], c["nely"], c["dt"], c["gamma"], c["beta"]) == (
        1.0e4, 0.2, 2.0, 4444.0, 0.4, 1.0, 1.0e-4, 10.0, 10.0, 10, 0.02, 0.6, 0.3025)


def test_terzaghi_ramp_limits():
    # ramp -> step; long after the ramp the ramp response lags the step by ~ Tvr/2
    assert abs(R.terzaghi_ramp(0.4, 0.2, 1e-6) - R.terzaghi_series(0.4, 0.2)) < 1e-5
    a = R.terzaghi_ramp(0.5, 0.3, 0.02)
    b = R.terzaghi_series(0.5, 0.3 - 0.01)
    assert abs(a - b) < 2e-4
    k = R.terzaghi_constants(1e4, 0.2, 4444.0, 0.4, 1e-4, 10.0, 10.0)
    assert math.isclose(R.terzaghi_settle_ramp(50.0, 0.02, k["s0"], k["sinf"]), k["sinf"],
                        rel_tol=1e-12)


def test_mandel_closed_form_selfcheck():
    out = R.closed_selfcheck()
    assert out["mandel_c_vs_cd_form"] < 1e-14
    assert out["mandel_series_vs_laplace_max_over_p0"] < 1e-8
    assert out["mandel_p_0plus_uniform_err"] < 1e-3


def test_mandel_peak_meets_the_kickoff_floor():
    mc = R.mandel_constants(**R.MANDEL_SOIL)
    pk, tv = R.mandel_peak(mc, nt=1500)
    assert pk >= 1.05 and 0.05 < tv < 0.09


def test_cases_file_reference_blocks_point_at_records():
    text = (HERE / "mpm_twophase_cases.txt").read_text(encoding="utf-8")
    for line in text.splitlines():
        if line.startswith("data ") and "<dt>" not in line:
            assert (R.REPO / line.split()[1]).exists(), line


def test_records_t2_terzaghi():
    d = _json("t2_terzaghi_fe")
    _provenance_ok(d)
    for run in d["runs"]:
        assert [r["Tv"] for r in run["rows"]] == list(R.TV_T2)
        for r in run["rows"][1:]:
            assert r["p_err_gimp_points"] < 0.02 and r["U_err"] < 0.02


def test_records_t3_mandel():
    d = _json("t3_mandel_fe")
    _provenance_ok(d)
    assert "UNVERIFIED" in d["closed_form_note"]
    assert d["selfcheck"]["mandel_series_vs_laplace_max_over_p0"] < 1e-8
    pk = d["closed_form_peak_ramp"]["p_over_p0"]
    assert pk >= 1.05 and math.isclose(d["peak_gate"], 1 + 0.5 * (pk - 1), rel_tol=1e-15)
    assert len(d["gated"]) == 2
    for run in d["gated"]:
        assert run["t_ramp"] == R.MANDEL["t_ramp"] and run["dt"] == 0.02
        assert [r["tv"] for r in run["rows"]] == list(R.MANDEL_TVS)
        assert all(b == 0.02 for b in run["band"])          # FE error < 0.01 everywhere


def test_records_t4_t5():
    h = _json("t4_hydro_fe")
    _provenance_ok(h)
    for run in h["gated"] + h["reported"]:
        assert run["p_err_over_gwH"] < 1e-6
        # sigma' is element-constant: the point error is the h/4 (GIMP) and h/(2 sqrt 3) offsets
        hL = 1.0 / run["nely"]
        assert math.isclose(run["sig_err_gimp_over_rhopgH"], hL / 4, rel_tol=1e-5)
        assert math.isclose(run["sig_err_gauss_over_rhopgH"], hL / (2 * math.sqrt(3)), rel_tol=1e-5)
    b = _json("t5_b4_fe")
    _provenance_ok(b)
    cb = {json.dumps(r["stab"]): r["CB_lap"] for r in b["runs"]}
    assert cb['"off"'] >= 0.5 and cb['["auto", 0.25]'] <= 0.3
