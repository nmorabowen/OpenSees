"""WP-232f (232f-V): checks the REFERENCES of tests/mpm_implicit_cases.txt -- never MPM.

* the cases file parses, covers every P1 gate row, carries the pre-registered tolerances
  (.claude/skills/ladruno-mpm/references/wp-briefs.md, "232f") and no TBD;
* the closed forms reproduce (Prandtl-Reissner N_q, Howland, the bar period);
* the elastica ODE is converged;
* the stored fork FE runs (tests/data/mpm/implicit/*.json) carry their provenance, converged in
  their own meshes, and the cases file's expected values are exactly what they give;
* the plate seeding rule reproduces the case's checksums and writes a valid Contract-1 file.

numpy only, except the Contract-1 write (h5py; skipped without it). No fork binary needed: the FE
legs are re-run by `python3.12 tests/mpm_implicit_reference.py fe_<leg>`.
"""
from __future__ import annotations

import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import mpm_implicit_reference as R  # noqa: E402

pytestmark = pytest.mark.zone_a

CASES = R.parse_cases()
DATA = R.DATA


def _json(name):
    return json.loads((DATA / f"{name}.json").read_text(encoding="utf-8"))


def _expects(case, obs):
    return [t for k, t in CASES[case]["lines"] if k == "expect" and t[0] == obs]


def _tol(tokens):
    i = tokens.index("tol")
    return float(tokens[i + 1]), float(tokens[i + 2])


def _eq_value(tokens):
    return [float(x) for x in tokens[tokens.index("eq") + 1:tokens.index("tol")]]


# ----------------------------------------------------------------------------------------------
# The cases file
# ----------------------------------------------------------------------------------------------


def test_no_tbd_remains():
    text = R.CASES.read_text(encoding="utf-8")
    body = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    assert not [ln for ln in body if re.search(r"\bTBD\b", ln)]


def test_every_run_case_states_its_reference_and_tolerance():
    for name, c in CASES.items():
        kind = c["kv"]["kind"][0][0]
        assert kind in ("run", "protocol", "refusal"), name
        if kind == "run":
            assert c["kv"].get("ref"), f"{name}: no 'ref'"
            assert c["kv"]["mode"][0][0] in ("transient", "quasistatic"), name
            assert c["kv"].get("expect"), f"{name}: no expectation"
            for t in c["kv"]["expect"]:
                assert "tol" in t or "le" in t or "ge" in t, f"{name}: {t}"
        if kind == "refusal":
            assert ["refused"] in c["kv"]["expect"], name


GATE_ROWS = {
    # gate-table row -> (case, observable, rel tolerance or bound)
    "cantilever_qs": ("cantilever_qs_r3", "tip_disp_y", 0.05),
    "cantilever_dyn_peak": ("cantilever_dyn_r3", "first_peak_tip_disp_y", 0.05),
    "cantilever_dyn_period": ("cantilever_dyn_r3", "first_period_tip", 0.05),
    "footing": ("footing_mildcone_r3", "cond_pressure_max", 0.10),
    "plate": ("plate_hole_dw050", "peak_sthth_over_snet", 0.05),
    "bar": ("bar_implicit_2d", "period_comx", 0.01),
}


@pytest.mark.parametrize("row", sorted(GATE_ROWS))
def test_gate_tolerances_are_the_pre_registered_ones(row):
    case, obs, rel = GATE_ROWS[row]
    assert CASES[case]["kv"].get("status", [["gate"]])[0][0] == "gate"
    (tok,) = [t for t in _expects(case, obs) if t[-1] != "report"]
    assert _tol(tok) == (0.0, rel)


def test_footing_gate_is_the_exact_prandtl_band():
    q = R.footing_closed()["q_exact"]
    for r in ("r1", "r2", "r3"):
        (tok,) = _expects(f"footing_mildcone_{r}", "cond_pressure_max")
        assert _eq_value(tok) == [pytest.approx(q, abs=5e-4)]
        assert _tol(tok) == (0.0, 0.10)


def test_column_is_the_slice1_small_case_leg():
    """Geometry, material and seeding equal tests/mpm_cases.txt column_selfweight_small_2d_gimp;
    the frozen bound 0.03 is kept (Amendment 4); the 20 % case stays a strict xfail."""
    s1kv = dict(_slice1_lines("column_selfweight_small_2d_gimp"))
    me = CASES["column_selfweight_small_qs"]["kv"]
    for k in ("ndm", "basis", "origin", "h", "cells", "material", "body", "gravity"):
        assert me[k][0] == s1kv[k], k
    assert _expects("column_selfweight_small_qs", "selfweight_error")[0][-2:] == ["le", "0.03"]
    assert CASES["column_selfweight_qs"]["kv"]["status"][0][0] == "xfail"


def _slice1_lines(name):
    out, on = [], False
    for raw in (HERE / "mpm_cases.txt").read_text(encoding="utf-8").splitlines():
        t = raw.split()
        if t[:2] == ["case", name]:
            on = True
            continue
        if on and t[:1] == ["end"]:
            break
        if on and t and not t[0].startswith("#"):
            out.append((t[0], t[1:]))
    return out


def test_protocol_tolerances():
    assert _tol([t for k, t in CASES["retry_bar"]["lines"] if k == "compare"][0]) == (0.0, 1e-12)
    assert _tol([t for k, t in CASES["reseed_bar_forced"]["lines"] if k == "compare"][0]) == \
        (0.0, 1e-12)
    for c in ("restart_bar", "restart_footing"):
        assert _tol([t for k, t in CASES[c]["lines"] if k == "compare"][0]) == (0.0, 1e-10)
        assert CASES[c]["kv"]["engines"][0] == ["tcl", "python"]
    # the retry reference ends at the same time as the retried run
    sched = [list(map(float, t)) for t in CASES["retry_bar_reference"]["kv"]["schedule"]]
    assert sum(d * n for d, n in sched) == pytest.approx(0.02 * 200, abs=1e-12)


def test_refusals_cover_the_gate_row():
    names = {n for n, c in CASES.items() if c["kv"]["kind"][0][0] == "refusal"}
    assert {"refuse_integrator_centraldifference", "refuse_ladrunohht",
            "refuse_ladrunogeneralizedalpha", "refuse_variabletransient", "refuse_arclength",
            "refuse_dispcontrol_gridnode"} <= names


# ----------------------------------------------------------------------------------------------
# Closed forms
# ----------------------------------------------------------------------------------------------


def test_prandtl_reissner_mild_cone():
    f = R.footing_closed()
    assert f["phi_ps"] == pytest.approx(27.47, abs=5e-3)
    assert f["Nq"] == pytest.approx(13.89, abs=5e-3)
    assert f["q_exact"] == pytest.approx(138.9, abs=0.01)
    # Prandtl's N_q at phi = 30 deg is 18.40 (textbook value): the formula itself
    assert R.prandtl_nq(30.0) == pytest.approx(18.401, abs=1e-3)
    # the apex regularizer's cohesion is a stated, small model-form term
    assert f["cNc_fraction"] < 0.03


def test_howland_and_heywood_agree_at_dw_050():
    assert R.howland_ktn(0.5) == pytest.approx(2.16)
    assert R.howland_ktn(0.0) == pytest.approx(3.0)            # Kirsch limit
    assert abs(R.heywood_ktn(0.5) / R.howland_ktn(0.5) - 1) < 0.02


def test_bar_period_and_its_discretization_gap():
    assert R.bar_closed()["period"] == 10.0
    for lumped in (True, False):
        d = R.bar_discrete(lumped=lumped)
        assert abs(d["period_newmark"] / 10.0 - 1) < 2e-4      # far inside the 1 % gate


def test_cantilever_linear_parameters():
    assert R.cantilever_linear(6e7)["qL3_EI"] == pytest.approx(2.0)
    assert R.cantilever_linear(1.2e8)["qL3_EI"] == pytest.approx(1.0)
    assert R.cantilever_linear(1.2e8)["period1_linear"] == pytest.approx(1.787, abs=1e-3)


# ----------------------------------------------------------------------------------------------
# The elastica
# ----------------------------------------------------------------------------------------------


def test_elastica_converged_and_linear_limit():
    a, b = R.elastica_tip(2.0, 500), R.elastica_tip(2.0, 2000)
    assert abs(a["tip_down_over_L"] - b["tip_down_over_L"]) < 1e-9
    assert R.elastica_tip(1e-3, 500)["tip_down_over_L"] == pytest.approx(1e-3 / 8, rel=1e-4)
    assert b["tip_down_over_L"] == pytest.approx(0.238507, abs=1e-6)


# ----------------------------------------------------------------------------------------------
# Stored FE legs: provenance, own-mesh convergence, and the cases file's values
# ----------------------------------------------------------------------------------------------
LEGS = ("fe_cantilever_qs", "fe_cantilever_dyn", "fe_footing", "fe_plate")


@pytest.mark.parametrize("leg", LEGS)
def test_fe_leg_provenance(leg):
    d = _json(leg)
    assert re.fullmatch(r"[0-9a-f]{40}", d["build"]), d["build"]
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", d["date"])
    assert d["deck"].startswith("tests/mpm_implicit_reference.py::")
    assert d["output"] == f"tests/data/mpm/implicit/{leg}.json"


def test_cantilever_qs_fe_converged_and_is_the_case_value():
    d = _json("fe_cantilever_qs")
    seq = [x["tip"]["r3"] for x in d["hypo"]]
    diffs = np.abs(np.diff(seq))
    assert np.all(diffs[1:] < diffs[:-1])                     # monotone convergence
    assert 1.7 < d["order_hypo"] < 2.3                         # shear locking: O(h^2)
    ref = d["reference"]
    assert abs(ref["r3"] - seq[-1]) / abs(ref["r3"]) < 1e-3   # extrapolation is small
    # an independent kinematic family (EICR corotational) agrees to 1 %
    assert abs(d["corot"][-1]["tip"]["r3"] / ref["r3"] - 1) < 0.01
    # and the elastica (inextensible, shear-rigid beam) to 1 %
    assert abs(d["elastica_tip"] / d["hypo"][-1]["tip_at_L"] - 1) < 0.01
    for r, x in (("r1", 9.875), ("r2", 9.9375), ("r3", 9.96875)):
        (tok,) = _expects(f"cantilever_qs_{r}", "tip_disp_y")
        assert float(tok[1]) == x
        assert _eq_value(tok) == [pytest.approx(ref[r], rel=1e-9)]


def test_cantilever_dyn_fe_converged_and_is_the_case_value():
    d = _json("fe_cantilever_dyn")
    cons = d["consistent"]
    for key in ("first_peak", "first_period"):
        seq = [x["tip"]["r3"][key] for x in cons]
        diffs = np.abs(np.diff(seq))
        assert np.all(diffs[1:] < diffs[:-1]), key
        assert abs(d["reference"]["r3"][key] / seq[-1] - 1) < 3e-3, key
    # mass lumping and halving dt move the answer far less than the 5 % gate
    c16 = [x for x in cons if x["ny"] == 16][0]["tip"]["r3"]
    for alt in (d["lumped"][0]["tip"]["r3"], d["dt_half"][0]["tip"]["r3"]):
        for key in ("first_peak", "first_period"):
            assert abs(alt[key] / c16[key] - 1) < 5e-3
    for r in ("r1", "r2", "r3"):
        (pk,) = _expects(f"cantilever_dyn_{r}", "first_peak_tip_disp_y")
        (pt,) = _expects(f"cantilever_dyn_{r}", "first_period_tip")
        assert _eq_value(pk) == [pytest.approx(d["reference"][r]["first_peak"], rel=1e-5)]
        assert _eq_value(pt) == [pytest.approx(d["reference"][r]["first_period"], rel=1e-5)]


def test_footing_fe_on_the_case_domain():
    """Context, not the gate: the FE control (B-bar) on the MPM case's own domain, staging and
    rough footing reproduces ADR-79's monotone, non-converging sequence (non-associated flow past
    the localization threshold): in the band at 4 and 8 across B, BELOW it (0.882) at 16 -- a
    recorded finding the cases file states. Roughness moves the coarse answer ~3 %."""
    d = _json("fe_footing")
    q = d["closed"]["q_exact"]
    rough = sorted(d["rough"], key=lambda x: x["n_across_b"])
    assert [x["n_across_b"] for x in rough] == [4, 8, 16]
    for x in rough:
        assert x["s_end"] == pytest.approx(0.2, abs=1e-9)
    r = [x["q_max"] / q for x in rough]
    assert all(b < a for a, b in zip(r, r[1:]))
    assert 0.90 <= r[0] <= 1.10 and 0.90 <= r[1] <= 1.10
    assert r[2] == pytest.approx(0.8821, abs=5e-4)
    assert abs(d["smooth_4"]["q_max"] / rough[0]["q_max"] - 1) < 0.04


def test_footing_locking_probe_is_recorded():
    """The same H8 without B-bar (full integration, as 2 x 2 particles per cell integrate):
    isochoric plastic flow (psi = 0) locks, the over-prediction falls with refinement but starts
    far outside the band. Recorded so a mode-A miss on this row is diagnosable (no anti-locking
    until ADR-232 §9 phase 4)."""
    d = _json("fe_footing")
    q = d["closed"]["q_exact"]
    std = sorted(d["std_no_antilocking"], key=lambda x: x["n_across_b"])
    qs = [x["q_max"] / q for x in std]
    assert qs[0] > 1.10
    assert all(b < a for a, b in zip(qs, qs[1:]))


def test_plate_fe_converged_and_within_one_percent_of_howland():
    d = _json("fe_plate")
    k = [x["Ktn"] for x in d["levels"]]
    assert abs(k[-1] / k[-2] - 1) < 1e-3
    assert abs(k[-1] / R.howland_ktn(0.5) - 1) < 0.01
    # the finite length and the displacement end condition are worth < 0.2 %
    assert abs(d["long_traction"]["Ktn"] / k[-2] - 1) < 2e-3
    (tok,) = _expects("plate_hole_dw050", "peak_sthth_over_snet")
    assert _eq_value(tok) == [R.howland_ktn(0.5)]


# ----------------------------------------------------------------------------------------------
# The plate's particle set
# ----------------------------------------------------------------------------------------------


def test_plate_seeding_checksums_are_the_case_values():
    c = R.plate_seed_checksums()
    assert c["volume_sum"] == pytest.approx(c["exact_area"], abs=1e-12)
    assert c["max_lp_over_h"] <= 0.5
    (n,) = _expects("plate_hole_dw050", "seed_count")
    (v,) = _expects("plate_hole_dw050", "seed_volume_sum")
    (m,) = _expects("plate_hole_dw050", "seed_moment")
    assert _eq_value(n) == [c["count"]]
    assert _eq_value(v) == [pytest.approx(c["volume_sum"], abs=1e-13)]
    assert _eq_value(m) == [pytest.approx(c["moment"][0], abs=1e-12),
                            pytest.approx(c["moment"][1], abs=1e-12)]


def test_plate_particle_set_is_a_valid_contract1_file(tmp_path):
    pytest.importorskip("h5py")
    import mpm_contracts as MC  # noqa: E402
    p = tmp_path / "plate_hole_dw050_h32.h5"
    R.write_plate_contract1(p)
    warns: list = []
    assert MC.validate_particle_set(str(p), warns) == []


def test_holeplate_cells_tile_the_region_at_any_resolution():
    for nt, nr in ((8, 4), (32, 15), (64, 30)):
        pos, vol, lp = R.holeplate_particles(1.0, 2.0, nt, nr)
        assert vol.sum() == pytest.approx(4.0 - math.pi / 4, abs=1e-12)
        assert np.all(np.hypot(pos[:, 0], pos[:, 1]) > 1.0)
        assert np.all(pos >= lp - 1e-15) and np.all(pos + lp <= 2.0 + 1e-15)
