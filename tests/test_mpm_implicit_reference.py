"""WP-232f (232f-V): checks the REFERENCES of tests/mpm_implicit_cases.txt -- never MPM.

* the cases file parses, covers every P1 gate row, carries the pre-registered tolerances
  (.claude/skills/ladruno-mpm/references/wp-briefs.md, "232f") and no TBD;
* the closed forms reproduce (Prandtl-Reissner N_q, Howland, the bar period);
* the elastica ODE is converged;
* the stored fork FE runs (tests/data/mpm/implicit/*.json) carry their provenance, converged in
  their own meshes, and the cases file's expected values are exactly what they give;
* the plate seeding rule reproduces the case's checksums and writes a valid Contract-1 file;
* WP-232j (232j-V), gate L3: the (2 + pi) k closed form, the FE record fe_footing_j2.json's
  shape and admission block, and the L3 block of tests/mpm_antilocking_cases.txt.

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


# ----------------------------------------------------------------------------------------------
# WP-232j (232j-V), gate L3: the von Mises strip footing -- closed form and the FE record.
# Never MPM: the MPM side of L3 is 232j-C's test.
# ----------------------------------------------------------------------------------------------
AL_CASES_PATH = HERE / "mpm_antilocking_cases.txt"
J2_RUNGS = {"r1": 4, "r2": 8, "r3": 16}


def test_footing_j2_closed_form():
    c = R.footing_j2_closed()
    assert c["k"] == pytest.approx(100.0, rel=1e-14)               # k = sigma_y / sqrt(3)
    assert c["q_u"] == pytest.approx((2.0 + math.pi) * 100.0, rel=1e-14)
    assert round(c["q_u"], 2) == 514.16
    assert c["band_factor"] == [0.97, 1.08]                        # wp-briefs.md "232j" row L3
    assert c["band"] == pytest.approx([498.7345, 555.2920], abs=1e-4)
    K, G = R.j2_moduli()
    assert K == pytest.approx(1.0e5 / 1.2) and G == pytest.approx(1.0e5 / 2.6)
    assert 9 * K * G / (3 * K + G) == pytest.approx(1.0e5)           # E back
    assert (3 * K - 2 * G) / (2 * (3 * K + G)) == pytest.approx(0.3)  # nu back


def test_q_lim_is_the_window_mean():
    s = np.arange(1, 51) * 0.002
    q = 1000.0 * s
    out = R.q_lim([[a, b] for a, b in zip(s, q)])
    assert out["n"] == 11
    assert out["mean"] == pytest.approx(40.0) and out["trapz"] == pytest.approx(40.0)
    # a run that stops before the window closes has no q_lim
    assert R.q_lim([[a, b] for a, b in zip(s[:20], q[:20])])["mean"] is None


def _j2_record():
    d = _json("fe_footing_j2")
    return d, d["legs"]


def test_footing_j2_record_provenance_and_setup():
    d, legs = _j2_record()
    assert re.fullmatch(r"[0-9a-f]{40}", d["build"]), d["build"]
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", d["date"])
    assert d["deck"] == "tests/mpm_implicit_reference.py::leg_footing_j2"
    assert d["output"] == "tests/data/mpm/implicit/fe_footing_j2.json"
    st = d["setup"]
    assert (st["q0"], st["weight"], st["rough"]) == (0.0, 0.0, True)
    assert (st["ds"], st["s_max"], st["s_lim"]) == (0.002, 0.1, [0.03, 0.05])
    assert st["ladder"] == [["KrylovNewton", 1e-6, 40], ["NewtonLineSearch", 1e-6, 60],
                            ["KrylovNewton", 1e-5, 80]]
    m = d["material"]
    assert m["command"][0] == "LadrunoJ2" and m["command"][4:6] == ["-iso", "voce"]
    assert m["command"][6] == pytest.approx(100.0 * math.sqrt(3.0))
    assert m["command"][7:] == [0.0, 0.0, 0.0]                       # perfectly plastic
    assert d["closed"]["q_u"] == pytest.approx(R.footing_j2_closed()["q_u"], rel=1e-14)
    # r1 is always present in both forms; every leg is from the record's build
    for form in ("bbar", "std"):
        assert "4" in legs[form], form
    for form, by_n in legs.items():
        assert form in ("bbar", "std")
        for n, r in by_n.items():
            assert r["build"] == d["build"], (form, n)
            assert r["material"] == "j2" and r["form"] == form and r["n_across_b"] == int(n)


def test_footing_j2_legs_shape_and_q_lim():
    d, legs = _j2_record()
    qu = d["closed"]["q_u"]
    for form, by_n in legs.items():
        for n, r in by_n.items():
            c = np.array(r["curve"])
            assert np.all(np.diff(c[:, 0]) > 0), (form, n)           # s increases
            assert r["reached_s_max"] and r["s_end"] == pytest.approx(0.1, abs=1e-9)
            assert r["steps"] == len(c)
            assert isinstance(r["failed_attempts"], int) and isinstance(r["relaxed"], int)
            assert r["wall_s"] >= 0
            ql = R.q_lim(r["curve"])
            assert r["q_lim"] == pytest.approx(ql["mean"], rel=1e-12)
            assert r["q_lim_n"] >= 11
            assert r["q_max"] == pytest.approx(c[:, 1].max(), rel=1e-12)
            assert r["q_lim_over_exact"] == pytest.approx(r["q_lim"] / qu, rel=1e-12)
    # full integration locks under isochoric flow: std above B-bar at every rung run
    for n in legs["bbar"]:
        if n in legs["std"]:
            assert legs["std"][n]["q_lim"] > legs["bbar"][n]["q_lim"], n


def test_footing_j2_admission_block_reads_the_band():
    """The admission verdict is recomputed from the legs; the band is the pre-registered one."""
    d, legs = _j2_record()
    a = d["admission"]
    assert a == R.footing_j2_admission(legs)
    assert a["band"] == pytest.approx(R.footing_j2_closed()["band"])
    if "16" in legs.get("bbar", {}):
        q = legs["bbar"]["16"]["q_lim"]
        assert a["admitted"] == (a["band"][0] <= q <= a["band"][1])


def test_antilocking_cases_l3_block():
    cases = R.parse_cases(AL_CASES_PATH)
    c = R.footing_j2_closed()
    K, G = R.j2_moduli()
    _, legs = _j2_record()
    for r, n in J2_RUNGS.items():
        name = f"footing_j2_{r}"
        kv = cases[name]["kv"]
        mild = CASES[f"footing_mildcone_{r}"]["kv"]
        # same domain, grid, seeding, footing and staging as footing_mildcone
        for key in ("ndm", "basis", "origin", "h", "cells", "body", "fix", "gridvel", "series",
                    "dt"):
            assert kv[key] == mild[key], (name, key)
        assert kv["mode"] == [["quasistatic"]]
        assert "pforce" not in kv and "gravity" not in kv         # weightless, no surcharge
        assert float(kv["h"][0][0]) == pytest.approx(2.0 / n)
        (mat,) = kv["material"]
        assert mat[1] == "ladrunoj2"
        assert [float(x) for x in mat[2:]] == pytest.approx([K, G, 100.0 * math.sqrt(3.0)])
        # 10 empty staging steps + 50 pushes of 2 mm = s 0.1 m; the window t 25..35 = s .03..05
        assert int(kv["steps"][0][0]) == 60
        assert kv["ladder"] == [["KrylovNewton", "1e-6", "40", "NewtonLineSearch", "1e-6", "60",
                                 "KrylovNewton", "1e-5", "80", "relaxed", "3", "halve"]]
        win = [t for t in kv["expect"] if t[0] == "cond_pressure_window"]
        assert sorted(t[6] for t in win) == ["ge", "le"]
        for t in win:
            assert t[1:6] == ["footing", "y", "1.0", "25.0", "35.0"]
            assert 0.002 * (float(t[4]) - 10.0) == pytest.approx(0.03)
            assert 0.002 * (float(t[5]) - 10.0) == pytest.approx(0.05)
            band = c["band"][0] if t[6] == "ge" else c["band"][1]
            assert float(t[7]) == pytest.approx(band, abs=1e-3)
        refs = " ".join(" ".join(t) for t in kv["ref"])
        assert "514.159" in refs
        if str(n) in legs.get("bbar", {}):
            assert f"{legs['bbar'][str(n)]['q_lim']:.2f}" in refs, name
            assert "PENDING" not in refs, name
        else:
            assert "PENDING" in refs, name
    # the r3 refs state the record's admission verdict
    a = _json("fe_footing_j2")["admission"]
    if "admitted" in a:
        r3 = " ".join(" ".join(t) for t in cases["footing_j2_r3"]["kv"]["ref"])
        assert ("ADMITTED" in r3) == a["admitted"]
        assert ("non-discriminating at r3" in r3) == bool(a.get("non_discriminating_at_r3"))
    # the gate is at r3; r1 and r2 report
    assert cases["footing_j2_r1"]["kv"]["status"] == [["report"]]
    assert cases["footing_j2_r2"]["kv"]["status"] == [["report"]]
    assert "status" not in cases["footing_j2_r3"]["kv"]
