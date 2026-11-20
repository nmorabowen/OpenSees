"""WP-237d: the study driver declares the apex face Z of contract amendment 2 (interface 19.2.1).

The worker runs in process on the 54-tet ``mini_pad`` miniature of ``upper3d_fixed_phi30`` from the reference tests
(an Upper whose cones and rows force conditions onto the apex at the wall corners and the footing rim), substituted for
the named case through ``limit_reference.case_problem``. Each test takes a few seconds with Clarabel.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
TESTS = HERE.parents[2] / "tests"
for path in (str(HERE), str(TESTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

clarabel = pytest.importorskip("clarabel")

import test_limit_contract as contract     # noqa: E402
import limit_reference as lr               # noqa: E402
import solver_study as study               # noqa: E402
from test_limit_reference import MINI_PAD  # noqa: E402


@pytest.fixture
def mini_pad_case(monkeypatch):
    case = contract.parse_case_text(MINI_PAD)["mini_pad"]
    real = lr.case_problem
    monkeypatch.setattr(lr, "case_problem",
                        lambda name, mesh=0: contract.build_case_problem(case, mesh) if name == "mini_pad"
                        else real(name, mesh))
    return "mini_pad"


def run_row(tmp_path, case, formulation, overrides=None, tag=""):
    study.worker(case, 0, formulation, "clarabel-qdldl", tmp_path, 300.0, overrides, tag)
    stem = study.stem_of(case, 0, formulation, "clarabel-qdldl", tag)
    return json.loads((tmp_path / f"{stem}.row.json").read_text(encoding="utf-8"))


def stored_result(tmp_path, row):
    result = tmp_path / Path(row["certificate_file"].replace("\\", "/"))
    result = result.with_name(result.name[: -len(".certificate.json")] + ".result.json")
    return contract.load_document(result)


def test_jump_upper_declares_the_face_and_certifies(tmp_path, mini_pad_case):
    row = run_row(tmp_path, mini_pad_case, lr.JUMP_ID)
    assert "error" not in row, row.get("error")
    assert row["apex_zero_enabled"] and row["apex_zero_declared"] and row["result_apex_zero"]
    assert row["n_apex_zero"] == row["n_apex_zero_flow"] + row["n_apex_zero_jump"] > 0
    assert row["n_apex_zero_jump"] > 0
    info = row["apex_zero_info"]
    assert info["n_apex_zero"] == row["n_apex_zero"] and info["stuck_off_apex"] == 0
    assert row["apex_zero_wall_s"] >= 0.0 and "apex_zero" in row["phases"]
    assert row["interior"] == "supplied"
    cert = row["certificate"]
    assert cert["classification"] == "certified_bound", cert
    assert row["rho_apex"] is not None and row["rho_apex"] <= contract.TOLERANCES["TOL_CONE_MARGIN"]
    assert row["method_apex_zero"]["count"] == row["n_apex_zero"]
    assert row["gate_admissible"] is True, row.get("gate_refusal")
    # Z is a result array, so the stored hash covers it (19.2.1 item 2).
    doc, arrays = stored_result(tmp_path, row)
    assert arrays["apex_zero"].shape == (row["n_apex_zero"], 2)
    assert doc["result_sha256"] == row["result_sha256"]
    assert doc["result_sha256"] == contract.document_hash(doc, arrays, "result_sha256")


def test_no_apex_zero_reproduces_the_old_result(tmp_path, mini_pad_case):
    row = run_row(tmp_path, mini_pad_case, lr.JUMP_ID, {"driver_apex_zero": False}, "noz")
    assert "error" not in row, row.get("error")
    assert row["apex_zero_enabled"] is False and row["apex_zero_declared"] is False and row["n_apex_zero"] == 0
    assert "apex_zero_info" not in row and row["result_apex_zero"] is False
    assert row["certificate"]["classification"] == "estimate"
    assert "interior_field_invalid" in row["certificate"]["reasons"]
    assert row["rho_apex"] is None
    _, arrays = stored_result(tmp_path, row)
    assert "apex_zero" not in arrays


def test_finder_failure_is_recorded_and_the_row_goes_on(tmp_path, mini_pad_case, monkeypatch):
    def refuse(program, time_limit=None):
        raise contract.ContractError("no_interior_field", "forced in a test")
    monkeypatch.setattr(lr, "find_apex_zero", refuse)
    row = run_row(tmp_path, mini_pad_case, lr.JUMP_ID)
    assert "error" not in row, row.get("error")
    assert row["apex_zero_declared"] is False and row["apex_zero_info"]["code"] == "no_interior_field"
    assert row["result_apex_zero"] is False
    assert row["certificate"]["classification"] in ("estimate", "certified_bound")


def test_lower_never_declares_a_face(tmp_path, mini_pad_case, monkeypatch):
    def forbidden(program, time_limit=None):
        raise AssertionError("find_apex_zero called on a Lower")
    monkeypatch.setattr(lr, "find_apex_zero", forbidden)
    row = run_row(tmp_path, mini_pad_case, lr.LOWER_ID)
    assert "error" not in row, row.get("error")
    assert "apex_zero_declared" not in row and "result_apex_zero" in row and row["result_apex_zero"] is False


def test_cli_flag_sets_the_driver_key_and_a_tag(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(study, "supervise", lambda *args: seen.setdefault("args", args))
    study.main(["job", "--case", "c", "--mesh", "0", "--formulation", lr.JUMP_ID, "--solver", "clarabel-faer",
                "--data", str(tmp_path), "--no-apex-zero"])
    overrides, tag = seen["args"][7], seen["args"][8]
    assert overrides == {"driver_apex_zero": False} and tag == "noz"
    assert study.PRESETS["noz"] == {"driver_apex_zero": False}
    assert "driver_apex_zero" in study.DRIVER_KEYS
