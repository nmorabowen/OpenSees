"""Pure-Python unit tests for WP-260w0 `summarize.py` sanity & brief §6 gates.

Run:  py -3.12 -m pytest Ladruno_scripts/wp260_profiles/test_summarize_gates.py -v
"""
from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SUMMARIZE = HERE / "summarize.py"
SMOKE_SNAP = HERE / "results" / "smoke_7e91aaf1c"
QUIET_SNAP = HERE / "results" / "smoke_quiet2x3x3_7e91aaf1c"


def _run_summarize(snap_dir: Path, out_md: Path) -> str:
    cp = subprocess.run(
        [sys.executable, str(SUMMARIZE), "--from-json", str(snap_dir), "--out-results", str(out_md)],
        capture_output=True, text=True, check=True,
    )
    assert cp.returncode == 0
    return out_md.read_text(encoding="utf-8")


def _make_mutated_snap(tmp_path: Path, src: Path, mutate_fn) -> Path:
    dst = tmp_path / "snap"
    dst.mkdir(parents=True, exist_ok=True)
    m = json.loads((src / "manifest.json").read_text(encoding="utf-8"))
    t = json.loads((src / "trees.json").read_text(encoding="utf-8"))
    mutate_fn(m, t)
    (dst / "manifest.json").write_text(json.dumps(m), encoding="utf-8")
    (dst / "trees.json").write_text(json.dumps(t), encoding="utf-8")
    return dst


def test_reproduces_committed_smoke_and_quiet_snapshots(tmp_path):
    txt_smoke = _run_summarize(SMOKE_SNAP, tmp_path / "smoke.md")
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_smoke
    assert "`d2_cpml_wave[smoke]`: " in txt_smoke and "one-off `Linear -factorOnce`" in txt_smoke
    assert "upd_constr" in txt_smoke
    assert "| `d4_mpm_modeA[smoke]` |" in txt_smoke and "| NO |" in txt_smoke

    txt_q = _run_summarize(QUIET_SNAP, tmp_path / "quiet.md")
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_q and "| YES |" in txt_q
    assert "| `d3a_explicit_cdl_cpml[smoke]` |" in txt_q and "| YES |" in txt_q
    assert "NOT RANKABLE" not in txt_q


def test_refuses_deep_overhead_below_095(tmp_path):
    def mut(m, _t):
        for r in m["records"]:
            if r["deck"] == "d1b_brick_concrete3d" and r["mode"] == "1t_coarse":
                r["steady_wall_s"] = 10.0  # makes deep_ovh < 0.1

    snap = _make_mutated_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "deep_ovh=" in txt and "<0.95" in txt
    assert "`d1b_brick_concrete3d[smoke]`: deep_ovh=" in txt


def test_refuses_scope_call_count_mismatch(tmp_path):
    def mut(m, t):
        for r in m["records"]:
            if r["deck"] == "d1b_brick_concrete3d" and r["mode"] == "nt_coarse":
                step = t[r["run_id"]]["tree"]["children"]["step"]
                step["children"]["solveCurrentStep"]["children"]["update"]["calls"] += 1
                break

    snap = _make_mutated_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "call-count mismatch" in txt


def test_refuses_residual_above_05_percent(tmp_path):
    def mut(m, t):
        for r in m["records"]:
            if r["deck"] == "d3a_explicit_cdl_cpml" and r["mode"] == "1t_deep":
                step = t[r["run_id"]]["tree"]["children"]["step"]
                step["wall_ns"] = int(step["wall_ns"] * 1.03)  # injects +3% residual

    snap = _make_mutated_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "|residual|=" in txt and ">0.5%" in txt


def test_refuses_unresolved_remainder_above_5_percent(tmp_path):
    def mut(m, t):
        for r in m["records"]:
            if r["deck"] == "d3a_explicit_cdl_cpml":
                r["constraints"] = "Plain"  # makes ~23% enforceSPs time unresolved!

    snap = _make_mutated_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "unresolved newStep_other=" in txt and ">5%" in txt
