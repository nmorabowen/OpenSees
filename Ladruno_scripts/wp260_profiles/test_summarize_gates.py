"""Pure-Python unit tests for WP-260w0 `summarize.py` sanity & brief §6 gates.

Run:  py -3.12 -m pytest Ladruno_scripts/wp260_profiles/test_summarize_gates.py -v
"""
from __future__ import annotations

import copy
import json
import subprocess
import sys
from pathlib import Path

from typing import List, Optional

HERE = Path(__file__).resolve().parent
SUMMARIZE = HERE / "summarize.py"
SMOKE_SNAP = HERE / "results" / "smoke_7e91aaf1c"
QUIET_SNAP = HERE / "results" / "smoke_quiet2x3x3_7e91aaf1c"


def _run_summarize(snap_dir: Path, out_md: Path, extra_args: Optional[List[str]] = None) -> str:
    cmd = [sys.executable, str(SUMMARIZE), "--from-json", str(snap_dir), "--out-results", str(out_md)]
    if extra_args:
        cmd.extend(extra_args)
    cp = subprocess.run(cmd, capture_output=True, text=True, check=True)
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


def _make_clean_snap(tmp_path: Path, src: Path, mutate_fn=None) -> Path:
    def wrap(m, t):
        for r in m["records"]:
            r["cpu_load_before_pct"] = 2.0
            r["cpu_load_after_pct"] = 2.0
            r["cpu_load_during_ext_pct"] = 1.0
        if mutate_fn is not None:
            mutate_fn(m, t)

    return _make_mutated_snap(tmp_path, src, wrap)


def test_committed_noisy_snapshots_refused_at_default_max_load(tmp_path):
    txt_smoke = _run_summarize(SMOKE_SNAP, tmp_path / "smoke_def.md")
    assert "Host CPU load gate (`--max-load`):" in txt_smoke and "`<= 10.0%`" in txt_smoke
    assert "`max_load`" in txt_smoke
    assert "cpu_load_before=" in txt_smoke and ">10.0%" in txt_smoke
    assert "| NO |" in txt_smoke and "| YES |" not in txt_smoke

    txt_q = _run_summarize(QUIET_SNAP, tmp_path / "quiet_def.md")
    assert "`<= 10.0%`" in txt_q
    assert "cpu_load_before=" in txt_q and ">10.0%" in txt_q
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_q and "| NO |" in txt_q
    assert "| `d3a_explicit_cdl_cpml[smoke]` |" in txt_q and "| YES |" not in txt_q


def test_reproduces_committed_smoke_and_quiet_snapshots_accounting(tmp_path):
    txt_smoke = _run_summarize(SMOKE_SNAP, tmp_path / "smoke.md", ["--max-load", "100"])
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_smoke
    assert "`d2_cpml_wave[smoke]`: " in txt_smoke and "one-off `Linear -factorOnce`" in txt_smoke
    assert "upd_constr" in txt_smoke
    assert "| `d4_mpm_modeA[smoke]` |" in txt_smoke and "| NO |" in txt_smoke

    txt_q = _run_summarize(QUIET_SNAP, tmp_path / "quiet.md", ["--max-load", "100"])
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_q and "| YES |" in txt_q
    assert "| `d3a_explicit_cdl_cpml[smoke]` |" in txt_q and "| YES |" in txt_q
    assert "NOT RANKABLE" not in txt_q


def test_clean_record_passes_and_high_load_before_refused(tmp_path):
    clean = _make_clean_snap(tmp_path / "clean", QUIET_SNAP)
    txt_ok = _run_summarize(clean, tmp_path / "ok.md")
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_ok and "| YES |" in txt_ok
    assert "| `d3a_explicit_cdl_cpml[smoke]` |" in txt_ok and "| NO |" not in txt_ok
    assert "NOT RANKABLE" not in txt_ok

    def mut_hi(m, _t):
        for r in m["records"]:
            if r["deck"] == "d1b_brick_concrete3d" and r["mode"] == "1t_coarse":
                r["cpu_load_before_pct"] = 12.5
                break

    bad = _make_clean_snap(tmp_path / "bad_before", QUIET_SNAP, mut_hi)
    txt_bad = _run_summarize(bad, tmp_path / "bad_before.md")
    assert "cpu_load_before=12.5%>10.0%" in txt_bad
    assert "| `d1b_brick_concrete3d[smoke]` |" in txt_bad and "12.5%" in txt_bad


def test_refuses_high_load_after(tmp_path):
    def mut(m, _t):
        for r in m["records"]:
            if r["deck"] == "d3a_explicit_cdl_cpml" and r["mode"] == "nt_coarse":
                r["cpu_load_after_pct"] = 14.2
                break

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "cpu_load_after=14.2%>10.0%" in txt
    assert "| `d3a_explicit_cdl_cpml[smoke]` |" in txt and "14.2%" in txt


def test_refuses_high_load_during_ext(tmp_path):
    def mut(m, _t):
        for r in m["records"]:
            if r["deck"] == "d1b_brick_concrete3d" and r["mode"] == "1t_deep":
                r["cpu_load_during_ext_pct"] = 11.4
                break

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "cpu_load_during_ext=11.4%>10.0%" in txt


def test_custom_max_load_flag(tmp_path):
    def mut(m, _t):
        for r in m["records"]:
            if r["deck"] == "d1b_brick_concrete3d" and r["mode"] == "1t_coarse":
                r["cpu_load_before_pct"] = 6.5
                break

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
    txt_pass = _run_summarize(snap, tmp_path / "pass.md")
    assert "NOT RANKABLE" not in txt_pass

    txt_fail = _run_summarize(snap, tmp_path / "fail.md", ["--max-load", "5.0"])
    assert "`<= 5.0%`" in txt_fail
    assert "cpu_load_before=6.5%>5.0%" in txt_fail


def test_refuses_deep_overhead_below_095(tmp_path):
    def mut(m, _t):
        for r in m["records"]:
            if r["deck"] == "d1b_brick_concrete3d" and r["mode"] == "1t_coarse":
                r["steady_wall_s"] = 10.0  # makes deep_ovh < 0.1

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
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

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "call-count mismatch" in txt


def test_refuses_residual_above_05_percent(tmp_path):
    def mut(m, t):
        for r in m["records"]:
            if r["deck"] == "d3a_explicit_cdl_cpml" and r["mode"] == "1t_deep":
                step = t[r["run_id"]]["tree"]["children"]["step"]
                step["wall_ns"] = int(step["wall_ns"] * 1.03)  # injects +3% residual

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "|residual|=" in txt and ">0.5%" in txt


def test_refuses_unresolved_remainder_above_5_percent(tmp_path):
    def mut(m, t):
        for r in m["records"]:
            if r["deck"] == "d3a_explicit_cdl_cpml":
                r["constraints"] = "Plain"  # makes ~23% enforceSPs time unresolved!

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut)
    txt = _run_summarize(snap, tmp_path / "out.md")
    assert "unresolved newStep_other=" in txt and ">5%" in txt


def test_coverage_gaps_does_not_subtract_non_elem_phases(tmp_path):
    """Regression guard: _compute_coverage_gaps must NOT subtract non-element phases
    (soe.trisolve, soe.factor, commit.elements) from elem_kernel when computing
    uncovered element shares."""
    txt_smoke = _run_summarize(SMOKE_SNAP, tmp_path / "smoke.md", ["--max-load", "100"])
    # On d3c_explicit_sms_cons[smoke], elem_kernel=80.4% (LadrunoBrick) and 260b targets
    # soe.trisolve (10.4%), NOT elem_kernel -> uncovered elem_kernel must be 80.4%, not 69.9%.
    assert "`d3c_explicit_sms_cons[smoke]` uncovered `elem_kernel` ≈ 80.4%" in txt_smoke
    # On d3a[smoke], if LadrunoBrick (33002) is 15% while commit.elements is ~19%,
    # commit.elements (targeted by 260b) must NOT mask uncovered LadrunoBrick!
    def mut_d3a(m, t):
        for r in m["records"]:
            if r["deck"] == "d3a_explicit_cdl_cpml" and r["mode"] == "1t_deep":
                step = t[r["run_id"]]["tree"]["children"]["step"]
                eu = step["children"]["solveCurrentStep"]["children"]["update"]["children"]["elem.update"]
                for row in eu["elem_by_type"]:
                    if row[0] == 33002:
                        row[2] = int(step["wall_ns"] * 0.15)
                        eu["wall_ns"] += row[2]
                        step["children"]["solveCurrentStep"]["children"]["update"]["wall_ns"] += row[2]
                        step["children"]["solveCurrentStep"]["wall_ns"] += row[2]
                        step["wall_ns"] += row[2]

    snap = _make_clean_snap(tmp_path, QUIET_SNAP, mut_d3a)
    txt = _run_summarize(snap, tmp_path / "d3a_gap.md")
    assert "`d3a_explicit_cdl_cpml[smoke]` uncovered `elem_kernel`" in txt
    assert "dominant `LadrunoBrick`" in txt
