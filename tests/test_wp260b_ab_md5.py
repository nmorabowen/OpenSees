"""Zone-A tests for WP-260b follow-up (#1167 items 1 & 3) in wp260b_ab_md5.py.

Verifies:
  1. Item 1: `wp260b_ab_md5.py` exits non-zero with a clear message when
     `before_sha == after_sha` (same `dist/bin` passed for both sides), and
     enforces `--expect-before` / `--expect-after` (>= 9 char prefix match).
  2. Item 3: `t7_quad_all`, `t7_shell_callers`, and `brick_formuri_damage`
     each emit >= 200 finite non-zero values (with `brick_formuri_damage`
     asserting Concrete3D `damage > 0.5` on `formUri(0)`).
"""
from __future__ import annotations

import importlib.util
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

ROOT = Path(__file__).resolve().parents[1]
HARNESS = ROOT / "Ladruno_scripts" / "wp260b_ab_md5.py"
BIN_DIR = os.path.dirname(os.path.abspath(ops.__file__))


def _load_harness_module():
    spec = importlib.util.spec_from_file_location("wp260b_ab_md5", HARNESS)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run_cli(*extra_args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(HARNESS), *extra_args],
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=60,
        encoding="utf-8",
        errors="replace",
    )


def test_same_bin_dir_exits_nonzero_with_clear_message():
    build_sha = str(ops.ladrunoBuild()).strip()
    cp = _run_cli("--before-bin", BIN_DIR, "--after-bin", BIN_DIR)
    assert cp.returncode != 0, f"expected non-zero exit on same build, got 0:\n{cp.stdout}\n{cp.stderr}"
    assert "same ladrunoBuild SHA" in cp.stderr
    assert build_sha in cp.stderr

    cp_short = _run_cli("--before", BIN_DIR, "--after", BIN_DIR)
    assert cp_short.returncode != 0
    assert "same ladrunoBuild SHA" in cp_short.stderr


def test_expect_before_and_expect_after_prefix_checks():
    actual = str(ops.ladrunoBuild()).strip()
    assert len(actual) >= 9, f"unexpected short ladrunoBuild SHA: {actual!r}"
    good_prefix = actual[:9]
    bad_prefix = "fffffffff" if not actual.lower().startswith("fffffffff") else "000000000"

    # Rejects < 9 characters on both flags
    cp_short_b = _run_cli("--before-bin", BIN_DIR, "--after-bin", BIN_DIR, "--expect-before", actual[:8])
    assert cp_short_b.returncode != 0
    assert "--expect-before requires at least 9 characters" in cp_short_b.stderr

    cp_short_a = _run_cli("--before-bin", BIN_DIR, "--after-bin", BIN_DIR, "--expect-after", actual[:8])
    assert cp_short_a.returncode != 0
    assert "--expect-after requires at least 9 characters" in cp_short_a.stderr

    # Rejects mismatched --expect-before
    cp_bad_b = _run_cli("--before-bin", BIN_DIR, "--after-bin", BIN_DIR, "--expect-before", bad_prefix)
    assert cp_bad_b.returncode != 0
    assert "does not match --expect-before prefix" in cp_bad_b.stderr

    # Accepts matching 9-char --expect-before and rejects mismatched --expect-after
    cp_bad_a = _run_cli(
        "--before-bin", BIN_DIR, "--after-bin", BIN_DIR,
        "--expect-before", good_prefix, "--expect-after", bad_prefix,
    )
    assert cp_bad_a.returncode != 0
    assert "does not match --expect-after prefix" in cp_bad_a.stderr

    # Matching both 9-char prefixes passes prefix checks and then hits same-SHA guard
    cp_same = _run_cli(
        "--before-bin", BIN_DIR, "--after-bin", BIN_DIR,
        "--expect-before", good_prefix.upper(), "--expect-after", good_prefix,
    )
    assert cp_same.returncode != 0
    assert "same ladrunoBuild SHA" in cp_same.stderr


def test_enlarged_cases_emit_at_least_200_values():
    mod = _load_harness_module()
    expected_min = {
        "t7_quad_all": 200,
        "t7_shell_callers": 200,
        "brick_formuri_damage": 200,
    }
    for case, min_n in expected_min.items():
        res = mod._run_one(BIN_DIR, case)
        vals = res["vals"]
        assert len(vals) >= min_n, f"{case}: expected >= {min_n} values, got {len(vals)}"
        assert all(math.isfinite(v) for v in vals), f"{case}: non-finite values detected"
        assert max(abs(v) for v in vals) > 0.0, f"{case}: all-zero output"
        assert len(res["md5"]) == 32
