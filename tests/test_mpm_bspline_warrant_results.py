"""WP-232i: the committed warrant run (tests/data/mpm/bspline/warrant.json, written by
tests/mpm_bspline_warrant.py on the 232i build) against the pre-registered gates.

A gated row that misses stays a miss (a strict xfail would hide nothing either way: the gate is
read from tests/mpm_bspline_cases.txt, never from the JSON). Skipped until the warrant run has
been committed. Not named test_ladruno_mpm_*: it reads a file, not the binary, so it does not
belong in the mutation suite.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_bspline_harness as B  # noqa: E402

pytestmark = pytest.mark.zone_a

DATA = Path(__file__).resolve().parent / "data" / "mpm" / "bspline"
ROWS = B.parse_rows()
GATED = ["G2", "G3", "G4", "G5-G2", "G5-G3", "G5-G4"]


def _warrant():
    p = DATA / "warrant.json"
    if not p.exists():
        pytest.skip("the warrant run is not committed yet (tests/mpm_bspline_warrant.py)")
    return json.loads(p.read_text())


@pytest.mark.parametrize("name", GATED)
def test_gated_row(name):
    w = _warrant()
    r = w["rows"].get(name)
    if r is None:
        pytest.skip(f"{name} not in the warrant run")
    row = ROWS[name]
    assert r["rungs"][:4] == B.as_list(row["ladder"]), r["rungs"]
    if r.get("validity_change") is not None:
        assert r["validity_change"] < float(row["validity"]), r["validity_change"]
    ok, why = B.verdict(row, r)
    assert ok, f"{name}: {why}; e = {r['e']}"


@pytest.mark.parametrize("mode", ["zero", "scale"])
def test_g6_mutants_turn_g2_or_g3_red(mode):
    p = DATA / f"g6_{mode}.json"
    if not p.exists():
        pytest.skip(f"the {mode.upper()} mutant run is not committed yet")
    d = json.loads(p.read_text())
    assert d["G6_killed"], d.get("red")
    for n, r in d["rows"].items():
        ok, _ = B.verdict(ROWS[n], r)
        assert ok == r["verdict"]
