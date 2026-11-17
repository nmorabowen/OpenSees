"""WP-232i: the committed warrant run (tests/data/mpm/bspline/warrant.json, written by
tests/mpm_bspline_warrant.py on the 232i build) against the pre-registered gates.

The gate is read from tests/mpm_bspline_cases.txt, never from the JSON, and the order p and the
monotonicity are recomputed here from the stored e and h (critic A3), so an edited verdict or p in
the file cannot pass. A gated MMS row missing from the file fails; every row must come from one
build. Skipped only until the warrant run has been committed. Not named test_ladruno_mpm_*: it
reads a file, not the binary, so it does not belong in the mutation suite.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_bspline_harness as B  # noqa: E402
import mpm_mms as M  # noqa: E402

pytestmark = pytest.mark.zone_a

DATA = Path(__file__).resolve().parent / "data" / "mpm" / "bspline"
ROWS = B.parse_rows()
GATED = ["G2", "G3", "G4", "G5-G2", "G5-G3", "G5-G4"]


def _warrant():
    p = DATA / "warrant.json"
    if not p.exists():
        pytest.skip("the warrant run is not committed yet (tests/mpm_bspline_warrant.py)")
    return json.loads(p.read_text())


def _recomputed(r):
    es, hs = r["e"], r["h"]
    assert len(es) == len(hs) == len(r["rungs"])
    assert all(abs(h - 1.0 / n) < 1e-15 for h, n in zip(hs, r["rungs"]))
    p = M.order(hs, es) if len(es) >= 3 else None
    nonincr = all(es[i + 1] <= es[i] for i in range(len(es) - 1))
    return dict(r, p=p, nonincr=nonincr)


def test_one_build_for_every_row():
    w = _warrant()
    builds = {r.get("ladrunoBuild") for r in w["rows"].values()}
    assert builds == {w["ladrunoBuild"]}, builds
    assert w.get("ladrunoMutation") in ("none", "none/unknown"), w.get("ladrunoMutation")


@pytest.mark.parametrize("name", GATED)
def test_gated_row(name):
    w = _warrant()
    assert name in w["rows"], f"{name} is missing from the warrant run"
    r = _recomputed(w["rows"][name])
    row = ROWS[name]
    assert r["rungs"][:4] == B.as_list(row["ladder"]), r["rungs"]
    if row.get("validity") is not None:
        assert r.get("validity_change") is not None and r["validity_change"] < float(row["validity"]), \
            r.get("validity_change")
        assert r.get("validity_reached") is not False
    ok, why = B.verdict(row, r)
    assert ok, f"{name}: {why}; e = {r['e']}"


@pytest.mark.parametrize("mode", ["zero", "scale"])
def test_g6_mutants_turn_g2_or_g3_red(mode):
    p = DATA / f"g6_{mode}.json"
    if not p.exists():
        pytest.skip(f"the {mode.upper()} mutant run is not committed yet")
    d = json.loads(p.read_text())
    assert d.get("ladrunoMutation", "").upper().endswith(mode.upper()), d.get("ladrunoMutation")
    red = [n for n, r in d["rows"].items() if B.verdict(ROWS[n], _recomputed(r))[0] is not True]
    assert set(red) & {"G2", "G3"}, red
