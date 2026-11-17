"""WP-232j: the warrant records (tests/data/mpm/antilocking/*.json) against the pre-registered gates.

tests/mpm_antilocking_warrant.py writes them on the build the orchestrator makes; this file
re-reads them on every Zone-A run (< 1 s) so the committed warrant cannot drift from its gate.
A record that is not there yet is skipped LOUDLY with its name, never passed. Gates are the
rows of wp-briefs.md "232j", as copied into tests/mpm_antilocking_cases.txt.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_antilocking_harness as A  # noqa: E402

pytestmark = pytest.mark.zone_a

DATA = Path(__file__).resolve().parent / "data" / "mpm" / "antilocking"


def _load(name):
    p = DATA / f"{name}.json"
    if not p.exists():
        pytest.skip(f"warrant record {p.name} not written yet (tests/mpm_antilocking_warrant.py)")
    return json.loads(p.read_text(encoding="utf-8"))


@pytest.mark.parametrize("name", ["footing_mildcone_r1_al", "footing_mildcone_r2_al",
                                  "footing_mildcone_r3_al"])
def test_l4_record(name):
    """L4: q_max within +-10 % of the FE B-bar q_max at that rung, and s = 0.2 m reached."""
    rec = _load(name)
    c = A.case_by_name(name)
    assert rec["antilocking"] == "fbar" and rec["ladrunoBuild"], rec.get("ladrunoBuild")
    for key in ("failed_attempts", "relaxed_steps", "halvings", "al_averaged", "al_fallback",
                "zone_fallback_share", "wall_s"):
        assert key in rec, key                         # K6 counts and F1 shares are recorded
    q = rec["q"]
    t = rec["t"]
    gate = [e for e in c["expects"] if e["obs"] == "cond_pressure_max" and not e["report"]]
    qmax = max(q)
    for e in gate:
        assert A.judge(e, qmax)[0] == "pass", (name, e["line"], qmax, rec["zone_fallback_share"])
    assert rec["stopped"] is None and abs(t[-1] - 110.0) <= 1e-9, (name, rec["stopped"], t[-1])


def test_l3_r1_record():
    """L3 at r1 is reported (its gate is at r3, which waits for the admission, freeze F2)."""
    rec = _load("footing_j2_r1")
    assert rec["ladrunoBuild"] and "zone_fallback_share" in rec
    print(json.dumps({k: rec[k] for k in ("observables", "failed_attempts", "relaxed_steps",
                                          "halvings", "zone_fallback_share")}, indent=1))


@pytest.mark.parametrize("name", ["mms_G2_al", "mms_G4_al"])
def test_l5_record(name):
    """L5: the MMS orders kept with anti-locking on: p >= 1.8 and e non-increasing."""
    rec = _load(name)
    assert rec["antilocking"] == "fbar" and rec["ladrunoBuild"]
    e = rec["e"]
    assert all(b <= a for a, b in zip(e, e[1:])), e
    assert rec["p"] is not None and rec["p"] >= 1.8, rec["p"]
