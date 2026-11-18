"""WP-232k S3 from the committed records (tests/data/mpm/soil/s3_mpm_<mat>_r<r>_<al>.json against
232k-V's s3_fe_<mat>_n<4|8|16>.json). No binary needed.

Gate (pre-registered): q at s/B = 0.01, 0.02, 0.05 within +-10 % of the FE B-bar record at that
rung, and the run reaches s = 0.05 B (a stop is a miss). HS: r1 and r2 gated (anti-locking on),
r3 and the anti-locking-off legs reported. SANISAND: waits for its FE deck.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.zone_a
DATA = Path(__file__).resolve().parent / "data" / "mpm" / "soil"
FE_N = {1: 4, 2: 8, 3: 16}


def _pair(mat, r, al):
    m = DATA / f"s3_mpm_{mat}_r{r}_{al}.json"
    f = DATA / f"s3_fe_{mat}_n{FE_N[r]}.json"
    if not m.exists() or not f.exists():
        pytest.skip(f"{m.name} or {f.name} not committed")
    return json.loads(m.read_text(encoding="utf-8")), json.loads(f.read_text(encoding="utf-8"))


@pytest.mark.parametrize("r", [1, 2])
def test_s3_hs_gated(r):
    mpm, fe = _pair("hs", r, "fbar")
    assert mpm["reached"], mpm.get("stopped")
    for sb in ("0.01", "0.02", "0.05"):
        q, qr = mpm["q_at"][sb], fe["q_at"][sb]
        assert abs(q / qr - 1.0) <= 0.10, (sb, q, qr)


@pytest.mark.parametrize("r,al", [(3, "fbar"), (1, "off"), (2, "off")])
def test_s3_hs_reported(r, al):
    """Reported legs: the record exists and carries q at the three settlements it reached."""
    mpm, fe = _pair("hs", r, al)
    assert "q_at" in mpm and "failed_attempts" in mpm
    for sb, q in mpm["q_at"].items():
        assert q > 0.0
