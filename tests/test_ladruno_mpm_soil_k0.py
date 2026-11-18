"""WP-232k K0 and K0p from the committed records (tests/data/mpm/soil/k0_meter_*.json; runs on
esmeralda and the Windows box, each record carrying its build). No binary needed.

Rule (pre-registered): A = m(SANISAND) / m(DP) <= 2, or m(SANISAND) <= 0.01, per mode and basis;
a gated leg with a post-flip floor event is not admitted (D23).
- Mode A (records k0_meter_A_{gimp,bspline2}.json): pass on both bases; no post-flip floor event.
- Mode F uncapped (b94c9d75d): not admitted -- every translated soil leg had post-flip floor
  events, every twin none.
- Mode F capped (D23, 15b109ad8): bspline2 passes; gimp is the recorded miss (owner: SANISAND in
  mode F on B-spline bases only).
- C1 (critic): the admitted passes engaged no SANISAND plasticity (max |d alpha| = 0); the failing
  gimp leg did. K0p (pre-registered) is the plastic meter; its records are asserted once present.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

pytestmark = pytest.mark.zone_a
DATA = Path(__file__).resolve().parent / "data" / "mpm" / "soil"


def _load(name):
    p = DATA / name
    if not p.exists():
        pytest.skip(f"record {name} not committed")
    return json.loads(p.read_text(encoding="utf-8"))


def _gate(ms, md):
    return (ms / md) <= 2.0 or ms <= 0.01


def _leg(basis, mode, leg, capped=False):
    if mode == "A":
        return _load(f"k0_meter_A_{basis}.json")["legs"][leg]
    if capped:
        return _load(f"k0_meter_F_{basis}_capped_{leg}.json")["legs"][leg]
    if leg in ("sanisand", "hs"):
        return _load(f"k0_meter_F_{basis}_{leg}.json")["legs"][leg]
    return _load(f"k0_meter_F_{basis}_elastic_dp_dp_mildcone.json")["legs"][leg]


@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
def test_k0_mode_a_passes(basis):
    s, d = _leg(basis, "A", "sanisand"), _leg(basis, "A", "dp")
    assert _gate(s["translated"]["meter"], d["translated"]["meter"])
    for leg in ("elastic", "dp", "sanisand", "hs", "dp_mildcone"):
        assert _leg(basis, "A", leg)["translated"]["floor_hits"] == 0
    assert s["state_drift"]["dalpha_max"] == 0.0          # C1: no plasticity engaged


@pytest.mark.parametrize("basis", ["gimp", "bspline2"])
def test_k0_mode_f_uncapped_not_admitted(basis):
    """D23: post-flip floor events in every translated soil leg, none in the twins."""
    for leg in ("sanisand", "hs"):
        r = _leg(basis, "F", leg)
        assert r["translated"]["floor_hits"] > 0 and r["twin"]["floor_hits"] == 0


def test_k0_mode_f_capped_bspline2_passes():
    s, d = _leg("bspline2", "F", "sanisand", True), _leg("bspline2", "F", "dp", True)
    assert s["translated"]["floor_hits"] == 0
    assert _gate(s["translated"]["meter"], d["translated"]["meter"])
    assert s["state_drift"]["dalpha_max"] == 0.0          # C1


@pytest.mark.xfail(strict=True, reason=(
    "recorded miss (K0 verdict 2026-11-18, build 15b109ad8): mode F capped gimp, SANISAND 0.18557"
    " (13 post-flip floor events) against DP 0.02335, A = 7.95. Diagnosis: GIMP's two-rows-per-"
    "cell stress pattern under explicit crossing noise, with alpha moved in 70 of 160 particles"
    " (max 0.58); plasticity is not separated from a p-dependent elastic instability. Owner:"
    " SANISAND in mode F on B-spline bases only (mpm analyze refuses gimp / linear)."))
def test_k0_mode_f_capped_gimp():
    s, d = _leg("gimp", "F", "sanisand", True), _leg("gimp", "F", "dp", True)
    assert s["translated"]["floor_hits"] == 0
    assert _gate(s["translated"]["meter"], d["translated"]["meter"])


def test_k0_matched_dp_does_not_flatter():
    """B-row note: the gated DP (SANISAND's stage-0 K, G) meters at or below dp_mildcone."""
    for basis in ("gimp", "bspline2"):
        a = _leg(basis, "A", "dp")["translated"]["meter"]
        b = _leg(basis, "A", "dp_mildcone")["translated"]["meter"]
        assert a <= b


K0P = [("A", "gimp"), ("A", "bspline2"), ("F", "bspline2"), ("F", "bspline3")]


@pytest.mark.parametrize("mode,basis", K0P)
def test_k0p(mode, basis):
    """K0p (pre-registered): valid (>= 50 % of SANISAND particles moved alpha in the pre-shear),
    no post-flip floor event, A_p <= 2 or m <= 0.01."""
    s = _load(f"k0_meter_{mode}_{basis}_k0p_sanisand.json")["legs"]["sanisand"]["translated"]
    d = _load(f"k0_meter_{mode}_{basis}_k0p_dp.json")["legs"]["dp"]["translated"]
    assert s["k0p"]["valid"] and s["k0p"]["alpha_moved_fraction"] >= 0.5
    assert s["meter"] is not None and d["meter"] is not None
    assert s["floor_hits"] == 0
    assert _gate(s["meter"], d["meter"])
