"""WP-232l T0 (i, ii): route E against 232l-V's numpy oracle (tests/data/mpm/twophase/t0_oracle.json).

The oracle is written from tests/mpm_twophase_interface.md §1 and §7a alone; the band is
rel 1e-9 per field (max |C++ - oracle| / max |oracle|), case t0_modeA of
tests/mpm_twophase_cases.txt. (i): `mpm debug particleBlocks` for particles 0, 5, 15 at the
initial committed state. (ii): three Newmark steps, every particle field after each.
"""
from __future__ import annotations

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_twophase_harness as H  # noqa: E402

pytestmark = pytest.mark.zone_a


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


def test_t0_i_particle_blocks():
    ops = _ops()
    rec = H.load(H.T0_REC)
    with tempfile.TemporaryDirectory() as tmp:
        rows = H.run_t0_i(ops, rec, tmp)
    for r in rows:
        assert "error" not in r, r
        assert r["max"] <= H.T0_TOL, r


def test_t0_ii_three_steps():
    ops = _ops()
    rows = H.run_t0_ii(ops, H.load(H.T0_REC))
    for r in rows:
        assert r["max"] <= H.T0_TOL, r
