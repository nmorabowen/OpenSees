"""WP-135: PDMY01/02/03 ground inside analyze() on a wild Newton iterate.

``setSubStrainRate()`` (PressureDependMultiYield{,02,03}.cpp) sized the
substep loop of ``getStress()`` as ``|d eps| / 1e-5`` (octahedral shear;
1e-4 in PDMY01) or ``d eps_v / 1e-5``, with no cap, converted to ``int``.
A Newton iterate with ``|du| ~ 1e1 .. 1e5`` (KrylovNewton on a step whose
tangent no longer matches the residual) asked for ~1e6 .. 1e9 substeps per
Gauss point per call: analyze() "hung" (WP-133's two-element deck ran over
15 min). WP-135:

* ``setTrialStrain`` / ``setTrialStrainIncr`` return
  ``LADRUNO_MATERIAL_REFUSED`` when the increment needs more than 1e5
  substeps (stage 1 only), so a host that propagates it (quad) fails the
  step and the analysis can cut it;
* ``setSubStrainRate`` caps the loop at 1e5, so a host that SWALLOWS the
  code (SSPquad, Brick, ...) still returns in bounded time.

Below the cap the vanilla expressions run unchanged.

Gates:
  H1  WP-133's two-element pair past the crossing returns in bounded time
      (was: stalled on step 20), with the refusal printed.
  H2  PDMY01 on the plain one-element dense deck (no brake) returns in
      bounded time (was: stalled on step 16).
  H3  an absurd one-step increment (|eps| = 1e3) is refused fast by all
      three PDMY classes in a propagating host (quad: analyze fails), and
      returns in bounded time in a swallowing host (SSPquad).
  H4  with the refusal, a step-cutting driver can take the pair past the
      crossing.
  B1  byte-identity where the cap is not hit: PDMY01/02 decks equal the
      baseline captured from the UNMODIFIED build (WP-133 HEAD d960d7e4f,
      which leaves PDMY01/02 untouched). PDMY03 byte-identity is WP-133's
      own G1 (tests/test_wp133_pdmy03_cs_params.py), re-run unchanged.

Every scenario that ground before WP-135 runs in a SUBPROCESS under a
wall-clock timeout, so a regression fails instead of hanging the session.
"""
import json
import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

_HERE = Path(__file__).resolve().parent
_WT = _HERE.parent
_DIST = str(_WT / "dist" / "bin")
if not os.path.isfile(os.path.join(_DIST, "opensees.pyd")):
    pytest.skip(f"worktree engine not built: {_DIST}", allow_module_level=True)

from _engine import bind_worktree_engine  # noqa: E402
ops = bind_worktree_engine(_DIST)

sys.path.insert(0, str(_HERE))
import wp135_pdmy_decks as W  # noqa: E402

pytestmark = [pytest.mark.zone_a]

REFUSAL = "substeps; refusing it (cut the step)"
PDMY = ["PressureDependMultiYield", "PressureDependMultiYield02",
        "PressureDependMultiYield03"]


def _scenario(name, timeout=120.0):
    """Run one scenario in a child interpreter; never hangs the session."""
    t0 = time.perf_counter()
    try:
        r = subprocess.run(
            [sys.executable, "-S", str(_HERE / "wp135_pdmy_decks.py"), _DIST, "run", name],
            capture_output=True, text=True, timeout=timeout, stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        pytest.fail(f"{name}: still inside analyze() after {timeout:.0f} s "
                    "(the pre-WP-135 substep grind)")
    out = (r.stdout or "") + (r.stderr or "")
    line = [ln for ln in out.splitlines() if ln.startswith("RESULT ")]
    assert r.returncode == 0 and line, out[-3000:]
    return json.loads(line[-1][7:]), out, time.perf_counter() - t0


# ---------------------------------------------------------------- H1
def test_h1_two_element_pair_past_crossing_returns():
    res, out, dt = _scenario("pair")
    # the refusal fails the crossing step instead of grinding through it
    assert REFUSAL in out, out[-2000:]
    assert res["n1"] == res["n2"] < W.D.NSTEP
    assert 15 <= res["n2"] <= 25, res  # WP-133: the brake fires at step ~20
    assert dt < 60.0, dt


# ---------------------------------------------------------------- H2
def test_h2_pdmy01_one_element_deck_returns():
    res, out, dt = _scenario("pdmy01_full")
    assert REFUSAL in out and "PressureDependMultiYield::setTrialStrain" in out
    assert res["n"] == W.PDMY01_GRIND_STEP - 1, res
    assert dt < 60.0, dt


# ---------------------------------------------------------------- H3
@pytest.mark.parametrize("mat", PDMY)
def test_h3_absurd_increment_refused_by_propagating_host(mat):
    res, out, dt = _scenario(f"absurd:{mat}:quad", timeout=60.0)
    assert res["n"] == 0, res  # analyze(1) failed: the step was refused
    assert f"{mat}::setTrialStrain - trial strain increment" in out, out[-2000:]
    assert dt < 30.0, dt


@pytest.mark.parametrize("mat", PDMY)
def test_h3_absurd_increment_bounded_in_swallowing_host(mat):
    # SSPquad drops setTrialStrain's return code; the loop cap alone must
    # keep getStress() bounded (1e5 substeps per call).
    res, out, dt = _scenario(f"absurd:{mat}:SSPquad", timeout=120.0)
    assert REFUSAL in out, out[-2000:]
    assert dt < 90.0, dt


# ---------------------------------------------------------------- H4
def test_h4_step_cutting_takes_the_pair_past_the_crossing():
    res, out, dt = _scenario("pair_cut", timeout=180.0)
    assert res["cuts"] >= 1, res
    assert res["n1"] == res["n2"] == W.D.NSTEP, res


# ---------------------------------------------------------------- B1
def test_b1_pdmy01_02_byte_identical_to_unmodified_build():
    base = json.loads((_HERE / "wp135_pdmy_byteid_baseline.json").read_text())
    assert set(base) == set(W.CASES)
    for name, fn in W.CASES.items():
        got = W.to_hex(fn(ops))
        assert got["n"] == base[name]["n"], name
        assert got == base[name], f"{name}: not byte-identical to the unmodified build"
