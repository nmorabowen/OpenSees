"""WP-127 (TIMs F20(a) + F21): the per-instance ModifiedEuler census
(`substepStats`) and the material-point state replay
(`ladrunoSANISANDReplay`, helper `Ladruno_scripts/sanisand_replay.py`).

What is pinned here:

1. BYTE IDENTITY.  The counters and the (null-by-default) trace hook live in
   the vanilla `ManzariDafalias` integrator.  Five material-point decks
   (`wp127_sanisand_byteid.py`: vanilla ManzariDafalias, LadrunoSANISAND 3D,
   the TIMs campaign set in plane strain incl. a huge-increment leg and a
   reversal leg) must reproduce, bit for bit (`float.hex`), the numbers the
   PRE-WP-127 binary produced (`wp127_sanisand_byteid_baseline.json`).
2. THE CENSUS.  `substepStats` has 17 documented columns; the cumulative ones
   are monotone, the census closes (every substep attempt ends in exactly one
   outcome), and the counters SURVIVE revertToLastCommit (a failed analyze)
   but are zeroed by revertToStart (`reset`).
3. THE REPLAY.  A zero increment returns the loaded state; both sign
   conventions and both wrappers agree; a replayed analysis step reproduces
   the analysis' own next committed stress; the trace closes against the
   census; non-deviatoric alpha is projected with a warning.

MEASURED WALL TIME: ~5 s for the file on the dev box.
"""
import json
import math
import os
import sys

import pytest

from _testbed import ops
import test_ladruno_sanisand as sani
import wp127_sanisand_byteid as byteid

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "Ladruno_scripts"))
import sanisand_replay as sr  # noqa: E402

pytestmark = [pytest.mark.zone_a]

_NSTATS = 17
_CUM = slice(0, 13)          # cumulative columns (since revertToStart)
(UPD, MECALLS, SUB, ACC, REJ, FORCED, CLAMP, REJLOWP, ABANDON, CAP, ENTRY,
 PNRESET, MAXONE, LSUB, LFORCED, LABANDON, LCAP) = range(_NSTATS)


def _stats(ele=1, gp=1):
    s = list(ops.eleResponse(ele, "material", gp, "substepStats"))
    assert len(s) == _NSTATS, s
    return s


def _census_closes(s):
    return s[SUB] == (s[ACC] + s[REJ] + s[FORCED] + s[REJLOWP] + s[ABANDON]
                      + s[CAP])


# ---------------------------------------------------------------------------
#  1. byte identity against the pre-WP-127 binary
# ---------------------------------------------------------------------------

def test_counters_are_byte_identical():
    with open(byteid.BASELINE) as fh:
        ref = json.load(fh)["decks"]
    got = byteid.run_all()
    assert sorted(got) == sorted(ref)
    for name in ref:
        assert len(got[name]) == len(ref[name]), name
        for k, (a, b) in enumerate(zip(got[name], ref[name])):
            assert a == b, f"deck {name} row {k}: first differing entry " \
                f"{next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)}"


# ---------------------------------------------------------------------------
#  2. the census
# ---------------------------------------------------------------------------

def test_substep_stats_monotone_and_closed():
    incs = ([(-2e-4, 2e-4)] * 10 + [(2e-4, -2e-4)] * 20 + [(-2e-4, 2e-4)] * 20)
    byteid._build_ps(byteid._CAMPAIGN, byteid._CAMPAIGN_OPTS, 10, 3.0e-6, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    assert _stats()[SUB] == 0          # elastic stage never enters ModifiedEuler
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    prev = _stats()
    for step in range(len(incs)):
        assert ops.analyze(1) == 0, step
        s = _stats()
        for i in range(13):
            assert s[i] >= prev[i], (step, i, s, prev)
        assert _census_closes(s), s
        assert s[MAXONE] >= s[LSUB]
        prev = s
    assert s[SUB] > s[MECALLS] > 0, s   # the leg really substepped
    assert s[UPD] >= s[MECALLS]
    # every Gauss point is its own instance: GP 1 and GP 3 of a homogeneous
    # quad count the same work, but separately
    assert _stats(gp=3)[SUB] == s[SUB]
    # revertToStart zeroes the census
    ops.reset()
    assert _stats() == [0.0] * _NSTATS


def _free_quad(max_substeps):
    """A LadrunoQuad-free, vanilla `quad` with GENUINE free DOFs (top edge
    LOADED), so a material refusal fails Newton and analyze() reverts."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j, (x, y) in enumerate(sani._XY):
        ops.node(j + 1, x, y)
    opts = list(byteid._CAMPAIGN_OPTS)
    opts[opts.index("-maxSubsteps") + 1] = max_substeps
    ops.nDMaterial("LadrunoSANISAND", 1, *byteid._CAMPAIGN, *opts)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    for j, (x, y) in enumerate(sani._XY):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j, (x, y) in enumerate(sani._XY):   # confining pressure 10 kPa
        ops.load(j + 1, -5.0 if x == 1. else 0.0, -5.0 if y == 1. else 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.1)
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.loadConst("-time", 0.0)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(sani._XY):   # deviatoric push on the top edge
        if y == 1.:
            ops.load(j + 1, 0.0, -5.0)


def test_counters_survive_a_failed_analyze():
    """F20(a)'s point: after a FAILED analyze (revertToLastCommit plus the
    zero-increment settle pass) the census must still hold what the failed
    step cost.  The legacy `substeps` response reads 0 here."""
    _free_quad(max_substeps=2)
    ops.integrator("LoadControl", 1.0)
    rc = ops.analyze(1)
    assert rc < 0, "the capped push was meant to fail"
    s = _stats()
    assert s[CAP] >= 1 and s[SUB] >= 3, s
    assert s[LCAP] == 1.0 and s[LSUB] >= 3, s
    assert _census_closes(s), s
    legacy = list(ops.eleResponse(1, "material", 1, "substeps"))
    assert legacy[0] < s[SUB], (legacy, s)


# ---------------------------------------------------------------------------
#  3. the replay
# ---------------------------------------------------------------------------

@pytest.fixture
def campaign():
    ops.wipe()
    sr.define_campaign_material(ops, 1)
    yield 1
    ops.wipe()


def _rows(path, n=None):
    rows = sr.read_ring_csv(path)
    return rows if n is None else rows[:n]


def test_finding_a_attachment_is_compression_positive():
    for path in sr.RING_CSVS:
        for r in sr.read_ring_csv(path):
            p_csv, plus, minus = sr.check_sign_convention(r)
            assert abs(p_csv - plus) <= 1e-5 * max(1.0, p_csv), (r["element"], r["gp"])
            assert abs(p_csv - minus) > 1e-3


def test_zero_increment_returns_the_loaded_state(campaign):
    r = _rows(sr.RING_CSVS[0])[5]
    for mat_type in ("3D", "PlaneStrain"):
        res = sr.replay_row(ops, campaign, r, [0.0] * 6, mat_type=mat_type)
        assert res["rc"] == 0
        assert res["path"] == "elastic"
        assert res["stats"]["substeps"] == 0
        dev = lambda v: [v[i] - (sum(v[:3]) / 3.0 if i < 3 else 0.0) for i in range(6)]
        for key, want in (("sigma", r["sigma"]), ("alpha", dev(r["alpha"])),
                          ("alpha_in", dev(r["alpha_in"])), ("z", dev(r["z"]))):
            got = res[key]
            scale = max(1.0, max(abs(x) for x in want))
            assert max(abs(a - b) for a, b in zip(got, want)) <= 1e-12 * scale, key
        assert abs(res["e"] - r["e"]) <= 1e-12
        assert abs(res["p"] - r["p_kPa"]) <= 1e-5 * max(1.0, r["p_kPa"])


def test_conventions_and_wrappers_agree(campaign):
    r = _rows(sr.RING_CSVS[1])[10]
    de = sr.probes(1e-5)["shear"]
    a = sr.replay_row(ops, campaign, r, de)
    b = sr.replay(ops, campaign, [-x for x in r["sigma"]], r["alpha"],
                  r["alpha_in"], r["z"], r["e"], [-x for x in de],
                  "tensionPositive")
    c = sr.replay_row(ops, campaign, r, de, mat_type="PlaneStrain")
    for x, y, z in zip(a["sigma"], b["sigma"], c["sigma"]):
        assert x == -y and x == z
    assert a["stats"] == c["stats"]
    assert a["trace"] == c["trace"]


def test_convention_is_required(campaign):
    with pytest.raises(Exception):
        ops.ladrunoSANISANDReplay(campaign, "-sigma", *([1.0] * 6),
                                  "-alpha", *([0.0] * 6), "-alphaIn", *([0.0] * 6),
                                  "-fabric", *([0.0] * 6), "-voidRatio", 0.7,
                                  "-dStrain", *([0.0] * 6))


def test_trace_closes_against_the_census_and_is_bounded(campaign):
    r = _rows(sr.RING_CSVS[0])[0]          # the worst b8 point, el 1950 gp 3
    res = sr.replay_row(ops, campaign, r, sr.probes(1e-4)["shear"])
    s = res["stats"]
    assert len(res["trace"]) == s["substeps"] and res["trace_dropped"] == 0
    by = {}
    for t in res["trace"]:
        by[t["code"]] = by.get(t["code"], 0) + 1
        assert t["atDTmin"] == (t["dT"] == 1e-6)
    assert by.get(0, 0) == s["accepted"] and by.get(1, 0) == s["rejectedErr"]
    assert by.get(2, 0) + by.get(3, 0) == s["forcedAtDTmin"]
    assert by.get(3, 0) == s["forcedClampMc"]
    # bounded: a cap of 3 keeps 3 records and counts the rest as dropped
    if s["substeps"] > 3:
        small = sr.replay_row(ops, campaign, r, sr.probes(1e-4)["shear"], trace=3)
        assert len(small["trace"]) == 3
        assert small["trace_dropped"] == s["substeps"] - 3
        assert small["sigma"] == res["sigma"]      # the trace changes nothing
    off = sr.replay_row(ops, campaign, r, sr.probes(1e-4)["shear"], trace=0)
    assert off["trace"] == [] and off["sigma"] == res["sigma"]


def test_trace_of_alpha_is_projected(campaign, capfd):
    rows = [r for r in _rows(sr.RING_CSVS[0])
            if abs(sum(r["alpha"][:3])) > 1e-4]
    assert rows, "b8 row 1859/2 carries tr(alpha) = 2.3e-3"
    r = rows[0]
    res = sr.replay_row(ops, campaign, r, [0.0] * 6)
    assert abs(res["tr_alpha0"] - sum(r["alpha"][:3])) < 1e-12
    assert abs(sum(res["alpha"][:3])) < 1e-12
    assert "tr(alpha)" in capfd.readouterr().err


def test_replay_reproduces_an_analysis_step():
    """Take the committed state of a real 3D analysis at step k, replay the
    strain increment of step k+1, and get the analysis' own stress at k+1
    (to round-off: the replay re-derives the strain from the void ratio)."""
    incs = byteid._iso_dev(40, 5.0e-3, 0.5)
    byteid._build_3d("LadrunoSANISAND", byteid._PARAMS, (), 10, 3.0e-6, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(20):
        assert ops.analyze(1) == 0

    def grab():
        g = lambda name: list(ops.eleResponse(1, "material", 1, name))
        return dict(sig=g("stress"), eps=g("strain"), alpha=g("alpha"),
                    ain=g("alpha_in"), z=g("fabric"), e=g("state")[24])
    prev_eps = None
    ops.analyze(1)
    k0 = grab()
    prev_eps = k0["eps"]
    assert ops.analyze(1) == 0
    k = grab()
    assert ops.analyze(1) == 0
    k1 = grab()
    dnorm = math.sqrt(sum((a - b) ** 2 for a, b in zip(k["eps"][:3], prev_eps[:3]))
                      + 0.5 * sum((a - b) ** 2 for a, b in zip(k["eps"][3:], prev_eps[3:])))
    de = [a - b for a, b in zip(k1["eps"], k["eps"])]
    res = sr.replay(ops, 1, k["sig"], k["alpha"], k["ain"], k["z"], k["e"], de,
                    "tensionPositive", prev_incr_norm=dnorm)
    assert res["rc"] == 0 and res["path"] in ("plastic", "elasticToPlastic")
    scale = max(abs(x) for x in k1["sig"])
    diff = max(abs(a - b) for a, b in zip(res["sigma"], k1["sig"]))
    assert diff <= 1e-9 * scale, (diff, res["sigma"], k1["sig"])
