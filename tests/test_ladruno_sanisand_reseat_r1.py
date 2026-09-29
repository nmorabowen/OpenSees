"""WP-151 (R1): the opt-in regularization of DM04's alpha_in re-seat singularity
in SAS-ME (IntScheme 129): -sasHFloor c_A, -sasReseatHyst c_rev, -sasSoftCap kappa.

What is pinned (Ladruno_implementation/151_sanisand_reseat_singularity.md):
  (a)  OFF is byte-identical: the flags absent (prototype 11) and given as 0 (12)
       reproduce, bit for bit on win32, 643 replays recorded on a build WITHOUT
       WP-151 (the WP-138 wall fan, the 80 ring rows x 4 probes, the WP-128
       reproducer) -- tests/data/wp151_sasme_byteid_baseline.json;
  (b)  negative control: flags OFF, the wall fan IS refused (102 of 320 on
       win32) -- the gate can fail;
  (c)  R1 ON (c_A 1, c_rev 1, kappa 0.5; and without the cap): 0 of 320 refused,
       and the end stress agrees with the WP-151 oracle (exact Radau on the
       modified model) to SAS-ME's tolerance;
  (d)  each piece alone does NOT clear the fan (the floor alone chatters into
       -maxSubsteps; the hysteresis alone keeps h = 1e10 in its band): the
       ablation the memo reports, and why the flags are recommended together;
  (e)  the three census columns count (and stay 0 when OFF);
  (f)  the parser: bad values and non-129 schemes are refused;
  (g)  a monotonic undrained triaxial chain moves by < 1e-3 q_max with R1 ON;
  (h)  the three options cross the datastore wire (a skeleton WITHOUT them is
       restored WITH them), which also guards the size-keyed datastore trap;
  (i)  after the #868 merge, every option family at once (CPPM on one point,
       WP-129 SAS + R1 on another) survives a database round trip by value, and
       the continuation is bit-identical;
  (j)  the two SAS values (i) cannot continue on (errorVars stress, alphaInMode
       bracket) cross the wire by value.

Runtime: ~1 min (C++ replays).
"""
import json
import math
import os
import sys

import pytest

from _testbed import ops
import wp151_reseat_tools as T
import sanisand_replay as sr


def _nrm(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2 + 2.0 * (v[3] ** 2 + v[4] ** 2 + v[5] ** 2))


@pytest.fixture(scope="module")
def runs():
    ops.wipe()
    T.define(ops)
    out = {}
    for tag in T.PROTOS:
        rows, raw = {}, {}
        for kind, key, st, de in T.jobs():
            o = T.replay(ops, tag, st, de)
            rows[key] = T.row(o)
            raw[key] = (kind, st, o)
        out[tag] = (rows, raw)
    yield out
    ops.wipe()


def _fan(run):
    return {k: v for k, v in run[1].items() if v[0] == "fan"}


def _refused(run, kind="fan"):
    return sum(1 for (knd, st, o) in run[1].values() if knd == kind and o["rc"] != 0)


# ----------------------------------------------------------------------- (a)
@pytest.mark.parametrize("tag", [11, 12])
def test_off_is_byte_identical(runs, tag):
    with open(T.BASELINE) as fh:
        base = json.load(fh)["rows"]
    rows = runs[tag][0]
    assert set(rows) == set(base)
    bad = [k for k in base if not T.rows_equal(rows[k], base[k])]
    assert not bad, f"prototype {tag}: {len(bad)} of {len(base)} replays moved, e.g. {bad[:5]}"


# ----------------------------------------------------------------------- (b)
def test_negative_control_wall_fan_refused_today(runs):
    n = _refused(runs[11])
    if sys.platform == "win32":
        assert n == 102, n
    else:
        assert 90 <= n <= 115, n      # round-off at singular states moves a few
    codes = [o["sas"]["lastRefuseCode"] for (k, st, o) in _fan(runs[11]).values() if o["rc"] != 0]
    assert sum(1 for c in codes if int(c) == 5) >= 0.9 * len(codes), codes   # loadingNonPosH


# ----------------------------------------------------------------------- (c)
@pytest.mark.parametrize("tag,ora", [(13, "T1B1S"), (16, "T1B1")])
def test_R1_integrates_the_wall_fan_and_matches_the_oracle(runs, tag, ora):
    assert _refused(runs[tag]) == 0
    with open(T.ORACLE_FAN) as fh:
        ref = json.load(fh)["fan"]
    rels = []
    for key, (kind, st, o) in _fan(runs[tag]).items():
        r = ref[key][ora]
        assert r["status"] == "ok", (key, r["status"])
        dc = [a - b for a, b in zip(o["sigma"], st["sigma"])]
        rels.append(_nrm([a - b for a, b in zip(dc, r["dsig"])]) / max(_nrm(r["dsig"]), 1e-12))
    rels.sort()
    # SAS-ME at TolR 1e-4 against the exact integration of the same model
    # (measured: median 2.8e-4, max 3.3e-3)
    assert rels[len(rels) // 2] < 1.0e-3, rels[len(rels) // 2]
    assert rels[-1] < 1.0e-2, rels[-1]


def test_R1_leaves_the_inadmissible_entries_refused(runs):
    # b8 1950/2-3 (rho_alpha 6.8/7.3) are refused at entry with or without R1
    assert _refused(runs[13], "ring") == _refused(runs[11], "ring") == 8


# ----------------------------------------------------------------------- (d)
def test_each_piece_alone_does_not_clear_the_fan(runs):
    floor_only, hyst_only = _refused(runs[14]), _refused(runs[15])
    assert floor_only >= 50, floor_only     # measured 87 (Zeno chatter -> -maxSubsteps)
    assert hyst_only >= 90, hyst_only       # measured 102 (h = 1e10 in the band)
    caps = [o["sas"]["lastRefuseCode"] for (k, st, o) in _fan(runs[14]).values() if o["rc"] != 0]
    assert sum(1 for c in caps if int(c) == 9) >= 0.8 * len(caps), caps       # maxSubsteps


# ----------------------------------------------------------------------- (e)
def test_census_columns(runs):
    def tot(tag, col):
        return sum(o["sas"][col] for (k, st, o) in _fan(runs[tag]).values())
    for col in ("hFloored", "hSoftCapped", "reseatHeld"):
        assert tot(11, col) == 0 and tot(12, col) == 0, col
    assert tot(13, "hFloored") > 0 and tot(13, "reseatHeld") > 0
    assert tot(13, "hSoftCapped") > 0      # the cap binds on the wall states (oracle: min H/X 0.17 without it)
    assert tot(16, "hSoftCapped") == 0     # no cap asked for
    assert len(sr.SAS_NAMES) == 44         # WP-151's three, then WP-152's four + the review's four


# ----------------------------------------------------------------------- (f)
@pytest.mark.parametrize("flags", [("-sasHFloor", -1.0), ("-sasReseatHyst", -0.5),
                                   ("-sasSoftCap", 1.0), ("-sasSoftCap", -0.1)])
def test_parser_refuses_bad_values(flags):
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.P, *T.EB_OPTS, *flags)
    ops.wipe()


def test_parser_refuses_the_flags_on_other_schemes():
    ops.wipe()
    for flag in ("-sasHFloor", "-sasReseatHyst", "-sasSoftCap"):
        with pytest.raises(Exception):
            ops.nDMaterial("LadrunoSANISAND", 1, *T.P, 1, 0, 1, 1e-7, 1e-4, flag, 0.5)
    ops.nDMaterial("LadrunoSANISAND", 1, *T.P, *T.EB_OPTS, *T.R1_ON)
    ops.wipe()


@pytest.mark.parametrize("mode,flags,ok", [
    ("bracket", ("-sasReseatHyst", 1.0), False),   # no re-seat happens in bracket: inert
    ("stale", ("-sasReseatHyst", 1.0), False),     # nor in stale (ModifiedEuler's alpha_in)
    ("stale", ("-sasHFloor", 1.0), False),         # stale keeps ModifiedEuler's h: inert
    ("bracket", ("-sasHFloor", 1.0), True),        # the floor replaces the 1e10 bracket
    ("stale", ("-sasSoftCap", 0.5), True),         # the cap acts in every mode
    ("bracket", ("-sasReseatHyst", 0.0), True),    # given as 0 = OFF, not inert
])
def test_parser_refuses_R1_flags_the_alpha_in_mode_makes_inert(mode, flags, ok):
    """#868's rule (a flag nothing reads is refused, not ignored) inside SAS-ME:
    the re-seat threshold needs -sasAlphaIn reseat, the floor anything but stale."""
    ops.wipe()
    args = ("LadrunoSANISAND", 1, *T.P, *T.EB_OPTS, "-sasAlphaIn", mode, *flags)
    if ok:
        ops.nDMaterial(*args)
    else:
        with pytest.raises(Exception):
            ops.nDMaterial(*args)
    ops.wipe()


# ----------------------------------------------------------------------- (g)
def _undrained_tc(tag, n=150, de=2.0e-4):
    st = dict(sigma=[100.0, 100.0, 100.0, 0.0, 0.0, 0.0], alpha=[0.0] * 6, alpha_in=[0.0] * 6,
              z=[0.0] * 6, e=T.P[2])
    q = []
    for _ in range(n):
        o = T.replay(ops, tag, st, [de, -0.5 * de, -0.5 * de, 0.0, 0.0, 0.0])
        assert o["rc"] == 0
        st = dict(sigma=list(o["sigma"]), alpha=list(o["alpha"]), alpha_in=list(o["alpha_in"]),
                  z=list(o["z"]), e=o["e"])      # undrained: e does not change
        q.append(o["q"])
    return q


def test_R1_leaves_a_monotonic_undrained_triaxial_unchanged():
    ops.wipe()
    T.define(ops, tags=[11, 13])
    q0, q1 = _undrained_tc(11), _undrained_tc(13)
    ops.wipe()
    dq = max(abs(a - b) for a, b in zip(q0, q1)) / max(abs(x) for x in q0)
    assert dq < 1.0e-3, dq


# ----------------------------------------------------------------------- (h)
def test_R1_options_cross_the_datastore_wire():
    """The three options travel in recvSelf, not in the deck: the skeleton is
    REBUILT WITHOUT them before the restore, so only the wire can bring them
    back. The restored stress also guards the size-keyed datastore trap
    (LEDGER_quirks, WP-151): a Ladruno block the size of the base's Vector(97)
    overwrites the base state and the restore comes back elsewhere."""
    import tempfile
    import test_ladruno_sanisand as S
    base = (129, 0, 1, 1.0e-7, 1.0e-4, "-Pmin", 0.0101, "-Presidual", 0.0)
    tag = 151
    with tempfile.TemporaryDirectory(prefix="ladruno_wp151_", ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, "r1_rt")
        S._build("LadrunoSANISAND", tag, base + T.R1_ON)
        S._elastic_leg(tag)
        ops.updateMaterialStage("-material", tag, "-stage", 1)
        for step in range(5):
            assert ops.analyze(1) == 0, step
        mid = S._stress()
        opt_saved = list(ops.eleResponse(1, "material", 1, "sasOptions"))
        assert opt_saved[6:9] == [1.0, 1.0, 0.5], opt_saved
        try:
            ops.database("File", dbpath)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f"database() unsupported in this build: {exc}")
        ops.save(1)
        S._build("LadrunoSANISAND", tag, base)          # skeleton WITHOUT the flags
        assert list(ops.eleResponse(1, "material", 1, "sasOptions"))[6:9] == [0.0, 0.0, 0.0]
        ops.database("File", dbpath)
        ops.restore(1)
        after = S._stress()
        opt_after = list(ops.eleResponse(1, "material", 1, "sasOptions"))
        ops.wipe()
    assert opt_after == opt_saved, (opt_saved, opt_after)
    assert S._reldiff(mid, after) <= 1.0e-12, ("restore did not reproduce the saved state", mid, after)


# ----------------------------------------------------------------------- (i)
def test_every_option_family_crosses_the_wire_at_once_after_the_868_merge():
    """#893 review scope (2), after the #868 merge: the SAS options block is 9
    wide (LWIRE_SAS_OPT_N 6 -> 9: WP-129's six, then WP-151's three) and the
    WP-151 layout tag is LWIRE_TAG, the last entry. #868's own two-block test,
    with WP-151 added: point 1 carries NON-default CPPM options (IntScheme 2),
    point 2 NON-default WP-129 SAS options AND R1 (IntScheme 129; alphaInMode
    stays `reseat`, else the R1 options would be inert in the continuation).
    Save, restore into a skeleton with every option at its DEFAULT, and require
    both option responses back by value, the census widths, and the next two
    steps bit-identical to a run that never went through the database.
    errorVars stays `full` here: `-sasErrorVars stress` (WP-129's attribution
    switch, not for production) leaves a committed state (step <= 12) that step
    13 refuses at its start, with or without WP-151 (checked on a ladruno-HEAD
    build), so its
    transport is value-checked in (j), without a continuation."""
    import tempfile
    import test_ladruno_sanisand_cppm_newton as t130
    cppm = (2, 2, 1, 1e-7, 1e-7, "-cppmTangent", "vanilla", "-cppmOnFail", "refuse",
            "-cppmHalvings", 5, "-cppmLineSearch", "on", "-cppmStart", "explicit")
    sas = (129, 0, 1, 1e-7, 1e-4, "-Presidual", 0.0, "-errFloor", 3.0, "-alphaBoundTol", 0.2,
           "-alphaProject", 1, "-alphaEntryTol", 3.0,
           "-sasHFloor", 0.5, "-sasReseatHyst", 2.0, "-sasSoftCap", 0.25,
           "-sasTensionCutoff", 0.3, 0.9)          # WP-152's two options too (no point separates here)
    cppm_def = (2, 2, 1, 1e-7, 1e-7)
    sas_def = (129, 0, 1, 1e-7, 1e-4, "-Presidual", 0.0)
    g = lambda e, name: list(ops.eleResponse(e, "material", 1, name))

    def advance(n):
        for _ in range(n):
            assert ops.analyze(1) == 0

    def stage(s):
        ops.updateMaterialStage("-material", 1, "-stage", s)
        ops.updateMaterialStage("-material", 2, "-stage", s)

    t130._two_cube_model(cppm, sas)
    stage(0)
    advance(10)
    stage(1)
    advance(4)
    cppm_saved, sas_opts_saved = g(1, "cppmOptions"), g(2, "sasOptions")
    widths = (len(g(1, "substepStats")), len(g(2, "sasStats")))
    # #868's wire order, then WP-151's three
    assert sas_opts_saved == [3.0, 0.2, 1.0, 0.0, 0.0, 3.0, 0.5, 2.0, 0.25, 0.3, 0.9, 4.5], sas_opts_saved   # p0max = 5 p_contact
    assert cppm_saved[:6] == [1.0, 5.0, 1.0, 0.0, 1.0, 0.0], cppm_saved
    with tempfile.TemporaryDirectory(prefix="ladruno_wp151_all_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "db")
        try:
            ops.database("File", db)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f"database() unsupported in this build: {exc}")
        ops.save(1)
        advance(2)
        ref = [g(1, "stress"), g(2, "stress")]

        t130._two_cube_model(cppm_def, sas_def)        # every option at its DEFAULT
        assert g(2, "sasOptions") == [-1.0, 0.1, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        assert g(1, "cppmOptions") != cppm_saved
        ops.database("File", db)
        ops.restore(1)
        assert g(1, "cppmOptions") == cppm_saved
        assert g(2, "sasOptions") == sas_opts_saved    # all nine, by value
        assert (len(g(1, "substepStats")), len(g(2, "sasStats"))) == widths
        stage(1)
        advance(2)
        got = [g(1, "stress"), g(2, "stress")]
        ops.wipe()
    assert got == ref, (ref, got)


# ----------------------------------------------------------------------- (j)
def test_the_remaining_sas_options_cross_the_wire_by_value():
    """(i)'s complement: the two SAS option values (i) cannot set on a path it
    continues -- errorVars `stress` (see (i)) and alphaInMode `bracket` (which
    refuses -sasReseatHyst) -- with the floor and the cap, saved after the
    elastic stage and restored into a default skeleton: all nine by value."""
    import tempfile
    import test_ladruno_sanisand_cppm_newton as t130
    cppm_def = (2, 2, 1, 1e-7, 1e-7)
    sas_def = (129, 0, 1, 1e-7, 1e-4, "-Presidual", 0.0)
    sas = sas_def + ("-sasErrorVars", "stress", "-sasAlphaIn", "bracket",
                     "-sasHFloor", 1.5, "-sasSoftCap", 0.3)
    g = lambda name: list(ops.eleResponse(2, "material", 1, name))
    t130._two_cube_model(cppm_def, sas)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    ops.updateMaterialStage("-material", 2, "-stage", 0)
    for _ in range(5):
        assert ops.analyze(1) == 0
    saved = g("sasOptions")
    assert saved == [-1.0, 0.1, 0.0, 1.0, 1.0, 2.0, 1.5, 0.0, 0.3, 0.0, 0.0, 0.0], saved
    with tempfile.TemporaryDirectory(prefix="ladruno_wp151_j_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "db")
        try:
            ops.database("File", db)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f"database() unsupported in this build: {exc}")
        ops.save(1)
        t130._two_cube_model(cppm_def, sas_def)
        assert g("sasOptions") == [-1.0, 0.1, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        ops.database("File", db)
        ops.restore(1)
        got = g("sasOptions")
        ops.wipe()
    assert got == saved, (saved, got)
