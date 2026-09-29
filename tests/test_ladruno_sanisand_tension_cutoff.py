"""WP-152: the SAS-ME tension cutoff (separation), `-sasTensionCutoff p_sep p_contact`.

What is pinned (Ladruno_implementation/152_sanisand_tension_cutoff.md):
  (a) it changes ONLY the qualifying refusals: on the 643 WP-151 replays every accepted
      update is bit-identical with the cutoff given, every refusal that is not E1/E2
      still refuses with the same code, and every E1/E2 one separates (p = p_min, no
      shear, alpha = 0);
  (b) the element paths (isotropic, triaxial extension, 4 open/close cycles) follow the
      state machine defined around the exact WP-134 oracle (tests/data/
      wp152_oracle_paths.json): the entry and re-contact steps exactly, the separated
      stress exactly, the normal phases to SAS-ME's tolerance, the census;
  (c) net work over the cycles is >= 0 and matches the oracle's;
  (d) the separation state crosses the database wire (a restored point stays
      separated and the continuation is bit-identical);
  (e) a separated point's tangentEP is C_e at p_min;
  (f) the parser.
Runtime: ~10 s."""
import json
import math
import os
import tempfile

import pytest

from _testbed import ops
import sanisand_replay as sr
import wp151_reseat_tools as T151
import wp152_cutoff_tools as T

PMIN = 0.0101
FIX = json.load(open(T.FIXTURE))


def _rel(a, b, floor=1.0):
    n = max(math.sqrt(sum(x * x for x in b)), floor)
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))) / n


# ----------------------------------------------------------------------- (a)
@pytest.fixture(scope="module")
def replays():
    base = json.load(open(os.path.join(T.HERE, "data", "wp151_sasme_byteid_baseline.json")))["rows"]
    ops.wipe()
    ops.nDMaterial("LadrunoSANISAND", 21, *T151.P, *T151.EB_OPTS, *T.CUTOFF)
    cur = {}
    for (_k, key, st, de) in T151.jobs():
        o = T151.replay(ops, 21, st, de)
        cur[key] = (T151.row(o), o, st)
    ops.wipe()
    return base, cur


def test_only_qualifying_refusals_change(replays):
    base, cur = replays
    assert len(cur) == len(base) == 643
    n_same = n_sep = n_ref = 0
    for key, ref in base.items():
        row, o, st = cur[key]
        p0 = sum(st["sigma"][0:3]) / 3.0            # compression positive, p_r = 0
        rc0, code0 = ref[0], float.fromhex(ref[27])
        if rc0 == 0:                                # accepted before: bit-identical now
            assert T151.rows_equal(row, ref), key
            n_same += 1
            continue
        e1 = code0 == 6 or (code0 == 3 and p0 <= 0.0)
        e2 = code0 in (4, 9) and p0 < T.CUTOFF[1]
        if e1 or e2:                                # separated: accepted, p = p_min, no shear
            assert int(o["rc"]) == 0, (key, code0)
            s = list(o["sigma"])
            assert all(abs(x - PMIN) < 1e-12 for x in s[0:3]) and all(abs(x) < 1e-12 for x in s[3:6]), (key, s)
            assert all(abs(x) < 1e-15 for x in o["alpha"]), key
            n_sep += 1
        else:                                       # still refused, same code
            assert int(o["rc"]) != 0 and float(o["sas"]["lastRefuseCode"]) == code0, (key, code0, o["rc"])
            n_ref += 1
    assert n_same > 400 and n_ref > 90, (n_same, n_sep, n_ref)   # the NonPosH wall fan still refuses


# ----------------------------------------------------------------------- (b), (c)
def _events(records):
    """(step, 'enter' / 'recontact') and the entry causes, from the COMMITTED census
    (counted once per committed transition)."""
    ev, causes, prev = [], [], [0.0] * 40
    for r in records:
        s = r["sas"]
        for col, cause in ((36, "E1"), (37, "E2")):
            if s[col] > prev[col]:
                assert s[col] == prev[col] + 1.0, (r["k"], col, s[col], prev[col])
                ev.append((r["k"], "enter"))
                causes.append(cause)
        if s[38] > prev[38]:
            assert s[38] == prev[38] + 1.0, (r["k"], s[38], prev[38])
            ev.append((r["k"], "recontact"))
        prev = s
    return ev, causes


@pytest.fixture(scope="module")
def element_runs():
    out = {}
    for name, fx in FIX["paths"].items():
        incs = T.paths()[name]
        assert incs == fx["increments"], f"fixture increments for {name} drifted: regenerate it"
        flip, rec = T.run(incs, T.BASE + T.R1 + T.CUTOFF)
        out[name] = (flip, rec)
    ops.wipe()
    return out


@pytest.mark.parametrize("name", ["iso", "te", "cyc"])
def test_element_path_follows_the_oracle_state_machine(element_runs, name):
    flip, rec = element_runs[name]
    fx = FIX["paths"][name]
    assert _rel(flip, FIX["flip"]["sigma"]) < 1e-12, "the post-flip state moved: regenerate the fixture"
    assert all(r["rc"] == 0 for r in rec) and len(rec) == len(fx["records"])
    ev, causes = _events(rec)
    # the oracle's entry is E1 (its p_floor stop); the C++ may reach its accuracy/cost limit
    # (E2, codes 4/9 at p0 < p_sep) in that SAME step before the exact tension point
    assert ev == [(k, "enter" if e.startswith("enter") else e) for k, e in fx["events"]]
    assert set(causes) <= {"E1", "E2"}
    if name in ("iso", "cyc"):              # elastic isotropic unloading: pure tension
        assert set(causes) == {"E1"}, causes
    worst = 0.0
    for r, o in zip(rec, fx["records"]):
        if o["mode"] == "S" and o["event"] != "recontact":
            assert all(abs(x - PMIN) < 1e-12 for x in r["sigma"][0:3]), (r["k"], r["sigma"])
            assert all(abs(x) < 1e-12 for x in r["sigma"][3:6]), (r["k"], r["sigma"])
            assert r["sas"][39] in (0.0, 1.0)
        else:
            worst = max(worst, _rel(r["sigma"], o["sigma"]))
    assert worst < 5e-3, worst
    last = rec[-1]["sas"]
    assert last[36] + last[37] == sum(1 for e in fx["events"] if e[1] == "enter_tension")
    assert last[38] == sum(1 for e in fx["events"] if e[1] == "recontact")
    assert last[39] == (1.0 if fx["records"][-1]["mode"] == "S" else 0.0)


def test_cycles_net_work_nonnegative_and_as_the_oracle(element_runs):
    flip, rec = element_runs["cyc"]
    incs = T.paths()["cyc"]
    prev, W = flip, 0.0
    for r, de in zip(rec, incs):
        W += sum(0.5 * (a + b) * d for a, b, d in zip(r["sigma"], prev, de))
        prev = r["sigma"]
    assert W >= 0.0
    assert abs(W - FIX["paths"]["cyc"]["work"]) < 5e-3 * max(abs(W), 1e-6) + 1e-7, (W, FIX["paths"]["cyc"]["work"])


# ----------------------------------------------------------------------- (d)
def test_separation_state_crosses_the_database_wire():
    incs = T.paths()["iso"]
    opts = T.BASE + T.R1 + T.CUTOFF
    k_save = 20                                           # separated (entry at 4, re-contact at 96)
    with tempfile.TemporaryDirectory(prefix="ladruno_wp152_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "db")
        T.build(incs, opts)
        T.stage0()
        for _ in range(k_save + 1):
            assert ops.analyze(1) == 0
        assert T.mresp("sasStats")[39] == 1.0
        try:
            ops.database("File", db)
        except Exception as exc:                          # noqa: BLE001
            pytest.skip(f"database() unsupported in this build: {exc}")
        ops.save(1)
        for _ in range(k_save + 1, len(incs)):
            assert ops.analyze(1) == 0
        ref = (T.sig_comp(), T.mresp("sasStats"))
        T.build(incs, opts)                               # a fresh skeleton, NORMAL
        assert T.mresp("sasStats")[39] == 0.0
        ops.database("File", db)
        ops.restore(1)
        assert T.mresp("sasStats")[39] == 1.0             # still separated after the restore
        ops.updateMaterialStage("-material", 1, "-stage", 1)
        for _ in range(k_save + 1, len(incs)):
            assert ops.analyze(1) == 0
        got = (T.sig_comp(), T.mresp("sasStats"))
        ops.wipe()
    assert got[0] == ref[0], (ref[0], got[0])
    assert got[1][36:40] == ref[1][36:40]


# ----------------------------------------------------------------------- (e)
def test_tangentEP_of_a_separated_point_is_Ce_at_p_min():
    incs = T.paths()["iso"]
    T.build(incs, T.BASE + T.R1 + T.CUTOFF)
    T.stage0()
    for _ in range(21):
        assert ops.analyze(1) == 0
    assert T.mresp("sasStats")[39] == 1.0
    C = T.mresp("tangentEP")
    ops.wipe()
    c11, c12, c44 = C[0], C[1], C[21]
    assert abs(C[7] - c11) < 1e-9 * c11 and abs(C[14] - c11) < 1e-9 * c11      # isotropic
    G = c44 if abs(c11 - c12 - 2 * c44) < 1e-9 * c11 else (c11 - c12) / 2.0   # Voigt shear convention
    assert abs((c11 - c12) / 2.0 - G) < 1e-6 * G
    # G at p_min; G uses e_init (mUseCurrentVoidRatioInG is false in every constructor)
    P = T.PARAMS
    e = P[2]
    G_pmin = P[0] * P[8] * (2.97 - e) ** 2 / (1.0 + e) * math.sqrt(PMIN / P[8])
    assert abs(G - G_pmin) < 1e-6 * G_pmin, (G, G_pmin)


# ----------------------------------------------------------------------- (f)
@pytest.mark.parametrize("vals,ok", [
    ((0.5, 1.0), True), ((0.0, 1.0), True),
    ((-0.1, 1.0), False), ((1.0, 1.0), False), ((1.0, 0.5), False),
    ((0.0, 0.005), False),                       # p_contact must exceed p_min (0.0101)
])
def test_parser(vals, ok):
    ops.wipe()
    args = ("LadrunoSANISAND", 1, *T.PARAMS, *T.BASE, "-sasTensionCutoff", *vals)
    if ok:
        ops.nDMaterial(*args)
    else:
        with pytest.raises(Exception):
            ops.nDMaterial(*args)
    ops.wipe()


def test_parser_refuses_other_schemes_and_implex():
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, 1, 0, 1, 1e-7, 1e-4, "-sasTensionCutoff", 0.5, 1.0)
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, *T.BASE, "-implex", "-sasTensionCutoff", 0.5, 1.0)
    ops.wipe()


def test_census_names_and_options_response():
    assert sr.SAS_NAMES[36:40] == ["sepEntriesTension", "sepEntriesLowP", "sepExits", "sepActive"]
    ops.wipe()
    T.build([[0, 0, 0, 0, 0, 0]], T.BASE + T.CUTOFF)
    assert len(T.mresp("sasStats")) == 40
    assert T.mresp("sasOptions")[9:11] == [0.5, 1.0]
    ops.wipe()
