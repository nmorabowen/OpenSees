"""WP-230 (TIMs F25): the SAS-ME low-pressure CLAMP, the comparison mode to WP-152's
separation: `-sasLowP clamp p_floor p_sep p_contact` (+ `-sasLowPAlpha leave|project`).

A point that loses its confinement (WP-152's triggers E1 / E2, with their gates) is held at
p = p_floor and KEEPS its committed deviator, frozen, until it re-contacts at p_contact.

What is pinned (Ladruno_implementation/230_sanisand_lowp_clamp.md):
  (a)  ORACLE FIRST: five element paths follow the state machine of the oracle driver
       (Ladruno_files/testbed/sanisand_lowp_clamp/clamp_oracle.py, tests/data/
       wp230_oracle_paths.json), under both alpha rules: entry and re-contact at the same
       steps; while clamped the deviator is EXACTLY the one committed before the entry and
       p follows p_floor + K(p_contact) max(g, 0); the normal phase agrees to SAS-ME's
       tolerance; the census counts each transition once;
  (b)  the two measures of what a clamped state costs -- f in yield-cone radii
       (`clampMaxFRatio`) and the stress ratio over the bounding ratio (`clampMaxBound`) --
       agree with the oracle; under `project` f is 0 and the stress path is the SAME as
       under `leave` (alpha does not enter a clamped point's stress);
  (c)  clamp against separation on the same path: same entry step, and the clamped point
       carries the deviator the separated one drops;
  (d)  OFF is byte-identical: the WP-152 legs (iso / te / cyc with -sasTensionCutoff) give,
       bit for bit on win32 (1e-6 elsewhere), what a build WITHOUT this WP gave (tests/data/
       wp230_cutoff_byteid_baseline.json), and so does the cutoff absent;
  (e)  E2 (a cost refusal at committed p0 < p_sep under a non-compressing increment) clamps,
       and a clamp entered from p0 > p_floor holds a stress OUTSIDE the bounding surface
       (clampMaxBound > 1) -- the case the A/B is about, which the exact oracle cannot reach;
       held under compression it refuses, as the cutoff does;
  (f)  the parser: exclusive with -sasTensionCutoff, the ranges, the inert alpha rule,
       other schemes, and the cutoff's own refusals (-Presidual, -pRe, no h floor, -implex);
  (g)  the options and the clamped state cross the datastore wire, and the continuation
       is bit-identical;
  (h)  a clamped point's tangent is the elastic stiffness at its state.

Runtime: ~20 s.
"""
import json
import math
import os
import tempfile

import pytest

from _testbed import ops
import wp151_reseat_tools as T151
import wp152_cutoff_tools as T
import wp230_clamp_tools as C

FIX = json.load(open(C.FIXTURE))
PF, PSEP, PCON = FIX["p_floor"], FIX["p_sep"], FIX["p_contact"]
CLAMP = ("-sasLowP", "clamp", PF, PSEP, PCON)
RULES = {"leave": (), "project": ("-sasLowPAlpha", "project")}
E1, E2, EXITS, ACTIVE, MAXFR, MAXBD = 36, 37, 38, 39, 44, 45


def _rel(a, b, floor=1.0):
    n = max(math.sqrt(sum(x * x for x in b)), floor)
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b))) / n


def _dev(s):
    p = (s[0] + s[1] + s[2]) / 3.0
    return [s[0] - p, s[1] - p, s[2] - p, s[3], s[4], s[5]]


def _events(records):
    ev, causes, prev = [], [], [0.0] * 46
    for r in records:
        s = r["sas"]
        for col, cause in ((E1, "E1"), (E2, "E2")):
            if s[col] > prev[col]:
                assert s[col] == prev[col] + 1.0, (r["k"], col)
                ev.append((r["k"], "enter"))
                causes.append(cause)
        if s[EXITS] > prev[EXITS]:
            assert s[EXITS] == prev[EXITS] + 1.0, r["k"]
            ev.append((r["k"], "recontact"))
        prev = s
    return ev, causes


@pytest.fixture(scope="module")
def runs():
    out = {}
    for key, fx in FIX["paths"].items():
        name, rule = key.split("|")
        flip, rec = T.run(fx["increments"], T.BASE + T.R1 + CLAMP + RULES[rule])
        out[key] = (flip, rec)
    ops.wipe()
    return out


# ----------------------------------------------------------------------- (a)
@pytest.mark.parametrize("key", sorted(FIX["paths"]))
def test_element_path_follows_the_oracle_state_machine(runs, key):
    flip, rec = runs[key]
    fx = FIX["paths"][key]
    assert _rel(flip, FIX["flip"]["sigma"]) < 1e-12, "the post-flip state moved: regenerate the fixture"
    assert all(r["rc"] == 0 for r in rec) and len(rec) == len(fx["records"])
    ev, causes = _events(rec)
    assert ev == [(k, "enter" if e.startswith("enter") else e) for k, e in fx["events"]]
    assert set(causes) <= {"E1", "E2"}
    worst, kept = 0.0, None
    for i, (r, o) in enumerate(zip(rec, fx["records"])):
        if o["mode"] == "C":
            if kept is None:                 # the deviator committed just before the entry
                kept = _dev(rec[i - 1]["sigma"] if i else flip)
            d = _dev(r["sigma"])
            assert max(abs(a - b) for a, b in zip(d, kept)) < 1e-13, (r["k"], d, kept)   # FROZEN, exactly
            p = sum(r["sigma"][0:3]) / 3.0
            if o["g"] <= 0.0:
                assert abs(p - PF) < 1e-13, (r["k"], p)
            else:
                assert abs(p - o["p"]) < 1e-9 * max(o["p"], 1.0), (r["k"], p, o["p"])
            assert r["sas"][ACTIVE] == 1.0
        else:
            if o["event"] == "recontact":
                # the exit state keeps the deviator too; p_re as the oracle's
                d = _dev(r["sigma"])
                assert max(abs(a - b) for a, b in zip(d, kept)) < 1e-13, (r["k"], d, kept)
                assert abs(sum(r["sigma"][0:3]) / 3.0 - o["p"]) < 1e-9 * max(o["p"], 1.0)
                kept = None
            # the WP-152 tolerance: 5e-3 relative above 4 kPa, 0.02 kPa absolute below
            worst = max(worst, _rel(r["sigma"], o["sigma"], floor=4.0))
            assert r["sas"][ACTIVE] == 0.0
    assert worst < 5e-3, worst
    last = rec[-1]["sas"]
    assert last[E1] + last[E2] == sum(1 for e in fx["events"] if e[1] == "enter_tension")
    assert last[EXITS] == sum(1 for e in fx["events"] if e[1] == "recontact")


# ----------------------------------------------------------------------- (b)
@pytest.mark.parametrize("name", ["te", "dev_open", "te_dev"])
def test_the_cost_measures_agree_with_the_oracle(runs, name):
    lv, pj = runs[f"{name}|leave"][1][-1]["sas"], runs[f"{name}|project"][1][-1]["sas"]
    ol, op = FIX["paths"][f"{name}|leave"], FIX["paths"][f"{name}|project"]
    assert ol["max_f_ratio"] > 10.0                                  # the stale alpha: many cone radii
    assert lv[MAXFR] == pytest.approx(ol["max_f_ratio"], rel=0.1), (lv[MAXFR], ol["max_f_ratio"])
    assert abs(pj[MAXFR]) < 1e-9 and abs(op["max_f_ratio"]) < 1e-9   # projected: ON the cone
    for got in (lv[MAXBD], pj[MAXBD]):                               # the stress, whatever alpha is
        assert got == pytest.approx(ol["max_bound"], rel=0.1, abs=0.02), (got, ol["max_bound"])
    assert ol["max_bound"] < 1.0                                     # an E1 entry RAISES p: inside the bounding surface


@pytest.mark.parametrize("name", ["iso", "te", "cyc", "dev_open", "te_dev"])
def test_the_alpha_rule_does_not_change_a_clamped_points_stress(runs, name):
    a, b = runs[f"{name}|leave"][1], runs[f"{name}|project"][1]
    fx = FIX["paths"][f"{name}|leave"]["records"]
    k_exit = [r["k"] for r in fx if r["event"] == "recontact"][0]
    for x, y in zip(a[:k_exit + 1], b[:k_exit + 1]):                 # up to and including the first exit
        assert [v.hex() for v in x["sigma"]] == [v.hex() for v in y["sigma"]], x["k"]
    # after re-contact both restart from the cone around the same stress
    assert _rel(a[-1]["sigma"], b[-1]["sigma"], floor=4.0) < 5e-3


def test_iso_path_has_no_deviator_to_keep(runs):
    """Isotropic unloading from an isotropic state: the clamp holds p_floor I, as the cutoff
    holds p_min I -- and f is at the bottom of the cone (-1 cone radius)."""
    _, rec = runs["iso|leave"]
    cl = [r for r, o in zip(rec, FIX["paths"]["iso|leave"]["records"]) if o["mode"] == "C" and o["g"] <= 0.0]
    assert cl and all(max(abs(x - PF) for x in r["sigma"][0:3]) < 1e-13 for r in cl)
    assert all(max(abs(x) for x in _dev(r["sigma"])) < 1e-13 for r in cl)


# ----------------------------------------------------------------------- (c)
def test_clamp_against_separation_on_the_same_path(runs):
    incs = FIX["paths"]["dev_open|leave"]["increments"]
    _, cl = runs["dev_open|leave"]
    _, sp = T.run(incs, T.BASE + T.R1 + ("-sasTensionCutoff", PSEP, PCON))
    ops.wipe()
    ev_c, _ = _events(cl)
    ev_s, _ = _events([dict(r, sas=r["sas"] + [0.0, 0.0][:max(0, 46 - len(r["sas"]))]) for r in sp])
    assert ev_c[0] == ev_s[0]                                        # the same trigger, the same step
    k = ev_c[0][0] + 5                                               # both open
    q_c = math.sqrt(1.5 * sum(x * x for x in _dev(cl[k]["sigma"])))
    q_s = math.sqrt(1.5 * sum(x * x for x in _dev(sp[k]["sigma"])))
    assert q_s < 1e-12 and q_c > 1e-3, (q_c, q_s)                    # the one thing that changes
    assert abs(sum(cl[k]["sigma"][0:3]) / 3.0 - PF) < 1e-13
    assert abs(sum(sp[k]["sigma"][0:3]) / 3.0 - 0.0101) < 1e-13


# ----------------------------------------------------------------------- (d)
@pytest.mark.parametrize("leg", ["cutoff", "none"])
def test_off_is_byte_identical(leg):
    if not os.path.isfile(C.BASELINE):
        pytest.skip("no baseline recorded (tests/wp230_clamp_tools.py --record on a build without WP-230)")
    base = json.load(open(C.BASELINE))
    got = C.byteid_rows(leg)
    ops.wipe()
    assert set(got) == set(base["rows"][leg])
    for name in got:
        ref = base["rows"][leg][name]
        assert len(got[name]) == len(ref), (leg, name, len(got[name]), len(ref))
        # EXACT on win32 (the baseline's platform); the fork's 1e-6 floor elsewhere
        bad = [k for k, (a, b) in enumerate(zip(got[name], ref)) if not T151.rows_equal(a, b)]
        assert not bad, f"{leg}/{name}: {len(bad)} of {len(ref)} steps moved, first at step {bad[0]}"


# ----------------------------------------------------------------------- (e)
CAP1 = T.with_max_substeps(T.BASE, 1)
EV0_LOW = 2.0e-6                           # stage 0 to p0 of a few tenths of a kPa, below p_sep


E_K0 = (2.0e-6 / 1.2, 0.2e-6 / 1.2, 0.2e-6 / 1.2)   # stage 0: eps_v = 2e-6 with a deviator (q/p about 0.6)


@pytest.mark.parametrize("pf,outside", [(0.02, True), (PSEP, False)])
def test_E2_clamps_and_can_hold_a_stress_outside_the_bounding_surface(pf, outside):
    """A cost refusal (code 9, -maxSubsteps 1) at committed p0 < p_sep under isochoric shear, from
    a K0-like state that carries a deviator. The clamp keeps that deviator and moves p from p0 to
    p_floor, so the stress ratio is multiplied by p0 / p_floor:
      p_floor = 0.02 << p0: far past the bounding surface, the state no alpha can make
                            admissible (clampMaxBound > 1) -- what the exact oracle cannot reach;
      p_floor = p_sep > p0: p is RAISED, the stress ratio falls, the stress stays inside
                            (clampMaxBound < 1). With p_sep = p_floor an E2 entry is always this."""
    opts = CAP1 + T.R1 + ("-sasLowP", "clamp", pf, PSEP, PCON)
    incs = [[3.0e-4, -1.5e-4, -1.5e-4, 0, 0, 0]] * 3
    flip, rec = C.run_aniso(incs, opts, E_K0)
    S = [list(ops.eleResponse(1, "material", gp, "sasStats")) for gp in range(1, 9)]
    ops.wipe()
    p0 = sum(flip[0:3]) / 3.0
    q0 = math.sqrt(1.5 * sum(x * x for x in _dev(flip)))
    assert 0.02 < p0 < PSEP and 0.3 < q0 / p0 < 1.0, (p0, q0 / p0)
    assert all(r["rc"] == 0 for r in rec)
    ev, causes = _events(rec)
    assert ev[0] == (0, "enter") and causes[0] == "E2", (ev, causes)
    after = rec[0]["sigma"]
    assert max(abs(a - b) for a, b in zip(_dev(flip), _dev(after))) < 1e-13      # the deviator is kept
    assert abs(sum(after[0:3]) / 3.0 - pf) < 1e-13
    assert all(s[ACTIVE] == 1.0 and s[40] in (4.0, 9.0) for s in S), [s[36:46] for s in S]
    eta = q0 / pf
    if outside:
        assert eta > 5.0, eta                                      # a stress ratio no sand state has
        assert all(s[MAXBD] > 2.0 for s in S), [s[MAXBD] for s in S]
    else:
        assert eta < q0 / p0
        assert all(0.0 < s[MAXBD] < 1.0 for s in S), [s[MAXBD] for s in S]
    assert all(s[MAXFR] > 0.0 for s in S)                          # alpha left as committed: outside the cone
    # isochoric shear of a clamped point: frozen (to round-off: each update rebuilds
    # p I + dev(sigma_committed)), and never a re-contact
    for r in rec[1:]:
        assert max(abs(a - b) for a, b in zip(r["sigma"], after)) < 1e-14, r["k"]
        assert r["sas"][EXITS] == 0.0


def test_E2_is_held_under_compression():
    opts = CAP1 + T.R1 + CLAMP
    flip, rec = T.run([[3.0e-4, 0, 0, 0, 0, 0]], opts, ev0=EV0_LOW)
    s = T.mresp("sasStats")
    ops.wipe()
    assert 0.0 < sum(flip[0:3]) / 3.0 < PSEP
    assert rec[0]["rc"] != 0 and s[ACTIVE] == 0.0 and s[43] >= 1.0 and s[E2] == 0.0, (rec[0]["rc"], s[36:46])


# ----------------------------------------------------------------------- (f)
@pytest.mark.parametrize("flags", [
    ("-sasLowP", "clamp", 0.1, 0.5, 1.0, "-sasTensionCutoff", 0.5, 1.0),     # two modes
    ("-sasTensionCutoff", 0.5, 1.0, "-sasLowP", "clamp", 0.1, 0.5, 1.0),
    ("-sasLowP", "clamp", 0.0, 0.5, 1.0),                                    # p_floor > 0
    ("-sasLowP", "clamp", 1.0, 0.5, 1.0),                                    # p_floor < p_contact
    ("-sasLowP", "clamp", 0.1, 1.0, 1.0),                                    # p_sep < p_contact
    ("-sasLowP", "clamp", 0.005, 0.5, 1.0),                                  # p_floor >= p_min (0.0101)
    ("-sasLowP", "separate", 0.1, 0.5, 1.0),                                 # unknown mode
    ("-sasLowP", "clamp", 0.1, 0.5),                                         # three values
    ("-sasLowPAlpha", "project"),                                            # inert without the clamp
    ("-sasTensionCutoff", 0.5, 1.0, "-sasLowPAlpha", "project"),
    ("-sasLowP", "clamp", 0.1, 0.5, 1.0, "-sasLowPAlpha", "cone"),           # unknown rule
])
def test_parser_refuses(flags):
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, *T.BASE, *T.R1, *flags)
    ops.wipe()


@pytest.mark.parametrize("base,extra", [
    ((1, 0, 1, 1e-7, 1e-4), ()),                                             # not SAS-ME
    (T.BASE[:-4] + ("-Presidual", 0.5, "-honorTolR", 0), T.R1),              # p_r != 0
    (T.BASE, ()),                                                            # no h floor
    (T.BASE + ("-pRe", 1.0), T.R1),                                          # -pRe
])
def test_parser_keeps_the_cutoffs_refusals(base, extra):
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, *base, *extra, *CLAMP)
    ops.wipe()


def test_parser_accepts_and_reports():
    ops.wipe()
    T.build([[0, 0, 0, 0, 0, 0]], T.BASE + T.R1 + CLAMP + RULES["project"])
    o = T.mresp("sasOptions")
    s = T.mresp("sasStats")
    ops.wipe()
    assert len(o) == 15 and o[9:15] == [PSEP, PCON, PCON, 1.0, PF, 1.0], o
    assert len(s) == 46 and s[MAXFR] == 0.0 and s[MAXBD] == 0.0
    T.build([[0, 0, 0, 0, 0, 0]], T.BASE + T.R1 + ("-sasTensionCutoff", PSEP, PCON))
    o = T.mresp("sasOptions")
    ops.wipe()
    assert o[9:15] == [PSEP, PCON, PCON, 0.0, 0.0, 0.0], o


# ----------------------------------------------------------------------- (g)
@pytest.mark.parametrize("rule", ["leave", "project"])
def test_clamp_crosses_the_database_wire(rule):
    fx = FIX["paths"][f"dev_open|{rule}"]
    incs = fx["increments"]
    opts = T.BASE + T.R1 + CLAMP + RULES[rule]
    k_save = fx["events"][0][0] + 10                                 # clamped
    with tempfile.TemporaryDirectory(prefix="ladruno_wp230_", ignore_cleanup_errors=True) as td:
        db = os.path.join(td, "db")
        T.build(incs, opts)
        T.stage0()
        for _ in range(k_save + 1):
            assert ops.analyze(1) == 0
        assert T.mresp("sasStats")[ACTIVE] == 1.0
        opt_saved = T.mresp("sasOptions")
        try:
            ops.database("File", db)
        except Exception as exc:                                     # noqa: BLE001
            pytest.skip(f"database() unsupported in this build: {exc}")
        ops.save(1)
        for _ in range(k_save + 1, len(incs)):
            assert ops.analyze(1) == 0
        ref = (T.sig_comp(), T.mresp("sasStats"))
        T.build(incs, T.BASE + T.R1)                                 # a skeleton WITHOUT the clamp
        assert T.mresp("sasOptions")[12] == 0.0
        ops.database("File", db)
        ops.restore(1)
        assert T.mresp("sasOptions") == opt_saved
        assert T.mresp("sasStats")[ACTIVE] == 1.0
        ops.wipe()
    # the continuation from a restored state needs the Path series of the saved model: rebuild
    # WITH the options and replay to the end -- the reference the restore must agree with
    T.build(incs, opts)
    T.stage0()
    for _ in range(len(incs)):
        assert ops.analyze(1) == 0
    again = (T.sig_comp(), T.mresp("sasStats"))
    ops.wipe()
    assert [v.hex() for v in again[0]] == [v.hex() for v in ref[0]]
    assert again[1] == ref[1]


# ----------------------------------------------------------------------- (h)
def test_tangentEP_of_a_clamped_point_is_Ce_at_its_state():
    fx = FIX["paths"]["dev_open|leave"]
    incs = fx["increments"]
    T.build(incs, T.BASE + T.R1 + CLAMP)
    T.stage0()
    for _ in range(fx["events"][0][0] + 3):
        assert ops.analyze(1) == 0
    assert T.mresp("sasStats")[ACTIVE] == 1.0
    Cep = T.mresp("tangentEP")
    ops.wipe()
    assert len(Cep) == 36
    M = [Cep[6 * i:6 * i + 6] for i in range(6)]
    # isotropic elastic: M[0][0] - M[0][1] = 2G = 2 M[3][3]; symmetric; positive
    assert M[3][3] > 0.0 and M[0][1] > 0.0
    assert M[0][0] - M[0][1] == pytest.approx(2.0 * M[3][3], rel=1e-12)
    assert all(abs(M[i][j] - M[j][i]) <= 1e-12 * M[0][0] for i in range(6) for j in range(6))
    assert all(abs(M[i][j]) <= 1e-12 * M[0][0] for i in range(3) for j in range(3, 6))
