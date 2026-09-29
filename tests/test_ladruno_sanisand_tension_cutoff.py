"""WP-152: the SAS-ME tension cutoff (separation), `-sasTensionCutoff p_sep p_contact`.

What is pinned (Ladruno_implementation/152_sanisand_tension_cutoff.md):
  (a) it changes ONLY the qualifying refusals: on the 643 WP-151 replays, R1 + cutoff
      against R1 alone (the cutoff requires -sasHFloor), every accepted update is
      bit-identical, every refusal that is not E1/E2 still refuses with the same code,
      and every E1/E2 one separates (p = p_min, no shear, alpha = 0); E1 needs
      p0 <= p0max, E2 a non-compressing increment;
  (b) the element paths (isotropic, triaxial extension, 4 open/close cycles) follow the
      state machine defined around the exact WP-134 oracle (tests/data/
      wp152_oracle_paths.json): the entry and re-contact steps exactly, the separated
      stress exactly, the normal phases to SAS-ME's tolerance, the census;
  (c) net work over the cycles is >= 0 and matches the oracle's;
  (d) the separation state crosses the database wire (a restored point stays
      separated and the continuation is bit-identical);
  (e) a separated point's tangentEP is C_e at p_min;
  (f) the parser (+ the review #2/#4/#7 refusals);
  (g) review #1-#3, #5, #9: the E1 bound (a well-confined point opened in one increment
      refuses, counted), E2 held under compression (a code 9 at p0 < p_sep refuses when
      the increment compresses and separates when it opens), the masked code in
      sepLastCode, Newton on a free-node column under gravity through entry and
      re-contact, a plane-strain smoke, isochoric shear while separated, and an
      InitialStateAnalysis revert that keeps the separation (a plain reset clears it).
Runtime: ~15 s."""
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
P0MAX = 5.0 * T.CUTOFF[2]                   # the default E1 bound, 5 p_contact


@pytest.fixture(scope="module")
def replays():
    ops.wipe()
    ops.nDMaterial("LadrunoSANISAND", 20, *T151.P, *T151.EB_OPTS, *T.R1)
    ops.nDMaterial("LadrunoSANISAND", 21, *T151.P, *T151.EB_OPTS, *T.R1, *T.CUTOFF)
    base, cur = {}, {}
    for (_k, key, st, de) in T151.jobs():
        base[key] = T151.row(T151.replay(ops, 20, st, de))
        o = T151.replay(ops, 21, st, de)
        cur[key] = (T151.row(o), o, st, de)
    ops.wipe()
    return base, cur


def test_only_qualifying_refusals_change(replays):
    base, cur = replays
    assert len(cur) == len(base) == 643
    n_same = n_sep = n_ref = n_held = 0
    for key, ref in base.items():
        row, o, st, de = cur[key]
        p0 = sum(st["sigma"][0:3]) / 3.0            # compression positive, p_r = 0
        rc0, code0 = ref[0], float.fromhex(ref[27])
        if rc0 == 0:                                # accepted by R1 alone: bit-identical now
            assert T151.rows_equal(row, ref), key
            n_same += 1
            continue
        e1q = code0 == 6 or (code0 == 3 and p0 <= 0.0)
        e2q = code0 in (4, 9) and p0 < T.CUTOFF[1]
        e1 = e1q and p0 <= P0MAX
        e2 = e2q and not (sum(de[0:3]) > 0.0)
        if e1 or e2:                                # separated: accepted, p = p_min, no shear
            assert int(o["rc"]) == 0, (key, code0)
            s = list(o["sigma"])
            assert all(abs(x - PMIN) < 1e-12 for x in s[0:3]) and all(abs(x) < 1e-12 for x in s[3:6]), (key, s)
            assert all(abs(x) < 1e-15 for x in o["alpha"]), key
            n_sep += 1
        else:                                       # still refused, same code
            assert int(o["rc"]) != 0 and float(o["sas"]["lastRefuseCode"]) == code0, (key, code0, o["rc"])
            n_held += int(e1q or e2q)
            n_ref += 1
    print(f"replays: {n_same} accepted identical, {n_sep} separated, {n_ref} refused ({n_held} held back)")
    assert n_same > 400, (n_same, n_sep, n_ref)


# ----------------------------------------------------------------------- (b), (c)
def _events(records):
    """(step, 'enter' / 'recontact') and the entry causes, from the COMMITTED census
    (counted once per committed transition)."""
    ev, causes, prev = [], [], [0.0] * 44
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
            # SAS-ME controls the stress error against sigma_ref = 1 kPa (-errFloor), so below
            # p ~ 1 kPa its NORMAL phase agrees with the exact oracle in ABSOLUTE terms. Measured
            # on the extension path: <= 0.0125 kPa in the deviator at p 0.7..0.007 kPa, before
            # the separation (pre-existing, not the cutoff; after re-contact <= 4e-6). Hence the
            # 4 kPa floor: 5e-3 relative above, 0.02 kPa absolute below.
            worst = max(worst, _rel(r["sigma"], o["sigma"], floor=4.0))
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
    ((0.5, 1.0, "-sasSepMaxP0", 2.0), True),
    ((0.5, 1.0, "-sasSepMaxP0", 0.2), False),    # p0max below p_sep
    ((0.5, 1.0, "-sasSepMaxP0", 0.0), False),
])
def test_parser(vals, ok):
    ops.wipe()
    args = ("LadrunoSANISAND", 1, *T.PARAMS, *T.BASE, *T.R1, "-sasTensionCutoff", *vals)
    if ok:
        ops.nDMaterial(*args)
    else:
        with pytest.raises(Exception):
            ops.nDMaterial(*args)
    ops.wipe()


def test_parser_refuses_other_schemes_and_implex():
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, 1, 0, 1, 1e-7, 1e-4, *T.R1, "-sasTensionCutoff", 0.5, 1.0)
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, *T.BASE, *T.R1, "-implex", "-sasTensionCutoff", 0.5, 1.0)
    ops.wipe()


def _with_presidual(opts, pr):
    o = list(opts)
    o[o.index("-Presidual") + 1] = pr
    return tuple(o)


@pytest.mark.parametrize("args", [
    pytest.param(T.BASE + T.CUTOFF, id="no-hFloor"),                                  # review #7
    pytest.param(T.BASE + ("-sasHFloor", 0.0) + T.CUTOFF, id="hFloor-0"),
    pytest.param(_with_presidual(T.BASE, 0.5) + T.R1 + T.CUTOFF, id="Presidual"),     # review #4
    pytest.param(T.BASE + T.R1 + ("-sasSepMaxP0", 2.0), id="p0max-without-cutoff"),
])
def test_parser_review_refusals(args):
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *T.PARAMS, *args)
    ops.wipe()


def test_census_names_and_options_response():
    assert sr.SAS_NAMES[36:44] == ["sepEntriesTension", "sepEntriesLowP", "sepExits", "sepActive",
                                   "sepLastCode", "sepMaxP0", "sepHeldHighP", "sepHeldCompressing"]
    ops.wipe()
    T.build([[0, 0, 0, 0, 0, 0]], T.BASE + T.R1 + T.CUTOFF)
    assert len(T.mresp("sasStats")) == 44
    assert T.mresp("sasOptions")[9:12] == [0.5, 1.0, 5.0]        # p0max defaults to 5 p_contact
    T.build([[0, 0, 0, 0, 0, 0]], T.BASE + T.R1 + T.CUTOFF + ("-sasSepMaxP0", 2.5))
    assert T.mresp("sasOptions")[11] == 2.5
    ops.wipe()


# ----------------------------------------------------------------------- (g) review
SEP_TENSION, SEP_LOWP, SEP_EXITS, SEP_ACTIVE = 36, 37, 38, 39
SEP_LAST_CODE, SEP_MAX_P0, SEP_HELD_HIGHP, SEP_HELD_COMP = 40, 41, 42, 43
OPEN_ISO = [[-4.0e-3, -4.0e-3, -4.0e-3, 0, 0, 0]]        # one increment far past p = 0


def _is_pmin(sig, tol=1e-12):
    return all(abs(x - PMIN) < tol for x in sig[0:3]) and all(abs(x) < tol for x in sig[3:6])


def test_E1_separates_below_the_bound_and_names_the_masked_code():
    """review #2/#3: from p0 ~ 2 kPa (< 5 p_contact) one isotropic increment far past p = 0
    separates on E1; the masked refusal (code 6) and p0 are in the census."""
    flip, rec = T.run(OPEN_ISO, T.BASE + T.R1 + T.CUTOFF)
    ops.wipe()
    assert rec[0]["rc"] == 0 and _is_pmin(rec[0]["sigma"])
    s = rec[0]["sas"]
    assert (s[SEP_TENSION], s[SEP_LOWP], s[SEP_ACTIVE], s[SEP_LAST_CODE]) == (1.0, 0.0, 1.0, 6.0), s[36:44]
    assert abs(s[SEP_MAX_P0] - sum(flip[0:3]) / 3.0) < 1e-9, (s[SEP_MAX_P0], flip)
    assert s[SEP_HELD_HIGHP] == 0.0


def test_E1_refuses_above_the_bound():
    """review #2: the same increment with -sasSepMaxP0 below p0 is a step to cut: it
    REFUSES with its own code 6, nothing separates, and the hold is counted."""
    flip, rec = T.run(OPEN_ISO, T.BASE + T.R1 + T.CUTOFF + ("-sasSepMaxP0", 1.0))
    s = T.mresp("sasStats")
    ops.wipe()
    assert sum(flip[0:3]) / 3.0 > 1.0
    assert rec[0]["rc"] != 0
    assert s[SEP_HELD_HIGHP] >= 1.0 and s[SEP_TENSION] == 0.0 and s[SEP_ACTIVE] == 0.0, s[36:44]
    assert s[sr.SAS_NAMES.index("lastRefuseCode")] == 6.0


EV0_LOW = 2.0e-6                   # stage 0 to p0 ~ 0.4 kPa < p_sep = 0.5
CAP1 = T.with_max_substeps(T.BASE, 1)


@pytest.mark.parametrize("de,separates", [
    pytest.param([3.0e-4, 0, 0, 0, 0, 0], False, id="compressing"),
    pytest.param([-3.0e-4, 0, 0, 0, 0, 0], True, id="opening"),
])
def test_E2_is_held_under_compression(de, separates):
    """review #1: a cost refusal (code 9, -maxSubsteps 1) at committed p0 < p_sep separates
    only under a non-compressing increment; under compression it refuses, counted."""
    flip, rec = T.run([de], CAP1 + T.R1 + T.CUTOFF, ev0=EV0_LOW)
    s = T.mresp("sasStats")
    ops.wipe()
    p0 = sum(flip[0:3]) / 3.0
    assert 0.0 < p0 < T.CUTOFF[1], p0
    if separates:
        assert rec[0]["rc"] == 0 and _is_pmin(rec[0]["sigma"])
        assert s[SEP_LAST_CODE] in (4.0, 6.0, 9.0) and s[SEP_ACTIVE] == 1.0, s[36:44]
        assert s[SEP_HELD_COMP] == 0.0
    else:
        assert rec[0]["rc"] != 0
        assert s[SEP_HELD_COMP] >= 1.0 and s[SEP_ACTIVE] == 0.0 and s[SEP_LOWP] == 0.0, s[36:44]
        assert s[sr.SAS_NAMES.index("lastRefuseCode")] in (4.0, 9.0)


def test_newton_column_under_gravity_separates_and_recontacts():
    """review #8/#9: a free-node column solved by Newton with a FORCE test. Opening the top
    separates the upper brick while the lower one carries the weight; closing re-contacts.
    Every step converges; the lower brick never separates; while the upper one is separated
    the lower one's sigma_zz is the weight plus p_min (equilibrium of the free nodes)."""
    top = [-2.9e-5 * (k + 1) / T.N0 for k in range(T.N0)]
    u = top[-1]
    for _ in range(40):
        u += 1.0e-5
        top.append(u)
    for _ in range(60):
        u -= 1.0e-5
        top.append(u)
    T.build_column(top, T.BASE + T.R1 + T.CUTOFF)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(T.N0):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    seen_sep = False
    for k in range(len(top) - T.N0):
        assert ops.analyze(1) == 0, k
        up, lo = T.eresp(2, "sasStats"), T.eresp(1, "sasStats")
        assert lo[SEP_ACTIVE] == 0.0 and lo[SEP_TENSION] + lo[SEP_LOWP] == 0.0, (k, lo[36:44])
        if up[SEP_ACTIVE] == 1.0:
            seen_sep = True
            szz_lo = -T.eresp(1, "stress")[2]
            assert abs(szz_lo - (T.COL_W + PMIN)) < 1e-6, (k, szz_lo)
    up = T.eresp(2, "sasStats")
    ops.wipe()
    assert seen_sep and up[SEP_TENSION] + up[SEP_LOWP] >= 1.0 and up[SEP_EXITS] >= 1.0, up[36:44]
    assert up[SEP_ACTIVE] == 0.0


def test_plane_strain_smoke():
    """review #9: the plane-strain view (quad, PlaneStrain) separates on isotropic opening
    and re-contacts on closing, sitting at p_min in between."""
    e0 = -T.EV0 / 2.0
    hist = [[e0 * (k + 1) / T.N0] * 2 for k in range(T.N0)]
    cur = list(hist[-1])
    for d in [1.0e-5] * 50 + [-1.0e-5] * 70:
        cur = [cur[0] + d, cur[1] + d]
        hist.append(list(cur))
    T.build_plane(hist, T.BASE + T.R1 + T.CUTOFF)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(T.N0):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    was_sep = False
    for k in range(len(hist) - T.N0):
        assert ops.analyze(1) == 0, k
        s = list(ops.eleResponse(1, "material", 1, "sasStats"))
        if s[SEP_ACTIVE] == 1.0:
            was_sep = True
            sig = list(ops.eleResponse(1, "material", 1, "stress"))
            assert all(abs(-x - PMIN) < 1e-12 for x in sig[0:2]) and abs(sig[2]) < 1e-12, (k, sig)
    ops.wipe()
    assert was_sep and s[SEP_TENSION] == 1.0 and s[SEP_EXITS] == 1.0 and s[SEP_ACTIVE] == 0.0, s[36:44]


def test_isochoric_shear_while_separated_never_recontacts():
    """review #6 (documented behaviour): re-contact is volumetric only -- constant-volume
    distortion of a separated point keeps it at p_min; the volumetric closing that
    follows re-contacts it."""
    D = T.D
    incs = [[-D, -D, -D, 0, 0, 0]] * 20 + [[D, -D, 0, 0, 0, 0]] * 50 + [[D, D, D, 0, 0, 0]] * 70
    flip, rec = T.run(incs, T.BASE + T.R1 + T.CUTOFF)
    ops.wipe()
    assert all(r["rc"] == 0 for r in rec) and len(rec) == len(incs)
    for r in rec[20:70]:
        assert r["sas"][SEP_ACTIVE] == 1.0 and _is_pmin(r["sigma"]), r["k"]
        assert r["sas"][SEP_EXITS] == 0.0
    assert rec[-1]["sas"][SEP_EXITS] == 1.0 and rec[-1]["sas"][SEP_ACTIVE] == 0.0


def test_initial_state_analysis_revert_keeps_the_separation():
    """review #5: InitialStateAnalysis off calls revertToStart with the flag still set,
    which keeps the stress -- so it keeps the separation that produced it. A plain
    reset (revertToStart outside ISA) re-initialises the point NORMAL."""
    incs = T.paths()["iso"]
    T.build(incs, T.BASE + T.R1 + T.CUTOFF)
    T.stage0()
    ops.InitialStateAnalysis("on")
    for _ in range(21):
        assert ops.analyze(1) == 0
    assert T.mresp("sasStats")[SEP_ACTIVE] == 1.0
    ops.InitialStateAnalysis("off")
    assert T.mresp("sasStats")[SEP_ACTIVE] == 1.0
    assert _is_pmin(T.sig_comp())
    ops.reset()
    assert T.mresp("sasStats")[SEP_ACTIVE] == 0.0
    ops.wipe()
