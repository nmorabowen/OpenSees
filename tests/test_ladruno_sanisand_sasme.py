"""WP-129: SAS-ME (IntScheme 129), the Sloan-Abbo-Sheng-style SANISAND integrator.

What is pinned here (test letters follow the WP-129 brief, T1-T5 WP-128 sec. 6):

  (a)  every EXISTING IntScheme is byte-identical to the unmodified WP-127
       binary (hex floats, wp129_sanisand_byteid.py + its baseline JSON);
  (b)  on benign states SAS-ME agrees with WP-134's independent oracle
       (`uw_model`, exact integration) -- tests/data/wp129_uw_model_reference.json.
       RungeKutta45 is NOT a reference (WP-128: dT_min 1e-3 hard-coded,
       Mc-clamp force-accept, no drift correction; and its stages 3/4 never
       compute dAlpha3/dAlpha4, so its alpha weights sum to 301/336);
  T1   the WP-128 reproducer grid: no alpha/alpha^b > 1 + kappa with rc 0,
       in both alpha_in modes; with BOTH fixes ablated the escape comes back
       (the gate can fail); with only the error ablated G alone keeps alpha in;
  T2   WP-128's vertUnload / extShear reversal chains: alpha stays inside, no
       f > 1e-6 returned as success, refusals counted;
  T3/(c) the 80 TIMs ring rows x 8 probes: no f > 1e-6 with rc 0, no alpha
       outside with rc 0, and b8 1950/2-3 (inadmissible) REFUSED with the named
       entry code; the substep table against today's ModifiedEuler (d);
  T4   1950/3 shear+ 1e-5 is refused, not returned with f = 0.0139;
  (e)  a refusal reaches analyze (rc < 0) in a one-element deck, and the
       sasStats census survives the failed step;
  (f)  `tangentEP` matches a one-sided finite difference at a plastic point.

Runtime: ~2-4 min (the ring and the chains are C++ replays).
"""
import json
import math
import os
import sys

import pytest

from _testbed import ops
import wp129_sasme_tools as W
import sanisand_replay as sr

_HERE = os.path.dirname(os.path.abspath(__file__))
KAPPA = W.KAPPA


def _define_all():
    W.define_prototypes(ops)
    ops.nDMaterial("LadrunoSANISAND", 7, *W.P, *W.sas_opts(1e-7))    # tight, for (b)
    ops.nDMaterial("LadrunoSANISAND", 8, *W.P, *W.sas_opts(1e-10))   # tighter, for (f)


@pytest.fixture(scope="module")
def protos():
    _define_all()
    yield
    ops.wipe()


def _sas_sub(o):
    return int(o["sas"]["substeps"]) if o["sas"] else 0


# ---------------------------------------------------------------------- (a)
def _rows_equal(cur, ref):
    """EXACT on the baseline's platform (win32, MSVC); elsewhere the fork's
    1e-6 cross-platform floor on floats (test_adr97_p4_inertness.py:146 --
    GCC/libm differ from MSVC in the last bits), non-floats (rc) exact."""
    if sys.platform == "win32":
        return cur == ref
    if len(cur) != len(ref):
        return False
    for rc_, rr in zip(cur, ref):
        if len(rc_) != len(rr) or rc_[0] != rr[0]:
            return False
        xs = [float.fromhex(x) for x in rc_[1:]]
        ys = [float.fromhex(y) for y in rr[1:]]
        scale = max([abs(y) for y in ys] + [1.0])
        if any(abs(x - y) > 1e-6 * scale for x, y in zip(xs, ys)):
            return False
    return True


def test_existing_schemes_byte_identical():
    """Every existing IntScheme vs the unmodified WP-127 binary (Windows).
    Bit-exact on win32; the 1e-6 cross-platform floor elsewhere."""
    import wp129_sanisand_byteid as B
    with open(B.BASELINE) as fh:
        ref = json.load(fh)["decks"]
    cur = B.run_all()
    assert set(cur) == set(ref)
    bad = []
    for name in ref:
        if name in B.NONDETERMINISTIC:
            # IntScheme 4: vanilla MaxEnergyInc passes UNINITIALISED `nG, nK`
            # into ForwardEuler once it sub-steps, so its plastic rows differ
            # run to run in the SAME process on the unmodified binary too
            # (LEDGER_quirks, WP-129). Only its elastic-stage rows are pinned.
            n0 = B.NONDETERMINISTIC[name]
            if not _rows_equal(cur[name][:n0], ref[name][:n0]):
                bad.append(f"{name}: elastic-stage rows differ")
            continue
        if not _rows_equal(cur[name], ref[name]):
            nrow = sum(1 for x, y in zip(cur[name], ref[name]) if not _rows_equal([x], [y]))
            bad.append(f"{name}: {nrow} rows differ")
    assert not bad, "existing schemes moved: " + "; ".join(bad)


# ---------------------------------------------------------------- parser
def test_sas_flags_refused_on_other_schemes():
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *W.P, 1, 0, 1, 1e-7, 1e-4, "-errFloor", 1.0)
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *W.P, 129, 0, 1, 1e-7, 1e-4, "-implex")
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *W.P, 129, 0, 1, 1e-7, 1e-4, "-sasAlphaIn", "sideways")
    ops.nDMaterial("LadrunoSANISAND", 1, *W.P, 129, 0, 1, 1e-7, 1e-4, "-errFloor", 0.5,
                   "-alphaBoundTol", 0.05, "-alphaProject", 1, "-sasAlphaIn", "bracket")
    ops.wipe()


# ---------------------------------------------------------------------- (b)
# The oracle: WP-134's independent reference integrator (draft PR #872), preset
# `uw_model` (DM04 + the UW constitutive additions U1-U5, the paper's alpha_in
# rule, continuous moduli, SciPy Radau rtol 1e-10 with event detection). The
# 3.12 runner has no scipy, so the oracle's answers are a fixture generated by
# Ladruno_files/testbed/wp129_sasme/gen_uw_model_reference.py (CPython 3.11).
REF = os.path.join(_HERE, "data", "wp129_uw_model_reference.json")


def _ref_cases(kind):
    with open(REF) as fh:
        return [c for c in json.load(fh)["cases"] if c["kind"] == kind]


def _rel(a, b):
    return W.norm([x - y for x, y in zip(a, b)]) / max(W.norm(b), 1e-12)


@pytest.fixture(scope="module")
def tight(protos):
    return 7


def test_benign_agreement_with_oracle(tight):
    """K0 states at 20/50/100 kPa, three directions, 1e-5 and 1e-4: SAS-ME at
    TolR 1e-7 lands on the oracle's stress to BENIGN_TOL (relative)."""
    worst = []
    for c in _ref_cases("benign"):
        assert c["ref"]["status"] == "ok"
        new, o = W.step(ops, tight, c, c["dstrain"])
        assert o["rc"] == 0, c
        worst.append((_rel(new["sigma"], c["ref"]["sigma"]), c["p0"], c["dir"], c["delta"]))
    worst.sort(reverse=True)
    assert worst[0][0] < BENIGN_TOL, worst[:3]


BENIGN_TOL = 5.0e-8   # measured 6.9e-9 at TolR 1e-7


def test_reproducer_matches_oracle(tight):
    c = _ref_cases("reproducer")[0]
    new, o = W.step(ops, tight, c, c["dstrain"])
    assert o["rc"] == 0
    rho = W.alpha_over_b(new["sigma"], new["alpha"], new["e"])
    # oracle: eta 0.534, rho_alpha 0.251
    assert abs(rho - c["ref"]["rho_alpha_end"]) < 0.02 * c["ref"]["rho_alpha_end"], (rho, c["ref"])
    assert abs(W.eta(new["sigma"]) - c["ref"]["eta_end"]) < 0.02 * c["ref"]["eta_end"]


# ---------------------------------------------------------------------- T1
def _t1_grid():
    for ps in (0.0101, 0.1, 1.0, 5.0):
        for delta in (1e-7, 1e-6, 1e-5, 3e-5, 1e-4, 3e-4):
            yield ps, delta


def _t1(tag):
    rows = []
    for ps, delta in _t1_grid():
        st = dict(sigma=[ps, ps, ps, 0, 0, 0], alpha=[0.0] * 6, alpha_in=[0.0] * 6,
                  z=[0.0] * 6, e=0.697787979641054)
        new, o = W.step(ops, tag, st, [0.0, delta, 0.0, 0.0, 0.0, 0.0])
        rows.append(dict(ps=ps, delta=delta, rc=o["rc"], sub=_sas_sub(o),
                         ab=W.alpha_over_b(new["sigma"], new["alpha"], new["e"]),
                         ab_n=W.alpha_over_b_n(new["sigma"], new["alpha"], new["e"]),
                         f=o["f_after"]))
    return rows


@pytest.mark.parametrize("tag", [W.TAG_SAS, W.TAG_SAS_BR])
def test_T1_reproducer_grid(protos, tag):
    rows = _t1(tag)
    bad = [r for r in rows if r["rc"] == 0 and not (r["ab"] <= 1.0 + KAPPA)]
    assert not bad, bad
    badf = [r for r in rows if r["rc"] == 0 and r["f"] > 1e-6]
    assert not badf, badf
    # the smallest reproducer is integrated, not refused, and lands near the
    # alpha-aware answer (WP-128: 0.266 for the port, 0.25 for RK45)
    r = next(r for r in rows if r["ps"] == 0.0101 and r["delta"] == 1e-4)
    assert r["rc"] == 0 and r["ab_n"] < 0.6, r


def test_T1_negative_control_today_escapes(protos):
    """The gate can fail: today's ModifiedEuler on the same grid escapes
    (WP-128: 5.14 at p_s 0.0101, delta 1e-4). With BOTH G and E ablated,
    SAS-ME still stays inside there (0.25): the stage moduli (U9) and the
    true-gradient loading classification (F/U10) remove this escape on their
    own -- WP-134's finding that F is load-bearing, not a compounder."""
    r = next(r for r in _t1(W.TAG_ME) if r["ps"] == 0.0101 and r["delta"] == 1e-4)
    assert r["rc"] == 0 and r["ab_n"] > 2.0, r
    r = next(r for r in _t1(W.TAG_SAS_ABL) if r["ps"] == 0.0101 and r["delta"] == 1e-4)
    assert r["rc"] == 0 and r["ab_n"] < 1.0, r


def test_T1_G_fix_alone_keeps_alpha_inside(protos):
    rows = _t1(W.TAG_SAS_NOE)
    bad = [r for r in rows if r["rc"] == 0 and r["ab_n"] > 1.0]
    assert not bad, bad


# ---------------------------------------------------------------------- T2
@pytest.mark.parametrize("pname", ["vertUnload", "extShear"])
@pytest.mark.parametrize("tag", [W.TAG_SAS, W.TAG_SAS_BR])
def test_T2_reversal_chains(protos, pname, tag):
    st0 = W.k0_state(2.0)
    incs = W.incs_for(W.PATHS[pname], 1e-4, 20)
    h = W.run_chain(ops, tag, st0, incs)
    ok = [x for x in h if x["rc"] == 0]
    assert ok, "every increment refused"
    assert max(x["ab"] for x in ok) <= 1.0 + KAPPA
    assert not [x for x in ok if x["f"] > 1e-6]
    ref = [x for x in h if x["rc"] != 0]
    # every refusal is counted with a named code
    for x in ref:
        assert x["sas"]["refusals"] == 1 and x["sas"]["lastRefuseCode"] > 0


# ------------------------------------------------------------------ T3/(c)
def _ring_rows():
    rows = []
    for path in sr.RING_CSVS:
        for r in sr.read_ring_csv(path):
            r["set"] = os.path.basename(path)
            rows.append(r)
    return rows


def test_T3_ring_states(protos):
    rows = _ring_rows()
    assert len(rows) == 80
    table = []
    for r in rows:
        for delta in (1e-6, 1e-5):
            for pname, de in W.ring_probes(delta).items():
                nc, oc = W.step(ops, W.TAG_SAS, r, de)
                nm, om = W.step(ops, W.TAG_ME, r, de)
                table.append(dict(el=r["element"], gp=r["gp"], set=r["set"], probe=pname,
                                  delta=delta, rc=oc["rc"], sub=_sas_sub(oc),
                                  code=int(oc["sas"]["lastRefuseCode"]),
                                  f=oc["f_after"],
                                  ab=W.alpha_over_b(nc["sigma"], nc["alpha"], nc["e"]),
                                  me_rc=om["rc"], me_sub=int(om["stats"]["substeps"]),
                                  me_f=om["f_after"], me_clamp=int(om["stats"]["forcedClampMc"])))
    out = os.environ.get("WP129_RING_JSON")
    if out:
        with open(out, "w") as fh:
            json.dump(table, fh)
    ok = [t for t in table if t["rc"] == 0]
    assert not [t for t in ok if t["f"] > 1e-6], "f > 1e-6 returned as success"
    assert not [t for t in ok if not (t["ab"] <= 1.0 + KAPPA)], "alpha outside with rc 0"
    bad_rows = [t for t in table if t["set"] == "ring_points_b8.csv" and t["el"] == 1950
                and t["gp"] in (2, 3)]
    assert len(bad_rows) == 16
    assert all(t["rc"] != 0 and t["code"] == 2 for t in bad_rows), bad_rows


def test_T3_ring_matches_oracle(protos):
    """The 78 admissible ring rows x 8 probes against WP-134's uw_model: the
    oracle has 0 escapes and |f| <= 1.7e-7 at exit; SAS-ME must have none
    either, and land near it where both integrate."""
    rows = []
    for c in _ref_cases("ring"):
        new, o = W.step(ops, W.TAG_SAS, c, c["dstrain"])
        rows.append(dict(c=c, rc=o["rc"], code=int(o["sas"]["lastRefuseCode"]), f=o["f_after"],
                         rho=W.alpha_over_b(new["sigma"], new["alpha"], new["e"]),
                         d=_rel(new["sigma"], c["ref"]["sigma"]) if c["ref"]["status"] == "ok" else None))
    ok = [r for r in rows if r["rc"] == 0]
    assert not [r for r in ok if r["f"] > 1e-6]
    rho0 = {(r["c"]["mesh"], r["c"]["element"], r["c"]["gp"]): r["c"]["ref"]["max_rho_alpha"] for r in rows}
    assert not [r for r in ok if r["rho"] > max(1.0, rho0[(r["c"]["mesh"], r["c"]["element"], r["c"]["gp"])]) + 1e-3]
    both = [r for r in ok if r["d"] is not None]
    ds = sorted(r["d"] for r in both)
    out = os.environ.get("WP129_ORACLE_JSON")
    if out:
        with open(out, "w") as fh:
            json.dump([dict(mesh=r["c"]["mesh"], el=r["c"]["element"], gp=r["c"]["gp"],
                            probe=r["c"]["probe"], delta=r["c"]["delta"], rc=r["rc"],
                            code=r["code"], f=r["f"], rho=r["rho"], d=r["d"],
                            ref_status=r["c"]["ref"]["status"]) for r in rows], fh)
    # measured at TolR 1e-4: median 8e-6, p95 5e-5
    assert ds[len(ds) // 2] < 1e-4, ds[len(ds) // 2]
    assert ds[int(0.95 * len(ds))] < 5e-4, ds[int(0.95 * len(ds))]
    assert ds[-1] < 5e-3, ds[-1]


def test_T4_1950_3_shear_plus_refused(protos):
    r = next(x for x in sr.read_ring_csv(sr.RING_CSVS[0]) if x["element"] == 1950 and x["gp"] == 3)
    _, om = W.step(ops, W.TAG_ME, r, [0, 0, 0, 1e-5, 0, 0])
    assert om["rc"] == 0 and om["f_after"] > 1e-3      # today: f = 0.0139 as success
    _, oc = W.step(ops, W.TAG_SAS, r, [0, 0, 0, 1e-5, 0, 0])
    assert oc["rc"] != 0 and int(oc["sas"]["refStartAlpha"]) == 1


# ---------------------------------------------------------------------- (e)
def _quad_deck(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    xy = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
    for j, (x, y) in enumerate(xy):
        ops.node(j + 1, x, y)
    ops.nDMaterial("LadrunoSANISAND", 1, *W.P, *opts)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    for j, (x, y) in enumerate(xy):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j, (x, y) in enumerate(xy):
        if x == 1.:
            ops.sp(j + 1, 1, -1.0e-4)
        if y == 1.:
            ops.sp(j + 1, 2, -1.0e-4)
    ops.constraints("Transformation")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20, 0)
    ops.algorithm("Newton")


def test_refusal_reaches_analyze():
    _quad_deck(W.sas_opts(1e-4)[:5] + ("-maxSubsteps", 1))
    ops.integrator("LoadControl", 0.01)
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    rcs = [ops.analyze(1) for _ in range(20)]
    assert min(rcs) < 0, rcs
    s = dict(zip(sr.SAS_NAMES, ops.eleResponse(1, "material", 1, "sasStats")))
    assert s["refusals"] >= 1 and s["refCap"] >= 1, s     # survives the failed step
    _define_all()   # the module's prototypes for the tests that follow


# ---------------------------------------------------------------------- (f)
def test_tangentEP_matches_finite_difference(protos):
    d = [0.3, 1.0, 0.0, 0.2, 0.0, 0.0]
    st = W.k0_state(50.0)
    prev = 0.0
    for _ in range(6):                     # drive to a plastic state on the cone
        de = [2e-5 * x for x in d]
        st, o = W.step(ops, 8, st, de, prev)
        assert o["rc"] == 0
        prev = W.ncov(de)
    _, o0 = W.step(ops, 8, st, [0.0] * 6, prev)
    Cep = o0["tangent_ep"]
    worst = 0.0
    for dd in (d, [0.2, 1.0, 0.1, 0.0, 0.0, 0.0], [0.3, 1.0, 0.0, 0.5, 0.0, 0.0]):
        h = 1e-9
        de = [h * x for x in dd]
        new, o = W.step(ops, 8, st, de, prev)
        assert o["rc"] == 0 and o["sas"]["elastic"] == 0   # a PLASTIC probe
        fd = [(a - b) / h for a, b in zip(new["sigma"], st["sigma"])]
        an = [sum(Cep[i][j] * dd[j] for j in range(6)) for i in range(6)]
        worst = max(worst, W.norm([a - b for a, b in zip(fd, an)]) / W.norm(an))
    assert worst < 1e-4, worst


# ------------------------------------------------ review of #871 (round 1)
_REVIEW_OPTS = (129, 0, 1, 1e-7, 1e-4, "-Presidual", 0.0, "-Pmin", 0.0101)


def _discard_deck(eletype):
    """The reviewer's p1_discard.py: confine, flip, then unload into tension so
    SAS-ME refuses; `eletype` DISCARDS (SSPquad, stdBrick) or forwards (quad)."""
    ops.wipe()
    if eletype == "stdBrick":
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        xy = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
        for k in range(2):
            for j, (x, y) in enumerate(xy):
                ops.node(4 * k + j + 1, x, y, float(k))
        ops.nDMaterial("LadrunoSANISAND", 1, *W.P, *_REVIEW_OPTS)
        ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
        for k in range(2):
            for j, (x, y) in enumerate(xy):
                ops.fix(4 * k + j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0, 1 if k == 0 else 0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for k in range(2):
            for j, (x, y) in enumerate(xy):
                n = 4 * k + j + 1
                if x == 1.:
                    ops.sp(n, 1, -1.0e-3)
                if y == 1.:
                    ops.sp(n, 2, -1.0e-3)
                if k == 1:
                    ops.sp(n, 3, -2.0e-3)
    else:
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        xy = [(0., 0.), (1., 0.), (1., 1.), (0., 1.)]
        for j, (x, y) in enumerate(xy):
            ops.node(j + 1, x, y)
        ops.nDMaterial("LadrunoSANISAND", 1, *W.P, *_REVIEW_OPTS)
        if eletype == "quad":
            ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
        else:
            ops.element("SSPquad", 1, 1, 2, 3, 4, 1, "PlaneStrain", 1.0)
        for j, (x, y) in enumerate(xy):
            ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
        ops.timeSeries("Linear", 1)
        ops.pattern("Plain", 1, 1)
        for j, (x, y) in enumerate(xy):
            if x == 1.:
                ops.sp(j + 1, 1, -1.0e-3)
            if y == 1.:
                ops.sp(j + 1, 2, -2.0e-3)
    ops.constraints("Transformation"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-12, 20, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.25); ops.analysis("Static")
    for _ in range(4):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    ops.integrator("LoadControl", -0.25)
    def mresp(name):
        if eletype == "SSPquad":   # one material, argv handed to it as-is
            return list(ops.eleResponse(1, name))
        return list(ops.eleResponse(1, "material", 1, name))
    rows = []
    for _ in range(12):
        rc = ops.analyze(1)
        rows.append((rc, mresp("strain"), mresp("stress")))
    s = dict(zip(sr.SAS_NAMES, mresp("sasStats")))
    return rows, s


@pytest.mark.parametrize("eletype", ["SSPquad", "stdBrick", "quad"])
def test_refusal_under_discarding_element_is_not_committed(eletype):
    """Review #871 item 1: before the fix SSPquad returned analyze() == 0 for 9
    steps after the first refusal with the strain climbing at frozen stress.
    Now the first step that meets a refusal fails (discarders: the WP-99 commit
    abort; forwarders: the trial return), and no later step commits a strain the
    material did not integrate."""
    rows, s = _discard_deck(eletype)
    first = next((k for k, r in enumerate(rows) if r[0] < 0), None)
    assert first is not None, [r[0] for r in rows]
    assert s["refusals"] >= 1
    ok_before = [r for r in rows[:first] if r[0] == 0]
    # nothing after the first refusal is reported converged
    assert all(r[0] < 0 for r in rows[first:]), [r[0] for r in rows]
    if ok_before:
        eps_last = ok_before[-1][1]
        # the strain the material reports never moves past its last good commit
        for r in rows[first:]:
            assert max(abs(a - b) for a, b in zip(r[1], eps_last)) < 1e-12, (r[1], eps_last)
    _define_all()


def test_refusal_warning_is_per_instance(capfd):
    """Review item 2: no process-wide warning budget. Every fresh instance (a
    replay's private copy) warns once, however many refused before it."""
    _define_all()
    r = next(x for x in sr.read_ring_csv(sr.RING_CSVS[0]) if x["element"] == 1950 and x["gp"] == 3)
    capfd.readouterr()
    for _ in range(12):
        W.step(ops, W.TAG_SAS, r, [0, 0, 0, 1e-5, 0, 0])
    err = capfd.readouterr().err
    assert err.count("update REFUSED") == 12, err[-500:]


def test_refused_update_state_diagnostics(protos):
    """Review item 7: after a refusal the `last` columns are NaN (no valid end
    state) and the void ratio is the committed one."""
    r = next(x for x in sr.read_ring_csv(sr.RING_CSVS[0]) if x["element"] == 1950 and x["gp"] == 3)
    new, o = W.step(ops, W.TAG_SAS, r, [1e-5, 1e-5, 0, 0, 0, 0])
    assert o["rc"] != 0
    assert math.isnan(o["sas"]["lastF"]) and math.isnan(o["sas"]["lastAlphaRatio"])
    assert abs(o["e"] - r["e"]) < 1e-12


def test_parser_and_runtime_refusals(capfd):
    """Review items 5, 6, 9."""
    ops.wipe()
    for bad in ((385, 0, 1, 1e-7, 1e-4),
                (129, 0, 1, 1e-7, 1e-4, "-errFloor", float("inf")),
                (129, 0, 1, 1e-7, 1e-4, "-reversalTol", 1e-8),
                (129, 0, 1, 1e-7, 1e-4, "-reversalRel", 0.1)):
        with pytest.raises(Exception):
            ops.nDMaterial("LadrunoSANISAND", 1, *W.P, *bad)
    # -reversal* are live under the bracket rule (the P2-5 guard runs there)
    ops.nDMaterial("LadrunoSANISAND", 1, *W.P, 129, 0, 1, 1e-7, 1e-4,
                   "-sasAlphaIn", "bracket", "-reversalTol", 1e-8)
    ops.wipe()
    capfd.readouterr()
    ops.nDMaterial("LadrunoSANISAND", 1, *W.P, 129, 0, 1, 1e-7, 1e-4, "-honorTolR", 1)
    assert "INERT: IntScheme 129" in capfd.readouterr().err
    # runtime IntegrationScheme: out of range, and 129 under -implex
    for opts, newv in (((1, 0, 1, 1e-7, 1e-4), 385.0),
                       ((1, 0, 1, 1e-7, 1e-4, "-maxSubsteps", 1000, "-implex"), 129.0)):
        _quad_deck(opts)
        ops.integrator("LoadControl", 0.01)
        ops.analysis("Static")
        ops.parameter(1, "element", 1, "IntegrationScheme", 1)
        ops.updateParameter(1, newv)
        ops.updateMaterialStage("-material", 1, "-stage", 1)
        for _ in range(3):
            ops.analyze(1)
        s = dict(zip(sr.SAS_NAMES, ops.eleResponse(1, "material", 1, "sasStats")))
        assert s["updates"] == 0, (opts, newv, s)   # SAS-ME never ran
    _define_all()


def test_no_hold_skip_census_under_sas(capfd):
    """Review item 5: under SAS-ME's paper alpha_in rule no P2-5 reversal test
    is skipped on a hold, so implexGuards[5] must not count one."""
    _quad_deck(W.sas_opts(1e-4)[:5])
    ops.integrator("LoadControl", 0.01)
    ops.analysis("Static")
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in range(3):
        ops.analyze(1)
    g0 = ops.eleResponse(1, "material", 1, "implexGuards")[5]
    ops.integrator("LoadControl", 0.0)
    for _ in range(3):
        ops.analyze(1)
    assert ops.eleResponse(1, "material", 1, "implexGuards")[5] == g0
    _define_all()



# ------------------------------------------------ review of #871 (numerics)
def test_elastic_path_is_exact(protos):
    """Numerics item 1: the elastic predictor is the closed form (sqrt p linear
    in dv), so it lands on the oracle to its own tolerance whatever TolR."""
    worst = 0.0
    for c in _ref_cases("elastic"):
        new, o = W.step(ops, W.TAG_SAS, c, c["dstrain"])
        if c["ref"]["status"] != "ok":
            # the oracle reaches p = 0 (p_floor): SAS-ME must refuse, not invent
            assert o["rc"] != 0
            continue
        assert o["rc"] == 0
        if o["sas"]["elastic"] == 1:
            worst = max(worst, _rel(new["sigma"], c["ref"]["sigma"]))
    assert worst < 1e-8, worst


def test_convergence_with_tolR_against_oracle(protos):
    """Numerics items 1/3/5: the error against the oracle falls with TolR --
    no floor from the elastic part, the intersection or the drift correction."""
    errs = {1e-4: [], 1e-7: []}
    for tol, tag in ((1e-4, W.TAG_SAS), (1e-7, 7)):
        for c in _ref_cases("conv"):
            if c["ref"]["status"] != "ok":
                continue
            new, o = W.step(ops, tag, c, c["dstrain"])
            errs[tol].append(_rel(new["sigma"], c["ref"]["sigma"]) if o["rc"] == 0 else None)
    pairs = [(a, b) for a, b in zip(errs[1e-4], errs[1e-7]) if a is not None and b is not None]
    assert len(pairs) >= 0.8 * len(errs[1e-4])
    b = sorted(x[1] for x in pairs)
    assert b[len(b) // 2] < CONV_MED and b[-1] < CONV_MAX, (b[len(b) // 2], b[-1])
    a = sorted(x[0] for x in pairs)
    assert b[len(b) // 2] < 0.1 * a[len(a) // 2]


CONV_MED, CONV_MAX = 1.0e-7, 5.0e-7   # measured 1.1e-8 / 6.3e-8 at TolR 1e-7


def test_psi_driven_exceedance_is_not_a_dead_end(protos):
    """Numerics item 2 (reviewer's p5.py): proportional elastic compression from
    rho_alpha 0.999 -- psi shrinks the bounding surface around a fixed alpha.
    Every increment integrates (no refusal) and follows the oracle's chain."""
    # the oracle's chain runs to 150 MPa; rho_alpha passes 1 + kappa_entry = 3
    # only beyond ~30 MPa (it is 2.04 at 12.5 MPa), far outside the model's
    # range -- the entry threshold's justification. Integrate to 10 MPa.
    chain = [c for c in _ref_cases("deadend") if c["ref"]["p_end"] <= 1.0e4]
    assert len(chain) >= 50
    st = {k: chain[0][k] for k in ("sigma", "alpha", "alpha_in", "z", "e")}
    prev = 0.0
    maxrho = 0.0
    for c in chain:
        new, o = W.step(ops, W.TAG_SAS, st, c["dstrain"], prev)
        assert o["rc"] == 0, (c["step"], o["sas"]["lastRefuseCode"])
        assert _rel(new["sigma"], c["ref"]["sigma"]) < 1e-7, c["step"]
        maxrho = max(maxrho, W.alpha_over_b(new["sigma"], new["alpha"], new["e"]))
        st, prev = new, W.ncov(c["dstrain"])
    assert maxrho > 1.1 + 1e-3   # the path really does leave 1 + kappa


def test_step_factor_cap(protos):
    """Numerics item 5: accepted substeps grow by at most 1.1."""
    c = next(x for x in _ref_cases("conv") if x["p0"] == 20.0 and x["dir"] == "act" and x["delta"] == 1e-3)
    _, o = W.step(ops, W.TAG_SAS, c, c["dstrain"], trace=100000)
    acc = [r for r in o["trace"] if r["outcome"] == "accept"]
    ratios = [b["dT"] / a["dT"] for a, b in zip(acc, acc[1:]) if b["T"] + b["dT"] < 1.0 - 1e-12]
    assert ratios and max(ratios) <= 1.1 * (1 + 1e-12), max(ratios)


def test_alpha_project_keeps_f_and_lands_on_the_surface(protos):
    """-alphaProject 1 on the inadmissible 1950/3: projected (counted), f stays
    on the cone, rho_alpha is brought to the bounding surface."""
    r = next(x for x in sr.read_ring_csv(sr.RING_CSVS[0]) if x["element"] == 1950 and x["gp"] == 3)
    new, o = W.step(ops, W.TAG_SAS_PRJ, r, [0, 0, 0, 1e-7, 0, 0])
    assert o["rc"] == 0 and o["sas"]["alphaProjected"] >= 1
    assert o["f_after"] <= 1e-6
    assert W.alpha_over_b(new["sigma"], new["alpha"], new["e"]) <= 1.0 + 1e-3


def test_ring_oracle_gap_closed(protos):
    """Numerics item 4: b16 element 5496 (the start-of-increment on-surface
    alpha_in rule) now lands on the oracle."""
    cs = [c for c in _ref_cases("ring") if c["mesh"] == "b16" and c["element"] == 5496]
    worst = 0.0
    for c in cs:
        new, o = W.step(ops, W.TAG_SAS, c, c["dstrain"])
        if o["rc"] == 0 and c["ref"]["status"] == "ok":
            worst = max(worst, _rel(new["sigma"], c["ref"]["sigma"]))
    assert worst < 1e-3, worst

