"""WP-129: SAS-ME (IntScheme 129), the Sloan-Abbo-Sheng-style SANISAND integrator.

What is pinned here (test letters follow the WP-129 brief, T1-T5 WP-128 sec. 6):

  (a)  every EXISTING IntScheme is byte-identical to the unmodified WP-127
       binary (hex floats, wp129_sanisand_byteid.py + its baseline JSON);
  (b)  on benign states SAS-ME agrees with an INDEPENDENT integration of the
       same rate equations -- WP-128's validated Python port with alpha and z
       in its error and the alpha_in bracket (tests/_testbed/sanisand_md_port.py)
       -- at a tight tolerance. RungeKutta45 is NOT a reference (WP-128: dT_min
       1e-3 hard-coded, Mc-clamp force-accept, no drift correction; and its
       stages 3/4 never compute dAlpha3/dAlpha4, so its alpha weights sum to
       301/336). A hook takes WP-134's reference values when they land;
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

import pytest

from _testbed import ops
import wp129_sasme_tools as W
import sanisand_replay as sr

_HERE = os.path.dirname(os.path.abspath(__file__))
KAPPA = W.KAPPA


@pytest.fixture(scope="module")
def protos():
    W.define_prototypes(ops)
    # attribution prototypes need the alpha check out of the way, so that the
    # thing measured is the ablated mechanism, not the backstop
    ops.nDMaterial("LadrunoSANISAND", W.TAG_SAS_ABL, *W.P,
                   *W.sas_opts(1e-4, extra=("-sasAlphaIn", "stale", "-sasErrorVars", "stress",
                                            "-alphaBoundTol", 1.0e6)))
    ops.nDMaterial("LadrunoSANISAND", W.TAG_SAS_NOE, *W.P,
                   *W.sas_opts(1e-4, extra=("-sasErrorVars", "stress", "-alphaBoundTol", 1.0e6)))
    ops.nDMaterial("LadrunoSANISAND", 7, *W.P, *W.sas_opts(1e-8))    # tight, for (b)
    ops.nDMaterial("LadrunoSANISAND", 8, *W.P, *W.sas_opts(1e-10))   # tighter, for (f)
    yield
    ops.wipe()


def _sas_sub(o):
    return int(o["sas"]["substeps"]) if o["sas"] else 0


# ---------------------------------------------------------------------- (a)
def test_existing_schemes_byte_identical():
    import wp129_sanisand_byteid as B
    with open(B.BASELINE) as fh:
        ref = json.load(fh)["decks"]
    cur = B.run_all()
    assert set(cur) == set(ref)
    bad = []
    for name in ref:
        if cur[name] != ref[name]:
            nrow = sum(1 for a, b in zip(cur[name], ref[name]) if a != b)
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
def _port_alpha_aware(tol):
    from _testbed import sanisand_md_port as MP
    return MP.Material(W.P, TolE=tol, alpha_err=True, fabric_err=True,
                       h_mode="macaulay", maxSubsteps=0)


@pytest.mark.parametrize("p0,dname,d", [
    (20.0, "active", [0.3, 1.0, 0.0, 0.0, 0.0, 0.0]),
    (50.0, "active", [0.3, 1.0, 0.0, 0.0, 0.0, 0.0]),
    (50.0, "shear", [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
    (20.0, "passive", [1.0, -0.3, 0.0, 0.0, 0.0, 0.0]),
])
def test_benign_agreement_with_alpha_aware_port(protos, p0, dname, d):
    port = _port_alpha_aware(1e-8)
    st_c = W.k0_state(p0)
    st_p = W.k0_state(p0)
    prev = 0.0
    worst = 0.0
    for k in range(8):
        de = [1e-5 * x for x in d]
        new_c, o = W.step(ops, 7, st_c, de, prev)
        assert o["rc"] == 0
        op = port.update(st_p["sigma"], st_p["alpha"], st_p["alpha_in"], st_p["z"],
                         st_p["e"], de, prev_incr_norm=prev)
        assert op["rc"] == 0
        ds = W.norm([a - b for a, b in zip(new_c["sigma"], op["sigma"])]) / W.norm(op["sigma"])
        da = W.norm([a - b for a, b in zip(new_c["alpha"], op["alpha"])])
        worst = max(worst, ds, da)
        # each integrator continues from its OWN state
        st_c = new_c
        st_p = dict(sigma=op["sigma"], alpha=op["alpha"], alpha_in=op["alpha_in"],
                    z=op["z"], e=op["e"])
        prev = W.ncov(de)
    assert worst < 2e-6, f"SAS-ME vs alpha-aware port: worst rel. difference {worst:.2e}"


WP134_REF = os.path.join(_HERE, "data", "wp134_sanisand_reference.json")


@pytest.mark.skipif(not os.path.exists(WP134_REF),
                    reason="hook: WP-134's independent reference values not delivered yet")
def test_benign_agreement_with_wp134_reference(protos):
    """Hook for WP-134. Expected file layout: a list of cases
    {sigma, alpha, alpha_in, z, e, dstrain, sigma_ref, alpha_ref, tol}
    (internal compression-positive convention)."""
    with open(WP134_REF) as fh:
        cases = json.load(fh)
    for c in cases:
        _, o = W.step(ops, 7, dict(sigma=c["sigma"], alpha=c["alpha"], alpha_in=c["alpha_in"],
                                   z=c["z"], e=c["e"]), c["dstrain"])
        assert o["rc"] == 0
        err = W.norm([a - b for a, b in zip(o["sigma"], c["sigma_ref"])]) / W.norm(c["sigma_ref"])
        assert err < c.get("tol", 1e-5)


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


def test_T1_gate_can_fail_with_both_fixes_ablated(protos):
    rows = _t1(W.TAG_SAS_ABL)
    r = next(r for r in rows if r["ps"] == 0.0101 and r["delta"] == 1e-4)
    assert r["rc"] == 0 and r["ab_n"] > 2.0, r   # the ME escape (5.14) is back


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
    ops.wipe()


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
        h = 1e-7
        de = [h * x for x in dd]
        new, o = W.step(ops, 8, st, de, prev)
        assert o["rc"] == 0 and o["sas"]["elastic"] == 0   # a PLASTIC probe
        fd = [(a - b) / h for a, b in zip(new["sigma"], st["sigma"])]
        an = [sum(Cep[i][j] * dd[j] for j in range(6)) for i in range(6)]
        worst = max(worst, W.norm([a - b for a, b in zip(fd, an)]) / W.norm(an))
    assert worst < 1e-3, worst
