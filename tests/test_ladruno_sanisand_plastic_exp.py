"""WP-229 (TIMs F26): `-sasPlasticExp nh`, a separate pressure exponent for the plastic
modulus of LadrunoSANISAND under SAS-ME (IntScheme 129):

    b0 = G0 h0 (1 - ch e) (p / p_at)^-nh          (DM04: nh = 0.5)

What is pinned (Ladruno_implementation/229_sanisand_plastic_exponent.md):
  (a)  `-sasPlasticExp 0.5` is bit-identical to the flag absent: three undrained
       triaxial chains, R1 off and R1 on, every committed state by float.hex.
       (The flag absent against a build WITHOUT this WP is the 643-replay baseline of
       test_ladruno_sanisand_reseat_r1.py, which this WP leaves green.)
  (b)  nh = 0.25 agrees with the WP-134 oracle carrying the same exponent
       (tests/data/wp229_oracle_nh.json, written by
       Ladruno_files/testbed/sanisand_plastic_exp/gate_nh_oracle.py --fixture) to SAS-ME's
       tolerance, and so does nh = 0.5 -- at p0 = 10, 100 and 1000 kPa;
  (c)  the option DOES something, in the direction the formula says: below p_at a
       smaller nh softens (q lower at the same strain), above p_at it stiffens --
       the negative control for (b);
  (d)  one plastic stage scales as the formula says: at a fixed state the plastic
       corrector moves with nh exactly as Kp ~ (p/p_at)^(0.5 - nh) requires (two exponents
       fix the one unknown, two more are predicted), at 5, 20 and 400 kPa;
  (e)  the parser: nh outside [0, 1] and any scheme but 129 are refused;
  (f)  the exponent is the 13th `sasOptions` value and crosses the datastore wire
       (a skeleton WITHOUT the flag is restored WITH it).

Runtime: a few seconds (C++ replays).
"""
import json
import math
import os

import pytest

from _testbed import ops
import sanisand_replay as sr

_HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(_HERE, "data", "wp229_oracle_nh.json")
R1_ON = ("-sasHFloor", 1.0, "-sasReseatHyst", 1.0, "-sasSoftCap", 0.5)
OPTS = (129, 0, 1, 1.0e-7, 1.0e-4, "-Pmin", 0.0101, "-Presidual", 0.0, "-maxSubsteps", 20000,
        "-flipAlphaIn", "init")


def _fx():
    with open(FIXTURE) as fh:
        return json.load(fh)


def _chain(flags, e0, p0, da, n, tolr=1.0e-4):
    """Undrained triaxial chain of committed replays; returns (q, p, hex rows)."""
    fx = _fx()
    P = list(fx["params"])
    P[2] = float(e0)
    ops.wipe()
    opts = list(OPTS)
    opts[4] = tolr
    ops.nDMaterial("LadrunoSANISAND", 1, *P, *opts, *flags)
    st = dict(sigma=[p0] * 3 + [0.0] * 3, alpha=[0.0] * 6, alpha_in=[0.0] * 6, z=[0.0] * 6, e=float(e0))
    q, p, rows = [], [], []
    for k in range(n):
        o = sr.replay(ops, 1, st["sigma"], st["alpha"], st["alpha_in"], st["z"], st["e"],
                      [da, -0.5 * da, -0.5 * da, 0.0, 0.0, 0.0], "compressionPositive", trace=0)
        assert o["rc"] == 0, (flags, e0, p0, k, o["rc"])
        st = dict(sigma=list(o["sigma"]), alpha=list(o["alpha"]), alpha_in=list(o["alpha_in"]),
                  z=list(o["z"]), e=o["e"])
        q.append(o["q"])
        p.append(o["p"])
        rows.append([float(x).hex() for x in st["sigma"] + st["alpha"] + st["alpha_in"] + st["z"]])
    ops.wipe()
    return q, p, rows


PATHS = [(0.833, 10.0), (0.833, 100.0), (0.735, 1000.0)]


# ----------------------------------------------------------------------- (a)
@pytest.mark.parametrize("r1", [(), R1_ON], ids=["R1off", "R1on"])
@pytest.mark.parametrize("e0,p0", PATHS)
def test_half_is_bit_identical_to_the_flag_absent(r1, e0, p0):
    fx = _fx()
    a = _chain(r1, e0, p0, fx["da"], fx["n"])
    b = _chain(r1 + ("-sasPlasticExp", 0.5), e0, p0, fx["da"], fx["n"])
    assert a[2] == b[2]
    assert [x.hex() for x in a[0]] == [x.hex() for x in b[0]]


# ----------------------------------------------------------------------- (b)
@pytest.mark.parametrize("nh", [0.5, 0.25])
@pytest.mark.parametrize("e0,p0", PATHS)
def test_cxx_matches_the_oracle_with_the_same_exponent(nh, e0, p0):
    fx = _fx()
    ref = fx["paths"][f"e{e0}_p{p0:g}_nh{nh}"]
    q, p, _ = _chain(("-sasPlasticExp", nh), e0, p0, fx["da"], fx["n"])
    qmax = max(ref["q"])
    dq = max(abs(a - b) for a, b in zip(q, ref["q"])) / qmax
    dp = max(abs(a - b) for a, b in zip(p, ref["p"])) / max(ref["p"])
    assert dq < 3.0e-3 and dp < 3.0e-3, (nh, e0, p0, dq, dp)


# ----------------------------------------------------------------------- (c)
def test_the_exponent_moves_the_response_the_way_the_formula_says():
    """b0(nh)/b0(0.5) = (p/p_at)^(0.5 - nh): < 1 below p_at for nh < 0.5 (softer), > 1 above.
    And the oracle's two exponents are far apart there, so (b) cannot pass by accident."""
    fx = _fx()
    lo5, lo25 = fx["paths"]["e0.833_p10_nh0.5"], fx["paths"]["e0.833_p10_nh0.25"]
    hi5, hi25 = fx["paths"]["e0.735_p1000_nh0.5"], fx["paths"]["e0.735_p1000_nh0.25"]
    assert lo25["q"][-1] < 0.9 * lo5["q"][-1]
    assert hi25["q"][-1] > 1.1 * hi5["q"][-1]
    q5, _, _ = _chain((), 0.833, 10.0, fx["da"], fx["n"])
    q25, _, _ = _chain(("-sasPlasticExp", 0.25), 0.833, 10.0, fx["da"], fx["n"])
    assert q25[-1] < 0.9 * q5[-1]
    q5, _, _ = _chain((), 0.735, 1000.0, fx["da"], fx["n"])
    q25, _, _ = _chain(("-sasPlasticExp", 0.25), 0.735, 1000.0, fx["da"], fx["n"])
    assert q25[-1] > 1.1 * q5[-1]


# ----------------------------------------------------------------------- (d)
@pytest.mark.parametrize("p0", [5.0, 20.0, 400.0])
def test_one_stage_scales_by_the_b0_ratio(p0):
    """At ONE state on the yield surface, materials that differ only in nh differ only in
    Kp = (2/3) p h b:n, h = b0 / ((alpha - alpha_in):n).  For a tiny plastic increment
        d_sigma = C:d_eps - lambda C:R,   lambda = N / (Kp + X),   Kp(nh) = Kp(0.5) c(nh),
    with c(nh) = (p/p_at)^(0.5 - nh).  The elastic part and N, X, C:R are common, so
        d_sigma(0.5) - d_sigma(nh)  =  (N / X) C:R  [ 1/(k c(nh) + 1) - 1/(k + 1) ],   k = Kp(0.5)/X.
    Checked: the differences are PARALLEL for every exponent (one direction, C:R), and their
    magnitudes follow the bracket -- k is solved from the 0.75 and 0.25 runs and then has to
    predict the 0.0 and 1.0 runs."""
    fx = _fx()
    P = list(fx["params"])
    e0 = 0.8
    P[2] = e0
    m, patm = P[9], P[8]
    eta = 0.6
    s = eta * p0 * math.sqrt(2.0 / 3.0)
    n = [2.0 / math.sqrt(6.0), -1.0 / math.sqrt(6.0), -1.0 / math.sqrt(6.0), 0.0, 0.0, 0.0]
    sig = [p0 + s * n[0], p0 + s * n[1], p0 + s * n[2], 0.0, 0.0, 0.0]
    al = [(s / p0 - math.sqrt(2.0 / 3.0) * m) * x for x in n]      # f = 0, n = the axial direction
    ain = [0.5 * x for x in al]                                    # (alpha - alpha_in):n > 0: h finite
    de = [1.0e-8, -0.5e-8, -0.5e-8, 0.0, 0.0, 0.0]

    def dsig(x):
        ops.wipe()
        ops.nDMaterial("LadrunoSANISAND", 1, *P, 129, 0, 1, 1.0e-7, 1.0e-4, "-Pmin", 0.0101,
                       "-Presidual", 0.0, "-flipAlphaIn", "init",
                       *(() if x is None else ("-sasPlasticExp", x)))
        o = sr.replay(ops, 1, sig, al, ain, [0.0] * 6, e0, de, "compressionPositive", trace=0)
        ops.wipe()
        assert o["rc"] == 0 and o["path"] == "plastic", (o["rc"], o["path"])
        return [o["sigma"][i] - sig[i] for i in range(3)]

    ref = dsig(None)
    diff = {x: [a_ - b_ for a_, b_ in zip(ref, dsig(x))] for x in (0.75, 0.25, 0.0, 1.0)}
    nrm = {x: math.sqrt(sum(v * v for v in d)) for x, d in diff.items()}
    assert min(nrm.values()) > 1.0e-3 * math.sqrt(sum(v * v for v in ref)), (nrm, ref)   # not vacuous
    u = [v / nrm[0.75] for v in diff[0.75]]
    amp = {}
    for x, d in diff.items():
        amp[x] = sum(a_ * b_ for a_, b_ in zip(d, u))             # signed amplitude along C:R
        perp = math.sqrt(max(nrm[x] ** 2 - amp[x] ** 2, 0.0))
        assert perp < 1.0e-5 * nrm[x], (x, perp, nrm[x])           # parallel

    def c(x):
        return (p0 / patm) ** (0.5 - x)

    def br(k, x):
        return 1.0 / (k * c(x) + 1.0) - 1.0 / (k + 1.0)

    target = amp[0.25] / amp[0.75]
    lo, hi = 1.0e-8, 1.0e8                                         # br(k,.25)/br(k,.75) is monotonic in k
    flo = br(lo, 0.25) / br(lo, 0.75) - target
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        fm = br(mid, 0.25) / br(mid, 0.75) - target
        if fm * flo > 0.0:
            lo, flo = mid, fm
        else:
            hi = mid
    k = math.sqrt(lo * hi)
    assert 1.0e-6 < k < 1.0e6, k
    for x in (0.0, 1.0):
        assert amp[x] / amp[0.75] == pytest.approx(br(k, x) / br(k, 0.75), rel=2.0e-4), (x, k)


# ----------------------------------------------------------------------- (e)
@pytest.mark.parametrize("v", [-0.1, 1.5, float("nan")])
def test_parser_refuses_bad_values(v):
    fx = _fx()
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *fx["params"], *OPTS, "-sasPlasticExp", v)
    ops.wipe()


@pytest.mark.parametrize("scheme", [1, 2, 45])
def test_parser_refuses_the_flag_on_other_schemes(scheme):
    fx = _fx()
    ops.wipe()
    with pytest.raises(Exception):
        ops.nDMaterial("LadrunoSANISAND", 1, *fx["params"], scheme, 0, 1, 1e-7, 1e-4, "-sasPlasticExp", 0.25)
    ops.nDMaterial("LadrunoSANISAND", 1, *fx["params"], *OPTS, "-sasPlasticExp", 0.25)
    ops.wipe()


# ----------------------------------------------------------------------- (f)
def test_the_exponent_crosses_the_datastore_wire():
    import tempfile
    import test_ladruno_sanisand as S
    base = (129, 0, 1, 1.0e-7, 1.0e-4, "-Pmin", 0.0101, "-Presidual", 0.0)
    tag = 229
    with tempfile.TemporaryDirectory(prefix="ladruno_wp229_", ignore_cleanup_errors=True) as td:
        dbpath = os.path.join(td, "nh_rt")
        S._build("LadrunoSANISAND", tag, base + ("-sasPlasticExp", 0.25))
        S._elastic_leg(tag)
        ops.updateMaterialStage("-material", tag, "-stage", 1)
        for step in range(5):
            assert ops.analyze(1) == 0, step
        mid = S._stress()
        saved = list(ops.eleResponse(1, "material", 1, "sasOptions"))
        assert len(saved) == 13 and saved[12] == 0.25, saved
        try:
            ops.database("File", dbpath)
        except Exception as exc:                       # noqa: BLE001
            pytest.skip(f"database() unsupported in this build: {exc}")
        ops.save(1)
        S._build("LadrunoSANISAND", tag, base)          # skeleton WITHOUT the flag
        assert list(ops.eleResponse(1, "material", 1, "sasOptions"))[12] == 0.5
        ops.database("File", dbpath)
        ops.restore(1)
        after = S._stress()
        got = list(ops.eleResponse(1, "material", 1, "sasOptions"))
        # the restored point must also CONTINUE with nh = 0.25
        assert ops.analyze(1) == 0
        cont = S._stress()
        ops.wipe()
    assert got == saved, (saved, got)
    assert S._reldiff(mid, after) <= 1.0e-12, (mid, after)
    S._build("LadrunoSANISAND", tag, base + ("-sasPlasticExp", 0.25))
    S._elastic_leg(tag)
    ops.updateMaterialStage("-material", tag, "-stage", 1)
    for step in range(6):
        assert ops.analyze(1) == 0, step
    ref = S._stress()
    ops.wipe()
    assert S._reldiff(ref, cont) <= 1.0e-12, (ref, cont)
