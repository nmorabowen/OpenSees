"""WP-130 (TIMs F18(c) + F18(d)): BackwardEuler_CPPM under a global Newton,
and the ModifiedEuler -> CPPM per-point fallback.

What is pinned here:

1. BYTE IDENTITY at the defaults.  Six IntScheme-2 decks
   (`wp130_sanisand_byteid.py`: vanilla ManzariDafalias, LadrunoSANISAND 3D
   incl. a huge-increment leg that drives the CPPM's halving ladder into its
   explicit fallback, the TIMs campaign set in plane strain incl. a reversal
   and a free-DOF global-Newton deck whose iteration counts are pinned too)
   reproduce, bit for bit, what the PRE-WP-130 binary produced.  WP-127's
   IntScheme-1 decks are pinned by its own test, unchanged.
2. THE PARSER refuses every flag combination that would be inert or is not
   qualified (-implex).
3. THE CENSUS: `substepStats` has 28 columns; at the defaults the vanilla
   silent explicit fallback is now COUNTED (F12 5.2/5.3 said it was invisible).
4. F18(c): with `-cppmOnFail refuse` a trial iterate the CPPM cannot return
   reaches analyze() as a failure (rc < 0) in bounded time, the refusal is
   counted, and nothing is integrated explicitly.
5. F18(d): the one-element fallback test -- a ModifiedEuler cap hit that
   refuses the step without `-meFallback` is carried by the CPPM with it, and
   agrees with an uncapped integration.

MEASURED WALL TIME: see the WP-130 PR (the byte-identity decks dominate).
"""
import json
import math
import sys
import time

import pytest

from _testbed import ops
import test_ladruno_sanisand as sani
import wp127_sanisand_byteid as b127
import wp130_sanisand_byteid as b130

pytestmark = [pytest.mark.zone_a]

_NSTATS = 28
(CPPM_CALLS, CPPM_NFAIL, CPPM_HALV, CPPM_EXPL, CPPM_LOWP, CPPM_REF,
 ME_FB, ME_FB_OK, LAST_CPPM_REF, GUESS_TRIES, GUESS_OK) = range(17, 28)
CAP = 9
SUB = 2


def _stats(ele=1, gp=1):
    s = list(ops.eleResponse(ele, "material", gp, "substepStats"))
    assert len(s) == _NSTATS, s
    return s


# ---------------------------------------------------------------------------
#  1. byte identity against the pre-WP-130 binary
# ---------------------------------------------------------------------------

def _as_float(v):
    """Hex-float strings (float.hex) back to floats; anything else unchanged."""
    if isinstance(v, str) and ("0x" in v or v in ("inf", "-inf", "nan")):
        return float.fromhex(v)
    return v


def _compare(got, ref, what):
    """Exact on win32 (the baselines' platform); elsewhere the fork's 1e-6
    cross-platform floor on floats (test_adr97_p4_inertness.py convention, as
    WP-127: GCC/libm differ from MSVC in the last bits), non-float entries
    (rc, Newton iteration counts) still exact."""
    assert sorted(got) == sorted(ref), what
    for name in ref:
        assert len(got[name]) == len(ref[name]), (what, name)
        if sys.platform == "win32":
            for k, (a, b) in enumerate(zip(got[name], ref[name])):
                assert a == b, f"{what}: deck {name} row {k}: first differing entry " \
                    f"{next(i for i, (x, y) in enumerate(zip(a, b)) if x != y)}"
            continue
        scale = max((abs(_as_float(x)) for row in ref[name] for x in row
                     if isinstance(_as_float(x), float)
                     and math.isfinite(_as_float(x))), default=1.0)
        tol = 1e-6 * max(scale, 1.0)
        for k, (a, b) in enumerate(zip(got[name], ref[name])):
            assert len(a) == len(b), f"{what}: deck {name} row {k}: length"
            for i, (x, y) in enumerate(zip(a, b)):
                fx, fy = _as_float(x), _as_float(y)
                if isinstance(fx, float) and isinstance(fy, float):
                    assert abs(fx - fy) <= tol or (fx != fx and fy != fy), (
                        f"{what}: deck {name} row {k} entry {i}: {fx!r} vs {fy!r} beyond "
                        f"the 1e-6 cross-platform floor ({tol:.3e}) on {sys.platform}")
                else:
                    assert fx == fy, f"{what}: deck {name} row {k} entry {i}: {x!r} vs {y!r}"


def test_optout_reproduces_the_pre_wp130_binary():
    """`-cppmTangent vanilla` (and no option at all on vanilla ManzariDafalias)
    reproduces the PRE-WP-130 binary on every IntScheme-2 deck: the seams, the
    stack-local LU, the de-static NewtonIter and the census move nothing."""
    with open(b130.BASELINE) as fh:
        ref = json.load(fh)["decks"]
    _compare(b130.run_all(("-cppmTangent", "vanilla")), ref, "opt-out vs pre-WP-130")


def test_default_fixed_tangent_regression_pin():
    """RE-BASELINED DELIBERATELY (owner decision, WP-130): `-cppmTangent fixed`
    is the LadrunoSANISAND default, because the vanilla CPPM TanType-2 tangent
    is MINUS the derivative of its own return map (NewtonSol `Cep = -1.0 *
    CSigma`; finite difference ||-T - D_fd||/||D_fd|| = 1.24e-3 against
    ||T - D_fd||/||D_fd|| = 2.0 -- LEDGER_quirks "IntScheme 2's TanType-2
    tangent is MINUS", Ladruno_files/testbed/hypo_bearing/wp130_f18c/
    q_tangent_fd.txt, and the two FD tests below). So the LadrunoSANISAND
    IntScheme-2 + TanType-2 decks change. This pins the new default against
    `wp130_sanisand_byteid_fixed_baseline.json` (written by the WP-130 build)
    and checks the change is EXACTLY the tangent sign: vanilla ManzariDafalias
    is untouched, and on every zero-free-DOF deck the stress / strain / state
    entries equal the pre-WP-130 baseline and only the tangent entries flip
    sign. (The free-DOF decks change throughout: the tangent steers Newton.)"""
    with open(b130.FIXED_BASELINE) as fh:
        ref_fixed = json.load(fh)["decks"]
    with open(b130.BASELINE) as fh:
        ref_pre = json.load(fh)["decks"]
    got = b130.run_all()
    _compare(got, ref_fixed, "default vs fixed baseline")
    if sys.platform != "win32":
        return
    assert got["md3d_s2"] == ref_pre["md3d_s2"]
    ntan = {"ls3d_s2": 36, "ls3d_s2_big": 36, "ls_ps_s2": 9, "ls_ps_s2_cyc": 9}
    for name, n in ntan.items():
        changed = 0
        for a, b in zip(got[name], ref_pre[name]):
            assert a[:-n] == b[:-n], (name, "a non-tangent entry moved")
            ta = [_as_float(x) for x in a[-n:]]
            tb = [_as_float(x) for x in b[-n:]]
            if ta != tb:
                changed += 1
                for x, y in zip(ta, tb):
                    assert x == -y, (name, "not a pure sign flip", x, y)
        assert changed > 0, (name, "no row changed: is the default really fixed?")


# ---------------------------------------------------------------------------
#  2. the parser
# ---------------------------------------------------------------------------

def _mat(*opts):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    try:
        ops.nDMaterial("LadrunoSANISAND", 1, *sani._PARAMS, *opts)
    except Exception:
        return False
    try:
        return 1 in ops.getNDMaterialTags() if hasattr(ops, "getNDMaterialTags") else True
    except Exception:
        return True


@pytest.mark.parametrize("opts", [
    (1, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "refuse"),               # scheme 1: inert
    (2, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "maybe"),                # bad token
    (2, 2, 1, 1e-7, 1e-7, "-cppmHalvings", 10),                   # out of range
    (2, 2, 1, 1e-7, 1e-7, "-cppmHalvings", -1),
    (1, 2, 1, 1e-7, 1e-7, "-cppmLineSearch", "on"),               # scheme 1, no fallback
    (1, 2, 1, 1e-7, 1e-7, "-meFallback", "cppm"),                 # no -maxSubsteps
    (2, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-meFallback", "cppm"),   # scheme 2
    (2, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-implex", "-cppmOnFail", "refuse"),
    (1, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-implex", "-meFallback", "cppm"),
    (1, 2, 1, 1e-7, 1e-7, "-cppmTangent", "fixed"),                # scheme 1, no fallback
    (2, 2, 1, 1e-7, 1e-7, "-cppmTangent", "right"),               # bad token
])
def test_parser_refuses_inert_or_unqualified(opts):
    with pytest.raises(Exception):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        ops.nDMaterial("LadrunoSANISAND", 1, *sani._PARAMS, *opts)


@pytest.mark.parametrize("opts", [
    (2, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "refuse", "-cppmHalvings", 0),
    (2, 2, 1, 1e-7, 1e-7, "-cppmOnFail", "explicit", "-cppmLineSearch", "on"),
    (1, 2, 1, 1e-7, 1e-7, "-maxSubsteps", 50, "-meFallback", "cppm",
     "-cppmHalvings", 3, "-cppmLineSearch", "on"),
])
def test_parser_accepts(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoSANISAND", 1, *sani._PARAMS, *opts)


# ---------------------------------------------------------------------------
#  3. the census makes vanilla's silent explicit fallback visible
# ---------------------------------------------------------------------------

def _free_push(extra, lateral=50.0, push=20.0, tangent="vanilla"):
    """The byte-id free quad (100 kPa, loaded edges) under IntScheme 2 +
    `extra`, one push step of 0.1*push kPa.  Returns (rc, wall s, census).
    `tangent` defaults to the VANILLA (sign-flipped) CPPM tangent on purpose:
    it is what produces F12's off-path trial iterates, i.e. the scenario the
    refusal machinery exists for. The LadrunoSANISAND default is `fixed`."""
    b130.build_free_quad(b130._CAMPAIGN_S2 + ("-cppmTangent", tangent) + tuple(extra),
                         lateral=lateral)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    for j, (x, y) in enumerate(sani._XY):
        if y == 1.:
            ops.load(j + 1, 0.0, -push)
    ops.integrator("LoadControl", 0.1)
    t = time.time()
    rc = ops.analyze(1)
    return rc, time.time() - t, _stats()


def test_default_counts_the_silent_explicit_fallback():
    """F12 5.3: a CPPM failure used to be invisible in every channel.  At the
    DEFAULTS (vanilla control flow, byte-identical) the census now shows it.
    Measured on build 428328adc: rc -3 after 31 Newton iterations in ~5.6 s,
    224 local-Newton failures, 440 half-increments, 4 silent explicit
    fallbacks, 0 refusals."""
    rc, wall, s = _free_push(())
    assert rc < 0
    assert s[CPPM_CALLS] > 0
    assert s[CPPM_NFAIL] > 0 and s[CPPM_HALV] > 0
    assert s[CPPM_EXPL] > 0, s          # the silent fallback, now counted
    assert s[CPPM_REF] == 0 and s[LAST_CPPM_REF] == 0


def test_fixed_tangent_carries_the_same_push():
    """The same deck and step with the LadrunoSANISAND DEFAULT tangent (fixed
    sign) and vanilla control flow otherwise: the global Newton converges."""
    rc, wall, s = _free_push((), tangent="fixed")
    assert rc == 0, (rc, s)


# ---------------------------------------------------------------------------
#  4. F18(c): refuse at once, and the refusal reaches analyze()
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("extra", [
    ("-cppmOnFail", "refuse", "-cppmHalvings", 0),
    ("-cppmOnFail", "refuse", "-cppmHalvings", 0, "-cppmStart", "explicit",
     "-cppmLineSearch", "on"),
])
def test_refusal_reaches_analyze_fast(extra):
    """Same deck and step as above: the refusal propagates (LADRUNO_MATERIAL_
    REFUSED -> quad -> Domain::update -> analyze rc < 0) on a trial iterate
    the CPPM cannot return, with nothing integrated explicitly and no halving.
    Measured 8-22 ms against ~5.6 s at the defaults."""
    rc_def, wall_def, _ = _free_push(())
    rc, wall, s = _free_push(extra)
    assert rc < 0
    assert s[CPPM_REF] == 1 and s[LAST_CPPM_REF] == 1, s
    assert s[CPPM_EXPL] == 0 and s[CPPM_HALV] == 0, s
    assert s[1] == 0, s                 # ModifiedEuler never ran
    assert wall < 1.0 and wall < 0.2 * wall_def, (wall, wall_def)


def test_refused_update_is_not_sticky():
    """The refusal flag is per update: after the failed step a much smaller
    load step from the same committed state integrates (nothing latched)."""
    rc, _, s0 = _free_push(("-cppmOnFail", "refuse", "-cppmHalvings", 0))
    assert rc < 0 and s0[CPPM_REF] == 1
    ops.integrator("LoadControl", 1.0e-4)
    assert ops.analyze(1) == 0
    s = _stats()
    # no new refusal; LAST_CPPM_REF keeps describing the last PLASTIC update
    # (a small step here can be elastic, which by design does not reset it)
    assert s[CPPM_REF] == s0[CPPM_REF], s


# ---------------------------------------------------------------------------
#  5. F18(d): the per-point ModifiedEuler -> CPPM fallback, one element
# ---------------------------------------------------------------------------

def _fallback_deck(extra, econf=3.0e-4, de=2.0e-3, steps=10):
    """Zero-free-DOF plane-strain quad, TIMs campaign set, IntScheme 1,
    confinement to p ~ 131 kPa (econf 3e-4), then `steps` steps of `de`
    axial compression / lateral extension. Returns (rcs, stress, census)."""
    opts = (1, 2, 1, 1.0e-7, 1.0e-7, "-Presidual", 0.0, "-Pmin", 0.0101,
            "-flipAlphaIn", "init") + tuple(extra)
    incs = b127._iso_dev(steps, de, 1.0)
    b127._build_ps(b127._CAMPAIGN, opts, 10, econf, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    rcs = []
    for _ in incs:
        rcs.append(ops.analyze(1))
        if rcs[-1] != 0:
            break
    return rcs, list(ops.eleResponse(1, "material", 1, "stress")), _stats()


def test_me_fallback_carries_a_capped_point():
    """THE one-element test (F10b(b)).  Measured on build 428328adc: uncapped
    ModifiedEuler takes up to 123 substeps in one update on this leg;
    `-maxSubsteps 20` refuses step 1; with `-meFallback cppm` every capped
    update is returned by the CPPM, all 10 steps commit, and the stress agrees
    with the uncapped integration to 1.3 % (the two integrators differ by
    their own discretisation error at this increment, F12 section 2)."""
    rc_ref, sig_ref, s_ref = _fallback_deck(())
    assert rc_ref == [0] * 10 and s_ref[CAP] == 0
    rc_cap, _, s_cap = _fallback_deck(("-maxSubsteps", 20))
    assert rc_cap[-1] < 0 and s_cap[CAP] >= 1 and s_cap[ME_FB] == 0
    rc_fb, sig_fb, s_fb = _fallback_deck(("-maxSubsteps", 20, "-meFallback", "cppm"))
    assert rc_fb == [0] * 10, rc_fb
    assert s_fb[ME_FB] >= 10 and s_fb[ME_FB_OK] == s_fb[ME_FB], s_fb
    assert s_fb[CAP] == s_fb[ME_FB]          # every cap hit was handed over
    assert s_fb[CPPM_REF] == 0 and s_fb[CPPM_EXPL] == 0 and s_fb[LAST_CPPM_REF] == 0
    for a, b in zip(sig_fb[:2], sig_ref[:2]):
        assert abs(a - b) <= 0.03 * abs(b), (sig_fb, sig_ref)


def test_me_fallback_refuses_when_the_cppm_fails_too():
    """The fallback CPPM never integrates explicitly (that would re-enter the
    ModifiedEuler that just failed): where it cannot return the increment the
    update is REFUSED.  Measured: at p ~ 44 kPa and 5e-3 steps the first
    plastic step fails both ways."""
    rcs, _, s = _fallback_deck(("-maxSubsteps", 20, "-meFallback", "cppm", "-cppmHalvings", 0),
                               econf=1.0e-4, de=5.0e-3)
    assert rcs[-1] < 0
    assert s[ME_FB] >= 1 and s[ME_FB_OK] < s[ME_FB] and s[CPPM_REF] >= 1, s
    assert s[CPPM_EXPL] == 0 and s[CPPM_LOWP] == 0, s


# ---------------------------------------------------------------------------
#  6. the CPPM's TanType-2 tangent: sign (vanilla) and -cppmTangent fixed
# ---------------------------------------------------------------------------

def _fd_tangent(extra, h=1.0e-8):
    """3D cube, IntScheme 2, TanType 2, 5 plastic steps; the element's
    `tangent` (= getTangent) after step 6 against a central finite
    difference of the return map from the committed state at step 5
    (ladrunoSANISANDReplay, which reproduces the analysis' step to ~1e-15).
    Returns (||T - D||/||D||, ||-T - D||/||D||, replay error)."""
    import math
    import os
    import sys
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    "Ladruno_scripts"))
    import sanisand_replay as sr

    incs = b127._iso_dev(40, 5.0e-3, 0.5)
    b127._build_3d("LadrunoSANISAND", b127._PARAMS, (2, 2, 1, 1e-7, 1e-7) + tuple(extra),
                   10, 3.0e-6, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)

    def grab():
        g = lambda name: list(ops.eleResponse(1, "material", 1, name))
        return dict(sig=g("stress"), eps=g("strain"), alpha=g("alpha"), ain=g("alpha_in"),
                    z=g("fabric"), e=g("state")[24], tan=g("tangent"))
    for _ in range(4):
        assert ops.analyze(1) == 0
    prev = grab()
    assert ops.analyze(1) == 0
    k = grab()
    assert ops.analyze(1) == 0
    k1 = grab()
    de = [a - b for a, b in zip(k1["eps"], k["eps"])]
    dn = math.sqrt(sum((a - b) ** 2 for a, b in zip(k["eps"][:3], prev["eps"][:3]))
                   + 0.5 * sum((a - b) ** 2 for a, b in zip(k["eps"][3:], prev["eps"][3:])))
    run = lambda d: sr.replay(ops, 1, k["sig"], k["alpha"], k["ain"], k["z"], k["e"], d,
                              "tensionPositive", prev_incr_norm=dn)
    base = run(de)
    rep = max(abs(a - b) for a, b in zip(base["sigma"], k1["sig"])) / max(map(abs, k1["sig"]))
    D = [[0.0] * 6 for _ in range(6)]
    for j in range(6):
        dp = list(de); dp[j] += h
        dm = list(de); dm[j] -= h
        rp, rm = run(dp), run(dm)
        for i in range(6):
            D[i][j] = (rp["sigma"][i] - rm["sigma"][i]) / (2 * h)
    T = [[k1["tan"][6 * i + j] for j in range(6)] for i in range(6)]
    nD = math.sqrt(sum(x * x for r in D for x in r))
    e_plus = math.sqrt(sum((T[i][j] - D[i][j]) ** 2 for i in range(6) for j in range(6))) / nD
    e_minus = math.sqrt(sum((T[i][j] + D[i][j]) ** 2 for i in range(6) for j in range(6))) / nD
    return e_plus, e_minus, rep


def test_vanilla_cppm_tangent_has_the_wrong_sign():
    """Pins the vanilla defect, reachable only through the opt-out: the
    TanType-2 tangent IntScheme 2 hands the element is MINUS the derivative of
    its own return map. Measured: ||-T - D_fd||/||D_fd|| = 1.24e-3 (one
    iterate stale), ||T - D_fd||/||D_fd|| = 2.0."""
    e_plus, e_minus, rep = _fd_tangent(("-cppmTangent", "vanilla"))
    assert rep < 1e-9
    assert e_minus < 1e-2 and e_plus > 1.9, (e_plus, e_minus)


def test_cppm_tangent_default_is_the_algorithmic_tangent():
    """The LadrunoSANISAND DEFAULT (owner decision) is the fixed sign."""
    e_plus, e_minus, rep = _fd_tangent(())
    assert rep < 1e-9
    assert e_plus < 1e-2 and e_minus > 1.9, (e_plus, e_minus)


def _ps_final(extra, s1=0.0, s2=0.0):
    """Zero-free-DOF plane-strain quad (campaign set, IntScheme 2, TanType 2,
    p ~ 131 kPa), 6 steps of 2e-3; the LAST increment's lateral / axial
    displacement increments are perturbed by s1 / s2 (strain units:
    d eps11 = -s1, d eps22 = -s2 on the unit square). Returns the element's
    final (tension-positive) stress and 3x3 tangent."""
    opts = (2, 2, 1, 1e-7, 1e-7, "-Presidual", 0.0, "-Pmin", 0.0101,
            "-flipAlphaIn", "init") + tuple(extra)
    incs = list(b127._iso_dev(6, 1.2e-2, 1.0))
    dl, da = incs[-1]
    incs[-1] = (dl + s1, da + s2)
    b127._build_ps(b127._CAMPAIGN, opts, 10, 3.0e-4, incs)
    ops.updateMaterialStage("-material", 1, "-stage", 0)
    for _ in range(10):
        assert ops.analyze(1) == 0
    ops.updateMaterialStage("-material", 1, "-stage", 1)
    for _ in incs:
        assert ops.analyze(1) == 0
    return (list(ops.eleResponse(1, "material", 1, "stress")),
            list(ops.eleResponse(1, "material", 1, "tangent")))


def _fd_tangent_ps(extra, h=1.0e-8):
    """The same check through the PLANE-STRAIN wrapper (open item 5). sigma_33
    is not exposed by the wrapper's `stress`, so the replay cannot be seeded;
    instead the whole deterministic run is repeated with the last increment
    perturbed, which differentiates the SAME return map from the SAME committed
    state. Compares the 2x2 normal block (e11, e22) of the element's 3x3
    `tangent` with the FD. Returns (e_plus, e_minus)."""
    sig0, tan0 = _ps_final(extra)
    T = [[tan0[3 * i + j] for j in range(2)] for i in range(2)]
    D = [[0.0, 0.0], [0.0, 0.0]]
    for j in range(2):
        sp = _ps_final(extra, *((h, 0.0) if j == 0 else (0.0, h)))[0]
        sm = _ps_final(extra, *((-h, 0.0) if j == 0 else (0.0, -h)))[0]
        for i in range(2):
            D[i][j] = (sp[i] - sm[i]) / (-2.0 * h)     # d eps_jj = -s_j
    nD = math.sqrt(sum(x * x for r in D for x in r))
    e_plus = math.sqrt(sum((T[i][j] - D[i][j]) ** 2 for i in range(2) for j in range(2))) / nD
    e_minus = math.sqrt(sum((T[i][j] + D[i][j]) ** 2 for i in range(2) for j in range(2))) / nD
    return e_plus, e_minus


def test_planestrain_wrapper_tangent_sign():
    """Open item 5 closed by measurement: the plane-strain wrapper hands out
    the same object -- vanilla sign flipped, default fixed."""
    ep_v, em_v = _fd_tangent_ps(("-cppmTangent", "vanilla"))
    ep_f, em_f = _fd_tangent_ps(())
    assert em_v < 1e-2 and ep_v > 1.9, (ep_v, em_v)
    assert ep_f < 1e-2 and em_f > 1.9, (ep_f, em_f)


# ---------------------------------------------------------------------------
#  7. a DISCARDING element: the CPPM refusal must not commit (WP-99 channel)
# ---------------------------------------------------------------------------

def test_cppm_refusal_under_a_discarding_element_does_not_commit():
    """SSPquad DISCARDS the material's setTrialStrain return code, so the
    trial-time refusal is invisible to it and Newton can 'converge' on a
    refused (unintegrated) state. LadrunoSANISAND::commitState's plain path
    now declares the CPPM refusal to Domain::commit() (WP-99's channel) and
    latches: analyze < 0, and the committed strain does not move."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j, (x, y) in enumerate(sani._XY):
        ops.node(j + 1, x, y)
    # the VANILLA tangent: it produces off-path iterates, i.e. refusals (see
    # _free_push); the commit-path check itself does not depend on it
    ops.nDMaterial("LadrunoSANISAND", 1, *b127._CAMPAIGN, *b130._CAMPAIGN_S2,
                   "-cppmTangent", "vanilla", "-cppmOnFail", "refuse", "-cppmHalvings", 0)
    ops.element("SSPquad", 1, 1, 2, 3, 4, 1, "PlaneStrain", 1.0)
    for j, (x, y) in enumerate(sani._XY):
        ops.fix(j + 1, 1 if x == 0. else 0, 1 if y == 0. else 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j, (x, y) in enumerate(sani._XY):
        ops.load(j + 1, -50.0 if x == 1. else 0.0, -50.0 if y == 1. else 0.0)
    ops.constraints("Plain")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.test("NormDispIncr", 1.0e-10, 30, 0)
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
    for j, (x, y) in enumerate(sani._XY):
        if y == 1.:
            ops.load(j + 1, 0.0, -20.0)
    ops.integrator("LoadControl", 0.1)
    rcs, refused = [], False
    for _ in range(5):
        eps_before = list(ops.eleResponse(1, "material", 1, "strain"))
        rc = ops.analyze(1)
        rcs.append(rc)
        if _stats()[CPPM_REF] > 0:
            refused = True
        if rc < 0:
            break
    assert refused, ("the deck never produced a CPPM refusal", rcs)
    assert rcs[-1] < 0, rcs
    # it was the COMMIT that refused (the per-instance commit latch, slot 4 of
    # `implexRefusals`, is set only by a commit-time refusal; -implex is off)
    assert list(ops.eleResponse(1, "material", 1, "implexRefusals"))[4] == 1.0
    # the refused step committed nothing: the material strain is the last
    # committed one, and further steps are refused too (latch), no drift
    assert list(ops.eleResponse(1, "material", 1, "strain")) == eps_before
    assert ops.analyze(1) < 0
    assert list(ops.eleResponse(1, "material", 1, "strain")) == eps_before
