"""Lemaitre coupled ductile damage — quantitative textbook validation (Zone-A).

This complements the qualitative self-consistency battery in test_lemaitre_damage.py
with the *quantitative* oracle cells of the ADR §6 table
(Ladruno_implementation/15_lemaitre_ductile_damage_adr.md):

  Layer A — uniaxial constitutive oracle (ADR V1)
    A1  LadrunoUniaxialJ2 -damage lemaitre  vs an INDEPENDENT from-the-equations
        backward-Euler reference integrator (_testbed/lemaitre_ref.py): sigma(eps)
        AND D(eps) match to ~1e-8 over a combined-hardening + AF + damage push.
    A1c same, with a load reversal (exercises the AF kinematic + damage together).
    A2  closed form: perfectly plastic + linear damage (s=1, pD=0) =>
        D(eps) = s0^2/(2 E r) * (eps - s0/E), an exact line — checked vs the code.

  Layer B — triaxiality carrier R_v (ADR V3, the multiaxial differentiator)
    B1  recover the effective stress from the 3D brick and verify the WIRED damage
        energy release rate equals the closed form Y = (a q^2 + b sigmaH^2)/(2E),
        a=(2/3)(1+nu), b=3(1-2nu), at BOTH a free-uniaxial (eta~1/3) and a laterally
        confined (eta>1/3) state — i.e. R_v(eta) is exactly the documented function.
    B2  R_v == 1 at uniaxial tension; R_v > 1 (and faster damage) under confinement.

The single-brick drivers reproduce a homogeneous uniaxial-stress state (1/8-symmetry
determinate restraints, free lateral contraction) or a laterally confined state, as
in test_lemaitre_damage.py.
"""
import math

import pytest

from _testbed import ops
from _testbed.lemaitre_ref import Lemaitre1D, integrate_rate

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

_E, _NU = 1000.0, 0.3


def _kg(E, nu):
    return E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))


_K, _G = _kg(_E, _NU)


# ==========================================================================
# 1D truss driver (drives the uniaxial material directly)
# ==========================================================================
def _build_truss(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 1.0)
    mat_fn(1)
    ops.fix(1, 1)
    ops.element("Truss", 1, 1, 2, 1.0, 1)   # area 1, length 1 => strain = disp
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run_truss(mat_fn, eps_hist):
    """Drive an arbitrary monotone-or-reversing strain history (list of total
    strains). Returns list of (eps, sigma, D, p)."""
    _build_truss(mat_fn)
    out = []
    prev = 0.0
    for eps in eps_hist:
        ops.integrator("DisplacementControl", 2, 1, eps - prev)
        if ops.analyze(1) != 0:
            break
        prev = eps
        sig = ops.eleResponse(1, "material", 1, "stress")
        D = ops.eleResponse(1, "material", 1, "damage")
        p = ops.eleResponse(1, "material", 1, "equivalentPlasticStrain")
        out.append((eps, sig[0] if sig else 0.0,
                    D[0] if D else 0.0, p[0] if p else 0.0))
    return out


def _ramp(eps_target, nsteps):
    return [eps_target * i / nsteps for i in range(1, nsteps + 1)]


# shared combined-hardening params (multi-backstress to exercise the full law)
_ISO = ("-iso", "voce", 5.0, 4.0, 30.0, 10.0)        # s0,Qinf,b,Hiso
_KIN = ("-kin", 2, 400.0, 80.0, 150.0, 12.0)          # N=2: (C1,g1),(C2,g2)
_DMG = (1.5, 1.0, 0.01, 0.95)                          # r,s,pD,Dc


def _mat1d(tag, dmg=None):
    args = ["LadrunoUniaxialJ2", tag, _E, *_ISO, *_KIN]
    if dmg:
        args += ["-damage", "lemaitre", *dmg]
    ops.uniaxialMaterial(*args)


def _ref(dmg=None):
    kw = dict(E=_E, s0=5.0, Qinf=4.0, b=30.0, Hiso=10.0,
              C=[400.0, 150.0], g=[80.0, 12.0])
    if dmg:
        kw.update(r=dmg[0], s=dmg[1], pD=dmg[2], Dc=dmg[3])
    return Lemaitre1D(**kw)


# ==========================================================================
# A1 — uniaxial oracle: code vs independent reference integrator (monotone)
# ==========================================================================
def test_A1_oracle_monotone():
    hist = _ramp(0.12, 60)
    code = _run_truss(lambda t: _mat1d(t, dmg=_DMG), hist)
    ref = _ref(_DMG)
    assert len(code) == len(hist)
    maxes = 0.0
    maxd = 0.0
    grew = False
    for (eps, sig, D, _p) in code:
        rs, rD, _rp, _rse = ref.step(eps)
        maxes = max(maxes, abs(sig - rs) / (abs(rs) + 1.0))
        maxd = max(maxd, abs(D - rD))
        grew = grew or D > 1e-3
    assert grew, "damage never accumulated in the validation push"
    assert maxes < 1e-8, f"stress vs reference: max rel err {maxes:.2e}"
    assert maxd < 1e-8, f"damage vs reference: max abs err {maxd:.2e}"


# ==========================================================================
# A3 — CROSS-METHOD oracle: the implicit return map must converge to an INDEPENDENT
#      forward-Euler explicit integration of the same rate equations (different AF
#      discretization, single forward plastic step, no Newton). This catches return-map
#      implementation errors that the same-scheme A1 port cannot. (Physical-equation
#      validation vs an external code, e.g. LS-DYNA *MAT_105, remains deferred.)
# ==========================================================================
def _rate_ref(eps_target, nsub):
    return integrate_rate([eps_target], E=_E, s0=5.0, Qinf=4.0, b=30.0, Hiso=10.0,
                          C=[400.0, 150.0], g=[80.0, 12.0],
                          r=_DMG[0], s=_DMG[1], pD=_DMG[2], Dc=_DMG[3], nsub=nsub)


def test_A3_implicit_converges_to_explicit_rate():
    sig_ref, D_ref = _rate_ref(0.10, 20000)          # fine explicit ≈ exact
    errs = []
    for n in (25, 100, 400):
        res = _run_truss(lambda t: _mat1d(t, dmg=_DMG), _ramp(0.10, n))
        sig, D = res[-1][1], res[-1][2]
        errs.append((abs(sig - sig_ref) / (abs(sig_ref) + 1.0), abs(D - D_ref)))
    assert D_ref > 1e-3, "rate reference accumulated no damage"
    # the implicit material converges to the independent explicit solution as Δε→0
    assert errs[-1][0] < errs[0][0] and errs[-1][1] < errs[0][1], \
        f"not converging to explicit rate ref: stress {[e[0] for e in errs]}, D {[e[1] for e in errs]}"
    assert errs[-1][0] < 0.01, f"refined stress vs explicit rate ref {errs[-1][0]:.2%} (>1%)"
    assert errs[-1][1] < 0.01, f"refined damage vs explicit rate ref {errs[-1][1]:.3f} (>0.01)"


# ==========================================================================
# A1c — same, with a load reversal (AF kinematic + damage together)
# ==========================================================================
def test_A1c_oracle_with_reversal():
    hist = _ramp(0.06, 30) + _ramp(-0.02, 20)[::1]   # push then partial unload/reverse
    # build a strictly sequential history: 0->0.06 then 0.06->-0.02
    hist = [0.06 * i / 30 for i in range(1, 31)] + \
           [0.06 + (-0.08) * j / 40 for j in range(1, 41)]   # reverse down to -0.02
    code = _run_truss(lambda t: _mat1d(t, dmg=_DMG), hist)
    ref = _ref(_DMG)
    maxes = maxd = 0.0
    for (eps, sig, D, _p) in code:
        rs, rD, _rp, _rse = ref.step(eps)
        maxes = max(maxes, abs(sig - rs) / (abs(rs) + 1.0))
        maxd = max(maxd, abs(D - rD))
    assert maxes < 1e-8, f"reversal stress vs reference: {maxes:.2e}"
    assert maxd < 1e-8, f"reversal damage vs reference: {maxd:.2e}"


# ==========================================================================
# A2 — closed form: perfectly plastic + linear damage (s=1, pD=0)
#      sigma~ == s0 on the yield surface => Y == s0^2/(2E) => D(p)= s0^2/(2Er) p,
#      with p = eps - s0/E in the plastic regime.  D(eps) is an exact line.
# ==========================================================================
def test_A2_perfectly_plastic_linear_damage_closed_form():
    s0, r = 5.0, 0.5
    Dc = 0.99

    def mat(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, _E,
                             "-iso", "voce", s0, 0.0, 0.0, 0.0,   # perfectly plastic
                             "-damage", "lemaitre", r, 1.0, 0.0, Dc)

    hist = _ramp(0.05, 100)
    code = _run_truss(mat, hist)
    slope = s0 * s0 / (2.0 * _E * r)        # dD/dp  (== dD/deps in the plastic regime)
    eps_y = s0 / _E
    checked = 0
    for (eps, _sig, D, p) in code:
        if D <= 0.0 or D >= Dc:             # skip elastic + capped tail
            continue
        # closed form in p (robust to the discrete first-plastic step)
        D_pred = slope * p
        assert abs(D - D_pred) <= 1e-6 + 1e-6 * D_pred, \
            f"D(p): code {D:.6e} vs closed form {D_pred:.6e} at p={p:.4e}"
        # and the eps-space line D = slope*(eps - eps_y) once fully plastic
        assert abs(p - (eps - eps_y)) < 1e-6, f"p != eps - s0/E ({p} vs {eps-eps_y})"
        checked += 1
    assert checked > 10, f"too few damaging steps validated ({checked})"


# ==========================================================================
# 3D single-brick driver (uniaxial stress free / laterally confined)
# ==========================================================================
def _build_brick(mat_fn, confine=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *map(float, c))
    mat_fn(1)
    if confine:
        ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 1, 1); ops.fix(4, 1, 1, 1)
        ops.fix(5, 1, 1, 1); ops.fix(6, 0, 1, 1); ops.fix(7, 0, 1, 1); ops.fix(8, 1, 1, 1)
    else:
        ops.fix(1, 1, 1, 1); ops.fix(2, 0, 1, 1); ops.fix(3, 0, 0, 1); ops.fix(4, 1, 0, 1)
        ops.fix(5, 1, 1, 0); ops.fix(6, 0, 1, 0); ops.fix(8, 1, 0, 0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, 0.25, 0.0, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-10, 100, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _run_brick(mat_fn, eps_target, nsteps, confine=False):
    """Returns list of (eps, sig6, D, p) at each step."""
    _build_brick(mat_fn, confine)
    out = []
    d = eps_target / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        if ops.analyze(1) != 0:
            break
        sig = list(ops.eleResponse(1, "stresses"))[0:6]
        D = ops.eleResponse(1, "material", 1, "damage")[0]
        p = ops.eleResponse(1, "material", 1, "equivalentPlasticStrain")[0]
        out.append((ops.nodeDisp(2, 1), sig, D, p))
    return out


def _mat3d(tag, dmg):
    ops.nDMaterial("LadrunoJ2", tag, _K, _G, *_ISO, *_KIN, "-damage", "lemaitre", *dmg)


def _Y_formula(sig6, D, E, nu):
    """Closed-form damage energy release rate from the EFFECTIVE stress recovered
    from the nominal stress vector and the current damage."""
    se = [s / (1.0 - D) for s in sig6]            # effective stress sigma~ = sigma/(1-D)
    sH = (se[0] + se[1] + se[2]) / 3.0
    dev = [se[0] - sH, se[1] - sH, se[2] - sH, se[3], se[4], se[5]]
    J2 = 0.5 * (dev[0]**2 + dev[1]**2 + dev[2]**2) + dev[3]**2 + dev[4]**2 + dev[5]**2
    q = math.sqrt(3.0 * J2)                        # von Mises (effective)
    a = (2.0 / 3.0) * (1.0 + nu)
    b = 3.0 * (1.0 - 2.0 * nu)
    Y = (a * q * q + b * sH * sH) / (2.0 * E)
    eta = sH / q if q > 1e-12 else 0.0
    Rv = a + b * eta * eta
    return Y, q, sH, eta, Rv


# ==========================================================================
# B1 — the WIRED damage energy release rate Y matches the closed form
#      Y = (a q^2 + b sigmaH^2)/(2E) at the actually-realized stress state, for
#      both a free-uniaxial and a confined (triaxial) brick.  Y_implied is backed
#      out of the committed damage increment: dD = (Y/r)^s dp_dmg.
# ==========================================================================
@pytest.mark.parametrize("confine,label", [(False, "free"), (True, "confined")])
def test_B1_Y_release_rate_matches_closed_form(confine, label):
    dmg = (3.0, 1.0, 0.0, 0.98)                    # s=1 => Y_implied = r*dD/dp_dmg
    r = dmg[0]
    res = _run_brick(lambda t: _mat3d(t, dmg), 0.06, 60, confine=confine)
    assert res and res[-1][2] > 1e-3, f"{label}: no damage accumulated"
    checks = 0
    worst = 0.0
    for (e0, _s0, D0, p0), (e1, s1, D1, p1) in zip(res, res[1:]):
        dD = D1 - D0
        dp = p1 - p0
        if dD <= 1e-9 or dp <= 1e-9 or D1 >= dmg[3] - 1e-6:
            continue
        Y_implied = r * (dD / dp)                  # s == 1
        Y_form, q, sH, eta, Rv = _Y_formula(s1, D1, _E, _NU)
        rel = abs(Y_implied - Y_form) / (abs(Y_form) + 1e-12)
        worst = max(worst, rel)
        checks += 1
    assert checks > 5, f"{label}: too few damaging increments ({checks})"
    # backward-Euler dD uses Y at end-of-step; the finite ratio dD/dp introduces a
    # small O(dp) discretization gap vs the pointwise closed form -> few-% tolerance
    assert worst < 0.05, f"{label}: Y(stress) vs closed form max rel err {worst:.2%}"


# ==========================================================================
# B3 — INDEPENDENT uniaxial hand-value: at a free-uniaxial brick the damage energy
#      release rate must reduce to Y = sigma~^2/(2E) EXACTLY and INDEPENDENT of nu
#      (R_v(1/3)==1). Asserting this bare closed form (not the code's a*q^2+b*sigmaH^2
#      form) pins the coefficient identity a + b/9 == 1; a wrong (a,b) pair would pass
#      B1 — which uses the code's own _Y_formula — but fail here. Run at two nu values.
# ==========================================================================
@pytest.mark.parametrize("nu", [0.2, 0.4])
def test_B3_uniaxial_Y_hand_value_nu_independent(nu):
    dmg = (3.0, 1.0, 0.0, 0.98)
    r = dmg[0]
    K, G = _kg(_E, nu)

    def mat(tag):
        ops.nDMaterial("LadrunoJ2", tag, K, G, *_ISO, *_KIN, "-damage", "lemaitre", *dmg)

    res = _run_brick(mat, 0.06, 60, confine=False)
    assert res and res[-1][2] > 1e-3, f"nu={nu}: no damage"
    worst = 0.0
    checks = 0
    for (e0, _s0, D0, p0), (e1, s1, D1, p1) in zip(res, res[1:]):
        dD, dp = D1 - D0, p1 - p0
        if dD <= 1e-9 or dp <= 1e-9 or D1 >= dmg[3] - 1e-6:
            continue
        sigEff_axial = s1[0] / (1.0 - D1)          # uniaxial stress => sigma_xx is the only stress
        Y_hand = sigEff_axial * sigEff_axial / (2.0 * _E)   # bare closed form, nu-FREE
        Y_implied = r * dD / dp                     # s == 1
        worst = max(worst, abs(Y_implied - Y_hand) / (abs(Y_hand) + 1e-12))
        checks += 1
    assert checks > 5, f"nu={nu}: too few damaging increments"
    assert worst < 0.05, \
        f"nu={nu}: wired Y deviates from the nu-independent uniaxial closed form sigma~^2/2E " \
        f"by {worst:.2%} (=> R_v(1/3) != 1, i.e. a+b/9 != 1)"


# ==========================================================================
# B2 — triaxiality signature: R_v == 1 at uniaxial; R_v > 1 confined; the
#      confined (higher-triaxiality) state damages faster at equal axial strain.
# ==========================================================================
def test_B2_triaxiality_Rv_signature():
    dmg = (3.0, 1.0, 0.0, 0.98)
    free = _run_brick(lambda t: _mat3d(t, dmg), 0.06, 60, confine=False)
    conf = _run_brick(lambda t: _mat3d(t, dmg), 0.06, 60, confine=True)

    # R_v at a free-uniaxial damaging step must be ~1 (eta ~ 1/3)
    fz = next(((s, D) for (_e, s, D, _p) in free if D > 1e-3), None)
    assert fz is not None, "free brick never damaged"
    _Yf, _qf, _sHf, eta_f, Rv_f = _Y_formula(fz[0], fz[1], _E, _NU)
    assert abs(eta_f - 1.0 / 3.0) < 0.02, f"free triaxiality {eta_f:.3f} != 1/3"
    assert abs(Rv_f - 1.0) < 0.02, f"R_v at uniaxial {Rv_f:.3f} != 1"

    # confined state has higher triaxiality and R_v > 1
    cz = next(((s, D) for (_e, s, D, _p) in conf if D > 1e-3), None)
    assert cz is not None, "confined brick never damaged"
    _Yc, _qc, _sHc, eta_c, Rv_c = _Y_formula(cz[0], cz[1], _E, _NU)
    assert eta_c > eta_f + 0.05, f"confined triaxiality {eta_c:.3f} !> free {eta_f:.3f}"
    assert Rv_c > 1.0 + 1e-3, f"confined R_v {Rv_c:.3f} !> 1"

    # and damages faster at the same applied axial strain
    assert conf[-1][2] > free[-1][2] + 1e-6, \
        f"confined D {conf[-1][2]:.4f} !> free D {free[-1][2]:.4f}"
