"""WP concrete3d-hang-diagnosis #877 follow-up: DEAD-POINT treatment at the material point, through the real nDMaterial
(kernel + wrapper + element), not the numpy oracle.

A point whose committed tensile damage has reached ``omegaDead`` (default 0.998, a residual strength fraction of 2e-3) is a
fully open crack: its tensile effective stress is carried ELASTICALLY (no plastic flow, so kappa_p stops growing) and the
return map runs on the compressive remainder only. Before this, kappa_p and sig_eff ran away on such a point (K&R coarse
kappa_p 3.4e4 / sig_eff 813 MPa; the G5 band 3.2e4 / 779 MPa) until the return map could not integrate it and the analysis
aborted. The oracle-level statement of the same contract is test_dead_point_tension_cutoff_gate in
test_ladrunoConcrete3D_material.py; this one drives a unit cube in UNIFORM uniaxial-strain tension so every Gauss point is
the same material point, well past the point of death, and then unloads.

Also pins the ``-deadThreshold`` flag (0.998 default; >= 1 disables the treatment and restores the runaway).
"""
import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

_NODES = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
          5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}
_E, _NU, _FC, _FT = 30000.0, 0.2, 30.0, 3.0
_EPS_UNIT = 1.0e-3            # load factor 1 == exx of 1e-3


def _build(threshold_args=()):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for n, c in _NODES.items():
        ops.node(n, *c)
    ops.nDMaterial("LadrunoConcrete3D", 1, _E, _NU, _FC, _FT, 0.1, 5.0, "-lch", 50.0, *threshold_args)
    for n in _NODES:                               # uniaxial STRAIN: uy = uz = 0 everywhere; x = 0 face also ux = 0
        ops.fix(n, 1 if n in (1, 4, 5, 8) else 0, 1, 1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.sp(n, 1, _EPS_UNIT)
    ops.system("UmfPack")                          # the tangent is non-symmetric
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-12, 20)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _point():
    r = lambda k: list(ops.eleResponse(1, "material", 1, k))
    return dict(strain=r("strain"), sig=r("stress"), sigEff=r("effectiveStress"), kp=r("kappaP")[0],
                dmg=r("damage"), fails=r("returnFailures")[0])


def _ramp(target, dlam, stop=None):
    """LoadControl towards `target` (load factor) in steps of dlam; `stop(point)` ends the ramp early. Returns the point
    list [(lam, point)]."""
    ops.integrator("LoadControl", dlam)
    out = []
    lam = ops.getTime()
    while (lam < target - 1e-12) if dlam > 0 else (lam > target + 1e-12):
        assert ops.analyze(1) == 0, f"analyze failed at load factor {lam:.4f}"
        lam = ops.getTime()
        p = _point()
        out.append((lam, p))
        if stop is not None and stop(p):
            break
    return out


def _run_to_ten_times_death(threshold_args=()):
    _build(threshold_args)
    up = _ramp(1.0e3, 0.02, stop=lambda p: p["dmg"][0] >= 0.998)
    lam_dead, p_dead = up[-1]
    assert p_dead["dmg"][0] >= 0.998, "the point never died"
    more = _ramp(10.0 * lam_dead, 0.1)
    return lam_dead, p_dead, more


def test_dead_point_default_freezes_tension_channel_and_unloads_elastically():
    lam_dead, p_dead, more = _run_to_ten_times_death()
    ctx = 1.0 - _NU
    e_prime = _E * ctx / ((1.0 + _NU) * (1.0 - 2.0 * _NU))          # uniaxial-strain (oedometer) stiffness
    exx_dead = p_dead["strain"][0]
    epsp_dead = exx_dead - p_dead["sigEff"][0] / e_prime                # frozen plastic strain (axial, uniaxial strain)
    assert p_dead["kp"] > 1.0                                           # a hardened, cracked point (the runaway regime)
    # (1) kappa_p frozen, zero return-map failures, sig_eff bounded by the elastic response on the frozen plastic strain
    for lam, p in more:
        assert p["kp"] == pytest.approx(p_dead["kp"], rel=1e-12, abs=1e-12), f"kappa_p moved on a dead point (lam {lam})"
        assert p["fails"] == 0.0
        exx = p["strain"][0]
        assert abs(p["sigEff"][0]) <= e_prime * abs(exx - epsp_dead) * (1.0 + 1e-6) + 1e-9
    last = more[-1][1]
    assert abs(last["sigEff"][0]) > 100.0 * _FT                          # the old runaway regime (sig_eff >> ft) ...
    # (2) ... where the nominal stress is at the FLOOR level, not a spurious fraction of ft
    assert abs(last["sig"][0]) <= 2.0e-6 * abs(last["sigEff"][0]) + 1e-12
    assert abs(last["sig"][0]) < 1.0e-3 * _FT
    assert last["dmg"][0] >= 1.0 - 1.0e-5
    # (3) unload to zero strain: no tension is ever carried; kappa_p stays frozen; zero failures
    lam_top = ops.getTime()
    down = _ramp(0.0, -0.5)
    assert ops.getTime() == pytest.approx(0.0, abs=1e-9)
    for lam, p in down:
        assert p["kp"] == pytest.approx(p_dead["kp"], rel=1e-12, abs=1e-12)
        assert p["fails"] == 0.0
        if p["strain"][0] > epsp_dead:
            assert max(p["sig"][:3]) <= 1.0e-3 * _FT                     # nominal tension stays at the floor
    # sig_eff returns to zero at eps = eps_p and is the elastic crack-closure compression E'(0 - eps_p) at zero strain
    zero = down[-1][1]
    assert zero["sigEff"][0] == pytest.approx(-e_prime * epsp_dead, rel=1e-6, abs=1e-6)
    assert lam_top > lam_dead


def test_dead_threshold_flag_disables_the_treatment():
    """-deadThreshold 2 (>= 1) restores the pre-fix behaviour on the SAME path: kappa_p keeps running on the dead point.
    Pins that the flag is parsed, carried by the copy the element runs, and honoured by the kernel."""
    _build(("-deadThreshold", 2.0))
    up = _ramp(1.0e3, 0.02, stop=lambda p: p["dmg"][0] >= 0.998)
    lam_dead, p_dead = up[-1]
    more = _ramp(10.0 * lam_dead, 0.1)
    assert more[-1][1]["kp"] > 1.5 * p_dead["kp"], "kappa_p should keep growing with the treatment disabled"
