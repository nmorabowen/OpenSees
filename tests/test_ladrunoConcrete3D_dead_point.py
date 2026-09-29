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


# --------------------------------------------------------------------------- #
# BeamFiber view (review M2): the confined-fibre view (driveConfinedFiber) had no dead-point treatment.
# A single unit-area NDFiber in a zeroLengthSection, axial strain imposed by DisplacementControl (lateral block condensed
# against the passive hoop). Fibre-material state is read through the section:
#     eleResponse(1, "section", "fiber", 0, 0, matTag, "damage" | "kappaP" | "returnFailures").
# --------------------------------------------------------------------------- #
def _fiber_build(Gc=5.0, lch=50.0, hoop=1200.0):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    ops.node(1, 0.0, 0.0, 0.0)
    ops.node(2, 0.0, 0.0, 0.0)
    ops.fix(1, 1, 1, 1, 1, 1, 1)
    ops.fix(2, 0, 1, 1, 1, 1, 1)                 # free axial dof only => pure axial fibre strain
    ops.nDMaterial("LadrunoConcrete3D", 1, _E, _NU, _FC, _FT, 0.1, Gc, "-lch", lch, "-hoop", hoop)
    ops.section("NDFiber", 1)
    ops.fiber(0.0, 0.0, 1.0, 1)                  # unit area => P == axial nominal stress
    ops.element("zeroLengthSection", 1, 1, 2, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0, 0, 0, 0, 0, 0)
    ops.system("FullGeneral")                    # the confined tangent is non-symmetric
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-8, 200, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")


def _fiber_state():
    r = lambda k: list(ops.eleResponse(1, "section", "fiber", 0.0, 0.0, 1, k))
    ops.eleResponse(1, "forces")
    return dict(eps=ops.nodeDisp(2, 1), P=float(ops.eleResponse(1, "section", "force")[0]),
                dmg=r("damage"), kp=r("kappaP")[0], fails=r("returnFailures")[0])


def _fiber_step(deps):
    ops.integrator("DisplacementControl", 2, 1, deps)
    assert ops.analyze(1) == 0, f"fibre step refused/failed at eps={ops.nodeDisp(2, 1):.5f} (deps {deps})"
    return _fiber_state()


def test_beamfiber_view_tension_death_then_parallel_compression():
    """Tension past omega_t = 0.998 (no kappa_p runaway, no refusal, nominal at the floor), unloading, then compression
    PARALLEL to the crack under the hoop: the strut keeps its plasticity, hardening and omega_c."""
    _fiber_build()
    st = _fiber_step(0.0)
    for _ in range(400):                          # to the death of the tensile channel
        st = _fiber_step(5.0e-5)
        if st["dmg"][0] >= 0.998:
            break
    assert st["dmg"][0] >= 0.998, "the fibre never reached omega_t >= 0.998"
    kp_dead, eps_dead = st["kp"], st["eps"]
    assert kp_dead > 1.0                           # a hardened cracked point: the regime that used to run away
    while st["eps"] < 10.0 * eps_dead:             # 10x further in tension
        st = _fiber_step(2.5e-4)
        assert st["kp"] == pytest.approx(kp_dead, rel=1e-12, abs=1e-12), "kappa_p moved on a tension-dead fibre"
        assert st["fails"] == 0.0
    assert abs(st["P"]) < 1.0e-3 * _FT             # nominal at the floor level, not a spurious fraction of ft
    assert st["dmg"][0] >= 1.0 - 1.0e-5
    while st["eps"] > 0.0:                         # unload to zero strain: no tension is carried
        st = _fiber_step(-1.0e-3)
        assert st["kp"] == pytest.approx(kp_dead, rel=1e-12, abs=1e-12) and st["fails"] == 0.0
        if st["eps"] > 0.0:
            assert st["P"] <= 1.0e-3 * _FT
    Pmin = 0.0
    for _ in range(200):                           # compression parallel to the crack
        st = _fiber_step(-2.0e-5)
        Pmin = min(Pmin, st["P"])
        assert st["fails"] == 0.0
        if st["eps"] < -3.0e-3:
            break
    # Measured (this material, hoop 1200, lch 50; strain to -3e-3): virgin fibre P_min = -39.8 (confined, still rising); cracked
    # fibre WITH the treatment -16.98 at eps = -9.2e-4, then softening (omega_c 0.94); cracked fibre with the treatment disabled
    # (-deadThreshold 2) only -5.63 (kappa_p runs to 1.5e4 and the tensile flow poisons the strut). The strut is a genuine, much
    # weaker-than-virgin cracked-concrete strut, 3x the legacy one; bound 0.5 fc.
    assert Pmin < -0.5 * _FC, f"the cracked fibre's compressive strut collapsed (P_min = {Pmin:.2f})"
    assert st["dmg"][1] > 0.0                      # omega_c evolves normally on the live compressive channel


def test_beamfiber_view_crushing_freezes_the_fibre():
    """Compression past omega_c = 0.998 (steep compressive law): kappa_p frozen, both damages at the floor, nominal at the
    floor level, no refusal; unloading elastic."""
    _fiber_build(Gc=0.075, lch=50.0, hoop=0.0)
    st = _fiber_step(0.0)
    for _ in range(600):
        st = _fiber_step(-5.0e-5)
        if st["dmg"][1] >= 0.998:
            break
    assert st["dmg"][1] >= 0.998, "the fibre never reached omega_c >= 0.998"
    kp_dead, eps_dead = st["kp"], st["eps"]
    for _ in range(40):
        st = _fiber_step(-5.0e-4)
        assert st["kp"] == pytest.approx(kp_dead, rel=1e-12, abs=1e-12) and st["fails"] == 0.0
    assert abs(st["P"]) < 1.0e-3 * _FC
    assert st["dmg"][0] >= 1.0 - 1.0e-5 and st["dmg"][1] >= 1.0 - 1.0e-5
    while st["eps"] < eps_dead:                    # unload (elastic, still the floor-level nominal)
        st = _fiber_step(1.0e-3)
        assert st["kp"] == pytest.approx(kp_dead, rel=1e-12, abs=1e-12) and st["fails"] == 0.0
        assert abs(st["P"]) < 1.0e-3 * _FC
