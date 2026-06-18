"""LadrunoRCConcrete Phase-3 tension stiffening (``-tensStiff {vc|cm}``, classTag
33015) — Zone-A battery, the material AS INTEGRATED IN OPENSEES.

Tension stiffening is a stress FLOOR on the principal tensile axis: between cracks
bonded reinforcement holds the average concrete tension above the bare
fracture-energy softening curve,
    sigma_ts = ft / (1 + sqrt(c * e1))             (vc, MCFT / Bentz)
    sigma_ts = alpha * ft / (1 + sqrt(500 * e1))   (cm, Collins-Mitchell)
active ONLY post-crack (e1 >= eps_cr); default off => baseline-identical. The
closed-form gate is proven first in the numpy oracle (tests/_testbed/rc_shell_ref.py
run_T1) and a standalone g++ build of the kernel; here we verify the C++ material
end-to-end on a uniaxial-stress ``stdBrick`` where the live principal tensile axis
is x, so e1 == eps_xx and the pinned stress sigma_xx == sigma_ts(eps_xx) exactly.
"""
import pytest

from _testbed import ops
from _testbed.rc_shell_ref import tens_stiff
from _testbed.roundtrip import database_roundtrip

pytestmark = [pytest.mark.zone_a]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}

E, NU, KC = 30000.0, 0.2, 2.0 / 3.0
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]
FT = max(TS)            # tension-backbone peak ft = 3
EPS_CR = FT / E         # 1e-4
C_DEFAULT = 500.0


def _rc(tag, tens_stiff_mode=None, c=None, alpha=None):
    args = ["LadrunoRCConcrete", tag, E, NU,
            "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if tens_stiff_mode is not None:
        args += ["-tensStiff", tens_stiff_mode]
        if c is not None:
            args += ["-tensStiffC", c]
        if alpha is not None:
            args += ["-tensStiffAlpha", alpha]
    ops.nDMaterial(*args)


def _build(mat_fn):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, co in _CUBE.items():
        ops.node(t, *co)
    mat_fn(1)
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 1)
    ops.fix(3, 0, 0, 1)
    ops.fix(4, 1, 0, 1)
    ops.fix(5, 1, 1, 0)
    ops.fix(6, 0, 1, 0)
    ops.fix(8, 1, 0, 0)
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


def _run(mat_fn, eps_target, nsteps):
    """Monotonic uniaxial-stress tension ramp of eps_xx. Returns [(eps_xx, sig_xx)]."""
    _build(mat_fn)
    out = []
    d = eps_target / nsteps
    ops.integrator("DisplacementControl", 2, 1, d)
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, f"analyze failed heading to eps={eps_target}"
        ops.eleResponse(1, "forces")
        sig = list(ops.eleResponse(1, "stresses"))[0:6]
        out.append((ops.nodeDisp(2, 1), sig[0]))
    return out


def _sig_ts(e1, mode=1, c=C_DEFAULT, alpha=1.0):
    s, _ = tens_stiff(e1, mode, FT, c, alpha)
    return s


# --------------------------------------------------------------------------
def test_pre_crack_untouched():
    """For eps_xx < eps_cr the elastic branch is identical with TS on vs off."""
    off = _run(lambda t: _rc(t), 8.0e-5, 8)
    vc = _run(lambda t: _rc(t, "vc"), 8.0e-5, 8)
    for (e0, s0), (e1, s1) in zip(off, vc):
        assert abs(s1 - s0) <= 1.0e-7 * (abs(s0) + 1.0), (e0, s0, s1)


def test_floor_closed_form():
    """Post-crack, where the floor binds, sigma_xx == sigma_ts(eps_xx) exactly,
    and never drops below the bare softening curve."""
    off = _run(lambda t: _rc(t), 4.0e-3, 200)
    vc = _run(lambda t: _rc(t, "vc"), 4.0e-3, 200)
    nbind = 0
    for (e, s_off), (_, s_on) in zip(off, vc):
        if e < EPS_CR:
            continue
        assert s_on >= s_off - 1.0e-6, (e, s_off, s_on)   # never below bare
        sts = _sig_ts(e)
        if sts > s_off + 1.0e-6:                           # floor binds here
            nbind += 1
            assert abs(s_on - sts) <= 5.0e-3 * (abs(sts) + 1.0), (e, s_on, sts)
    assert nbind > 20, f"floor never exercised ({nbind} binding steps)"


def test_above_bare_softening():
    """Far post-crack the bare curve has softened toward ~0 while TS holds a
    finite floor: TS-on sigma_xx is materially above TS-off."""
    e_far = 3.0e-3
    off = _run(lambda t: _rc(t), e_far, 150)[-1]
    vc = _run(lambda t: _rc(t, "vc"), e_far, 150)[-1]
    assert vc[1] > off[1] + 0.1, (off, vc)
    assert abs(vc[1] - _sig_ts(e_far)) <= 5.0e-3 * (_sig_ts(e_far) + 1.0), vc


def test_vc_cm_default_equivalence():
    """vc default (c=500) and cm default (alpha=1, fixed 500) are the same curve."""
    vc = _run(lambda t: _rc(t, "vc"), 3.0e-3, 120)
    cm = _run(lambda t: _rc(t, "cm"), 3.0e-3, 120)
    for (_, s_vc), (_, s_cm) in zip(vc, cm):
        assert abs(s_vc - s_cm) <= 1.0e-6 * (abs(s_vc) + 1.0), (s_vc, s_cm)


def test_c_knob_changes_floor():
    """A larger c lowers the tension-stiffening floor (steeper decay)."""
    e_far = 3.0e-3
    lo = _run(lambda t: _rc(t, "vc", c=200.0), e_far, 150)[-1]
    hi = _run(lambda t: _rc(t, "vc", c=2000.0), e_far, 150)[-1]
    assert lo[1] > hi[1] + 0.05, (lo, hi)
    assert abs(lo[1] - _sig_ts(e_far, c=200.0)) <= 5.0e-3 * (_sig_ts(e_far, c=200.0) + 1.0)
    assert abs(hi[1] - _sig_ts(e_far, c=2000.0)) <= 5.0e-3 * (_sig_ts(e_far, c=2000.0) + 1.0)


def _homog(mat_fn, exx, eyy, ezz, nsteps=120):
    """Impose a HOMOGENEOUS strain (no shear) on the unit cube by prescribing every
    nodal disp u=(exx*x, eyy*y, ezz*z). Returns the centroid-GP stress (6-vec)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, co in _CUBE.items():
        ops.node(t, *co)
    mat_fn(1)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)
        ops.sp(t, 2, eyy * y)
        ops.sp(t, 3, ezz * z)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-8, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for _ in range(nsteps):
        assert ops.analyze(1) == 0, "homogeneous-strain analyze failed"
        ops.eleResponse(1, "forces")
    return list(ops.eleResponse(1, "stresses"))[0:6]


def test_equibiaxial_floor_both_normals():
    """DEGENERATE (equibiaxial) state: with exx==eyy the principal axis is ill-defined;
    the floor must pin BOTH in-plane normals to sigma_ts(e1) (e1 == exx). The earlier
    rank-1 (0.5,0.5) reuse reached only halfway (self-consistency coeff 0.5)."""
    e = 1.5e-3
    sig = _homog(lambda t: _rc(t, "vc"), e, e, 0.0)
    sts = _sig_ts(e)
    assert abs(sig[0] - sts) <= 1.0e-2 * (sts + 1.0), (sig[0], sts)
    assert abs(sig[1] - sts) <= 1.0e-2 * (sts + 1.0), (sig[1], sts)


# --------------------------------------------------------------------------
# PlateFiber / shell view — the real deployment path (ASDShellQ4 + LayeredShell,
# getCopy("PlateFiber") + the sigma_33=0 inner Newton + condensation).
_QUAD = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}
_H = 0.1
_NL = 4


def _build_shell(tens_stiff_mode=None, interlock=False):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for t, (x, y) in _QUAD.items():
        ops.node(t, x, y, 0.0)
    args = ["LadrunoRCConcrete", 1, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
            "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC]
    if tens_stiff_mode is not None:
        args += ["-tensStiff", tens_stiff_mode]
    if interlock:
        args += ["-interlock"]
    ops.nDMaterial(*args)
    layers = []
    for _ in range(_NL):
        layers += [1, _H / _NL]
    ops.section("LayeredShell", 1, _NL, *layers)
    ops.element("ASDShellQ4", 1, 1, 2, 3, 4, 1)


def _shell_tension(tens_stiff_mode, exx_max, nsteps, interlock=False):
    _build_shell(tens_stiff_mode, interlock)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for t, (x, y) in _QUAD.items():
        ops.sp(t, 1, exx_max * x)
        ops.sp(t, 2, 0.0); ops.sp(t, 3, 0.0)
        ops.sp(t, 4, 0.0); ops.sp(t, 5, 0.0); ops.sp(t, 6, 0.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e12, 1.0e12)
    ops.test("NormDispIncr", 1.0e-7, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    nxx = []
    for k in range(nsteps):
        assert ops.analyze(1) == 0, "shell tension analyze failed"
        ops.eleResponse(1, "forces")
        nxx.append((exx_max * (k + 1) / nsteps, list(ops.eleResponse(1, "stresses"))[0]))
    return nxx


def test_plate_fiber_shell_floor():
    """End-to-end in ASDShellQ4 + LayeredShell (PlateFiber view): membrane tension
    resultant Nxx floors at sigma_ts(eps_xx)*h with TS on, vs softening toward ~0 with
    TS off. Exercises the sigma_33=0 inner Newton + condensation with the TS term."""
    off = _shell_tension(None, 2.5e-3, 120)
    on = _shell_tension("vc", 2.5e-3, 120)
    e_end, nxx_off = off[-1]
    _, nxx_on = on[-1]
    assert nxx_on > 1.5 * nxx_off, (nxx_off, nxx_on)            # floor materially above bare
    assert nxx_on == pytest.approx(_sig_ts(e_end) * _H, rel=0.08), (nxx_on, _sig_ts(e_end) * _H)


def test_with_interlock_monotonic_converges():
    """TS + fixed-crack interlock under PROPORTIONAL (non-rotating) membrane tension:
    the combined run converges and the membrane resultant still floors at sigma_ts*h
    (interlock only bounds crack-plane SHEAR, ~0 under proportional tension)."""
    on = _shell_tension("vc", 2.5e-3, 120, interlock=True)
    e_end, nxx_on = on[-1]
    assert nxx_on == pytest.approx(_sig_ts(e_end) * _H, rel=0.10), (nxx_on, _sig_ts(e_end) * _H)


def test_cm_alpha_curve():
    """cm mode with alpha != 1 scales the floor by alpha (vs the oracle)."""
    e_far = 3.0e-3
    sig = _run(lambda t: _rc(t, "cm", alpha=0.7), e_far, 150)[-1]
    target = _sig_ts(e_far, mode=2, alpha=0.7)
    assert abs(sig[1] - target) <= 5.0e-3 * (target + 1.0), (sig, target)


def test_load_unload_floor_tracks_live_e1():
    """TS v1 is a MONOTONIC backbone floor on the LIVE e1 (no eps1max memory). The analysis
    must stay robust through a load->unload, and the floored stress tracks sigma_ts(live e1)
    throughout — including RE-INFLATING on unload (sigma_ts grows as e1 shrinks). This pins
    the documented monotonic-scope behavior (not a hysteretic cyclic model)."""
    _build(lambda t: _rc(t, "vc"))
    ops.integrator("DisplacementControl", 2, 1, 2.0e-3 / 50)
    for _ in range(50):                       # load to eps=2e-3 (floor binding)
        assert ops.analyze(1) == 0
    ops.eleResponse(1, "forces")
    s_peak = list(ops.eleResponse(1, "stresses"))[0]
    assert s_peak == pytest.approx(_sig_ts(2.0e-3), rel=0.05), (s_peak, _sig_ts(2.0e-3))
    ops.integrator("DisplacementControl", 2, 1, -1.5e-3 / 50)
    for _ in range(50):                       # unload to eps=0.5e-3
        assert ops.analyze(1) == 0
    ops.eleResponse(1, "forces")
    e_lo = ops.nodeDisp(2, 1)
    s_lo = list(ops.eleResponse(1, "stresses"))[0]
    # documented v1 behavior: floor tracks the LIVE e1, so it re-inflates on unload
    assert s_lo == pytest.approx(_sig_ts(e_lo), rel=0.05), (s_lo, _sig_ts(e_lo))
    assert s_lo > s_peak, (s_peak, s_lo)      # re-inflated (sigma_ts decreasing in e1)


def test_parser_rejects_bad_args():
    """The new validation paths must reject: -tensStiffC<=0 (sqrt-NaN guard) and an unknown
    -tensStiff mode. A dropped guard would re-introduce post-crack NaN and stay green."""
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    base = ["LadrunoRCConcrete", 1, E, NU, "-Ce", *CE, "-Cs", *CS, "-Te", *TE, "-Ts", *TS]
    with pytest.raises(Exception):
        ops.nDMaterial(*base, "-tensStiff", "vc", "-tensStiffC", 0.0)
    with pytest.raises(Exception):
        ops.nDMaterial(*[*base[:1], 2, *base[2:]], "-tensStiff", "vc", "-tensStiffC", -5.0)
    with pytest.raises(Exception):
        ops.nDMaterial(*[*base[:1], 3, *base[2:]], "-tensStiff", "bogus")


def test_serialization_roundtrip():
    """schema v4 send/recv balance: drive a TS-on material past cracking, round-trip
    the committed state through the FE_Datastore, and confirm the pinned stress
    survives bit-for-bit."""
    def build():
        _build(lambda t: _rc(t, "vc"))
        ops.integrator("DisplacementControl", 2, 1, 1.0e-3)
        assert ops.analyze(3) == 0
        ops.eleResponse(1, "forces")

    database_roundtrip(
        build, probe_nodes=[2], ndf=3,
        probe_fn=lambda: list(ops.eleResponse(1, "stresses"))[0:6],
    )


def test_serialization_nondefault_params():
    """The v4 wire must carry tensStiffC / tensStiffAlpha / ftPeak, not just survive at defaults.
    Round-trip a NON-default vc c=2000 (and a cm alpha=0.6): if a slot is dropped/misordered,
    recvSelf reverts to the ctor default (c=500 / alpha=1) and the restored stress diverges from
    the rebuilt parser value, so before!=after fails. The probe targets the floored stress, which
    depends on those values."""
    for mode, kw in (("vc", {"c": 2000.0}), ("cm", {"alpha": 0.6})):
        def build(mode=mode, kw=kw):
            _build(lambda t: _rc(t, mode, **kw))
            ops.integrator("DisplacementControl", 2, 1, 1.0e-3)
            assert ops.analyze(3) == 0
            ops.eleResponse(1, "forces")
        database_roundtrip(
            build, probe_nodes=[2], ndf=3, dbname=f"rc_ts_{mode}",
            probe_fn=lambda: list(ops.eleResponse(1, "stresses"))[0:6],
        )


def test_implex_with_tensstiff():
    """TS under IMPL-EX (`-implex`): the analysis converges and the floor is still applied —
    the floored stress tracks sigma_ts(ε1) (within the one-step extrapolation lag) and stays
    materially above the bare-off softening run."""
    def _rc_implex(t):
        ops.nDMaterial("LadrunoRCConcrete", t, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                       "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC, "-tensStiff", "vc", "-implex")
    _build(_rc_implex)
    ops.integrator("DisplacementControl", 2, 1, 3.0e-3 / 150)
    for _ in range(150):
        assert ops.analyze(1) == 0, "implex+TS analyze failed"
    ops.eleResponse(1, "forces")
    e = ops.nodeDisp(2, 1)
    s_on = list(ops.eleResponse(1, "stresses"))[0]
    s_off = _run(lambda t: _rc(t), 3.0e-3, 150)[-1][1]
    assert s_on > 1.5 * s_off, (s_off, s_on)                       # floor applied under implex
    assert s_on == pytest.approx(_sig_ts(e), rel=0.06), (s_on, _sig_ts(e))


def test_tensstiff_interlock_nonproportional_robust():
    """TS (live p1) + fixed-crack interlock (frozen crack) under a NON-proportional path
    (tension to crack, then add shear so the principal axis rotates). This drives the documented
    TS→interlock boundary: the interlock (frozen crack) reworks the in-plane stress AFTER the TS
    floor (live p1), so the floor is NOT cleanly preserved here — exactly why combined TS+interlock
    is scoped to PROPORTIONAL loading. The robustness claim is what matters: the solve must stay
    well-behaved (every step converges, finite stress) through the boundary."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, co in _CUBE.items():
        ops.node(t, *co)
    ops.nDMaterial("LadrunoRCConcrete", 1, E, NU, "-Ce", *CE, "-Cs", *CS, "-Cd", *CD,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD, "-Kc", KC,
                   "-tensStiff", "vc", "-interlock", "-crackSpacing", 50.0)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    exx, gmax = 2.0e-3, 1.0e-3
    times = [0.0, 1.0, 2.0]
    ops.timeSeries("Path", 1, "-time", *times, "-values", 0.0, 1.0, 1.0)   # tension then hold
    ops.timeSeries("Path", 2, "-time", *times, "-values", 0.0, 0.0, 1.0)   # shear in stage 2
    ops.pattern("Plain", 1, 1)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 1, exx * x)                                   # u_x = exx*X
        ops.sp(t, 3, 0.0)
    ops.pattern("Plain", 2, 2)
    for t, (x, y, z) in _CUBE.items():
        ops.sp(t, 2, gmax * x)                                  # u_y = gamma*X (engineering shear)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Penalty", 1.0e14, 1.0e14)
    ops.test("NormDispIncr", 1.0e-7, 100, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / 80)
    ops.analysis("Static")
    for _ in range(160):                                        # stage 1 (tension) + stage 2 (shear)
        assert ops.analyze(1) == 0, "non-proportional TS+interlock analyze failed"
    ops.eleResponse(1, "forces")
    sig = list(ops.eleResponse(1, "stresses"))[0:6]
    import math
    assert all(math.isfinite(x) for x in sig), sig            # robust: no NaN/inf at the boundary
    assert abs(sig[3]) > 1.0e-3, sig                          # shear actually developed (axes rotated)
    # stress stays bounded (no runaway) — combined non-proportional TS+interlock is well-behaved
    # even though the floor is not exactly preserved (the documented proportional-only boundary).
    assert max(abs(x) for x in sig) < 5.0 * FT, sig
