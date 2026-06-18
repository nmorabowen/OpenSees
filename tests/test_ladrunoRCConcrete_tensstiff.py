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
