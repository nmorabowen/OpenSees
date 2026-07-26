"""LadrunoCohesiveHinge (Ladruno discrete cohesive moment-rotation hinge,
uniaxial classTag 33003) -- Zone-A material battery.

The Tier-2 building block for the LadrunoDispBeamColumn embedded strong-
discontinuity hinge (ADR 32 sec.Tier-2). A UniaxialMaterial whose "strain" is the
rotation jump [[theta]] and whose "stress" is the cohesive moment M([[theta]]).
The hinge carries a *fracture energy* Gf directly -- there is no characteristic
length to calibrate (the whole point vs. Tier-1 lch smearing).

Law: near-rigid penalty pre-peak (hinge closed) up to M = Mc, then softening
(exponential default / linear option) calibrated so the monotonic envelope
integrates to EXACTLY Gf. Irreversible secant-to-origin damage.

Correctness anchors (the cheap solver-independent gates ADR 32 says to build
FIRST, before any element wiring or path-follower):

  * **energy gate** int M d[[theta]] == Gf.  LINEAR softening reaches M=0 at a
    finite jump so the area is exact to machine precision; EXP converges to Gf as
    the jump grows (tail truncation only).
  * **peak == Mc**: the cohesive law activates at M = Mc and never exceeds it.
  * **irreversibility**: unload returns to ~0 moment at zero jump (secant to
    origin); reload shows no strength recovery.
  * **near-rigid pre-peak**: initial tangent == the penalty Kpen; the hinge stays
    essentially closed (jump ~ 0) until M reaches Mc.
  * **guarded penalty floor**: a -penalty below Mc^2/(2Gf) is clamped (warned) and
    the material still dissipates exactly Gf.
  * **symmetry**: negative jump mirrors positive (single isotropic damage var).
  * **state cycle**: commit / revert / sendSelf-recvSelf preserve the law + history.

Plan: Ladruno_implementation/32_ladruno_dispbeamcolumn_regularization_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]


# --------------------------------------------------------------------------- #
#  driver: a 1-DOF zeroLength with the hinge material; the DOF IS the jump,
#  the element force IS the cohesive moment M.  Displacement-controlled so the
#  softening branch (negative tangent) is traced robustly.
# --------------------------------------------------------------------------- #
def _build(mat_args):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0)
    ops.node(2, 0.0)
    ops.fix(1, 1)
    ops.uniaxialMaterial("LadrunoCohesiveHinge", 1, *mat_args)
    ops.element("zeroLength", 1, 1, 2, "-mat", 1, "-dir", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(2, 1.0)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Plain")
    ops.test("NormDispIncr", 1.0e-12, 50)
    ops.algorithm("Newton")


def _moment():
    return ops.eleResponse(1, "material", 1, "stress")[0]


def _jump():
    return ops.nodeDisp(2, 1)


def _mat_energy():
    return ops.eleResponse(1, "material", 1, "energy")[0]


def _push(targets):
    """Drive node-2 displacement through the ordered list of jump `targets`
    (each linear segment substepped) and return the sampled (jump, moment) path
    plus the final material-reported energy.  Returning the path lets the test
    integrate int M d[[theta]] in python independently of the material's own
    getEnergy."""
    path = [(_jump(), _moment())]
    cur = _jump()
    for tgt, nsub in targets:
        incr = (tgt - cur) / nsub
        ops.integrator("DisplacementControl", 2, 1, incr)
        ops.analysis("Static")
        for _ in range(nsub):
            assert ops.analyze(1) == 0, "displacement-controlled step failed"
            path.append((_jump(), _moment()))
        cur = tgt
    return path, _mat_energy()


def _trap(path):
    """trapezoidal int M d[[theta]] over a (jump, moment) path."""
    e = 0.0
    for (k0, m0), (k1, m1) in zip(path, path[1:]):
        e += 0.5 * (m0 + m1) * (k1 - k0)
    return e


# analytic calibration mirrors LadrunoCohesiveHinge::computeDerived
def _kpen_default(Mc, Gf, ratio=1000.0):
    return ratio * Mc * Mc / (2.0 * Gf)


def _kappa0(Mc, Gf, Kpen):
    return Mc / Kpen


def _kappaf(Mc, Gf, Kpen):
    Esoft = Gf - Mc * Mc / (2.0 * Kpen)
    return 2.0 * Esoft / Mc


_Mc = 10.0
_Gf = 4.0


# --------------------------------------------------------------------------- #
#  parsing / construction
# --------------------------------------------------------------------------- #
def test_builds_exp_and_linear_and_options():
    for args in (
        (_Mc, _Gf),                                  # default: exp, default penalty
        (_Mc, _Gf, "-linear"),
        (_Mc, _Gf, "-exp"),
        (_Mc, _Gf, "-penaltyRatio", 500.0),
        (_Mc, _Gf, "-penalty", 1.0e4, "-linear"),
    ):
        _build(args)  # builds + a trivial solve setup; no exception == pass
    assert True


def test_rejects_nonpositive_Mc_Gf():
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)
    for bad in ((-1.0, _Gf), (_Mc, 0.0), (0.0, 0.0)):
        with pytest.raises(Exception):
            ops.uniaxialMaterial("LadrunoCohesiveHinge", 99, *bad)


# --------------------------------------------------------------------------- #
#  THE energy gate -- the cheap solver-independent capstone (ADR 32)
# --------------------------------------------------------------------------- #
def test_energy_gate_linear_exact():
    """LINEAR softening: drive to the fully-broken jump exactly; the area under
    the piecewise-linear envelope is Gf to machine precision.  The push targets
    land ON the kink (kappa0) and the break (kappa0+kappaf) so trapezoidal
    integration of each linear segment is exact."""
    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    kf = _kappaf(_Mc, _Gf, Kpen)
    _build((_Mc, _Gf, "-linear"))
    path, e_mat = _push([(k0, 200), (k0 + kf, 800)])

    e_path = _trap(path)
    assert math.isclose(e_path, _Gf, rel_tol=1.0e-9), (e_path, _Gf)
    # the material's own getEnergy agrees
    assert math.isclose(e_mat, _Gf, rel_tol=1.0e-9), (e_mat, _Gf)


def test_energy_gate_exp_converges():
    """EXP softening: drive to many softening characteristic lengths; the area
    converges to Gf (the only error is the truncated exponential tail)."""
    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    Esoft = _Gf - _Mc * _Mc / (2.0 * Kpen)
    a = _Mc / Esoft
    k_end = k0 + 25.0 / a          # tail moment ~ Mc*e^-25 -> negligible
    _build((_Mc, _Gf, "-exp"))
    path, e_mat = _push([(k0, 200), (k_end, 2000)])

    e_path = _trap(path)
    assert math.isclose(e_path, _Gf, rel_tol=2.0e-3), (e_path, _Gf)
    assert math.isclose(e_mat, _Gf, rel_tol=2.0e-3), (e_mat, _Gf)


# --------------------------------------------------------------------------- #
#  peak == Mc, activation, near-rigid pre-peak
# --------------------------------------------------------------------------- #
def test_peak_equals_Mc_and_near_rigid_preeak():
    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    kf = _kappaf(_Mc, _Gf, Kpen)
    _build((_Mc, _Gf, "-linear"))
    # Stop just SHY of full exhaustion: at kappa >= kappaf the cohesive moment
    # and tangent are both exactly zero, so the 1-DOF static tangent is
    # legitimately singular. Pre-ADR-76 the LAPACK solver SWALLOWED that
    # singular solve (rc 0) and the push limped through it; post-fix it
    # correctly refuses. Every assertion below (peak == Mc at kappa0, near-
    # rigid pre-peak) is decided long before 0.999*kappaf.
    path, _ = _push([(k0 + 0.999 * kf, 1000)])

    moments = [m for _, m in path]
    peak = max(moments)
    # cohesive law activates at Mc and never exceeds it
    assert math.isclose(peak, _Mc, rel_tol=1.0e-6), peak
    # the jump at peak is the (tiny) activation jump kappa0 -> near-rigid pre-peak
    k_at_peak = next(k for k, m in path if math.isclose(m, peak, rel_tol=1e-9))
    assert math.isclose(k_at_peak, k0, rel_tol=1.0e-6), (k_at_peak, k0)
    # near-rigid: activation jump is 1/penaltyRatio of the post-peak span
    assert k0 < 0.01 * kf


def test_initial_tangent_is_penalty():
    Kpen = _kpen_default(_Mc, _Gf)
    _build((_Mc, _Gf, "-exp"))
    # ops.getEleTangentStiff not exposed for materials; probe via a tiny pre-peak
    # step: M/jump in the elastic penalty regime == Kpen
    k0 = _kappa0(_Mc, _Gf, Kpen)
    path, _ = _push([(0.1 * k0, 5)])
    k, m = path[-1]
    assert math.isclose(m / k, Kpen, rel_tol=1.0e-6), (m / k, Kpen)


# --------------------------------------------------------------------------- #
#  irreversibility: secant-to-origin unload, no strength recovery on reload
# --------------------------------------------------------------------------- #
def test_irreversible_secant_unload_and_no_recovery():
    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    kf = _kappaf(_Mc, _Gf, Kpen)
    k_damaged = k0 + 0.5 * kf       # half-way down the softening branch
    _build((_Mc, _Gf, "-linear"))

    # load into softening, capture the envelope moment there
    path1, _ = _push([(k_damaged, 500)])
    m_env = path1[-1][1]
    assert 0.0 < m_env < _Mc

    # unload to zero jump -> moment returns to ~0 (secant to origin, no residual)
    path2, _ = _push([(0.0, 500)])
    k_zero, m_zero = path2[-1]
    assert abs(k_zero) < 1.0e-12
    assert abs(m_zero) < 1.0e-9 * _Mc

    # reload -> follows the secant; peak on reload must NOT exceed the damaged
    # envelope moment (no strength recovery)
    path3, _ = _push([(k_damaged, 500)])
    assert max(m for _, m in path3) <= m_env * (1.0 + 1.0e-6)


# --------------------------------------------------------------------------- #
#  symmetry: negative jump mirrors positive
# --------------------------------------------------------------------------- #
def test_symmetric_in_jump_sign():
    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    kf = _kappaf(_Mc, _Gf, Kpen)
    k_probe = k0 + 0.3 * kf

    _build((_Mc, _Gf, "-linear"))
    pos, _ = _push([(k_probe, 400)])
    m_pos = pos[-1][1]

    _build((_Mc, _Gf, "-linear"))
    neg, _ = _push([(-k_probe, 400)])
    m_neg = neg[-1][1]

    assert math.isclose(m_pos, -m_neg, rel_tol=1.0e-9), (m_pos, m_neg)


# --------------------------------------------------------------------------- #
#  guarded penalty floor: below-floor penalty is clamped, Gf still exact
# --------------------------------------------------------------------------- #
def test_guarded_penalty_floor_still_dissipates_Gf():
    # floor = Mc^2/(2Gf); ask for half the floor -> must be clamped, not blow up
    floor = _Mc * _Mc / (2.0 * _Gf)
    _build((_Mc, _Gf, "-linear", "-penalty", 0.5 * floor))
    # with the clamp the penalty becomes the default 1000*floor; drive to break
    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    kf = _kappaf(_Mc, _Gf, Kpen)
    path, e_mat = _push([(k0, 200), (k0 + kf, 800)])
    assert math.isclose(_trap(path), _Gf, rel_tol=1.0e-9)
    assert math.isclose(e_mat, _Gf, rel_tol=1.0e-9)


# --------------------------------------------------------------------------- #
#  state cycle: commit advances history, sendSelf/recvSelf roundtrip preserves it
# --------------------------------------------------------------------------- #
def test_database_roundtrip_preserves_law():
    from _testbed.roundtrip import database_roundtrip

    Kpen = _kpen_default(_Mc, _Gf)
    k0 = _kappa0(_Mc, _Gf, Kpen)
    kf = _kappaf(_Mc, _Gf, Kpen)
    target = k0 + 0.4 * kf            # committed partway down the softening branch

    def build():
        _build((_Mc, _Gf, "-linear"))
        ops.integrator("DisplacementControl", 2, 1, target / 300.0)
        ops.analysis("Static")
        for _ in range(300):
            assert ops.analyze(1) == 0

    # probe_fn reads the committed cohesive HISTORY (moment, max-jump, energy) --
    # state that only survives if sendSelf/recvSelf round-trip CkappaMax/Cwork.
    database_roundtrip(
        build, [2], 1,
        probe_fn=lambda: [
            ops.eleResponse(1, "material", 1, "stress")[0],
            ops.eleResponse(1, "material", 1, "kappaMax")[0],
            ops.eleResponse(1, "material", 1, "energy")[0],
        ],
    )
