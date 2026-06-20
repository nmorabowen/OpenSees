"""LadrunoCohesiveHingeBiaxial (coupled Mz-My cohesive interaction hinge, ND
classTag 33004) -- material-point battery for the CONSISTENT off-radial tangent
and the Benzeggagh-Kenane mode-mix (ADR 34 items #1 + #2).

The order-2 nDMaterial carries the rotation-jump vector alpha = [alpha_z, alpha_y]
as "strain" and the cohesive moment vector M = [M_z, M_y] as "stress", with a 2x2
tangent dM/dalpha. We drive the bare material through a `zeroLengthND` (the free
node's two translational DOFs ARE alpha), prescribe alpha exactly, and read the raw
2x2 tangent (material response "tangent").

Two coupled facts being pinned here:

  * With the DEFAULT linear fracture-energy mix (eta=1), the exact-reduction
    calibration Kpen_i = ratio*Mc_i^2/(2 Gf_i) forces S = (2/ratio) Gf_mix, so the
    softening rate A = S/Esoft = 2/(ratio-1) is mix-INDEPENDENT: the damage D(r)
    does not depend on the mode mix and the v1 "frozen-mix" tangent is ALREADY the
    exact off-radial Jacobian.
  * With a Benzeggagh-Kenane mix (eta != 1), Gf_mix = Gf_z + (Gf_y-Gf_z) w_y^eta is
    no longer proportional to S, so A and hence D become mix-DEPENDENT. Now the
    consistent tangent's mode-mix term dD/dw . dw/dalpha is LIVE -- the frozen-mix
    tangent is wrong off-radial, and the hardened tangent matches the true Jacobian.

Gates:
  * stress matches the closed form at mixed states (eta = 0.5, 1, 2);
  * CONSISTENT TANGENT == central FD of M(alpha) off-radial, for eta in {0.5,1,2};
  * the term BITES: for eta != 1 the frozen-mix tangent deviates from FD while the
    full (consistent) tangent matches -- so the mode-mix term is necessary;
  * mix-dependence of D flips on with eta (invariant at eta=1, varies at eta!=1);
  * pure-axis reduction stays exact for any eta (Gf_z / Gf_y, diagonal tangent);
  * B-K dissipation: a 45deg radial break dissipates Gf_z + (Gf_y-Gf_z)*(1/2)^eta;
  * pre-peak closed tangent is the diagonal penalty.

Plan: Ladruno_implementation/34_ladruno_cohesive_hinge_biaxial_adr.md.
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# law parameters (Gfz != Gfy so the mode-mix sensitivity is non-trivial)
_McZ, _GfZ = 50.0, 20.0
_McY, _GfY = 35.0, 12.0
_RATIO = 1000.0


# --------------------------------------------------------------------------- #
#  closed-form reference (mirrors LadrunoCohesiveHingeBiaxial on a loading path:
#  r_max = r(alpha) for a monotone radial drive, so M is a pure function of alpha)
# --------------------------------------------------------------------------- #
def _derived(ratio=_RATIO):
    Kpz = ratio * _McZ * _McZ / (2.0 * _GfZ)
    Kpy = ratio * _McY * _McY / (2.0 * _GfY)
    return Kpz, Kpy, _McZ / Kpz, _McY / Kpy            # Kpz, Kpy, a0z, a0y


def _law(az, ay, shape, eta, ratio=_RATIO):
    """returns ez, ey, r, wz, wy, S, Esoft, Teff, dTeff, D on the loading surface."""
    Kpz, Kpy, a0z, a0y = _derived(ratio)
    ez, ey = az / a0z, ay / a0y
    r = math.hypot(ez, ey)
    if r > 1.0e-300:
        wz, wy = ez * ez / (r * r), ey * ey / (r * r)
    else:
        wz, wy = 1.0, 0.0
    cz, cy = _McZ * _McZ / Kpz, _McY * _McY / Kpy
    S = wz * cz + wy * cy
    Gfmix = _GfZ + (_GfY - _GfZ) * (wy ** eta)
    Esoft = Gfmix - 0.5 * S
    if r <= 1.0:
        Teff, dTeff = S * r, S
    else:
        s = r - 1.0
        if shape == "linear":
            rf = 2.0 * Esoft / S
            if s >= rf:
                Teff, dTeff = 0.0, 0.0
            else:
                Teff, dTeff = S * (1.0 - s / rf), -S / rf
        else:
            A = S / Esoft
            e = math.exp(-A * s)
            Teff, dTeff = S * e, -A * S * e
    D = 0.0
    if r > 1.0e-300 and S > 0.0:
        D = min(max(1.0 - Teff / (S * r), 0.0), 1.0)
    return ez, ey, r, wz, wy, S, Esoft, Teff, dTeff, D


def _M_analytic(az, ay, shape="exp", eta=1.0, ratio=_RATIO):
    Kpz, Kpy, _, _ = _derived(ratio)
    _, _, _, _, _, _, _, _, _, D = _law(az, ay, shape, eta, ratio)
    return (1.0 - D) * Kpz * az, (1.0 - D) * Kpy * ay


def _tangent_analytic(az, ay, shape="exp", eta=1.0, include_mix=True, ratio=_RATIO):
    """raw 2x2 [K00,K01,K10,K11]; include_mix=False reproduces the v1 frozen-mix tangent."""
    Kpz, Kpy, a0z, a0y = _derived(ratio)
    ez, ey, r, wz, wy, S, Esoft, Teff, dTeff, D = _law(az, ay, shape, eta, ratio)
    K = [(1.0 - D) * Kpz, 0.0, 0.0, (1.0 - D) * Kpy]
    if r > 1.0 and S > 0.0:
        s = r - 1.0
        dDdr = (Teff - dTeff * r) / (S * r * r)
        dDaz = dDdr * az / (a0z * a0z * r)
        dDay = dDdr * ay / (a0y * a0y * r)
        if include_mix and eta != 1.0 and wz > 0.0 and wy > 0.0:
            cz, cy = _McZ * _McZ / Kpz, _McY * _McY / Kpy
            Sp = cz - cy
            Gp = eta * (_GfZ - _GfY) * (wy ** (eta - 1.0))
            Ep = Gp - 0.5 * Sp
            Ap = (Sp * Esoft - S * Ep) / (Esoft * Esoft)
            dDdw = (s / (2.0 * r)) * Ap if shape == "linear" else (s / r) * (Teff / S) * Ap
            r4 = r ** 4
            dDaz += dDdw * (2.0 * ez * ey * ey / (r4 * a0z))
            dDay += dDdw * (-2.0 * ez * ez * ey / (r4 * a0y))
        K[0] -= Kpz * az * dDaz
        K[1] -= Kpz * az * dDay
        K[2] -= Kpy * ay * dDaz
        K[3] -= Kpy * ay * dDay
    return K


# --------------------------------------------------------------------------- #
#  driver: a zeroLengthND material point (ndm 2 / ndf 3 -> numDOF 6; the order-2
#  ND material maps the two translational DOFs to alpha = [alpha_z, alpha_y]).  A
#  parallel zeroLength spring keeps the unused rotational DOF non-singular (it stays
#  0); alpha is prescribed exactly via sp + LoadControl (one monotone step ->
#  loading, r_max = r(alpha)).
# --------------------------------------------------------------------------- #
def _drive(az, ay, shape="-exp", eta=1.0, ratio=_RATIO):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)
    ops.fix(1, 1, 1, 1)
    args = [79, _McZ, _GfZ, _McY, _GfY, shape, "-penaltyRatio", ratio]
    if eta != 1.0:
        args += ["-bk", eta]
    ops.nDMaterial("LadrunoCohesiveHingeBiaxial", *args)
    ops.element("zeroLengthND", 1, 1, 2, 79)
    ops.uniaxialMaterial("Elastic", 91, 1.0e6)
    ops.element("zeroLength", 2, 1, 2, "-mat", 91, "-dir", 3)   # rotational DOF stiffness (stays 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.sp(2, 1, az)
    ops.sp(2, 2, ay)
    ops.system("FullGeneral")
    ops.numberer("Plain")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1.0e-13, 20)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, "material-point drive failed"
    Mz = ops.eleResponse(1, "material", "momentZ")[0]
    My = ops.eleResponse(1, "material", "momentY")[0]
    K = ops.eleResponse(1, "material", "tangent")        # [K00, K01, K10, K11]
    return Mz, My, K


# --------------------------------------------------------------------------- #
#  stress matches the closed form at mixed states (validates the FD reference)
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("eta", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("shape", ["-exp", "-linear"])
@pytest.mark.parametrize("ez,ey", [(1.2, 0.8), (0.7, 1.4), (1.0, 1.0)])
def test_stress_matches_closed_form(eta, shape, ez, ey):
    _, _, a0z, a0y = _derived()
    az, ay = ez * a0z, ey * a0y
    Mz, My, _ = _drive(az, ay, shape, eta)
    mz, my = _M_analytic(az, ay, shape.lstrip("-"), eta)
    assert math.isclose(Mz, mz, rel_tol=1e-9, abs_tol=1e-12), (Mz, mz)
    assert math.isclose(My, my, rel_tol=1e-9, abs_tol=1e-12), (My, my)


# --------------------------------------------------------------------------- #
#  THE consistent-tangent gate: raw 2x2 == central FD of M(alpha) off-radial, AND
#  == the analytic Jacobian.  For eta != 1 the mode-mix term is live (verified by
#  the discrimination assert: the frozen-mix tangent deviates from FD).
# --------------------------------------------------------------------------- #
# Run at a SOFTER penalty (ratio=10): at the near-rigid default (1000) the damage is
# dominated by the mix-independent secant geometry (1 - 1/r), so the mode-mix term is
# real but numerically tiny (~1e-4); ratio=10 makes it a clear ~1e-2 so the gate bites.
@pytest.mark.parametrize("eta", [0.5, 1.0, 2.0])
@pytest.mark.parametrize("shape", ["-exp", "-linear"])
@pytest.mark.parametrize("ez,ey", [(1.2, 0.9), (0.9, 1.4), (1.6, 0.7)])
def test_consistent_tangent_matches_fd_offradial(eta, shape, ez, ey):
    ratio = 10.0
    _, _, a0z, a0y = _derived(ratio)
    az, ay = ez * a0z, ey * a0y
    sh = shape.lstrip("-")
    _, _, K = _drive(az, ay, shape, eta, ratio)

    # central FD of M(alpha): each endpoint is the loading-surface value M(alpha)
    hz, hy = 1e-6 * a0z, 1e-6 * a0y
    J = [
        (_M_analytic(az + hz, ay, sh, eta, ratio)[0] - _M_analytic(az - hz, ay, sh, eta, ratio)[0]) / (2 * hz),
        (_M_analytic(az, ay + hy, sh, eta, ratio)[0] - _M_analytic(az, ay - hy, sh, eta, ratio)[0]) / (2 * hy),
        (_M_analytic(az + hz, ay, sh, eta, ratio)[1] - _M_analytic(az - hz, ay, sh, eta, ratio)[1]) / (2 * hz),
        (_M_analytic(az, ay + hy, sh, eta, ratio)[1] - _M_analytic(az, ay - hy, sh, eta, ratio)[1]) / (2 * hy),
    ]
    scale = max(abs(x) for x in J) + 1.0

    # C++ analytic tangent == FD (exact to FD precision), and == my independent analytic tangent
    Ka = _tangent_analytic(az, ay, sh, eta, include_mix=True, ratio=ratio)
    for k in range(4):
        assert abs(K[k] - J[k]) <= 1e-6 * scale, ("cpp-vs-fd", k, K[k], J[k])
        assert abs(K[k] - Ka[k]) <= 1e-7 * scale, ("cpp-vs-analytic", k, K[k], Ka[k])

    # the off-diagonals are genuinely live here
    assert abs(J[1]) > 1e-3 * scale and abs(J[2]) > 1e-3 * scale, J

    # discrimination: the frozen-mix (v1) tangent equals FD only when eta == 1; for a
    # nonlinear mix it is measurably WRONG off-radial -> the mode-mix term is necessary.
    Kf = _tangent_analytic(az, ay, sh, eta, include_mix=False, ratio=ratio)
    gap = max(abs(Kf[k] - J[k]) for k in range(4))
    if eta == 1.0:
        assert gap <= 1e-6 * scale, ("eta=1 frozen should match", gap)
    else:
        assert gap > 1e-3 * scale, ("mode-mix term must bite for eta!=1", gap, Kf, J)


# --------------------------------------------------------------------------- #
#  mix-dependence of D flips on with eta: invariant at eta=1, varies at eta!=1
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("shape", ["-exp", "-linear"])
def test_damage_mix_invariance_only_at_eta1(shape):
    ratio = 10.0
    _, _, a0z, a0y = _derived(ratio)
    rr = 2.2

    def dmg(ez, ey, eta):
        _drive(ez * a0z, ey * a0y, shape, eta, ratio)
        return ops.eleResponse(1, "material", "damage")[0]

    # eta = 1: damage depends only on r (mix-invariant -- the v1 mix-cancellation holds
    # for ANY penalty ratio because S stays proportional to Gf_mix)
    d_z = dmg(rr, 0.0, 1.0)
    d_45 = dmg(rr / math.sqrt(2), rr / math.sqrt(2), 1.0)
    assert abs(d_z - d_45) <= 1e-9, (d_z, d_45)

    # eta = 2: the nonlinear B-K mix makes damage vary with the mix at fixed r
    e_z = dmg(rr, 0.0, 2.0)
    e_45 = dmg(rr / math.sqrt(2), rr / math.sqrt(2), 2.0)
    assert abs(e_z - e_45) > 1e-3, (e_z, e_45)


# --------------------------------------------------------------------------- #
#  pure-axis reduction stays exact for any eta; tangent is diagonal there
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("eta", [0.5, 2.0])
@pytest.mark.parametrize("shape", ["-exp", "-linear"])
def test_pure_axis_unaffected_by_eta(eta, shape):
    _, _, a0z, a0y = _derived()
    sh = shape.lstrip("-")
    az = 1.4 * a0z                                        # pure strong-axis, softening
    Mz, My, K = _drive(az, 0.0, shape, eta)
    mz, _ = _M_analytic(az, 0.0, sh, eta)                 # Gf_mix = Gf_z regardless of eta
    assert math.isclose(Mz, mz, rel_tol=1e-9), (Mz, mz)
    assert abs(My) <= 1e-12, My
    diag = abs(K[0]) + abs(K[3]) + 1.0
    assert abs(K[1]) <= 1e-6 * diag and abs(K[2]) <= 1e-6 * diag, K


# --------------------------------------------------------------------------- #
#  B-K dissipation on a 45deg radial break: Gf_mix = Gf_z + (Gf_y-Gf_z)*(1/2)^eta
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("eta", [0.5, 1.0, 2.0])
def test_bk_dissipation_on_45deg(eta):
    # drive a 45deg-in-normalized-space radial path well past the break, sampling the
    # path work; the dissipated energy is the B-K mixed-mode fracture energy.
    _, _, a0z, a0y = _derived()
    Kpz, Kpy = _derived()[0], _derived()[1]
    # break radius: LINEAR reaches 0 at r = 1 + rf; use a generous end well beyond it.
    wz = wy = 0.5
    S = wz * _McZ * _McZ / Kpz + wy * _McY * _McY / Kpy
    Gfmix = _GfZ + (_GfY - _GfZ) * (wy ** eta)
    rf = 2.0 * (Gfmix - 0.5 * S) / S
    r_end = 1.0 + 1.05 * rf
    n = 4000
    # integrate int M . dalpha along the radial path az=ez*a0z, ay=ey*a0y, ez=ey=r/sqrt2
    e_path = 0.0
    prev = None
    for i in range(n + 1):
        r = r_end * i / n
        ez = ey = r / math.sqrt(2.0)
        az, ay = ez * a0z, ey * a0y
        mz, my = _M_analytic(az, ay, "linear", eta)
        if prev is not None:
            paz, pay, pmz, pmy = prev
            e_path += 0.5 * ((pmz + mz) * (az - paz) + (pmy + my) * (ay - pay))
        prev = (az, ay, mz, my)
    assert math.isclose(e_path, Gfmix, rel_tol=1e-3), (e_path, Gfmix, eta)


# --------------------------------------------------------------------------- #
#  pre-peak (closed hinge): tangent is the diagonal penalty, no coupling
# --------------------------------------------------------------------------- #
def test_prepeak_tangent_is_diagonal_penalty():
    Kpz, Kpy, a0z, a0y = _derived()
    az, ay = 0.4 * a0z, 0.3 * a0y                         # r < 1 -> closed
    _, _, K = _drive(az, ay, "-exp", 2.0)
    assert math.isclose(K[0], Kpz, rel_tol=1e-9), (K[0], Kpz)
    assert math.isclose(K[3], Kpy, rel_tol=1e-9), (K[3], Kpy)
    assert abs(K[1]) <= 1e-9 * Kpz and abs(K[2]) <= 1e-9 * Kpy, K
