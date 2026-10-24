"""WP-205: ASDShellQ4/T3 geometric stiffness.

(a) -corotational: the stock tangent had NO geometric stiffness (EICR rotation
    gradient G zeroed under the polar-aligned frame, and the K_GP/K_GR sign was
    the Kratos one), so eigen after gravity + loadConst showed no softening and
    a P-Delta step needed ~34 Newton iterations.
(b) -pDelta (new): linear kinematics + the initial-stress stiffness of the
    current membrane forces, for linear eigen/RSA about gravity (ETABS-style).

Model: 30 m cantilever wall strip (B = 1, t = 0.3), tip mass, axial P = 0.5 Pcr.
Oracles: ShellNLDKGQ/ShellNLDKGT on the same mesh, the Euler period ratio
1/sqrt(1 - P/Pcr) = sqrt(2), the static amplification ~ 1/(1 - P/Pcr) = 2, and
P = 0 (nothing may change).
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

E, NU, T, B, H, M = 25.0e6, 0.2, 0.3, 1.0, 30.0, 100.0
PCR = math.pi ** 2 * E * (B * T ** 3 / 12.0) / (4.0 * H ** 2)


def _run(el, opt, pfac, tri=False, nx=2, nz=20, lateral_dir=2):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    tg = lambda i, k: 1 + i + k * (nx + 1)
    inplane = lateral_dir == 1   # sway in the wall's own plane: hold it planar
    for k in range(nz + 1):
        for i in range(nx + 1):
            ops.node(tg(i, k), i * B / nx, 0.0, k * H / nz)
            if k == 0:
                ops.fix(tg(i, k), 1, 1, 1, 1, 1, 1)
            elif inplane:
                ops.fix(tg(i, k), 0, 1, 0, 1, 0, 1)
    ops.section('ElasticMembranePlateSection', 1, E, NU, T, 0.0)
    e = 1
    for k in range(nz):
        for i in range(nx):
            a, b, c, d = tg(i, k), tg(i + 1, k), tg(i + 1, k + 1), tg(i, k + 1)
            if tri:
                ops.element(el, e, a, b, c, 1, *opt); e += 1
                ops.element(el, e, a, c, d, 1, *opt); e += 1
            else:
                ops.element(el, e, a, b, c, d, 1, *opt); e += 1
    top = [tg(i, nz) for i in range(nx + 1)]
    for n in top:
        ops.mass(n, M / len(top), 0.0 if inplane else M / len(top), 0, 0, 0, 0)
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for n in top:
        ops.load(n, 0, 0, -pfac * PCR / len(top), 0, 0, 0)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.test('NormDispIncr', 1.0e-10, 100)
    ops.algorithm('Newton'); ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.loadConst('-time', 0.0)
    lam = ops.eigen('-fullGenLapack', 1)[0]
    period = 2.0 * math.pi / math.sqrt(lam)
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    f = [0.0] * 6
    f[lateral_dir - 1] = 1.0 / len(top)
    for n in top:
        ops.load(n, *f)
    assert ops.analyze(1) == 0
    iters = ops.testIter()
    return period, iters, ops.nodeDisp(top[0], lateral_dir)


CASES = [('ASDShellQ4', ['-corotational'], False, 'ShellNLDKGQ'),
         ('ASDShellQ4', ['-pDelta'], False, 'ShellNLDKGQ'),
         ('ASDShellT3', ['-corotational'], True, 'ShellNLDKGT'),
         ('ASDShellT3', ['-pDelta'], True, 'ShellNLDKGT')]


@pytest.mark.parametrize("el,opt,tri,ref", CASES)
def test_softening_matches_reference_and_theory(el, opt, tri, ref):
    t0, _, u0 = _run(el, opt, 0.0, tri)
    tp, it, up = _run(el, opt, 0.5, tri)
    r0, _, _ = _run(ref, [], 0.0, tri)
    rp, _, _ = _run(ref, [], 0.5, tri)
    ratio, ref_ratio = tp / t0, rp / r0
    assert ratio == pytest.approx(ref_ratio, rel=2.0e-3), (ratio, ref_ratio)
    assert ratio == pytest.approx(math.sqrt(2.0), rel=1.0e-2), ratio
    # second-order static amplification of the lateral tip displacement
    assert up / u0 == pytest.approx(2.0, rel=2.0e-2), up / u0
    # a consistent tangent: Newton converges like the reference element
    assert it <= 8, it


@pytest.mark.parametrize("el,tri", [('ASDShellQ4', False), ('ASDShellT3', True)])
def test_zero_axial_load_unchanged(el, tri):
    lin = _run(el, [], 0.0, tri)
    for opt in (['-corotational'], ['-pDelta']):
        t, _, u = _run(el, opt, 0.0, tri)
        assert t == pytest.approx(lin[0], rel=1e-7)   # round-off: 1.5e-9 on Linux LAPACK
        assert u == pytest.approx(lin[2], rel=1e-4)


def test_pdelta_inplane_sway_matches_corotational_and_theory():
    """a wall swaying in its own plane (out-of-plane held): the u/v part of
    K_sigma. ShellNLDKGQ is NOT a reference here: it carries no in-plane
    geometric stiffness (its period ratio is 0.9999 on this case)."""
    pcr_in = math.pi ** 2 * E * (T * B ** 3 / 12.0) / (4.0 * H ** 2)
    pfac = 0.3 * pcr_in / PCR
    t0, _, _ = _run('ASDShellQ4', ['-pDelta'], 0.0, lateral_dir=1, nx=4)
    tp, it, _ = _run('ASDShellQ4', ['-pDelta'], pfac, lateral_dir=1, nx=4)
    c0, _, _ = _run('ASDShellQ4', ['-corotational'], 0.0, lateral_dir=1, nx=4)
    cp, _, _ = _run('ASDShellQ4', ['-corotational'], pfac, lateral_dir=1, nx=4)
    assert tp / t0 == pytest.approx(cp / c0, rel=1e-3)
    assert tp / t0 == pytest.approx(1.0 / math.sqrt(0.7), rel=5e-3)
    assert it <= 8


def test_linear_kinematics_has_no_softening():
    """without -pDelta / -corotational nothing changes with P (the baseline)."""
    t0, _, _ = _run('ASDShellQ4', [], 0.0)
    tp, _, _ = _run('ASDShellQ4', [], 0.5)
    assert tp == pytest.approx(t0, rel=1e-3)


def test_pdelta_and_corotational_are_exclusive():
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], start=1):
        ops.node(i, x, y, 0.0)
    ops.section('ElasticMembranePlateSection', 1, E, NU, T, 0.0)
    with pytest.raises(Exception):
        ops.element('ASDShellQ4', 1, 1, 2, 3, 4, 1, '-corotational', '-pDelta')


@pytest.mark.parametrize("el,tri", [('ASDShellQ4', False), ('ASDShellT3', True)])
def test_pdelta_wide_wall_frozen_preload_after_loadconst(el, tri):
    """WP-209: after loadConst, -pDelta holds K_sigma at the committed preload
    membrane forces N0 instead of re-evaluating N0 + dN_L on the lateral step.
    Without freezing N0:
      * Term A (K_sigma(dN_L) * u_0) injects a spurious P/(E*A) = 2.18% horizontal
        base-shear error on a wide wall under Linear (or shifts u_L under Newton);
      * Term B (K_sigma(dN_L) * du_L) injects a quadratic vertical reaction error
        sum(FZ) != P_grav and requires multiple Newton iterations."""
    bw, hw, tw = 4.0, 9.0, 0.25
    ew, nuw = 30.0e6, 0.2
    p_grav = 10000.0
    f_lat = 100.0
    nx, nz = 4, 9

    def _solve_wall(algo):
        ops.wipe()
        ops.model('basic', '-ndm', 3, '-ndf', 6)
        tg = lambda i, k: 1 + i + k * (nx + 1)
        for k in range(nz + 1):
            for i in range(nx + 1):
                ops.node(tg(i, k), i * bw / nx, 0.0, k * hw / nz)
                if k == 0:
                    ops.fix(tg(i, k), 1, 1, 1, 1, 1, 1)
        ops.section('ElasticMembranePlateSection', 1, ew, nuw, tw, 0.0)
        e = 1
        for k in range(nz):
            for i in range(nx):
                a, b, c, d = tg(i, k), tg(i + 1, k), tg(i + 1, k + 1), tg(i, k + 1)
                if tri:
                    ops.element(el, e, a, b, c, 1, '-pDelta'); e += 1
                    ops.element(el, e, a, c, d, 1, '-pDelta'); e += 1
                else:
                    ops.element(el, e, a, b, c, d, 1, '-pDelta'); e += 1
        top = [tg(i, nz) for i in range(nx + 1)]
        base = [tg(i, 0) for i in range(nx + 1)]
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)
        for i, n in enumerate(top):
            w = 0.5 if i in (0, nx) else 1.0
            ops.load(n, 0.0, 0.0, -p_grav * (w / nx), 0.0, 0.0, 0.0)
        ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
        ops.integrator('LoadControl', 1.0)
        ops.test('NormUnbalance', 1.0e-8, 50)
        ops.algorithm('Newton'); ops.analysis('Static')
        assert ops.analyze(1) == 0
        ops.loadConst('-time', 0.0)

        ops.timeSeries('Linear', 2)
        ops.pattern('Plain', 2, 2)
        for i, n in enumerate(top):
            w = 0.5 if i in (0, nx) else 1.0
            ops.load(n, f_lat * (w / nx), 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.algorithm(algo)
        assert ops.analyze(1) == 0
        iters = ops.testIter() if algo == 'Newton' else 1
        ops.reactions()
        rx = sum(ops.nodeReaction(n, 1) for n in base)
        rz = sum(ops.nodeReaction(n, 3) for n in base)
        ux = ops.nodeDisp(top[-1], 1)
        return rx, rz, ux, iters

    rx_lin, rz_lin, ux_lin, _ = _solve_wall('Linear')
    rx_newt, rz_newt, ux_newt, it_newt = _solve_wall('Newton')

    assert it_newt == 1
    assert rx_lin == pytest.approx(-f_lat, rel=1e-12)
    assert rx_newt == pytest.approx(-f_lat, rel=1e-12)
    assert rz_lin == pytest.approx(p_grav, abs=1e-10)
    assert rz_newt == pytest.approx(p_grav, abs=1e-10)
    assert ux_lin == pytest.approx(ux_newt, rel=1e-12)
@pytest.mark.parametrize('el,tri', [('ASDShellQ4', False), ('ASDShellT3', True)])
def test_asdshell_pdelta_loadconst_survives_database_save_restore(el, tri, tmp_path):
    """WP-209: Domain::sendSelf/recvSelf persists m_loadConstActive alongside
    ASDShellQ4/T3::sendSelf/recvSelf (m_pdelta_N, m_pdelta_frozen), so a lateral
    step after database restore(1) keeps the frozen preload geometric stiffness."""
    bw, hw, tw = 4.0, 9.0, 0.25
    ew, nuw = 30.0e6, 0.2
    p_grav = 10000.0
    f_lat = 100.0
    nx, nz = 2, 3

    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    tg = lambda i, k: 1 + i + k * (nx + 1)
    for k in range(nz + 1):
        for i in range(nx + 1):
            ops.node(tg(i, k), i * bw / nx, 0.0, k * hw / nz)
            if k == 0:
                ops.fix(tg(i, k), 1, 1, 1, 1, 1, 1)
    ops.section('ElasticMembranePlateSection', 1, ew, nuw, tw, 0.0)
    e = 1
    for k in range(nz):
        for i in range(nx):
            a, b, c, d = tg(i, k), tg(i + 1, k), tg(i + 1, k + 1), tg(i, k + 1)
            if tri:
                ops.element(el, e, a, b, c, 1, '-pDelta'); e += 1
                ops.element(el, e, a, c, d, 1, '-pDelta'); e += 1
            else:
                ops.element(el, e, a, b, c, d, 1, '-pDelta'); e += 1
    top = [tg(i, nz) for i in range(nx + 1)]
    base = [tg(i, 0) for i in range(nx + 1)]
    ops.timeSeries('Linear', 1)
    ops.pattern('Plain', 1, 1)
    for i, n in enumerate(top):
        w = 0.5 if i in (0, nx) else 1.0
        ops.load(n, 0.0, 0.0, -p_grav * (w / nx), 0.0, 0.0, 0.0)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.test('NormUnbalance', 1.0e-8, 50)
    ops.algorithm('Newton'); ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.loadConst('-time', 0.0)

    db_prefix = str(tmp_path / f"db_{el}")
    ops.database('File', db_prefix)
    ops.save(1)
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 6)
    ops.database('File', db_prefix)
    ops.restore(1)

    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    for i, n in enumerate(top):
        w = 0.5 if i in (0, nx) else 1.0
        ops.load(n, f_lat * (w / nx), 0.0, 0.0, 0.0, 0.0, 0.0)
    ops.system('UmfPack'); ops.numberer('RCM'); ops.constraints('Plain')
    ops.integrator('LoadControl', 1.0)
    ops.algorithm('Linear'); ops.analysis('Static')
    assert ops.analyze(1) == 0
    ops.reactions()
    rx = sum(ops.nodeReaction(n, 1) for n in base)
    rz = sum(ops.nodeReaction(n, 3) for n in base)
    assert rx == pytest.approx(-f_lat, rel=1e-6)
    assert rz == pytest.approx(p_grav, rel=1e-6)
