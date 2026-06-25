"""LadrunoRigidBody (ADR-58, ELE_TAG 33015) — Zone-A no-regression gate.

A 6-DOF rigid-body Element (LS-DYNA *CONSTRAINED_NODAL_RIGID_BODY analogue): a
zero-stiffness Element owning a private internal CoM Node, slaves tied by rigid-link
MP_Constraints, rotation integrated OFF the global solve in a side channel (spatial
angular momentum L + orientation quaternion, midpoint exp-map in commitState).

Ported from the standalone gates in Ladruno_scripts/rigidbody_tests/ (which keep the
tighter tolerances / longer runs); here trimmed for CI speed.

  P1   test_ballistic            free body falls as z=-1/2 g t^2; slaves follow exactly
  P2S1 test_dzhanibekov          torque-free spin: ||L|| conserved + intermediate-axis flip
  P2S2 test_slave_tracking       slaves track the FINITE rotation u_i=(R-I)d_i^0
  P2S2 test_gather_sign          a toe/contact reaction drives the spin (correct sign)
  P2S2 test_housner_rocking      free-rocking quarter-period vs Housner analytic
"""
import math

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

G = 9.81


def _quat_to_R(q):
    """rotation matrix from quaternion (qx,qy,qz,qs)."""
    x, y, z, w = q
    return [
        [1 - 2*(y*y + z*z), 2*(x*y - z*w),     2*(x*z + y*w)],
        [2*(x*y + z*w),     1 - 2*(x*x + z*z), 2*(y*z - x*w)],
        [2*(x*z - y*w),     2*(y*z + x*w),     1 - 2*(x*x + y*y)],
    ]


# ----------------------------------------------------------------------------- P1
def test_ballistic():
    """A free rigid body under gravity falls ballistically (a=-g, mass-independent)
    and every slave tracks the CoM exactly (R=I)."""
    m, nbody = 2.0, 3
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    ops.node(1, 0.0, 0.0, 10.0); ops.mass(1, m, m, m)
    ops.node(2, 1.0, 0.0, 10.0); ops.mass(2, m, m, m)
    ops.node(3, 0.0, 1.0, 10.0); ops.mass(3, m, m, m)
    ops.element('LadrunoRigidBody', 100, 3, 1, 2, 3)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for t in (1, 2, 3):
        ops.load(t, 0.0, 0.0, -m * G)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    dt, nsteps = 1.0e-4, 2000
    ops.analyze(nsteps, dt)
    t = nsteps * dt
    z_exact = -0.5 * G * t * t
    z_com = ops.eleResponse(100, 'comDisp')[2]
    assert abs(z_com - z_exact) / abs(z_exact) < 1e-3
    assert abs(ops.nodeDisp(1, 3) - z_com) < 1e-9          # exact follower (R=I)


# --------------------------------------------------------------------------- P2 S1
def test_dzhanibekov():
    """Torque-free asymmetric top: spatial angular momentum ||L|| is conserved to
    machine precision (momentum-form integrator) and the intermediate-axis (tennis-
    racket) instability fires — the body-frame omega_y reverses sign."""
    a, b, c = 2.0, 1.5, 1.0          # a>b>c -> intermediate principal axis = y
    pts = [(a, 0, 0), (-a, 0, 0), (0, b, 0), (0, -b, 0), (0, 0, c), (0, 0, -c)]
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for i, (x, y, z) in enumerate(pts, start=1):
        ops.node(i, float(x), float(y), float(z))
        ops.mass(i, 1.0, 1.0, 1.0)
    ops.element('LadrunoRigidBody', 100, len(pts), *range(1, len(pts) + 1),
                '-omega', 0.12, 6.0, 0.12)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    L0 = ops.eleResponse(100, 'angularMom')
    L0n = math.sqrt(sum(v * v for v in L0))
    dt, nsteps = 1.0e-3, 15000
    maxdev, signs = 0.0, []
    for s in range(nsteps):
        ops.analyze(1, dt)
        if s % 50 == 0:
            Ln = math.sqrt(sum(v * v for v in ops.eleResponse(100, 'angularMom')))
            maxdev = max(maxdev, abs(Ln - L0n) / L0n)
            signs.append(1 if ops.eleResponse(100, 'omegaBody')[1] >= 0 else -1)
    flips = sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])
    assert maxdev < 1e-10, f"||L|| not conserved (rel-dev {maxdev:.2e})"
    assert flips >= 1, "no intermediate-axis flip observed"


# --------------------------------------------------------------------------- P2 S2
def test_slave_tracking():
    """Finite-rotation slave following: each slave's displacement equals (R-I)d_i^0
    for the orientation it was imposed at (same-step), to machine precision."""
    a, b, c = 2.0, 1.5, 1.0
    pts = [(a, 0, 0), (-a, 0, 0), (0, b, 0), (0, -b, 0), (0, 0, c), (0, 0, -c)]
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for i, (x, y, z) in enumerate(pts, start=1):
        ops.node(i, float(x), float(y), float(z))
        ops.mass(i, 1.0, 1.0, 1.0)
    tags = list(range(1, len(pts) + 1))
    ops.element('LadrunoRigidBody', 100, len(pts), *tags, '-omega', 0.12, 6.0, 0.12)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    d0 = [list(p) for p in pts]                 # CoM = origin -> d_i^0 = x_i^0
    maxerr = 0.0
    for s in range(3000):
        q_before = ops.eleResponse(100, 'orientation')   # q used to impose this step
        ops.analyze(1, 1.0e-3)
        if s % 20 == 0:
            R = _quat_to_R(q_before)
            for i, tag in enumerate(tags):
                u = [ops.nodeDisp(tag, k + 1) for k in range(3)]
                exp = [sum(R[r][cc] * d0[i][cc] for cc in range(3)) - d0[i][r]
                       for r in range(3)]
                maxerr = max(maxerr, math.sqrt(sum((u[k] - exp[k])**2 for k in range(3))))
    assert maxerr < 1e-9, f"slaves do not track finite rotation (err {maxerr:.2e})"


def test_gather_sign():
    """A toe/contact reaction drives the spin: a torsional oscillator (two y-springs at
    opposite slaves) DECELERATES a +z spin (restoring couple), with no off-axis leak."""
    body = {1: (1, 0, 0), 2: (-1, 0, 0), 3: (0, 1, 0), 4: (0, -1, 0)}
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for t, (x, y, z) in body.items():
        ops.node(t, float(x), float(y), float(z))
        ops.mass(t, 1.0, 1.0, 1.0)
    ops.node(99, 1.0, 0.0, 0.0);  ops.fix(99, 1, 1, 1)
    ops.node(98, -1.0, 0.0, 0.0); ops.fix(98, 1, 1, 1)
    ops.uniaxialMaterial('Elastic', 1, 20.0)
    ops.element('zeroLength', 200, 99, 1, '-mat', 1, '-dir', 2)
    ops.element('zeroLength', 201, 98, 2, '-mat', 1, '-dir', 2)
    Omega = 1.0
    ops.element('LadrunoRigidBody', 100, 4, 1, 2, 3, 4, '-omega', 0.0, 0.0, Omega)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    Lz0 = 4.0 * Omega                            # Iz = 4
    Lz_hist = []
    for s in range(5000):                        # ~ a quarter+ of the torsional period
        ops.analyze(1, 2.0e-4)
        Lz_hist.append(ops.eleResponse(100, 'angularMom')[2])
    L = ops.eleResponse(100, 'angularMom')
    assert min(Lz_hist) < 0.5 * Lz0, "spin never decelerated — gather not driving rotation"
    assert max(Lz_hist) < 1.05 * Lz0, "spin ACCELERATED — gather sign flipped"
    assert abs(L[0]) < 0.1 * Lz0 and abs(L[1]) < 0.1 * Lz0, "off-axis momentum leak"


def test_housner_rocking():
    """Housner free-rocking quarter-period: a block released from theta0<alpha rocking
    about a stiff-pinned toe matches I0*thetaddot = -W*Rdist*sin(alpha-theta), i.e.
    T_quarter = (1/p) acosh(1/(1-theta0/alpha)), p=sqrt(W*Rdist/I0). Exercises slaving +
    moment gather + side-channel SO(3) and their consistency end-to-end."""
    b, h, M = 0.2, 1.0, 1.0
    Rdist = math.hypot(b, h)
    alpha = math.atan2(b, h)
    I0 = M * (b*b + h*h) + M * Rdist * Rdist      # I_cm + M*Rdist^2 for 4 corner masses
    p = math.sqrt(M * G * Rdist / I0)
    theta0 = 0.5 * alpha
    T_quarter = (1.0 / p) * math.acosh(1.0 / (1.0 - theta0 / alpha))

    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    upright = {1: (0.0, 0.0, 0.0), 2: (-2*b, 0.0, 0.0),
               3: (0.0, 0.0, 2*h), 4: (-2*b, 0.0, 2*h)}
    c, s = math.cos(theta0), math.sin(theta0)     # tilt by theta0 about the toe (+y)
    for t, (x, y, z) in upright.items():
        ops.node(t, x * c + z * s, y, -x * s + z * c)
        ops.mass(t, M/4, M/4, M/4)
    ops.node(99, 0.0, 0.0, 0.0); ops.fix(99, 1, 1, 1)
    ops.uniaxialMaterial('Elastic', 1, 3.0e4)     # stiff elastic pin (toe penetration ~3e-4)
    ops.element('zeroLength', 200, 99, 1, '-mat', 1, 1, '-dir', 1, 3)
    ops.element('LadrunoRigidBody', 100, 4, 1, 2, 3, 4)
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    for t in (1, 2, 3, 4):
        ops.load(t, 0.0, 0.0, -(M/4) * G)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    # dt=2e-5 keeps the CI run ~10s; the penalty/lag bias ~k*dt gives ~-2.5% here
    # (the standalone gate runs dt=1e-5 -> -1.2% at the tighter 3% tol).
    dt = 2.0e-5
    psi_prev, t_cross = 0.0, None
    for step in range(int(1.5 * T_quarter / dt)):
        ops.analyze(1, dt)
        q = ops.eleResponse(100, 'orientation')
        psi = 2.0 * math.atan2(q[1], q[3])        # body rotation about +y from start
        if psi <= -theta0:                        # reached upright
            frac = (psi_prev + theta0) / (psi_prev - psi) if psi != psi_prev else 1.0
            t_cross = (step + frac) * dt
            break
        psi_prev = psi
    assert t_cross is not None, "block never reached upright — rocking did not occur"
    rel = abs(t_cross - T_quarter) / T_quarter
    assert rel < 0.04, f"Housner quarter-period off by {rel*100:.1f}% (> 4%)"
