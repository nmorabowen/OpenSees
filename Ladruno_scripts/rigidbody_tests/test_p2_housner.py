"""
ADR-58 P2 Stage-2 GATE — Housner rocking block (the rocking-foundation MVP).

A rigid block rocking about a bottom toe is, for one swing, a physical pendulum:
    I0 * thetaddot = -W * Rdist * sin(alpha - theta)
with theta the tilt from upright, alpha = atan(b/h) the block slenderness, Rdist the
CoM-to-toe distance, and I0 the rotational inertia about the toe. Released FROM REST
at theta0 < alpha, the quarter period (theta0 -> upright) is Housner's
    T_quarter = (1/p) * acosh( 1 / (1 - theta0/alpha) ),   p = sqrt(W*Rdist/I0).

This end-to-end exercises the WHOLE Stage-2 chain at once:
  - finite-rotation slave following (the toe slave must track R*d^0 so the toe pin
    sees the rotated geometry),
  - the CoM moment gather (the pin reaction's moment about the CoM drives the spin),
  - the side-channel SO(3) integrator (rotation about the CoM with I_cm),
  - AND their consistency: the parallel-axis term M*Rdist^2 of I0 emerges from the
    CoM-TRANSLATION channel (the CoM swings on an arc about the pinned toe). If any
    piece or the translation/rotation coupling is wrong, the period misses.

HARNESS (ADR D5 fences the contact engine out): the toe is PINNED by a stiff elastic
zeroLength (dirs x,z) to a coincident fixed ground node — a pin joint, and an ELEMENT
so the incident-cache gather reads its reaction. The lifting corner is free. The pin
is the instantaneous rocking pivot; this reproduces Housner's EOM exactly for the
first quarter swing (no impact -> no restitution modelling needed; the penalty-impact
restitution is NOT reliably reproducible and is intentionally not gated here).

    python3.12 Ladruno_scripts/rigidbody_tests/test_p2_housner.py
"""
import os
import sys
import math

_dist = os.path.join(os.path.dirname(__file__), "..", "..", "dist", "bin")
if os.path.isfile(os.path.join(_dist, "opensees.pyd")):
    sys.path.insert(0, _dist)
    os.add_dll_directory(os.path.abspath(_dist))
    import opensees as ops
else:
    import openseespy.opensees as ops

os.environ["LADRUNO_OPENSEES_QUIET"] = "1"

# ---- block geometry / analytic Housner quarter period -----------------------
b, h = 0.2, 1.0                 # half-width, half-height (slender -> clean rocking)
M, g = 1.0, 9.81
Rdist = math.hypot(b, h)
alpha = math.atan2(b, h)
# 4 corner masses M/4 -> I_cm(y) = M*(b^2+h^2); I0 about toe = I_cm + M*Rdist^2
I_cm = M * (b*b + h*h)
I0 = I_cm + M * Rdist * Rdist
p = math.sqrt(M * g * Rdist / I0)
theta0 = 0.5 * alpha
T_quarter = (1.0 / p) * math.acosh(1.0 / (1.0 - theta0 / alpha))
print(f"alpha={alpha:.4f} rad  Rdist={Rdist:.4f}  I0={I0:.4f}  p={p:.4f}")
print(f"theta0={theta0:.4f} rad  analytic quarter-period T={T_quarter:.5f} s")

# ---- model: block tilted by theta0 about the toe O=(0,0,0) -------------------
ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# upright corners relative to the toe O (bottom-right): BR is the pivot at O
upright = {1: (0.0, 0.0, 0.0),       # BR = toe (pivot)
           2: (-2*b, 0.0, 0.0),      # BL
           3: (0.0, 0.0, 2*h),       # TR
           4: (-2*b, 0.0, 2*h)}      # TL
# tilt by +theta0 about +y (x' = x c + z s, z' = -x s + z c): CoM moves +x toward O
c, s = math.cos(theta0), math.sin(theta0)
for t, (x, y, z) in upright.items():
    xt = x * c + z * s
    zt = -x * s + z * c
    ops.node(t, xt, y, zt)
    ops.mass(t, M/4, M/4, M/4)

# coincident fixed ground node at the (tilted) toe and a stiff elastic pin (x,z).
# Penalty choice: the rotation (side channel) and CoM translation (global solve) are
# coupled through the pin with a one-step lag, so an over-stiff pin injects a parasitic
# restoring (period biased fast) that shrinks ~k*dt as dt->0 (verified convergence
# sweep). k=3e4 keeps the toe rigid for the physics (static penetration W/k ~ 3e-4,
# i.e. 3e-4 of the half-height) while dt=1e-5 lands the quarter period within ~1-2% of
# the Housner analytic.
ops.node(99, 0.0, 0.0, 0.0)
ops.fix(99, 1, 1, 1)
kpin = 3.0e4
ops.uniaxialMaterial('Elastic', 1, kpin)
ops.element('zeroLength', 200, 99, 1, '-mat', 1, 1, '-dir', 1, 3)   # pin toe in x and z

# the rigid body over the 4 corners (pin element already in the domain -> cached)
ops.element('LadrunoRigidBody', 100, 4, 1, 2, 3, 4)

# gravity at the corners (gathered to the CoM by the MPs; its moment about the toe is
# the rocking restoring, transmitted through the pin reaction)
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
for t in (1, 2, 3, 4):
    ops.load(t, 0.0, 0.0, -(M/4) * g)

ops.constraints('Transformation')
ops.numberer('Plain')
ops.system('Diagonal')
ops.test('NormDispIncr', 1e-12, 10, 0)
ops.algorithm('Linear')
ops.integrator('CentralDifference')
ops.analysis('Transient')

# ---- run: detect the first upright crossing (body rotated by theta0 from start) --
dt = 1.0e-5
nmax = int(1.5 * T_quarter / dt)
psi_prev = 0.0
t_cross = None
for step in range(nmax):
    ops.analyze(1, dt)
    q = ops.eleResponse(100, 'orientation')           # (qx,qy,qz,qs)
    psi = 2.0 * math.atan2(q[1], q[3])                # body rotation about +y from start
    # upright = body has un-tilted by theta0 (psi reaches -theta0)
    if t_cross is None and psi <= -theta0:
        # linear-interpolate the crossing time within this step
        frac = (psi_prev - (-theta0)) / (psi_prev - psi) if psi != psi_prev else 1.0
        t_cross = (step + frac) * dt
        break
    psi_prev = psi

assert t_cross is not None, "block never reached upright — rocking did not occur"
rel = abs(t_cross - T_quarter) / T_quarter
print(f"measured quarter-period = {t_cross:.5f} s   analytic = {T_quarter:.5f} s   "
      f"rel-err = {rel*100:.2f}%")
assert rel < 0.03, f"HOUSNER GATE FAILED: quarter-period rel-err {rel*100:.2f}% (> 3%)"
print("HOUSNER GATE PASS: free-rocking quarter-period within 3% of Housner analytic "
      "(slaving + moment gather + side-channel SO(3) consistent end-to-end)")
