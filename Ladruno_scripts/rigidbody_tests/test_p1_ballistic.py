"""
ADR-58 P1 smoke test — LadrunoRigidBody (zero-stiffness Element, D1(b)).

Gate (ADR 58 §7 P1):
  1. BALLISTIC: a free rigid body under gravity falls as z(t) = z0 - 1/2 g t^2,
     independent of the body mass (a = F/m = -g). This exercises the condensed
     CoM mass on the private internal node + the diagonal (translational) solve.
  2. EXACT FOLLOWER: with R = I in P1, every slaved node displacement equals the
     CoM displacement to machine precision (u_i = u_R).

Run with the fork's Python against the freshly-built module in dist/bin:
    python3.12 Ladruno_scripts/rigidbody_tests/test_p1_ballistic.py
"""
import os
import sys

# prefer the fork build in dist/bin (imports as `opensees`); fall back to a
# pip-installed openseespy if the local build is not present
_dist = os.path.join(os.path.dirname(__file__), "..", "..", "dist", "bin")
if os.path.isfile(os.path.join(_dist, "opensees.pyd")):
    sys.path.insert(0, _dist)
    os.add_dll_directory(_dist)
    import opensees as ops
else:
    import openseespy.opensees as ops

g = 9.81
m = 2.0            # mass per slave node
nbody = 3          # 3 non-collinear slaves -> non-singular inertia
W = nbody * m * g  # total weight

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# three non-collinear material points that make up the rigid body, each massed
ops.node(1, 0.0, 0.0, 10.0); ops.mass(1, m, m, m)
ops.node(2, 1.0, 0.0, 10.0); ops.mass(2, m, m, m)
ops.node(3, 0.0, 1.0, 10.0); ops.mass(3, m, m, m)

# the rigid body: condenses the slaves' mass to a private internal CoM node
ops.element('LadrunoRigidBody', 100, 3, 1, 2, 3)

# gravity as a constant load on the slaves (condensed to the CoM through the MPs)
ops.timeSeries('Constant', 1)
ops.pattern('Plain', 1, 1)
ops.load(1, 0.0, 0.0, -m * g)
ops.load(2, 0.0, 0.0, -m * g)
ops.load(3, 0.0, 0.0, -m * g)

ops.constraints('Transformation')          # eliminate the rigid-link MP DOFs
ops.numberer('Plain')
ops.system('Diagonal')                      # D3: diagonal translational mass
ops.test('NormDispIncr', 1e-12, 10, 0)
ops.algorithm('Linear')
ops.integrator('CentralDifference')
ops.analysis('Transient')

dt = 1.0e-4
nsteps = 2000
ops.analyze(nsteps, dt)
t = nsteps * dt

z_exact = -0.5 * g * t * t
z_com = ops.eleResponse(100, 'comDisp')[2]
z_s1 = ops.nodeDisp(1, 3)

rel = abs(z_com - z_exact) / abs(z_exact)
foll = abs(z_s1 - z_com)
print(f"t={t:.4f}  z_exact={z_exact:.6e}  z_com={z_com:.6e}  z_slave1={z_s1:.6e}")
print(f"  ballistic rel-err = {rel:.3e}   follower abs-err = {foll:.3e}")

assert rel < 1e-3,  f"BALLISTIC gate failed: rel-err {rel:.3e}"
assert foll < 1e-9, f"FOLLOWER gate failed: abs-err {foll:.3e}"
print("P1 GATE PASS: ballistic free-fall + exact node-follower")
