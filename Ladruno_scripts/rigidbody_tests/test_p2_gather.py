"""
ADR-58 P2 Stage-2 PROBE — CoM moment gather (gatherCoMTorque), sign + drive.

A torsional oscillator. The body spins about +z (omega0 = (0,0,Omega)). Two stiff
y-springs at opposite slaves form a RESTORING couple:
  - slave 1 at (+1,0,0): a +z spin moves it +y; a y-spring then pushes it -y
    -> moment r x f = (1,0,0) x (0,-f,0) = (0,0,-f)         (restoring, -z)
  - slave 2 at (-1,0,0): a +z spin moves it -y; the y-spring pushes it +y
    -> moment (-1,0,0) x (0,+f,0) = (0,0,-f)                (restoring, -z)
Net force = 0 (a pure couple -> CoM does not translate). So the +z spin must be
DECELERATED by the gathered moment: L_z falls below its initial value and the body
torsionally oscillates. With the gather broken (m=0) L_z would stay constant (free
spin). A flipped gather sign would ACCELERATE the spin (L_z grows) -> caught here.

Note: zeroLength force comes from relative DISPLACEMENT (the ground node is placed AT
the slave, not offset); the body's rotation supplies that displacement.

    python3.12 Ladruno_scripts/rigidbody_tests/test_p2_gather.py
"""
import os
import sys

_dist = os.path.join(os.path.dirname(__file__), "..", "..", "dist", "bin")
if os.path.isfile(os.path.join(_dist, "opensees.pyd")):
    sys.path.insert(0, _dist)
    os.add_dll_directory(os.path.abspath(_dist))
    import opensees as ops
else:
    import openseespy.opensees as ops

os.environ["LADRUNO_OPENSEES_QUIET"] = "1"

ops.wipe()
ops.model('basic', '-ndm', 3, '-ndf', 3)

# 4 unit masses in the xy-plane, CoM at origin; Iz = sum m (x^2+y^2) = 4
body = {1: (1, 0, 0), 2: (-1, 0, 0), 3: (0, 1, 0), 4: (0, -1, 0)}
for t, (x, y, z) in body.items():
    ops.node(t, float(x), float(y), float(z))
    ops.mass(t, 1.0, 1.0, 1.0)

# ground nodes co-located with slaves 1 and 2 (zeroLength needs them coincident),
# fixed; stiff y-springs resist the slaves' y-motion -> a restoring couple about z.
ops.node(99, 1.0, 0.0, 0.0);  ops.fix(99, 1, 1, 1)
ops.node(98, -1.0, 0.0, 0.0); ops.fix(98, 1, 1, 1)
k = 20.0
ops.uniaxialMaterial('Elastic', 1, k)
ops.element('zeroLength', 200, 99, 1, '-mat', 1, '-dir', 2)   # dir 2 = global y
ops.element('zeroLength', 201, 98, 2, '-mat', 1, '-dir', 2)

# the rigid body (springs 200/201 already in the domain -> incident cache sees them)
Omega = 1.0
ops.element('LadrunoRigidBody', 100, 4, 1, 2, 3, 4, '-omega', 0.0, 0.0, Omega)

ops.constraints('Transformation')
ops.numberer('Plain')
ops.system('Diagonal')
ops.test('NormDispIncr', 1e-12, 10, 0)
ops.algorithm('Linear')
ops.integrator('CentralDifference')
ops.analysis('Transient')

Iz = 4.0
Lz0 = Iz * Omega
dt = 2.0e-4   # torsional period ~ 2*pi/sqrt(2k/Iz) ~ 2 s -> run ~3 s to see recovery
Lz_hist = []
for s in range(15000):
    ops.analyze(1, dt)
    Lz_hist.append(ops.eleResponse(100, 'angularMom')[2])

L = ops.eleResponse(100, 'angularMom')
Lz_min = min(Lz_hist)
Lz_max = max(Lz_hist)
print(f"Lz0 = {Lz0:.4f}   final L = ({L[0]:.3e}, {L[1]:.3e}, {L[2]:.3e})")
print(f"Lz range over run: [{Lz_min:.4f}, {Lz_max:.4f}]")

# restoring couple => Lz must drop well below its initial value (gather drives rotation,
# correct sign). A broken gather keeps Lz==Lz0; a flipped sign pushes Lz above Lz0.
imin = Lz_hist.index(Lz_min)
recovery = max(Lz_hist[imin:])           # does Lz swing back up after the first dip?
assert Lz_min < 0.5 * Lz0,  f"GATHER FAILED: Lz never decelerated (min {Lz_min:.3f}, L0 {Lz0:.3f})"
assert Lz_max < 1.05 * Lz0, f"GATHER SIGN FLIPPED: Lz grew above L0 (max {Lz_max:.3f})"
assert abs(L[0]) < 0.1 * Lz0 and abs(L[1]) < 0.1 * Lz0, \
    f"off-axis momentum leaked: L=({L[0]:.3e},{L[1]:.3e},{L[2]:.3e})"
# oscillation: after the first minimum, the restoring couple swings Lz back up
assert recovery > 0.5 * Lz0, f"no torsional recovery after the dip (peak {recovery:.3f})"
print("GATHER PROBE PASS: spring reaction drives the spin, correct restoring sign, "
      f"no off-axis leak, clean torsional oscillation (Lz in [{Lz_min:.2f}, {Lz_max:.2f}])")
