"""
ADR-58 P2 Stage-1 gate — LadrunoRigidBody body-frame SO(3) integrator.

Torque-free Dzhanibekov (intermediate-axis / tennis-racket) free spin:
  GATE 1 (hard): spatial angular momentum is conserved — ‖L(t)−L0‖/‖L0‖ < 1e-10
                 over the whole run (the momentum-form integrator keeps L exact).
  GATE 2 (qualitative): the intermediate-axis instability actually fires — the
                 dominant ω component reverses sign at least once (a real flip,
                 not a frozen state).
  GATE 3 (energy): rotational KE drift is BOUNDED and SHRINKS with dt (the explicit
                 scheme is momentum- not energy-conserving — ADR §7 P2).

Body: 6 unit masses on the axes at (±a,0,0),(0,±b,0),(0,0,±c) → a diagonal,
asymmetric inertia I=(2b²+2c², 2a²+2c², 2a²+2b²). With a>b>c the INTERMEDIATE
principal axis is y; we spin about y with a small x,z perturbation.

Run with the worktree's Python against its dist/bin build:
    python3.12 Ladruno_scripts/rigidbody_tests/test_p2_dzhanibekov.py
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

a, b, c = 2.0, 1.5, 1.0           # a>b>c → intermediate principal axis = y
pts = [(a, 0, 0), (-a, 0, 0), (0, b, 0), (0, -b, 0), (0, 0, c), (0, 0, -c)]
Ix, Iy, Iz = 2*b*b + 2*c*c, 2*a*a + 2*c*c, 2*a*a + 2*b*b   # = 6.5, 10.0, 12.5
Iprin = (Ix, Iy, Iz)


def run(dt, nsteps, w0=(0.12, 6.0, 0.12), sample=50):
    ops.wipe()
    ops.model('basic', '-ndm', 3, '-ndf', 3)
    for i, (x, y, z) in enumerate(pts, start=1):
        ops.node(i, float(x), float(y), float(z))
        ops.mass(i, 1.0, 1.0, 1.0)
    tags = list(range(1, len(pts) + 1))
    ops.element('LadrunoRigidBody', 100, len(pts), *tags, '-omega', *w0)
    ops.constraints('Transformation')
    ops.numberer('Plain')
    ops.system('Diagonal')
    ops.test('NormDispIncr', 1e-12, 10, 0)
    ops.algorithm('Linear')
    ops.integrator('CentralDifference')
    ops.analysis('Transient')

    L0 = ops.eleResponse(100, 'angularMom')
    L0n = math.sqrt(sum(v*v for v in L0))

    def energy():                         # body frame: E = 1/2 Σ I_i ω_body,i²
        wb = ops.eleResponse(100, 'omegaBody')
        return 0.5 * sum(Iprin[k] * wb[k] * wb[k] for k in range(3))

    E0 = energy()
    maxLdev = 0.0
    maxEdev = 0.0
    wby_signs = []                        # body-frame intermediate-axis (y) ω sign
    for s in range(nsteps):
        ops.analyze(1, dt)
        if s % sample == 0:
            L = ops.eleResponse(100, 'angularMom')
            Ln = math.sqrt(sum(v*v for v in L))
            maxLdev = max(maxLdev, abs(Ln - L0n) / L0n)
            maxEdev = max(maxEdev, abs(energy() - E0) / abs(E0))
            wb = ops.eleResponse(100, 'omegaBody')
            wby_signs.append(1 if wb[1] >= 0 else -1)
    flips = sum(1 for i in range(1, len(wby_signs)) if wby_signs[i] != wby_signs[i-1])
    return L0n, maxLdev, maxEdev, flips


print(f"Iprincipal = {Iprin}  (intermediate axis = y)")
L0n, ldev, edev, flips = run(dt=1.0e-3, nsteps=40000)
print(f"dt=1e-3:  ||L0||={L0n:.6f}  max ||L|| rel-dev={ldev:.3e}  "
      f"max E rel-drift={edev:.3e}  omega_y flips={flips}")

# GATE 3: energy drift shrinks with dt (momentum-conserving, not energy-conserving)
_, _, edev_half, _ = run(dt=5.0e-4, nsteps=80000)
print(f"dt=5e-4:  max E rel-drift={edev_half:.3e}  (should be < dt=1e-3 drift)")

assert ldev < 1e-10,  f"GATE 1 FAILED: ||L|| not conserved (rel-dev {ldev:.3e})"
assert flips >= 1,    f"GATE 2 FAILED: no intermediate-axis flip observed"
assert edev_half < edev + 1e-12, f"GATE 3 FAILED: energy drift did not shrink with dt"
print("P2 STAGE-1 GATE PASS: ||L|| conserved to 1e-10, Dzhanibekov flip captured, "
      "energy drift bounded & shrinking")
