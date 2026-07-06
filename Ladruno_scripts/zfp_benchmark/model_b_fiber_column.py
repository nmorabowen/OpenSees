"""ZFP benchmark input B — fiber-section element history (ADR 08 config: "a
fiber-section element history").

3D RC cantilever column (5 forceBeamColumn elements, Concrete02+Steel02 fiber
section, 5 Lobatto IPs each), growing-amplitude cyclic lateral load into the
inelastic range, 2000 steps. Fiber stress/strain histories are the payload —
many components, hysteretic (kinked) time content: the pessimistic case for a
transform codec.

  python model_b_fiber_column.py <out_dir>
"""
import os
import sys

import numpy as np
import opensees as ops

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
target = os.path.join(OUT, "bench_b_fiber_column.ladruno")
if os.path.exists(target):
    os.remove(target)

NELE = 8
LCOL = 3000.0                 # mm
B = D = 500.0                 # square column
COVER = 40.0
DT, NSTEP = 0.005, 2000

ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 6)
for k in range(NELE + 1):
    ops.node(1 + k, 0.0, 0.0, k * LCOL / NELE)
    if k > 0:
        ops.mass(1 + k, 0.5, 0.5, 0.5, 1e-6, 1e-6, 1e-6)  # tonne
ops.fix(1, 1, 1, 1, 1, 1, 1)

# materials: unconfined+confined concrete, reinforcing steel (N-mm)
ops.uniaxialMaterial("Concrete02", 1, -30.0, -0.002, -6.0, -0.02, 0.1, 3.0, 1500.0)
ops.uniaxialMaterial("Concrete02", 2, -36.0, -0.0035, -18.0, -0.04, 0.1, 3.5, 1500.0)
ops.uniaxialMaterial("Steel02", 3, 420.0, 200000.0, 0.01, 18.0, 0.925, 0.15)

secT = 10
ops.section("Fiber", secT, "-GJ", 1.0e12)
c = COVER
# confined core + unconfined cover patches
ops.patch("rect", 2, 12, 12, -(B / 2 - c), -(D / 2 - c), B / 2 - c, D / 2 - c)
ops.patch("rect", 1, 14, 2, -B / 2, -D / 2, B / 2, -(D / 2 - c))
ops.patch("rect", 1, 14, 2, -B / 2, D / 2 - c, B / 2, D / 2)
ops.patch("rect", 1, 2, 10, -B / 2, -(D / 2 - c), -(B / 2 - c), D / 2 - c)
ops.patch("rect", 1, 2, 10, B / 2 - c, -(D / 2 - c), B / 2, D / 2 - c)
# 12 x φ20 bars (As = 314 mm^2)
AS = 314.0
ops.layer("straight", 3, 4, AS, -(B / 2 - c), -(D / 2 - c), B / 2 - c, -(D / 2 - c))
ops.layer("straight", 3, 4, AS, -(B / 2 - c), D / 2 - c, B / 2 - c, D / 2 - c)
ops.layer("straight", 3, 2, AS, -(B / 2 - c), -(D / 2 - c) / 3, -(B / 2 - c), (D / 2 - c) / 3)
ops.layer("straight", 3, 2, AS, B / 2 - c, -(D / 2 - c) / 3, B / 2 - c, (D / 2 - c) / 3)

ops.geomTransf("PDelta", 1, 0.0, 1.0, 0.0)
# displacement-based formulation: forceBeamColumn's internal state
# determination diverges (dW ~1e8) under concrete softening + P-Delta near
# peak push regardless of -iter budget; dispBeamColumn has no element-level
# Newton, produces the same section/fiber result tree, and 8 elements give
# adequate curvature resolution for a compression benchmark
ops.beamIntegration("Legendre", 1, secT, 3)
for k in range(NELE):
    ops.element("dispBeamColumn", 1 + k, 1 + k, 2 + k, 1, 1)

# gravity (kept in a constant pattern), then cyclic lateral history
ops.timeSeries("Constant", 9)
ops.pattern("Plain", 9, 9)
ops.load(NELE + 1, 0.0, 0.0, -1.5e6, 0.0, 0.0, 0.0)  # 1500 kN axial

t = np.arange(NSTEP + 1) * DT
amp = np.minimum(t / (0.6 * t[-1]), 1.0)              # ramp to full amplitude
# +-170 kN cyclic: past yield (check_inputs.py gate) without the collapse
# branch (+-200 kN failed on the first full-cycle reversal; +-140 kN stayed
# elastic, max fiber strain 1.55e-3 < eps_y)
sig = 1.7e5 * amp * np.sin(2 * np.pi * 0.5 * t)
ops.timeSeries("Path", 1, "-dt", DT, "-values", *sig.tolist())
ops.pattern("Plain", 1, 1)
ops.load(NELE + 1, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

ops.recorder("ladruno", target,
             "-N", "displacement", "acceleration",
             "-E", "localForce", "section.force", "section.deformation",
             "section.fiber.stress", "section.fiber.strain",
             "-precision", "f64", "-T", "dt", 0.0)

ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormDispIncr", 1.0e-6, 300)  # 1e-6 mm on a 3 m column; 1e-7 stalls
ops.algorithm("NewtonLineSearch")      # in the concrete-softening regime
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")

# adaptive time-stepping: halve dt on failure (recorder records per commit, so
# a locally refined dt is exactly what a production run would produce anyway)
ok = 0
tend = NSTEP * DT
tcur, step = 0.0, DT
while tcur < tend - 1e-9:
    rc = ops.analyze(1, step)
    if rc != 0:  # stiff-spot fallback before cutting dt
        ops.algorithm("KrylovNewton")
        rc = ops.analyze(1, step)
        ops.algorithm("NewtonLineSearch")
    if rc == 0:
        tcur += step
        step = min(DT, step * 2.0)
    else:
        step *= 0.5
        if step < DT / 64.0:
            ok = rc
            break
ops.wipe()

print("MODEL_B", "ok" if ok == 0 else f"FAILED rc={ok} at t={tcur:.3f}",
      target, os.path.getsize(target) if os.path.exists(target) else -1)
if ok != 0:
    sys.exit(1)
