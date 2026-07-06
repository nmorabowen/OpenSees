"""ZFP benchmark input C — large nodal field (ADR 08 config: "a large nodal
field").

3D elastic soil block (16x16x8 stdBrick grid, 2601 nodes), Ricker-pulse base
excitation (UniformExcitation), 600 steps: wave propagation = the
spatially-correlated field ZFP's block transform is designed for. This is the
optimistic case, and the one the component-separated/time-major layout
prototype (config 5) is judged on.

  python model_c_brick_field.py <out_dir>
"""
import os
import sys

import numpy as np
import opensees as ops

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
target = os.path.join(OUT, "bench_c_brick_field.ladruno")
if os.path.exists(target):
    os.remove(target)

NX, NY, NZ = 16, 16, 8
LX, LY, LZ = 32000.0, 32000.0, 16000.0  # mm (32 m x 32 m x 16 m block)
dx, dy, dz = LX / NX, LY / NY, LZ / NZ

E, NU, RHO = 30000.0, 0.3, 2.4e-9   # stiff soil / soft rock, N-mm-tonne
DT, NSTEP = 0.001, 600


def nid(i, j, k):
    return 1 + i + j * (NX + 1) + k * (NX + 1) * (NY + 1)


ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
for k in range(NZ + 1):
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.node(nid(i, j, k), i * dx, j * dy, k * dz)
for j in range(NY + 1):
    for i in range(NX + 1):
        ops.fix(nid(i, j, 0), 1, 1, 1)

ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
etag = 1
for k in range(NZ):
    for j in range(NY):
        for i in range(NX):
            ops.element("stdBrick", etag,
                        nid(i, j, k), nid(i + 1, j, k),
                        nid(i + 1, j + 1, k), nid(i, j + 1, k),
                        nid(i, j, k + 1), nid(i + 1, j, k + 1),
                        nid(i + 1, j + 1, k + 1), nid(i, j + 1, k + 1), 1)
            etag += 1

# Ricker wavelet base acceleration, fp = 4 Hz
t = np.arange(NSTEP + 1) * DT
fp, t0 = 4.0, 0.25
a = (np.pi * fp * (t - t0)) ** 2
acc = 2000.0 * (1.0 - 2.0 * a) * np.exp(-a)   # mm/s^2 peak ~0.2 g
ops.timeSeries("Path", 1, "-dt", DT, "-values", *acc.tolist())
ops.pattern("UniformExcitation", 1, 1, "-accel", 1)

ops.recorder("ladruno", target,
             "-N", "displacement", "velocity", "acceleration",
             "-precision", "f64", "-T", "dt", 0.0)

ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
ok = ops.analyze(NSTEP, DT)
ops.wipe()

print("MODEL_C", "ok" if ok == 0 else f"FAILED rc={ok}",
      target, os.path.getsize(target) if os.path.exists(target) else -1)
if ok != 0:
    sys.exit(1)
