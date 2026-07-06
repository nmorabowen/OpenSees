"""ZFP benchmark input A — nonlinear transient, many steps (ADR 08 config: "a
nonlinear transient with many steps").

2D J2-plasticity cantilever plate (30x8 quads, plane strain), multi-frequency
tip load driving partial yielding, 3000 steps. Written with the INSTALLED
Ladruno opensees (venv boot .pth) at default lossless f64 -> the oracle file
the re-encode spike compresses.

  python model_a_nonlinear2d.py <out_dir>
"""
import os
import sys

import numpy as np
import opensees as ops

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
os.makedirs(OUT, exist_ok=True)
target = os.path.join(OUT, "bench_a_nonlinear2d.ladruno")
if os.path.exists(target):
    os.remove(target)

NX, NY = 30, 8
L, H = 3000.0, 800.0          # mm
THICK = 100.0
dx, dy = L / NX, H / NY

E, NU = 200000.0, 0.3          # steel, N-mm
K = E / (3.0 * (1.0 - 2.0 * NU))
G = E / (2.0 * (1.0 + NU))
SIG0, SIGINF, DELTA, HISO = 250.0, 400.0, 10.0, 1000.0
RHO = 7.85e-9                  # tonne/mm^3

DT, NSTEP = 0.002, 3000


def nid(i, j):
    return 1 + i + j * (NX + 1)


ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
for j in range(NY + 1):
    for i in range(NX + 1):
        ops.node(nid(i, j), i * dx, j * dy)
for j in range(NY + 1):
    ops.fix(nid(0, j), 1, 1)

ops.nDMaterial("J2Plasticity", 1, K, G, SIG0, SIGINF, DELTA, HISO)
ops.nDMaterial("PlaneStrain", 2, 1)
etag = 1
for j in range(NY):
    for i in range(NX):
        ops.element("quad", etag, nid(i, j), nid(i + 1, j),
                    nid(i + 1, j + 1), nid(i, j + 1),
                    THICK, "PlaneStrain", 2, 0.0, RHO)
        etag += 1

# broadband tip load: 120 random-phase tones, flat 0.5-25 Hz (seeded) — a few
# discrete sines would be unrealistically codec-friendly vs real earthquake
# input (adversarial-review finding #7); this is the honest transient content
t = np.arange(NSTEP + 1) * DT
rng = np.random.default_rng(42)
freqs = np.linspace(0.5, 25.0, 120)
phases = rng.uniform(0.0, 2.0 * np.pi, freqs.size)
sig = np.zeros_like(t)
for f0, ph in zip(freqs, phases):
    sig += np.sin(2.0 * np.pi * f0 * t + ph)
sig *= np.minimum(t / 1.0, 1.0)              # 1 s ramp-in envelope
sig *= 1.5e5 / np.abs(sig).max()  # N per node -> base stress past SIG0 at peaks
ops.timeSeries("Path", 1, "-dt", DT, "-values", *sig.tolist())
ops.pattern("Plain", 1, 1)
for j in range(NY + 1):
    ops.load(nid(NX, j), 0.0, -1.0)

ops.recorder("ladruno", target,
             "-N", "displacement", "velocity", "acceleration", "reactionForce",
             "-E", "stress", "strain",
             "-precision", "f64", "-T", "dt", 0.0)

ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("UmfPack")
ops.test("NormDispIncr", 1.0e-7, 100)
ops.algorithm("Newton")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
ok = ops.analyze(NSTEP, DT)
ops.wipe()  # flush + close recorder

print("MODEL_A", "ok" if ok == 0 else f"FAILED rc={ok}",
      target, os.path.getsize(target) if os.path.exists(target) else -1)
if ok != 0:
    sys.exit(1)
