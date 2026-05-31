"""Energy-balance result type — model runner (D8 verification).

Runs ONE dynamic model that emits the new ON_DOMAIN / ON_REGIONS energy channel
from `recorder ladruno -G energy <regionTag>`, plus the plain-text
`EnergyBalance` sidecar on the SAME model so the two can be cross-checked (they
share EnergyBalanceKernel.h, so the whole-model values must agree).

Physics: an undamped SDOF axial spring-mass with a SUDDENLY-APPLIED constant
load. The system is conservative, so all external work goes to kinetic + strain
energy every step:  ULW = KE + IE  (DW = 0)  =>  RES = ULW-(KE+IE+DW) ~ 0.
A clean closed energy balance (small |ERR%|) is the acceptance check.

Run with the BUILD python (no boot .pth), pointing sys.path at a copy of the
fresh OpenSeesPy.dll -> opensees.pyd (same recipe as parity_model.py):
    python energy_model.py <dir_with_opensees_pyd> <out_dir>
"""

from __future__ import annotations

import math
import os
import sys

TMP = sys.argv[1]                      # dir holding the fresh opensees.pyd
OUT = sys.argv[2] if len(sys.argv) > 2 else TMP
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"

os.add_dll_directory(DIST)             # MKL etc. (HDF5 is static-linked)
sys.path.insert(0, TMP)

import opensees as ops  # noqa: E402

out_ladruno = os.path.join(OUT, "energy.ladruno")
out_sidecar = os.path.join(OUT, "energy_sidecar.txt")
for p in (out_ladruno, out_sidecar):
    if os.path.exists(p):
        os.remove(p)

# ---- model: undamped axial spring-mass (SDOF in x) ----
E = 1.0e3
A = 1.0
L = 1.0
m = 10.0
P = 1.0
k = E * A / L
w = math.sqrt(k / m)
T = 2.0 * math.pi / w
dt = T / 100.0
nsteps = 200          # two natural periods

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.node(1, 0.0, 0.0)
ops.node(2, L, 0.0)
ops.fix(1, 1, 1)
ops.fix(2, 0, 1)                       # node 2 free in x only
ops.uniaxialMaterial("Elastic", 1, E)
ops.element("truss", 1, 1, 2, A, 1)
ops.mass(2, m, m)                      # y mass is inert (dof fixed)

# region 1 = element 1 (MeshRegion derives its connected nodes 1,2), so
# ON_REGIONS has real content and SET_1 self-describes both nodes+elements
ops.region(1, "-ele", 1)

# new energy channel: whole-model (ON_DOMAIN) + region 1 (ON_REGIONS)
ops.recorder("ladruno", out_ladruno, "-G", "energy", 1, "-T", "dt", 0.0)
# plain-text sidecar on the same model (shared-kernel cross-check)
ops.recorder("EnergyBalance", "-file", out_sidecar, "-time", "-region", 1)

# suddenly-applied constant load (conservative forced response)
ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
ops.load(2, P, 0.0)

ops.constraints("Transformation")
ops.numberer("RCM")
ops.system("BandSPD")
ops.test("NormDispIncr", 1.0e-12, 10)
ops.algorithm("Linear")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")
ops.analyze(nsteps, dt)
ops.wipe()                            # flush + close both recorders

print("LADRUNO:", out_ladruno, os.path.exists(out_ladruno))
print("SIDECAR:", out_sidecar, os.path.exists(out_sidecar))
print(f"PARAMS w={w:.6g} T={T:.6g} dt={dt:.6g} nsteps={nsteps}")
print("ENERGY_MODEL OK")
