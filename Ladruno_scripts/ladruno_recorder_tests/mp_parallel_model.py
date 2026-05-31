"""Parallel (openseesmp / interpreter-per-rank) MPCO_Ladruno output gate -- model.

Each MPI rank builds its own piece of a 3-truss frame (the canonical OpenSeesMP
manual-partition example) and adds the SAME recorder command. The point of the
test: before the parallel fix every rank wrote the one file "<name>.ladruno" and
they corrupted each other; now each rank must write its own
"mp_out.part-<rank>.ladruno" (the per-partition contract), with INFO/PARTITION_ID
== its rank and INFO/NUM_PARTITIONS == np.

Run under mpiexec with the MP build python:
    mpiexec -n <N> python mp_parallel_model.py <dir_with_openseesmp_pyd> <mkl_dir> [out_dir]
"""

from __future__ import annotations

import os
import sys

MODDIR = sys.argv[1]                       # dir holding openseesmp.pyd
MKLDIR = sys.argv[2] if len(sys.argv) > 2 else MODDIR
OUT = sys.argv[3] if len(sys.argv) > 3 else MODDIR

# openseesmp.pyd is dynamically linked against Intel MPI (impi.dll) + libfabric +
# the Intel OpenMP runtime. Python resolves a .pyd's dependencies from the .pyd dir
# and add_dll_directory dirs (NOT PATH), and libfabric additionally needs a
# process-local PATH (matches the installer boot module). Register every oneAPI bin
# that exists so the import succeeds under mpiexec.
_IMPI = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"
_LIBFABRIC = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"
_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
for d in (MODDIR, MKLDIR, _IMPI, _LIBFABRIC, _MKL, _ICOMP):
    if d and os.path.isdir(d):
        os.add_dll_directory(d)
os.environ["PATH"] = os.pathsep.join(
    [p for p in (_IMPI, _LIBFABRIC, MODDIR, MKLDIR, os.environ.get("PATH", "")) if p])
sys.path.insert(0, MODDIR)

import openseesmp as ops  # noqa: E402

pid = ops.getPID()
nproc = ops.getNP()
ops.start()
if nproc < 2:
    print("WARN: need >=2 ranks to exercise partitioned output")

ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 3000.0)

if pid == 0:
    ops.node(1, 0.0, 0.0)
    ops.node(4, 72.0, 96.0)
    ops.fix(1, 1, 1)
    ops.mass(4, 100.0, 100.0)
    ops.element("Truss", 1, 1, 4, 10.0, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 100.0, -50.0)
else:
    # workers each hold one diagonal truss to the shared apex node 4
    base = 1 + pid  # distinct node tags per rank
    ops.node(base * 10 + 2, 144.0 + pid, 0.0)
    ops.node(4, 72.0, 96.0)
    ops.fix(base * 10 + 2, 1, 1)
    ops.mass(4, 100.0, 100.0)
    ops.element("Truss", 100 + pid, base * 10 + 2, 4, 5.0, 1)

out = os.path.join(OUT, "mp_out.ladruno")
ops.recorder("mpcoLadruno", out, "-N", "displacement", "-kind", "transient",
             "-T", "dt", 0.0)

ops.constraints("Transformation")
ops.numberer("ParallelPlain")
ops.test("NormDispIncr", 1e-6, 6, 0)
ops.algorithm("Linear")
ops.system("MPIDiagonal")
ops.integrator("ExplicitDifference")
ops.analysis("Transient")
for _ in range(10):
    ops.analyze(1, 1e-6)
ops.wipe()  # flush + close the recorder file

print(f"RANK {pid}/{nproc} DONE")
