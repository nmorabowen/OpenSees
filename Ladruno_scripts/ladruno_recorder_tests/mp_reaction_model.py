"""Partitioned-reaction gate (WP-126) -- model. The WP-126 reproduction: a SUPPORT node
shared by two partitions. Checker: mp_reaction_check.py.

Two-truss 2D model, statically determinate, total reaction known by hand:

  node 1 (0,0)  fixed (x, y)           <- SHARED by both ranks: the reaction we test
  node 2 (0,1)  fixed x; Fy = -10      truss 1: 1-2 (vertical)      -> rank 0
  node 3 (1,1)  fixed x; Fy = -20      truss 2: 1-3 (45 degrees)    -> rank 1

  total reaction at node 1 (serial):  (Rx, Ry) = (-20, +30) up to sign convention of Rx
  expected partials:  rank 0 -> (0, 10)   rank 1 -> (-20, 20)

Per rank it prints ops.nodeReaction(1) and writes (a) a Ladruno streaming recorder
(reactionForce), (b) a Ladruno -envelope recorder, (c) a vanilla Node recorder, into
OUT/<tag>.part-<rank>.* (Ladruno appends .part-N itself when partitioned).

    serial:  python mp_reaction_model.py serial <dir_with_opensees_pyd> <out_dir>
    MP:      mpiexec -n 2 python mp_reaction_model.py mp <dir_with_openseesmp_pyd> <out_dir>
"""
from __future__ import annotations

import os
import sys

MODE, MODDIR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
_IMPI = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"
_LIBFABRIC = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"
_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
for d in (MODDIR, _IMPI, _LIBFABRIC, _MKL, _ICOMP):
    if os.path.isdir(d):
        os.add_dll_directory(d)
os.environ["PATH"] = os.pathsep.join([_IMPI, _LIBFABRIC, MODDIR, os.environ.get("PATH", "")])
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
sys.path.insert(0, MODDIR)

if MODE == "mp":
    import openseesmp as ops
    pid, nproc = ops.getPID(), ops.getNP()
    ops.start()
else:
    import opensees as ops
    pid, nproc = 0, 1

ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 2)
ops.uniaxialMaterial("Elastic", 1, 1.0e6)
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)

ops.node(1, 0.0, 0.0)          # the shared support
ops.fix(1, 1, 1)
if pid == 0:                   # rank 0 (or the serial run) holds truss 1
    ops.node(2, 0.0, 1.0)
    ops.fix(2, 1, 0)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.load(2, 0.0, -10.0)
if pid == 1 or nproc == 1:     # rank 1 (or the serial run) holds truss 2
    ops.node(3, 1.0, 1.0)
    ops.fix(3, 1, 0)
    ops.element("Truss", 2, 1, 3, 1.0, 1)
    ops.load(3, 0.0, -20.0)

tag = f"{MODE}_np{nproc}"
ops.recorder("ladruno", os.path.join(OUT, f"{tag}_stream.ladruno"), "-N", "displacement", "reactionForce")
ops.recorder("ladruno", os.path.join(OUT, f"{tag}_env.ladruno"), "-N", "displacement", "reactionForce", "-envelope")
ops.recorder("Node", "-file", os.path.join(OUT, f"{tag}_node.rank{pid}.txt"), "-node", 1,
             "-dof", 1, 2, "reaction")

ops.constraints("Transformation")
if MODE == "mp":
    ops.numberer("ParallelPlain")
    ops.system("Mumps")
else:
    ops.numberer("Plain")
    ops.system("FullGeneral")
ops.test("NormDispIncr", 1.0e-12, 10)
ops.algorithm("Newton")
ops.integrator("LoadControl", 0.5)
ops.analysis("Static")
for step in range(2):          # two load steps: 50 %, 100 %
    rc = ops.analyze(1)
    ops.reactions()
    r = ops.nodeReaction(1)
    print(f"RANK {pid}/{nproc} step {step + 1} rc={rc} nodeReaction(1) = ({r[0]:+.12g}, {r[1]:+.12g})",
          flush=True)
ops.wipe()                     # flush + close recorders
print(f"RANK {pid}/{nproc} DONE", flush=True)
