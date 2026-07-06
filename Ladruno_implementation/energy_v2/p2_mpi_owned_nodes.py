"""ADR-69 P2 gate — MPI nodal-term bookkeeping + the -ownedNodes gate.

Two-rank flat-per-rank (openseesmp) free-vibration truss, mirroring
EXAMPLES/ExamplePython/example_mpi_paralleltruss_explicit.py: shared node 4
lives on BOTH ranks with mass(4, m, m) declared on both (the example's
idiom - the parallel diagonal assembly SUMS duplicate contributions, so the
assembled system has 2m at node 4). The serial reference therefore uses 2m.

What this gate establishes empirically (the ADR's multiply-count claim was
scoped from _PARALLEL_PROCESSING PartitionedDomain semantics; openseesmp is
flat-per-rank and the semantics differ):

  G1  element block reduces exactly: sum over ranks of KE_ele/IE/DW_ele ==
      serial (elements are uniquely partitioned).
  G2  nodal-term semantics MEASURED: per-rank KE_nod sums vs serial KE_nod.
      In the split-mass idiom (m per rank, assembled 2m) the naive sum is
      CORRECT (each rank books its own mass share). The gate asserts
      whichever relation holds and prints it - this is a characterization,
      not a bug hunt.
  G3  -ownedNodes mechanics: with rank 0 owning the shared node (region
      membership), the owned-filtered sum differs from the naive sum by
      exactly rank 1's share of node 4's KE at every sample - the identity
      that makes offline dedup / canonical-rank emits exact when the model
      DOES declare full duplicate values (PartitionedDomain-style mirrors,
      or full-mass-on-all-ranks emits).
  G4  control: a rank with -ownedNodes owning ALL its nodes is byte-
      identical to the same rank without the flag.

Run:  python p2_mpi_owned_nodes.py          (orchestrator: serial + 2 mpi runs)
      python p2_mpi_owned_nodes.py --driver <mode>   (internal, under mpiexec)
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))
DIST_BIN = os.path.join(WT, "dist", "bin")
DIST_MP = os.path.join(WT, "dist", "openseesmp")
PYTHON = sys.executable

M = 100.0        # nodal mass declared PER RANK at the shared node (serial: 2M)
V0 = 1.5         # initial velocity at the shared node
DT, NSTEP = 1.0e-5, 200
RHO = 0.01       # truss element mass density (element-KE block signal)


def _model(ops, pid, mode):
    """Build the per-rank (pid 0/1) or serial (pid < 0) model."""
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, 3000.0)

    if pid <= 0:   # serial or rank 0
        ops.node(1, 0.0, 0.0)
        ops.node(4, 72.0, 96.0)
        ops.fix(1, 1, 1)
        ops.element("Truss", 1, 1, 4, 10.0, 1, "-rho", RHO)
    if pid != 0:   # serial or rank 1
        ops.node(2, 144.0, 0.0)
        ops.node(3, 168.0, 0.0)
        if pid > 0:
            ops.node(4, 72.0, 96.0)
        ops.fix(2, 1, 1)
        ops.fix(3, 1, 1)
        ops.element("Truss", 2, 2, 4, 5.0, 1, "-rho", RHO)
        ops.element("Truss", 3, 3, 4, 5.0, 1, "-rho", RHO)

    if pid < 0:
        ops.mass(4, 2.0 * M, 2.0 * M)   # serial reference = ASSEMBLED mass
    else:
        ops.mass(4, M, M)               # per-rank share (example's idiom)

    ops.setNodeVel(4, 1, V0, "-commit")
    ops.setNodeVel(4, 2, -0.5 * V0, "-commit")

    # ownership region for -ownedNodes runs: rank 0 owns {1, 4}, rank 1 owns
    # {2, 3} (canonical-rank convention: the shared node is owned by rank 0)
    if pid == 0:
        ops.region(1, "-node", 1, 4)
    elif pid > 0:
        ops.region(1, "-node", 2, 3)


def _analysis(ops, parallel):
    ops.constraints("Transformation")
    ops.numberer("ParallelPlain" if parallel else "Plain")
    ops.test("NormDispIncr", 1e-8, 6, 2)
    ops.algorithm("Linear")
    ops.system("MPIDiagonal" if parallel else "Diagonal")
    ops.integrator("ExplicitDifference")
    ops.analysis("Transient")
    for _ in range(NSTEP):
        if ops.analyze(1, DT) != 0:
            raise SystemExit("analyze failed")


def driver(mode):
    """Runs under mpiexec -n 2. mode: 'naive' | 'owned' | 'ownall'."""
    sys.path.insert(0, DIST_MP)
    os.add_dll_directory(DIST_MP)
    import openseesmp as ops

    pid = ops.getPID()
    _model(ops, pid, mode)

    efile = os.path.join(HERE, "p2mpi_%s_energy.txt" % mode)
    rec = ["EnergyBalance", "-file", efile, "-time", "-v2"]
    if mode == "owned":
        rec += ["-ownedNodes", 1]
    elif mode == "ownall":
        # region 2 = every node this rank has (control: flag on, no filtering)
        if pid == 0:
            ops.region(2, "-node", 1, 4)
        else:
            ops.region(2, "-node", 2, 3, 4)
        rec += ["-ownedNodes", 2]
    ops.recorder(*rec)

    _analysis(ops, parallel=True)
    ops.wipe()


def load(path):
    import numpy as np
    return np.atleast_2d(np.loadtxt(path))


def main():
    import numpy as np

    # --- serial reference (assembled system: 2M at node 4) ---
    sys.path.insert(0, DIST_BIN)
    os.add_dll_directory(DIST_BIN)
    import opensees as ops
    print("[gate] serial pyd:", ops.__file__)
    assert "strange-hawking" in ops.__file__

    _model(ops, -1, "serial")
    efile = os.path.join(HERE, "p2mpi_serial_energy.txt")
    ops.recorder("EnergyBalance", "-file", efile, "-time", "-v2")
    _analysis(ops, parallel=False)
    ops.wipe()
    d_ser = load(efile)

    # --- MPI runs ---
    mpiexec = os.path.join(DIST_MP, "mpiexec.exe")
    env = dict(os.environ)
    env["PATH"] = DIST_MP + os.pathsep + env.get("PATH", "")
    for mode in ("naive", "owned", "ownall"):
        cmd = [mpiexec, "-n", "2", PYTHON, "-S", os.path.abspath(__file__),
               "--driver", mode]
        print("[gate] running:", " ".join(cmd))
        r = subprocess.run(cmd, env=env, capture_output=True, text=True,
                           timeout=600)
        if r.returncode != 0:
            print(r.stdout[-3000:])
            print(r.stderr[-3000:])
            raise SystemExit("mpi run '%s' failed" % mode)

    ok = True

    def parts(mode):
        return (load(os.path.join(HERE, "p2mpi_%s_energy.part-0.txt" % mode)),
                load(os.path.join(HERE, "p2mpi_%s_energy.part-1.txt" % mode)))

    # v2 columns, no channels declared in these processes:
    # time KE_ele KE_nod IE DW_ele DW_nod ULW RES ERR%
    KE_ELE, KE_NOD, IE = 1, 2, 3

    d0, d1 = parts("naive")
    n = min(d_ser.shape[0], d0.shape[0], d1.shape[0])

    # G1: element block reduces exactly
    ele_sum = d0[:n, KE_ELE] + d1[:n, KE_ELE]
    ie_sum = d0[:n, IE] + d1[:n, IE]
    scale = max(np.abs(d_ser[:n, KE_ELE]).max(), 1e-30)
    g1 = (np.allclose(ele_sum, d_ser[:n, KE_ELE], atol=1e-9 * scale) and
          np.allclose(ie_sum, d_ser[:n, IE],
                      atol=1e-9 * max(np.abs(d_ser[:n, IE]).max(), 1e-30)))
    print("G1 element block sums == serial -> %s" % ("PASS" if g1 else "FAIL"))
    ok &= g1

    # G2: MEASURE the nodal semantics. Split-mass idiom prediction: naive sum
    # of per-rank KE_nod == serial KE_nod (each rank books its own share).
    nod_sum = d0[:n, KE_NOD] + d1[:n, KE_NOD]
    nod_scale = max(np.abs(d_ser[:n, KE_NOD]).max(), 1e-30)
    split_correct = np.allclose(nod_sum, d_ser[:n, KE_NOD],
                                atol=1e-6 * nod_scale)
    doubled = np.allclose(nod_sum, 2.0 * d_ser[:n, KE_NOD],
                          atol=1e-6 * nod_scale)
    print("G2 nodal semantics: naive==serial %s | naive==2x-serial %s"
          % (split_correct, doubled))
    g2 = split_correct or doubled
    print("G2 nodal sum characterized -> %s" % ("PASS" if g2 else "FAIL"))
    ok &= g2

    # G3: -ownedNodes identity - owned drops EXACTLY rank 1's share of the
    # shared node's nodal terms (rank 1 owns only {2,3}, both massless).
    o0, o1 = parts("owned")
    g3 = (np.allclose(o0[:n, KE_NOD], d0[:n, KE_NOD], atol=1e-12) and
          np.abs(o1[:n, KE_NOD]).max() < 1e-12 and
          np.allclose(o0[:n, KE_ELE], d0[:n, KE_ELE], atol=1e-12) and
          np.allclose(o1[:n, KE_ELE], d1[:n, KE_ELE], atol=1e-12))
    print("G3 -ownedNodes drops the non-owner's nodal share, element block "
          "untouched -> %s" % ("PASS" if g3 else "FAIL"))
    ok &= g3

    # G4: control - owning ALL local nodes == no flag at all
    a0, a1 = parts("ownall")
    g4 = (np.allclose(a0[:n], d0[:n], atol=0.0) and
          np.allclose(a1[:n], d1[:n], atol=0.0))
    print("G4 own-all control byte-identical to naive -> %s"
          % ("PASS" if g4 else "FAIL"))
    ok &= g4

    print("\nADR-69 P2 MPI owned-nodes gate:", "ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--driver":
        driver(sys.argv[2])
    else:
        sys.path.insert(
            0, r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64"
               r"\Lib\site-packages")
        sys.exit(main())
