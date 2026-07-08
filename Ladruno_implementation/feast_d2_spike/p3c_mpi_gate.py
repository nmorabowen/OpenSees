"""ADR-43 P3c-MPI gate — distributed FEAST-RCI inner solve (L3-only, R0).

Under openseesmp with np>1, `eigen -feast ... -rci` routes its contour solves
through LadrunoDistBlockZKernel: the 2n block-real (zM-K) system factored/solved
by DISTRIBUTED MUMPS (SYM=2) across the ranks, the dfeast_srci outer loop run
REPLICATED on every rank, the solution broadcast so the RCI stays in lockstep.

The differential (all in ONE mpiexec run, identical replicated frame model):
  G1  distributed `-rci` spectrum == the driver-path `-feast` (dfeast_scsrgv,
      per-rank serial) spectrum: lambda rel <= 1e-8, MAC >= 0.999 per mode.
      The driver path is the in-run oracle (serial -rci == driver was pinned by
      the serial gate #530); so G1 certifies the distributed inner solve.
  G2  LOCKSTEP: every rank returns the identical `-rci` spectrum (a divergent
      replicated RCI would deadlock or return per-rank-different eigenvalues).
  G3  Phi^T M Phi == I for the distributed modes (mass-normalization survives
      the distributed solve + broadcast).
Run at np=2 AND np=4 (different triplet partitions of the same 2n system).

Run:  python p3c_mpi_gate.py                 (orchestrator: np=2 and np=4)
      python p3c_mpi_gate.py --driver        (internal, under mpiexec)
"""
import json
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))
DIST_MP = os.path.join(WT, "dist", "openseesmp")
PYTHON = sys.executable

TWO_PI = 2.0 * math.pi
NX, NY = 3, 10                      # frame grid (~132 DOF), matches test_feastEigen
E, A, IZ, MASS = 30.0e9, 0.09, 6.75e-4, 5.0e3


def _hz(lam):
    return math.sqrt(lam) / TWO_PI


def _build_frame(ops):
    """Replicated 2D elastic frame — every rank builds the FULL model (the
    replicated-interpreter assumption the L3 distribution rests on)."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 3)
    tag = 0
    for j in range(NY + 1):
        for i in range(NX + 1):
            tag += 1
            ops.node(tag, i * 5.0, j * 3.0)
            if j == 0:
                ops.fix(tag, 1, 1, 1)
            else:
                ops.mass(tag, MASS, MASS, 0.0)
    ops.geomTransf("Linear", 1)
    ele = 0

    def nid(i, j):
        return j * (NX + 1) + i + 1

    for j in range(1, NY + 1):
        for i in range(NX + 1):
            ele += 1
            ops.element("elasticBeamColumn", ele, nid(i, j - 1), nid(i, j),
                        A, E, IZ, 1)
    for j in range(1, NY + 1):
        for i in range(NX):
            ele += 1
            ops.element("elasticBeamColumn", ele, nid(i, j), nid(i + 1, j),
                        A, E, IZ, 1)


def _free_nodes():
    return [j * (NX + 1) + i + 1 for j in range(1, NY + 1) for i in range(NX + 1)]


def _shapes(ops, nmodes):
    """Stacked (ux, uy) per free node for each mode -> list[list[float]]."""
    nodes = _free_nodes()
    out = []
    for k in range(1, nmodes + 1):
        v = []
        for n in nodes:
            ev = ops.nodeEigenvector(n, k)
            v.extend(ev[:2])
        out.append(v)
    return out


def driver(outdir):
    """Under mpiexec: driver-path -feast and distributed -rci on one model.
    Each rank writes its result to outdir/rank_<world>.json (files, not stdout,
    so the large mode-shape payloads can't interleave across ranks)."""
    sys.path.insert(0, DIST_MP)
    os.add_dll_directory(DIST_MP)
    os.environ["LADRUNO_OPENSEES_QUIET"] = "1"
    # Determinism + no thread oversubscription under MPI (ADR R2): the
    # replicated dfeast_srci must issue the SAME ijob sequence on every rank
    # for the collective inner solves to line up — single-threaded MKL LAPACK
    # keeps the internal reduced-eig arithmetic bit-identical across ranks.
    os.environ["MKL_NUM_THREADS"] = "1"
    os.environ["OMP_NUM_THREADS"] = "1"
    import openseesmp as ops

    world = ops.getPID()
    npr = ops.getNP()

    # band enclosing the first 8 modes (from a per-rank ARPACK reference)
    _build_frame(ops)
    lam_ref = list(ops.eigen(10))
    f_hi = 0.5 * (_hz(lam_ref[7]) + _hz(lam_ref[8]))

    # driver path (dfeast_scsrgv, per-rank serial) — the in-run oracle
    _build_frame(ops)
    lam_drv = list(ops.eigen("-feast", 0.0, f_hi))
    phi_drv = _shapes(ops, len(lam_drv))

    # distributed -rci (np>1 -> LadrunoDistBlockZKernel across the ranks)
    _build_frame(ops)
    lam_rci = list(ops.eigen("-feast", 0.0, f_hi, "-rci"))

    out = dict(world=world, npr=npr, f_hi=f_hi, lam_rci=lam_rci)
    # only the host ships the heavy mode-shape payload (G1 MAC + G3 Gram);
    # the other ranks ship just lam_rci for the G2 lockstep check
    if world == 0:
        out.update(lam_drv=lam_drv, phi_drv=phi_drv,
                   phi_rci=_shapes(ops, len(lam_rci)))
    with open(os.path.join(outdir, f"rank_{world}.json"), "w") as fh:
        json.dump(out, fh)


def _mac(a, b):
    import numpy as np
    a = np.asarray(a); b = np.asarray(b)
    return float((a @ b) ** 2 / ((a @ a) * (b @ b)))


def _mpi_run(nranks):
    import tempfile
    mpiexec = os.path.join(DIST_MP, "mpiexec.exe")
    env = dict(os.environ)
    env["PATH"] = DIST_MP + os.pathsep + env.get("PATH", "")
    env["MKL_NUM_THREADS"] = "1"
    env["OMP_NUM_THREADS"] = "1"
    outdir = tempfile.mkdtemp(prefix=f"p3c_n{nranks}_")
    cmd = [mpiexec, "-n", str(nranks), PYTHON, "-S",
           os.path.abspath(__file__), "--driver", outdir]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=600)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"mpi run n={nranks} failed rc={r.returncode}")
    res = []
    for w in range(nranks):
        fp = os.path.join(outdir, f"rank_{w}.json")
        if not os.path.exists(fp):
            print(r.stdout); print(r.stderr, file=sys.stderr)
            raise RuntimeError(f"rank {w} produced no result file")
        with open(fp) as fh:
            res.append(json.load(fh))
    res.sort(key=lambda d: d["world"])
    return res


def _check(nranks):
    import numpy as np
    res = _mpi_run(nranks)
    r0 = next(r for r in res if r["world"] == 0)
    lam_drv = np.asarray(r0["lam_drv"])
    lam_rci = np.asarray(r0["lam_rci"])
    assert lam_rci.shape == lam_drv.shape and lam_rci.size >= 1, \
        f"n={nranks}: shape mismatch drv={lam_drv.shape} rci={lam_rci.shape}"

    # G1: distributed -rci == driver path (lambda + MAC)
    np.testing.assert_allclose(lam_rci, lam_drv, rtol=1e-8,
                               err_msg=f"n={nranks} G1 lambda")
    for k in range(lam_rci.size):
        mac = _mac(r0["phi_rci"][k], r0["phi_drv"][k])
        assert mac >= 0.999, f"n={nranks} G1 MAC mode {k+1}: {mac:.6f}"

    # G2: lockstep — every rank returns the identical -rci spectrum
    for r in res:
        if r["world"] == 0:
            continue
        np.testing.assert_allclose(
            np.asarray(r["lam_rci"]), lam_rci, rtol=1e-12,
            err_msg=f"n={nranks} G2 lockstep: rank {r['world']} diverged")

    # G3: Phi^T M Phi == I (nodal masses only -> exact reduced Gram)
    phi = np.asarray(r0["phi_rci"])         # p x (2*nfree)
    gram = MASS * (phi @ phi.T)
    np.testing.assert_allclose(gram, np.eye(lam_rci.size), atol=1e-8,
                               err_msg=f"n={nranks} G3 Phi^T M Phi")

    print(f"  n={nranks}: G1+G2+G3 PASS — {lam_rci.size} modes, "
          f"max |drv-rci|/|drv| = "
          f"{np.max(np.abs(lam_rci-lam_drv)/np.abs(lam_drv)):.2e}, "
          f"all {nranks} ranks agree")


def main():
    if "--driver" in sys.argv:
        driver(sys.argv[sys.argv.index("--driver") + 1])
        return
    print("ADR-43 P3c-MPI gate: distributed FEAST-RCI inner solve")
    for nranks in (2, 4):
        _check(nranks)
    print("ALL GATES PASS")


if __name__ == "__main__":
    main()
