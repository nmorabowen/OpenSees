"""ADR-43 L2 justification profiler — strong-scaling of the L3-only distributed
FEAST-RCI inner solve on a 3D solid box.

THE QUESTION (ADR 43 R0 / umbrella ADR 45): L2 (quadrature-parallel: split the
world into N_q node-groups, run the L3 kernel in each, Allreduce the projector)
is worth building ONLY if the L3-only path SATURATES below the rank budget — i.e.
adding ranks to ONE distributed solve stops helping while N_q(~<=8, ~x4 after
half-contour) independent solves sit idle.

METHOD: time the distributed `eigen -feast 0 fmax -rci` call (which under np>1
routes every contour solve through LadrunoDistBlockZKernel = distributed dmumps
SYM=2 on the 2n block system) at np = 1,2,4,8,... on ONE fixed replicated model.
  - T(np) that keeps dropping ~linearly to the max np  => L3 NOT saturated
    => L2 NOT justified yet (spend ranks on L3).
  - T(np) that FLATTENS at some R* < budget  => headroom L2 could capture:
    a ~min(N_q_eff, budget/R*)x throughput multiplier L3 can't deliver.

np=1 is the serial PARDISO kernel (byte-identical fallback), the L3=1 baseline.

Run:  python l2_profile.py [ne] [np_csv]      orchestrator
        e.g.  python l2_profile.py 30 1,2,4,8
      python l2_profile.py --driver <ne> <want_modes> <outdir>   (under mpiexec)
"""
import json
import os
import re
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))
# Module dir (where openseesmp.pyd/.so lives) + launcher are platform-specific:
#   Windows: dist/openseesmp/{openseesmp.pyd, mpiexec.exe}
#   Linux (esmeralda): build/Release/openseesmp.so, mpirun on PATH
# Override via env: L2_MODDIR (module dir), L2_MPIEXEC (launcher), L2_LDPATH
# (extra LD_LIBRARY_PATH dirs for the Linux run, ':'-joined).
MODDIR = os.environ.get("L2_MODDIR", os.path.join(WT, "dist", "openseesmp"))
MPIEXEC = os.environ.get(
    "L2_MPIEXEC",
    os.path.join(MODDIR, "mpiexec.exe") if os.name == "nt" else "mpirun")
PYTHON = sys.executable

sys.path.insert(0, HERE)


def driver(ne, want_modes, outdir):
    """Under mpiexec: build the box, calibrate the band, TIME the distributed
    -rci solve, write rank_<world>.json."""
    sys.path.insert(0, MODDIR)
    if hasattr(os, "add_dll_directory"):      # Windows-only DLL search dir
        os.add_dll_directory(MODDIR)
    os.environ["LADRUNO_OPENSEES_QUIET"] = "1"
    import openseesmp as ops
    import solid_box as sb

    world = ops.getPID()
    npr = ops.getNP()

    # deterministic band (same on every replicated rank). L2_FMAX overrides the
    # ARPACK calibration (reuse a coarse-mesh band on a large mesh — the lowest
    # modes are ~mesh-converged, and ARPACK on a 1M-DOF solid is itself costly).
    fmax_env = os.environ.get("L2_FMAX")
    if fmax_env:
        fmax = float(fmax_env)
    else:
        fmax, _ = sb.calibrate_fmax(ops, ne, want_modes)

    # fresh model, then TIME only the distributed -rci eigen call.
    # Pass an explicit -m0 (L2_M0) to SKIP the auto-seed: without it the solver
    # runs a full SERIAL dfeast_scsrgv stochastic estimate (complex PARDISO,
    # replicated per rank, NOT distributed) before the RCI loop — a fixed serial
    # cost that would contaminate the L3 strong-scaling measurement. A fixed m0
    # isolates the distributed-factorization scaling that L2 would multiply.
    m0 = int(os.environ.get("L2_M0", "0"))
    ndof, nnode, nele = sb.build_box(ops, ne)
    t0 = time.perf_counter()
    if m0 > 0:
        lam = list(ops.eigen("-feast", 0.0, fmax, "-m0", m0, "-rci"))
    else:
        lam = list(ops.eigen("-feast", 0.0, fmax, "-rci"))
    t1 = time.perf_counter()

    out = dict(world=world, npr=npr, ne=ne, ndof=ndof, nnode=nnode, nele=nele,
               fmax=fmax, nmodes=len(lam), t_rci=t1 - t0,
               lam=lam if world == 0 else None)
    with open(os.path.join(outdir, f"rank_{world}.json"), "w") as fh:
        json.dump(out, fh)


def _run(nranks, ne, want_modes, mkl_threads):
    import tempfile
    env = dict(os.environ)
    env["MKL_NUM_THREADS"] = str(mkl_threads)
    env["OMP_NUM_THREADS"] = str(mkl_threads)
    env["LADRUNO_FEAST_MPI_DEBUG"] = "1"
    if os.name == "nt":
        env["PATH"] = MODDIR + os.pathsep + env.get("PATH", "")
    else:
        # Linux: openseesmp.so needs its runtime deps on LD_LIBRARY_PATH
        # (OpenMPI, netlib scalapack, MKL). Pass extras via L2_LDPATH.
        ld = [MODDIR] + os.environ.get("L2_LDPATH", "").split(os.pathsep)
        ld.append(env.get("LD_LIBRARY_PATH", ""))
        env["LD_LIBRARY_PATH"] = os.pathsep.join(p for p in ld if p)
        # single-dynamic-lib MKL: pin LP64 (MKL_INT==int, matches the code) +
        # sequential threading (avoid OpenMP oversubscription under MPI).
        env.setdefault("MKL_INTERFACE_LAYER", "LP64")
        env.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
    outdir = tempfile.mkdtemp(prefix=f"l2_ne{ne}_n{nranks}_",
                              dir=os.environ.get("L2_TMPDIR") or None)
    mpi_flags = os.environ.get("L2_MPI_FLAGS", "").split()
    cmd = [MPIEXEC, "-n", str(nranks), *mpi_flags, PYTHON, "-S",
           os.path.abspath(__file__), "--driver", str(ne), str(want_modes), outdir]
    t0 = time.perf_counter()
    # esmeralda wall time is effectively unlimited (internally-controlled cluster,
    # partition MaxTime=INFINITE) -> default the per-point cap high so large models
    # don't self-kill; override with L2_TIMEOUT (seconds).
    _to = int(os.environ.get("L2_TIMEOUT", "14400"))
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, timeout=_to)
    wall = time.perf_counter() - t0
    if r.returncode != 0:
        print(r.stdout); print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"mpi run n={nranks} ne={ne} failed rc={r.returncode}")

    # count distributed solves (factorizations) hosted per rank + prove distribution
    solves = {}
    for m in re.finditer(r"LADRUNO_FEAST_MPI rank=(\d+) np=(\d+) nBlk=(\d+) "
                         r"nzLoc=(\d+)", r.stderr):
        rk, npr_dbg, nblk, nz = (int(m.group(1)), int(m.group(2)),
                                 int(m.group(3)), int(m.group(4)))
        solves.setdefault(rk, []).append(nz)
    distributed = (nranks == 1) or (set(solves.keys()) == set(range(nranks)))
    nsolves = len(solves.get(0, [])) if solves else 0

    with open(os.path.join(outdir, "rank_0.json")) as fh:
        r0 = json.load(fh)
    return dict(np=nranks, t_rci=r0["t_rci"], wall=wall, ndof=r0["ndof"],
                nmodes=r0["nmodes"], fmax=r0["fmax"], nnode=r0["nnode"],
                nele=r0["nele"], nsolves=nsolves, distributed=distributed)


def main():
    if "--driver" in sys.argv:
        i = sys.argv.index("--driver")
        driver(int(sys.argv[i + 1]), int(sys.argv[i + 2]), sys.argv[i + 3])
        return

    ne = int(sys.argv[1]) if len(sys.argv) > 1 else 30
    nps = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2
                            else "1,2,4,8".split(","))]
    want_modes = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    mkl_threads = int(os.environ.get("L2_MKL_THREADS", "1"))

    print(f"ADR-43 L2 profiler — solid box ne={ne}^3, band holds ~{want_modes} "
          f"modes, MKL_THREADS={mkl_threads}")
    print(f"  np list: {nps}\n")
    rows = []
    for npr in nps:
        print(f"  running np={npr} ...", flush=True)
        rows.append(_run(npr, ne, want_modes, mkl_threads))

    base = rows[0]["t_rci"]
    hdr = base and rows[0]
    print("\n=== RESULT ===")
    print(f"  model: ne={ne}^3  ndof={rows[0]['ndof']:,}  nele={rows[0]['nele']:,}"
          f"  nmodes={rows[0]['nmodes']}  fmax={rows[0]['fmax']:.4g} Hz"
          f"  nsolves(rank0)={rows[0]['nsolves']}")
    print(f"  {'np':>3} {'t_rci(s)':>10} {'speedup':>8} {'effic':>7} "
          f"{'t/solve(s)':>11} {'dist':>5}")
    for row in rows:
        sp = base / row["t_rci"] if row["t_rci"] else float("nan")
        eff = sp / row["np"]
        tps = row["t_rci"] / row["nsolves"] if row["nsolves"] else float("nan")
        print(f"  {row['np']:>3} {row['t_rci']:>10.3f} {sp:>8.2f} {eff:>7.2f} "
              f"{tps:>11.4f} {str(row['distributed']):>5}")

    outp = os.path.join(HERE, f"scaling_ne{ne}.json")
    with open(outp, "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"\n  wrote {outp}")

    # --- L2 ceiling from the L3 curve (derivation in README) ---
    # In the factorization-dominated limit, running N_q_eff contour solves
    # CONCURRENTLY on P/G-rank groups (L2) vs SEQUENTIALLY on P ranks (L3) gives
    #   speedup_L2/L3  ~=  E(P/G) / E(P)          (G = N_q_eff groups)
    # E = strong-scaling efficiency. Flat E (perfect L3) => 1 => L2 useless;
    # E collapsing over the last factor-of-G of ranks => that's the L2 recovery.
    eff = {}
    for row in rows:
        sp = base / row["t_rci"] if row["t_rci"] else 0.0
        eff[row["np"]] = sp / row["np"]
    P = rows[-1]["np"]
    print(f"\n  L3 strong-scaling efficiency at np={P}: {eff[P]:.0%}")
    for G in (4, 8):                      # N_q_eff after half-contour (4), full (8)
        Pg = P // G
        if Pg in eff and eff[P] > 0:
            ceil_ = eff[Pg] / eff[P]
            print(f"  L2 ceiling estimate (G={G} groups, P/G={Pg} ranks/group): "
                  f"E({Pg})/E({P}) = {eff[Pg]:.2f}/{eff[P]:.2f} = {ceil_:.2f}x")
    print("\n  READING: L2-over-L3 speedup ceiling is E(P/G)/E(P). >~2x on a real "
          "rank budget => L2 worth the (largest novel-math) build; ~1x => not. "
          "This LOCAL np<=8 run is a first estimate; the real budget is esmeralda "
          "np=16-64 (push P until E(P) collapses, then E(P/G)/E(P) is the L2 win).")


if __name__ == "__main__":
    main()
