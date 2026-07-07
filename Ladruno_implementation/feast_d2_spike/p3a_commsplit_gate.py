"""ADR-43 P3a gate — MumpsParallelSolver honors an MPI_Comm_split sub-communicator.

This is the D2 spike's originally-scoped toy, promoted to a permanent unit
test: 4 openseesmp ranks split into 2 sub-communicator groups via
`system('Mumps', '-commSplit', color)`. Each group solves a DIFFERENT static
model concurrently through its own MumpsParallelSOE:

  group 0 (world ranks 0,1): 2-rank flat-per-rank truss, load P0
  group 1 (world ranks 2,3): same topology, load 3*P0 and stiffer material

Gates:
  G1  each group's solution matches a SERIAL reference of its own model
      (rel <= 1e-10) — the sub-comm factorization is correct.
  G2  the two groups' answers differ by the expected factor — no cross-talk:
      a WORLD-comm regression would deadlock or corrupt one group's solve.
  G3  control: the same 4-rank run with the plain WORLD path (no -commSplit,
      single 4-rank model) still solves — the default is untouched.

Run:  python p3a_commsplit_gate.py            (orchestrator: serial + mpi runs)
      python p3a_commsplit_gate.py --driver <mode>   (internal, under mpiexec)
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))
DIST_MP = os.path.join(WT, "dist", "openseesmp")
DIST_BIN = os.path.join(WT, "dist", "bin")
PYTHON = sys.executable

# NOTE: P/E must differ between the groups (displacement scales as P/E for
# this linear model) or G2's "answers differ" check is vacuous.
E_SOFT, E_STIFF = 3000.0, 9000.0
P_SOFT, P_STIFF = 10.0, 45.0
AREA = 10.0


def _build_group_model(ops, local_rank, E, P):
    """Flat-per-rank 2-rank truss (or serial when local_rank < 0).

    Shared loaded node 4; rank 0 owns element 1, rank 1 owns elements 2+3.
    """
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, E)

    if local_rank <= 0:      # serial or group-local rank 0
        ops.node(1, 0.0, 0.0)
        ops.node(4, 72.0, 96.0)
        ops.fix(1, 1, 1)
        ops.element("Truss", 1, 1, 4, AREA, 1)
    if local_rank != 0:      # serial or group-local rank 1
        ops.node(2, 144.0, 0.0)
        ops.node(3, 168.0, 0.0)
        if local_rank > 0:
            ops.node(4, 72.0, 96.0)
        ops.fix(2, 1, 1)
        ops.fix(3, 1, 1)
        ops.element("Truss", 2, 2, 4, AREA / 2.0, 1)
        ops.element("Truss", 3, 3, 4, AREA / 2.0, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if local_rank <= 0:      # load applied once (on the group host / serial)
        ops.load(4, P, -P)


def _analyze_static(ops, comm_split_color=None):
    """MP static solve. ParallelPlain numbering is WORLD-collective (it runs
    on the interpreter's channels, independent of the Mumps communicator);
    consistent global numbering across the two groups holds because both
    groups build byte-identical topologies (same node/element tags and
    fixities — only E and P differ, which numbering never sees)."""
    ops.numberer("ParallelPlain")
    ops.constraints("Plain")
    if comm_split_color is not None:
        ops.system("Mumps", "-commSplit", comm_split_color)
    else:
        ops.system("Mumps")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    return ops.analyze(1)


def driver_split():
    """Runs under mpiexec -n 4: ranks 0,1 -> group 0; ranks 2,3 -> group 1."""
    sys.path.insert(0, DIST_MP)
    os.add_dll_directory(DIST_MP)
    import openseesmp as ops

    world = ops.getPID()
    color = 0 if world < 2 else 1
    local_rank = world % 2
    E = E_SOFT if color == 0 else E_STIFF
    P = P_SOFT if color == 0 else P_STIFF

    _build_group_model(ops, local_rank, E, P)
    ok = _analyze_static(ops, comm_split_color=color)
    ux, uy = ops.nodeDisp(4, 1), ops.nodeDisp(4, 2)
    print(f"RESULT {json.dumps(dict(world=world, color=color, ok=ok, ux=ux, uy=uy))}",
          flush=True)


def driver_world():
    """Control: plain 4-rank WORLD-comm Mumps solve (no -commSplit).

    Flat-per-rank: rank 0 hosts element 1 + load; ranks 1..3 each add one
    supported leg to the shared node.
    """
    sys.path.insert(0, DIST_MP)
    os.add_dll_directory(DIST_MP)
    import openseesmp as ops

    world = ops.getPID()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.uniaxialMaterial("Elastic", 1, E_SOFT)
    ops.node(4, 72.0, 96.0)
    if world == 0:
        ops.node(1, 0.0, 0.0)
        ops.fix(1, 1, 1)
        ops.element("Truss", 1, 1, 4, AREA, 1)
    else:
        ops.node(10 + world, 144.0 + world, 0.0)
        ops.fix(10 + world, 1, 1)
        ops.element("Truss", 10 + world, 10 + world, 4, AREA / 3.0, 1)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if world == 0:
        ops.load(4, P_SOFT, -P_SOFT)

    ok = _analyze_static(ops, comm_split_color=None)
    ux, uy = ops.nodeDisp(4, 1), ops.nodeDisp(4, 2)
    print(f"RESULT {json.dumps(dict(world=world, ok=ok, ux=ux, uy=uy))}",
          flush=True)


def _mpi_run(mode, nranks):
    mpiexec = os.path.join(DIST_MP, "mpiexec.exe")
    env = dict(os.environ)
    env["PATH"] = DIST_MP + os.pathsep + env.get("PATH", "")
    cmd = [mpiexec, "-n", str(nranks), PYTHON, "-S",
           os.path.abspath(__file__), "--driver", mode]
    r = subprocess.run(cmd, env=env, capture_output=True, text=True,
                       timeout=300)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr, file=sys.stderr)
        raise RuntimeError(f"mpi run '{mode}' failed rc={r.returncode}")
    results = [json.loads(line.split("RESULT ", 1)[1])
               for line in r.stdout.splitlines() if line.startswith("RESULT ")]
    return results


def main():
    if "--driver" in sys.argv:
        mode = sys.argv[sys.argv.index("--driver") + 1]
        # -S runs need manual site wiring for os.add_dll_directory users
        {"split": driver_split, "world": driver_world}[mode]()
        return

    # serial oracles (run in subprocesses so the two module loads don't clash)
    def _oracle(E, P):
        code = (f"import sys; sys.path.insert(0, {DIST_BIN!r}); "
                f"import os; os.add_dll_directory({DIST_BIN!r}); "
                f"import opensees as ops; "
                f"ops.model('basic','-ndm',2,'-ndf',2); "
                f"ops.uniaxialMaterial('Elastic',1,{E}); "
                f"ops.node(1,0.,0.); ops.node(4,72.,96.); "
                f"ops.node(2,144.,0.); ops.node(3,168.,0.); "
                f"ops.fix(1,1,1); ops.fix(2,1,1); ops.fix(3,1,1); "
                f"ops.element('Truss',1,1,4,{AREA},1); "
                f"ops.element('Truss',2,2,4,{AREA/2.},1); "
                f"ops.element('Truss',3,3,4,{AREA/2.},1); "
                f"ops.timeSeries('Linear',1); ops.pattern('Plain',1,1); "
                f"ops.load(4,{P},-{P}); "
                f"ops.numberer('Plain'); ops.constraints('Plain'); "
                f"ops.system('FullGeneral'); "
                f"ops.integrator('LoadControl',1.); ops.algorithm('Linear'); "
                f"ops.analysis('Static'); assert ops.analyze(1)==0; "
                f"import json; print('RESULT', json.dumps("
                f"[ops.nodeDisp(4,1), ops.nodeDisp(4,2)]))")
        # -S: skip site/.pth processing (the boot .pth preloads a stale pyd
        # from another worktree — see project_opensees_test_env quirk)
        r = subprocess.run([PYTHON, "-S", "-c", code], capture_output=True,
                           text=True, timeout=120)
        assert r.returncode == 0, r.stderr
        line = [l for l in r.stdout.splitlines() if l.startswith("RESULT")][0]
        return json.loads(line.split("RESULT ", 1)[1])

    ref0 = _oracle(E_SOFT, P_SOFT)
    ref1 = _oracle(E_STIFF, P_STIFF)
    print(f"serial refs: group0={ref0}  group1={ref1}")

    # G1 + G2: the split run
    res = _mpi_run("split", 4)
    assert len(res) == 4, f"expected 4 rank results, got {len(res)}"
    for r in res:
        assert r["ok"] == 0, f"rank {r['world']} analyze failed"
        ref = ref0 if r["color"] == 0 else ref1
        for got, want, dof in ((r["ux"], ref[0], "ux"), (r["uy"], ref[1], "uy")):
            rel = abs(got - want) / max(abs(want), 1e-30)
            assert rel <= 1e-10, (f"G1 FAIL rank {r['world']} color {r['color']} "
                                  f"{dof}: got {got}, want {want}, rel {rel:.2e}")
    g0 = [r for r in res if r["color"] == 0][0]
    g1 = [r for r in res if r["color"] == 1][0]
    assert abs(g0["ux"] - g1["ux"]) > 1e-12, "G2 FAIL: groups returned identical ux"
    print(f"G1+G2 PASS: group0 ux={g0['ux']:.6e}, group1 ux={g1['ux']:.6e} "
          f"(both match their own serial oracle, and differ from each other)")

    # G3: control — WORLD path untouched
    res_w = _mpi_run("world", 4)
    assert len(res_w) == 4 and all(r["ok"] == 0 for r in res_w), "G3 FAIL"
    ux_set = {round(r["ux"], 14) for r in res_w}
    assert len(ux_set) == 1, f"G3 FAIL: WORLD ranks disagree on ux: {ux_set}"
    print(f"G3 PASS: plain WORLD-comm Mumps run solves, all ranks agree "
          f"(ux={res_w[0]['ux']:.6e})")
    print("ALL GATES PASS")


if __name__ == "__main__":
    main()
