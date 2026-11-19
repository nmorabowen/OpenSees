"""WP-240 -- LadrunoDynamicRelaxation in OpenSeesMP: one global kinetic energy
and one global residual.

Ported from jaabell/OpenSees 3c2f585b5 (José Abell). Each process used to
compute the kinetic energy, the residual and max|u| of its own equations, so
the processes zeroed their velocities at different kinetic-energy peaks and an
interface DOF saw only its own elements' share. The reductions are now global
(`Ladruno::drAllreduceSum/Max`, pointed at `MPI_Allreduce` by the two parallel
interpreters), a run on several processes must use `system MPIDiagonal`, and a
distributed diagonal system keeps its first M* (no refresh).

Model: a bar of 8 LadrunoBrick elements of LadrunoJ2 along x, split by x over
the processes, loaded past yield (tests/tcl/wp240_dr_mp.tcl, and the same
model in Python below).

  P-1  OpenSeesMP.exe on 1, 2 and 4 processes: residual, kinetic energy and
       max|u| of every step equal the serial OpenSees.exe run with
       -noAutoRefresh, identical on every process; same nodal displacements.
  P-2  the same for openseesmp.pyd against opensees.pyd.
  P-3  `ladrunoDR relax`: every process returns the same status and step
       count, the serial one.
  P-4  several processes with `system Diagonal` fail at the first step and say
       why; one process may use it.
  P-5  without -noAutoRefresh the run says the refresh is off and gives the
       same numbers.

Runs the engines as subprocesses (`python -S` for the modules: a boot .pth can
preload a stale engine, see test_ladruno_up_mp_smoke.py).
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "tcl" / "wp240_dr_mp.tcl"
DIST = REPO / "dist"
_EXE = ".exe" if os.name == "nt" else ""
SER_EXE = DIST / "bin" / ("OpenSees" + _EXE)
MP_EXE = DIST / "bin" / ("OpenSeesMP" + _EXE)
MPIEXEC = DIST / "openseesmp" / ("mpiexec" + _EXE)
SER_PYD = DIST / "bin" / "opensees.pyd"
MP_PYD = DIST / "openseesmp" / "openseesmp.pyd"

needs_exe = pytest.mark.skipif(
    not (SER_EXE.exists() and MP_EXE.exists() and MPIEXEC.exists()),
    reason="dist OpenSees / OpenSeesMP / mpiexec not built")
needs_pyd = pytest.mark.skipif(
    not (SER_PYD.exists() and MP_PYD.exists() and MPIEXEC.exists()),
    reason="dist opensees.pyd / openseesmp.pyd / mpiexec not built")

_RTOL = 1.0e-9
_NOTE = "MPIDiagonal keeps the first M*"
_GUARD = "use `system MPIDiagonal`"

# The Tcl deck, in Python. No numpy and no pytest: it runs under `python -S`.
_PY_DRIVER = r'''
import os, sys
ENGINE, MODDIR, OUT, MODE, SYSTEM, REFRESH = sys.argv[1:7]
if ENGINE == "mp":
    _ONEAPI = r"C:\Program Files (x86)\Intel\oneAPI"
    _DIRS = [MODDIR, _ONEAPI + r"\mpi\latest\bin",
             _ONEAPI + r"\mpi\latest\opt\mpi\libfabric\bin",
             _ONEAPI + r"\mkl\latest\bin", _ONEAPI + r"\compiler\latest\bin"]
    for d in _DIRS:
        if os.path.isdir(d) and hasattr(os, "add_dll_directory"):
            os.add_dll_directory(d)
    os.environ["PATH"] = os.pathsep.join(
        [d for d in _DIRS if os.path.isdir(d)] + [os.environ.get("PATH", "")])
    sys.path.insert(0, MODDIR)
    import openseesmp as ops
    pid, np_ = ops.getPID(), ops.getNP()
    ops.start()
else:
    if hasattr(os, "add_dll_directory"):
        os.add_dll_directory(MODDIR)
    os.environ["PATH"] = MODDIR + os.pathsep + os.environ.get("PATH", "")
    for k in ("PMI_RANK", "PMI_SIZE", "OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE"):
        os.environ.pop(k, None)
    sys.path.insert(0, MODDIR)
    import opensees as ops
    pid, np_ = 0, 1
if SYSTEM == "auto":
    SYSTEM = "MPIDiagonal" if ENGINE == "mp" else "Diagonal"

NEL = 8
E, nu, sig0 = 2.0e5, 0.3, 250.0
K, G = E / (3.0 * (1.0 - 2.0 * nu)), E / (2.0 * (1.0 + nu))
P = 0.5 * sig0 if MODE == "relax" else 1.1 * sig0

def ntag(i, j, k):
    return 1 + 4 * i + 2 * j + k

ops.model("basic", "-ndm", 3, "-ndf", 3)
ops.nDMaterial("LadrunoJ2", 1, K, G, "-iso", "voce", sig0, 0.0, 0.0, 200.0)
mine = [e for e in range(NEL) if (e * np_) // NEL == pid]
held = []
for e in mine:
    for i in (e, e + 1):
        if i in held:
            continue
        held.append(i)
        for j in (0, 1):
            for k in (0, 1):
                ops.node(ntag(i, j, k), float(i), float(j), float(k))
for e in mine:
    a, b = e, e + 1
    ops.element("LadrunoBrick", e + 1,
                ntag(a, 0, 0), ntag(b, 0, 0), ntag(b, 1, 0), ntag(a, 1, 0),
                ntag(a, 0, 1), ntag(b, 0, 1), ntag(b, 1, 1), ntag(a, 1, 1), 1)
if 0 in held:
    for j in (0, 1):
        for k in (0, 1):
            ops.fix(ntag(0, j, k), 1, 1, 1)
ops.timeSeries("Constant", 1)
ops.pattern("Plain", 1, 1)
if 7 in mine:
    for j in (0, 1):
        for k in (0, 1):
            ops.load(ntag(8, j, k), P / 4.0, 0.0, 0.0)
if 3 in mine:
    for j in (0, 1):
        for k in (0, 1):
            ops.load(ntag(4, j, k), 0.0, 0.0, -0.01 * P / 4.0)

ops.constraints("Transformation")
ops.numberer("ParallelPlain" if ENGINE == "mp" else "Plain")
ops.system(SYSTEM)
ops.test("NormUnbalance", 1.0e30, 1, 0)
ops.algorithm("Linear")
if REFRESH == "1":
    ops.integrator("LadrunoDynamicRelaxation")
else:
    ops.integrator("LadrunoDynamicRelaxation", "-noAutoRefresh")
ops.analysis("Transient")

rc, result = 0, []
if MODE == "relax":
    try:
        result = ops.ladrunoDR("relax", "-tolResidual", 1.0e-6 * P, "-maxSteps", 20000)
    except Exception:
        rc = -1
else:
    with open("%s.hist%d" % (OUT, pid), "w") as fh:
        for s in range(1, 401):
            rc = ops.analyze(1, 1.0)
            if rc != 0:
                break
            fh.write("%d %.17g %.17g %.17g\n" % (
                s, ops.ladrunoDR("residualNorm"), ops.ladrunoDR("kineticEnergy"),
                ops.ladrunoDR("maxDisplacement")))
with open("%s.node%d" % (OUT, pid), "w") as fh:
    for i in sorted(held):
        for j in (0, 1):
            for k in (0, 1):
                t = ntag(i, j, k)
                fh.write("%d %.17g %.17g %.17g\n" % ((t,) + tuple(ops.nodeDisp(t))))
with open("%s.rc%d" % (OUT, pid), "w") as fh:
    fh.write("rc %d\n" % rc)
    fh.write("result " + " ".join("%.17g" % v for v in result) + "\n")
'''


def _launch(cmd, env, cwd):
    try:
        p = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                           encoding="utf-8", errors="replace", timeout=600, env=env,
                           cwd=str(cwd))
    except subprocess.TimeoutExpired as exc:
        pytest.fail("run did not finish (a hung collective?): %s" % exc)
    return p.stdout + p.stderr


def _mpi_prefix(np, env):
    env["PATH"] = str(MPIEXEC.parent) + os.pathsep + env.get("PATH", "")
    return [str(MPIEXEC)] + (["-localonly"] if os.name == "nt" else []) + ["-n", str(np)]


def _run(tmp, engine, np, mode="hist", system=None, refresh=False):
    """engine: 'tcl' (OpenSees.exe / OpenSeesMP.exe) or 'py' (opensees.pyd /
    openseesmp.pyd). np = 0 is the serial engine; np >= 1 the MP one."""
    tmp = Path(tmp)
    name = "%s_%s_np%d_%s_%d" % (engine, mode, np, system or "auto", int(refresh))
    out = tmp / name
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    if engine == "tcl":
        env["WP240_OUT"] = str(out)
        env["WP240_MODE"] = mode
        env["WP240_REFRESH"] = "1" if refresh else "0"
        if system:
            env["WP240_SYSTEM"] = system
        tcl = DIST / "lib" / "tcl8.6"
        if "TCL_LIBRARY" not in env and tcl.exists():
            env["TCL_LIBRARY"] = str(tcl)
        cmd = [str(SER_EXE), str(DECK)] if np == 0 else \
            _mpi_prefix(np, env) + [str(MP_EXE), str(DECK)]
    else:
        drv = tmp / "wp240_driver.py"
        if not drv.exists():
            drv.write_text(_PY_DRIVER)
        tail = [str(out), mode, system or "auto", "1" if refresh else "0"]
        if np == 0:
            cmd = [sys.executable, "-S", str(drv), "serial", str(SER_PYD.parent)] + tail
        else:
            cmd = _mpi_prefix(np, env) + [sys.executable, "-S", str(drv), "mp",
                                          str(MP_PYD.parent)] + tail
    log = _launch(cmd, env, tmp)
    nproc = max(np, 1)
    ranks = []
    for r in range(nproc):
        f = Path("%s.rc%d" % (out, r))
        assert f.exists(), "%s: process %d wrote no result\n%s" % (name, r, log[-3000:])
        lines = f.read_text().splitlines()
        d = {"rc": int(lines[0].split()[1]),
             "result": [float(x) for x in lines[1].split()[1:]],
             "hist": [], "node": {}}
        h = Path("%s.hist%d" % (out, r))
        if h.exists():
            d["hist"] = [[float(x) for x in ln.split()] for ln in h.read_text().splitlines()]
        for ln in Path("%s.node%d" % (out, r)).read_text().splitlines():
            t = ln.split()
            d["node"][int(t[0])] = [float(x) for x in t[1:]]
        ranks.append(d)
    return {"log": log, "ranks": ranks, "name": name}


def _close(a, b, scale):
    return abs(a - b) <= _RTOL * scale


def _assert_same_history(run, ref):
    href = ref["ranks"][0]["hist"]
    assert len(href) == 400, (ref["name"], len(href), ref["log"][-2000:])
    # the model is a probe: the residual moved and the bar went past yield
    assert href[0][1] > 10.0 * href[-1][1] > 0.0
    scale = [max(abs(row[c]) for row in href) for c in range(4)]
    for r, d in enumerate(run["ranks"]):
        assert d["rc"] == 0, (run["name"], r, run["log"][-2000:])
        assert len(d["hist"]) == 400, (run["name"], r, len(d["hist"]))
        # global values: every process reports the same numbers
        assert d["hist"] == run["ranks"][0]["hist"], (run["name"], r)
        for row, rrow in zip(d["hist"], href):
            for c in (1, 2, 3):
                assert _close(row[c], rrow[c], scale[c]), (
                    run["name"], r, int(row[0]), c, row[c], rrow[c])
    umax = scale[3]
    seen = set()
    for d in run["ranks"]:
        for tag, u in d["node"].items():
            seen.add(tag)
            for a, b in zip(u, ref["ranks"][0]["node"][tag]):
                assert _close(a, b, umax), (run["name"], tag, u)
    assert seen == set(ref["ranks"][0]["node"]), run["name"]


# ------------------------------------------------------------------- P-1 / P-2

@needs_exe
@pytest.mark.parametrize("np", [1, 2, 4])
def test_P1_OpenSeesMP_exe_matches_serial(tmp_path, np):
    ref = _run(tmp_path, "tcl", 0)
    run = _run(tmp_path, "tcl", np)
    _assert_same_history(run, ref)


@needs_pyd
@pytest.mark.parametrize("np", [1, 2, 4])
def test_P2_openseesmp_pyd_matches_serial(tmp_path, np):
    ref = _run(tmp_path, "py", 0)
    run = _run(tmp_path, "py", np)
    _assert_same_history(run, ref)


@needs_exe
@needs_pyd
def test_P2_the_two_serial_engines_agree(tmp_path):
    a = _run(tmp_path, "tcl", 0)
    b = _run(tmp_path, "py", 0)
    _assert_same_history(a, b)


# ------------------------------------------------------------------------- P-3

def _assert_same_relax(run, ref):
    rref = ref["ranks"][0]
    assert rref["rc"] == 0 and len(rref["result"]) == 5, (ref["name"], ref["log"][-2000:])
    assert int(rref["result"][0]) == 0 and 1 < rref["result"][1] < 20000, rref["result"]
    for r, d in enumerate(run["ranks"]):
        assert d["rc"] == 0 and len(d["result"]) == 5, (run["name"], r, run["log"][-2000:])
        assert d["result"] == run["ranks"][0]["result"], (run["name"], r)
        assert d["result"][:2] == rref["result"][:2], (run["name"], d["result"], rref["result"])
        assert _close(d["result"][4], rref["result"][4], rref["result"][4])


@needs_exe
@pytest.mark.parametrize("np", [2, 4])
def test_P3_relax_in_OpenSeesMP_exe(tmp_path, np):
    _assert_same_relax(_run(tmp_path, "tcl", np, mode="relax"),
                       _run(tmp_path, "tcl", 0, mode="relax"))


@needs_pyd
@pytest.mark.parametrize("np", [2, 4])
def test_P3_relax_in_openseesmp_pyd(tmp_path, np):
    _assert_same_relax(_run(tmp_path, "py", np, mode="relax"),
                       _run(tmp_path, "py", 0, mode="relax"))


# ------------------------------------------------------------------------- P-4

def _assert_guard(tmp_path, engine):
    run = _run(tmp_path, engine, 2, system="Diagonal")
    for r, d in enumerate(run["ranks"]):
        assert d["rc"] != 0, (run["name"], r)
        assert d["hist"] == [], (run["name"], r)
    assert _GUARD in run["log"], run["log"][-2000:]
    # one process is not a parallel run: `system Diagonal` is allowed
    one = _run(tmp_path, engine, 1, system="Diagonal")
    _assert_same_history(one, _run(tmp_path, engine, 0))
    assert _GUARD not in one["log"]


@needs_exe
def test_P4_several_processes_need_MPIDiagonal_exe(tmp_path):
    _assert_guard(tmp_path, "tcl")


@needs_pyd
def test_P4_several_processes_need_MPIDiagonal_pyd(tmp_path):
    _assert_guard(tmp_path, "py")


# ------------------------------------------------------------------------- P-5

@needs_pyd
def test_P5_a_distributed_system_keeps_the_first_mass(tmp_path):
    ref = _run(tmp_path, "py", 0)
    run = _run(tmp_path, "py", 2, refresh=True)
    assert _NOTE in run["log"], run["log"][-2000:]
    _assert_same_history(run, ref)
    quiet = _run(tmp_path, "py", 2)
    assert _NOTE not in quiet["log"]
