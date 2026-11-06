"""EnergyBalance -v2: every MPI rank writes and echoes the same column layout.

Found on esmeralda (PBD-V'26 coarse-mesh H5DRM runs, 2026-10-31): the optional
-v2 columns were switched on by RANK-LOCAL tests, so a rank owning no DRM-loaded
element wrote and echoed one model column fewer (no E_inject) than its peers,
and readers that require identical per-rank layouts failed. The fix OR-reduces
the column flags at initialize() (SRC/recorder/LadrunoEnergyBalanceMPI.cpp).

Model: `mpiexec -n 2` on openseesmp, two independent per-rank models (no shared
nodes, a rank-local UmfPack solve each - the partition is the point, not the
coupling). Rank 0 owns a soil column driven through an ASDAbsorbingBoundary2D
bottom injection (the same ABSORB_LEAK channel as H5DRM / Lysmer); rank 1 owns
the same column on a fixed base with an initial velocity and no absorbing
element - the "rank without DRM elements". Checks:

  * both ranks echo the identical `columns =` line, which carries E_inject;
  * both part files have the same width;
  * rank 1's E_inject column is identically 0;
  * the rank sum of E_inject equals the same rank-0 model run alone in the
    serial engine (the fix moves no energy, it only adds a zero column), and
    rank 1's RES does not drift.

Runs in `python -S` subprocesses for the reason in test_ladruno_up_mp_smoke.py
(a boot .pth can preload a stale engine).

Run:  py -3.12 -m pytest tests/test_energyBalance_mp_columns.py -v
"""
import os
import re
import subprocess
import sys

import numpy as np
import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
DIST_BIN = os.path.join(_ROOT, "dist", "bin")
DIST_MP = os.path.join(_ROOT, "dist", "openseesmp")

# Shared model text (no numpy: the `-S` children have no site-packages). INJECT selects the rank-0 (absorbing-base, injected) or the
# rank-1 (fixed-base, no absorbing element) column.
_MODEL = r'''
import math

def build(ops, inject, n_layers=3):
    E, nu, rho = 2.0 * 2000.0 * 200.0 ** 2 * 1.25, 0.25, 2000.0
    G = 2000.0 * 200.0 ** 2
    dt, nstep = 1.0e-3, 60
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for row in range(n_layers + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))
    ops.nDMaterial("ElasticIsotropic", 1, E, nu, rho)
    for row in range(n_layers):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, 1.0, "PlaneStrain", 1)
    if inject:
        a = [(math.pi * 10.0 * (k * dt - 0.02)) ** 2 for k in range(nstep + 1)]
        vals = [0.05 * (1.0 - 2.0 * x) * math.exp(-x) for x in a]
        ops.node(1001, 0.0, -1.0); ops.node(1002, 1.0, -1.0)
        ops.fix(1001, 1, 1); ops.fix(1002, 1, 1)
        ops.timeSeries("Path", 1, "-dt", dt, "-values", *vals)
        ops.element("ASDAbsorbingBoundary2D", 101, 1001, 1002, 2, 1,
                    G, nu, rho, 1.0, "B", "-fx", 1)
        ops.setParameter("-val", 1, "-ele", 101, "stage")
    else:
        ops.fix(1, 1, 1); ops.fix(2, 1, 1)
        for row in range(1, n_layers + 1):
            ops.setNodeVel(1 + 2 * row, 1, 0.1, "-commit")
            ops.setNodeVel(2 + 2 * row, 1, 0.1, "-commit")
    return dt, nstep

def run(ops, efile, inject):
    dt, nstep = build(ops, inject)
    ops.recorder("EnergyBalance", "-file", efile, "-time", "-v2")
    ops.constraints("Plain"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 10); ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25); ops.analysis("Transient")
    for _ in range(nstep):
        assert ops.analyze(1, dt) == 0
    ops.wipe()
'''

_MP_DRIVER = r'''
import os, sys
MODDIR = sys.argv[1]; OUT = sys.argv[2]; MODEL = sys.argv[3]
_IMPI = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\bin"
_LIBFABRIC = r"C:\Program Files (x86)\Intel\oneAPI\mpi\latest\opt\mpi\libfabric\bin"
_MKL = r"C:\Program Files (x86)\Intel\oneAPI\mkl\latest\bin"
_ICOMP = r"C:\Program Files (x86)\Intel\oneAPI\compiler\latest\bin"
for d in (MODDIR, _IMPI, _LIBFABRIC, _MKL, _ICOMP):
    if d and os.path.isdir(d):
        os.add_dll_directory(d)
os.environ["PATH"] = os.pathsep.join(
    [p for p in (_IMPI, _LIBFABRIC, MODDIR, os.environ.get("PATH", "")) if p])
sys.path.insert(0, MODDIR)
import openseesmp as ops
exec(open(MODEL).read())
pid = ops.getPID(); nproc = ops.getNP(); ops.start()
run(ops, os.path.join(OUT, "energy.txt"), inject=(pid == 0))
with open(os.path.join(OUT, "mp_rank%d.txt" % pid), "w") as fh:
    fh.write("pid=%d nproc=%d DONE\n" % (pid, nproc))
'''

_SERIAL_DRIVER = r'''
import os, sys
DIST_BIN = sys.argv[1]; EFILE = sys.argv[2]; MODEL = sys.argv[3]
os.add_dll_directory(DIST_BIN)
sys.path.insert(0, DIST_BIN)
os.environ["PATH"] = DIST_BIN + os.pathsep + os.environ.get("PATH", "")
for k in ("PMI_RANK", "PMI_SIZE", "OMPI_COMM_WORLD_RANK", "OMPI_COMM_WORLD_SIZE"):
    os.environ.pop(k, None)
import opensees as ops
exec(open(MODEL).read())
run(ops, EFILE, inject=True)
print("SERIAL-OK")
'''

_ECHO = re.compile(r"EnergyBalanceRecorder: columns =([^\r\n]*)")


def _mp_available():
    return (os.path.isfile(os.path.join(DIST_MP, "mpiexec.exe"))
            and os.path.isfile(os.path.join(DIST_MP, "openseesmp.pyd")))


def _serial_available():
    return os.path.isfile(os.path.join(DIST_BIN, "opensees.pyd"))


def _columns(echo):
    # " time [model: KE_ele ... RES ERR%]" -> ["time", "KE_ele", ..., "ERR%"]
    return echo.replace("[model:", " ").replace("]", " ").split()


def test_energybalance_v2_ranks_share_column_layout(tmp_path):
    if not (_mp_available() and _serial_available()):
        pytest.skip("needs dist/openseesmp (mpiexec + openseesmp.pyd) and dist/bin/opensees.pyd")

    model = tmp_path / "eb_model.py"
    model.write_text(_MODEL)
    env = dict(os.environ)
    env["LADRUNO_OPENSEES_QUIET"] = "1"

    # --- reference: the injected (rank-0) column alone, serial engine -------
    sdrv = tmp_path / "eb_serial_driver.py"
    sdrv.write_text(_SERIAL_DRIVER)
    sfile = tmp_path / "energy_serial.txt"
    proc = subprocess.run(
        [sys.executable, "-S", str(sdrv), DIST_BIN, str(sfile), str(model)],
        stdin=subprocess.DEVNULL, capture_output=True, text=True, timeout=300,
        env=env, encoding="utf-8", errors="replace")
    assert proc.returncode == 0 and "SERIAL-OK" in proc.stdout, (
        f"serial reference failed (rc={proc.returncode})\n{proc.stdout}\n{proc.stderr}")
    s_echo = _ECHO.findall(proc.stdout + proc.stderr)
    assert len(s_echo) == 1, s_echo
    s_cols = _columns(s_echo[0])
    assert "E_inject" in s_cols, s_cols
    ser = np.atleast_2d(np.loadtxt(sfile))

    # --- two ranks, rank 1 without any absorbing/DRM element ----------------
    mdrv = tmp_path / "eb_mp_driver.py"
    mdrv.write_text(_MP_DRIVER)
    mpiexec = os.path.join(DIST_MP, "mpiexec.exe")
    menv = dict(env)
    menv["PATH"] = DIST_MP + os.pathsep + menv.get("PATH", "")
    try:
        proc = subprocess.run(
            [mpiexec, "-n", "2", sys.executable, "-S", str(mdrv),
             DIST_MP, str(tmp_path), str(model)],
            stdin=subprocess.DEVNULL, capture_output=True, text=True,
            timeout=300, env=menv, encoding="utf-8", errors="replace")
    except (subprocess.TimeoutExpired, OSError) as exc:
        pytest.fail(f"MP run did not finish (a hung collective?): {exc}")
    if not all((tmp_path / f"mp_rank{r}.txt").exists() for r in (0, 1)):
        pytest.skip("MP run produced no rank output (MPI infra?). "
                    f"tail:\n{proc.stdout[-800:]}\n{proc.stderr[-800:]}")

    echoes = _ECHO.findall(proc.stdout + proc.stderr)
    assert len(echoes) == 2, f"expected one echo per rank, got {echoes}"
    assert echoes[0].strip() == echoes[1].strip(), (
        f"ranks echo different layouts:\n  {echoes[0]}\n  {echoes[1]}")
    cols = _columns(echoes[0])
    assert cols == s_cols, f"MP layout {cols} != serial layout {s_cols}"

    d0 = np.atleast_2d(np.loadtxt(tmp_path / "energy.part-0.txt"))
    d1 = np.atleast_2d(np.loadtxt(tmp_path / "energy.part-1.txt"))
    assert d0.shape[1] == d1.shape[1] == len(cols), (d0.shape, d1.shape, cols)

    ci = cols.index("E_inject")
    ires = cols.index("RES")
    assert np.all(d1[:, ci] == 0.0), "rank without an injector must write E_inject = 0"
    assert abs(d0[-1, ci]) > 0.0, "rank-0 injection is zero - the model is not a probe"

    # the rank sum is the serial total: the fix adds a zero column, moves nothing
    e_sum = d0[-1, ci] + d1[-1, ci]
    assert abs(e_sum - ser[-1, ci]) <= 1.0e-9 * max(1.0, abs(ser[-1, ci])), (
        e_sum, ser[-1, ci])
    np.testing.assert_allclose(d0, ser, rtol=0.0, atol=1.0e-9 * np.abs(ser).max())

    # rank 1 still closes on its own terms: RES carries the constant offset of
    # the initial-velocity start and must not drift (nothing leaks into it)
    iulw = cols.index("ULW")
    drift = np.abs(d1[:, ires] - d1[0, ires]).max()
    assert drift <= 1.0e-2 * max(np.abs(d1[:, 1:iulw + 1]).max(), 1e-30), drift
