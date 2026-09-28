"""WP-132 (TIMs F22) — `system Pardiso -deterministic` / `-cbwr <BRANCH>`.

Threaded MKL PARDISO is not byte-reproducible run to run (75c Trap 7: a 14^3
Lane-B model at 4 threads split 5/5 over 10 runs, ~1 ULP; size-dependent), and
on a path-dependent nonlinear deck that last bit is amplified — TIMs §1.6 saw a
bearing wall move 30 % between two identical runs. The mode turns on MKL's
Conditional Numerical Reproducibility (mkl_cbwr_set + PARDISO iparm(34)).

Every case runs in a CHILD process: the CNR mode is process-wide and sticky, so
setting it inside the pytest process would leak into every later test. The
child imports the SAME opensees.pyd the parent loaded (directory passed in the
environment, child started with -S) so it cannot pick up a stale build through
a boot .pth.

  test_notice_and_iparm34            the one-time notice reports the branch MKL
                                     reads back, iparm(34) = thread count, ACTIVE
  test_no_flag_no_notice             default path: no notice (byte-identical)
  test_cbwr_explicit_branch          -cbwr COMPATIBLE is what MKL reports
  test_env_branch_kept               MKL_CBWR=COMPATIBLE + bare -deterministic:
                                     the launcher's branch is kept
  test_unknown_cbwr_degrades         -cbwr BOGUS warns and falls back to AUTO
  test_refused_after_lapack_warns_not_active
                                     eigen (LAPACK) first -> rc -8, relaunch hint,
                                     notice says CNR NOT ACTIVE
  test_prior_pardiso_solve_platform_rule
                                     a PARDISO solve first does NOT block it on
                                     Windows; on Linux it does (rc -8)
  test_mode_on_byte_identical        [slow] ~22k-DOF J2 push at 8 threads, N
                                     runs with the mode on: ONE distinct result
  test_mode_off_reported             [slow] the same N runs with the mode off:
                                     the distinct-result count is REPORTED, not
                                     asserted (drift is size/CPU-dependent, and a
                                     reproducible run with the mode off is not a
                                     failure of anything)

Measured 2026-09-27 (AMD Ryzen AI 7 PRO 350, 16 threads, MKL_NUM_THREADS=8,
N_RUNS=5, box shared with other builds): mode ON 1 distinct result / 5 runs,
706 s; mode OFF 5 distinct / 5 runs, 2 distinct curves, max relative spread
1.5e-16 (1 ULP) in the load factor, 1159 s. The slow pair is ~31 min.
"""
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from _testbed import ops  # noqa: E402  (the parent's pyd — its dir goes to the child)

pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="PARDISO requires MKL (Windows/oneAPI build)"),
]

PYD_DIR = str(Path(ops.__file__).resolve().parent)
N_RUNS = 5
THREADS = "8"

# ---------------------------------------------------------------------------
# The child: a J2 block pushed laterally under DisplacementControl. NX=NY=24,
# NZ=12 stdBricks -> 8125 nodes, 22 500 free DOFs (a few-thousand-DOF deck can
# be reproducible even with the mode off — Trap 7 is size-dependent).
# Prints ONE json line: sha256 of every nodal displacement after every step,
# plus the (lambda, u_top) curve and the iteration counts.
# ---------------------------------------------------------------------------
CHILD = r'''
import hashlib, json, os, struct, sys
d = os.environ["WP132_PYD_DIR"]
getattr(os, "add_dll_directory", lambda _d: None)(d)  # Windows-only API
sys.path.insert(0, d)
import opensees as ops
assert os.path.dirname(os.path.abspath(ops.__file__)) == os.path.abspath(d), ops.__file__

nx, ny, nz = [int(v) for v in os.environ["WP132_MESH"].split(",")]
nsteps = int(os.environ["WP132_STEPS"])
L, H = 1000.0, 500.0
def nid(i, j, k): return 1 + i + (nx + 1) * (j + (ny + 1) * k)
ops.wipe()
ops.model("basic", "-ndm", 3, "-ndf", 3)
# E=200000, nu=0.3 -> K, G; sigma_y 250, saturation 300, linear H 1000
ops.nDMaterial("J2Plasticity", 1, 166666.7, 76923.1, 250.0, 300.0, 10.0, 1000.0)
for k in range(nz + 1):
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(nid(i, j, k), i * L / nx, j * L / ny, k * H / nz)
for j in range(ny + 1):
    for i in range(nx + 1):
        ops.fix(nid(i, j, 0), 1, 1, 1)
t = 1
for k in range(nz):
    for j in range(ny):
        for i in range(nx):
            ops.element("stdBrick", t, nid(i, j, k), nid(i+1, j, k),
                        nid(i+1, j+1, k), nid(i, j+1, k), nid(i, j, k+1),
                        nid(i+1, j, k+1), nid(i+1, j+1, k+1), nid(i, j+1, k+1), 1)
            t += 1
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
for j in range(ny + 1):
    for i in range(nx + 1):
        ops.load(nid(i, j, nz), 1.0, 0.1, 0.0)
ctrl = nid(nx // 2, ny // 2, nz)
ops.constraints("Plain")
ops.numberer("RCM")
if os.environ.get("WP132_PARDISO_FIRST"):
    ops.system("Pardiso")
    ops.test("NormDispIncr", 1.0e-8, 40); ops.algorithm("Newton")
    ops.integrator("LoadControl", 0.0); ops.analysis("Static")
    assert ops.analyze(1) == 0
    ops.wipeAnalysis()
if os.environ.get("WP132_EIGEN_FIRST"):
    # a LAPACK solve BEFORE the request: MKL's dispatch is then initialized
    ops.mass(ctrl, 1.0, 1.0, 1.0)
    ops.eigen("-fullGenLapack", 1)
ops.system(*json.loads(os.environ["WP132_SYSTEM"]))
ops.test("NormDispIncr", 1.0e-8, 40)
ops.algorithm("Newton")
ops.integrator("DisplacementControl", ctrl, 1, float(os.environ["WP132_DU"]))
ops.analysis("Static")
h = hashlib.sha256()
curve, iters = [], []
tags = ops.getNodeTags()
for s in range(nsteps):
    ok = ops.analyze(1)
    if ok != 0:
        print(json.dumps({"ok": ok, "step": s})); sys.exit(0)
    iters.append(ops.testIter())
    for n in tags:
        h.update(struct.pack("<3d", *ops.nodeDisp(n)))
    curve.append([float.hex(ops.getLoadFactor(1)),
                  float.hex(ops.nodeDisp(ctrl, 1))])
print(json.dumps({"ok": 0, "sha": h.hexdigest(), "curve": curve, "iters": iters}))
'''


def _child(system_args, *, mesh=(4, 4, 2), steps=1, du=0.05, extra_env=None,
           threads=THREADS, timeout=900):
    env = dict(os.environ)
    env.pop("MKL_CBWR", None)          # every case states its own CNR setting
    env.update({
        "WP132_PYD_DIR": PYD_DIR,
        "WP132_MESH": ",".join(str(v) for v in mesh),
        "WP132_STEPS": str(steps),
        "WP132_DU": repr(du),
        "WP132_SYSTEM": json.dumps(list(system_args)),
        "MKL_NUM_THREADS": threads,
        "LADRUNO_OPENSEES_QUIET": "1",
    })
    if extra_env:
        env.update(extra_env)
    proc = subprocess.run([sys.executable, "-S", "-c", CHILD],
                          stdin=subprocess.DEVNULL, env=env, capture_output=True,
                          text=True, timeout=timeout)
    assert proc.returncode == 0, proc.stderr[-3000:]
    lines = [ln for ln in proc.stdout.splitlines() if ln.startswith("{")]
    assert lines, proc.stdout[-2000:] + proc.stderr[-2000:]
    out = json.loads(lines[-1])
    assert out["ok"] == 0, out
    return out, proc.stderr


def _notice(stderr):
    ln = [x for x in stderr.splitlines() if x.startswith("PARDISO deterministic mode:")]
    return ln


# ---------------------------------------------------------------- fast cases
def test_notice_and_iparm34():
    _, err = _child(["Pardiso", "-deterministic"])
    ln = _notice(err)
    assert len(ln) == 1, err            # ONE notice per solver object
    assert "CNR ACTIVE" in ln[0], ln[0]
    # AUTO; the "(-> BRANCH)" suffix appears only where MKL resolves one
    # (Intel). On AMD mkl_cbwr_get_auto_branch() returns AUTO itself.
    assert "branch AUTO" in ln[0] and "(-> AUTO)" not in ln[0], ln[0]
    assert "iparm(34)=" + THREADS + " thread(s)" in ln[0], ln[0]


def test_no_flag_no_notice():
    _, err = _child(["Pardiso"])
    assert not _notice(err), err


def test_cbwr_explicit_branch():
    _, err = _child(["Pardiso", "-cbwr", "COMPATIBLE"])
    ln = _notice(err)
    assert len(ln) == 1 and "branch COMPATIBLE," in ln[0] and "CNR ACTIVE" in ln[0], err


def test_env_branch_kept():
    _, err = _child(["Pardiso", "-deterministic"],
                    extra_env={"MKL_CBWR": "COMPATIBLE"})
    ln = _notice(err)
    assert len(ln) == 1 and "branch COMPATIBLE," in ln[0], err


def test_unknown_cbwr_degrades():
    _, err = _child(["Pardiso", "-cbwr", "BOGUS"])
    assert "unknown -cbwr branch BOGUS" in err, err
    ln = _notice(err)
    assert len(ln) == 1 and "branch AUTO" in ln[0] and "CNR ACTIVE" in ln[0], err


def test_refused_after_lapack_warns_not_active():
    # MKL_CBWR_ERR_MODE_CHANGE_FAILURE (-8): measured trigger = an earlier
    # LAPACK eigen solve. Python degrades (the Tcl parser returns TCL_ERROR).
    _, err = _child(["Pardiso", "-deterministic"],
                    extra_env={"WP132_EIGEN_FIRST": "1"})
    assert "failed (rc -8" in err and "MKL_CBWR=AUTO" in err, err
    ln = _notice(err)
    assert len(ln) == 1 and "CNR NOT ACTIVE" in ln[0], err


def test_prior_pardiso_solve_platform_rule():
    # Whether a PARDISO solve (no flag) earlier in the SAME process stops the
    # mode from being set afterwards depends on the platform. Both halves were
    # measured:
    #   Windows, oneMKL 2025.1 (mkl_intel_thread): it does NOT block.
    #   Linux (esmeralda, 2026-09-28): it DOES, rc -8, with oneMKL 2024.2 and
    #   2025.1 and with the gnu_thread and sequential layers alike. An earlier
    #   UmfPack solve blocks it too there, because UMFPACK calls MKL's BLAS.
    _, err = _child(["Pardiso", "-deterministic"],
                    extra_env={"WP132_PARDISO_FIRST": "1"})
    ln = _notice(err)
    assert len(ln) == 1, err
    if sys.platform == "win32":
        assert "CNR ACTIVE" in ln[0], err
    else:
        assert "failed (rc -8" in err and "CNR NOT ACTIVE" in ln[0], err


# ---------------------------------------------------------------- the F22 gate
BIG = dict(mesh=(24, 24, 12), steps=6, du=1.0)


def _spread(curves):
    """Largest relative spread across runs, over every step, of the load
    factor and of the control displacement: max|x - x_run0| / |x_run0|."""
    worst = [0.0, 0.0]
    for c in curves[1:]:
        for step, ref in zip(c, curves[0]):
            for q in (0, 1):
                a, b = float.fromhex(step[q]), float.fromhex(ref[q])
                if b != 0.0:
                    worst[q] = max(worst[q], abs(a - b) / abs(b))
    return worst


def _repeat(system_args):
    shas, curves, iters = [], [], []
    for _ in range(N_RUNS):
        out, _ = _child(system_args, **BIG)
        shas.append(out["sha"])
        curves.append(tuple(map(tuple, out["curve"])))
        iters.append(tuple(out["iters"]))
    return shas, curves, iters


@pytest.mark.slow
def test_mode_on_byte_identical():
    shas, curves, iters = _repeat(["Pardiso", "-deterministic"])
    print(f"\n[WP-132] mode ON : {len(set(shas))} distinct result(s) over "
          f"{N_RUNS} runs at MKL_NUM_THREADS={THREADS}; iters {iters[0]}")
    assert len(set(shas)) == 1, shas
    assert len(set(curves)) == 1, curves


@pytest.mark.slow
def test_mode_off_reported():
    shas, curves, iters = _repeat(["Pardiso"])
    print(f"\n[WP-132] mode OFF: {len(set(shas))} distinct result(s) over "
          f"{N_RUNS} runs at MKL_NUM_THREADS={THREADS}; "
          f"{len(set(curves))} distinct curve(s); iters {sorted(set(iters))}; "
          f"max rel. spread vs run 0: load factor {_spread(curves)[0]:.3e}, "
          f"control disp {_spread(curves)[1]:.3e}")
    # Reported, not asserted: see the module docstring.
