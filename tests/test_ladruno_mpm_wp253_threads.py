"""WP-253 (#1097, program #1075): mode F's threaded step is bit-identical to 1 thread.

Ladruno_implementation/253_mpm_performance.md §2 ("Determinism") and §5: `ladrunoThreads n`
threads mode F's transfer loops (supports, P2G as a node-centric gather in the serial summation
order, the grid update, G2P, the second projection, the particle update) and -- when every
particle material is on the WP-107 / WP-192 allowlist (Update AND Form) -- the material calls.
Every particle array the engine keeps (x, v, V, F, R, eps, sigHat, sig) after K steps must be
EXACTLY equal (float64 ==, no tolerance) at 2 and 4 threads to the 1-thread run, over:

* a 3D gimp elastic column on a frictional floor (the #1097 headline model, small);
* 3D bspline2 with `fix` faces (slaved ghosts and their images, the serial fold order);
* 2D plane strain with `mpm antiLocking fbar` (the corner term inside the gather; the 2D
  elastic material is NOT Form-safe -- class-wide statics -- so its calls stay serial);
* DruckerPrager (not allowlisted: serial material calls, threaded transfer);
* FLIP < 1 with local damping; `mpm gridvel`; a particle body force;
* failure paths: A7 (a particle carried out of the grid) and a commit refusal (W17) print the
  same error at n threads as at 1.

A loop-structure guard (`mpm debug threads` -> [transfer-loop, material-loop] threads of the
last step) asserts that 3D elastic threads both loops and DP / 2D elastic only the transfer.
One Tcl parity test (classic Tcl at 2 threads == classic Tcl at 1 thread, and == Python).
Particle counts and steps are set at run time from the model; no number is frozen.
Skipped when the binary has no OpenMP (nothing to thread).
"""
from __future__ import annotations

import json
import math
import os
import sys
import tempfile
from pathlib import Path

import pytest

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parent
sys.path.insert(0, str(TESTS))
import ladruno_mpm_implicit_harness as H  # noqa: E402
from _testbed.subprocess_run import run_python_script  # noqa: E402

_PROBE = (
    "import sys; sys.path.insert(0, %r)\n"
    "from _testbed import ops\n"
    "print('LADRUNO_PROBE_HW', ops.ladrunoThreads(1 << 20))\n"
) % str(TESTS)


def _openmp():
    try:
        code, out = run_python_script(_PROBE)
    except Exception:
        return True, ""
    if code != 0 or "LADRUNO_PROBE_HW" not in out:
        return True, ""
    if "built WITHOUT LADRUNO_OPENMP" in out:
        return False, "binary built with LADRUNO_OPENMP=OFF: nothing to thread"
    return True, ""


_OPENMP_IN, _SKIP_WHY = _openmp()
pytestmark = [pytest.mark.zone_a, pytest.mark.skipif(not _OPENMP_IN, reason=_SKIP_WHY)]

FIELDS = ("x", "v", "V", "F", "R", "epsHat", "sigHat", "sig")


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


@pytest.fixture(autouse=True)
def _one_thread_after():
    yield
    try:
        from _testbed import ops
        ops.ladrunoThreads(1)
    except Exception:
        pass


def _threads():
    ops = _ops()
    hw = int(ops.ladrunoThreads(1024))
    ops.ladrunoThreads(1)
    ns = sorted({n for n in (2, 4, min(8, hw)) if 2 <= n <= hw})
    if not ns:
        pytest.skip("this machine has one hardware thread")
    return ns


def _threads_ends():
    """The failure children: the smallest and the largest thread count (each is a process)."""
    ns = _threads()
    return sorted({ns[0], ns[-1]})


def _l(x):
    return [float(v) for v in x] if isinstance(x, (list, tuple)) else [float(x)]


# --------------------------------------------------------------------------- #
#  Models (command lists, so the Tcl test runs the same deck)                  #
# --------------------------------------------------------------------------- #
DP = ("DruckerPrager", 1, 1.6666666666666666e7, 7.6923076923076922e6, 50.0, 0.21, 0.0, 0.0,
      0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def col3d(basis="gimp", mat="elastic", flip=None, damping=None, antilocking=None, small=False):
    """The #1097 column (small): 3D, quarter symmetry, a frictional floor."""
    m = ("ElasticIsotropic", 1, 2.0e7, 0.3) if mat == "elastic" else DP
    cmds = [("wipe", []), ("model", ["basic", "-ndm", 3, "-ndf", 3]), ("nDMaterial", list(m)),
            ("mpm", ["grid", 232, "-origin", 0.0, 0.0, 0.0, "-h", 0.25, "-cells", 6, 4, 10,
                     "-basis", basis]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 0.0, 0.75, 0.5, 1.5, "-ppc", 2, "-mat", 1,
                     "-rho", 2000.0]),
            ("mpm", ["fix", "xmin", "x"]), ("mpm", ["fix", "ymin", "y"])]
    if basis == "gimp":
        cmds.append(("mpm", ["wall", "zmin", 0.3]))
    else:
        cmds.append(("mpm", ["fix", "zmin", "xyz"]))         # images on three faces
    cmds.append(("mpm", ["gravity", 0.0, 0.0, -10.0]))
    if flip is not None:
        cmds.append(("mpm", ["flip", flip]))
    if damping is not None:
        cmds.append(("mpm", ["damping", "local", damping]))
    if antilocking:
        cmds.append(("mpm", ["antiLocking", antilocking]))
    if small:
        cmds.append(("mpm", ["kinematics", "small"]))
    return cmds


def col2d(basis="gimp", antilocking="fbar", gridvel=False, mat="elastic"):
    m = ("ElasticIsotropic", 1, 2.0e7, 0.45) if mat == "elastic" else DP
    cmds = [("wipe", []), ("model", ["basic", "-ndm", 2, "-ndf", 2]), ("nDMaterial", list(m)),
            ("mpm", ["grid", 232, "-origin", 0.0, 0.0, "-h", 0.25, "-cells", 10, 12,
                     "-basis", basis]),
            ("mpm", ["body", 1, "box", 0.0, 0.0, 1.0, 2.0, "-ppc", 2, "-mat", 1, "-rho", 2000.0]),
            ("mpm", ["fix", "xmin", "x"]), ("mpm", ["wall", "ymin", 0.3]),
            ("mpm", ["gravity", 0.0, -10.0])]
    if antilocking:
        cmds.append(("mpm", ["antiLocking", antilocking]))
    if gridvel:
        cmds.append(("mpm", ["gridvel", "push", "box", 0.9, 1.0, 3.0, 3.5, "v", -0.2, "free"]))
    return cmds


def _body_force(ops):
    n = int(ops.mpm("info", "numParticles"))
    nd = len(_l(ops.mpm("particleData", "x"))) // n
    bf = [0.0] * (nd * n)
    for p in range(n):
        bf[nd * p] = 30.0 * math.cos(0.37 * p)
    ops.mpm("particleBodyForce", *bf)


def steps_at_quarter_dtcrit(ops, k):
    """k single steps, each at 0.25 x the dt_crit of its own start (the particles accelerate,
    so a dt fixed at step 0 can exceed s x dt_crit later). The same dt sequence at every thread
    count, because dt_crit is itself bit-identical."""
    for _ in range(k):
        ops.mpm("analyze", 1, 0.25 * float(ops.mpm("criticalTimeStep")))


def run(cmds, nthreads, post=None, steps=None, deep=False):
    """Build, set the thread count, run K steps at 0.25 dt_crit (under the -deep profiler when
    asked); the particle arrays, the loop-structure guard, and the step count."""
    ops = _ops()
    for name, args in cmds:
        getattr(ops, name)(*H._py_args(args))
    if post is not None:
        post(ops)
    ops.ladrunoThreads(nthreads)
    n = int(ops.mpm("info", "numParticles"))
    k = steps if steps is not None else max(12, min(40, 20000 // max(n, 1)))
    x0 = _l(ops.mpm("particleData", "x"))
    if deep:
        ops.profiler("reset")
        ops.profiler("start", "-deep")
    try:
        steps_at_quarter_dtcrit(ops, k)
    finally:
        if deep:
            ops.profiler("stop")
            ops.profiler("reset")
    data = {f: _l(ops.mpm("particleData", f)) for f in FIELDS}
    data["moved"] = max(abs(a - b) for a, b in zip(data["x"], x0))
    thr = [int(v) for v in _l(ops.mpm("debug", "threads"))]
    ops.ladrunoThreads(1)
    return data, thr, n


CASES = {
    "3d_gimp_elastic_wall": dict(cmds=lambda: col3d("gimp"), both=True),
    "3d_bspline2_fix_images": dict(cmds=lambda: col3d("bspline2"), both=True),
    "2d_plane_strain_fbar": dict(cmds=lambda: col2d("gimp", "fbar"), both=False),
    "2d_bspline2_fbar": dict(cmds=lambda: col2d("bspline2", "fbar"), both=False),
    "3d_druckerprager": dict(cmds=lambda: col3d("gimp", mat="dp"), both=False),
    "3d_flip_damping": dict(cmds=lambda: col3d("gimp", flip=0.7, damping=0.1), both=True),
    "2d_gridvel": dict(cmds=lambda: col2d("gimp", None, gridvel=True), both=False),
    "3d_body_force": dict(cmds=lambda: col3d("bspline2"), both=True, post=_body_force),
    "3d_gimp_elastic_fbar": dict(cmds=lambda: col3d("gimp", antilocking="fbar"), both=True),
    "3d_small_kinematics": dict(cmds=lambda: col3d("gimp", small=True), both=True),
    # the -deep profiler: threaded material calls fold their walls per particle (afterParallel)
    "3d_elastic_deep_profiled": dict(cmds=lambda: col3d("gimp"), both=True, deep=True),
}


@pytest.mark.parametrize("name", sorted(CASES))
def test_threads_bit_identical(name):
    """Every particle array after K steps: float64 == at 2 and 4 threads vs 1 thread."""
    c = CASES[name]
    ref, thr1, n = run(c["cmds"](), 1, c.get("post"))      # the reference is never profiled
    assert thr1 == [1, 1], thr1
    assert n > 50, n
    assert ref["moved"] > 0.0, "the model must move for the comparison to mean anything"
    for nt in _threads():
        got, thr, _ = run(c["cmds"](), nt, c.get("post"), deep=c.get("deep", False))
        assert thr == [nt, nt if c["both"] else 1], (name, nt, thr)
        for f in FIELDS:
            assert len(got[f]) == len(ref[f]), (name, nt, f)
            bad = [i for i, (a, b) in enumerate(zip(got[f], ref[f])) if a != b]
            assert not bad, f"{name}: {f} differs at {nt} threads in {len(bad)} entries (first {bad[0]})"


def test_loop_structure_guard():
    """`mpm debug threads`: 3D elastic threads the transfer AND the material loop; DP and the 2D
    elastic material (class-wide statics, not Form-safe) thread only the transfer loops."""
    nt = _threads()[-1]
    assert run(col3d("gimp"), nt, steps=1)[1] == [nt, nt]
    assert run(col3d("gimp", mat="dp"), nt, steps=1)[1] == [nt, 1]
    assert run(col2d("gimp", None), nt, steps=1)[1] == [nt, 1]
    assert run(col3d("gimp"), 1, steps=1)[1] == [1, 1]


# --------------------------------------------------------------------------- #
#  Failure paths -- in a child process, so the error text is captured          #
# --------------------------------------------------------------------------- #
_CHILD = r'''
import json, sys
sys.path.insert(0, sys.argv[1])
import ladruno_mpm_implicit_harness as H
from _testbed import ops
cmds, nt, kind = json.loads(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
for name, args in cmds:
    getattr(ops, name)(*H._py_args(args))
ops.ladrunoThreads(nt)
def step():
    ops.mpm("analyze", 1, 0.25 * float(ops.mpm("criticalTimeStep")))
def state(tag):
    print(tag, json.dumps({f: [float(v) for v in ops.mpm("particleData", f)]
                           for f in ("x", "v", "F", "sig", "sigHat", "epsHat")}))
def failing(f):
    try:
        f()
        print("NO-FAILURE")
    except Exception:
        pass
ids = [int(v) for v in ops.mpm("particleData", "id")]
n = len(ids)
# two particles far apart in storage order: in different blocks at every thread count
chosen = [ids[(2 * n) // 3], ids[n // 3]]
print("LOWEST", min(chosen))
if kind == "a7":
    ops.mpm("velocity", 1, "uniform", 40.0, 0.0, 0.0)
    step(); step()
    failing(lambda: [step() for _ in range(400)])
    state("STATE")
elif kind == "inject":
    step(); step()
    ops.mpm("debug", "injectRefusal")
    failing(step)
    state("STATE")
elif kind == "failTrial":
    step(); step()
    ops.mpm("debug", "failTrial", *chosen)
    failing(step)
    state("STATE")                       # the committed state, untouched
    ops.mpm("debug", "failTrial", "off")
    step(); step()                       # the retried steps
    state("RETRY")
elif kind == "failCommit":
    step(); step()
    ops.mpm("debug", "failCommit", *chosen)
    failing(step)
    state("STATE")                       # W17: particles not moved
    ops.mpm("debug", "failCommit", "off")
    ops.mpm("reset")
    step(); step()
    state("RETRY")
print("THREADS", [int(v) for v in ops.mpm("debug", "threads")])
'''


def _child(cmds, nt, kind):
    code, out = run_python_script(_CHILD, argv=(str(TESTS), json.dumps(cmds), str(nt), kind))
    assert code == 0, out[-3000:]
    assert "NO-FAILURE" not in out, out[-3000:]
    lines = out.splitlines()
    errors = [ln.strip() for ln in lines if "ERROR mpm analyze" in ln]
    tagged = {ln.split(" ", 1)[0]: ln.split(" ", 1)[1] for ln in lines
              if ln.split(" ", 1)[0] in ("STATE", "RETRY", "LOWEST", "THREADS")}
    return errors, tagged


_EXPECT = {
    "a7": "out of the grid",
    "inject": "REFUSED this commit",
    "failTrial": "failed to integrate the trial strain",
    "failCommit": "2 particle material point(s) REFUSED this commit",
}


@pytest.mark.parametrize("mat", ["elastic", "dp"])     # threaded material calls / the split path
@pytest.mark.parametrize("kind", sorted(_EXPECT))
def test_failure_text_is_the_serial_one(kind, mat):
    '''A7, a W17 refusal (injected; and declared by two chosen particles' commits) and two
    chosen particles' trial failures: at n threads the same error line naming the LOWEST id, the
    same particle state after the failed step, and the same state after a retry, as at 1 thread
    -- on both material paths (3D elastic threads its material calls, DruckerPrager does not).'''
    cmds = col3d("gimp", mat=mat)
    err1, t1 = _child(cmds, 1, kind)
    assert err1 and _EXPECT[kind] in err1[-1], (kind, mat, err1)
    lowest = t1["LOWEST"]
    if kind == "failTrial":
        assert f"the material of particle {lowest} failed" in err1[-1], err1
        assert f"(particle id {lowest}," in err1[-1], err1
    if kind == "failCommit":
        assert f"first refusing: particle id {lowest}," in err1[-1], err1
    for nt in _threads_ends():
        errn, tn = _child(cmds, nt, kind)
        assert errn == err1, (kind, mat, nt, errn, err1)
        for tag in ("STATE", "RETRY"):
            assert tn.get(tag) == t1.get(tag), (kind, mat, nt, tag)
        if kind != "a7":
            # the failed / retried steps ran threaded (a7's last step is the refused one too)
            want = [nt, nt if mat == "elastic" else 1]
            assert json.loads(tn["THREADS"]) == want, (kind, mat, nt, tn["THREADS"])


# --------------------------------------------------------------------------- #
#  Classic Tcl                                                                 #
# --------------------------------------------------------------------------- #
def _tcl_exe():
    exe = os.environ.get("LADRUNO_TCL_EXE")
    if exe:
        return exe if Path(exe).exists() else None
    cand = REPO / "dist" / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
    return str(cand) if cand.exists() else None


def _tcl(cmds, nt, k):
    exe = _tcl_exe()
    if exe is None:
        pytest.skip("classic OpenSees exe not found (build it, or set LADRUNO_TCL_EXE)")
    lines = [name + " " + " ".join(H._tcl_arg(a) for a in args) for name, args in cmds]
    lines += [f"ladrunoThreads {nt}",
              f"for {{set i 0}} {{$i < {k}}} {{incr i}} {{ mpm analyze 1 [expr {{0.25 * [mpm criticalTimeStep]}}] }}"]
    lines += [f'puts "P_{f} [mpm particleData {f}]"' for f in FIELDS]
    lines += ['puts "THREADS [mpm debug threads]"', 'puts "DECK-DONE"']
    with tempfile.TemporaryDirectory(ignore_cleanup_errors=True) as td:
        out = H._run_deck("\n".join(lines) + "\n", exe, td, f"wp253_t{nt}", 600)
    res = {}
    for line in out.splitlines():
        t = line.split()
        if t and (t[0].startswith("P_") or t[0] == "THREADS"):
            res[t[0]] = [float(v) for v in t[1:]]
    return res


def test_tcl_threads_parity():
    """Classic Tcl: 2 threads == 1 thread exactly (the same print path), and == Python's
    1-thread run to the print precision (as the other MPM Tcl parity tests)."""
    nt = _threads()[0]
    cmds = col3d("gimp")
    ops = _ops()
    for name, args in cmds:
        getattr(ops, name)(*H._py_args(args))
    k = 15
    steps_at_quarter_dtcrit(ops, k)
    py = {f: _l(ops.mpm("particleData", f)) for f in FIELDS}
    t1, tn = _tcl(cmds, 1, k), _tcl(cmds, nt, k)
    assert tn["THREADS"] == [nt, nt], tn["THREADS"]
    for f in FIELDS:
        assert tn["P_" + f] == t1["P_" + f], f
        scale = max(abs(v) for v in py[f]) or 1.0
        assert len(t1["P_" + f]) == len(py[f])
        assert max(abs(a - b) for a, b in zip(t1["P_" + f], py[f])) <= 1e-14 * scale, f
