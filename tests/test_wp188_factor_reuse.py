"""WP-188: factorization reuse across time steps (`algorithm KrylovNewton -reuseFactor m`).

ADR: Ladruno_implementation/188_ladruno_factor_reuse_adr.md.

The deck is tests/wp188_factor_reuse/frame_gm.tcl: a 3-story, 2-bay RC frame
(fiber forceBeamColumn, Concrete01/Steel02, P-Delta columns) or its elastic twin,
gravity by LoadControl and then 10 s of a synthetic ground motion (Newmark,
dt = 0.01, 1000 steps), run by the classic OpenSees.exe and by OpenSeesMP.exe.
The reuse statistics are parsed from the summary line the algorithm prints
(`print -algorithm` per rank, and its destructor on rank 0).

What is pinned:
  * option off == the vanilla path: no -reuseFactor and `-reuseFactor 0` give the
    same bytes; with LADRUNO_WP188_BASE_DIST pointing at a dist/ built from the
    base commit, the option-absent outputs are md5-equal to it as well;
  * linear model: reuse equals full Newton to round-off and forms exactly
    1 + ceil((N - 1) / (m + 1)) factors (step 1, then one per m + 1 steps; the
    first step's dt is unknown to the algorithm, so step 2 forms too);
  * nonlinear model: converges, > 90 % of the steps start on a kept factor, the
    EDPs stay within the test-tolerance scale of full Newton; a dt change forces a
    factorization; the iteration and rate triggers fire; a forced non-convergence
    (test max 7 iterations, no in-step refresh) is undone and retried;
  * the option parsers fail closed (classic Tcl and openseespy);
  * OpenSeesMP on 2 and 4 ranks (MUMPS): completes, every rank reports the same
    statistics, results match the no-reuse MP run within the tolerance scale, and
    the reduction is live (a rank that votes "form" makes every rank form).

Measured wall time: about 40 s for the whole module on the 24-core desktop
(serial runs ~0.7 s each, MP runs ~2 s each).
"""
import hashlib
import math
import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "wp188_factor_reuse" / "frame_gm.tcl"
DIST = REPO / "dist"
SER_EXE = DIST / "bin" / ("OpenSees.exe" if os.name == "nt" else "OpenSees")
MP_EXE = DIST / "bin" / ("OpenSeesMP.exe" if os.name == "nt" else "OpenSeesMP")
MPIEXEC = DIST / "openseesmp" / ("mpiexec.exe" if os.name == "nt" else "mpiexec")

needs_serial = pytest.mark.skipif(not SER_EXE.exists(), reason="dist OpenSees not built")
needs_mp = pytest.mark.skipif(not (MP_EXE.exists() and MPIEXEC.exists()),
                              reason="dist OpenSeesMP / mpiexec not built")

SUM = re.compile(
    r"factor reuse \(WP-188\) ([\w ]+?)(?: \[rank (\d+)\])?: steps (\d+), started on a kept factor (\d+) "
    r"\(([\d.e+-]+) %\), factorizations (\d+) = start (\d+) \+ iter (\d+) \+ rate (\d+) \+ krylov (\d+) "
    r"\+ retry (\d+), retries (\d+) \(failed (\d+)\), iterations (\d+) \(on kept factors (\d+)\)")
KEYS = ("steps", "reused", "pct", "fact", "start", "iter", "rate", "krylov", "retry",
        "retries", "rfailed", "its", "its_kept")


def _summaries(text):
    out = []
    for m in SUM.finditer(text):
        d = {"when": m.group(1), "rank": None if m.group(2) is None else int(m.group(2))}
        for k, v in zip(KEYS, m.groups()[2:]):
            d[k] = float(v) if k == "pct" else int(v)
        out.append(d)
    return out


def _run(tmp, name, dist=DIST, np=1, extra_env=None, **opts):
    tmp = Path(tmp)
    env = dict(os.environ)
    env["WP188_OUT"] = str(tmp / name)
    env["LADRUNO_OPENSEES_QUIET"] = "1"
    for k, v in opts.items():
        env["WP188_" + k.upper()] = str(v)
    if extra_env:
        env.update(extra_env)
    tcl = Path(dist) / "lib" / "tcl8.6"
    if "TCL_LIBRARY" not in env and tcl.exists():
        env["TCL_LIBRARY"] = str(tcl)
    if np == 1:
        cmd = [str(Path(dist) / "bin" / SER_EXE.name), str(DECK)]
    else:
        mpi = Path(dist) / "openseesmp" / MPIEXEC.name
        env["PATH"] = str(mpi.parent) + os.pathsep + env.get("PATH", "")
        cmd = [str(mpi)] + (["-localonly"] if os.name == "nt" else []) + \
              ["-n", str(np), str(Path(dist) / "bin" / MP_EXE.name), str(DECK)]
    p = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600, env=env, cwd=str(tmp))
    log = p.stdout + p.stderr
    edp_f = tmp / (name + ".edp")
    assert edp_f.exists(), f"{name}: no EDP file (rc={p.returncode})\n{log[-3000:]}"
    toks = edp_f.read_text().split()
    edp = {toks[i]: float(toks[i + 1]) for i in range(0, len(toks), 2)}
    roof_f = tmp / (name + ".roof")
    roof = [float(x) for x in roof_f.read_text().split()]
    md5 = hashlib.md5(edp_f.read_bytes() + roof_f.read_bytes()).hexdigest()
    algo = []
    for r in range(np):
        f = tmp / (name + ".algo%d" % r)
        s = _summaries(f.read_text(errors="replace")) if f.exists() else []
        algo.append(s[-1] if s else None)
    return {"rc": p.returncode, "log": log, "edp": edp, "roof": roof, "md5": md5, "algo": algo}


def _hist_diff(a, b):
    return max(abs(x - y) for x, y in zip(a["roof"], b["roof"])) / b["edp"]["maxRoof"]


def _edp_diff(a, b):
    ref = b["edp"]
    d = [abs(a["edp"][k] - ref[k]) / abs(ref[k]) for k in ("maxRoof", "maxIdr")]
    d.append(abs(a["edp"]["resid"] - ref["resid"]) / ref["maxRoof"])
    return max(d)


# ----------------------------------------------------------------------------- serial

@needs_serial
def test_option_off_is_the_vanilla_path(tmp_path):
    a = _run(tmp_path, "absent")
    z = _run(tmp_path, "zero", reuse=0)
    assert a["edp"]["fail"] == 0
    assert a["md5"] == z["md5"]
    assert a["algo"][0] is None and z["algo"][0] is None      # no reuse summary when off
    base = os.environ.get("LADRUNO_WP188_BASE_DIST")
    if base:
        for opts in ({}, {"tol": 1e-4}, {"model": "elastic"}, {"halve": 300}):
            b = _run(tmp_path, "b%d" % len(opts), dist=base, **opts)
            n = _run(tmp_path, "n%d" % len(opts), **opts)
            assert b["md5"] == n["md5"], f"option-absent output differs from the base build for {opts}"


@needs_serial
@pytest.mark.parametrize("m", [4, 1000])
def test_linear_reuse_equals_full_newton(tmp_path, m):
    full = _run(tmp_path, "full", model="elastic", tol=1e-10)
    reuse = _run(tmp_path, "reuse", model="elastic", tol=1e-10, reuse=m)
    assert reuse["edp"]["fail"] == 0
    assert _hist_diff(reuse, full) < 1e-12
    s = reuse["algo"][0]
    n = int(reuse["edp"]["steps"])
    assert s["steps"] == n
    assert s["start"] == 1 + math.ceil((n - 1) / (m + 1))
    assert s["fact"] == s["start"] and s["retries"] == 0


@needs_serial
def test_nonlinear_frame_reuse(tmp_path):
    full = _run(tmp_path, "full", tol=1e-6, maxdim=20)
    reuse = _run(tmp_path, "reuse", tol=1e-6, maxdim=20, reuse=30)
    assert reuse["edp"]["fail"] == 0 and reuse["edp"]["steps"] == 1000
    s = reuse["algo"][0]
    assert s["pct"] > 90.0
    assert s["fact"] < 0.1 * s["steps"]
    assert s["retries"] == 0
    # converged to another point of the same tolerance ball: differences are of
    # the order of the test tolerance (1e-6 m on a 0.2 m roof), far above round-off
    d = max(_edp_diff(reuse, full), _hist_diff(reuse, full))
    assert 1e-12 < d < 1e-4, d


@needs_serial
def test_dt_change_forces_a_factorization(tmp_path):
    steady = _run(tmp_path, "steady", tol=1e-6, maxdim=20, reuse=1000)
    halved = _run(tmp_path, "halved", tol=1e-6, maxdim=20, reuse=1000, halve=300)
    assert steady["algo"][0]["start"] == 2
    # halving forms once, restoring dt forms once more
    assert halved["algo"][0]["start"] == 4
    assert halved["edp"]["fail"] == 0


@needs_serial
def test_iteration_and_rate_triggers(tmp_path):
    it = _run(tmp_path, "iter", tol=1e-6, maxdim=20, reuse=100, iter=3)
    rt = _run(tmp_path, "rate", tol=1e-6, maxdim=20, reuse=100, rate=0.3)
    assert it["edp"]["fail"] == 0 and rt["edp"]["fail"] == 0
    assert it["algo"][0]["iter"] > 0 and it["algo"][0]["rate"] == 0
    assert rt["algo"][0]["rate"] > 0 and rt["algo"][0]["iter"] == 0


@needs_serial
def test_forced_nonconvergence_is_retried(tmp_path):
    # full Newton needs at most 6 iterations per step here; with max 7 and no
    # in-step refresh (-reuseIter 100 -reuseRate 0) many kept-factor steps fail
    full = _run(tmp_path, "full", tol=1e-8, maxit=7, maxdim=20)
    reuse = _run(tmp_path, "reuse", tol=1e-8, maxit=7, maxdim=20, reuse=100, iter=100, rate=0)
    assert full["edp"]["fail"] == 0
    assert reuse["edp"]["fail"] == 0 and reuse["edp"]["steps"] == 1000
    s = reuse["algo"][0]
    assert s["retries"] > 0 and s["retry"] == s["retries"] and s["rfailed"] == 0
    assert "undoing it and retrying with a fresh tangent" in reuse["log"]
    assert _hist_diff(reuse, full) < 1e-5


_REFUSE_DECK = r"""
model BasicBuilder -ndm 1 -ndf 1
node 1 0.0
node 2 1.0
fix 1 1
uniaxialMaterial Elastic 1 100.0
element truss 1 1 2 1.0 1
test NormDispIncr 1e-8 10
set r [catch {algorithm KrylovNewton %s} msg]
puts "RESULT $r"
"""


@needs_serial
@pytest.mark.parametrize("args,refused", [
    ("-reuseFactor 5", 0),
    ("-reuseFactor 5 -reuseIter 4 -reuseRate 0.5 -reuseReport 100", 0),
    ("-reuseIter 4", 1),
    ("-reuseFactor -1", 1),
    ("-reuseFactor", 1),
    ("-reuseFactor 5 -reuseRate x", 1),
    ("-reuseFactor 5 -reuseIter 0", 1),
])
def test_classic_parser_fails_closed(tmp_path, args, refused):
    deck = tmp_path / "refuse.tcl"
    deck.write_text(_REFUSE_DECK % args)
    env = dict(os.environ)
    tcl = DIST / "lib" / "tcl8.6"
    if "TCL_LIBRARY" not in env and tcl.exists():
        env["TCL_LIBRARY"] = str(tcl)
    p = subprocess.run([str(SER_EXE), str(deck)], stdin=subprocess.DEVNULL, capture_output=True,
                       text=True, encoding="utf-8", errors="replace", timeout=60, env=env)
    m = re.search(r"RESULT (\d)", p.stdout + p.stderr)
    assert m, p.stdout + p.stderr
    assert int(m.group(1)) == refused, p.stdout + p.stderr


def test_openseespy_parser_and_linear_reuse():
    """The interpreter path (OPS_KrylovNewton: openseespy / openseesmp)."""
    from _testbed import ops

    def build():
        ops.wipe()
        ops.model("basic", "-ndm", 1, "-ndf", 1)
        ops.node(1, 0.0); ops.node(2, 1.0)
        ops.fix(1, 1)
        ops.mass(2, 1.0)
        ops.uniaxialMaterial("Elastic", 1, 100.0)
        ops.element("truss", 1, 1, 2, 1.0, 1)
        ops.timeSeries("Sine", 1, 0.0, 10.0, 0.5)
        ops.pattern("Plain", 1, 1)
        ops.load(2, 1.0)
        ops.constraints("Plain"); ops.numberer("Plain"); ops.system("BandGeneral")
        ops.test("NormDispIncr", 1e-12, 20)

    def drive(*algo):
        build()
        ops.algorithm("KrylovNewton", *algo)
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
        out = []
        for _ in range(50):
            assert ops.analyze(1, 0.01) == 0
            out.append(ops.nodeDisp(2, 1))
        return out

    full = drive()
    reuse = drive("-reuseFactor", 10)
    assert max(abs(a - b) for a, b in zip(full, reuse)) < 1e-14
    for bad in (("-reuseIter", 3), ("-reuseFactor", -2), ("-reuseFactor", 3, "-reuseRate", -1.0)):
        build()
        with pytest.raises(Exception):
            ops.algorithm("KrylovNewton", *bad)
    ops.wipe()


# ----------------------------------------------------------------------------- MPI

@needs_mp
@pytest.mark.parametrize("np", [2, 4])
def test_mp_reuse_ranks_agree(tmp_path, np):
    full = _run(tmp_path, "full", np=np, tol=1e-6, maxdim=20)
    reuse = _run(tmp_path, "reuse", np=np, tol=1e-6, maxdim=20, reuse=30)
    assert full["edp"]["fail"] == 0 and reuse["edp"]["fail"] == 0
    assert reuse["edp"]["steps"] == 1000
    sums = reuse["algo"]
    assert all(s is not None for s in sums), "a rank wrote no reuse summary"
    strip = [tuple((k, v) for k, v in s.items() if k not in ("rank", "when")) for s in sums]
    assert len(set(strip)) == 1, f"ranks disagree: {sums}"
    assert sums[0]["pct"] > 90.0
    assert _hist_diff(reuse, full) < 1e-4


@needs_mp
def test_mp_retry_and_reduction_is_live(tmp_path):
    # forced non-convergence under MUMPS on 2 ranks
    retry = _run(tmp_path, "retry", np=2, tol=1e-8, maxit=7, maxdim=20, reuse=100, iter=100, rate=0)
    assert retry["edp"]["fail"] == 0
    s0, s1 = retry["algo"]
    assert s0["retries"] > 0 and s0["rfailed"] == 0
    assert {k: v for k, v in s0.items() if k != "rank"} == {k: v for k, v in s1.items() if k != "rank"}
    # rank 1 votes "form" at every step start: rank 0 can only learn that through
    # the MPI_Allreduce, so rank 0 forming at every step proves the reduction is
    # compiled into OpenSeesMP (not the dead OBJECT-library #ifdef)
    vote = _run(tmp_path, "vote", np=2, tol=1e-6, maxdim=20, reuse=30,
                extra_env={"LADRUNO_WP188_TEST_VOTE_RANK": "1"})
    v0, v1 = vote["algo"]
    assert vote["edp"]["fail"] == 0
    assert v0["reused"] == 0 and v1["reused"] == 0
    assert v0["start"] == v0["steps"] == 1000
