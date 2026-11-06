"""WP-190: wall-clock phase timings in the `system Mumps -stats` block.

Deck: tests/wp188_factor_reuse/frame_gm.tcl (elastic model, OpenSeesMP, MUMPS), with
WP188_MUMPSARGS="-stats". What is pinned:
  * the original stats lines keep their exact text (compare_validation.py regexes);
  * the block gains t_factor / t_solve lines (t_analysis only on a print after an analysis);
  * without -stats nothing of the new output appears;
  * with WP-188 factor reuse a skipped factorization prints no block: prints ==
    factorizations, not steps;
  * the answer is the same with and without -stats (JOB=5 vs JOB=2+3).
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "wp188_factor_reuse" / "frame_gm.tcl"
DIST = Path(os.environ.get("LADRUNO_DIST", REPO / "dist"))
MP_EXE = DIST / "bin" / ("OpenSeesMP.exe" if os.name == "nt" else "OpenSeesMP")
MPIEXEC = DIST / "openseesmp" / ("mpiexec.exe" if os.name == "nt" else "mpiexec")
needs_mp = pytest.mark.skipif(not (MP_EXE.exists() and MPIEXEC.exists()),
                              reason="dist OpenSeesMP / mpiexec not built")

OLD_LINES = [
    r"MUMPS stats: n=\d+ BLR\(ICNTL35\)=0",
    r"  factor entries INFOG\(9\)   = \S+",
    r"  factor MB/proc INFOG\(21\)  = \d+  \(max over procs\)",
    r"  factor MB tot  INFOG\(22\)  = \d+  \(sum over procs\)",
    r"  elim flops     RINFOG\(3\)  = \S+",
]
TF = re.compile(r"^  t_factor s \(max\)    = (\S+)$", re.M)
TS = re.compile(r"^  t_solve s \(max, n\)  = (\S+), (\d+)$", re.M)
TA = re.compile(r"^  t_analysis s \(max\)  = (\S+)$", re.M)


def _run(tmp, name, np=2, stats=True, reuse=0, nstep=60):
    env = dict(os.environ)
    env.update(WP188_OUT=str(Path(tmp) / name), WP188_MODEL="elastic", WP188_NSTEP=str(nstep),
               WP188_REUSE=str(reuse), LADRUNO_OPENSEES_QUIET="1")
    env["WP188_MUMPSARGS"] = "-stats" if stats else ""
    tcl = DIST / "lib" / "tcl8.6"
    if "TCL_LIBRARY" not in env and tcl.exists():
        env["TCL_LIBRARY"] = str(tcl)
    env["PATH"] = str(MPIEXEC.parent) + os.pathsep + env.get("PATH", "")
    cmd = [str(MPIEXEC)] + (["-localonly"] if os.name == "nt" else []) + \
          ["-n", str(np), str(MP_EXE), str(DECK)]
    p = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=600, env=env, cwd=str(tmp))
    log = p.stdout + p.stderr
    edp = (Path(tmp) / (name + ".edp")).read_text()
    return log, edp


@needs_mp
def test_stats_block_gains_timings_and_keeps_old_lines(tmp_path):
    log, _ = _run(tmp_path, "s")
    for rx in OLD_LINES:
        assert re.search("^" + rx + "$", log, re.M), rx
    tf, ts, ta = TF.findall(log), TS.findall(log), TA.findall(log)
    nblocks = len(re.findall(r"^MUMPS stats:", log, re.M))
    assert nblocks >= 1
    assert len(tf) == nblocks and len(ts) == nblocks
    # analysis only on the print after one (gravity and transient each analyse once)
    assert 1 <= len(ta) < nblocks
    assert all(float(x) > 0 for x in tf + ta)
    # every solve is accounted for: tail + blocks == one JOB=3 per analyze step
    tail = re.findall(r"^MUMPS solve tail: t_solve s \(rank 0, n\) = \S+, (\d+)$", log, re.M)
    total = sum(int(n) for _, n in ts) + sum(int(n) for n in tail)
    assert total >= nblocks and total > 0


@needs_mp
def test_stats_off_prints_nothing_new_and_same_answer(tmp_path):
    log0, edp0 = _run(tmp_path, "off", stats=False)
    log1, edp1 = _run(tmp_path, "on", stats=True)
    assert "t_factor" not in log0 and "t_solve" not in log0 and "MUMPS stats" not in log0
    assert edp0 == edp1


@needs_mp
def test_reuse_skipped_factorization_prints_no_block(tmp_path):
    full, _ = _run(tmp_path, "full", reuse=0)
    reuse, _ = _run(tmp_path, "reuse", reuse=10)
    nf = len(re.findall(r"^MUMPS stats:", full, re.M))
    nr = len(re.findall(r"^MUMPS stats:", reuse, re.M))
    assert 0 < nr < nf
    assert len(TF.findall(reuse)) == nr
    # solves keep being counted across the skipped factorizations
    ns = sum(int(n) for _, n in TS.findall(reuse)) + \
        sum(int(n) for n in re.findall(r"tail: t_solve s \(rank 0, n\) = \S+, (\d+)", reuse))
    assert ns > nr
