"""WP-194: mixed-precision MUMPS, `system Mumps ... -precision single [-refine k]`.

Deck: tests/wp194_mumps_mixed/mixed_frame.tcl, which runs the WP-188 RC frame
(tests/wp188_factor_reuse/frame_gm.tcl: fiber forceBeamColumn, Concrete01/Steel02,
KrylovNewton, NormDispIncr 1e-8) on OpenSeesMP with 2 ranks and appends
$WP194_MUMPSARGS to its `system Mumps`. What is pinned:
  * default off: no option and `-precision double` run the same path (same
    iterations; histories equal up to the 2-rank MUMPS run-to-run last-bit noise);
  * a build without single-precision MUMPS refuses `-precision single` (no silent
    double factorization); bad values and `-refine` without single are refused;
  * with SMUMPS: the nonlinear frame converges with `-precision single` to the
    same ConvergenceTest, the roof history agrees with double to well inside the
    tolerance-driven noise, and the iteration count is printed (single vs double);
  * `-refine k` drives the double residual down (the `-stats` residual lines);
  * it coexists with `-BLR` and with WP-188 `-reuseFactor`.
"""
import os
import re
import subprocess
from pathlib import Path

import pytest

pytestmark = [pytest.mark.zone_a]

REPO = Path(__file__).resolve().parents[1]
DECK = REPO / "tests" / "wp194_mumps_mixed" / "mixed_frame.tcl"
DIST = Path(os.environ.get("LADRUNO_DIST", REPO / "dist"))
MP_EXE = Path(os.environ.get("WP194_MP_EXE", DIST / "bin" / ("OpenSeesMP.exe" if os.name == "nt" else "OpenSeesMP")))
MPIEXEC = Path(os.environ.get("WP194_MPIEXEC", DIST / "openseesmp" / ("mpiexec.exe" if os.name == "nt" else "mpiexec")))
needs_mp = pytest.mark.skipif(not (MP_EXE.exists() and MPIEXEC.exists()),
                              reason="dist OpenSeesMP / mpiexec not built")

REFUSED = "needs single-precision MUMPS"
NSTEP = 300


def _run(tmp, name, mumpsargs="", model="fiber", reuse=0, nstep=NSTEP, np=2):
    env = dict(os.environ)
    env.update(WP188_OUT=str(Path(tmp) / name), WP188_MODEL=model, WP188_NSTEP=str(nstep),
               WP188_REUSE=str(reuse), WP194_MUMPSARGS=mumpsargs, LADRUNO_OPENSEES_QUIET="1",
               OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    tcl = DIST / "lib" / "tcl8.6"
    if "TCL_LIBRARY" not in env and tcl.exists():
        env["TCL_LIBRARY"] = str(tcl)
    env["PATH"] = str(MPIEXEC.parent) + os.pathsep + env.get("PATH", "")
    cmd = [str(MPIEXEC)] + (["-localonly"] if os.name == "nt" else []) + \
          ["-n", str(np), str(MP_EXE), str(DECK)]
    p = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True,
                       encoding="utf-8", errors="replace", timeout=900, env=env, cwd=str(tmp))
    log = p.stdout + p.stderr
    out = Path(tmp) / name
    res = {"log": log}
    if (out.with_suffix(".edp")).exists():
        edp = out.with_suffix(".edp").read_text()
        res["edp"] = edp
        res["fail"] = int(re.search(r"fail (\d+)", edp).group(1))
        res["steps"] = int(re.search(r"steps (\d+)", edp).group(1))
        res["roof"] = [float(x) for x in out.with_suffix(".roof").read_text().split()]
    if (out.with_suffix(".iter")).exists():
        res["iters"] = int(out.with_suffix(".iter").read_text().split()[1])
    return res


@pytest.fixture(scope="module")
def single_available(tmp_path_factory):
    if not (MP_EXE.exists() and MPIEXEC.exists()):
        pytest.skip("dist OpenSeesMP / mpiexec not built")
    r = _run(tmp_path_factory.mktemp("probe"), "probe", "-precision single", model="elastic", nstep=2)
    return REFUSED not in r["log"]


@needs_mp
def test_default_off_matches_explicit_double(tmp_path):
    # Same code path (precision 0 never touches the WP-194 code), but a 2-rank
    # MUMPS run is not byte-reproducible run to run (about 1 run in 6-12 differs in
    # the last bits, same binary, MKL_CBWR=AUTO,STRICT included; quirk
    # WP-194-mumps-two-rank-runs-are-not-byte-reproducible). So: same iterations,
    # same history to 1e-12 relative. The byte-identity check against the
    # pre-WP-194 binary is in ADR-194 section 7 (distributions, not n=1).
    a = _run(tmp_path, "none", "")
    b = _run(tmp_path, "dbl", "-precision double")
    assert a["fail"] == 0 and b["fail"] == 0
    assert a["iters"] == b["iters"]
    peak = max(abs(x) for x in a["roof"])
    assert max(abs(x - y) for x, y in zip(a["roof"], b["roof"])) <= 1e-12 * peak
    assert "SMUMPS" not in a["log"] and "SMUMPS" not in b["log"]


@needs_mp
@pytest.mark.parametrize("args", ["-precision quad", "-refine 2", "-precision double -refine 1",
                                  "-refine -1 -precision single"])
def test_bad_options_are_refused(tmp_path, args):
    r = _run(tmp_path, "bad", args, model="elastic", nstep=2)
    assert "edp" not in r, "the deck ran past `system` with " + args
    assert "Mumps ERROR" in r["log"]


@needs_mp
def test_missing_smumps_refuses(tmp_path, single_available):
    if single_available:
        pytest.skip("this build has single-precision MUMPS")
    r = _run(tmp_path, "nosingle", "-precision single", model="elastic", nstep=2)
    assert REFUSED in r["log"]
    assert "edp" not in r


@needs_mp
@pytest.mark.parametrize("alias", ["-precision single", "-mixed"])
def test_single_converges_like_double(tmp_path, single_available, alias):
    if not single_available:
        pytest.skip("no single-precision MUMPS in this build")
    d = _run(tmp_path, "d", "")
    s = _run(tmp_path, "s", alias)
    assert d["fail"] == 0 and s["fail"] == 0
    assert s["steps"] == d["steps"] == NSTEP
    peak = max(abs(x) for x in d["roof"])
    diff = max(abs(x - y) for x, y in zip(s["roof"], d["roof"]))
    print(f"\nWP-194 iterations double {d['iters']} single {s['iters']} "
          f"(+{100.0 * (s['iters'] - d['iters']) / d['iters']:.1f} %), "
          f"max |roof_s - roof_d| = {diff:.3e} (peak {peak:.3e})")
    # NormDispIncr 1e-8 stops each step somewhere inside the same tolerance ball
    assert diff <= 1e-5 * max(peak, 1e-3)
    assert s["iters"] <= 1.5 * d["iters"] + 5


@needs_mp
def test_refine_reduces_double_residual_and_stats_report_precision(tmp_path, single_available):
    if not single_available:
        pytest.skip("no single-precision MUMPS in this build")
    r = _run(tmp_path, "ref", "-precision single -refine 2 -stats", model="elastic", nstep=20)
    assert r["fail"] == 0
    log = r["log"]
    assert re.search(r"^MUMPS stats: n=\d+ BLR\(ICNTL35\)=0$", log, re.M)
    assert re.search(r"^  precision           = single \(SMUMPS factor, double residual\), refine 2$", log, re.M)
    hist = re.findall(r"^  mixed: \|\|b-Ax\|\|inf/\|\|b\|\|inf  step (\d) = (\S+)$", log, re.M)
    assert hist, "no residual history"
    by_step = {}
    for k, v in hist:
        by_step.setdefault(int(k), []).append(float(v))
    assert set(by_step) == {0, 1, 2}
    for r0, r2 in zip(by_step[0], by_step[2]):
        assert r2 < r0 or r0 < 1e-14


@needs_mp
def test_single_with_blr(tmp_path, single_available):
    if not single_available:
        pytest.skip("no single-precision MUMPS in this build")
    r = _run(tmp_path, "blr", "-precision single -BLR 1e-5 -stats", nstep=100)
    assert r["fail"] == 0 and r["steps"] == 100
    assert re.search(r"^MUMPS stats: n=\d+ BLR\(ICNTL35\)=1 eps\(CNTL7\)=\S+$", r["log"], re.M)


@needs_mp
def test_single_with_factor_reuse(tmp_path, single_available):
    if not single_available:
        pytest.skip("no single-precision MUMPS in this build")
    d = _run(tmp_path, "dr", "", reuse=10)
    s = _run(tmp_path, "sr", "-precision single -stats", reuse=10)
    assert d["fail"] == 0 and s["fail"] == 0 and s["steps"] == NSTEP
    # the kept factor is the float one: fewer factorization blocks than analyze steps
    nfact = len(re.findall(r"^MUMPS stats:", s["log"], re.M))
    assert 0 < nfact < NSTEP
    peak = max(abs(x) for x in d["roof"])
    diff = max(abs(x - y) for x, y in zip(s["roof"], d["roof"]))
    print(f"\nWP-194 reuse=10: iterations double {d['iters']} single {s['iters']}, "
          f"factorizations (single) {nfact}, max |roof_s - roof_d| = {diff:.3e}")
    assert diff <= 1e-5 * max(peak, 1e-3)
