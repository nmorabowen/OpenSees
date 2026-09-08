"""`system Pardiso -stats` / `-pardisoStats` (ADR-75 P1k) — the MUMPS-shaped
factorization-stats block.

TIMs PM-01 D26 wants factorization memory and fill in every desktop leg's log.
The serial `MumpsSolver` this fork ships is never compiled in (see
Ladruno_implementation/LEDGER_vanilla_files.md and the ADR-75 P5 docs), so MKL
PARDISO (`system Pardiso`) is the only desktop route that can produce that
number, and it is now printed in the SAME labelled-block shape as the shipped
MUMPS `-stats` (MumpsParallelSolver.cpp:295-319): a one-line header followed by
indented "label iparm(N)  = value" lines, checkable directly against the MKL
Developer Guide's PARDISO iparm table.

  test_stats_block_present_with_flag       -stats prints every required label
  test_pardiso_stats_alias_matches_stats   -pardisoStats is a synonym for -stats
  test_stats_block_absent_without_flag     no flag => no "PARDISO stats:" line
  test_stats_prints_every_factorization    NOT latched to "once per pattern" --
                                            a second refactorization of the SAME
                                            sparsity pattern logs its own block

Model: a single small elastic 3-D brick cantilever (2x2x2 stdBrick elements)
solved in one LoadControl step -- correctness of the mechanics is not the
point here (test_pardiso_solver.py already covers PARDISO answer-parity on a
real plastic kernel); this file only exercises the -stats REPORTING path.
"""
import os
import sys

import pytest

os.environ.setdefault("MKL_NUM_THREADS", "1")

from _testbed import ops  # noqa: E402

# Same gating as test_pardiso_solver.py: PARDISO is MKL-only, Windows/oneAPI
# build; Zone-A's Ubuntu reference-LAPACK runner never compiles it.
pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="PARDISO requires MKL (Windows/oneAPI build)"),
]

NX = 2
LEL = 100.0
E, NU = 200000.0, 0.3


def _nid(i, j, k):
    return 1 + i + (NX + 1) * (j + (NX + 1) * k)


def _build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    for k in range(NX + 1):
        for j in range(NX + 1):
            for i in range(NX + 1):
                ops.node(_nid(i, j, k), i * LEL, j * LEL, k * LEL)
    for j in range(NX + 1):
        for i in range(NX + 1):
            ops.fix(_nid(i, j, 0), 1, 1, 1)
    tag = 1
    for k in range(NX):
        for j in range(NX):
            for i in range(NX):
                ops.element(
                    "stdBrick", tag,
                    _nid(i, j, k), _nid(i + 1, j, k),
                    _nid(i + 1, j + 1, k), _nid(i, j + 1, k),
                    _nid(i, j, k + 1), _nid(i + 1, j, k + 1),
                    _nid(i + 1, j + 1, k + 1), _nid(i, j + 1, k + 1),
                    1)
                tag += 1
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(NX + 1):
        for i in range(NX + 1):
            ops.load(_nid(i, j, NX), 1.0e4, 0.0, -1.0e4)
    ops.constraints("Plain")
    ops.numberer("RCM")


def _run(system_args, nsteps=1):
    """Build the cantilever and run `nsteps` LoadControl steps under
    `system_args`. Each step re-assembles A and (with a changed tangent, or
    just re-solved) triggers PARDISO's `factored == false` path, i.e. a real
    phase-22 refactorization -- an ELASTIC model still refactorizes step to
    step because the SOE is rebuilt/re-zeroed and reloaded each analyze()."""
    _build()
    ops.system(*system_args)
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / nsteps)
    ops.analysis("Static")
    for s in range(nsteps):
        assert ops.analyze(1) == 0, (
            f"analyze failed at step {s + 1} with system {system_args}")


def _drain(capfd):
    """Collect console output; skip cleanly if the build has no MKL/PARDISO
    (same rationale as test_pardiso_solver.py's `_drain`: a build without MKL
    still parses `system Pardiso` and silently falls back to ProfileSPD, which
    would make every assertion here pass or fail vacuously)."""
    cap = capfd.readouterr()
    text = cap.out + cap.err
    if "has no MKL" in text:
        pytest.skip("this build was configured without MKL/PARDISO")
    return text


REQUIRED_LABELS = (
    "factor entries iparm(18)",
    "peak memory KB iparm(15)",
    "perm memory KB iparm(16)",
    "fact memory KB iparm(17)",
    "factor Mflops  iparm(19)",
)


def _assert_block(text):
    assert "PARDISO stats:" in text, text
    assert "matrixType=" in text
    assert "threads=" in text
    for label in REQUIRED_LABELS:
        assert label in text, f"missing '{label}' in:\n{text}"


def _factor_entries(text):
    """Pull the integer value off the 'factor entries iparm(18)  = N' line."""
    for line in text.splitlines():
        if "factor entries iparm(18)" in line:
            return int(line.split("=")[-1].strip())
    raise AssertionError(f"no 'factor entries iparm(18)' line in:\n{text}")


# ---- tests ------------------------------------------------------------------

def test_stats_block_present_with_flag(capfd):
    """`-stats` prints the full MUMPS-shaped block, with every required label,
    and the reported factor-entry count (nnz in L+U) is a positive integer."""
    _run(["Pardiso", "-stats"])
    text = _drain(capfd)
    _assert_block(text)
    n_entries = _factor_entries(text)
    assert isinstance(n_entries, int) and n_entries > 0, (
        f"factor entries iparm(18) should be a positive integer, got {n_entries}")


def test_pardiso_stats_alias_matches_stats(capfd):
    """`-pardisoStats` is a bare-flag synonym for `-stats` (mirrors the
    `-stats`/`-mumpsStats` pair `system Mumps` already offers) -- same block,
    same labels."""
    _run(["Pardiso", "-pardisoStats"])
    text = _drain(capfd)
    _assert_block(text)


def test_stats_block_absent_without_flag(capfd):
    """No `-stats`/`-pardisoStats` => no "PARDISO stats:" line at all -- the
    feature must be strictly opt-in (it costs extra MKL analysis time)."""
    _run(["Pardiso"])
    text = _drain(capfd)
    assert "PARDISO stats:" not in text, text


def test_stats_prints_every_factorization(capfd):
    """P1k (superseding P1d's "once per sparsity pattern" latch): a SECOND
    refactorization of the identical sparsity pattern -- forced here by two
    successive analyze() calls under LoadControl, both of which reassemble A
    and hit `factored == false` -- must log its OWN block, not be swallowed
    by a per-pattern latch. TIMs PM-01 D26 wants every desktop factorization
    event in the log, matching how the shipped MUMPS `-stats` prints inside
    `if (factored == false)` on every job=5 call."""
    _run(["Pardiso", "-stats"], nsteps=2)
    text = _drain(capfd)
    assert text.count("PARDISO stats:") >= 2, (
        f"expected at least 2 stats blocks (one per refactorization), got:\n{text}")
