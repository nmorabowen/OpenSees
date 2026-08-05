"""`system Pardiso` (ADR-75 P1/P1b/P1d/P1e/P1f) — CI regression battery.

Until now every PARDISO correctness claim lived in offline bench scripts
(`Ladruno_files/testbed/perf/phase1/`) that no CI lane runs; a regression in
the solver's state machine — the `factored` gate, the `factorsCurrent` CGS
staleness guard, the half-storage CSR build — would only be caught by someone
re-running the benches by hand. This battery promotes the P1d smoke's checks
into pytest on the fork's own kernels (LadrunoBrick + LadrunoJ2, plasticity
engaged so the tangent actually changes between iterations).

  test_pardiso_matches_umfpack        unsym PARDISO == UmfPack on a plastic run
  test_symmetric_modes_match          -matrixType 1/2 exactness vs full CSR (P1d)
  test_modified_newton_reuse          the phase-22-skip path under an unchanged A
  test_krylov_matches_direct          -krylov CGS answers == direct (P1e)
  test_krylov_refused_on_indefinite   -krylov + -matrixType 2 warns, runs direct
  test_stats_prints_once_per_pattern  -stats fires once, labels the mtype (P1d)
  test_profiler_brackets_present      the phase-split profiler brackets report
  test_bad_options_degrade_not_null   parse failure degrades, never ProfileSPD

NOT covered here (documented gap, not an oversight): the addA asymmetry
detector needs a genuinely unsymmetric element tangent inside the first three
assemblies (contact / non-associated flow / LadrunoUP) — a heavier model than
a CI gate should carry; it stays exercised by the perf decks.

Tolerances are chosen to hold at ANY MKL thread count: threaded PARDISO is not
byte-reproducible run-to-run (P1f, ~1 ULP), so nothing here asserts bitwise
equality; converged-Newton parity at rtol 1e-9 is far above that noise floor
and far below any real solver defect.
"""
import os
import sys

import pytest

# Best-effort determinism: MKL reads this at its first init, so a battery that
# already ran an MKL kernel keeps its count — which is why the tolerances above
# do not depend on it. Must be set before the module is first imported.
os.environ.setdefault("MKL_NUM_THREADS", "1")

from _testbed import ops  # noqa: E402

# PARDISO IS MKL: linked on the Windows/oneAPI build only (Zone-A Ubuntu CI
# uses reference LAPACK and never compiles SRC/system_of_eqn/linearSOE/pardiso,
# same gating as FEAST). Enforced on the Windows box.
pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(sys.platform != "win32",
                       reason="PARDISO requires MKL (Windows/oneAPI build)"),
]

# ---- the shared model: LadrunoBrick + LadrunoJ2 cantilever ------------------
# Shrunk from phase1/p1d_smoke.py (nx=8 -> 4; ~300 free DOF) so the whole
# battery costs seconds. The lateral top load drives the base past yield
# (sigma_base ~ 1.5*s0 at full load), so the tangent changes between Newton
# iterations — without that, the factored-gate and -krylov paths under test
# would never actually be exercised.
E, NU = 200000.0, 0.3
KMOD = E / (3.0 * (1.0 - 2.0 * NU))
GMOD = E / (2.0 * (1.0 + NU))
S0, QINF, BEXP, HISO = 250.0, 0.0, 0.0, 2000.0
NX = 4
LEL = 100.0
NSTEPS = 5
TOL = 1.0e-8


def _nid(i, j, k):
    return 1 + i + (NX + 1) * (j + (NX + 1) * k)


def _run(system_args, algorithm="Newton"):
    """Build the cantilever, run NSTEPS load steps, return the tip (ux, uz)."""
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("LadrunoJ2", 1, KMOD, GMOD, "-iso", "voce",
                   S0, QINF, BEXP, HISO)
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
                    "LadrunoBrick", tag,
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
            ops.load(_nid(i, j, NX), 4.0e5, 0.0, -4.0e5)
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system(*system_args)
    ops.test("NormDispIncr", TOL, 30)
    ops.algorithm(algorithm)
    ops.integrator("LoadControl", 1.0 / NSTEPS)
    ops.analysis("Static")
    for s in range(NSTEPS):
        assert ops.analyze(1) == 0, (
            f"analyze failed at step {s + 1} with system {system_args}")
    tip = _nid(NX // 2, NX // 2, NX)
    return ops.nodeDisp(tip, 1), ops.nodeDisp(tip, 3)


def _drain(capfd):
    """Collect this test's console output; skip if the build lacks PARDISO.

    A Windows build configured without MKL still parses `system Pardiso` —
    the factory warns "this build has no MKL" and returns 0, and OpenSees
    quietly falls back to ProfileSPD, which would make every parity test here
    pass vacuously. Turn that into an explicit skip instead.
    """
    cap = capfd.readouterr()
    text = cap.out + cap.err
    if "has no MKL" in text:
        pytest.skip("this build was configured without MKL/PARDISO")
    return text


def _rel(a, b):
    return abs(a - b) / max(abs(b), 1e-30)


# ---- tests ------------------------------------------------------------------

def test_pardiso_matches_umfpack(capfd):
    """P1/P1b: the unsymmetric desktop path reproduces UmfPack through a
    5-step plastic run (both fully converged Newton at the same tolerance)."""
    ref = _run(["UmfPack"])
    got = _run(["Pardiso"])
    _drain(capfd)
    assert _rel(got[0], ref[0]) < 1e-9 and _rel(got[1], ref[1]) < 1e-9, (
        f"Pardiso {got} vs UmfPack {ref}")


def test_symmetric_modes_match(capfd):
    """P1d exactness: -matrixType 2 (LDL^T Bunch-Kaufman) and 1 (SPD) are
    EXACT, not approximate — same answer as the full unsymmetric CSR. The J2
    cantilever's consistent tangent is symmetric (associative flow) and stays
    positive definite (hardening), so both modes are legal here."""
    ref = _run(["Pardiso"])
    sym = _run(["Pardiso", "-matrixType", 2])
    spd = _run(["Pardiso", "-matrixType", 1])
    _drain(capfd)
    for name, got in (("-matrixType 2", sym), ("-matrixType 1", spd)):
        assert _rel(got[0], ref[0]) < 1e-9 and _rel(got[1], ref[1]) < 1e-9, (
            f"{name} {got} vs unsymmetric {ref}")


def test_modified_newton_reuse(capfd):
    """P1: ModifiedNewton holds A fixed across iterations, so all but the
    first solve of each step must take the phase-33 shortcut (the `factored`
    gate). A defect there answers with a stale factorization — which this
    parity against UmfPack under the SAME algorithm would catch."""
    ref = _run(["UmfPack"], algorithm="ModifiedNewton")
    got = _run(["Pardiso"], algorithm="ModifiedNewton")
    _drain(capfd)
    assert _rel(got[0], ref[0]) < 1e-9 and _rel(got[1], ref[1]) < 1e-9, (
        f"Pardiso {got} vs UmfPack {ref} (ModifiedNewton)")


def test_krylov_matches_direct(capfd):
    """P1e: -krylov reuses the retained L/U as a CGS preconditioner across a
    CHANGED tangent. A win leaves the stored factors one tangent behind — the
    `factorsCurrent` guard exists so a later solve cannot silently answer with
    the previous tangent. Parity is tolerance-limited (eps_CGS = 1e-6 feeds
    Newton, which then converges to TOL), hence the looser 1e-6 here; a
    staleness bug is orders of magnitude above that."""
    ref = _run(["Pardiso"])
    got = _run(["Pardiso", "-krylov", 6])
    text = _drain(capfd)
    assert "CGS solve (phase 23)" not in text, (
        "phase-23 hard error during the -krylov run")
    assert _rel(got[0], ref[0]) < 1e-6 and _rel(got[1], ref[1]) < 1e-6, (
        f"-krylov 6 {got} vs direct {ref}")


def test_krylov_refused_on_indefinite(capfd):
    """P1e: Intel documents no CGS mode for mtype -2, so -krylov with
    -matrixType 2 must WARN and continue with direct factorization — making
    the answer exactly the plain -matrixType 2 one."""
    ref = _run(["Pardiso", "-matrixType", 2])
    got = _run(["Pardiso", "-matrixType", 2, "-krylov", 6])
    text = _drain(capfd)
    assert "-krylov is not available" in text, (
        "expected the mtype -2 refusal warning")
    assert _rel(got[0], ref[0]) < 1e-9 and _rel(got[1], ref[1]) < 1e-9, (
        f"refused -krylov {got} vs plain -matrixType 2 {ref}")


def test_stats_prints_once_per_pattern(capfd):
    """P1d: -stats reports the peak-memory counters ONCE per sparsity pattern
    (not per Newton iteration — pattern-determined counters reprinted every
    solve would be noise), labelled with the storage actually in use."""
    _run(["Pardiso", "-stats"])
    text = _drain(capfd)
    assert text.count("PARDISO stats:") == 1, text
    assert "TOTAL PEAK" in text
    assert "unsymmetric, full CSR" in text

    _run(["Pardiso", "-matrixType", 2, "-stats"])
    text = capfd.readouterr()
    text = text.out + text.err
    assert text.count("PARDISO stats:") == 1, text
    assert "symmetric, upper-triangle CSR" in text


def test_profiler_brackets_present(tmp_path, capfd):
    """ADR-75 profiling: PARDISO reports its phase split under the fork's
    stack profiler with the same bracket names as UmfPack (soe.symbolic /
    soe.factor / soe.trisolve), plus its own soe.cgs (the -krylov phase-23
    regime), the dc.s.fill/dc.s.verify setSize sub-brackets, and the
    deep-gated soe.addA scatter bracket.

    ONE deep -krylov run touches all of them: symbolic once per pattern,
    factor+trisolve on the first solve (no factors retained yet), cgs on the
    later changed tangents, addA on every assembly. The rollup stores each
    bracket as an HDF5 GROUP NAME, so presence is verifiable on the raw bytes
    — no h5py dependency in a zone_a test.
    """
    ops.wipe()
    ops.profiler("start", "-deep")
    _run(["Pardiso", "-krylov", 6])
    ops.profiler("stop")
    h5 = tmp_path / "pardiso_prof.h5"
    ops.profiler("report", str(h5), "-run", "pardiso")
    ops.profiler("reset")
    _drain(capfd)

    data = h5.read_bytes()
    assert data[:4] == b"\x89HDF", "profiler report is not an HDF5 file"
    for name in (b"soe.symbolic", b"soe.factor", b"soe.trisolve", b"soe.cgs",
                 b"dc.s.fill", b"dc.s.verify", b"soe.addA"):
        assert name in data, (
            f"bracket {name.decode()} missing from the profiler report")

    # Negative control, double duty: (a) soe.addA is DEEP-gated — a coarse run
    # must not pay (or report) the per-element scatter bracket; (b) proves the
    # byte-scan above is not vacuous — a bracket that does not fire really is
    # absent from the report.
    ops.profiler("start")
    _run(["Pardiso"])
    ops.profiler("stop")
    h5c = tmp_path / "pardiso_prof_coarse.h5"
    ops.profiler("report", str(h5c), "-run", "coarse")
    ops.profiler("reset")
    _drain(capfd)

    data = h5c.read_bytes()
    assert b"soe.factor" in data
    assert b"soe.addA" not in data, (
        "deep-gated soe.addA leaked into a coarse profiler run")
    assert b"soe.cgs" not in data, (
        "soe.cgs reported on a run without -krylov")


def test_bad_options_degrade_not_null(capfd):
    """P1d's hard-won rule: a parse failure must DEGRADE, never return a null
    SOE — a null silently selects ProfileSPD (catastrophically slow on a 3D
    mesh, cost a 25-minute dark bench). This is the exact historical trigger:
    -matrixType passed as a STRING. The unknown-option warning is pinned in
    the same run."""
    ref = _run(["Pardiso"])
    got = _run(["Pardiso", "-matrixType", "2", "-bogus"])
    text = _drain(capfd)
    assert "failed to get -matrixType" in text, (
        "expected the string-value degrade warning")
    assert "unknown option -bogus" in text
    assert "ProfileSPD" not in text, (
        "the factory returned null and OpenSees fell back to ProfileSPD")
    assert _rel(got[0], ref[0]) < 1e-9 and _rel(got[1], ref[1]) < 1e-9, (
        f"degraded run {got} vs unsymmetric {ref}")
