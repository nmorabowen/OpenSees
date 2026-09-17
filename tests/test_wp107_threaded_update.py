"""WP-107 -- the threaded `Domain::update()` element loop (ADR-75b stage L3-1).

This file is the WARRANT for WP-107.  It exists because the red-team review of
PR #843 found the work package shipped with **no automated test at all** (B1):
the PR's "128 passed" was the pre-existing suite run once by hand at the default
1 thread, i.e. it exercised the serial fallback and nothing else.  Three
mutations would have survived it untouched:

  M1  delete the allowlist loop (`SRC/domain/domain/Domain.cpp`, the
      `ladrunoThreadSafeUpdate() == false` sweep),
  M2  flip `Element::ladrunoThreadSafeUpdate()`'s default from false to true
      (`SRC/element/Element.h`),
  M3  revert `LadrunoQuad::shp` / `::shpBar` from `thread_local` back to plain
      `static` (`SRC/element/ladrunoPlane/LadrunoQuad.h`).

WHICH TEST CATCHES WHICH -- the whole point of the file, so state it plainly:

  M1  caught by `test_refusal_*` and `test_reaudit_*`: with the sweep gone the
      un-audited element is never noticed, so no `WARNING ladrunoThreads:
      element <tag> ... not on the ... allowlist` line is emitted and the
      `assert ... in text` fails.  (They also assert the run was NOT announced
      as THREADED, so a mutation that keeps the message but threads anyway
      still fails.)
  M2  caught by the same two: every element then answers true, the refusal
      never fires, and the deck that must run serial runs threaded.
  M3  caught by `test_bit_identity_*`: `shp`/`shpBar` are one shared buffer
      written by `shapeFunction()`/`computeShapeBar()` and read by `formB()`,
      so on a shared `static` two threads interleave mid-element and the
      strains -- hence every nodal displacement and reaction -- stop matching
      the 1-thread field.  The deck is sized (800 `LadrunoQuad`, `-bbar`, 4
      threads, dynamic schedule) so the collision is not a lucky miss.

THE ORACLE IS THE FULL FIELD, per red-team S6.  The shipped bench harness
hashed two numbers per step, one of them (`settlement_m`) a PRESCRIBED `sp`
value that is identical by construction at every thread count, and the other a
SUM of reactions in which a per-element scramble can cancel.  Here the hash
covers every node's `nodeDisp` and `nodeReaction` at `repr()` precision plus one
element's stress vector, recorded **at every step** rather than only at the end
-- a transient that diverges and re-converges is then still caught.

Wall time: ~25 s (8 identity runs + 5 refusal decks + 4 child processes).

See `Ladruno_implementation/107_ladruno_openmp_element_loop.md`.
"""
import hashlib
import os
import sys

import pytest

# Same rule as ADR-95's probes: every test runs the binary built in THIS
# worktree.  `_testbed` does a bare `import opensees`, which an installed
# Ladruno .pth can answer; putting dist/bin first makes the local build win
# unaided (memory entry "Installed Ladruno .pth hijacks import opensees").
_DIST_BIN = os.environ.get("LADRUNO_DIST_BIN") or os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "dist", "bin")
if os.path.isdir(_DIST_BIN) and _DIST_BIN not in sys.path:
    sys.path.insert(0, _DIST_BIN)

from _testbed import ops                        # noqa: E402
from _testbed.subprocess_run import run_python_script   # noqa: E402


# ---------------------------------------------------------------------------
#  IS THE FEATURE COMPILED INTO THIS BINARY AT ALL?   (added by the PR #843
#  CI fix, 2026-09-16)
# ---------------------------------------------------------------------------
#  `LADRUNO_OPENMP` is a BUILD option.  In a binary built without it,
#  `ladrunoThreads(n)` stores n, warns, and the element loop stays serial --
#  so every threading assertion below fails for a reason that has nothing to
#  do with the code under test.  That is not hypothetical: it is exactly how
#  this file went red on Zone-A (run 35160539366, 10 failed), because CI
#  configures with a bare `cmake` and the option then defaulted OFF.
#
#  The option now defaults ON in CMakeLists.txt, so CI RUNS these tests.  This
#  probe exists for the other direction: a developer who deliberately builds
#  with `-DLADRUNO_OPENMP=OFF` (or `set LADRUNO_NO_OPENMP=1`) gets a clean
#  skip instead of 10 misleading failures.
#
#  THE PROBE ASKS THE BINARY, and it is biased to RUN rather than to skip.
#  There is no `ladrunoOpenMP` query verb -- the honest signal the binary
#  already emits is the `ladrunoThreads` warning itself, so that is what is
#  read.  Every ambiguous outcome (child crashed, no marker, a box with a
#  single hardware thread where the warning cannot fire) resolves to "run the
#  tests", because a file that silently skips on CI is worth less than one
#  that fails loudly: the skip must never become the way this WP stops being
#  gated.  A child process is used so the probe cannot perturb the in-process
#  thread count the tests below depend on.
_OPENMP_PROBE = (
    "import sys; sys.path.insert(0, %r)\n"
    "import opensees as ops\n"
    "print('LADRUNO_PROBE_HW', ops.ladrunoThreads(1 << 20))\n"
)
_OPENMP_ABSENT_MARKER = "built WITHOUT LADRUNO_OPENMP"


def _probe_openmp_compiled_in():
    """-> (compiled_in, why). `why` is the skip reason when compiled_in is False."""
    try:
        code, out = run_python_script(_OPENMP_PROBE % _DIST_BIN)
    except Exception:                      # the probe must never break collection
        return True, ""
    if code != 0 or "LADRUNO_PROBE_HW" not in out:
        return True, ""                    # ambiguous -> run, and fail loudly
    if _OPENMP_ABSENT_MARKER in out:
        return False, (
            "this binary was built with LADRUNO_OPENMP=OFF, so `ladrunoThreads`"
            " cannot thread anything and WP-107 has no threaded path to test."
            " Rebuild with -DLADRUNO_OPENMP=ON -- it is the CMakeLists default"
            " since the PR #843 CI fix, and `Ladruno_scripts\\build.bat` passes"
            " it explicitly. (Zone-A builds with it ON; these tests RUN there.)")
    return True, ""


_OPENMP_IN, _OPENMP_SKIP_WHY = _probe_openmp_compiled_in()

pytestmark = [
    pytest.mark.zone_a,
    pytest.mark.skipif(not _OPENMP_IN, reason=_OPENMP_SKIP_WHY),
]

E, NU, RHO = 30000.0, 0.25, 2.4
H = 0.25


def _threads_available():
    """The counts worth testing on THIS box.

    `ladruno_setNumThreads` clamps to hardware concurrency (red-team S7), so on
    a 2-core runner asking for 4 would silently become 2 and the "different
    thread counts agree" claim would be vacuous.  Ask the binary what it will
    actually honour.
    """
    hw = ops.ladrunoThreads(1024)     # clamped to hardware concurrency
    ops.ladrunoThreads(1)
    return [n for n in (1, 2, 4) if n <= hw]


# ---------------------------------------------------------------------------
#  the deck
# ---------------------------------------------------------------------------
def _build(nx, ny, form="bbar", mix="none", mode="static"):
    """A rectangular `LadrunoQuad` strip, self-weight + a skew edge load.

    `mix` injects exactly ONE un-audited object at element tag 7:
      "wrapper"    -- an `nDMaterial PlaneStrain` wrapper (material-level
                      refusal: the element is allowlisted, its material is not),
      "eas"        -- `-formulation eas` (element-level: formEAStrue runs on
                      ~12 shared function-scope statics and condenses through
                      Matrix inversion).
    "zerolength" adds an un-audited ELEMENT beside the allowlisted ones.
    """
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    nid = {}
    t = 1
    for i in range(nx + 1):
        for j in range(ny + 1):
            nid[(i, j)] = t
            ops.node(t, i * H, j * H)
            t += 1
    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    if mix == "wrapper":
        ops.nDMaterial("ElasticIsotropic", 2, E, NU, RHO)
        ops.nDMaterial("PlaneStrain", 3, 2)

    eid = 1
    for i in range(nx):
        for j in range(ny):
            mt, f = 1, form
            if mix == "wrapper" and eid == 7:
                mt = 3
            if mix == "eas" and eid == 7:
                f = "eas"
            ops.element("LadrunoQuad", eid, nid[(i, j)], nid[(i + 1, j)],
                        nid[(i + 1, j + 1)], nid[(i, j + 1)], mt,
                        "-thick", 1.0, "-type", "PlaneStrain",
                        "-formulation", f, "-rho", RHO, "-body", 0.0, -9.81)
            eid += 1

    for i in range(nx + 1):
        ops.fix(nid[(i, 0)], 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for i in range(nx + 1):
        ops.load(nid[(i, ny)], 3.7, -11.3)

    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-10, 50, 0)
    ops.algorithm("Newton")
    if mode == "static":
        ops.integrator("LoadControl", 0.1)
        ops.analysis("Static")
    else:
        ops.integrator("Newmark", 0.5, 0.25)
        ops.analysis("Transient")
    return nid


_ZL_FIXED, _ZL_FREE = 999998, 999999
_ZL_TAG = 900000


def _add_unaudited_zerolength():
    """A `zeroLength` -- un-audited, and NOT a LadrunoQuad, so it exercises the
    element half of the allowlist rather than the material half.

    BOTH nodes are fully restrained on purpose: the re-audit test removes this
    element again mid-run, and a node left with a free DOF and nothing attached
    to it makes the next factorization singular -- which would look like a
    WP-107 failure and is only a defect in the probe.
    """
    ops.uniaxialMaterial("Elastic", 99, 1.0e6)
    ops.node(_ZL_FIXED, 0.0, -1.0)
    ops.node(_ZL_FREE, 0.0, -1.0)
    ops.fix(_ZL_FIXED, 1, 1)
    ops.fix(_ZL_FREE, 1, 1)
    ops.element("zeroLength", _ZL_TAG, _ZL_FIXED, _ZL_FREE, "-mat", 99,
                "-dir", 1)


def _run_and_hash(nid, steps, mode="static", probe_ele=1):
    """md5 over the FULL field at EVERY step (red-team S6).

    Every node's displacement and reaction plus one element's stress vector,
    all at `repr()` precision -- i.e. every bit of the double.
    """
    h = hashlib.md5()
    rcs = []
    for _ in range(steps):
        rc = ops.analyze(1) if mode == "static" else ops.analyze(1, 0.02)
        rcs.append(rc)
        ops.reactions()
        for k in sorted(nid.values()):
            h.update(("%d|%s|%s" % (k, repr(ops.nodeDisp(k)),
                                    repr(ops.nodeReaction(k)))).encode())
        sig = ops.eleResponse(probe_ele, "stress")
        assert sig, (
            "eleResponse(%d,'stress') came back empty -- the per-step element "
            "term of the oracle would be a constant, which is exactly the "
            "defect red-team S6 raised about the shipped bench harness"
            % probe_ele)
        h.update(repr(list(sig)).encode())
        if rc != 0:
            break
    return h.hexdigest(), rcs


# ---------------------------------------------------------------------------
#  (a) full-field bit-identity at 1 / 2 / 4 threads
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mode", ["static", "transient"])
@pytest.mark.parametrize("form", ["std", "bbar"])
def test_bit_identity_across_thread_counts(mode, form, capfd):
    """MUTATION M3 (`shp`/`shpBar` back to `static`) dies here.

    Bit-identity is a CONSTRUCTION claim, not a tolerance one: `Domain::update`
    carries no floating-point reduction (ADR-75b section 2.1), so the threaded
    loop performs the identical arithmetic in the identical order and the md5
    must match exactly.  Any thread-count-dependent digit is a defect.
    """
    counts = _threads_available()
    assert 1 in counts
    if len(counts) < 2:
        pytest.skip("box reports a single hardware thread")

    digests, announced = {}, {}
    for n in counts:
        ops.ladrunoThreads(n)
        nid = _build(40, 20, form=form, mode=mode)
        digests[n], rcs = _run_and_hash(nid, 6, mode=mode)
        assert all(rc == 0 for rc in rcs), (n, rcs)
        announced[n] = capfd.readouterr()
    ops.ladrunoThreads(1)

    ref = digests[counts[0]]
    assert all(d == ref for d in digests.values()), digests

    # ... and the run really was threaded, so the identity above is evidence
    # about the threaded path and not about four serial runs.  (Red-team S3
    # un-latched this announcement: before the fix only the first model in the
    # process said anything, so this assertion could not have been written.)
    for n in counts:
        text = announced[n].err + announced[n].out
        if n == 1:
            assert "THREADED" not in text, text
        else:
            assert ("element update loop THREADED on %d threads" % n) in text, text


# ---------------------------------------------------------------------------
#  (b) refusal -- an un-audited object anywhere makes the WHOLE loop serial
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mix,expect_tag", [
    ("wrapper", 7),          # allowlisted element, un-audited MATERIAL
    ("eas", 7),              # un-audited FORMULATION of an allowlisted class
])
def test_refusal_names_the_tag_and_is_byte_identical_to_serial(mix, expect_tag,
                                                               capfd):
    """MUTATIONS M1 and M2 die here.

    All-or-nothing on purpose: a mixed parallel/serial split would still let an
    audited element run concurrently with an un-audited one that writes shared
    node state.
    """
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("box reports a single hardware thread")
    n = counts[-1]

    ops.ladrunoThreads(1)
    nid = _build(20, 10, mix=mix)
    serial, rc_s = _run_and_hash(nid, 4)
    capfd.readouterr()

    ops.ladrunoThreads(n)
    nid = _build(20, 10, mix=mix)
    threaded, rc_t = _run_and_hash(nid, 4)
    text = "".join(capfd.readouterr())
    ops.ladrunoThreads(1)

    assert "not on the" in text and "allowlist" in text, text
    assert ("element %d " % expect_tag) in text, text
    assert "THREADED" not in text, text
    assert threaded == serial, (mix, serial, threaded)
    assert rc_t == rc_s


def test_refusal_on_an_unaudited_element_class(capfd):
    """Same, but the offender is an un-audited ELEMENT (`zeroLength`) rather
    than an un-audited material -- the two halves of the allowlist."""
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("box reports a single hardware thread")
    n = counts[-1]

    ops.ladrunoThreads(1)
    nid = _build(20, 10)
    _add_unaudited_zerolength()
    serial, _ = _run_and_hash(nid, 4)
    capfd.readouterr()

    ops.ladrunoThreads(n)
    nid = _build(20, 10)
    _add_unaudited_zerolength()
    threaded, _ = _run_and_hash(nid, 4)
    text = "".join(capfd.readouterr())
    ops.ladrunoThreads(1)

    assert ("element %d " % _ZL_TAG) in text, text
    assert "allowlist" in text, text
    assert "THREADED" not in text, text
    assert threaded == serial


# ---------------------------------------------------------------------------
#  (c) re-audit -- red-team S3
# ---------------------------------------------------------------------------
def test_reaudit_after_an_unaudited_element_is_added_mid_run(capfd):
    """The audit re-runs every `Domain::update()`; S3 was that the MESSAGE was
    latched behind a process-wide `static bool`, so an element added after the
    first threaded step refused SILENTLY -- a run that quietly went serial and
    a run that stayed threaded looked identical.

    Sequence: threaded step -> add an un-audited element -> the next step must
    refuse ALOUD and name it -> remove it -> the loop must announce THREADED
    again.  Each transition is keyed on the element generation, which
    `addElement` / `removeElement` / `clearAll` bump.
    """
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("box reports a single hardware thread")
    n = counts[-1]

    ops.ladrunoThreads(n)
    nid = _build(20, 10)
    assert ops.analyze(1) == 0
    first = "".join(capfd.readouterr())
    assert ("element update loop THREADED on %d threads" % n) in first, first

    _add_unaudited_zerolength()          # generation++, un-audited
    assert ops.analyze(1) == 0
    after_add = "".join(capfd.readouterr())
    assert ("element %d " % _ZL_TAG) in after_add, after_add
    assert "allowlist" in after_add, after_add

    ops.remove("element", _ZL_TAG)       # generation++, clean again
    assert ops.analyze(1) == 0
    after_rm = "".join(capfd.readouterr())
    assert ("element update loop THREADED on %d threads" % n) in after_rm, after_rm
    ops.ladrunoThreads(1)


def test_wipe_makes_the_announcement_speak_again(capfd):
    """`wipe` is the other S3 case: the second model in a process was mute."""
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("box reports a single hardware thread")
    n = counts[-1]

    ops.ladrunoThreads(n)
    for round_ in range(2):
        nid = _build(20, 10)             # _build() starts with ops.wipe()
        assert ops.analyze(1) == 0
        text = "".join(capfd.readouterr())
        assert ("element update loop THREADED on %d threads" % n) in text, (
            round_, text)
    ops.ladrunoThreads(1)


def test_announcement_does_not_repeat_every_newton_iteration(capfd):
    """The converse of S3: un-latching must not turn the announcement into
    per-iteration spam.  One line per (model, outcome, thread count)."""
    counts = _threads_available()
    if len(counts) < 2:
        pytest.skip("box reports a single hardware thread")
    n = counts[-1]

    ops.ladrunoThreads(n)
    nid = _build(20, 10)
    for _ in range(5):
        assert ops.analyze(1) == 0
    text = "".join(capfd.readouterr())
    ops.ladrunoThreads(1)
    assert text.count("element update loop THREADED") == 1, text


# ---------------------------------------------------------------------------
#  (d) the knob itself -- red-team S7
# ---------------------------------------------------------------------------
def test_default_is_one_thread():
    """ADR-40's standing anti-goal is "OpenMP-by-default"."""
    code, out = run_python_script(
        "import sys; sys.path.insert(0, %r)\n"
        "import opensees as ops\n"
        "print('DEFAULT', ops.ladrunoThreads())\n" % _DIST_BIN)
    assert code == 0, out
    assert "DEFAULT 1" in out, out


@pytest.mark.parametrize("value,why", [
    ("abc", "not an integer"),
    ("2.7", "trailing text"),
    ("4x", "trailing text"),
    ("0", "< 1"),
    ("-3", "< 1"),
])
def test_env_var_rejects_bad_values_loudly(value, why):
    """Before red-team S7 every one of these ran SERIAL with NO MESSAGE AT ALL.

    A typo'd `LADRUNO_THREADS` in a job script is precisely how a bench lies,
    which is the anti-goal `SRC/utility/LadrunoThreads.h` exists to serve.  A
    child process is required: the env var is read once, lazily, per process.
    """
    prev = os.environ.get("LADRUNO_THREADS")
    os.environ["LADRUNO_THREADS"] = value
    try:
        code, out = run_python_script(
            "import sys; sys.path.insert(0, %r)\n"
            "import opensees as ops\n"
            "print('GOT', ops.ladrunoThreads())\n" % _DIST_BIN,
            merge_stderr=True)
    finally:
        if prev is None:
            os.environ.pop("LADRUNO_THREADS", None)
        else:
            os.environ["LADRUNO_THREADS"] = prev

    assert code == 0, out
    assert "GOT 1" in out, out
    assert "LADRUNO_THREADS" in out and "WARNING" in out, out


def test_env_var_above_hardware_concurrency_is_clamped_with_a_warning():
    """`LADRUNO_THREADS=99999` used to be honoured as 1024 threads -- measured
    by the red team at 5x SLOWER than serial, announced as a threaded run."""
    prev = os.environ.get("LADRUNO_THREADS")
    os.environ["LADRUNO_THREADS"] = "99999"
    try:
        code, out = run_python_script(
            "import sys; sys.path.insert(0, %r)\n"
            "import opensees as ops\n"
            "print('GOT', ops.ladrunoThreads())\n" % _DIST_BIN,
            merge_stderr=True)
    finally:
        if prev is None:
            os.environ.pop("LADRUNO_THREADS", None)
        else:
            os.environ["LADRUNO_THREADS"] = prev

    assert code == 0, out
    assert "CLAMPED" in out, out
    got = int(out.split("GOT")[1].split()[0])
    assert 1 <= got <= 1024, out
    assert got < 99999, out


def test_verb_clamps_and_reports(capfd):
    """The verb path clamps through the same helper as the env path, so both
    say the same thing in the same words."""
    assert ops.ladrunoThreads(0) == 1
    t0 = "".join(capfd.readouterr())
    assert "WARNING" in t0 and "< 1" in t0, t0

    hw = ops.ladrunoThreads(99999)
    t1 = "".join(capfd.readouterr())
    assert "CLAMPED" in t1, t1
    assert hw < 99999
    assert ops.ladrunoThreads(1) == 1
    capfd.readouterr()
