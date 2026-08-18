"""`printA` / `printB` must never KILL the interpreter, however the SOE got left unsized. CI gate.

THE REPORT (2026-08-18, found while debugging ADR-85 T1b)
    `ops.printA("-ret")` hard-crashed the Python interpreter on an elastic deck
    whenever `constraints LadrunoContact` was the active handler.

THE MECHANISM, and why it looked contact-specific
    `printA` ends in `LinearSOE::getA()`.  Only TWO classes in the whole
    `SRC/system_of_eqn/linearSOE` tree override that virtual -- `FullGenLinSOE`
    and `DiagonalSOE` -- and both used to `exit(-1)` when their `matA` wrapper
    was still null.  `matA` is allocated in `setSize()`, and `setSize()` is
    reached from `StaticAnalysis::domainChanged()` only AFTER
    `ConstraintHandler::handle()` has returned >= 0.  So any route that leaves
    the SOE unsized turns the standard tangent-extraction probe into
    whole-process death.  Every other `system` was already safe by accident:
    `LinearSOE::getA()` returns 0 by default and `printA` handles null.

    The contact handler is what makes that reachable in practice.  Under the
    upstream handlers `handle()` essentially never fails.  But the contact
    subsystem's entire ADR-78/ADR-85 abort discipline is BUILT on returning -1
    from `handle()` for a refused contact (`ladrunoContactFatal()`), so
    `domainChanged()` failing is a routine, designed outcome -- and the very
    next thing a user reaches for to debug the refusal is `printA`, which took
    the interpreter down with it.  Hence "printA crashes under LadrunoContact".

    NOT a contact bug, then: a pre-existing dense-SOE landmine that the contact
    handler is uniquely good at stepping on.  The `before_analyze` rows below
    fire under `Plain` too, and they are kept precisely to say so.

SYMPTOM THAT MAKES THIS EXPENSIVE TO BISECT -- identical to its sibling file
    `exit(-1)` is whole-process death with NO Python traceback and NO exception;
    the FATAL text goes to `opserr`, which the Python module redirects, so
    nothing prints.  `faulthandler` stays silent because it is a clean `exit()`,
    not a signal.  It looks exactly like a bug in the test file.  See
    tests/test_soe_zero_free_equations.py, which is the same `exit(-1)` family
    reached by the other door (a model with zero free equations, where
    `setSize()` DID run but skipped its wrapper block).

    That kinship is why every case here runs in a CHILD interpreter: on a
    regression the child dies and the parent reports an ordinary test failure
    instead of the whole pytest run being killed.  The two Windows traps in the
    subprocess bootstrap (an oversized PYTHONPATH, an inherited stdin handle)
    are commented at their call sites and were paid for by the sibling file.

WHAT EACH ROW GATES
    before_analyze : printA/printB with no analyze() at all -- the SOE has never
                     been sized.  Pre-fix: DIES for FullGeneral and Diagonal
                     under EVERY handler.  Post-fix: warns, returns empty.
    refused_contact: a contact the handler REFUSES (kn = 0) -> handle() returns
                     -1 -> domainChanged() fails -> analyze() != 0 -> the SOE is
                     unsized.  This is the reported crash's actual route, and
                     the only row that needs contact at all.
    ok             : the falsifier on the falsifiers.  A successful analyze must
                     still return the REAL tangent, and `LadrunoContact` with no
                     contact declared must return it bit-for-bit identically to
                     `Plain` (the handler's standing graph-neutrality claim).
                     Without this row, "printA returns empty" would pass
                     everything.

Upstream elements/material only (quad / stdBrick + ElasticIsotropic) => zone_a.
"""
import json
import os
import subprocess
import sys

import pytest

from _testbed import ops

pytestmark = [pytest.mark.zone_a]

# The child must load the SAME binary the parent resolved -- hand it that one
# directory explicitly. Do NOT forward the parent's whole `sys.path` as
# PYTHONPATH: under the full suite that string grows past Windows' 32 KB
# environment-block limit and every `CreateProcess` dies with
# `OSError: [WinError 50] The request is not supported` -- which looks exactly
# like the crash this file is meant to detect.
ENGINE_DIR = os.path.dirname(os.path.abspath(ops.__file__))

# EVERY `system` reachable in a serial desktop build. Parametrizing over all of
# them is load-bearing, not thoroughness theatre: the first cut of this gate
# covered only FullGeneral/Diagonal/UmfPack -- the three the `getA()` survey
# pointed at -- and it PASSED while six other systems still killed the
# interpreter through `printB`, and while `SparseSYM` still killed it through
# `printA` by a completely different mechanism (a null deref in `zeroA()`,
# reached from `formTangent()` before `getA()` is ever called). Measured
# pre-fix on an unsized SOE: printA died on 2 of 10, printB on 6 of 10.
SYSTEMS_ALL = [
    "FullGeneral",    # overrides getA() -- died in printA (exit)
    "Diagonal",       # overrides getA() -- died in printA (exit)
    "SparseSYM",      # died in printA via zeroA() null deref (0xC0000005)
    "BandGeneral",    # died in printB
    "BandSPD",        # died in printB
    "ProfileSPD",     # died in printB
    "SProfileSPD",    # died in printB
    "SparseGeneral",  # died in printB
    "Pardiso",        # died in printB
    "UmfPack",        # the control: never died, must stay that way
]

# Kept as documentation of the survey result that scoped the getA() half of the
# fix: these are the ONLY two classes in SRC/system_of_eqn/linearSOE that
# override getA() at all. Everything else inherits LinearSOE's safe `return 0`.
SYSTEMS_GETA = ["FullGeneral", "Diagonal"]

# The POSITIVE controls need a system that can actually solve the deck. That
# excludes `Diagonal`, for the reason its sibling file already records: it only
# ever solves the DIAGONAL of the tangent, so a Newton static on a real quad is
# Jacobi iteration and does not converge (measured: rc = -3). It still belongs in
# every unsized-SOE row above -- that is the path it can crash on, and the only
# one that matters for this gate.
SYSTEMS_SOLVABLE = [s for s in SYSTEMS_ALL if s != "Diagonal"]

# Only FullGeneral can be asked for a full dense tangent: every other system here
# either inherits `LinearSOE::getA()`'s `return 0` (so printA is legitimately
# empty even after a GOOD analyze -- upstream behaviour, not a bug) or reports
# only its diagonal.
SYSTEMS_DENSE_A = ["FullGeneral"]

HANDLERS = ["LadrunoContact", "Plain"]

CHILD = r'''
import json, os, sys

_D = %(ENGINE_DIR)r
if os.path.isdir(_D):
    os.environ["PATH"] = _D + os.pathsep + os.environ.get("PATH", "")
    # add_dll_directory is WINDOWS-ONLY. Probe it with getattr rather than
    # calling it inside try/except OSError: the absence is an AttributeError,
    # which that clause does not catch, so on Linux CI the child would die
    # before importing anything and every case would report "the crash is back".
    _add = getattr(os, "add_dll_directory", None)
    if _add is not None:
        try:
            _add(_D)
        except OSError:
            pass
    sys.path.insert(0, _D)

try:
    import opensees as ops
except ModuleNotFoundError:
    import openseespy.opensees as ops

HANDLER = sys.argv[1]
SYS = sys.argv[2]
SCEN = sys.argv[3]

E, NU = 1000.0, 0.25


def build_quad():
    """One 2D quad, bottom edge fixed. 4 free equations."""
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for t, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)], 1):
        ops.node(t, float(x), float(y))
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("quad", 1, 1, 2, 3, 4, 1.0, "PlaneStrain", 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 0.0, -1.0)


def build_brick_refused_contact():
    """3D cube + a contact the HANDLER refuses (kn = 0).

    The refusal is handle()-time on purpose: a declaration-time refusal would
    raise OpenSeesError in the parser and never leave the SOE unsized, which is
    the state under test. `contact ... 0.0` trips
    "segment contact needs kn > 0 ... ABORTING (ADR-78 P1)".
    """
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    crds = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    for t, (x, y, z) in enumerate(crds, 1):
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 0.0, 0.0, -1.0)
    ops.contactSurface(1, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(2, "-slave", 1, 2)
    ops.contact(1, 1, 2, 0.0, 0.0, 0.0)      # kn = 0 => refused in handle()


def setup():
    ops.constraints(HANDLER)
    ops.numberer("Plain")
    ops.system(SYS)
    ops.test("NormDispIncr", 1.0e-10, 20, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")


out = {"handler": HANDLER, "system": SYS, "scenario": SCEN}

if SCEN == "before_analyze":
    build_quad()
    setup()
    # No analyze() -> setSize() has never run -> matA/vectB are null.
    out["A"] = list(ops.printA("-ret"))
    out["B"] = list(ops.printB("-ret"))
elif SCEN == "refused_contact":
    build_brick_refused_contact()
    setup()
    try:
        out["rc"] = ops.analyze(1)
    except Exception as ex:
        out["rc"] = "exception: " + repr(ex)
    # The refusal left the SOE unsized. THIS is the reported crash.
    out["A"] = list(ops.printA("-ret"))
    out["B"] = list(ops.printB("-ret"))
else:  # "ok"
    build_quad()
    setup()
    out["rc"] = ops.analyze(1)
    out["A"] = list(ops.printA("-ret"))
    out["B"] = list(ops.printB("-ret"))

out["engine"] = os.path.dirname(os.path.abspath(ops.__file__))
print("RESULT " + json.dumps(out))
'''


def _run_child(handler, system, scenario):
    script = CHILD % {"ENGINE_DIR": ENGINE_DIR}
    proc = subprocess.run(
        [sys.executable, "-c", script, handler, system, scenario],
        # stdin=DEVNULL is load-bearing on Windows: with the inherited stdin,
        # `Popen` tries to `DuplicateHandle` whatever pytest's capture left
        # there and raises `OSError: [WinError 50]` -- but only once some OTHER
        # module in the suite has run, so the file passes standalone and every
        # case fails under `pytest tests/`. That failure is indistinguishable
        # here from the crash this file exists to detect.
        stdin=subprocess.DEVNULL,
        capture_output=True, text=True, timeout=300,
        encoding="utf-8", errors="replace",
        # the splash banner's UTF-8 glyphs raise UnicodeDecodeError in a
        # text-mode capture on a cp1252 console
        env=dict(os.environ, LADRUNO_OPENSEES_QUIET="1"),
    )
    assert proc.returncode == 0, (
        f"{handler}/{system}/{scenario}: child process DIED (exit "
        f"{proc.returncode}) -- printA/printB can kill the interpreter again. "
        f"An unsized SOE must warn and return empty, not exit().\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )
    for line in proc.stdout.splitlines():
        if line.startswith("RESULT "):
            out = json.loads(line[len("RESULT "):])
            # a child that quietly loaded a DIFFERENT binary would report on
            # something nobody built -- the classic stale-engine false pass
            assert os.path.normcase(out["engine"]) \
                == os.path.normcase(ENGINE_DIR), (
                f"{handler}/{system}/{scenario}: child loaded {out['engine']}, "
                f"parent is {ENGINE_DIR}"
            )
            return out
    raise AssertionError(
        f"{handler}/{system}/{scenario}: no RESULT line.\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )


@pytest.mark.parametrize("system", SYSTEMS_ALL)
@pytest.mark.parametrize("handler", HANDLERS)
def test_printa_before_analyze_survives(handler, system):
    """printA/printB with the SOE never sized: must return empty, not exit()."""
    out = _run_child(handler, system, "before_analyze")
    # The contract is "no tangent to report", expressed as an empty result --
    # NOT a fabricated matrix. Anything nonempty here would mean printA read a
    # matrix that was never assembled.
    assert out["A"] == [], out
    assert out["B"] == [], out


@pytest.mark.parametrize("system", SYSTEMS_ALL)
def test_printa_after_refused_contact_survives(system):
    """The reported crash: a refused contact fails handle(), then printA."""
    out = _run_child("LadrunoContact", system, "refused_contact")
    # analyze() MUST have failed -- if it ever starts succeeding, the deck no
    # longer reproduces the unsized-SOE state and this row is vacuous.
    assert out["rc"] != 0, (
        "the kn = 0 contact was NOT refused, so the SOE got sized after all "
        f"and this row no longer tests anything: {out}"
    )
    assert out["A"] == [], out
    assert out["B"] == [], out


@pytest.mark.parametrize("system", SYSTEMS_DENSE_A)
def test_printa_still_returns_the_real_tangent(system):
    """Falsifier on the falsifiers: a good analyze still yields a real matrix."""
    out = _run_child("LadrunoContact", system, "ok")
    assert out["rc"] == 0, out
    assert len(out["A"]) == 16, out            # dense 4x4 over 4 free equations
    assert any(v != 0.0 for v in out["A"]), out
    assert len(out["B"]) == 4, out


@pytest.mark.parametrize("system", SYSTEMS_SOLVABLE)
def test_printb_still_returns_the_real_rhs(system):
    """The `getX`/`getB` half needs its own falsifier, on EVERY solvable system.

    The empty-Vector fallback is a shared function-local static. If it ever
    shadowed the live path -- or if a future edit returned it unconditionally --
    every unsized row above would still pass while `printB` silently reported an
    empty right-hand side on a perfectly good model. This is the row that
    notices. 4 free equations => 4 entries, not all zero (the deck is loaded).
    """
    out = _run_child("LadrunoContact", system, "ok")
    assert out["rc"] == 0, out
    assert len(out["B"]) == 4, out
    assert any(v != 0.0 for v in out["B"]), out


def test_contact_handler_tangent_matches_plain():
    """`LadrunoContact` with NO contact declared is graph- and tangent-neutral.

    Pins the fix to "printA works again", not "printA returns something". A
    handler that injected a stray adapter or numbered differently would show up
    here as a changed tangent even though every row above still passed.
    """
    a = _run_child("LadrunoContact", "FullGeneral", "ok")
    b = _run_child("Plain", "FullGeneral", "ok")
    assert a["rc"] == 0 and b["rc"] == 0, (a, b)
    assert a["A"] == b["A"], (
        "the contact handler's tangent diverged from PlainHandler's on a deck "
        f"with no contact declared:\ncontact={a['A']}\nplain={b['A']}"
    )
    assert a["B"] == b["B"], (a["B"], b["B"])
