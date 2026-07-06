"""ladruno_logs — per-run analysis logging for the Ladruno fork.

Every analysis run should leave a self-contained record on disk. This module
gives that convention a blessed helper: an ``AnalysisLog`` context manager that

  1. creates a ``logs/`` folder next to the run (configurable),
  2. tees ALL terminal output into ``logs/log.log`` — the Python side
     (``print``/tracebacks) via a stdout/stderr tee, and the OpenSees C++ side
     (``opserr``: warnings, convergence failures, banner) via the interpreter's
     own ``logFile`` command in ``-append`` mode on the same file, and
  3. writes a model-info block (ndm/ndf, node/element counts, per-type element
     breakdown, equation count) plus an analysis run-time summary (per-phase
     wall times, total ``analyze()`` time, total wall time) into the log.

Usage (OpenSeesPy, local build or installed openseespy)::

    import opensees as ops
    from ladruno_logs import AnalysisLog

    with AnalysisLog(ops) as log:              # ./logs/log.log starts here
        build_model()
        log.model_info()                       # optional; auto-printed at exit
        with log.timed("gravity"):
            ops.analyze(10)
        with log.timed("pushover"):
            ok = log.analyze(400)              # timed passthrough to ops.analyze

Everything printed inside the block reaches BOTH the console and the log
(console behaviour is unchanged). An exception inside the block is written to
the log (full traceback) before it propagates, and the summary still records
the run as FAILED with the elapsed time — a crashed run leaves a diagnosable
log, not a truncated one.

Tcl runs get the C++ half natively — put ``file mkdir logs; logFile logs/log.log``
at the top of the script (that captures opserr AND ``puts`` goes to the console
only; use ``puts $logfd`` patterns or run under ``tee`` for the rest).

Pure stdlib + an ``ops`` handle. No OpenSees source change, no rebuild.
"""
from __future__ import annotations

import datetime
import io
import sys
import time
import traceback
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

_RULE = "-" * 66


class _Tee(io.TextIOBase):
    """Write-through to the real console stream AND the log file handle."""

    def __init__(self, console, fh):
        self._console = console
        self._fh = fh

    def write(self, s):
        n = self._console.write(s)
        try:
            self._fh.write(s)
            self._fh.flush()      # keep file current even if the run dies
        except ValueError:        # log already closed (interpreter teardown)
            pass
        return n

    def flush(self):
        self._console.flush()
        try:
            self._fh.flush()
        except ValueError:
            pass

    def isatty(self):
        return self._console.isatty()

    @property
    def encoding(self):
        return getattr(self._console, "encoding", "utf-8")


class AnalysisLog:
    """Context manager: tee terminal output to ``<log_dir>/<filename>`` and
    record model info + analysis run time for one OpenSees run.

    ops      : the OpenSees module handle (opensees / openseespy.opensees).
    log_dir  : folder for the log (created if missing). Default ``"logs"``.
    filename : log file name inside ``log_dir``. Default ``"log.log"``.
    append   : append to an existing log instead of starting fresh.
    tag      : optional free-text run label written into the header.
    """

    def __init__(self, ops, log_dir="logs", filename="log.log",
                 append=False, tag=None):
        self.ops = ops
        self.path = Path(log_dir) / filename
        self._append = append
        self._tag = tag
        self._fh = None
        self._old = None
        self._t0 = None
        self._phases = []           # [(label, seconds), ...]
        self._analyze_s = 0.0
        self._analyze_calls = 0
        self._analyze_steps = 0     # steps REQUESTED via log.analyze(n, ...)
        self._analyze_fails = 0
        self._info_printed = False

    # -- lifecycle -----------------------------------------------------------
    def __enter__(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(self.path, "a" if self._append else "w",
                        encoding="utf-8")
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        head = [_RULE, "LADRUNO analysis log"]
        if self._tag:
            head.append("run      : " + str(self._tag))
        head += ["started  : " + now,
                 "script   : " + (sys.argv[0] or "<interactive>"),
                 "cwd      : " + str(Path.cwd()),
                 _RULE, ""]
        self._fh.write("\n".join(head))
        self._fh.flush()
        self._old = (sys.stdout, sys.stderr)
        sys.stdout = _Tee(self._old[0], self._fh)
        sys.stderr = _Tee(self._old[1], self._fh)
        # OpenSees C++ stream (opserr) -> same file. -append so it does not
        # truncate the header; echo stays on so the console is unchanged.
        try:
            self.ops.logFile(str(self.path), "-append")
        except Exception as e:              # log still works Python-side
            print(f"[ladruno_logs] WARNING: ops.logFile failed ({e}); "
                  "C++-side output will reach the console only")
        self._t0 = time.perf_counter()
        return self

    def __exit__(self, exc_type, exc, tb):
        status = "OK"
        if exc_type is not None:
            status = f"FAILED ({exc_type.__name__}: {exc})"
            traceback.print_exception(exc_type, exc, tb, file=sys.stderr)
        if not self._info_printed:
            self.model_info()
        self._summary(status)
        sys.stdout, sys.stderr = self._old
        self._fh.close()
        return False                        # propagate any exception

    # -- model info ----------------------------------------------------------
    def model_info(self):
        """Print (console + log) a summary of the current Domain."""
        print("\n" + _RULE + "\nmodel info\n" + _RULE)
        try:
            nodes = self.ops.getNodeTags() or []
            eles = self.ops.getEleTags() or []
            ndm, ndf = self._scalar("getNDM"), self._scalar("getNDF")
            print(f"  ndm / ndf     : {ndm} / {ndf}")
            print(f"  nodes         : {len(nodes):,}")
            print(f"  elements      : {len(eles):,}")
            for name, n in self._ele_breakdown(eles):
                print(f"    {name:<28s} {n:>10,}")
            try:
                neq = self.ops.systemSize()
                print(f"  equations     : {int(neq):,}" if neq else
                      "  equations     : n/a (system not yet sized — "
                      "sized on first analyze)")
            except Exception:
                print("  equations     : n/a (no analysis defined)")
        except Exception as e:              # never let logging kill a run
            print(f"  [ladruno_logs] model query failed: {e}")
        print(_RULE)
        self._info_printed = True

    def _scalar(self, cmd):
        try:
            v = getattr(self.ops, cmd)()
            return v[0] if isinstance(v, (list, tuple)) else v
        except Exception:
            return "?"

    def _ele_breakdown(self, eles):
        """[(typeName, count)] by element class, largest first."""
        try:
            cts = self.ops.getEleClassTags() or []
            if len(cts) != len(eles):
                return []
            counts = Counter(cts)
            first = {}
            for tag, ct in zip(eles, cts):
                first.setdefault(ct, tag)
            rows = []
            for ct, n in counts.most_common():
                try:
                    name = self.ops.eleType(first[ct])
                except Exception:
                    name = f"classTag {ct}"
                rows.append((name, n))
            return rows
        except Exception:
            return []

    # -- timing --------------------------------------------------------------
    def analyze(self, *args):
        """Timed passthrough to ``ops.analyze``; returns its code."""
        t = time.perf_counter()
        rc = self.ops.analyze(*args)
        self._analyze_s += time.perf_counter() - t
        self._analyze_calls += 1
        if args and isinstance(args[0], int):
            self._analyze_steps += args[0]
        if rc != 0:
            self._analyze_fails += 1
        return rc

    @contextmanager
    def timed(self, label):
        """Record the wall time of a named phase in the summary."""
        t = time.perf_counter()
        try:
            yield
        finally:
            self._phases.append((label, time.perf_counter() - t))

    # -- summary -------------------------------------------------------------
    def _summary(self, status):
        total = time.perf_counter() - self._t0
        now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        row = lambda k, v: print(f"  {k:<21s}: {v}")
        print("\n" + _RULE + "\nrun summary\n" + _RULE)
        for label, s in self._phases:
            row(f"phase {label!r}", f"{s:10.3f} s")
        if self._analyze_calls:
            fail = (f", {self._analyze_fails} non-zero return(s)"
                    if self._analyze_fails else "")
            row("analyze() wall time", f"{self._analyze_s:10.3f} s  "
                f"({self._analyze_calls} call(s), "
                f"{self._analyze_steps} step(s) requested{fail})")
        row("total wall time", f"{total:10.3f} s")
        row("finished", now)
        row("status", status)
        print(_RULE)


# --------------------------------------------------------------------------
# self-test: tiny 2D truss pushover under the logger. Checks the log file
# carries (a) Python-side prints, (b) C++-side opserr output, (c) the model
# info block with the element breakdown, and (d) the run-time summary.
# --------------------------------------------------------------------------
def _selftest():
    import os
    import tempfile

    try:
        import opensees as ops
    except ModuleNotFoundError:
        dist = os.environ.get(
            "LADRUNO_DIST",
            r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin")
        if os.path.isdir(dist):
            os.add_dll_directory(dist)
            sys.path.insert(0, dist)
        import opensees as ops

    work = Path(tempfile.mkdtemp(prefix="ladruno_logs_"))
    os.chdir(work)

    with AnalysisLog(ops, tag="self-test truss") as log:
        ops.wipe()
        ops.model("basic", "-ndm", 2, "-ndf", 2)
        ops.node(1, 0.0, 0.0); ops.node(2, 144.0, 0.0)
        ops.node(3, 168.0, 0.0); ops.node(4, 72.0, 96.0)
        for n in (1, 2, 3):
            ops.fix(n, 1, 1)
        ops.uniaxialMaterial("Elastic", 1, 3000.0)
        ops.element("Truss", 1, 1, 4, 10.0, 1)
        ops.element("Truss", 2, 2, 4, 5.0, 1)
        ops.element("Truss", 3, 3, 4, 5.0, 1)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
        ops.load(4, 100.0, -50.0)
        ops.constraints("Plain"); ops.numberer("RCM")
        ops.system("BandSPD"); ops.test("NormDispIncr", 1e-8, 10)
        ops.algorithm("Newton"); ops.integrator("LoadControl", 0.1)
        ops.analysis("Static")
        log.model_info()
        with log.timed("pushover"):
            rc = log.analyze(10)
        print(f"analyze rc={rc}, ux(4)={ops.nodeDisp(4, 1):.6f}")
        # provoke a C++-side opserr warning (duplicate node tag) so the test
        # proves the opserr stream lands in the same file
        try:
            ops.node(4, 0.0, 0.0)
        except Exception:
            pass

    text = Path("logs/log.log").read_text(encoding="utf-8")
    checks = {
        "header": "LADRUNO analysis log" in text,
        "python print": "analyze rc=0" in text,
        "model info": "model info" in text and "nodes         : 4" in text,
        "ele breakdown": "Truss" in text,
        "opserr tee": "WARNING" in text,
        "phase timing": "phase 'pushover'" in text,
        "summary": "run summary" in text and "status               : OK" in text,
    }
    for k, ok in checks.items():
        print(f"  check {k:<14s}: {'PASS' if ok else 'FAIL'}")
    print("log at:", work / "logs" / "log.log")
    all_ok = all(checks.values())
    print("SELF-TEST:", "PASS" if all_ok else "FAIL")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(_selftest())
