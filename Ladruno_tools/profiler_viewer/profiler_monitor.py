"""
Ladruno Stack Profiler — live console monitor (no rebuild required).

Wraps an OpenSees ``analyze()`` loop and shows, in real time:

  * convergence health  — per-step iterations + achieved residual norm + dt
  * speed / ETA         — wall per step, steps/s, estimated time remaining
  * memory growth       — profiler peak bytes (delta since start)
  * phase breakdown     — cumulative formTangent / linearSolve / elem.* shares,
                          snapshotted every N steps

It uses ONLY what the shipped ``profiler`` command and standard openseespy
introspection already expose (``testIter`` / ``testNorms`` / ``getTime`` /
``systemSize`` / ``profiler('memory')`` / ``profiler('report')``), so it runs on
today's binary with no C++ change. The phase breakdown reads a throw-away
mid-run ``profiler('report')`` snapshot; it needs ``h5py`` in the analysis
interpreter and silently disables itself (with a one-line note) if absent.

Two ways to use it
------------------
Drive the loop yourself (full control over analyze args / substepping)::

    from profiler_monitor import ProfilerMonitor
    mon = ProfilerMonitor(ops, total=1000, phase_every=100)
    mon.start(deep=True, memory=True)
    for i in range(1000):
        code = ops.analyze(1, 0.01)          # Transient: (nsteps, dt)
        if not mon.step(code):               # False => the step did not converge
            break
    mon.stop("run.h5", run="caseA")          # stop + (optional) final report

Or let the monitor run the loop::

    ProfilerMonitor(ops, total=1000).run(1000, 0.01, deep=True, memory=True)

``ops`` is passed in (not imported) so this module is decoupled from however the
Ladruno ``opensees`` build is put on ``sys.path``.
"""

from __future__ import annotations

import os
import sys
import time
import math
import tempfile
import contextlib

# Windows + HDF5: the reader must not take a file lock on the writer's file.
os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

# Structural wrapper nodes that are just sums of their children — excluded from
# the "which phase is hot" line (we want the actual work leaves).
_WRAPPER_NODES = {"root", "step", "solveCurrentStep"}


def _fmt_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if abs(n) < 1024.0:
            return f"{n:5.1f}{unit}"
        n /= 1024.0
    return f"{n:5.1f}TB"


def _fmt_hms(secs: float) -> str:
    if not math.isfinite(secs) or secs < 0:
        return "--:--"
    secs = int(secs)
    h, rem = divmod(secs, 3600)
    m, s = divmod(rem, 60)
    return f"{h:d}:{m:02d}:{s:02d}" if h else f"{m:02d}:{s:02d}"


@contextlib.contextmanager
def _suppress_cstderr():
    """Redirect fd 2 to the null device for the duration of the block.

    The engine's ``mergedRollup()`` prints a benign 'called while still enabled'
    warning on a mid-run ``report`` (single-threaded runs are safe — the analysis
    is idle between ``analyze`` calls). Swallow just that fd so the live line
    stays clean; real Python-level exceptions still propagate.
    """
    try:
        saved = os.dup(2)
    except (OSError, ValueError):
        yield  # no real fd 2 (e.g. some notebooks); nothing to redirect
        return
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        sys.stderr.flush()
        os.dup2(devnull, 2)
        yield
    finally:
        sys.stderr.flush()
        os.dup2(saved, 2)
        os.close(devnull)
        os.close(saved)


class ProfilerMonitor:
    def __init__(self, ops, total: int | None = None, phase_every: int = 0,
                 stream=None, ema: float = 0.2):
        """
        ops         the openseespy module (or Tcl-like shim with the same verbs)
        total       total step count, for ETA + progress (optional)
        phase_every snapshot the phase breakdown every N steps (0 = never)
        stream      where the live line goes (default sys.stderr)
        ema         smoothing factor for the ms/step estimate (0..1)
        """
        self.ops = ops
        self.total = total
        self.phase_every = int(phase_every)
        self.stream = stream if stream is not None else sys.stderr
        self.ema = float(ema)

        self.n = 0                    # completed steps
        self._t0 = None               # wall at start()
        self._tlast = None            # wall at previous step
        self._ms_ema = None           # smoothed ms/step
        self._peak0 = 0               # peak bytes at start()
        self._peak = 0
        self._memory = False
        self._isatty = bool(getattr(self.stream, "isatty", lambda: False)())
        self._last_line_len = 0
        self._diverged = False
        self._phase_ok = self._phase_available()
        self._phase_note_shown = False

    # -- lifecycle ----------------------------------------------------------
    def start(self, deep: bool = False, memory: bool = False,
              per_step: bool = True):
        flags = []
        if deep:
            flags.append("-deep")
        if memory:
            flags.append("-memory")
        if per_step:
            flags.append("-perStep")
        self.ops.profiler("start", *flags)
        self._memory = memory
        self.n = 0
        self._t0 = self._tlast = time.perf_counter()
        self._ms_ema = None
        self._peak0 = self._peak = self._read_peak() if memory else 0
        self._diverged = False
        self._ndof_shown = False
        hdr = f"[profiler-monitor] armed{' -deep' if deep else ''}" \
              f"{' -memory' if memory else ''}"
        if self.phase_every and not self._phase_ok:
            hdr += "  (phase breakdown off: h5py not importable here)"
        print(hdr, file=self.stream, flush=True)
        return self

    def step(self, code: int | None = None) -> bool:
        """Call once after each ``analyze()``. Returns False if the step failed
        to converge (nonzero return code), so a caller can break the loop."""
        now = time.perf_counter()
        dwall = now - self._tlast
        self._tlast = now
        self.n += 1

        ms = dwall * 1e3
        self._ms_ema = ms if self._ms_ema is None else \
            self.ema * ms + (1 - self.ema) * self._ms_ema
        if self._memory:
            self._peak = self._read_peak()

        # nDOF is only known once the first analyze has formed the SOE.
        if not self._ndof_shown:
            nd = self._ndof()
            if nd > 0:
                print(f"[profiler-monitor] system formed: nDOF={nd}",
                      file=self.stream, flush=True)
                self.ndof = nd
                self._ndof_shown = True

        if code is not None and code != 0:
            self._diverged = True
            self._render(code=code)
            print("", file=self.stream, flush=True)  # break the \r line
            print(f"[profiler-monitor] step {self.n} did NOT converge "
                  f"(analyze returned {code})", file=self.stream, flush=True)
            return False

        self._render(code=0)
        if self.phase_every and self._phase_ok and self.n % self.phase_every == 0:
            self._phase_snapshot()
        return True

    def stop(self, h5: str | None = None, run: str = "run"):
        self.ops.profiler("stop")
        if not self._diverged:
            print("", file=self.stream, flush=True)  # finish the \r line
        wall = time.perf_counter() - self._t0 if self._t0 else 0.0
        avg = (wall * 1e3 / self.n) if self.n else 0.0
        summ = (f"[profiler-monitor] done: {self.n} step(s) in {_fmt_hms(wall)} "
                f"({avg:.2f} ms/step avg")
        if self._memory:
            summ += f", peak {_fmt_bytes(self._peak)} "
            summ += f"(+{_fmt_bytes(self._peak - self._peak0)})"
        summ += ")"
        print(summ, file=self.stream, flush=True)
        if h5:
            self.ops.profiler("report", h5, "-run", run)
            print(f"[profiler-monitor] report -> {h5} (run '{run}'); "
                  f"view with:  python launch.py {h5}", file=self.stream, flush=True)

    def run(self, nsteps: int, *analyze_args, deep: bool = False,
            memory: bool = False, h5: str | None = None, run: str = "run"):
        """Convenience: start, ``analyze(1, *analyze_args)`` × nsteps with the
        live monitor, then stop (+ optional report). Stops early on divergence."""
        if self.total is None:
            self.total = nsteps
        self.start(deep=deep, memory=memory)
        for _ in range(nsteps):
            code = self.ops.analyze(1, *analyze_args)
            if not self.step(code):
                break
        self.stop(h5=h5, run=run)
        return self.n

    # -- rendering ----------------------------------------------------------
    def _render(self, code: int):
        try:
            iters = int(self.ops.testIter())
        except Exception:
            iters = -1
        resid = self._resid()
        t = self._time()
        ms = self._ms_ema or 0.0
        sps = (1e3 / ms) if ms > 0 else 0.0

        parts = []
        if self.total:
            pct = 100.0 * self.n / self.total
            parts.append(f"{self.n:>6}/{self.total} {pct:4.0f}%")
        else:
            parts.append(f"step {self.n:>6}")
        parts.append(f"t={t:.4g}")
        parts.append(f"it={iters:>2}" if iters >= 0 else "it= ?")
        parts.append(f"|r|={resid:.1e}" if math.isfinite(resid) else "|r|=  n/a")
        parts.append(f"{ms:6.1f}ms/st")
        parts.append(f"{sps:5.0f}st/s")
        if self.total:
            eta = (self.total - self.n) * ms / 1e3
            parts.append(f"ETA {_fmt_hms(eta)}")
        if self._memory:
            parts.append(f"peak {_fmt_bytes(self._peak)}")

        line = " | ".join(parts)
        if self._isatty:
            pad = max(0, self._last_line_len - len(line))
            self.stream.write("\r" + line + " " * pad)
            self.stream.flush()
            self._last_line_len = len(line)
        else:
            # non-tty (piped/log/notebook): throttle to avoid a wall of lines
            if self.total is None or self.n == 1 or self.n == self.total \
                    or (self.phase_every and self.n % self.phase_every == 0):
                print(line, file=self.stream, flush=True)

    def _phase_snapshot(self):
        top = self._read_phase_shares()
        if not top:
            return
        cells = " | ".join(f"{name} {share:.0%}" for name, share in top)
        prefix = "\n" if self._isatty else ""
        print(f"{prefix}  phases @ step {self.n}: {cells}",
              file=self.stream, flush=True)
        self._last_line_len = 0  # next \r starts a fresh line

    # -- engine/openseespy access ------------------------------------------
    def _ndof(self) -> int:
        try:
            return int(self.ops.systemSize())
        except Exception:
            return -1

    def _read_peak(self) -> int:
        try:
            return int(self.ops.profiler("memory"))
        except Exception:
            return 0

    def _time(self) -> float:
        try:
            return float(self.ops.getTime())
        except Exception:
            return float("nan")

    def _resid(self) -> float:
        try:
            norms = self.ops.testNorms()
            if norms:
                return float(norms[-1])
        except Exception:
            pass
        return float("nan")

    # -- phase breakdown via a throw-away mid-run report --------------------
    def _phase_available(self) -> bool:
        def _try():
            import h5py  # noqa: F401
            from profiler_results import ProfilerResults  # noqa: F401
            return True
        try:
            return _try()
        except Exception:
            pass
        # A model run started with `python -S` (the common way to dodge a stale
        # boot-.pth pyd preload) has no site-packages on the path, so h5py is
        # absent. Add the site dir DIRECTLY (sys.path.append does NOT process
        # .pth files, so this does not re-trigger the preload) and retry.
        try:
            import sysconfig
            purelib = sysconfig.get_paths().get("purelib")
            if purelib and purelib not in sys.path:
                sys.path.append(purelib)
            return _try()
        except Exception:
            return False

    def _read_phase_shares(self, k: int = 4):
        if not self._phase_ok:
            return []
        from profiler_results import ProfilerResults
        tmp = os.path.join(tempfile.gettempdir(),
                           f"_ops_monitor_{os.getpid()}_{self.n}.h5")
        try:
            with _suppress_cstderr():
                self.ops.profiler("report", tmp, "-run", "snap")
            pr = ProfilerResults(tmp)
            root = pr.rollup("snap")
            pr.close()
        except Exception:
            if not self._phase_note_shown:
                print("\n[profiler-monitor] phase snapshot unavailable "
                      "(report/read failed); continuing without it.",
                      file=self.stream, flush=True)
                self._phase_note_shown = True
                self._phase_ok = False
            return []
        finally:
            with contextlib.suppress(OSError):
                os.remove(tmp)

        # Flatten to work leaves (exclude the pure-sum wrappers), rank by share.
        # Explicit stack walk (avoid recursion-depth surprises on deep trees).
        flat = {}
        work = [root]
        while work:
            nd = work.pop()
            name = nd["name"]
            if name not in _WRAPPER_NODES:
                # keep the largest share seen for a given phase name
                flat[name] = max(flat.get(name, 0.0), nd.get("share", 0.0))
            work.extend(nd.get("children", []))
        ranked = sorted(flat.items(), key=lambda kv: kv[1], reverse=True)
        return ranked[:k]
