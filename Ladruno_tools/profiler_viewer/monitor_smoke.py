"""
monitor_smoke.py  —  binary-free smoke test for ProfilerMonitor.

Runs anywhere h5py is present (no OpenSees build needed): a fake `ops` shim
implements just the verbs the monitor calls. `profiler('report')` writes a real
untimed-'root' rollup via the schema, so the phase-breakdown path is exercised
end to end (including the untimed-root backfill in ProfilerResults).

Run:  python monitor_smoke.py
"""
import io
import os
import tempfile

import h5py
from profiler_schema import write_run
from profiler_monitor import ProfilerMonitor


class FakeOps:
    """Minimal stand-in for the openseespy module."""
    def __init__(self, fail_at=None):
        self._t = 0.0
        self._peak = 0
        self._fail_at = fail_at
        self._nstep = 0
        self._armed = False

    def profiler(self, *args):
        verb = args[0]
        if verb == "start":
            self._armed = True
        elif verb == "stop":
            self._armed = False
        elif verb == "memory":
            return self._peak
        elif verb == "report":
            path = args[1]
            w = lambda ms: int(ms * 1e6)
            # untimed 'root' anchor over a 'step' subtree (real-engine shape)
            rollup = {
                "name": "root", "calls": 0, "wall_ns": 0,
                "wall_ns_min": 0, "wall_ns_max": 0, "cpu_ns": 0,
                "children": [{
                    "name": "step", "calls": self._nstep, "wall_ns": w(100),
                    "wall_ns_min": 0, "wall_ns_max": 0, "cpu_ns": w(100),
                    "children": [
                        {"name": "linearSolve", "calls": self._nstep, "wall_ns": w(70),
                         "wall_ns_min": 0, "wall_ns_max": 0, "cpu_ns": w(70), "children": []},
                        {"name": "formTangent", "calls": self._nstep, "wall_ns": w(30),
                         "wall_ns_min": 0, "wall_ns_max": 0, "cpu_ns": w(30), "children": []},
                    ],
                }],
            }
            with h5py.File(path, "w") as f:
                write_run(f, "snap", meta=dict(nDOF=20, nSteps=self._nstep,
                                               wall_ms_total=100.0), rollup=rollup)
        return 0

    def analyze(self, n, *args):
        self._nstep += 1
        self._t += (args[0] if args else 1.0)
        self._peak += 512
        if self._fail_at and self._nstep >= self._fail_at:
            return -3
        return 0

    def testIter(self):
        return 3

    def testNorms(self):
        return [1.0, 1e-4, 1e-9]

    def getTime(self):
        return self._t

    def systemSize(self):
        return 20


def test_happy_path():
    out = io.StringIO()
    ops = FakeOps()
    mon = ProfilerMonitor(ops, total=12, phase_every=4, stream=out)
    n = mon.run(12, 0.01, deep=True, memory=True,
                h5=os.path.join(tempfile.mkdtemp(), "m.h5"), run="A")
    text = out.getvalue()
    assert n == 12, f"expected 12 steps, got {n}"
    assert "system formed: nDOF=20" in text
    assert "phases @ step 4" in text and "linearSolve" in text, text
    assert "done: 12 step(s)" in text
    # peak grew (12 * 512 B) and is reported in the summary
    assert "peak" in text


def test_divergence_stops():
    out = io.StringIO()
    ops = FakeOps(fail_at=5)
    mon = ProfilerMonitor(ops, total=100, phase_every=0, stream=out)
    n = mon.run(100, 0.01, memory=False)
    text = out.getvalue()
    assert n == 5, f"expected stop at step 5, got {n}"
    assert "did NOT converge" in text, text


def main():
    test_happy_path()
    test_divergence_stops()
    print("OK — monitor smoke (happy path + phase snapshot + divergence stop)")


if __name__ == "__main__":
    main()
