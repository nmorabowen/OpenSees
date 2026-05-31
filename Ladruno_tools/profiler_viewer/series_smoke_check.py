"""P0#3 per-step series smoke -- checker (venv python, h5py).

Reads series_smoke.h5 via ProfilerResults.series() and asserts the per-step
hook populated the time history:
  * the run HAS a /series group (None means the hook never fired)
  * exactly N rows (one per committed step)
  * step index is strictly increasing (getCommitTag is monotonic)
  * dt == DT for every step (transient time step round-trips)
  * the Newton iteration count is recorded and positive on at least one step

    python series_smoke_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

from profiler_results import ProfilerResults

N = 20
DT = 0.01


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    path = os.path.join(out, "series_smoke.h5")
    problems = 0

    with ProfilerResults(path) as pr:
        s = pr.series("transient")

    def check(cond, msg):
        nonlocal problems
        print(f"  [{'OK' if cond else 'FAIL'}] {msg}")
        if not cond:
            problems += 1

    check(s is not None, "run has a /series group (per-step hook fired)")
    if s is None:
        print("\nSERIES_SMOKE_CHECK: 1 PROBLEM(S)")
        return 1

    step, dt, iters, t = s["step"], s["dt"], s["iters"], s["t"]
    print(f"  rows={len(step)}  step[0..]={step[:3]}…{step[-1:]}  "
          f"dt0={dt[0] if dt else None}  iters(max)={max(iters) if iters else None}")

    check(len(step) == N, f"{N} rows recorded (got {len(step)})")
    check(all(step[i] < step[i + 1] for i in range(len(step) - 1)),
          "step index strictly increasing (monotonic commit tag)")
    check(len(dt) == N and all(abs(d - DT) < 1e-12 for d in dt),
          f"every dt == {DT}")
    check(len(iters) == N and max(iters) >= 1,
          "Newton iteration count recorded (>=1 on some step)")
    check(len(t) == N and all(t[i] <= t[i + 1] + 1e-12 for i in range(len(t) - 1)),
          "pseudo-time non-decreasing")

    print("\n" + "=" * 52)
    if problems == 0:
        print("SERIES_SMOKE_CHECK: ALL PASS (per-step series populated)")
        return 0
    print(f"SERIES_SMOKE_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
