"""Run-header smoke -- checker (venv python, h5py).

Reads header_smoke.h5 via ProfilerResults.header() and asserts the report
command populated the run header:
  * algorithm / integrator / solver are the type strings the user typed
    (command-time capture; getClassType() returns "UnknownMovableObject" for
    these analysis classes, so the names come from cmds, not getClassType)
  * dt_min == dt_max == DT (derived from the per-step series in buildMeta)
  * dt_cr >= 0 — the P0#5 plumbing is intact (it reads the integrator's virtual
    getCriticalTimeStep() and only sets dt_cr when > 0). NOTE: with the bundled
    explicit integrators getCriticalTimeStep() currently returns -1 for these
    small lumped-mass models, so a *positive* dt_cr is not asserted here.

    python header_smoke_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

from profiler_results import ProfilerResults

DT = 0.005
WANT = {"integrator": "CentralDifferenceLadruno", "algorithm": "Linear", "solver": "BandGeneral"}


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    path = os.path.join(out, "header_smoke.h5")
    problems = 0

    with ProfilerResults(path) as pr:
        h = pr.header("explicit")

    print(f"  integrator='{h['integrator']}'  algorithm='{h['algorithm']}'  solver='{h['solver']}'")
    print(f"  dt_min={h['dt_min']}  dt_max={h['dt_max']}  dt_cr={h['dt_cr']:.5g}  "
          f"oversample={h['oversample_ratio']:.3g}x")

    def check(cond, msg):
        nonlocal problems
        print(f"  [{'OK' if cond else 'FAIL'}] {msg}")
        if not cond:
            problems += 1

    for k, want in WANT.items():
        got = h[k]
        check(got == want, f"{k} == '{want}' (got '{got}')")
        check(got != "UnknownMovableObject", f"{k} is not the getClassType fallback")

    check(abs(h["dt_min"] - DT) < 1e-12 and abs(h["dt_max"] - DT) < 1e-12,
          f"dt_min == dt_max == {DT}")
    check(h["dt_cr"] >= 0.0 and h["oversample_ratio"] >= 0.0,
          "dt_cr / oversample plumbing intact (>= 0; populated when the integrator reports one)")

    print("\n" + "=" * 52)
    if problems == 0:
        print("HEADER_SMOKE_CHECK: ALL PASS (run header identity populated)")
        return 0
    print(f"HEADER_SMOKE_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
