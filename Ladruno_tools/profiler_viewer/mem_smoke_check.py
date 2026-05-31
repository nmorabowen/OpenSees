"""P4 memory-counter smoke -- checker (venv python, h5py).

Reads mem_smoke.h5 via ProfilerResults.memory() and asserts the P4 per-type
Matrix/Vector counters were populated and are self-consistent:
  * matrix_live > 0 and vector_live > 0   (the analysis allocated both)
  * peak_bytes  > 0 and >= max(matrix_live, vector_live)
  * id_live == 0                           (ID counters deferred to a follow-up)

    python mem_smoke_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

from profiler_results import ProfilerResults


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    path = os.path.join(out, "mem_smoke.h5")
    problems = 0

    with ProfilerResults(path) as pr:
        mem = pr.memory("memcase")

    ml, vl, il, pk = (mem["matrix_live"], mem["vector_live"],
                      mem["id_live"], mem["peak_bytes"])
    print(f"  matrix_live={ml}  vector_live={vl}  id_live={il}  peak_bytes={pk}")

    def check(cond, msg):
        nonlocal problems
        if cond:
            print(f"  [OK] {msg}")
        else:
            print(f"  [FAIL] {msg}")
            problems += 1

    check(ml > 0, f"matrix_live > 0 (got {ml})")
    check(vl > 0, f"vector_live > 0 (got {vl})")
    check(pk > 0, f"peak_bytes > 0 (got {pk})")
    check(pk >= max(ml, vl), f"peak_bytes >= max(matrix_live, vector_live) ({pk} >= {max(ml, vl)})")
    check(il == 0, f"id_live == 0 (ID counters deferred) (got {il})")

    print("\n" + "=" * 52)
    if problems == 0:
        print("MEM_SMOKE_CHECK: ALL PASS (per-type Matrix/Vector counters populated)")
        return 0
    print(f"MEM_SMOKE_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
