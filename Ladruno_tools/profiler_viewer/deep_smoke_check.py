"""P3 deep per-element-type smoke -- checker (venv python, h5py).

Reads deep_smoke.h5 via ProfilerResults.rollup() and asserts the P3 per-element
buckets were populated:
  * an `elem.tangent` and an `elem.residual` node exist in the rollup tree
  * each carries a non-empty elem_by_type breakdown
  * the Truss classTag (ELE_TAG_Truss = 12) appears, and `count` is a positive
    whole multiple of N_TRUSS (every form-tangent / form-residual call evaluates
    all N_TRUSS elements, so count == N_TRUSS x n_calls; the exact multiplier
    depends on the Newton iteration count and is not hard-coded)
  * every row's wall_ns_per_elem >= 0 and fb_coupled is False (Truss is
    displacement-based)

    python deep_smoke_check.py <out_dir>
"""

from __future__ import annotations

import os
import sys

from profiler_results import ProfilerResults

ELE_TAG_TRUSS = 12
N_TRUSS = 20


def _find(node: dict, name: str):
    """Depth-first search for the first rollup node with the given name."""
    if node.get("name") == name:
        return node
    for ch in node.get("children", []):
        hit = _find(ch, name)
        if hit is not None:
            return hit
    return None


def main() -> int:
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    path = os.path.join(out, "deep_smoke.h5")
    problems = 0

    with ProfilerResults(path) as pr:
        root = pr.rollup("deepcase")

    def check(cond, msg):
        nonlocal problems
        if cond:
            print(f"  [OK] {msg}")
        else:
            print(f"  [FAIL] {msg}")
            problems += 1

    for phase in ("elem.tangent", "elem.residual"):
        node = _find(root, phase)
        check(node is not None, f"`{phase}` node present in rollup tree")
        if node is None:
            continue
        rows = node.get("elem_by_type", [])
        names = ", ".join(f"{r['classTag']}:{r['count']}" for r in rows)
        print(f"    {phase} elem_by_type ({len(rows)} class(es)): {names}")
        check(len(rows) > 0, f"`{phase}` has a non-empty elem_by_type breakdown")
        truss = next((r for r in rows if r["classTag"] == ELE_TAG_TRUSS), None)
        check(truss is not None, f"`{phase}` carries the Truss classTag ({ELE_TAG_TRUSS})")
        if truss is not None:
            c = truss["count"]
            check(c > 0 and c % N_TRUSS == 0,
                  f"`{phase}` Truss count is a positive multiple of {N_TRUSS} "
                  f"(got {c} == {N_TRUSS}x{c // N_TRUSS})")
            check(truss["wall_ns_per_elem"] >= 0.0,
                  f"`{phase}` Truss wall_ns_per_elem >= 0 (got {truss['wall_ns_per_elem']:.1f})")
            check(truss["fb_coupled"] is False,
                  f"`{phase}` Truss fb_coupled is False (displacement-based)")

    print("\n" + "=" * 52)
    if problems == 0:
        print("DEEP_SMOKE_CHECK: ALL PASS (per-element-type buckets populated)")
        return 0
    print(f"DEEP_SMOKE_CHECK: {problems} PROBLEM(S)")
    return 1


if __name__ == "__main__":
    sys.exit(main())
