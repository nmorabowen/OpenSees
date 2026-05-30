"""Multi-stage regression — check (the C1/C2 fix gate).

Asserts, for ms_ref.mpco (frozen oracle) and ms_test.ladruno (new recorder):

  1. STRUCTURE: the new .ladruno contains >= 2 MODEL_STAGE groups, and the SAME
     set of MODEL_STAGE stamps as the frozen .mpco. (Pre-fix, the .ladruno had
     only the FIRST stage -> this is the regression that proves the rebuild.)
  2. VALUE PARITY: every common nodal (node, result:comp, step) value matches to
     1e-12, across ALL stages. Reuses parity_check's reader/normalizer so the
     tolerance and result-name normalization are identical to the L1 gate.

Run with the venv python (has h5py):
    python multistage_check.py <ms_ref.mpco> <ms_test.ladruno>
"""

from __future__ import annotations

import sys

import h5py

import parity_check as pc  # ATOL/RTOL + normalize_nodal_h5 (all-stages reader)


def stage_stamps(path: str) -> set[int]:
    """Set of MODEL_STAGE[<id>] stamps present in an .mpco / .ladruno file."""
    out: set[int] = set()
    with h5py.File(path, "r") as f:
        for k in f:
            if k.startswith("MODEL_STAGE[") and k.endswith("]"):
                out.add(int(k[len("MODEL_STAGE["):-1]))
    return out


def main() -> int:
    ref_path, new_path = sys.argv[1], sys.argv[2]

    ref_stages = stage_stamps(ref_path)
    new_stages = stage_stamps(new_path)
    print(f"stages: ref={sorted(ref_stages)} new={sorted(new_stages)}")

    ok = True

    # --- 1. structure -----------------------------------------------------
    if len(new_stages) < 2:
        print(f"FAIL: new .ladruno has {len(new_stages)} MODEL_STAGE group(s); "
              f"expected >= 2 (multi-stage not written -> C1/C2 regression)")
        ok = False
    if ref_stages != new_stages:
        print(f"FAIL: MODEL_STAGE stamp sets differ "
              f"(only-ref={sorted(ref_stages - new_stages)}, "
              f"only-new={sorted(new_stages - ref_stages)})")
        ok = False
    elif len(new_stages) >= 2:
        print(f"OK structure: {len(new_stages)} matching MODEL_STAGE groups")

    # --- 2. value parity across ALL stages --------------------------------
    ref = pc.normalize_nodal_h5(ref_path)
    new = pc.normalize_nodal_h5(new_path)
    if not ref or not new:
        print(f"FAIL: empty nodal results (ref={len(ref)} new={len(new)})")
        return 1

    common = set(ref) & set(new)
    only_ref = set(ref) - set(new)
    only_new = set(new) - set(ref)
    mism = [(k, ref[k], new[k]) for k in sorted(common)
            if not (abs(ref[k] - new[k]) <= pc.ATOL + pc.RTOL * abs(ref[k]))]

    print(f"keys: ref={len(ref)} new={len(new)} common={len(common)}")
    if only_ref:
        print(f"  WARNING {len(only_ref)} keys only in ref (e.g. {sorted(only_ref)[:3]})")
    if only_new:
        print(f"  WARNING {len(only_new)} keys only in new (e.g. {sorted(only_new)[:3]})")

    if not common:
        print("FAIL: no overlapping nodal keys to compare")
        ok = False
    elif mism:
        print(f"PARITY FAIL: {len(mism)} mismatches (showing 10):")
        for (k, a, b) in mism[:10]:
            print(f"  {k}: ref={a!r} new={b!r} (d={a - b:.3e})")
        ok = False
    else:
        print(f"OK parity: {len(common)} nodal values match to {pc.ATOL:g} across all stages")

    print("MULTISTAGE_CHECK", "OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
