"""Parity gate — value diff (L1 of the test plan).

Reads ref.mpco (frozen oracle) and test.ladruno (new recorder), normalizes the
nodal results of BOTH to a canonical {(node, result:comp, step): value} dict, and
asserts every value in the intersection matches to 1e-12.

Both formats put nodal results at MODEL_STAGE[*]/RESULTS/ON_NODES/<name>/{ID,
COMPONENTS attr, DATA/STEP_k}, so one generic h5py reader handles both. Result-
group names are normalized (lowercase, alnum-only) so REACTION_FORCE ~ ReactionForce.

Run with the venv python (has h5py):
    python parity_check.py <ref.mpco> <test.ladruno>
"""

from __future__ import annotations

import sys

import h5py
import numpy as np

import ladruno_format as lf  # for the MPCO-style attr decoder (_attr)

ATOL = 1e-12
RTOL = 1e-9


def _norm(name: str) -> str:
    return "".join(c for c in name.lower() if c.isalnum())


def normalize_nodal_h5(path: str) -> dict[tuple, float]:
    """{(node_tag, 'normresult:comp', step): value} over all ON_NODES results,
    for either a legacy .mpco or a .ladruno file."""
    out: dict[tuple, float] = {}
    with h5py.File(path, "r") as f:
        for stage in [k for k in f if k.startswith("MODEL_STAGE")]:
            base = f.get(f"{stage}/RESULTS/ON_NODES")
            if base is None:
                continue
            for rname in base:
                grp = base[rname]
                if "ID" not in grp or "DATA" not in grp:
                    continue  # skip modal/irregular groups for this first cut
                ids = np.asarray(grp["ID"][...]).reshape(-1)  # ID may be [nN] or [nN×1]
                comps = [c for c in lf._attr(grp, "COMPONENTS").split(",") if c]
                rkey = _norm(rname)
                # handles both the chunked .ladruno DATA[T×nIds×nComp] layout and
                # the per-step .mpco DATA/STEP_<k> layout.
                for step, arr in lf.iter_step_slices(grp):
                    for r, tag in enumerate(ids):
                        for c, cn in enumerate(comps):
                            if c < arr.shape[1]:
                                out[(int(tag), f"{rkey}:{cn}", step)] = float(arr[r, c])
    return out


def main() -> int:
    ref_path, new_path = sys.argv[1], sys.argv[2]
    ref = normalize_nodal_h5(ref_path)
    new = normalize_nodal_h5(new_path)

    if not ref:
        print("FAIL: no nodal results read from the reference .mpco")
        return 1
    if not new:
        print("FAIL: no nodal results read from the new .ladruno")
        return 1

    common = set(ref) & set(new)
    only_ref = set(ref) - set(new)
    only_new = set(new) - set(ref)

    mism = []
    for k in sorted(common):
        a, b = ref[k], new[k]
        if not (abs(a - b) <= ATOL + RTOL * abs(a)):
            mism.append((k, a, b))

    print(f"keys: ref={len(ref)} new={len(new)} common={len(common)}")
    if only_ref:
        print(f"  WARNING {len(only_ref)} keys only in ref (e.g. {sorted(only_ref)[:3]})")
    if only_new:
        print(f"  WARNING {len(only_new)} keys only in new (e.g. {sorted(only_new)[:3]})")

    if mism:
        print(f"PARITY FAIL: {len(mism)} mismatches (showing 10):")
        for (k, a, b) in mism[:10]:
            print(f"  {k}: ref={a!r} new={b!r} (d={a - b:.3e})")
        return 1
    if not common:
        print("PARITY FAIL: no overlapping keys to compare (naming mismatch?)")
        return 1

    print(f"PARITY OK — {len(common)} nodal values match to {ATOL:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
