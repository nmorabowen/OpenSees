"""Element-result parity — value diff (L1, element channel).

Reads eref.mpco (frozen oracle) and etest.ladruno (new recorder) and asserts
every ELEMENT result value matches to 1e-12.

Both formats store element results at the SAME structural location:
    MODEL_STAGE[*]/RESULTS/ON_ELEMENTS/<result>/<bucket>/{ID[nElem,1], DATA/STEP_k[nElem,nCol]}
(.mpco adds a META/ subgroup, .ladruno a COLUMN_MAP/ subgroup — different metadata,
identical DATA). Because ladruno's ElementResultSource::evaluate() is a
byte-faithful port of the frozen recordResultsOnElements packing (same per-element
Response getData() order), column j of element tag t at step s must be identical.

We therefore key each value by (element_tag, result, column_index, step) — read off
each bucket's ID dataset — and never depend on the bucket NAME (the frozen bracket
carries an extra header-index field) or on decoding component names.

Run with the venv python (has h5py):
    python element_parity_check.py <eref.mpco> <etest.ladruno>
"""

from __future__ import annotations

import sys

import h5py
import numpy as np

import ladruno_format as lf  # MPCO-style attr decoder (_attr squeezes size-1 arrays)

ATOL = 1e-12
RTOL = 1e-9


def _norm(name: str) -> str:
    return "".join(c for c in name.lower() if c.isalnum())


def normalize_element_h5(path: str) -> dict[tuple, float]:
    """{(elem_tag, 'normresult', col_index, step): value} over all ON_ELEMENTS
    buckets, for either a legacy .mpco or a .ladruno file."""
    out: dict[tuple, float] = {}
    with h5py.File(path, "r") as f:
        for stage in [k for k in f if k.startswith("MODEL_STAGE")]:
            base = f.get(f"{stage}/RESULTS/ON_ELEMENTS")
            if base is None:
                continue
            for rname in base:                       # result display, e.g. "stress"
                rgrp = base[rname]
                if not isinstance(rgrp, h5py.Group):
                    continue
                rkey = _norm(rname)
                for bname in rgrp:                   # bucket, e.g. "31-FourNodeQuad[201:0:0]"
                    bgrp = rgrp[bname]
                    if not isinstance(bgrp, h5py.Group):
                        continue
                    if "ID" not in bgrp or "DATA" not in bgrp:
                        continue
                    ids = np.asarray(bgrp["ID"][...]).reshape(-1)
                    # chunked .ladruno DATA[T×nElem×nCol] or per-step .mpco DATA/STEP_<k>
                    for step, arr in lf.iter_step_slices(bgrp):
                        for r, tag in enumerate(ids):
                            for j in range(arr.shape[1]):
                                out[(int(tag), rkey, j, step)] = float(arr[r, j])
    return out


def main() -> int:
    ref_path, new_path = sys.argv[1], sys.argv[2]
    ref = normalize_element_h5(ref_path)
    new = normalize_element_h5(new_path)

    if not ref:
        print("FAIL: no element results read from the reference .mpco")
        return 1
    if not new:
        print("FAIL: no element results read from the new .ladruno")
        return 1

    common = set(ref) & set(new)
    only_ref = set(ref) - set(new)
    only_new = set(new) - set(ref)

    mism = [(k, ref[k], new[k]) for k in sorted(common)
            if not (abs(ref[k] - new[k]) <= ATOL + RTOL * abs(ref[k]))]

    print(f"keys: ref={len(ref)} new={len(new)} common={len(common)}")
    if only_ref:
        print(f"  WARNING {len(only_ref)} keys only in ref (e.g. {sorted(only_ref)[:3]})")
    if only_new:
        print(f"  WARNING {len(only_new)} keys only in new (e.g. {sorted(only_new)[:3]})")

    # Guard against a trivially-passing all-zero comparison.
    nz = sum(1 for k in common if abs(ref[k]) > 0.0)
    print(f"  nonzero ref values among common: {nz}/{len(common)}")

    if mism:
        print(f"ELEMENT PARITY FAIL: {len(mism)} mismatches (showing 10):")
        for (k, a, b) in mism[:10]:
            print(f"  {k}: ref={a!r} new={b!r} (d={a - b:.3e})")
        return 1
    if not common:
        print("ELEMENT PARITY FAIL: no overlapping keys (structure mismatch?)")
        return 1
    if nz == 0:
        print("ELEMENT PARITY FAIL: all compared values are zero (model produced no stress?)")
        return 1

    print(f"ELEMENT PARITY OK — {len(common)} element values match to {ATOL:g} ({nz} nonzero)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
