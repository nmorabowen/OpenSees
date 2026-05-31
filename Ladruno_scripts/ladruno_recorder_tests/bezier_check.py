"""BezierTri6 recorder gate -- checker.

Validates bezier_test.ladruno against schema v1 AND asserts the BezierTri6
self-description landed correctly (the whole point of the inverted-dispatch
contract):

  * schema validate() returns NO problems  (this alone proves GP_PARAM is
    [num_gp x len(ORDER)] and GP_WEIGHT len == num_gp -- the old recorder,
    which didn't register BezierTri6, produced a Custom group with no
    QUADRATURE and would fail here)
  * FAMILY == "bernstein"   (NOT the lagrange guess)
  * TOPOLOGY == "tri", PARAM_DOMAIN == "bary", RATIONAL == 0, NUM_CTRL == 6
  * ORDER == [2, 2]         (total degree, one entry per bary direction)
  * GEOMETRY attr == 101    (Triangle_6N)
  * GP_PARAM == the element's exact 3-pt area-coord rule, GP_WEIGHT == 1/6 each
  * CONNECTIVITY is [2 x 7]  (2 elements, 1 tag + 6 ctrl nodes)

Then a value-parity pass on the element-stress channel vs the frozen .mpco
oracle (both share the element engine, so per-(elem,col,step) stress must
match to 1e-12).

Run with the venv python (has h5py):
    python bezier_check.py [out_dir]
"""

from __future__ import annotations

import os
import sys

import h5py
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ladruno_format import LadrunoReader, iter_step_slices, validate  # noqa: E402

OUT = sys.argv[1] if len(sys.argv) > 1 else "."
LAD = os.path.join(OUT, "bezier_test.ladruno")
REF = os.path.join(OUT, "bezier_ref.mpco")

TOL = 1e-12
fails: list[str] = []


def check(cond: bool, msg: str) -> None:
    print(("  ok  " if cond else " FAIL ") + msg)
    if not cond:
        fails.append(msg)


# ---------------------------------------------------------------- schema
problems = validate(LAD)
check(not problems, f"schema validate(): {len(problems)} problem(s)")
for p in problems:
    print("        -", p)

# ---------------------------------------------------------------- basis
with LadrunoReader(LAD) as r:
    stages = r.stages()
    check(len(stages) >= 1, f"at least one MODEL_STAGE ({len(stages)})")
    stage = stages[0]
    groups = r.element_groups(stage)

    bez = {n: g for n, g in groups.items() if "BezierTri6" in n}
    check(len(bez) == 1, f"exactly one BezierTri6 element group ({len(bez)})")

    if bez:
        name, g = next(iter(bez.items()))
        print("        group:", name)
        check(g["family"] == "bernstein", f'FAMILY == bernstein (got {g["family"]!r})')
        check(g["topology"] == "tri", f'TOPOLOGY == tri (got {g["topology"]!r})')
        check(g["param_domain"] == "bary", f'PARAM_DOMAIN == bary (got {g["param_domain"]!r})')
        check(g["rational"] == 0, f'RATIONAL == 0 (got {g["rational"]})')
        check(g["num_ctrl"] == 6, f'NUM_CTRL == 6 (got {g["num_ctrl"]})')
        check(g["num_gp"] == 3, f'NUM_GP == 3 (got {g["num_gp"]})')
        check(tuple(g["order"]) == (2, 2), f'ORDER == (2,2) (got {tuple(g["order"])})')
        check(
            g["connectivity"].shape == (2, 7),
            f'CONNECTIVITY [2 x 7] (got {g["connectivity"].shape})',
        )

        # GEOMETRY attr (legacy enum retained) == Triangle_6N == 101
        with h5py.File(LAD, "r") as f:
            grp = f[f"{stage}/MODEL/ELEMENTS/{name}"]
            geom = int(np.atleast_1d(grp.attrs["GEOMETRY"]).reshape(-1)[0])
            cdim = int(np.atleast_1d(grp.attrs["CUSTOM_INTEGRATION_RULE_DIMENSION"]).reshape(-1)[0])
        check(geom == 101, f"GEOMETRY == 101 (Triangle_6N) (got {geom})")
        check(cdim == 2, f"CUSTOM_INTEGRATION_RULE_DIMENSION == 2 (got {cdim})")

        # exact 3-pt area-coord rule + weights (BezierTri6.cpp GP3_xi / GP3_w)
        exp_param = np.array(
            [[1 / 6, 1 / 6], [2 / 3, 1 / 6], [1 / 6, 2 / 3]], dtype=float
        )
        exp_w = np.full(3, 1 / 6, dtype=float)
        check(
            g["gp_param"].shape == (3, 2)
            and np.allclose(g["gp_param"], exp_param, atol=TOL),
            f'GP_PARAM == 3-pt area-coord rule (got shape {g["gp_param"].shape})',
        )
        check(
            g["gp_weight"].shape == (3,) and np.allclose(g["gp_weight"], exp_w, atol=TOL),
            f'GP_WEIGHT == [1/6,1/6,1/6] (got {np.asarray(g["gp_weight"]).ravel()})',
        )


# ------------------------------------------------------- element-stress parity
def read_elem_stress(path: str) -> dict[tuple, float]:
    """{(elem_tag, col, step): value} over ON_ELEMENTS/stress/* buckets.

    Layout is shared by .mpco and .ladruno: each bucket has ID [nElem] and
    DATA/STEP_k [nElem x ncols]. Column order is identical (same element
    engine getData), so we diff by raw column index -- no META/COLUMN_MAP."""
    out: dict[tuple, float] = {}
    with h5py.File(path, "r") as f:
        stages = [k for k in f if k.startswith("MODEL_STAGE")]
        for st in stages:
            base = f.get(f"{st}/RESULTS/ON_ELEMENTS/stress")
            if base is None:
                continue
            for bname in base:
                bucket = base[bname]
                if "ID" not in bucket or "DATA" not in bucket:
                    continue
                ids = bucket["ID"][...].reshape(-1)
                # chunked .ladruno DATA[T×nElem×nCol] or per-step .mpco DATA/STEP_<k>
                for k, arr in iter_step_slices(bucket):
                    for row, tag in enumerate(ids):
                        for col in range(arr.shape[1]):
                            out[(int(tag), col, k)] = float(arr[row, col])
    return out


if os.path.exists(REF):
    a = read_elem_stress(REF)
    b = read_elem_stress(LAD)
    common = set(a) & set(b)
    check(len(common) > 0, f"element-stress keys present in both ({len(common)})")
    nonzero = sum(1 for k in common if abs(a[k]) > 1e-9)
    worst = max((abs(a[k] - b[k]) for k in common), default=0.0)
    check(worst <= TOL, f"element-stress parity vs frozen <= {TOL} (worst {worst:.2e})")
    check(nonzero > 0, f"element stresses are non-trivial ({nonzero} nonzero)")
    print(f"        compared {len(common)} (elem,col,step) keys, {nonzero} nonzero")
else:
    print("  ..  frozen oracle bezier_ref.mpco absent; skipped element parity")

# ---------------------------------------------------------------- verdict
print()
if fails:
    print(f"BEZIER_CHECK FAILED ({len(fails)} problem(s))")
    sys.exit(1)
print("BEZIER_CHECK PASSED")
