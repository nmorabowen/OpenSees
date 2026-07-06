"""Sanity gate on the benchmark input files themselves (adversarial-review
NOTE #11): the fiber-column arm is the PESSIMISTIC case for a transform codec
only if it actually went inelastic — if steel never yielded, the aggregate is
carried by the two smooth models and the ADR gate is vacuous.

  python check_inputs.py <bench_b_fiber_column.ladruno> [others...]

Prints per-file dataset inventory; for the fiber file asserts
max |section.fiber.strain| > steel yield strain (420/200000 = 2.1e-3).
"""
from __future__ import annotations

import os
import sys

import h5py
import numpy as np

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

YIELD_STRAIN = 420.0 / 200000.0


def main():
    rc = 0
    for path in sys.argv[1:]:
        print(f"\n== {os.path.basename(path)} ==")
        fiber_strain_max = None
        with h5py.File(path, "r") as f:
            def visit(name, obj):
                nonlocal fiber_strain_max
                if (isinstance(obj, h5py.Dataset) and name.endswith("/DATA")
                        and obj.ndim == 3):
                    arr = obj[...]
                    amax = float(np.abs(arr).max())
                    print(f"  {name}: {obj.shape} raw {arr.nbytes/1e6:8.2f} MB"
                          f"  max|v| {amax:.3e}")
                    if "fiber.strain" in name.lower() or "fiberstrain" in name.lower():
                        fiber_strain_max = max(fiber_strain_max or 0.0, amax)
            f.visititems(visit)
        if "fiber" in os.path.basename(path):
            if fiber_strain_max is None:
                print("  !! no fiber-strain dataset found")
                rc = 1
            elif fiber_strain_max <= YIELD_STRAIN:
                print(f"  !! NOT YIELDED: max fiber strain {fiber_strain_max:.3e}"
                      f" <= eps_y {YIELD_STRAIN:.3e} — pessimistic arm is vacuous")
                rc = 1
            else:
                print(f"  yielded OK: max fiber strain {fiber_strain_max:.3e}"
                      f" > eps_y {YIELD_STRAIN:.3e}"
                      f" (ductility ~{fiber_strain_max/YIELD_STRAIN:.1f}x)")
    sys.exit(rc)


if __name__ == "__main__":
    main()
