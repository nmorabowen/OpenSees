"""Flag-order check — every option on the Tcl recorder line was honored.

Companion to flag_order_model.tcl, which puts `-G energy 1` in a NON-final
position (`... -G energy 1 -T nsteps 1`). The parse bug this guards against
desynced the option scan by one token in the classic Tcl exe, so a "passing"
parse could still silently drop a trailing option. Assert all three channels
landed: -N displacement (ON_NODES), -G energy (ON_DOMAIN), -G region 1
(ON_REGIONS).

Usage: python flag_order_check.py <flag_order.ladruno>
"""

from __future__ import annotations

import sys

import h5py


def _has_result(f, family, name):
    for stage in sorted(k for k in f if k.startswith("MODEL_STAGE")):
        grp = f.get(f"{stage}/RESULTS/{family}/{name}")
        if grp is not None and "DATA" in grp:
            return True
    return False


def main():
    path = sys.argv[1]
    ok = True
    with h5py.File(path, "r") as f:
        for family, name, flag in (
            ("ON_NODES", "DISPLACEMENT", "-N displacement"),
            ("ON_DOMAIN", "energyBalance", "-G energy"),
            ("ON_REGIONS", "energyBalance", "-G energy 1"),
        ):
            if _has_result(f, family, name):
                print(f"  OK   {family}/{name}  ({flag})")
            else:
                print(f"  FAIL {family}/{name} missing  ({flag})")
                ok = False
    if not ok:
        print("FLAG_ORDER CHECK FAILED")
        return 1
    print("FLAG_ORDER CHECK OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
