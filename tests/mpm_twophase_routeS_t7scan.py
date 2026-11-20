"""T7 diagnostic (WP-232l): the empirical stability factor, in multiples of the drained pencil, of
route S on T6's Terzaghi column; 2000 steps per factor. `mpm debug routeSLimit off` lifts the
run-time limit and its cap (debug only) so the scan measures above them.

Usage: python tests/mpm_twophase_routeS_t7scan.py <dist\\bin> <basis> [off|auto] [up|twin] [Kf]
"""
import os
import sys
from pathlib import Path

B = sys.argv[1]
if hasattr(os, "add_dll_directory"):
    os.add_dll_directory(B)
sys.path.insert(0, B)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import opensees as ops  # noqa: E402
import ladruno_mpm_twophase_cases_harness as H  # noqa: E402

basis = sys.argv[2]
stab = "off" if (len(sys.argv) > 3 and sys.argv[3] == "off") else ("auto", 0.25)
twin = len(sys.argv) > 4 and sys.argv[4] == "twin"


def spec():
    sp = H.t6_terzaghi_spec(basis, stab)
    if len(sys.argv) > 5:
        f = list(sp["fluid"])
        f[0] = float(sys.argv[5])
        sp["fluid"] = tuple(f)
    if twin:
        sp.pop("fluid")
    return sp


H.build(ops, spec())
dtd = H.as_list(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
print(ops.ladrunoBuild(), basis, "twin" if twin else "up", "dt_d", dtd, flush=True)
for fac in (0.3, 0.5, 0.6, 0.7, 0.8, 0.9):
    H.build(ops, spec())
    ops.mpm("debug", "routeSLimit", "off")
    dt = fac * dtd
    vmax = 0.0
    st = "stable"
    for s in range(2000):
        try:
            ops.mpm("analyze", 1, dt)
        except Exception:
            st = "refused at step %d" % s
            break
        vm = max(abs(x) for x in H.pdata(ops, "v"))
        vmax = max(vmax, vm)
        if vm > 1.0:
            st = "blew up at step %d" % (s + 1)
            break
    print("  factor %.2f dt %.5f: %s, vmax %.3e" % (fac, dt, st, vmax), flush=True)
