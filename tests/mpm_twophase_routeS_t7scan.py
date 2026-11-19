"""T7 diagnostic: the empirical stability factor (x drained pencil) of route S on T6's Terzaghi column
after interface 7c; 2000 steps; argv: basis [stab] [twin]"""
import os
import sys
B = r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\compassionate-mccarthy-18cf89\dist\bin"
os.add_dll_directory(B)
sys.path.insert(0, B)
sys.path.insert(0, r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\agent-adf46e5c8275479c0\tests")
import opensees as ops  # noqa: E402
import ladruno_mpm_twophase_cases_harness as H  # noqa: E402

basis = sys.argv[1]
stab = "off" if (len(sys.argv) > 2 and sys.argv[2] == "off") else ("auto", 0.25)
twin = len(sys.argv) > 3 and sys.argv[3] == "twin"


def spec(s=None):
    sp = H.t6_terzaghi_spec(basis, stab, safety=s)
    if len(sys.argv) > 4:
        f = list(sp["fluid"]); f[0] = float(sys.argv[4]); sp["fluid"] = tuple(f)
    if twin:
        sp.pop("fluid")
    return sp


H.build(ops, spec())
dtd = H.as_list(ops.mpm("criticalTimeStep", "drainedPencil"))[0]
print(basis, "twin" if twin else "up", "dt_d", dtd, flush=True)
for fac in (0.3, 0.5, 0.6, 0.7, 0.8, 0.9):
    H.build(ops, spec(1.0))
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
