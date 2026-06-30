#!/usr/bin/env python
# ADR-62 P1 regression — the kinematic mesh-tie on the REAL OpenSees solver.
#
# Proves on the shipped binary (CentralDifferenceLadruno + LadrunoProjection +
# system Diagonal) that a NON-CONFORMING mesh-tie, expressed as a weighted multi-master
# equationConstraint  u_s = Σ N_i u_{m,i}  and enforced KINEMATICALLY (no penalty), is:
#   (a) EXACT          — the tie holds to machine precision every step (zero penetration);
#   (b) Δt-NEUTRAL     — adding the tie does NOT change dt_cr (no penalty stiffness added),
#                        the whole point: a stiff penalty tie would collapse dt_cr.
#   (c) LOAD-CARRYING  — load transmits across the non-conforming interface (the bars act tied).
#
# Model: two axial (truss) bars meeting at a NON-CONFORMING interface (mirrors the P0 oracle).
#   master bar: nodes 1-2-3 (x=0,1,2)      slave bar: nodes 4-5 (x=1.3,2.3)
#   slave node 4 projects onto master facet [2,3] at ξ̄=0.3  ->  u_4 = 0.7 u_2 + 0.3 u_3.

import os, sys
try:
    sys.stdout.reconfigure(encoding="utf-8")   # Windows console: print the unicode glyphs
except Exception:
    pass

# --- bootstrap the local dist build (no oneAPI setvars needed to RUN) ---
DIST = r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin"
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
    sys.path.insert(0, DIST)
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
import opensees as ops

N2, N3 = 0.7, 0.3        # master facet shape weights at the slave's projection (ξ̄=0.3)
E, A, RHO = 1000.0, 1.0, 1.0   # element density -> element mass (the dt_cr kernel reads ELEMENT mass)


def build(tied: bool, with_load: bool):
    ops.wipe()
    ops.model("basic", "-ndm", 1, "-ndf", 1)         # pure 1D axial chain (clean dt_cr)
    for tag, x in [(1, 0.0), (2, 1.0), (3, 2.0), (4, 1.3), (5, 2.3)]:
        ops.node(tag, x)
    ops.fix(1, 1)                                     # anchor the master base
    ops.uniaxialMaterial("Elastic", 1, E)
    ops.element("Truss", 1, 1, 2, A, 1, "-rho", RHO)    # master bar (element mass via -rho)
    ops.element("Truss", 2, 2, 3, A, 1, "-rho", RHO)
    ops.element("Truss", 3, 4, 5, A, 1, "-rho", RHO)    # slave bar
    if tied:
        # kinematic tie u_4 = 0.7 u_2 + 0.3 u_3  ->  1*u_4 + (-0.7)*u_2 + (-0.3)*u_3 = 0
        ops.equationConstraint(4, 1, 1.0, 2, 1, -N2, 3, 1, -N3)
        ops.constraints("LadrunoProjection")
    else:
        ops.constraints("Plain")
    if with_load:
        ops.timeSeries("Constant", 1); ops.pattern("Plain", 1, 1)
        ops.load(3, 1.0)                   # pull the master free end
    ops.numberer("Plain"); ops.system("Diagonal")
    ops.test("NormDispIncr", 1e-12, 10); ops.algorithm("Linear")
    ops.integrator("CentralDifferenceLadruno", "-cfl")   # -cfl: populate the dt_cr query
    ops.analysis("Transient")


def main():
    print("ADR-62 P1 — kinematic mesh-tie on the REAL OpenSees solver\n")
    ok = True

    # (b) Δt-NEUTRALITY — dt_cr with the tie vs without it.
    # dt_cr populates at domainChanged, so take ONE safe step (dt0 << dt_cr~0.03) then query.
    dt0 = 1.0e-3
    build(tied=True, with_load=False);  ops.analyze(1, dt0);  dt_tied = ops.criticalTimeStep()
    build(tied=False, with_load=False); ops.analyze(1, dt0);  dt_free = ops.criticalTimeStep()
    print(f"(b) Δt-NEUTRAL:  dt_cr tied = {dt_tied:.6e}   untied = {dt_free:.6e}")
    print(f"                 ratio tied/untied = {dt_tied/dt_free:.6f}  "
          f"(kinematic tie adds NO stiffness -> dt_cr unchanged; a stiff penalty would collapse it)")
    ok &= abs(dt_tied/dt_free - 1.0) < 1e-9

    # (a) EXACT + (c) LOAD-CARRYING — run the tied bar under load, check the tie every step.
    build(tied=True, with_load=True)
    dt = 0.4 * dt_tied
    tie_err = 0.0
    for _ in range(400):
        assert ops.analyze(1, dt) == 0, "tied explicit analyze failed"
        u2, u3, u4 = ops.nodeDisp(2, 1), ops.nodeDisp(3, 1), ops.nodeDisp(4, 1)
        tie_err = max(tie_err, abs(u4 - (N2 * u2 + N3 * u3)))
    u5 = ops.nodeDisp(5, 1)
    print(f"\n(a) EXACT:       max |u_4 - (0.7 u_2 + 0.3 u_3)| = {tie_err:.3e}   over 400 steps")
    print(f"(c) LOAD-CARRYING: master free end u_3 = {u3:.4e}, slave end u_5 = {u5:.4e}  "
          f"(slave bar moves -> load crossed the non-conforming tie)")
    ok &= tie_err < 1e-10
    ok &= abs(u5) > 1e-9          # the slave actually responded (tie transmits)

    print("\n" + ("ALL PASS" if ok else "FAILURES PRESENT"))
    print("KEY: weighted multi-master tie enforced KINEMATICALLY on the shipped CentralDifference"
          "Ladruno+LadrunoProjection — EXACT, dt_cr-neutral, load-carrying. No penalty, no added mass.")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
