"""ADR-39 P2b-2 PROBE — does the MERGED P2b-1 code already handle a DEFORMABLE
master with the main-term (kn*BᵀB) tangent? (The ∂n/∂u block was deferred.)

Setup: one LadrunoBrick (unit cube), bottom face fixed, TOP face = master segment
(4 free nodes → deformable). A slave node sits just above the top face, pressed down
with P. Series response: contact penetration gap = P/kn, brick uniaxial compression
= P*L/(E*A). We check (a) static Newton CONVERGES with only kn*BᵀB, (b) the contact
force balances P, (c) the brick-top displacement matches the analytic series value.

If this converges + matches, deformable contact already works (residual exact, main
tangent enough) and the ∂n/∂u block is a convergence refinement, not a correctness gap.
"""
import os
import sys

D = "C:/Users/nmora/Github/OpenSees_Compile/OpenSees/.claude/worktrees/peaceful-bassi-ebd7ae/dist/bin"
os.add_dll_directory(D)
sys.path.insert(0, D)
import opensees as ops  # noqa: E402  (the built pyd, not the openseespy package)

KN = 1.0e6
E = 2.0e4      # brick Young's modulus (soft, so it deforms measurably)
NU = 0.0       # zero Poisson → clean uniaxial
P = 1.0e2
L = 1.0        # cube edge


def run():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    # unit cube nodes: bottom z=0 (1-4), top z=L (5-8)
    coords = {
        1: (0, 0, 0), 2: (L, 0, 0), 3: (L, L, 0), 4: (0, L, 0),
        5: (0, 0, L), 6: (L, 0, L), 7: (L, L, L), 8: (0, L, L),
    }
    for t, (x, y, z) in coords.items():
        ops.node(t, float(x), float(y), float(z))
    for t in (1, 2, 3, 4):
        ops.fix(t, 1, 1, 1)              # bottom fixed
    ops.nDMaterial("ElasticIsotropic", 1, E, NU)
    ops.element("LadrunoBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)

    # slave node just above the top face center, free in z only
    cx, cy = L / 2, L / 2
    ops.node(99, cx, cy, L - 1.0e-8)     # just penetrated into the top face
    ops.fix(99, 1, 1, 0)                  # free z

    # master = top quad (nodes 5,6,7,8, CCW seen from above → outward +z); slave = 99
    ops.contactSurface(10, "-master", 4, 5, 6, 7, 8)
    ops.contactSurface(20, "-slave", 99)
    ops.contact(1, 10, 20, KN, 0.0, 0.0, "-outward", 0.0, 0.0, 1.0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(99, 0.0, 0.0, -P)           # press the slave down into the brick

    ops.constraints("LadrunoContact")
    ops.numberer("Plain")
    ops.system("FullGeneral")
    ops.integrator("LoadControl", 1.0)
    ops.test("NormDispIncr", 1e-10, 50, 0)
    ops.algorithm("Newton")
    ops.analysis("Static")
    ok = ops.analyze(1)
    print(f"analyze returned {ok}")
    if ok != 0:
        print("DID NOT CONVERGE with kn*BᵀB main-term tangent")
        return
    z_slave = (L - 1.0e-8) + ops.nodeDisp(99, 3)
    z_top = L + ops.nodeDisp(5, 3)          # a top brick node z-displacement
    pen = z_top - z_slave                    # contact penetration (slave below top face)
    brick_comp = -(ops.nodeDisp(5, 3))       # top moved down by this
    # analytic: contact gap = P/kn ; brick uniaxial compression = P*L/(E*A), A=L^2=1
    exp_pen = P / KN
    exp_comp = P * L / (E * (L * L))
    print(f"contact penetration = {pen:.6e}  (expect P/kn = {exp_pen:.6e})")
    print(f"brick top compression = {brick_comp:.6e}  (expect PL/EA = {exp_comp:.6e})")
    print(f"slave total descent = {-ops.nodeDisp(99,3):.6e}  (expect pen+comp = {exp_pen+exp_comp:.6e})")
    pen_ok = abs(pen - exp_pen) / exp_pen < 0.02
    comp_ok = abs(brick_comp - exp_comp) / exp_comp < 0.05
    print(f"PEN match: {pen_ok} | COMP match: {comp_ok}")
    print("RESULT:", "DEFORMABLE CONTACT WORKS (main tangent converges)" if (pen_ok and comp_ok)
          else "forces/geometry off — investigate")


if __name__ == "__main__":
    run()
