"""RC squat shear-wall panel harness — Phase 2b.2c.4 (Tran-Wallace pinching validation, WIP).

This is the STARTING POINT for the one remaining Phase-2b.2c item: a meshed, reinforced,
non-homogeneous RC wall under cyclic in-plane shear, where principal-direction rotation
produces the pinched hysteresis *waist* that a homogeneous material point structurally cannot
(see ADR 19 -- "PINCHING is a PANEL phenomenon"). It is NOT a committed pytest gate: it lives
under _testbed and its filename does not start with `test_`, so pytest does not collect it.

Status (run locally 2026-06-17, OpenSeesPy build of this branch):
  * The model ASSEMBLES and RUNS: a structured NX x NY grid of ASDShellQ4 on a LayeredShell
    section = 4 LadrunoRCConcrete concrete layers (full cyclic stack: -beta -lublinerReduced
    -interlock -cyclic -implex) + 2 smeared PlateRebar(Steel02) layers (vertical + horizontal
    web steel). No gmsh -- a pure-ops grid, so it is environment-portable.
  * It produces a real cyclic shear response WITH HYSTERETIC DISSIPATION on the first drift
    cycle (V ~ +/-146 kN at 0.3 mm; closed-loop area > 0).
  * It WALLS on convergence at larger drift (~0.6 mm onward) -- the classic cyclic-softening
    RC-wall barrier. A robust multi-cycle solve to the drifts where the pinched waist is
    pronounced needs dedicated convergence engineering, which is exactly why this is the
    research-grade Zone-B validation the ADR flags (not a clean material slice).

To finish the validation (the deferred work):
  1. Robust solver: arc-length / LadrunoIndirectControl follower (the snap-back driver built
     for exactly this), or dynamic relaxation; finer DisplacementControl substeps; possibly
     element-by-element step-cut on the IMPL-EX error monitor.
  2. Calibrate geometry / reinforcement / backbone to a named specimen (Tran-Wallace RW-A20-P10
     or a PEER squat-wall) and assert pinching SHAPE + cumulative hysteretic ENERGY vs the
     measured loops -- the ADR's primary squat-wall gate.
  3. (Optional) re-mesh via gmsh/apeGmsh for a graded mesh + boundary elements; mark zone_b.

The cyclic MATERIAL physics this wall exercises is already complete + verified at the
material-point and single-shell-element scale (Phases 1..2b.2c.3): compression softening,
the v_ci,max interlock bound, cyclic friction-slip, X-cracking + wear, the -shearRetention
curves, IMPL-EX robustness, rigid-rotation objectivity, and unilateral crack closure.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from _testbed import ops

# SI units N, mm, MPa
E, NU = 30000.0, 0.2
CE = [0.0, 0.0007, 0.0020, 0.0100]
CS = [0.0, 24.0,   30.0,   5.0]
CD = [0.0, 0.0,    0.25,   1.0 - 5.0 / 45.0]
TE = [0.0, 0.0001, 0.0010]
TS = [0.0, 3.0,    0.5]
TD = [0.0, 0.0,    1.0 - 0.5 / 5.0]

L, H, TH = 1000.0, 1000.0, 120.0     # squat wall 1 m x 1 m, 120 mm thick (aspect ~1)
NX, NY = 4, 4
Fy, Es, BSTEEL = 420.0, 200000.0, 0.01
RHO_H, RHO_V = 0.006, 0.010          # web steel ratios


def build(implex=True):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nid = {}
    k = 1
    for j in range(NY + 1):
        for i in range(NX + 1):
            nid[(i, j)] = k
            ops.node(k, i * L / NX, j * H / NY, 0.0)
            k += 1
    cargs = ["LadrunoRCConcrete", 1, E, NU,
             "-Ce", *CE, "-Cs", *CS, "-Cd", *CD, "-Te", *TE, "-Ts", *TS, "-Td", *TD,
             "-beta", "-lublinerReduced", "-interlock", "-cyclic",
             "-crackSpacing", 80.0, "-agg", 16.0]
    if implex:
        cargs += ["-implex"]
    ops.nDMaterial(*cargs)
    ops.uniaxialMaterial("Steel02", 11, Fy, Es, BSTEEL, 18.0, 0.925, 0.15)
    ops.nDMaterial("PlateRebar", 21, 11, 90.0)   # vertical steel
    ops.nDMaterial("PlateRebar", 22, 11, 0.0)    # horizontal steel
    tc = TH * 0.94
    ops.section("LayeredShell", 1, 6,
                1, tc / 4, 1, tc / 4, 1, tc / 4, 1, tc / 4,
                21, RHO_V * TH, 22, RHO_H * TH)
    eid = 1
    for j in range(NY):
        for i in range(NX):
            ops.element("ASDShellQ4", eid,
                        nid[(i, j)], nid[(i + 1, j)], nid[(i + 1, j + 1)], nid[(i, j + 1)], 1)
            eid += 1
    base = [nid[(i, 0)] for i in range(NX + 1)]
    top = [nid[(i, NY)] for i in range(NX + 1)]
    for n in base:
        ops.fix(n, 1, 1, 1, 1, 1, 1)
    for (i, j), n in nid.items():
        if j > 0:
            ops.fix(n, 0, 0, 1, 1, 1, 0)        # keep it an in-plane 2D problem
    return nid, base, top


def run(amps=(0.3, 0.3, 0.6, 0.6, 1.2, 1.2, 2.5, 2.5), nstep=80):
    nid, base, top = build()
    ctrl = top[NX // 2]
    for n in top:
        if n != ctrl:
            ops.equalDOF(ctrl, n, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(ctrl, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-4, 200, 0)
    ops.algorithm("KrylovNewton")
    ops.analysis("Static")

    out = []
    cur = 0.0
    for a in amps:
        for target in (a, -a):
            d = (target - cur) / nstep
            ops.integrator("DisplacementControl", ctrl, 1, d)
            for _ in range(nstep):
                ok = ops.analyze(1)
                if ok != 0:
                    for alg in ("NewtonLineSearch", "ModifiedNewton", "Broyden"):
                        ops.algorithm(alg)
                        ok = ops.analyze(1)
                        ops.algorithm("KrylovNewton")
                        if ok == 0:
                            break
                    if ok != 0:
                        print(f"  NON-CONVERGE at target {target} mm (cur {cur:.3f})")
                        return out
                ops.reactions()
                out.append((ops.nodeDisp(ctrl, 1), -sum(ops.nodeReaction(n, 1) for n in base)))
            cur = target
        print(f"  amp {a} mm done, V={out[-1][1]:.1f} N")
    return out


if __name__ == "__main__":
    out = run()
    print(f"steps recorded: {len(out)}")
    if out:
        area = sum(0.5 * (out[i][1] + out[i - 1][1]) * (out[i][0] - out[i - 1][0])
                   for i in range(1, len(out)))
        vmax = max(abs(v) for _, v in out)
        dmax = max(abs(d) for d, _ in out)
        print(f"Vmax={vmax:.1f} N  dmax={dmax:.3f} mm  hysteretic area={area:.1f} N.mm")
