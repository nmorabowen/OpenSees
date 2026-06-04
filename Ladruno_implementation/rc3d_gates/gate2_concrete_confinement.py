"""RC-3D Pre-build GATE 2 — does ASDConcrete3D produce EMERGENT confinement?

The embedded-reinforcement ADR (20_ladruno_embedded_reinforcement_adr.md) assumes
confined column/joint strength self-emerges from the 3D triaxial concrete model
(no Mander-style fc inflation). The adversarial review flagged this as a BLOCKER,
claiming ASDConcrete3D cannot confine. This gate settles it empirically.

Test: a single 8-node brick under CONSTANT lateral confining pressure p (held via a
load pattern frozen with loadConst) + axial displacement control. Symmetry rollers on
the three x=0/y=0/z=0 faces; +x,+y faces loaded to -p and kept planar via equalDOF;
+z face pushed down. Peak axial stress = confined strength fcc(p). Compare to Mander
(unconfined-calibrated) and Richart fcc = fc + 4.1 p.

ACCEPTANCE (gate PASSES if all hold):
  * unconfined (p=0) peak == fc within 1%   (backbone is correctly calibrated)
  * fcc(p) within ~10% of Mander for p/fc in [0, 0.20]
  * peak strain eps_pk increases monotonically with p (ductility gain)

RESULT 2026-06-03 (build f6700775, stdBrick + ASDConcrete3D):
  p/fc   fcc_sim  Mander  sim/Mander  eps_pk
  0.00    30.00   30.00     1.00      0.0020
  0.05    37.24   39.30     0.95      0.0025
  0.10    44.50   46.95     0.95      0.0029
  0.15    51.71   53.47     0.97      0.0034
  0.20    58.96   59.16     1.00      0.0039
  => GATE 2 PASS. Emergent confinement is real (Lubliner Kc surface); within 5% of
     Mander across the range; ductility gain present. The "no dilation knob" caveat
     stands but does NOT prevent emergent confinement. NB: needed KrylovNewton +
     recursive step-halving to converge at high p (plain Newton diverged) — evidence
     for the adaptive-stepping tooling need.

Run:  pythoncore-3.12-64/python.exe gate2_concrete_confinement.py
"""
import os
import sys
import math


def _bootstrap_opensees():
    try:
        import opensees as ops  # local build already on path / CI copy
        return ops
    except ModuleNotFoundError:
        pass
    dist = os.environ.get(
        "LADRUNO_DIST",
        r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\dist\bin",
    )
    if os.path.isdir(dist):
        os.add_dll_directory(dist)
        sys.path.insert(0, dist)
    import opensees as ops
    return ops


ops = _bootstrap_opensees()

# --- material: a calibrated concrete backbone (fc=30, E=30000) --------------
E, NU, FC, FT = 30000.0, 0.2, 30.0, 3.0
# compression backbone: first point = elastic limit (0.4 fc, d=0); peak -fc at -0.002
CE = [-0.4 * FC / E, -0.002, -0.005, -0.015]
CS = [-0.4 * FC, -FC, -0.5 * FC, -0.2 * FC]
CD = [1.0 - abs(CS[i]) / (E * abs(CE[i])) for i in range(len(CE))]
TE = [FT / E, 1.0e-3]
TS = [FT, 0.3]
TD = [0.0, 1.0 - 0.3 / (E * 1.0e-3)]

_CUBE = {1: (0, 0, 0), 2: (1, 0, 0), 3: (1, 1, 0), 4: (0, 1, 0),
         5: (0, 0, 1), 6: (1, 0, 1), 7: (1, 1, 1), 8: (0, 1, 1)}


def mander(p):
    return FC * (-1.254 + 2.254 * math.sqrt(1 + 7.94 * p / FC) - 2 * p / FC)


def richart(p):
    return FC + 4.1 * p


def _adapt(node, dof, target, du0):
    """Drive a control DOF to `target` with recursive step-halving on divergence."""
    cur, du = 0.0, du0
    while abs(cur - target) > 1.0e-9:
        step = du if abs(target - cur) >= abs(du) else (target - cur)
        ops.integrator("DisplacementControl", node, dof, step)
        if ops.analyze(1) == 0:
            cur += step
            yield cur
        else:
            du *= 0.5
            if abs(du) < abs(du0) * 1.0e-3:
                return


def run_conf(p, umax=-12e-3, n0=1500):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for t, c in _CUBE.items():
        ops.node(t, *c)
    ops.nDMaterial("ASDConcrete3D", 1, E, NU,
                   "-Te", *TE, "-Ts", *TS, "-Td", *TD,
                   "-Ce", *CE, "-Cs", *CS, "-Cd", *CD)
    ops.element("stdBrick", 1, 1, 2, 3, 4, 5, 6, 7, 8, 1)
    for n in (1, 4, 5, 8):
        ops.fix(n, 1, 0, 0)          # x=0 -> ux=0
    for n in (1, 2, 5, 6):
        ops.fix(n, 0, 1, 0)          # y=0 -> uy=0
    for n in (1, 2, 3, 4):
        ops.fix(n, 0, 0, 1)          # z=0 -> uz=0
    ops.equalDOF(2, 3, 1); ops.equalDOF(2, 6, 1); ops.equalDOF(2, 7, 1)  # +x face ux
    ops.equalDOF(4, 3, 2); ops.equalDOF(4, 7, 2); ops.equalDOF(4, 8, 2)  # +y face uy
    ops.equalDOF(5, 6, 3); ops.equalDOF(5, 7, 3); ops.equalDOF(5, 8, 3)  # +z face uz
    # phase 1: ramp lateral confining pressure p
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for n in (2, 3, 6, 7):
        ops.load(n, -p / 4.0, 0, 0)
    for n in (3, 4, 7, 8):
        ops.load(n, 0, -p / 4.0, 0)
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("FullGeneral")
    ops.test("NormDispIncr", 1e-6, 500, 0); ops.algorithm("KrylovNewton")
    ops.integrator("LoadControl", 1.0 / 100); ops.analysis("Static")
    if p > 0 and ops.analyze(100) != 0:
        return None
    ops.loadConst("-time", 0.0)
    # phase 2: axial reference load + displacement control on +z
    ops.timeSeries("Linear", 2); ops.pattern("Plain", 2, 2)
    for n in (5, 6, 7, 8):
        ops.load(n, 0, 0, -0.25)
    fcc, epk = 0.0, 0.0
    for cur in _adapt(5, 3, umax, umax / n0):
        sax = -ops.eleResponse(1, "material", 1, "stress")[2]
        if sax > fcc:
            fcc, epk = sax, -cur
    return fcc, epk


def main():
    print(f"{'p':>5} {'p/fc':>5} {'fcc_sim':>8} {'Mander':>7} {'Richart':>8} "
          f"{'s/Mand':>7} {'eps_pk':>7}")
    ok = True
    last_epk = -1.0
    for p in (0.0, 1.5, 3.0, 4.5, 6.0):
        out = run_conf(p)
        if out is None:
            print(f"{p:5.1f}  phase-1 (confining) diverged"); ok = False; continue
        fcc, epk = out
        r = fcc / mander(p)
        print(f"{p:5.1f} {p / FC:5.2f} {fcc:8.2f} {mander(p):7.2f} "
              f"{richart(p):8.2f} {r:7.2f} {epk:7.4f}")
        if p == 0.0 and abs(fcc - FC) > 0.01 * FC:
            ok = False
        if not (0.90 <= r <= 1.10):
            ok = False
        if epk < last_epk:
            ok = False
        last_epk = epk
    print("\nGATE 2:", "PASS — emergent confinement confirmed" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
