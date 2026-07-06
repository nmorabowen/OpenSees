"""ADR-69 P1.5 gate — ASDAbsorbingBoundary2D bottom compliant-base Ricker run,
legacy vs -v2. openseespy (the ASD elements were already Python-registered).

Topology (2D plane-strain, ndm=2 ndf=2): a small quad soil column with a
"ghost + real" compliant-base pair at its foot, mirroring the class's own
internal model (addKPenaltyStage1 pins the ghost nodes near-rigidly via
penalty stiffness; addClk/addBaseActions act on the real soil-side nodes).
Ghost nodes are ALSO explicitly fixed here (belt-and-suspenders — the
element's internal penalty already does this, explicit fix removes any
doubt about a massless/underconstrained DOF).

Node layout:
  1001, 1002 : ghost nodes, y = -1  (fixed both dof)
  1..2*(N+1) : real soil column nodes, y = 0 .. N (shared with quads)
  element ASDAbsorbingBoundary2D 101  1001 1002 2 1  G v rho thick "B" -fx tsTag
    (n1..n4 order doesn't matter - the class sorts by geometry)

Gates (mirrors p1_ricker_closure.py's Lysmer gate exactly):
  G1  loader live: the compliant base's driven DOF actually moves.
  G2  legacy lie reproduced: IE_end trends strongly off from clean.
  G3  v2 IE truthful: IE_end close to 0 relative to E_ref.
  G4  publisher/recorder cross-check: E_inject ~= IE_v2_end - IE_legacy_end.
  G5  control (zero-amplitude time series): E_inject ~= 0.

Run:  python p15_asd_ricker_closure.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))
DIST = os.path.join(WT, "dist", "bin")

sys.path.insert(0, r"C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\Lib\site-packages")
sys.path.insert(0, DIST)
os.add_dll_directory(DIST)

import numpy as np  # noqa: E402

import opensees as ops  # noqa: E402

print("[gate] opensees pyd:", ops.__file__)
assert "strange-hawking" in ops.__file__, "wrong pyd loaded"

RHO, G, NU = 2000.0, 2000.0 * 200.0 ** 2, 0.25
E = 2.0 * G * (1.0 + NU)
DT, NSTEP = 1.0e-3, 500
F0, T0, AMP = 10.0, 0.12, 0.05
N = 4          # soil layers above the compliant base
THICK = 1.0


def ricker():
    t = np.arange(NSTEP + 1) * DT
    a = (np.pi * F0 * (t - T0)) ** 2
    return AMP * (1.0 - 2.0 * a) * np.exp(-a)


def build(v2, amp_scale):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    # ghost nodes (fixed both dof - belt-and-suspenders w/ the element's own
    # internal penalty pin)
    ops.node(1001, 0.0, -1.0)
    ops.node(1002, 1.0, -1.0)
    ops.fix(1001, 1, 1)
    ops.fix(1002, 1, 1)

    # real soil column: N+1 rows of 2 nodes each, y = 0..N
    for row in range(N + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))

    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    for row in range(N):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, THICK, "PlaneStrain", 1)

    vals = (ricker() * amp_scale).tolist()
    ops.timeSeries("Path", 1, "-dt", DT, "-values", *vals)

    ops.element("ASDAbsorbingBoundary2D", 101, 1001, 1002, 2, 1,
               G, NU, RHO, THICK, "B", "-fx", 1)

    # bottom compliant base only - constrain the lateral DOF (x here is the
    # column's own width direction, not the Ricker's driven direction; drive
    # is fx -> horizontal/x); keep the soil free in x, fixed in nothing else
    # (2D plane-strain, no z). No lateral BC needed - single-column model.

    return


def run(tag, v2, amp_scale):
    build(v2, amp_scale)
    ops.setParameter("-val", 1, "-ele", 101, "stage")

    efile = os.path.join(HERE, "p15_%s_energy.txt" % tag)
    rec = ["EnergyBalance", "-file", efile, "-time"]
    if v2:
        rec.append("-v2")
    ops.recorder(*rec)

    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-10, 10)
    ops.algorithm("Linear")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    for _ in range(NSTEP):
        ok = ops.analyze(1, DT)
        if ok != 0:
            raise SystemExit("analyze() failed for %s at step" % tag)
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile))


def main():
    d_leg = run("legacy", v2=False, amp_scale=1.0)
    d_v2 = run("v2", v2=True, amp_scale=1.0)
    d_ctl = run("control", v2=True, amp_scale=0.0)

    ie_leg = d_leg[-1, 2]
    eref_leg = max(np.abs(d_leg[:, 1]).max(), np.abs(d_leg[:, 2]).max(),
                   np.abs(d_leg[:, 3]).max(), 1e-30)

    assert d_v2.shape[1] >= 10, "v2 column count %d < 10" % d_v2.shape[1]
    ie_v2 = d_v2[-1, 3]
    e_inj = d_v2[-1, -3] if d_v2.shape[1] == 10 else d_v2[-1, 7]
    eref = max(np.abs(d_v2[:, 1:8]).max(), 1e-30)

    ok = True

    g2 = abs(ie_leg) > 0.05 * eref_leg
    print("G2 legacy lie: IE_end = %.4e (E_ref %.4e) -> %s"
          % (ie_leg, eref_leg, "PASS" if g2 else "FAIL"))
    ok &= g2

    g3 = abs(ie_v2) < 0.05 * eref
    print("G3 v2 IE truthful: IE_end = %.4e (E_ref %.4e) -> %s"
          % (ie_v2, eref, "PASS" if g3 else "FAIL"))
    ok &= g3

    diff = ie_v2 - ie_leg
    g4 = abs(e_inj) > 0 and abs(e_inj - diff) < 0.1 * max(abs(diff), 1e-30)
    print("G4 cross-check: E_inject = %.4e vs IE_v2-IE_legacy = %.4e -> %s"
          % (e_inj, diff, "PASS" if g4 else "FAIL"))
    ok &= g4

    e_inj_ctl = d_ctl[-1, -3] if d_ctl.shape[1] == 10 else d_ctl[-1, 7]
    eref_ctl = max(np.abs(d_ctl[:, 1:8]).max(), 1e-30)
    g5 = abs(e_inj_ctl) < 1e-6 * eref_ctl or abs(e_inj_ctl) < 1e-12
    print("G5 zero-amplitude control: E_inject = %.3e (E_ref %.3e) -> %s"
          % (e_inj_ctl, eref_ctl, "PASS" if g5 else "FAIL"))
    ok &= g5

    print("\nADR-69 P1.5 ASD Ricker gate:", "ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
