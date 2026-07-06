"""ADR-69 P1.6 gate — ASDAbsorbingBoundary2D LATERAL free-field transfer,
initial-velocity free vibration, legacy vs -v2. openseespy.

Topology (2D plane-strain, ndm=2 ndf=2): the p15 quad soil column at x in
[0,1], plus a free-field column of nodes at x = -LX and one lateral "L"
ASDAbsorbingBoundary2D per soil layer. The class sorts nodes by geometry
(SorterLeft: x asc, then y asc), so the outer pair (x=-LX) become the
FREE-FIELD nodes 0-1 (addKff/addMff/addCff column dynamics) and the inner
pair (x=0, shared with the quads) the SOIL nodes 2-3 that addRffToSoil
writes onto.

Fixes: soil base + FF BASE only. The upper FF nodes MUST stay free — stage-1
penalty machinery skips vertical boundaries entirely (addKPenaltyStage1
early-returns), so the FF column's dynamics live in addKff/addMff/addCff;
fixing them would zero the FF strain and kill addRffToSoil. (Deliberately
NOT the p15 ghost-fix pattern.)

Excitation: INITIAL-VELOCITY free vibration (both DOFs, all free nodes; the
transfer pairs cross-axis, so both axes must move for first-order work).
UniformExcitation is REJECTED for this gate, for two source-verified
reasons (first p16 attempt failed on both):
  1. ASDAbsorbingBoundary2D::addInertiaLoadToUnbalance is a deliberate
     no-op ("we don't need this!"), so the free-field masses never receive
     -M*ug'' and the FF column rides rigidly (zero strain, zero transfer).
  2. The quads DO implement it by storing -M*ug'' in their element load
     vector Q, and getResistingForce returns K*u - Q, so the seismic input
     work hides inside the recorder's IE (observed: IE = -DW exactly, RES
     accidentally closed) — the legacy-IE-lie signal would be dominated by
     input work the v2 rebucket deliberately does not touch.
With a pure initial-velocity drive there are no patterns, Q = 0, IE is
pure strain energy, and the only IE polluter left is the lateral transfer
itself. "L"/"R" take NO -fx/-fy args (parser only accepts the time-series
block for BND_BOTTOM).

Gates (mirrors p15_asd_ricker_closure.py):
  G1  free-field column live: an upper FF node's x-displacement moves.
  G2  legacy lie reproduced: IE_end trends strongly off from clean.
  G3  v2 IE truthful: IE_end close to 0 relative to E_ref.
  G4  publisher/recorder cross-check: E_inject ~= IE_v2_end - IE_legacy_end.
  G5  control (zero-amplitude): E_inject ~= 0.

Run:  python p16_asd_lateral_ricker_closure.py
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
DT, NSTEP = 1.0e-3, 800     # long enough for several transfer cycles
V0 = 0.5                    # initial velocity on all free nodes, both DOFs
N = 4          # soil layers
THICK = 1.0
LX = 1.0       # FF column horizontal offset (= tributary width)


def build(amp_scale):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)

    # soil column nodes: rows y = 0..N at x = 0 and x = 1 (as p15)
    for row in range(N + 1):
        ops.node(1 + 2 * row, 0.0, float(row))
        ops.node(2 + 2 * row, 1.0, float(row))
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)

    # free-field column: x = -LX, one node per row; ONLY the base fixed
    for row in range(N + 1):
        ops.node(1000 + row, -LX, float(row))
    ops.fix(1000, 1, 1)

    ops.nDMaterial("ElasticIsotropic", 1, E, NU, RHO)
    for row in range(N):
        n1, n2 = 1 + 2 * row, 2 + 2 * row
        n3, n4 = 2 + 2 * (row + 1), 1 + 2 * (row + 1)
        ops.element("quad", 200 + row, n1, n2, n3, n4, THICK, "PlaneStrain", 1)

    # one lateral ASD element per layer; node order arbitrary (geometry sort)
    tags = []
    for row in range(N):
        tag = 100 + row
        ops.element("ASDAbsorbingBoundary2D", tag,
                    1000 + row, 1000 + row + 1,   # FF pair (x=-LX)
                    1 + 2 * row, 1 + 2 * (row + 1),  # soil pair (x=0)
                    G, NU, RHO, THICK, "L")
        tags.append(tag)

    return tags


def run(tag, v2, amp_scale):
    ele_tags = build(amp_scale)
    for t in ele_tags:
        ops.setParameter("-val", 1, "-ele", t, "stage")

    # initial-velocity drive (no patterns; see module docstring): both DOFs
    # on every free node — soil interior/edge rows above the base and the
    # free-field column above its base
    v0 = V0 * amp_scale
    for row in range(1, N + 1):
        ops.setNodeVel(1 + 2 * row, 1, v0, "-commit")
        ops.setNodeVel(1 + 2 * row, 2, v0, "-commit")
        ops.setNodeVel(2 + 2 * row, 1, v0, "-commit")
        ops.setNodeVel(2 + 2 * row, 2, v0, "-commit")
        ops.setNodeVel(1000 + row, 1, v0, "-commit")
        ops.setNodeVel(1000 + row, 2, v0, "-commit")

    # mass-proportional Rayleigh so the FREE-FIELD column also decays: the
    # lateral addClk dashpot writes only SOIL rows (one-way coupling), so an
    # undamped FF column rings forever and its genuine spring energy keeps
    # IE_end from settling toward 0 (G3 would measure real elastic energy,
    # not a lie). addCff applies element Rayleigh to the FF mass, and the
    # recorder's DW picks all of it up via the element C matrices. alphaM
    # only — betaK on a dashpot-carrying element is the classic trap.
    ops.rayleigh(10.0, 0.0, 0.0, 0.0)

    efile = os.path.join(HERE, "p16_%s_energy.txt" % tag)
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
    ff_probe = 1000 + N   # top free-field node
    ff_max = 0.0
    for _ in range(NSTEP):
        ok = ops.analyze(1, DT)
        if ok != 0:
            raise SystemExit("analyze() failed for %s" % tag)
        ff_max = max(ff_max, abs(ops.nodeDisp(ff_probe, 1)))
    ops.wipe()
    return np.atleast_2d(np.loadtxt(efile)), ff_max


def main():
    d_leg, ff_leg = run("legacy", v2=False, amp_scale=1.0)
    d_v2, _ = run("v2", v2=True, amp_scale=1.0)
    d_ctl, _ = run("control", v2=True, amp_scale=0.0)

    ok = True

    g1 = ff_leg > 1e-8
    print("G1 free-field column live: max |u_x(FF top)| = %.4e -> %s"
          % (ff_leg, "PASS" if g1 else "FAIL"))
    ok &= g1

    ie_leg = d_leg[-1, 2]
    eref_leg = max(np.abs(d_leg[:, 1]).max(), np.abs(d_leg[:, 2]).max(),
                   np.abs(d_leg[:, 3]).max(), 1e-30)
    g2 = abs(ie_leg) > 0.05 * eref_leg
    print("G2 legacy lie: IE_end = %.4e (E_ref %.4e) -> %s"
          % (ie_leg, eref_leg, "PASS" if g2 else "FAIL"))
    ok &= g2

    assert d_v2.shape[1] >= 10, "v2 column count %d < 10" % d_v2.shape[1]
    ie_v2 = d_v2[-1, 3]
    # E_inject: front-anchored col 7 (right after ULW, with -time) — never
    # index from the back, E_lnvd may or may not be declared in-process
    e_inj = d_v2[-1, 7]
    eref = max(np.abs(d_v2[:, 1:8]).max(), 1e-30)
    g3 = abs(ie_v2) < 0.05 * eref
    print("G3 v2 IE truthful: IE_end = %.4e (E_ref %.4e) -> %s"
          % (ie_v2, eref, "PASS" if g3 else "FAIL"))
    ok &= g3

    diff = ie_v2 - ie_leg
    g4 = abs(e_inj) > 0 and abs(e_inj - diff) < 0.1 * max(abs(diff), 1e-30)
    print("G4 cross-check: E_inject = %.4e vs IE_v2-IE_legacy = %.4e -> %s"
          % (e_inj, diff, "PASS" if g4 else "FAIL"))
    ok &= g4

    e_inj_ctl = d_ctl[-1, 7]
    eref_ctl = max(np.abs(d_ctl[:, 1:8]).max(), 1e-30)
    g5 = abs(e_inj_ctl) < 1e-6 * eref_ctl or abs(e_inj_ctl) < 1e-12
    print("G5 zero-velocity control: E_inject = %.3e (E_ref %.3e) -> %s"
          % (e_inj_ctl, eref_ctl, "PASS" if g5 else "FAIL"))
    ok &= g5

    print("\nADR-69 P1.6 ASD lateral free-vibration gate:",
          "ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
