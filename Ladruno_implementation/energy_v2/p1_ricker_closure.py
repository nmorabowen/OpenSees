"""ADR-69 P1 gate — Lysmer compliant-base Ricker run, legacy vs -v2.

Requires the WORKTREE Tcl exe (dist/bin/OpenSees.exe) built with:
  * the eleLoad -lysmerVelocityLoader hook (previously unreachable — P0.5),
  * the LysmerTriangle ABSORB_LEAK publisher,
  * the EnergyBalanceRecorder -v2 channel layout.

Gates (chosen to be free of the F2 confound — under implicit Newmark the
Lysmer element's damping force never enters the residual, so the SOLVE is
energy-inconsistent and RES genuinely cannot close; RES is also invariant
under the E_inject rebucket by design):
  G1  loader is live: the column actually moves (mid-node |vx| > 0).
  G2  legacy lie reproduced: IE_end strongly negative (books -W_inject).
  G3  v2 IE truthful: IE_end >= -2% of E_ref (leak rebucketted out).
  G4  publisher/recorder cross-check: E_inject ~= IE_v2_end - IE_legacy_end
      (same leak measured through two independent data paths).
  G5  radiation-only control: E_inject ~= 0, IE >= 0 (no false positives).

Run:  python p1_ricker_closure.py
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WT = os.path.abspath(os.path.join(HERE, "..", ".."))
EXE = os.path.join(WT, "dist", "bin", "OpenSees.exe")

RHO, VS, NU = 2000.0, 200.0, 0.25
VP = VS * np.sqrt(2.0 * (1.0 - NU) / (1.0 - 2.0 * NU))
E = 2.0 * RHO * VS**2 * (1.0 + NU)
DT, NSTEP = 1.0e-3, 600
F0, T0, AMP = 10.0, 0.15, 0.1


def ricker():
    t = np.arange(NSTEP + 1) * DT
    a = (np.pi * F0 * (t - T0)) ** 2
    return AMP * (1.0 - 2.0 * a) * np.exp(-a)


def deck(path, efile, vfile, v2, with_loader):
    vals = " ".join("%.10e" % v for v in ricker())
    lines = ["model basic -ndm 3 -ndf 3",
             "nDMaterial ElasticIsotropic 1 %.8e %.2f %.1f" % (E, NU, RHO)]
    nid = 0
    for k in range(11):
        for (x, y) in ((0, 0), (1, 0), (1, 1), (0, 1)):
            nid += 1
            lines.append("node %d %.1f %.1f %.1f" % (nid, x, y, float(k)))
    for n in range(1, 45):
        lines.append("fix %d 0 1 1" % n)
    for k in range(10):
        b = 4 * k
        lines.append("element stdBrick %d %d %d %d %d %d %d %d %d 1"
                     % (k + 1, b + 1, b + 2, b + 3, b + 4, b + 5, b + 6, b + 7, b + 8))
    lines.append("element LysmerTriangle 101 1 2 3 %.1f %.8f %.1f" % (RHO, VP, VS))
    lines.append("element LysmerTriangle 102 1 3 4 %.1f %.8f %.1f" % (RHO, VP, VS))
    if with_loader:
        lines.append("timeSeries Path 1 -dt %.4e -values {%s}" % (DT, vals))
        lines.append("pattern Plain 1 1 { eleLoad -ele 101 102 -type "
                     "-lysmerVelocityLoader 1 }")
    else:
        # radiation-only control: initial velocity kick on the top face
        for n in (41, 42, 43, 44):
            lines.append("setNodeVel %d 1 0.5 -commit" % n)
    rec = "recorder EnergyBalance -file {%s} -time" % efile
    if v2:
        rec += " -v2"
    lines.append(rec)
    lines.append("recorder Node -file {%s} -time -node 21 -dof 1 vel" % vfile)
    lines += ["constraints Plain", "numberer RCM", "system UmfPack",
              "test NormDispIncr 1.0e-10 10", "algorithm Linear",
              "integrator Newmark 0.5 0.25", "analysis Transient",
              "analyze %d %.4e" % (NSTEP, DT), "wipe", "exit"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def run(tag, v2, with_loader):
    tcl = os.path.join(HERE, "p1_%s.tcl" % tag)
    ef = os.path.join(HERE, "p1_%s_energy.txt" % tag).replace("\\", "/")
    vf = os.path.join(HERE, "p1_%s_vel.txt" % tag).replace("\\", "/")
    deck(tcl, ef, vf, v2, with_loader)
    r = subprocess.run([EXE, tcl], capture_output=True, text=True, timeout=300)
    if not os.path.exists(ef):
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise SystemExit("run %s produced no energy file" % tag)
    return np.atleast_2d(np.loadtxt(ef)), np.atleast_2d(np.loadtxt(vf))


def main():
    assert os.path.exists(EXE), "worktree OpenSees.exe not built: %s" % EXE
    ok = True

    d_leg, v_leg = run("legacy", v2=False, with_loader=True)
    d_v2, v_v2 = run("v2", v2=True, with_loader=True)
    d_ctl, _ = run("control_v2", v2=True, with_loader=False)

    # legacy: time KE IE DW ULW RES ERR
    ie_leg = d_leg[-1, 2]
    eref_leg = max(np.abs(d_leg[:, 1]).max(), np.abs(d_leg[:, 2]).max(),
                   np.abs(d_leg[:, 3]).max(), 1e-30)
    # v2 (ABSORB_LEAK declared, no LNVD):
    #   time KE_ele KE_nod IE DW_ele DW_nod ULW E_inject RES ERR
    assert d_v2.shape[1] == 10, "v2 column count %d != 10" % d_v2.shape[1]
    ie_v2 = d_v2[-1, 3]
    e_inj = d_v2[-1, 7]
    eref = max(np.abs(d_v2[:, 1:8]).max(), 1e-30)

    g1 = np.abs(v_v2[:, 1]).max() > 1e-6
    print("G1 loader live: max|vx(mid)| = %.3e -> %s"
          % (np.abs(v_v2[:, 1]).max(), "PASS" if g1 else "FAIL"))

    g2 = ie_leg < -0.2 * eref_leg
    print("G2 legacy lie: IE_end = %.4e (E_ref %.4e) -> %s"
          % (ie_leg, eref_leg, "PASS" if g2 else "FAIL"))

    g3 = ie_v2 > -0.02 * eref
    print("G3 v2 IE truthful: IE_end = %.4e (E_ref %.4e) -> %s"
          % (ie_v2, eref, "PASS" if g3 else "FAIL"))

    diff = ie_v2 - ie_leg
    g4 = e_inj > 0 and abs(e_inj - diff) < 0.05 * max(abs(diff), 1e-30)
    print("G4 cross-check: E_inject = %.4e vs IE_v2-IE_legacy = %.4e -> %s"
          % (e_inj, diff, "PASS" if g4 else "FAIL"))

    ie_ctl = d_ctl[-1, 3]
    einj_ctl = d_ctl[-1, 7]
    eref_ctl = max(np.abs(d_ctl[:, 1:8]).max(), 1e-30)
    g5 = abs(einj_ctl) < 1e-3 * eref_ctl and ie_ctl > -1e-3 * eref_ctl
    print("G5 radiation control: E_inject = %.3e, IE_end = %.3e -> %s"
          % (einj_ctl, ie_ctl, "PASS" if g5 else "FAIL"))

    ok = g1 and g2 and g3 and g4 and g5
    print("\nADR-69 P1 Ricker gate:", "ALL PASS" if ok else "FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
