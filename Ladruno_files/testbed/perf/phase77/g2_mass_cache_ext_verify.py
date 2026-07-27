"""ADR-77 G2 extension — acceptance for the shared LadrunoMassCache rollout.

Five elements got the cache (BezierTet10, BezierTri6, LadrunoQuad, LadrunoLST,
LadrunoSolidShell). The claim is the same as the brick's: the cache is
unobservable in results — only in time. So per element:

  A/B: an implicit transient run (Newmark + Newton, material-rho mass) with the
       cache ON vs `-noMassCache` must give BYTE-EQUAL displacements. UmfPack,
       serial ⇒ deterministic, so equality is exact or it is a bug.

Plus one guard-liveness check through the SHARED helper (the brick's guards are
separate code, so the helper's own invalidation path needs its own witness):

  setNodeCoord mid-analysis on a LadrunoQuad deck — the coords guard must trip,
  re-form, and keep cached ≡ uncached byte-equal. (rho-liveness was proven on
  the brick; the helper's rho guard is the same compare loop as the coords one.)

Usage:  OPS_PYD=... python3.12 -S g2_mass_cache_ext_verify.py
"""
from __future__ import annotations

import os
import sys

sys.stdout.reconfigure(line_buffering=True)
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")

DIST = os.environ.get(
    "OPS_PYD",
    r"C:\Users\nmb\Documents\Github\OpenSees"
    r"\.claude\worktrees\pardisio-profiling-0a03b1\dist\bin",
)
sys.path.insert(0, DIST)
if os.path.isdir(DIST):
    os.add_dll_directory(DIST)
import opensees as ops  # noqa: E402

DT, NSTEPS, RHO = 1.0e-4, 3, 7.85e-9


def analysis():
    ops.constraints("Plain")
    ops.numberer("RCM")
    ops.system("UmfPack")
    ops.test("NormDispIncr", 1.0e-8, 50)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")


def run(nsteps=NSTEPS):
    for _ in range(nsteps):
        if ops.analyze(1, DT) != 0:
            raise RuntimeError("analyze failed")


def tet10(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3, RHO)
    # 4 vertices then mid-edges (1-2)(2-3)(1-3)(1-4)(3-4)(2-4)
    V = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1)}
    mids = [(1, 2), (2, 3), (1, 3), (1, 4), (3, 4), (2, 4)]
    for t, c in V.items():
        ops.node(t, *[float(x) for x in c])
    for k, (a, b) in enumerate(mids, start=5):
        ops.node(k, *[(V[a][i] + V[b][i]) / 2.0 for i in range(3)])
    for n in (1, 2, 3, 5, 6, 7):        # z=0 face
        ops.fix(n, 1, 1, 1)
    ops.element("BezierTet10", 1, *range(1, 11), 1, "-cMass", *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(4, 10.0, 5.0, -20.0)
    analysis()
    run()
    return tuple(ops.nodeDisp(4))


def tri6(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3, RHO)
    C = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0)}
    mids = [(1, 2), (2, 3), (3, 1)]
    for t, c in C.items():
        ops.node(t, *c)
    for k, (a, b) in enumerate(mids, start=4):
        ops.node(k, (C[a][0] + C[b][0]) / 2.0, (C[a][1] + C[b][1]) / 2.0)
    for n in (1, 2, 4):
        ops.fix(n, 1, 1)
    ops.element("BezierTri6", 1, 1, 2, 3, 4, 5, 6, 1.0, "PlaneStress", 1,
                "-cMass", *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 6.0, -5.0)
    analysis()
    run()
    return tuple(ops.nodeDisp(3))


def quad(opts, kick_coord=False):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3, RHO)
    for t, c in {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (1.0, 1.0), 4: (0.0, 1.0)}.items():
        ops.node(t, *c)
    ops.fix(1, 1, 1)
    ops.fix(2, 1, 1)
    ops.element("LadrunoQuad", 1, 1, 2, 3, 4, 1, "-thick", 1.0, *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 4.0, -6.0)
    analysis()
    if not kick_coord:
        run()
    else:
        # guard liveness through the SHARED helper: perturb a free node's
        # coordinates mid-analysis; the coords guard must trip and re-form.
        run(2)
        # signature is setNodeCoord(tag, dim, value), dim 1-based
        ops.setNodeCoord(4, 1, 0.05)
        ops.setNodeCoord(4, 2, 1.02)
        ops.domainChange()
        run(2)
    return tuple(ops.nodeDisp(3))


def lst(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3, RHO)
    C = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0)}
    mids = [(1, 2), (2, 3), (3, 1)]
    for t, c in C.items():
        ops.node(t, *c)
    for k, (a, b) in enumerate(mids, start=4):
        ops.node(k, (C[a][0] + C[b][0]) / 2.0, (C[a][1] + C[b][1]) / 2.0)
    for n in (1, 2, 4):
        ops.fix(n, 1, 1)
    ops.element("LadrunoLST", 1, 1, 2, 3, 4, 5, 6, 1, "-thick", 1.0, *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(3, 6.0, -5.0)
    analysis()
    run()
    return tuple(ops.nodeDisp(3))


def solidshell(opts):
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ops.nDMaterial("ElasticIsotropic", 1, 200000.0, 0.3, RHO)
    crds = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
            (0, 0, 0.1), (1, 0, 0.1), (1, 1, 0.1), (0, 1, 0.1)]
    for t, c in enumerate(crds, start=1):
        ops.node(t, *[float(x) for x in c])
    for n in (1, 2, 3, 4):
        ops.fix(n, 1, 1, 1)
    ops.element("LadrunoSolidShell", 1, *range(1, 9), 1, "-nz", 2, *opts)
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(7, 2.0, 1.0, -5.0)
    analysis()
    run()
    return tuple(ops.nodeDisp(7))


def main():
    fails = 0
    print("== A/B: cache on vs -noMassCache (must be byte-equal) ==")
    for name, fn in (("BezierTet10 (-cMass)", tet10),
                     ("BezierTri6 (-cMass)", tri6),
                     ("LadrunoQuad", quad),
                     ("LadrunoLST", lst),
                     ("LadrunoSolidShell", solidshell)):
        a = fn(())
        b = fn(("-noMassCache",))
        ok = a == b
        fails += 0 if ok else 1
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<22} on {a}\n"
              f"        {'':<22} off {b}" if not ok else
              f"  PASS  {name:<22} {a}")

    print("== coords-guard liveness (setNodeCoord mid-analysis, shared helper) ==")
    a = quad((), kick_coord=True)
    b = quad(("-noMassCache",), kick_coord=True)
    ok = a == b
    fails += 0 if ok else 1
    print(f"  {'PASS' if ok else 'FAIL'}  LadrunoQuad kick     on {a} off {b}")

    print(f"\n{'ALL CHECKS PASSED' if not fails else f'{fails} CHECK(S) FAILED'}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
