"""ADR-75 P1e follow-up — what does `-krylov` cost when CGS starts FAILING?

Every P1e measurement ran on a HARDENING tangent, where CGS won 340/340 solves.
That left the real question open: near a limit point, where the tangent goes
near-singular, PARDISO's fallback costs up to ~half a wasted iteration budget
AND the factorization. Is -krylov a net loss there?

Softening J2 (Hiso < 0) under LoadControl marches straight into a limit point,
so the last steps before divergence are exactly the ill-conditioned regime.
-matrixType 0 (mtype 11) is the only legal configuration: -matrixType 2 refuses
-krylov outright, and softening is precisely what makes a tangent indefinite.

Reports, per configuration: steps survived, Newton iterations, wall time, and
the full CGS decision trace (wins / fallbacks / cgs_error codes / perturbed
pivots). A DIFFERENT failure step between direct and -krylov is itself a
finding: it would mean the inexact solve changed the equilibrium path.
"""
import os
import re
import subprocess
import tempfile

DIST = (r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees"
        r"\opensees-pr-strategy-9c17b8\dist\bin")
PY = r"C:\Users\nmb\AppData\Local\Python\bin\python3.12"

DRIVER = r'''
import os, sys, time
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
DIST = r"{dist}"
sys.path.insert(0, DIST); os.add_dll_directory(DIST)
import opensees as ops
E, nu = 200000.0, 0.3
Kb = E/(3.0*(1.0-2.0*nu)); G = E/(2.0*(1.0+nu)); LEN=100.0
n = {n}; NS = {ns}
def nid(i,j,k): return 1 + i + (n+1)*(j + (n+1)*k)
ops.wipe(); ops.model("basic","-ndm",3,"-ndf",3)
# Hiso < 0 => LINEAR SOFTENING: sigma_y = s0 + Hiso*ep. The tangent loses
# positive definiteness once a band yields, and marches toward singular.
ops.nDMaterial("LadrunoJ2",1,Kb,G,"-iso","voce",250.0,0.0,0.0,{hiso})
for k in range(n+1):
    for j in range(n+1):
        for i in range(n+1): ops.node(nid(i,j,k), i*LEN, j*LEN, k*LEN)
for j in range(n+1):
    for i in range(n+1): ops.fix(nid(i,j,0),1,1,1)
e=1
for k in range(n):
    for j in range(n):
        for i in range(n):
            ops.element("LadrunoBrick",e, nid(i,j,k),nid(i+1,j,k),nid(i+1,j+1,k),
                        nid(i,j+1,k), nid(i,j,k+1),nid(i+1,j,k+1),
                        nid(i+1,j+1,k+1),nid(i,j+1,k+1),1); e+=1
ops.timeSeries("Linear",1); ops.pattern("Plain",1,1)
for j in range(n+1):
    for i in range(n+1): ops.load(nid(i,j,n), {load}, 0.0, -{load})
ops.constraints("Plain"); ops.numberer("RCM")
ops.system(*{sysargs})
ops.test("NormDispIncr",1e-7,40); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0/NS); ops.analysis("Static")
done=0; iters=0
t0=time.perf_counter()
for s in range(NS):
    if ops.analyze(1) != 0:
        break
    done += 1; iters += ops.testIter()
wall=time.perf_counter()-t0
print("RESULT steps=%d/%d iters=%d wall=%.4f ux=%.12e" % (
      done, NS, iters, wall, ops.nodeDisp(nid(n//2,n//2,n),1)))
'''


def run(sysargs, hiso, load, n=10, ns=40):
    src = DRIVER.format(dist=DIST, sysargs=repr(sysargs), hiso=hiso, n=n, ns=ns,
                        load=load)
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
        f.write(src)
        p = f.name
    out = subprocess.run([PY, "-S", p], capture_output=True, text=True)
    txt = out.stdout + out.stderr
    os.unlink(p)
    r = re.search(r"RESULT steps=(\d+)/(\d+) iters=(\d+) wall=([\d.]+) ux=(\S+)", txt)
    wins = [int(x) for x in re.findall(r"CGS converged in (\d+)", txt)]
    falls = re.findall(r"CGS gave up after (\d+) iteration\(s\) \(cgs_error (\d+)\)", txt)
    pert = len(re.findall(r"perturbed \d+ pivot", txt))
    if not r:
        return None
    return dict(steps=int(r.group(1)), ns=int(r.group(2)), iters=int(r.group(3)),
                wall=float(r.group(4)), ux=float(r.group(5)),
                wins=len(wins), maxit=max(wins) if wins else 0,
                falls=len(falls), errs=sorted({c for _, c in falls}), pert=pert)


print("Softening J2, -matrixType 0 (mtype 11), 4 threads, 10^3 (~3.6k DOF), "
      "LoadControl to divergence\n")
hdr = (f"{'load':>8s} {'Hiso':>7s} {'config':9s} {'steps':>7s} {'iters':>6s} "
       f"{'wall_s':>8s} {'CGSwin':>7s} {'maxit':>6s} {'FALLB':>6s} {'errs':>8s} "
       f"{'pert':>5s}  ux")
print(hdr)
print("-" * len(hdr))
# Sweep LOAD as well as softening slope: the first attempt kept every row at
# 40/40 steps and 101 iterations with ux barely moving from Hiso 0 to -6000,
# i.e. the plastic zone was too small for softening to bite at all. A limit
# point needs the structure driven well past first yield.
CASES = [(8.0e5, -2000.0), (8.0e5, -6000.0),
         (1.2e6, -2000.0), (1.2e6, -6000.0),
         (1.6e6, -6000.0), (1.6e6, -20000.0)]
for load, hiso in CASES:
    rows = {}
    for label, args in (("direct", ["Pardiso", "-stats"]),
                        ("krylov6", ["Pardiso", "-krylov", 6, "-stats"])):
        r = run(args, hiso, load)
        rows[label] = r
        if r is None:
            print(f"{load:8.1e} {hiso:7.0f} {label:9s}   <no RESULT line>")
            continue
        print(f"{load:8.1e} {hiso:7.0f} {label:9s} {r['steps']:3d}/{r['ns']:<3d} "
              f"{r['iters']:6d} {r['wall']:8.3f} {r['wins']:7d} {r['maxit']:6d} "
              f"{r['falls']:6d} {str(r['errs']):>8s} {r['pert']:5d}  {r['ux']:.9e}")
    d, k = rows.get("direct"), rows.get("krylov6")
    if d and k:
        note = []
        if d["steps"] != k["steps"]:
            note.append(f"*** DIFFERENT FAILURE STEP: direct {d['steps']} vs "
                        f"krylov {k['steps']} ***")
        if d["steps"] == k["steps"] and d["wall"] > 0:
            note.append(f"speedup {d['wall']/k['wall']:.3f}x")
            if abs(d["ux"]) > 0:
                note.append(f"ux rel diff {abs(k['ux']-d['ux'])/abs(d['ux']):.2e}")
        print(f"        -> {'  |  '.join(note)}\n")
