"""Does -krylov diverge from direct BEFORE the limit point, or only after?

The softening sweep showed ux rel diffs up to 20.8 and a different give-up
step. That is either (a) a solver defect, or (b) post-limit-point path chaos,
where LoadControl on a softening structure is ill-posed and ANY perturbation
sends the two runs onto different, equally meaningless branches.

Decisive test: stop at increasing step counts and compare. If the runs agree
while the response is still physical and only separate as the limit point is
approached, it is (b) — and the practical warning is about reproducibility in
the post-peak regime, not about correctness.
"""
import os
import re
import subprocess
import tempfile

DIST = (r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees"
        r"\opensees-pr-strategy-9c17b8\dist\bin")
PY = r"C:\Users\nmb\AppData\Local\Python\bin\python3.12"

DRIVER = r'''
import os, sys
os.environ.setdefault("MKL_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
DIST = r"{dist}"
sys.path.insert(0, DIST); os.add_dll_directory(DIST)
import opensees as ops
E, nu = 200000.0, 0.3
Kb = E/(3.0*(1.0-2.0*nu)); G = E/(2.0*(1.0+nu)); LEN=100.0
n = 10; NS = 40; STOP = {stop}
def nid(i,j,k): return 1 + i + (n+1)*(j + (n+1)*k)
ops.wipe(); ops.model("basic","-ndm",3,"-ndf",3)
ops.nDMaterial("LadrunoJ2",1,Kb,G,"-iso","voce",250.0,0.0,0.0,-2000.0)
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
    for i in range(n+1): ops.load(nid(i,j,n), 8.0e5, 0.0, -8.0e5)
ops.constraints("Plain"); ops.numberer("RCM")
ops.system(*{sysargs})
ops.test("NormDispIncr",1e-7,40); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0/NS); ops.analysis("Static")
done=0
for s in range(STOP):
    if ops.analyze(1) != 0: break
    done += 1
print("RESULT done=%d ux=%.14e" % (done, ops.nodeDisp(nid(n//2,n//2,n),1)))
'''


def run(sysargs, stop):
    src = DRIVER.format(dist=DIST, sysargs=repr(sysargs), stop=stop)
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
        f.write(src)
        p = f.name
    out = subprocess.run([PY, "-S", p], capture_output=True, text=True)
    os.unlink(p)
    m = re.search(r"RESULT done=(\d+) ux=(\S+)", out.stdout + out.stderr)
    return (int(m.group(1)), float(m.group(2))) if m else (None, None)


print("Softening J2 (Hiso -2000), load 8e5, -matrixType 0. Direct vs -krylov 6")
print("stopping at increasing step counts.\n")
print(f"{'stop':>5s} {'done':>5s} {'direct ux':>22s} {'krylov ux':>22s} {'rel diff':>10s}")
print("-" * 70)
for stop in (5, 10, 15, 20, 24, 28, 30, 32, 34, 35):
    dn, du = run(["Pardiso"], stop)
    kn, ku = run(["Pardiso", "-krylov", 6], stop)
    if du is None or ku is None:
        print(f"{stop:5d}   <no RESULT>")
        continue
    rel = abs(ku - du) / abs(du) if du else float("nan")
    flag = ""
    if dn != kn:
        flag = f"  *** step mismatch {dn} vs {kn} ***"
    elif rel > 1e-9:
        flag = "  <-- separates here"
    print(f"{stop:5d} {dn:5d} {du:22.14e} {ku:22.14e} {rel:10.2e}{flag}")
