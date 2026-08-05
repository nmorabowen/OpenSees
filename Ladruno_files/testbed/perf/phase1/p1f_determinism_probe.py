"""ADR-75 P1f — is threaded MKL PARDISO byte-reproducible RUN TO RUN?

ADR-75 P1 concluded it was, from ONE run per thread count. That design cannot
distinguish "deterministic" from "the same value came up twice" — with a 50/50
split a single-sample check passes half the time. This runs the SAME binary N
times at a FIXED thread count and counts distinct answers, which is the question
that was actually being asked.

Measured 2026-07-25 (14^3 Lane B, 8 steps, N=10):
    1 thread : 1 distinct / 10   reproducible
    4 threads: 2 distinct / 10   *** 5/5 split, ~1 ULP apart ***
An 8^3 model is reproducible even at 4 threads, so the effect is size-dependent.

HOW THIS CAME UP, because the debugging path is the reusable part: an A/B of the
P1f addA change reported "ux DIFFERS" at 4 threads. The tempting reading is "my
change broke exactness". The correct next move was to ask whether each binary
reproduces ITSELF — the pre-change binary showed the identical 5/5 split, which
proved the variation was MKL's and the change was exact (confirmed separately at
1 thread: 10/10 identical between old and new).

Run:
  python3.12 -S p1f_determinism_probe.py
  OPS_PYD          module dir under test (default dist/bin)
  OPS_PYD_OLD      optional second dir, to cross-check two builds
  OPS_PROBE_N      samples per configuration (default 10)
  OPS_PROBE_GRID   grid sizes to sweep (default 8,14)
"""
from __future__ import annotations
import os
import re
import subprocess
import sys
import tempfile

sys.stdout.reconfigure(line_buffering=True)

ROOT = (r"C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees"
        r"\opensees-pr-strategy-9c17b8")
NEW = os.environ.get("OPS_PYD", os.path.join(ROOT, "dist", "bin"))
OLD = os.environ.get("OPS_PYD_OLD")
N = int(os.environ.get("OPS_PROBE_N", "10"))
GRIDS = [int(g) for g in os.environ.get("OPS_PROBE_GRID", "8,14").split(",")]
PY = sys.executable

DRIVER = r'''
import os, sys
os.environ["MKL_NUM_THREADS"] = "{threads}"
os.environ["OMP_NUM_THREADS"] = "{threads}"
os.environ.setdefault("PMI_RANK", "0")
os.environ.setdefault("LADRUNO_OPENSEES_QUIET", "1")
DIST = r"{dist}"
sys.path.insert(0, DIST); os.add_dll_directory(DIST)
import opensees as ops
assert ops.__file__.lower().startswith(DIST.lower()), "WRONG PYD " + ops.__file__
E, nu = 200000.0, 0.3
Kb = E/(3.0*(1.0-2.0*nu)); G = E/(2.0*(1.0+nu)); LEN=100.0
n = {n}; NS = 8
def nid(i,j,k): return 1 + i + (n+1)*(j + (n+1)*k)
ops.wipe(); ops.model("basic","-ndm",3,"-ndf",3)
ops.nDMaterial("LadrunoJ2",1,Kb,G,"-iso","voce",250.0,0.0,0.0,2000.0)
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
    for i in range(n+1): ops.load(nid(i,j,n), 4.0e5, 0.0, -4.0e5)
ops.constraints("Plain"); ops.numberer("RCM")
ops.system(*{sysargs})
ops.test("NormDispIncr",1e-10,50); ops.algorithm("Newton")
ops.integrator("LoadControl", 1.0/NS); ops.analysis("Static")
for s in range(NS):
    if ops.analyze(1) != 0: print("FAILED", s); raise SystemExit(1)
print("RESULT ux=%.17e" % ops.nodeDisp(nid(n//2,n//2,n),1))
'''


def once(dist, n, threads, sysargs):
    src = DRIVER.format(dist=dist, n=n, threads=threads, sysargs=repr(sysargs))
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False) as f:
        f.write(src)
        p = f.name
    o = subprocess.run([PY, "-S", p], capture_output=True, text=True)
    os.unlink(p)
    m = re.search(r"RESULT ux=(\S+)", o.stdout + o.stderr)
    if not m:
        raise RuntimeError((o.stdout + o.stderr)[-400:])
    return m.group(1)


def probe(label, dist, n, threads, sysargs=("Pardiso",)):
    rs = [once(dist, n, threads, list(sysargs)) for _ in range(N)]
    uniq = sorted(set(rs))
    verdict = "reproducible" if len(uniq) == 1 else "*** NON-DETERMINISTIC ***"
    print(f"  {label:34s} {len(uniq)} distinct / {N} runs  {verdict}")
    for u in uniq:
        print(f"        {u}  (x{rs.count(u)})")
    return uniq


def main():
    print(f"P1f determinism probe | N={N} | grids={GRIDS}\n  NEW={NEW}"
          + (f"\n  OLD={OLD}" if OLD else ""))
    for n in GRIDS:
        ndof = 3 * (n + 1) ** 2 * n
        print(f"\n=== grid {n}^3 (~{ndof/1000:.1f}k DOF) ===")
        for threads in ("1", "4"):
            probe(f"NEW @ {threads} thread(s)", NEW, n, threads)
            if OLD:
                probe(f"OLD @ {threads} thread(s)", OLD, n, threads)


if __name__ == "__main__":
    main()
