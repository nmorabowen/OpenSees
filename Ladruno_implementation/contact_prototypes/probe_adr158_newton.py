"""Python Newton on the OpenSees contact residual (printB) for step K+1, with the assembled tangent
(printA) or a central-FD tangent. Springs + loads are added analytically (setNodeDisp does not
update() elements). Run: run_with_bin <bin> pynewton.py k=<K> tan=os|fd [model opts]"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import probe_adr158_mortar_tangent_fd as F
import opensees as ops

TAN = F.args.get("tan", "os")
NIT = int(F.args.get("nit", 30))
F.NOCONTACT[0] = True
F.build(); ops.test("NormUnbalance", 1e-7, 1, 0); ops.integrator("LoadControl", 0.0); ops.analyze(1)
nq = ops.systemSize(); KLIN = np.array(ops.printA("-ret")).reshape(nq, nq).T
F.NOCONTACT[0] = False
allslave, ridge = F.build()
if F.K == 0:
    ops.test("NormUnbalance", 1e-7, 1, 0); ops.integrator("LoadControl", 0.0)
    print("  init-only analyze rc", ops.analyze(1))   # fails after 1 it, reverts: initialises the model
for s in range(F.K):
    rc = ops.analyze(1)
    assert rc == 0, f"OS failed at step {s+1}"
neq = ops.systemSize()
eq = {}
for n in allslave:
    for d in range(3): eq[(n, d)] = ops.nodeDOFs(n)[d]
P = np.zeros(neq)
def ridge_or_surf():
    return [n for n in allslave if n < 3000]
for n in ridge_or_surf(): P[eq[(n, 2)]] = -F.PZ
lam_c = F.K / F.NST

uc = np.zeros(neq)
for n in allslave:
    for d in range(3): uc[eq[(n, d)]] = ops.nodeDisp(n, d + 1)


def setu(u):
    for (n, d), e in eq.items(): ops.setNodeDisp(n, d + 1, u[e], "-commit")


def B(u):
    setu(u); return np.array(ops.printB("-ret"))


c0 = B(uc) + lam_c * P - KLIN @ uc if F.K > 0 else np.zeros(neq)


def R(u, lam):
    return lam * P - KLIN @ u + B(u) - c0


def Kos(u):
    setu(u); return np.array(ops.printA("-ret")).reshape(neq, neq).T   # Matrix data is column-major


def Kfd(u, lam, h=1e-9):
    K = np.zeros((neq, neq))
    for e in range(neq):
        up = u.copy(); up[e] += h; um = u.copy(); um[e] -= h
        K[:, e] = -(R(up, lam) - R(um, lam)) / (2 * h)
    return K


try:
    import probe_adr158_pair_replica as FP
except ImportError:
    FP = None
mq, sq = F.GEOM["mq"], F.GEOM["sq"]
MFS = [mq[i:i + 4] for i in range(0, len(mq), 4)]; SFS = [sq[i:i + 4] for i in range(0, len(sq), 4)]
SKIP_SLIVER = float(F.args.get("skipsliver", 0))


def Knp(u, consistent):
    K = KLIN.copy()
    for sf in SFS:
        e = [eq[(t, d)] for t in sf for d in range(3)]
        us = u[e].reshape(4, 3); X0 = np.array([ops.nodeCoord(t) for t in sf])
        for mf in MFS:
            Xm = np.array([ops.nodeCoord(t) for t in mf])
            f, info = FP.force(X0 + us, Xm, F.EPS, F.EPS * F.ETR, F.MU, max(F.COH, 0), us)
            if info is None or info["area"] < SKIP_SLIVER: continue
            K[np.ix_(e, e)] += FP.kan(info, F.EPS, F.EPS * F.ETR, F.MU, max(F.COH, 0), consistent)
    return K


lam = (F.K + 1) / F.NST
u = uc.copy()
HIST = [u.copy()]
print(f"step {F.K+1} tan={TAN}  |R(uc)| at lam_c = {np.linalg.norm(R(uc, lam_c)):.2e}")
for it in range(NIT):
    r = R(u, lam)
    nr = np.linalg.norm(r)
    K = Kos(u) if TAN == "os" else (Kfd(u, lam) if TAN == "fd" else Knp(u, TAN == "np"))
    du = np.linalg.solve(K, r)
    if "chk" in F.args:
        Ko = Kos(u); Kn = Knp(u, F.args["chk"] == "1"); sc = np.abs(Kn - KLIN).max()
        print(f"     chk |Kos-Knp|/|Kc| = {np.abs(Ko - Kn).max() / sc:.2e}")
        if it == 1:
            i, j = np.unravel_index(np.abs(Ko - Kn).argmax(), Ko.shape)
            inv = {e_: k_ for k_, e_ in eq.items()}; ni, nj = inv[i][0], inv[j][0]
            ei = [eq[(ni, d)] for d in range(3)]; ej = [eq[(nj, d)] for d in range(3)]
            np.set_printoptions(precision=1, suppress=True, linewidth=150)
            print("     worst", inv[i], inv[j], ops.nodeCoord(ni), ops.nodeCoord(nj))
            print((Ko - KLIN)[np.ix_(ei, ej)]); print((Kn - KLIN)[np.ix_(ei, ej)])
    print(f"  it {it:2d} |R|={nr:.3e} |du|={np.linalg.norm(du):.3e}")
    if nr < 1e-7:
        break
    u = u + du
    HIST.append(u.copy())
np.save(F.args.get("save", "hist.npy"), np.array(HIST))
import json; json.dump({str(k): v for k, v in {f"{n},{d}": e for (n, d), e in eq.items()}.items()}, open("eqmap.json", "w"))
