"""ADR-158 (R0.7) FD probe -- the binary half of the WP-158 diagnosis: assembled mortar-friction tangent (printA) vs central FD of the residual (printB)
on the creased roof (ADR-157 tests b/c geometry). Run with the fork's opensees module first on sys.path: python probe_adr158_mortar_tangent_fd.py [opts]
Traps (LEDGER_quirks, WP-158): setNodeDisp needs -commit (else the other DOFs reset to committed);
printA -ret is column-major (transpose after reshape); zeroLength/brick forces are not update()d by
setNodeDisp, so their exact linear part is added analytically.
opts: split=0/1 mu=<float> coh=<float> steps=<n> k=<probe after k steps> h=<fd step> eps=<epsN>
"""
import math, sys
import numpy as np
import opensees as ops

args = dict(a.split("=") for a in sys.argv[1:])
SPLIT = int(args.get("split", 0))
MU = float(args.get("mu", 0.0))
COH = float(args.get("coh", 10.0))
EPS = float(args.get("eps", 1e6))
ETR = float(args.get("etr", 1.0))
DELTA = float(args.get("delta", 1e-3))
PZ = float(args.get("pz", 40.0))
NST = int(args.get("steps", 4))
K = int(args.get("k", 2))
H = float(args.get("h", 1e-9))
FLAGS = args.get("flags", "")
NX = int(args.get("nx", 4))
SOLID = float(args.get("solid", 0.0))   # brick E (0 = element-less slave)
NOCONTACT = [False]
MFREE = int(args.get("mfree", 0))
ALPHA = math.radians(float(args.get("alpha", 20.0)))
kd = 1.0e3
GEOM = {}


def build():
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    ta = math.tan(ALPHA)
    MW = float(args.get("mw", 1.0)); YW = float(args.get("yw", 0.0))
    xm, ym = [-MW, 0.0, MW], [0.0 - YW, 1.0 + YW]
    mt, t = {}, 1
    for j, y in enumerate(ym):
        for i, x in enumerate(xm):
            ops.node(t, x, y, -ta * abs(x)); mt[(i, j)] = t; t += 1
            if not MFREE: ops.fix(t - 1, 1, 1, 1)
    mq = []
    for i in range(2):
        mq += [mt[(i, 0)], mt[(i + 1, 0)], mt[(i + 1, 1)], mt[(i, 1)]]
    xs = [-1.0, -0.45, 0.0, 0.45, 1.0] if NX == 4 else list(np.linspace(-1, 1, NX + 1))
    ys = [0.0, 0.3, 0.75, 1.0]
    st, allslave, ridge = {}, [], []
    t = 101
    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            for cp in (("L", "R") if (SPLIT and abs(x) < 1e-14) else ("",)):
                ops.node(t, float(x), y, -ta * abs(x) - DELTA)
                st[(i, j, cp)] = t; allslave.append(t)
                if abs(x) < 1e-14:
                    ridge.append(t)
                t += 1

    def sn(i, j, side):
        return st.get((i, j, side), st.get((i, j, "")))
    sq = []
    for j in range(len(ys) - 1):
        for i in range(len(xs) - 1):
            side = "L" if xs[i + 1] <= 1e-14 else "R"
            sq += [sn(i, j, side), sn(i + 1, j, side), sn(i + 1, j + 1, side), sn(i, j + 1, side)]
    GEOM["mq"], GEOM["sq"] = list(mq), list(sq)
    ops.contactSurface(1, "-master", 4, *mq)
    ops.contactSurface(2, "-slave-segments", 4, *sq)
    opts = ["-mortar", "-epsN", EPS, "-epsT", EPS * ETR, "-outward", 0.0, 0.0, 1.0]
    if COH > 0:
        opts += ["-cohesion", COH]
    if MU > 0:
        opts += ["-mu", MU]
    opts += [f for f in FLAGS.split(",") if f]
    if not NOCONTACT[0]:
        ops.contact(1, 1, 2, *opts)
    surf = list(allslave)
    if MFREE:
        allslave.extend(sorted(set(mq)))
    if SOLID > 0:
        ops.nDMaterial("ElasticIsotropic", 7, SOLID, 0.25)
        bot = {}
        for s_ in surf:
            c = ops.nodeCoord(s_); bot[s_] = 3000 + s_
            ops.node(bot[s_], c[0], c[1], c[2] - 0.3); allslave.append(bot[s_])
        for q in range(0, len(sq), 4):
            a, b, c_, d = sq[q:q + 4]
            ops.element("stdBrick", 7000 + q, bot[a], bot[b], bot[c_], bot[d], a, b, c_, d, 7)
    ops.uniaxialMaterial("Elastic", 1, kd)
    for s in allslave:
        g = 50000 + s
        ops.node(g, *ops.nodeCoord(s)); ops.fix(g, 1, 1, 1)
        ops.element("zeroLength", 60000 + s, g, s, "-mat", 1, 1, 1, "-dir", 1, 2, 3)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1)
    for s in surf:
        ops.load(s, 0.0, 0.0, -PZ)
    ops.constraints("LadrunoContact"); ops.numberer("Plain"); ops.system("FullGeneral")
    ops.test("NormUnbalance", 1.0e-7, 60, 0); ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0 / NST); ops.analysis("Static")
    return allslave, ridge


def run_all():
    allslave, ridge = build()
    U, its = [], []
    for s in range(NST):
        rc = ops.analyze(1)
        its.append(ops.testIter())
        if rc != 0:
            its[-1] = -its[-1]
            break
        U.append({n: list(ops.nodeDisp(n)) for n in allslave})
    return U, its


def mats(allslave):
    neq = ops.systemSize()
    A = np.array(ops.printA("-ret")).reshape(neq, neq).T   # Matrix data is column-major
    return A, neq


def resid():
    return np.array(ops.printB("-ret"))


def probe(allslave, ridge, u_trial):
    for n in allslave:
        for d in range(3):
            ops.setNodeDisp(n, d + 1, u_trial[n][d], "-commit")
    A, neq = mats(allslave)
    R0 = resid()
    Kfd = np.zeros((neq, neq))
    eqmap = {}
    for n in allslave:
        eqs = ops.nodeDOFs(n)
        for d in range(3):
            e = eqs[d]
            eqmap[e] = (n, d)
            u = u_trial[n][d]
            ops.setNodeDisp(n, d + 1, u + H, "-commit"); Rp = resid()
            ops.setNodeDisp(n, d + 1, u - H, "-commit"); Rm = resid()
            ops.setNodeDisp(n, d + 1, u, "-commit")
            Kfd[:, e] = -(Rp - Rm) / (2 * H)
            Kfd[e, e] += kd   # zeroLength springs: setNodeDisp does not update() elements
    # re-evaluate the residual at the base to check idempotence
    R1 = resid()
    E = A - Kfd
    scale = np.abs(Kfd).max()
    print(f"  |K|max={scale:.4e}  |K-Kfd|max={np.abs(E).max():.4e} rel={np.abs(E).max()/scale:.3e}"
          f"  idem={np.abs(R1-R0).max():.2e}  asym(K)={np.abs(A-A.T).max()/scale:.2e}"
          f"  asym(Kfd)={np.abs(Kfd-Kfd.T).max()/scale:.2e}")
    # blocks by node
    rows = []
    for n in allslave:
        en = ops.nodeDOFs(n)[:3]
        blk_e = E[np.ix_(en, range(neq))]
        blk_k = Kfd[np.ix_(en, range(neq))]
        rows.append((np.abs(blk_e).max() / max(np.abs(blk_k).max(), 1e-300), n, n in ridge,
                     ops.nodeCoord(n)))
    rows.sort(reverse=True)
    for r in rows[:8]:
        print(f"    node {r[1]:4d} ridge={int(r[2])} xyz=({r[3][0]:+.3f},{r[3][1]:.2f},{r[3][2]:+.4f})"
              f"  row-rel-err {r[0]:.3e}")
    # Newton contraction estimate from the tangent error: spectral radius of I - K^-1 Kfd
    try:
        G = np.eye(neq) - np.linalg.solve(A, Kfd)
        print(f"  rho(I - K^-1 Kfd) = {max(abs(np.linalg.eigvals(G))):.3e}")
    except Exception as ex:
        print("  rho failed", ex)
    return A, Kfd


if __name__ == "__main__":
    U, its = run_all()
    print(f"split={SPLIT} mu={MU} coh={COH} eps={EPS} etr={ETR} flags={FLAGS} iterations per step: {its}")
    if len(U) < K + 1 and K > 0:
        print("  not enough converged steps to probe; probing the last converged + trial guess")
    allslave, ridge = build()
    for s in range(K):
        ops.analyze(1)
    tgt = U[min(K, len(U) - 1)] if U else None
    print(f" probe at committed step {K} (stick-like trial):")
    probe(allslave, ridge, {n: list(ops.nodeDisp(n)) for n in allslave})
    if tgt is not None:
        print(f" probe at trial = converged step {K + 1} (slip trial):")
        probe(allslave, ridge, tgt)
