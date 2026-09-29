"""Decompose the SAS-ME stage-1 plastic denominator
    H = Kp + 2G (B - C tr n^3) - K D qv,   Kp = 2/3 p h b:n,  qv = n:alpha + sqrt(2/3) m
at every Gauss point of a committed WP-138 footing checkpoint (read-only on the
orchestrator's analysis copy). Mirrors LadrunoSANISANDSasME.cpp:335-383 and
ManzariDafalias::GetStateDependent / GetElasticModuli (e_G = e_init, p_r = 0).
H at stage 1 depends ONLY on the committed state, so a GP with H <= 0 refuses every
loading (N > 0) increment -- the loadingNonPosH wall.
"""
import math, sys, glob, os
import numpy as np

# campaign set (footing_ab.py SAN)
G0, NU, EINIT, MC, CC, LAMC, E0, XI, PATM, MM = 264.32, 0.312885, 0.6944, 1.3309, 0.71, 0.027, 0.83, 0.45, 101.0, 0.005
H0, CH, NB, A0, ND, ZMAX, CZ = 1.3, 0.968, 3.5, 0.05, 5.75, 12.5, 1100.0
PMIN, SMALL = 0.0101, 1e-10
R23 = math.sqrt(2.0 / 3.0)
over = dict(a.split("=") for a in sys.argv[2:])
for k, v in over.items():
    globals()[k] = float(v)


def mat(v):
    return np.array([[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]])


def ddot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + 2*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5])


def decompose(npz):
    d = np.load(npz, allow_pickle=True)
    st, s6, e, psi = d["st"], d["s6"], d["st"][:, 24], d["psi"]
    S = -s6                                    # compression positive, tensor shear
    out = []
    for k in range(len(e)):
        Sk = S[k]; al = st[k, 6:12]; z = st[k, 12:18]; ain = st[k, 18:24]
        p = (Sk[0] + Sk[1] + Sk[2]) / 3.0
        p = max(p, SMALL)
        dev = Sk.copy(); dev[:3] -= p
        n = dev - p * al
        nn = math.sqrt(max(ddot(n, n), 0.0)); n = n / (nn if nn >= SMALL else 1.0)
        Nm = mat(n)
        c3 = float(np.clip(math.sqrt(6.0) * np.trace(Nm @ Nm @ Nm), -1, 1))
        g = 2 * CC / ((1 + CC) - (1 - CC) * c3)
        aB = g * MC * math.exp(-NB * psi[k]) - MM
        aD = g * MC * math.exp(ND * psi[k]) - MM
        b0 = G0 * H0 * (1 - CH * e[k]) / math.sqrt(p / PATM)
        b = R23 * aB * n - al
        dd = R23 * aD * n - al
        x = ddot(al - ain, n)
        sentinel = abs(x) < SMALL or x < SMALL      # SAS bracket: x < small -> 1e10
        h = 1e10 if sentinel else b0 / x
        A = A0 * (1 + max(ddot(z, n), 0.0))
        D = A * ddot(dd, n)
        if p < 0.05 * PATM:
            D *= 1.0 / (1.0 + math.exp(7.6349 - 7.2713 * 101.0 / PATM * p))
        B = 1.0 + 1.5 * (1 - CC) / CC * g * c3
        C = 3.0 * math.sqrt(1.5) * (1 - CC) / CC * g
        pn = max(p, PMIN)
        G = G0 * PATM * (2.97 - EINIT) ** 2 / (1 + EINIT) * math.sqrt(pn / PATM)
        K = 2.0 / 3.0 * (1 + NU) / (1 - 2 * NU) * G
        qv = ddot(n, al) + R23 * MM
        bn = ddot(b, n)
        Kp = 2.0 / 3.0 * p * h * bn
        T2 = 2 * G * (B - C * np.trace(Nm @ Nm @ Nm))
        T3 = -K * D * qv
        out.append((k, d["gx"][k], d["gy"][k], p, psi[k], x, h, bn, Kp, T2, T3, Kp + T2 + T3,
                    D, ddot(z, n), A, sentinel, G))
    return np.array(out, dtype=float), float(d["s_over_B"])


if __name__ == "__main__":
    for f in sys.argv[1].split(","):
        R, sb = decompose(f)
        H = R[:, 11]; T2 = R[:, 9]
        name = os.path.basename(os.path.dirname(os.path.dirname(f))) + "/" + os.path.basename(f)
        print(f"\n== {name}  s/B={sb:.4f}  GPs={len(R)}")
        print(f"   H<=0: {int((H<=0).sum())}   H/T2<0.05: {int((H/T2<0.05).sum())}   "
              f"sentinel h: {int(R[:,15].sum())}   Kp<0: {int((R[:,8]<0).sum())}   D>0 (contract): {int((R[:,12]>0).sum())}")
        idx = np.argsort(H / T2)[:8]
        print("   lowest H/(2G B'):  k      x       y       p     psi    (a-ain):n     h        b:n       Kp/T2   T3/T2   H/T2    D     z:n   A   sent")
        for i in idx:
            r = R[i]
            print(f"   {int(r[0]):6d} {r[1]:+7.3f} {r[2]:+7.3f} {r[3]:7.1f} {r[4]:+.3f} {r[5]:+.2e} {r[6]:9.2e} {r[7]:+.2e} "
                  f"{r[8]/r[9]:+8.3f} {r[10]/r[9]:+6.3f} {r[11]/r[9]:+7.3f} {r[12]:+.3f} {r[13]:+6.2f} {r[14]:.3f} {int(r[15])}")
