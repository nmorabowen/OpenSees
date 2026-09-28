"""Vectorized plane-strain acoustic-tensor scan of SANISAND's continuum elastoplastic
tangent D = Ce - (Ce:R)(x)(Q:Ce)/H at every committed GP (see acoustic.py for the
definitions). Rank-one update of the isotropic elastic acoustic tensor:
    A(m) = A_e(m) - a (x) b / H,  a = m.(Ce:R),  b = (Q:Ce).m
    det A / det A_e = 1 - b . A_e^{-1} a / H,  A_e^{-1} = (I - k m(x)m)/G,  k = (K+G/3)/(K+4G/3)
Reports the real (non-associated) tangent, an associated control (R -> Q), and the
ADR-90 V4 viscous blend (1-beta) Ce + beta D at a few beta.
"""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import h_decomp as hd

TH = np.linspace(0.0, math.pi, 721)
M = np.stack([np.cos(TH), np.sin(TH)], 1)            # (T,2)


def tens(v):                                         # (N,6) tensor-shear -> (N,3,3)
    t = np.zeros((len(v), 3, 3))
    t[:, 0, 0], t[:, 1, 1], t[:, 2, 2] = v[:, 0], v[:, 1], v[:, 2]
    t[:, 0, 1] = t[:, 1, 0] = v[:, 3]
    t[:, 1, 2] = t[:, 2, 1] = v[:, 4]
    t[:, 0, 2] = t[:, 2, 0] = v[:, 5]
    return t


def ratio(CeR, QCe, H, G, K):
    """min over theta of det A / det A_e ; returns (min, theta_deg)."""
    a = np.einsum("ti,nij->ntj", M, CeR[:, :2, :2])
    b = np.einsum("nkl,tl->ntk", QCe[:, :2, :2], M)
    k = ((K + G / 3.0) / (K + 4.0 * G / 3.0))[:, None]
    ba = np.einsum("ntk,ntk->nt", b, a)
    bm = np.einsum("ntk,tk->nt", b, M)
    am = np.einsum("ntk,tk->nt", a, M)
    q = (ba - k * bm * am) / G[:, None]
    r = 1.0 - q / H[:, None]
    i = np.argmin(r, axis=1)
    return r[np.arange(len(r)), i], np.degrees(TH[i])


def run(npz):
    R, sb = hd.decompose(npz)
    d = np.load(npz, allow_pickle=True)
    S = -d["s6"]; st = d["st"]
    p = np.maximum(S[:, :3].sum(1) / 3.0, hd.SMALL)
    dev = S.copy(); dev[:, :3] -= p[:, None]
    al = st[:, 6:12]
    nv = dev - p[:, None] * al
    nn = np.sqrt(nv[:, 0]**2 + nv[:, 1]**2 + nv[:, 2]**2 + 2 * (nv[:, 3]**2 + nv[:, 4]**2 + nv[:, 5]**2))
    nv = nv / np.where(nn >= hd.SMALL, nn, 1.0)[:, None]
    N = tens(nv); A = tens(al); I = np.eye(3)[None]
    N3 = np.einsum("nij,njk,nkl->nil", N, N, N)
    c3 = np.clip(math.sqrt(6) * np.trace(N3, axis1=1, axis2=2), -1, 1)
    g = 2 * hd.CC / ((1 + hd.CC) - (1 - hd.CC) * c3)
    B = 1 + 1.5 * (1 - hd.CC) / hd.CC * g * c3
    C = 3 * math.sqrt(1.5) * (1 - hd.CC) / hd.CC * g
    D = R[:, 12]; G = R[:, 16]; K = 2.0 / 3.0 * (1 + hd.NU) / (1 - 2 * hd.NU) * G
    NN = np.einsum("nij,njk->nik", N, N)
    Rt = B[:, None, None] * N - C[:, None, None] * (NN - I / 3.0) + (D / 3.0)[:, None, None] * I
    qv = np.einsum("nij,nij->n", N, A) + hd.R23 * hd.MM
    Qt = N - (qv / 3.0)[:, None, None] * I

    def ce(T):
        tr = np.trace(T, axis1=1, axis2=2)
        devT = T - (tr / 3.0)[:, None, None] * I
        return 2 * G[:, None, None] * devT + (K * tr)[:, None, None] * I

    CeR, CeQ = ce(Rt), ce(Qt)
    H = R[:, 11]; Kp = R[:, 8]
    Ha = Kp + np.einsum("nij,nij->n", Qt, CeQ)
    r, th = ratio(CeR, CeQ, H, G, K)
    ra, _ = ratio(CeQ, CeQ, Ha, G, K)
    # viscous blend (1-beta)Ce + beta*D = Ce - beta CeR (x) QCe / H  -> H/beta
    rb = {beta: ratio(CeR, CeQ, H / beta, G, K)[0] for beta in (0.5, 0.9, 0.99)}
    return dict(sb=sb, x=d["gx"], y=d["gy"], p=p, psi=d["psi"], H2G=H / (2 * G), r=r, th=th, c3=c3,
                ra=ra, rb=rb, Kp2G=Kp / (2 * G), D=D, f=d["f"])


if __name__ == "__main__":
    for f in sys.argv[1].split(","):
        o = run(f)
        tag = os.path.basename(os.path.dirname(os.path.dirname(f))) + "/" + os.path.basename(f)
        r, ra = o["r"], o["ra"]
        near = (np.abs(np.abs(o["x"]) - 0.80) < 0.20) & (o["y"] > -3.2)
        print(f"\n== {tag}  s/B {o['sb']:.4f}  GPs {len(r)}", flush=True)
        print(f"   LOSS OF ELLIPTICITY (det<=0): real {int((r<=0).sum())}  [edge zone {int((r[near]<=0).sum())}]"
              f"   associated control {int((ra<=0).sum())}   blend beta=0.5/0.9/0.99: "
              + " / ".join(str(int((o['rb'][b] <= 0).sum())) for b in (0.5, 0.9, 0.99)))
        print(f"   min det ratio: real {r.min():+.3e}  assoc {ra.min():+.3e};  real quantiles 0.1/1/5 %: "
              f"{np.quantile(r,0.001):+.4f} {np.quantile(r,0.01):+.4f} {np.quantile(r,0.05):+.4f}")
        bad = np.where(r <= 0)[0]
        if len(bad):
            print(f"   det<=0 GPs: x in [{o['x'][bad].min():+.2f},{o['x'][bad].max():+.2f}] "
                  f"y in [{o['y'][bad].min():+.2f},{o['y'][bad].max():+.2f}]  H/2G median {np.median(o['H2G'][bad]):+.3f} "
                  f"(min {o['H2G'][bad].min():+.3f})  Kp/2G median {np.median(o['Kp2G'][bad]):+.4f}  "
                  f"band angle median {np.median(o['th'][bad]):.1f} deg  D median {np.median(o['D'][bad]):+.3f}")
        for i in np.argsort(r)[:5]:
            print(f"     ({o['x'][i]:+.3f},{o['y'][i]:+.3f}) p {o['p'][i]:6.1f} psi {o['psi'][i]:+.3f} "
                  f"H/2G {o['H2G'][i]:+.3f} Kp/2G {o['Kp2G'][i]:+.4f} D {o['D'][i]:+.3f} det {r[i]:+.3e} @ {o['th'][i]:5.1f} "
                  f"assoc {ra[i]:+.3e}", flush=True)
