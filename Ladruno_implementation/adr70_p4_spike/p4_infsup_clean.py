import numpy as np
from p4_adversarial import t6_from_corners, build
from p4_adversarial import t6_dN  # noqa
np.set_printoptions(precision=4, suppress=True, linewidth=140)
from p4_adversarial import t6_N, GP3, W3

# Clean discrete inf-sup for P2-velocity / P1-continuous-per-macro-pressure.
# Homogeneous Dirichlet velocity on the outer boundary (K SPD on interior dofs),
# mean-zero pressure (remove global constant). beta_h^2 = min eig of
#   Qm^{-1/2} B Kint^{-1} B^T Qm^{-1/2}   restricted to mean-zero pressure.
def run(ncell, curved=False, seed=3):
    ecl=[];macro=[]
    for i in range(ncell):
        for j in range(ncell):
            a=(i,j);b=(i+1,j);c=(i+1,j+1);d=(i,j+1)
            ecl.append([a,b,c]);macro.append(i*ncell+j)
            ecl.append([a,c,d]);macro.append(i*ncell+j)
    mp=None
    if curved:
        rng=np.random.default_rng(seed)
        mp=[(rng.random((3,2))-0.5)*0.1 for _ in ecl]
    elems,dm,nd,nodes=build(ecl,mid_perturbs=mp)
    npat=ncell*ncell
    K=np.zeros((nd,nd)); B=np.zeros((3*npat,nd)); Q=np.zeros((3*npat,3*npat))
    for e,corners in enumerate(ecl):
        xy6=t6_from_corners(corners, None if mp is None else mp[e]); loc=dm[e]
        for (xi,eta),w in zip(GP3,W3):
            dN=t6_dN(xi,eta);J=xy6.T@dN;detJ=np.linalg.det(J);dNx=dN@np.linalg.inv(J)
            Gmat=np.zeros((4,12)); div=np.zeros(12)
            for aa in range(6):
                Gmat[0,2*aa]=dNx[aa,0];Gmat[1,2*aa]=dNx[aa,1]
                Gmat[2,2*aa+1]=dNx[aa,0];Gmat[3,2*aa+1]=dNx[aa,1]
                div[2*aa]=dNx[aa,0];div[2*aa+1]=dNx[aa,1]
            gx=t6_N(xi,eta)@xy6;q=np.array([1.,gx[0],gx[1]]);wd=w*detJ;m=macro[e]
            Gg=np.zeros((4,nd));dg=np.zeros(nd)
            for aa in range(6):
                ga=loc[2*aa]//2
                Gg[:,2*ga:2*ga+2]=Gmat[:,2*aa:2*aa+2]; dg[2*ga:2*ga+2]=div[2*aa:2*aa+2]
            K+=wd*(Gg.T@Gg)
            for r in range(3):
                B[3*m+r]+=wd*q[r]*dg
                for s in range(3): Q[3*m+r,3*m+s]+=wd*q[r]*q[s]
    # interior velocity dofs: nodes not on the outer boundary
    bnd=lambda p:(abs(p[0]-0)<1e-9 or abs(p[0]-ncell)<1e-9 or abs(p[1]-0)<1e-9 or abs(p[1]-ncell)<1e-9)
    free=[i for i,p in enumerate(nodes) if not bnd(p)]
    fdof=[2*i+d for i in free for d in (0,1)]
    Kii=K[np.ix_(fdof,fdof)]; Bi=B[:,fdof]
    # mean-zero pressure projector: remove global constant (vector of the const pressure)
    # const pressure coefficient vector: q=1 on every macro -> coeff (1,0,0) per macro
    c=np.zeros(3*npat)
    for m in range(npat): c[3*m]=1.0
    # Q-orthonormal complement of span{c}
    Qc=Q@c; cn=c/np.sqrt(c@Qc)
    # Schur: S = Bi Kii^{-1} Bi^T  (pressure space)
    Kii_i=np.linalg.inv(Kii)
    S=Bi@Kii_i@Bi.T
    # generalized eig S x = lam Q x, restrict to mean-zero (project out cn)
    # transform to standard: L=chol(Q); solve L^{-1} S L^{-T}
    L=np.linalg.cholesky(Q); Linv=np.linalg.inv(L)
    A=Linv@S@Linv.T; A=0.5*(A+A.T)
    # remove constant-pressure direction in transformed coords: d = L^T c (unit)
    d=L.T@c; d=d/np.linalg.norm(d)
    P=np.eye(3*npat)-np.outer(d,d)
    A=P@A@P
    ev=np.sort(np.linalg.eigvalsh(A).real)
    ev=ev[ev>1e-12*ev[-1]]   # drop the projected-out zero
    return ev.min(), len(ev), 3*npat
print("Clean inf-sup (Dirichlet velocity, mean-zero pressure) -- STRAIGHT meshes")
print(f"{'N':>3} {'pres_dof':>8} {'#pos_eig':>8} {'beta_h^2':>12} {'beta_h':>10}")
prev=None
for N in [2,3,4,5]:
    b2,npos,pd=run(N,curved=False)
    print(f"{N:3d} {pd:8d} {npos:8d} {b2:12.4e} {np.sqrt(b2):10.4e}"
          + ("" if prev is None else f"   ratio={b2/prev:.3f}"))
    prev=b2
print()
print("Same, CURVED meshes (10% midside perturbation):")
prev=None
for N in [2,3,4]:
    b2,npos,pd=run(N,curved=True)
    print(f"{N:3d} {pd:8d} {npos:8d} {b2:12.4e} {np.sqrt(b2):10.4e}"
          + ("" if prev is None else f"   ratio={b2/prev:.3f}"))
    prev=b2
