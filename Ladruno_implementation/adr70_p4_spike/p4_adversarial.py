# Independent adversarial re-derivation for ADR-70 P4 claims.
# I rebuild the T6 machinery from scratch (do NOT import the spike) so a shared
# bug can't hide.  I add: distorted/irregular pairs (break symmetry), a 4-elem
# patch (over-relaxation trend), a CURVED (non-affine) T6 (probe Claims B & D),
# and a pressure-mode / inf-sup probe on a structured mesh of pairs.
import numpy as np
np.set_printoptions(precision=4, suppress=True, linewidth=140)

# ---- T6 quadratic shape functions & derivatives (area coords) ----
def t6_N(xi, eta):
    L1, L2, L3 = 1-xi-eta, xi, eta
    return np.array([L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1),
                     4*L1*L2, 4*L2*L3, 4*L3*L1])

def t6_dN(xi, eta):
    L1, L2, L3 = 1-xi-eta, xi, eta
    dL = np.array([[-1.,-1.],[1.,0.],[0.,1.]])
    dN = np.zeros((6,2)); L=[L1,L2,L3]
    for i in range(3):
        dN[i] = (4*L[i]-1)*dL[i]
    dN[3] = 4*(L2*dL[0]+L1*dL[1])
    dN[4] = 4*(L3*dL[1]+L2*dL[2])
    dN[5] = 4*(L1*dL[2]+L3*dL[0])
    return dN

GP3 = [(2/3,1/6),(1/6,2/3),(1/6,1/6)]; W3=[1/6,1/6,1/6]
M = np.array([1.,1.,0.])

def elem_rows(xy6):
    """returns list of (B 3x12, b 12, w*detJ, gp_phys 2) per GP. Independent build."""
    out=[]
    for (xi,eta),w in zip(GP3,W3):
        dN=t6_dN(xi,eta); J=xy6.T@dN; detJ=np.linalg.det(J)
        if detJ<=0: raise RuntimeError("winding/degenerate detJ<=0")
        dNx=dN@np.linalg.inv(J)
        B=np.zeros((3,12))
        for a in range(6):
            B[0,2*a]=dNx[a,0]; B[1,2*a+1]=dNx[a,1]
            B[2,2*a]=dNx[a,1]; B[2,2*a+1]=dNx[a,0]
        gpx = t6_N(xi,eta)@xy6
        out.append((B,B[0]+B[1],w*detJ,gpx))
    return out

def t6_from_corners(corners, mid_perturb=None):
    c=np.asarray(corners,float)
    mids=np.array([(c[0]+c[1])/2,(c[1]+c[2])/2,(c[2]+c[0])/2])
    if mid_perturb is not None:
        mids = mids + np.asarray(mid_perturb,float)   # curved (non-affine) edges
    return np.vstack([c,mids])

def assemble_and_rank(elems, dofmaps, ndof, projector, tol=1e-10, return_S=False):
    """elems[e]=elem_rows list; dofmaps[e]=global dof list (len 12)."""
    all_gp=[]
    for e,rows in enumerate(elems):
        for (B,b,wd,gx) in rows:
            Bg=np.zeros((3,ndof)); bg=np.zeros(ndof)
            for a in range(6):
                ga=dofmaps[e][2*a]//2
                Bg[:,2*ga:2*ga+2]=B[:,2*a:2*a+2]
                bg[2*ga:2*ga+2]=b[2*a:2*a+2]
            all_gp.append((e,Bg,bg,wd,gx))
    btil=projector(all_gp)
    S=np.vstack([np.sqrt(wd)*(Bg+0.5*np.outer(M,btil[k]-bg))
                 for k,(e,Bg,bg,wd,gx) in enumerate(all_gp)])
    sv=np.linalg.svd(S,compute_uv=False)
    nz=ndof-int((sv>tol*sv[0]).sum())
    return (nz,S,sv) if return_S else nz

def proj_none(a):          return [bg for (_,_,bg,_,_) in a]
def proj_const_elem(a):
    out=[]
    for (e0,_,_,_,_) in a:
        num=sum(wd*bg for (e,_,bg,wd,_) in a if e==e0); den=sum(wd for (e,_,_,wd,_) in a if e==e0)
        out.append(num/den)
    return out
def proj_const_patch(a):
    num=sum(wd*bg for (_,_,bg,wd,_) in a); den=sum(wd for (_,_,_,wd,_) in a)
    return [num/den]*len(a)
def proj_P1_patch(a):
    P=np.array([[1.,gx[0],gx[1]] for (_,_,_,_,gx) in a])
    Wd=np.diag([wd for (_,_,_,wd,_) in a]); Bm=np.vstack([bg for (_,_,bg,_,_) in a])
    coef=np.linalg.solve(P.T@Wd@P, P.T@Wd@Bm); Bt=P@coef
    return [Bt[k] for k in range(len(a))]
def proj_P2_patch(a):   # richer projection: P2 over the patch (6 dof) -- to bracket
    P=np.array([[1.,gx[0],gx[1],gx[0]**2,gx[0]*gx[1],gx[1]**2] for (_,_,_,_,gx) in a])
    Wd=np.diag([wd for (_,_,_,wd,_) in a]); Bm=np.vstack([bg for (_,_,bg,_,_) in a])
    coef=np.linalg.lstsq(P.T@Wd@P, P.T@Wd@Bm, rcond=None)[0]; Bt=P@coef
    return [Bt[k] for k in range(len(a))]

# global node dedup builder
def build(elem_corner_list, mid_perturbs=None):
    nodes=[]
    def nid(p):
        for i,q in enumerate(nodes):
            if np.allclose(p,q,atol=1e-12): return i
        nodes.append(np.array(p,float)); return len(nodes)-1
    elems=[]; dofmaps=[]
    for i,corners in enumerate(elem_corner_list):
        mp = None if mid_perturbs is None else mid_perturbs[i]
        xy6=t6_from_corners(corners, mp)
        maps=[nid(p) for p in xy6]
        elems.append(elem_rows(xy6))
        dofmaps.append([2*n+d for n in maps for d in (0,1)])
    return elems, dofmaps, 2*len(nodes), nodes

print("="*78)
print("REPRO on symmetric unit-square split (independent rebuild)")
A,Bp,C,D=(0,0),(1,0),(1,1),(0,1)
elems,dm,nd,_=build([[A,Bp,C],[A,C,D]])
for nm,pj,exp in [("std",proj_none,3),("const-elem",proj_const_elem,4),
                  ("patch-const (Claim A)",proj_const_patch,5),
                  ("patch-P1 (Claim D)",proj_P1_patch,3),
                  ("patch-P2",proj_P2_patch,None)]:
    print(f"  {nm:24s}: zero modes = {assemble_and_rank(elems,dm,nd,pj)}  (ref {exp})")

print("="*78)
print("ATTACK 1: DISTORTED / IRREGULAR pair (symmetry cannot hide rank)")
# irregular quad split along a diagonal, no symmetry at all
P0,P1_,P2_,P3_=(0,0),(1.3,-0.2),(1.6,1.1),(-0.15,0.9)
elems,dm,nd,_=build([[P0,P1_,P2_],[P0,P2_,P3_]])
for nm,pj in [("std",proj_none),("patch-const",proj_const_patch),
              ("patch-P1",proj_P1_patch),("patch-P2",proj_P2_patch)]:
    print(f"  {nm:24s}: zero modes = {assemble_and_rank(elems,dm,nd,pj)}")

print("="*78)
print("ATTACK 1b: another distorted pair, opposite diagonal + skew")
Q0,Q1,Q2,Q3=(0.1,0.0),(2.0,0.3),(2.3,1.7),(0.4,1.2)
elems,dm,nd,_=build([[Q0,Q1,Q3],[Q1,Q2,Q3]])   # split along Q1-Q3
for nm,pj in [("std",proj_none),("patch-const",proj_const_patch),
              ("patch-P1",proj_P1_patch)]:
    print(f"  {nm:24s}: zero modes = {assemble_and_rank(elems,dm,nd,pj)}")

print("="*78)
print("ATTACK 2: 4-element patch (over-relaxation trend)")
# central node fan: square with center, 4 triangles sharing center vertex
cen=(0.5,0.5)
tris=[[A,Bp,cen],[Bp,C,cen],[C,D,cen],[D,A,cen]]
elems,dm,nd,_=build(tris)
for nm,pj in [("std",proj_none),("patch-const",proj_const_patch),
              ("patch-P1",proj_P1_patch),("patch-P2",proj_P2_patch)]:
    print(f"  {nm:24s}: zero modes = {assemble_and_rank(elems,dm,nd,pj)}  (dof={nd})")
print("  NOTE: RBM=3. patch-const excess = spurious count for a 4-tri patch.")

print("="*78)
print("ATTACK 3: CURVED (non-affine) T6 -- Claim B identity & Claim D rank")
# perturb midside nodes off the straight edges -> detJ varies over element
pert1=[[0.06,-0.05],[-0.04,0.07],[0.05,0.05]]
pert2=[[0.03,0.04],[-0.06,-0.03],[0.05,-0.04]]
elems,dm,nd,_=build([[A,Bp,C],[A,C,D]],mid_perturbs=[pert1,pert2])
# Claim B: element-local P1 projection identity at GPs, on a CURVED element
r1=elems[0]
P=np.array([[1.,gx[0],gx[1]] for (_,_,_,gx) in r1])
Wd=np.diag([wd for (_,_,wd,_) in r1]); Bm=np.vstack([b for (_,b,_,_) in r1])
coef=np.linalg.solve(P.T@Wd@P,P.T@Wd@Bm); Bt=P@coef
print(f"  Claim B on CURVED T6: max|b_tilde-b| at GPs = {np.max(np.abs(Bt-Bm)):.2e}")
for nm,pj in [("std",proj_none),("patch-const",proj_const_patch),
              ("patch-P1",proj_P1_patch),("patch-P2",proj_P2_patch)]:
    print(f"  {nm:24s}: zero modes = {assemble_and_rank(elems,dm,nd,pj)}")
print("  (curved-edge strains are RATIONAL, not linear -> 3-GP dev sampling may alias)")

print("="*78)
print("ATTACK 3b: explicit conformal mode injection on distorted pair")
# Build u = a*(z-z_p)^2 nodally on the distorted pair and confirm zero energy
# under patch-const, NONzero under patch-P1.  z_p = patch area centroid.
P0,P1_,P2_,P3_=(0,0),(1.3,-0.2),(1.6,1.1),(-0.15,0.9)
elems,dm,nd,nodes=build([[P0,P1_,P2_],[P0,P2_,P3_]])
# patch area centroid (area-weighted over the two triangles)
def tri_area_cent(c):
    c=np.asarray(c,float); ar=0.5*np.linalg.det(np.array([c[1]-c[0],c[2]-c[0]]))
    return ar, c.mean(0)
a1,g1=tri_area_cent([P0,P1_,P2_]); a2,g2=tri_area_cent([P0,P2_,P3_])
zc=(a1*complex(*g1)+a2*complex(*g2))/(a1+a2)
def conformal_disp(nodes, a, zc):
    u=np.zeros(2*len(nodes))
    for i,p in enumerate(nodes):
        z=complex(p[0],p[1]); w=a*(z-zc)**2
        u[2*i]=w.real; u[2*i+1]=w.imag
    return u
def energy(S,u):
    e=S@u; return float(e@e)
for tag,a in [("a=1 (real)",1+0j),("a=i (imag)",1j),("a=1+i",1+1j)]:
    u=conformal_disp(nodes,a,zc)
    _,Sc,_=assemble_and_rank(elems,dm,nd,proj_const_patch,return_S=True)
    _,Sp,_=assemble_and_rank(elems,dm,nd,proj_P1_patch,return_S=True)
    print(f"  conformal {tag:12s}: E(patch-const)={energy(Sc,u):.3e}  E(patch-P1)={energy(Sp,u):.3e}")

print("="*78)
print("ATTACK 4: pressure-mode / inf-sup probe on structured 2x2 mesh of pairs")
# Build 2x2 grid of unit squares, each split A,B,C / A,C,D -> 8 T6 triangles,
# 4 macro-patches. Equivalent mixed pair: P2 velocity / P1-cont-per-macro pressure.
# Assemble G[p,u]=int q*div(u).  Spurious pressure modes = nullspace of G^T beyond
# the single global-constant (which is spurious on a fully-free mesh).
def build_grid(ncell):
    # returns element corner list + which macro each pair belongs to
    ecl=[]; macro=[]
    for i in range(ncell):
        for j in range(ncell):
            x0,y0=i,j
            a=(x0,y0); b=(x0+1,y0); c=(x0+1,y0+1); d=(x0,y0+1)
            ecl.append([a,b,c]); ecl.append([a,c,d])
            mid=len(macro)//1
            macro.append(len(ecl)//2-1 if False else (i*ncell+j))
            macro.append(i*ncell+j)
    return ecl,macro
ncell=2
ecl,macro=build_grid(ncell)
elems,dm,nd,nodes=build(ecl)
npatch=ncell*ncell
# pressure dofs: 3 per macro (1,x,y continuous within macro), discontinuous across
# G row block for pressure of macro mp = sum over its 2 elems, GPs: q_basis * b (div)
Gp_rows=[]  # (npatch*3) x nd
Grows=np.zeros((npatch*3,nd))
# recompute gp lists with element->macro
gp_all=[]
idx=0
for e,rows in enumerate(elems):
    for (B,b,wd,gx) in rows:
        bg=np.zeros(nd)
        for a_ in range(6):
            ga=dm[e][2*a_]//2; bg[2*ga:2*ga+2]=b[2*a_:2*a_+2]
        gp_all.append((macro[e],bg,wd,gx))
for (mp,bg,wd,gx) in gp_all:
    q=np.array([1.,gx[0],gx[1]])
    for r in range(3):
        Grows[3*mp+r]+= wd*q[r]*bg
# inf-sup: singular values of G (pressure x velocity). Count pressure null modes.
svG=np.linalg.svd(Grows,compute_uv=False)
ptol=1e-9*svG[0]
p_null=int((svG<ptol).sum())
print(f"  mesh {ncell}x{ncell}: velocity dof={nd}, pressure dof={npatch*3}={Grows.shape[0]}")
print(f"  G singular values (sorted):\n   {np.round(svG,4)}")
print(f"  # pressure null modes (spurious incl. global const) = {p_null}")
print(f"  smallest NONZERO sing val (inf-sup proxy) = {svG[svG>ptol].min():.4e}")

# refine to 3x3 to see if inf-sup proxy shrinks toward 0 (checkerboard) or holds
ncell=3
ecl,macro=build_grid(ncell)
elems,dm,nd,nodes=build(ecl)
npatch=ncell*ncell
Grows=np.zeros((npatch*3,nd)); gp_all=[]
for e,rows in enumerate(elems):
    for (B,b,wd,gx) in rows:
        bg=np.zeros(nd)
        for a_ in range(6):
            ga=dm[e][2*a_]//2; bg[2*ga:2*ga+2]=b[2*a_:2*a_+2]
        gp_all.append((macro[e],bg,wd,gx))
for (mp,bg,wd,gx) in gp_all:
    q=np.array([1.,gx[0],gx[1]])
    for r in range(3): Grows[3*mp+r]+=wd*q[r]*bg
svG=np.linalg.svd(Grows,compute_uv=False); ptol=1e-9*svG[0]
print(f"  mesh 3x3: pressure dof={npatch*3}, #pressure null={int((svG<ptol).sum())}, "
      f"min nonzero sv={svG[svG>ptol].min():.4e}")
