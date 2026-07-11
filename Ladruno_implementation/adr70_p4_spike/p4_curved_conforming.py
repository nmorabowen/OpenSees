import numpy as np
from p4_adversarial import (elem_rows, assemble_and_rank, proj_none,
                            proj_const_patch, proj_P1_patch, proj_P2_patch)
np.set_printoptions(precision=4, suppress=True, linewidth=140)

# CONFORMING curved pair: build 9 global nodes explicitly, share the diagonal
# midside node between both triangles so geometry is C0 (no crack).
# corners: A B C D ; T1=[A,B,C], T2=[A,C,D]; shared edge A-C with ONE midnode.
A,B,C,D=np.array([0.,0]),np.array([1.,0]),np.array([1.,1]),np.array([0.,1])
def midcurve(p,q,off):  # curved midside = midpoint + offset
    return 0.5*(p+q)+np.array(off)

def build_conforming(offsets):
    """offsets: dict of edge-> (dx,dy) for the 5 unique edges:
       AB, BC, CA(shared), CD, DA.  CA is shared by both triangles."""
    mAB=midcurve(A,B,offsets['AB']); mBC=midcurve(B,C,offsets['BC'])
    mCA=midcurve(C,A,offsets['CA']); mCD=midcurve(C,D,offsets['CD'])
    mDA=midcurve(D,A,offsets['DA'])
    # T1 = [A,B,C] with mids (AB, BC, CA)
    xy1=np.vstack([A,B,C, mAB,mBC,mCA])
    # T2 = [A,C,D] with mids (AC=CA, CD, DA). NOTE edge A-C midnode == mCA (shared!)
    xy2=np.vstack([A,C,D, mCA,mCD,mDA])
    # global node table
    nodes=[A,B,C,D,mAB,mBC,mCA,mCD,mDA]
    idx={'A':0,'B':1,'C':2,'D':3,'mAB':4,'mBC':5,'mCA':6,'mCD':7,'mDA':8}
    # T1 nodes: A B C mAB mBC mCA  -> global 0 1 2 4 5 6
    map1=[0,1,2,4,5,6]
    # T2 nodes: A C D mCA mCD mDA  -> global 0 2 3 6 7 8
    map2=[0,2,3,6,7,8]
    e1=elem_rows(xy1); e2=elem_rows(xy2)
    dof1=[2*n+d for n in map1 for d in (0,1)]
    dof2=[2*n+d for n in map2 for d in (0,1)]
    return [e1,e2],[dof1,dof2],2*len(nodes)

print("CONFORMING curved pair (shared diagonal midnode identical in both tris)")
print(f"{'case':>28} | std | pc | pP1 | pP2")
# (1) straight reference
off0={k:(0,0) for k in ['AB','BC','CA','CD','DA']}
elems,dm,nd=build_conforming(off0)
print(f"{'straight (ref)':>28} | {assemble_and_rank(elems,dm,nd,proj_none):3d} |"
      f" {assemble_and_rank(elems,dm,nd,proj_const_patch):2d} |"
      f" {assemble_and_rank(elems,dm,nd,proj_P1_patch):3d} |"
      f" {assemble_and_rank(elems,dm,nd,proj_P2_patch):3d}")
# (2) curve ONLY the shared diagonal edge (bulges the interior interface)
for tag,off in [
    ("shared-edge CA curved only", {**off0,'CA':(0.08,0.08)}),
    ("all outer edges curved",     {**off0,'AB':(0,-0.07),'BC':(0.07,0),'CD':(0,0.07),'DA':(-0.07,0)}),
    ("all 5 edges curved",         {'AB':(0.03,-0.06),'BC':(0.06,0.02),'CA':(0.05,0.04),
                                    'CD':(-0.04,0.06),'DA':(-0.06,-0.03)}),
]:
    elems,dm,nd=build_conforming(off)
    r_std=assemble_and_rank(elems,dm,nd,proj_none)
    r_pc =assemble_and_rank(elems,dm,nd,proj_const_patch)
    r_p1 =assemble_and_rank(elems,dm,nd,proj_P1_patch)
    r_p2 =assemble_and_rank(elems,dm,nd,proj_P2_patch)
    print(f"{tag:>28} | {r_std:3d} | {r_pc:2d} | {r_p1:3d} | {r_p2:3d}")
print()
print("If pP1 stays 3 with conforming curvature -> my earlier '5' was a mesh CRACK")
print("(different shared-midnode perturbation per triangle), and Claim D SURVIVES.")
print("If pP1 = 5 even conforming -> Claim D genuinely fails on curved T6.")
