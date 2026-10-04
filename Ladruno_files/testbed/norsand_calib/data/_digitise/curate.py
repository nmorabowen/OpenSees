import numpy as np
def despike(pts, win=2, thr=0.10, protect_x=0.0):
    """pts: sorted [(x,y)]. drop points deviating more than thr from the median of their 2*win neighbours' linear interpolation."""
    pts=sorted(pts); keep=[]
    xs=np.array([p[0] for p in pts]); ys=np.array([p[1] for p in pts])
    for i,(x,y) in enumerate(pts):
        if x<protect_x: keep.append((x,y)); continue
        idx=[j for j in range(max(0,i-win),min(len(pts),i+win+1)) if j!=i]
        if len(idx)<3: keep.append((x,y)); continue
        # local linear fit through neighbours
        A=np.vstack([xs[idx],np.ones(len(idx))]).T
        m,c=np.linalg.lstsq(A,ys[idx],rcond=None)[0]
        if abs(y-(m*x+c))<=thr: keep.append((x,y))
    return keep
