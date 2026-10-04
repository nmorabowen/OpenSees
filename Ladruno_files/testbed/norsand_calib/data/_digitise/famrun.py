import json
from gp import *
from curate import despike
PALETTE=[(255,0,0),(0,150,0),(0,0,255),(200,0,200),(255,140,0),(0,170,170),(120,80,0),(90,90,90)]
def family(P, strip_seeds, chain_seeds, xs, ylim, ranges, xmerge=0.12, despike_thr=0.10, base_tol=0.09, max_hold=6, order=None, chain_area=(40,150), maxthick=17, extra=None):
    """returns {name:[(x,y,src)]} merging marker-chain points with column-strip points."""
    names=list(strip_seeds)
    S=track_family(P, strip_seeds, xs, ylim, base_tol=base_tol, max_hold=max_hold, maxthick=maxthick)
    bl=blobs_px(P,(P.x0,P.x1),ylim)
    C=chain_family(P, bl, chain_seeds, (P.x0,P.x1), ylim, single=chain_area, order=order or list(chain_seeds)) if chain_seeds else {}
    out={}
    for n in names:
        lo,hi=ranges.get(n,(P.x0,P.x1))
        ch=[]
        for cn,cv in C.items():
            if cn==n or cn.startswith(n+'#'):
                ch+=[(x,y,'chain') for (x,y,a,f) in cv if lo<=x<=hi]
        if extra and n in extra: ch+=[(x,y,'row') for (x,y) in extra[n] if lo<=x<=hi]
        st=[(x,y,'strip') for x,(y,s,u) in S[n].items() if s in('ok','seed') and lo<=x<=hi]
        pts=list(ch)
        for (x,y,s) in st:
            if all(abs(x-c[0])>xmerge for c in ch): pts.append((x,y,s))
        pts.sort()
        kept=despike([(x,y) for x,y,s in pts],thr=despike_thr)
        ks=set(kept)
        out[n]=[(x,y,s) for x,y,s in pts if (x,y) in ks]
    return out

def follow_isolated(P, xs, y0, tol, ylim, maxthick=40, halfw=2):
    """follow one isolated curve over the ascending abscissae xs starting from guess y0: nearest cluster within tol of the previous value"""
    out=[]; yprev=y0
    for x in xs:
        cl=[c for c in clusters2(P,x,ylim,halfw=halfw) if c[1]<=maxthick]
        if not cl: continue
        c=min(cl,key=lambda c:abs(c[0]-yprev))
        if abs(c[0]-yprev)<=tol:
            out.append((x,c[0])); yprev=c[0]
    return out
