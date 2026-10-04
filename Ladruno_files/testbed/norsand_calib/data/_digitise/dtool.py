"""Digitising helpers working on the native 1-bit page scans extracted from the PDFs (WP-144 P3 data pack, see data/README.md section 8).
Tool-kit only: scripts here were run from a scratch directory; paths are made relative to this folder (work/ is created by extract_native_images.py).
Python 3.11 + numpy, scipy, Pillow, PyMuPDF."""
import os, json
import numpy as np
from PIL import Image, ImageDraw, ImageFont
SP = os.environ.get("DIG_WORK") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "work")   # holds <paper>_native_<page>.png (extract_native_images.py) and the overlay PNGs
_cache = {}
def page(tag, n):
    k = (tag, n)
    if k not in _cache:
        _cache[k] = np.array(Image.open(f"{SP}/dig/{tag}_native_{n:02d}.png").convert("L")) < 128   # True = ink
    return _cache[k]
def prev2nat(tag, box):
    """preview(110dpi) box -> native box"""
    sx = 2848/909.0; sy = (4015 if tag=="ft" else 3984)/1286.4
    return (int(box[0]*sx), int(box[1]*sy), int(box[2]*sx), int(box[3]*sy))
def crop(tag, n, nbox):
    x0,y0,x1,y1 = nbox
    return page(tag, n)[y0:y1, x0:x1], (x0, y0)
def runs(v, minlen):
    """indices (start,end) of True-runs >= minlen in 1D bool array"""
    out=[]; s=None
    for i,x in enumerate(list(v)+[False]):
        if x and s is None: s=i
        if (not x) and s is not None:
            if i-s>=minlen: out.append((s,i-1))
            s=None
    return out
def lines(a, axis, minlen):
    """find long straight ink lines. axis=0 -> vertical lines (columns), axis=1 -> horizontal lines (rows).
    returns list of (center_index, start, end) merging adjacent rows/cols"""
    res=[]
    n = a.shape[1] if axis==0 else a.shape[0]
    cand=[]
    for i in range(n):
        v = a[:,i] if axis==0 else a[i,:]
        for (s,e) in runs(v,minlen): cand.append((i,s,e))
    # merge adjacent
    merged=[]
    for c in cand:
        if merged and c[0]-merged[-1][-1][0]<=1 and abs(c[1]-merged[-1][-1][1])<40: merged[-1].append(c)
        else: merged.append([c])
    return [(np.mean([c[0] for c in m]), min(c[1] for c in m), max(c[2] for c in m), len(m)) for m in merged]
def ticks(a, axis_kind, pos, rng, side, minlen=10, maxoff=40):
    """ticks along an axis line. axis_kind 'v': vertical axis line at column pos; ticks are horizontal ink stubs at rows.
    side=+1 stubs to the right of pos, -1 to the left. rng=(lo,hi) index range to scan. returns centers"""
    out=[]
    lo,hi=rng
    if axis_kind=='v':
        cols = range(int(pos)+side*3, int(pos)+side*(3+minlen), side) if side>0 else range(int(pos)-3, int(pos)-3-minlen, -1)
        sub = np.zeros(hi-lo, bool)+True
        for c in cols: sub &= a[lo:hi, c]
        idx = np.where(sub)[0]+lo
    else:
        rows = range(int(pos)+side*3, int(pos)+side*(3+minlen), side) if side>0 else range(int(pos)-3, int(pos)-3-minlen, -1)
        sub = np.zeros(hi-lo, bool)+True
        for r in rows: sub &= a[r, lo:hi]
        idx = np.where(sub)[0]+lo
    g=[]
    for i in idx:
        if g and i-g[-1][-1]<=2: g[-1].append(i)
        else: g.append([i])
    return [float(np.mean(x)) for x in g]
def save_crop(a, path, scale=1.0, marks=None):
    im = Image.fromarray((~a).astype(np.uint8)*255).convert("RGB")
    if scale!=1.0: im = im.resize((int(im.width*scale), int(im.height*scale)), Image.LANCZOS)
    im.save(path); return im

def deskew(a, amax=1.5, step=0.05, refine=0.01):
    """rotate bool ink image to maximise row/col projection sharpness. returns (rotated bool, angle_deg)"""
    im = Image.fromarray(a.astype(np.uint8)*255)
    def score(ang):
        r = np.array(im.rotate(ang, resample=Image.BILINEAR, fillcolor=0)) > 100
        cs = r.sum(0).astype(float); rs = r.sum(1).astype(float)
        return (cs**2).sum() + (rs**2).sum()
    best=(-1,0)
    for ang in np.arange(-amax, amax+1e-9, step):
        s=score(ang)
        if s>best[0]: best=(s,ang)
    a0=best[1]
    for ang in np.arange(a0-step, a0+step+1e-9, refine):
        s=score(ang)
        if s>best[0]: best=(s,ang)
    ang=best[1]
    r = np.array(im.rotate(ang, resample=Image.BILINEAR, fillcolor=0)) > 100
    return r, ang

def fit_vline(a, xr, yr, minrun=60):
    """Robust fit x = m*y + c of a near-vertical ink line. xr=(x0,x1) column window, yr=(y0,y1) rows."""
    pts=[]
    for c in range(xr[0], xr[1]):
        for s,e in runs(a[yr[0]:yr[1], c], minrun):
            pts.append(((s+e)/2.0+yr[0], c))
    if len(pts)<3: return None
    P=np.array(pts,float)
    keep=np.ones(len(P),bool)
    for _ in range(6):
        m,cc=np.polyfit(P[keep,0],P[keep,1],1)
        res=np.abs(P[:,1]-(m*P[:,0]+cc)); keep=res<2.0
        if keep.sum()<3: break
    return m,cc,int(keep.sum())

def fit_hline(a, yr, xr, minrun=60):
    """Robust fit y = m*x + c of a near-horizontal ink line."""
    pts=[]
    for r in range(yr[0], yr[1]):
        for s,e in runs(a[r, xr[0]:xr[1]], minrun):
            pts.append(((s+e)/2.0+xr[0], r))
    if len(pts)<3: return None
    P=np.array(pts,float)
    keep=np.ones(len(P),bool)
    for _ in range(6):
        m,cc=np.polyfit(P[keep,0],P[keep,1],1)
        res=np.abs(P[:,1]-(m*P[:,0]+cc)); keep=res<2.0
        if keep.sum()<3: break
    return m,cc,int(keep.sum())

def vticks(a, line, yr, side, k0=4, k1=12):
    """tick stubs (horizontal) along a near-vertical axis x=m*y+c. side=+1 right of the axis, -1 left. returns rows."""
    m,c=line[0],line[1]
    rows=[]
    for y in range(yr[0],yr[1]):
        x=int(round(m*y+c))
        ks=[x+side*k for k in range(k0,k1)]
        if all(a[y,k] for k in ks): rows.append(y)
    g=[]
    for y in rows:
        if g and y-g[-1][-1]<=2: g[-1].append(y)
        else: g.append([y])
    return [float(np.mean(v)) for v in g]

def hticks(a, line, xr, side, k0=4, k1=12):
    """tick stubs (vertical) along a near-horizontal axis y=m*x+c. side=+1 below the axis, -1 above."""
    m,c=line[0],line[1]
    cols=[]
    for x in range(xr[0],xr[1]):
        y=int(round(m*x+c))
        ks=[y+side*k for k in range(k0,k1)]
        if all(a[k,x] for k in ks): cols.append(x)
    g=[]
    for x in cols:
        if g and x-g[-1][-1]<=2: g[-1].append(x)
        else: g.append([x])
    return [float(np.mean(v)) for v in g]


class Panel:
    """Bilinear data<->pixel map for one plot panel in GLOBAL page-pixel coordinates.
    L=(m,c,..) fit x=m*row+c of the left frame/axis line (data x = x0); R same for the right line (data x = x1).
    yL = [(row, value),...] ticks on the left line; yR optional ticks on the right line in the SAME y units; else tilt from slope bt (bottom line dy/dx)."""
    def __init__(self, a, L, R, x0, x1, yL, yR=None, bt=0.0, name=""):
        self.a=a; self.L=L; self.R=R; self.x0=x0; self.x1=x1; self.name=name
        yl=np.array(yL,float); self.pL=np.polyfit(yl[:,1], yl[:,0], 1)   # row = p0*y+p1
        self.resL=float(np.max(np.abs(np.polyval(self.pL, yl[:,1])-yl[:,0])))
        if yR is not None:
            yr=np.array(yR,float); self.pR=np.polyfit(yr[:,1], yr[:,0], 1)
            self.resR=float(np.max(np.abs(np.polyval(self.pR, yr[:,1])-yr[:,0])))
        else:
            wid=(R[0]*0+R[1])-(L[1]); self.pR=np.array([self.pL[0], self.pL[1]+bt*wid]); self.resR=0.0
    def xL(self,r): return self.L[0]*r+self.L[1]
    def xR(self,r): return self.R[0]*r+self.R[1]
    def to_pix(self,x,y):
        u=(x-self.x0)/(self.x1-self.x0)
        rl=np.polyval(self.pL,y); rr=np.polyval(self.pR,y); r=rl+u*(rr-rl)
        c=self.xL(r)+u*(self.xR(r)-self.xL(r))
        return c,r
    def to_data(self,c,r):
        r=np.asarray(r,float); c=np.asarray(c,float)
        u=(c-self.xL(r))/(self.xR(r)-self.xL(r))
        for _ in range(3):
            # solve r = rl(y)+u*(rr(y)-rl(y)) for y (linear)
            A=(1-u)*self.pL[0]+u*self.pR[0]; B=(1-u)*self.pL[1]+u*self.pR[1]
            y=(r-B)/A
            u=(c-self.xL(r))/(self.xR(r)-self.xL(r))
        return self.x0+u*(self.x1-self.x0), y
    def column_clusters(self, x, ylim, halfw=1, gap=2):
        """ink clusters along the data-vertical line at abscissa x. returns list of (y_center, thickness_in_rows)"""
        r0=int(np.polyval(self.pL,ylim[1]))-3; r1=int(np.polyval(self.pL,ylim[0]))+3
        r0,r1=min(r0,r1),max(r0,r1)
        u=(x-self.x0)/(self.x1-self.x0)
        rows=[]
        for r in range(r0-20,r1+20):
            c=self.xL(r)+u*(self.xR(r)-self.xL(r))
            ci=int(round(c))
            if self.a[r,ci-halfw:ci+halfw+1].any(): rows.append(r)
        g=[]
        for r in rows:
            if g and r-g[-1][-1]<=gap: g[-1].append(r)
            else: g.append([r])
        out=[]
        for v in g:
            rc=float(np.mean(v)); c=self.xL(rc)+u*(self.xR(rc)-self.xL(rc))
            _,yy=self.to_data(c,rc); out.append((float(yy),len(v)))
        return out
    def grid_overlay(self, path, xs, ys, region, scale=1.0, pts=None, secondary=None):
        """draw data gridlines (xs,ys lists) + optional points {label:(color,[(x,y),...])}; region=(c0,r0,c1,r1) global px"""
        c0,r0,c1,r1=region
        img=Image.fromarray(((~self.a[r0:r1,c0:c1]).astype(np.uint8))*255).convert("RGB")
        d=ImageDraw.Draw(img)
        for x in xs:
            p=[self.to_pix(x,y) for y in np.linspace(ys[0],ys[-1],40)]
            d.line([(c-c0,r-r0) for c,r in p], fill=(120,170,255), width=1)
            c,r=self.to_pix(x,ys[-1]); d.text((c-c0+2,r-r0+2),f"{x:g}",fill=(0,0,200))
        for y in ys:
            p=[self.to_pix(x,y) for x in np.linspace(xs[0],xs[-1],40)]
            d.line([(c-c0,r-r0) for c,r in p], fill=(255,170,120), width=1)
            c,r=self.to_pix(xs[0],y); d.text((c-c0+2,r-r0-10),f"{y:g}",fill=(200,60,0))
        if pts:
            for lab,(col,P) in pts.items():
                for (x,y) in P:
                    c,r=self.to_pix(x,y); d.ellipse((c-c0-5,r-r0-5,c-c0+5,r-r0+5),outline=col,width=2)
                if P:
                    c,r=self.to_pix(*P[0]); d.text((c-c0+7,r-r0-12),lab,fill=col)
        if scale!=1.0: img=img.resize((int(img.width*scale),int(img.height*scale)),Image.LANCZOS)
        img.save(path); return img.size


def zoom_view(P, path, xr, yr, xstep, ystep, scale=2.0, pts=None, label_every=1, pad=10):
    """Gridded zoom of the data window xr=(x0,x1), yr=(y0,y1). scale = s or (sx,sy)."""
    sx,sy = (scale,scale) if not isinstance(scale,(tuple,list)) else scale
    corners=[P.to_pix(x,y) for x in xr for y in yr]
    c0=int(min(c for c,r in corners))-pad; c1=int(max(c for c,r in corners))+pad
    r0=int(min(r for c,r in corners))-pad; r1=int(max(r for c,r in corners))+pad
    img=Image.fromarray(((~P.a[r0:r1,c0:c1]).astype(np.uint8))*255).convert("RGB")
    img=img.resize((int(img.width*sx),int(img.height*sy)),Image.LANCZOS)
    d=ImageDraw.Draw(img)
    xs=np.arange(xr[0],xr[1]+1e-9,xstep); ys=np.arange(yr[0],yr[1]+1e-9,ystep)
    for i,x in enumerate(xs):
        p=[P.to_pix(x,y) for y in (yr[0],yr[1])]
        d.line([((c-c0)*sx,(r-r0)*sy) for c,r in p],fill=(90,140,255),width=1)
        if i%label_every==0: d.text(((p[1][0]-c0)*sx+2,(p[1][1]-r0)*sy+2),f"{x:g}",fill=(0,0,220))
    for j,y in enumerate(ys):
        p=[P.to_pix(x,y) for x in (xr[0],xr[1])]
        d.line([((c-c0)*sx,(r-r0)*sy) for c,r in p],fill=(255,150,90),width=1)
        d.text(((p[0][0]-c0)*sx+2,(p[0][1]-r0)*sy-11),f"{y:g}",fill=(210,70,0))
    if pts:
        for lab,(col,PP) in pts.items():
            for (x,y) in PP:
                c,r=P.to_pix(x,y); d.ellipse(((c-c0)*sx-6,(r-r0)*sy-6,(c-c0)*sx+6,(r-r0)*sy+6),outline=col,width=2)
    img.save(path); return img.size

def fill_small_holes(a, region, maxarea=140):
    from scipy import ndimage as ndi
    c0,r0,c1,r1=region
    sub=a[r0:r1,c0:c1]
    holes=ndi.binary_fill_holes(sub)&~sub
    lab,n=ndi.label(holes)
    if n:
        areas=ndi.sum(holes,lab,range(1,n+1))
        keep=np.isin(lab,[i+1 for i,ar in enumerate(areas) if ar<=maxarea])
    else: keep=holes
    F=a.copy(); F[r0:r1,c0:c1]=sub|keep
    return F


def clusters2(P, x, ylim, halfw=2, gap=2):
    """like Panel.column_clusters but returns (yc, thickness_rows, y_a, y_b) with y_a<y_b the data extents"""
    r0=int(np.polyval(P.pL,ylim[1]))-3; r1=int(np.polyval(P.pL,ylim[0]))+3
    r0,r1=min(r0,r1),max(r0,r1)
    u=(x-P.x0)/(P.x1-P.x0)
    rows=[]
    for r in range(r0,r1):
        c=int(round(P.xL(r)+u*(P.xR(r)-P.xL(r))))
        if P.a[r,c-halfw:c+halfw+1].any(): rows.append(r)
    g=[]
    for r in rows:
        if g and r-g[-1][-1]<=gap: g[-1].append(r)
        else: g.append([r])
    out=[]
    for v in g:
        rc=float(np.mean(v)); ca=lambda rr: P.xL(rr)+u*(P.xR(rr)-P.xL(rr))
        yc=P.to_data(ca(rc),rc)[1]; ya=P.to_data(ca(v[0]),v[0])[1]; yb=P.to_data(ca(v[-1]),v[-1])[1]
        out.append((float(yc),len(v),float(min(ya,yb)),float(max(ya,yb))))
    return out

def track(P, seed, xs, ylim, base_tol=0.09, slope_tol=1.0, maxthick=17, halfw=2, damp=0.8, max_hold=3):
    """Follow one curve from seed=(x0,y0). Returns dict x->(y,status,half_unc).
    ok: single-marker cluster accepted; held: curve lies inside a thicker (merged) cluster, prediction held (excluded from history);
    stop: tracking abandoned after max_hold consecutive non-ok points (remaining abscissae omitted)."""
    xs=sorted(set(float(x) for x in xs)); x0,y0=seed
    res={x0:(y0,'seed',0.0)}
    for direction in (+1,-1):
        seq=[x for x in xs if (x>x0 if direction>0 else x<x0)]
        if direction<0: seq=seq[::-1]
        hist=[(x0,y0)]; nhold=0
        for x in seq:
            if len(hist)>=2:
                (xa,ya),(xb,yb)=hist[-2],hist[-1]
                yp=yb+damp*(yb-ya)/(xb-xa)*(x-xb) if xb!=xa else yb
            else: yp=hist[-1][1]
            tol=base_tol+slope_tol*abs(yp-hist[-1][1])*0.5
            cl=clusters2(P,x,ylim,halfw=halfw)
            best=None
            for (yc,t,ya_,yb_) in cl:
                if t<=maxthick:
                    d=abs(yc-yp)
                    if d<=tol and (best is None or d<best[0]): best=(d,yc,'ok',0.0)
            if best is None:
                for (yc,t,ya_,yb_) in cl:
                    if t>maxthick and ya_-0.03<=yp<=yb_+0.03:
                        best=(0,yp,'held',(yb_-ya_)/2); break
            if best is None:
                nhold+=1
                if nhold>max_hold: break
                res[x]=(yp,'held',tol); continue
            _,y,st,half=best
            if st=='ok': res[x]=(y,'ok',0.03); hist.append((x,y)); nhold=0
            else:
                nhold+=1
                if nhold>max_hold: break
                res[x]=(y,'held',half)
    return dict(sorted(res.items()))


def marker_blobs(P, xr, yr, rad=3, minarea=40, maxarea=260):
    """centroids (data units) + area of marker-sized blobs of P.a (use the hole-filled image) inside the data window"""
    from scipy import ndimage as ndi
    corners=[P.to_pix(x,y) for x in xr for y in yr]
    c0=int(min(c for c,r in corners))-10; c1=int(max(c for c,r in corners))+10
    r0=int(min(r for c,r in corners))-10; r1=int(max(r for c,r in corners))+10
    sub=P.a[r0:r1,c0:c1]
    yy,xx=np.ogrid[-rad:rad+1,-rad:rad+1]; disk=(xx**2+yy**2)<=rad*rad
    op=ndi.binary_opening(sub,structure=disk)
    lab,n=ndi.label(op)
    out=[]
    for i in range(1,n+1):
        m=lab==i; ar=int(m.sum())
        cy,cx=ndi.center_of_mass(m)
        x,y=P.to_data(cx+c0,cy+r0)
        if xr[0]<=x<=xr[1] and yr[0]<=y<=yr[1]: out.append((float(x),float(y),ar))
    return sorted(out)


def track_family(P, seeds, xs, ylim, base_tol=0.10, maxthick=17, halfw=2, damp=0.85, max_hold=4, ylim_mask=None):
    """Track several curves at once from seeds {name:(x0,y0)} (all same x0). Greedy global matching curve<->cluster per abscissa,
    so a cluster can be claimed by one curve only.  Returns {name:{x:(y,status,unc)}}; status seed|ok|held."""
    names=list(seeds); x0=seeds[names[0]][0]
    slope0={n:(seeds[n][2] if len(seeds[n])>2 else None) for n in names}
    seeds={n:(seeds[n][0],seeds[n][1]) for n in names}
    xs=sorted(set(float(x) for x in xs))
    res={n:{x0:(seeds[n][1],'seed',0.0)} for n in names}
    for direction in (+1,-1):
        seq=[x for x in xs if (x>x0 if direction>0 else x<x0)]
        if direction<0: seq=seq[::-1]
        hist={n:([(x0-direction*0.5,seeds[n][1]-direction*0.5*slope0[n]),(x0,seeds[n][1])] if slope0[n] is not None else [(x0,seeds[n][1])]) for n in names}
        if direction<0:
            for n in names:
                if slope0[n] is not None: hist[n]=[(x0+0.5,seeds[n][1]+0.5*slope0[n]),(x0,seeds[n][1])]
        nhold={n:0 for n in names}; alive={n:True for n in names}
        for x in seq:
            pred={}
            for n in names:
                if not alive[n]: continue
                h=hist[n]
                if len(h)>=2 and h[-1][0]!=h[-2][0]:
                    (xa,ya),(xb,yb)=h[-2],h[-1]; pred[n]=yb+damp*(yb-ya)/(xb-xa)*(x-xb)
                else: pred[n]=h[-1][1]
            cl=clusters2(P,x,ylim,halfw=halfw)
            thin=[c for c in cl if c[1]<=maxthick]; thick=[c for c in cl if c[1]>maxthick]
            pairs=[]
            for n,yp in pred.items():
                h=hist[n]; tol=base_tol+0.5*abs(yp-h[-1][1])
                for j,c in enumerate(thin):
                    d=abs(c[0]-yp)
                    if d<=tol: pairs.append((d/tol,n,j))
            pairs.sort(); usedn=set(); usedj=set(); got={}
            for d,n,j in pairs:
                if n in usedn or j in usedj: continue
                usedn.add(n); usedj.add(j); got[n]=thin[j][0]
            for n,yp in pred.items():
                if n in got:
                    res[n][x]=(got[n],'ok',0.03); hist[n].append((x,got[n])); nhold[n]=0
                else:
                    inside=[c for c in thick if c[2]-0.03<=yp<=c[3]+0.03]
                    nhold[n]+=1
                    if nhold[n]>max_hold or not inside:
                        if nhold[n]>max_hold: alive[n]=False
                        continue
                    c=inside[0]; res[n][x]=(yp,'held',(c[3]-c[2])/2)
    return {n:dict(sorted(r.items())) for n,r in res.items()}


def blobs_px(P, xr, yr, rad=3, pad=10):
    """marker-sized blob centroids in GLOBAL pixel coords with areas, for the data window"""
    from scipy import ndimage as ndi
    corners=[P.to_pix(x,y) for x in xr for y in yr]
    c0=int(min(c for c,r in corners))-pad; c1=int(max(c for c,r in corners))+pad
    r0=int(min(r for c,r in corners))-pad; r1=int(max(r for c,r in corners))+pad
    sub=P.a[r0:r1,c0:c1]
    yy,xx=np.ogrid[-rad:rad+1,-rad:rad+1]; disk=(xx**2+yy**2)<=rad*rad
    op=ndi.binary_opening(sub,structure=disk)
    lab,n=ndi.label(op)
    out=[]
    for i in range(1,n+1):
        m=lab==i; ar=int(m.sum()); cy,cx=ndi.center_of_mass(m)
        ys,xs_=np.where(m)
        out.append(dict(c=cx+c0,r=cy+r0,area=ar,w=int(xs_.max()-xs_.min()+1),h=int(ys.max()-ys.min()+1)))
    return out

def chain_family(P, blobs, seeds, xr, yr, single=(40,150), maxdev=15.0, maxstep=70.0, area_pen=0.25, order=None):
    """seeds {name:(x,y)} data coords -> each snapped to nearest single blob. Chain both directions in pixel space.
    returns {name:[(x,y,area,flag)]} sorted by x. A blob is claimed by one curve only (processed in `order`)."""
    sing=[b for b in blobs if single[0]<=b['area']<=single[1]]
    used=set(); out={}
    names=order or list(seeds)
    for n in names:
        sx,sy=seeds[n]; sc,sr_=P.to_pix(sx,sy)
        cand=[(np.hypot(b['c']-sc,b['r']-sr_),i) for i,b in enumerate(sing) if i not in used]
        d,i0=min(cand)
        if d>20: print("seed",n,"no blob within 20px (nearest",round(d,1),")"); 
        chain=[i0]; used.add(i0); seedarea=sing[i0]['area']
        for direction in (+1,-1):
            cur=i0; prev=None
            hist=[i0]
            while True:
                bc=sing[cur]
                if prev is not None:
                    bp=sing[prev]; vec=np.array([bc['c']-bp['c'],bc['r']-bp['r']])
                else:
                    vec=np.array([direction*25.0,0.0])
                if np.hypot(*vec)<1: vec=np.array([direction*25.0,0.0])
                # predicted next position: same step length along direction of last step
                pred=np.array([bc['c'],bc['r']])+vec
                best=None
                for j,b in enumerate(sing):
                    if j in used: continue
                    dx=b['c']-bc['c']
                    if direction*dx<=2: continue
                    dist=np.hypot(b['c']-bc['c'],b['r']-bc['r'])
                    if dist>maxstep: continue
                    dev=np.hypot(b['c']-pred[0],b['r']-pred[1])
                    # allow for different step lengths: use perpendicular deviation from the extended line
                    u=vec/np.hypot(*vec); w=np.array([b['c']-bc['c'],b['r']-bc['r']])
                    perp=abs(u[0]*w[1]-u[1]*w[0]); along=u@w
                    if along<4: continue
                    cost=perp+area_pen*abs(b['area']-seedarea)/10.0
                    if perp<=maxdev and (best is None or cost<best[0]): best=(cost,j)
                if best is None:
                    # fallback: wider corridor, nearest blob ahead (bridges a stick-slip dip / one missing marker)
                    for j,b in enumerate(sing):
                        if j in used: continue
                        dx=b['c']-bc['c']
                        if direction*dx<=2: continue
                        w=np.array([b['c']-bc['c'],b['r']-bc['r']]); dist=np.hypot(*w)
                        if dist>maxstep: continue
                        u=vec/np.hypot(*vec); perp=abs(u[0]*w[1]-u[1]*w[0]); along=u@w
                        if along<4 or perp>2.4*maxdev: continue
                        cost=dist+area_pen*abs(b['area']-seedarea)/10.0
                        if best is None or cost<best[0]: best=(cost,j)
                if best is None: break
                j=best[1]; used.add(j); chain.append(j); prev,cur=cur,j
        pts=[]
        for j in chain:
            b=sing[j]; x,y=P.to_data(b['c'],b['r']); pts.append((float(x),float(y),b['area'],''))
        out[n]=sorted(pts)
    return out


def row_clusters(P, y, xlim, halfw=1, gap=2):
    """ink clusters along the data-horizontal line at ordinate y, over x in xlim. returns [(x_center, thickness_px)]"""
    us=np.linspace((xlim[0]-P.x0)/(P.x1-P.x0),(xlim[1]-P.x0)/(P.x1-P.x0),int(abs(P.to_pix(xlim[1],y)[0]-P.to_pix(xlim[0],y)[0]))+1)
    pts=[]
    for u in us:
        x=P.x0+u*(P.x1-P.x0); c,r=P.to_pix(x,y)
        ci=int(round(c)); ri=int(round(r))
        pts.append((x,P.a[ri-halfw:ri+halfw+1,ci].any()))
    out=[]; cur=[]
    for x,v in pts+[(None,False)]:
        if v: cur.append(x)
        else:
            if cur:
                out.append((float(np.mean(cur)),len(cur)))
                cur=[]
    # merge gaps<=gap
    merged=[]
    for c in out:
        if merged and abs(c[0]-merged[-1][0])<gap*0+0.0: pass
        merged.append(c)
    return merged

def row_track(P, y_levels, x_guess, xwin, tol, maxthick_px=60):
    """follow a steep curve upward through y_levels (ascending), predicting x by extrapolation. returns [(x,y)]"""
    res=[]; hist=[]
    xp=x_guess
    for y in y_levels:
        if len(hist)>=2 and hist[-1][1]!=hist[-2][1]:
            (xa,ya),(xb,yb)=hist[-2],hist[-1]; xp=xb+(xb-xa)/(yb-ya)*(y-yb)
        elif hist: xp=hist[-1][0]
        cl=row_clusters(P,y,(max(P.x0,xp-xwin),min(P.x1,xp+xwin)))
        cl=[c for c in cl if c[1]<=maxthick_px]
        if not cl: continue
        c=min(cl,key=lambda c:abs(c[0]-xp))
        if abs(c[0]-xp)>tol: continue
        res.append((c[0],y)); hist.append((c[0],y))
    return res
