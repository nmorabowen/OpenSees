from dtool import *
def make_panel(tag, pg, xL, xR, y0, y1, x0, x1, rowsL, rowsR, vals, fill=True, name=""):
    """rowsL/rowsR: tick rows (global px) on left/right line for values `vals` (same y units). fits axis lines near xL, xR."""
    a=page(tag,pg)
    L=fit_vline(a,(int(xL)-14,int(xL)+14),(y0,y1),40); R=fit_vline(a,(int(xR)-16,int(xR)+16),(y0,y1),40)
    if fill: a=fill_small_holes(a,(int(xL)+8,y0+8,int(xR)-8,y1-8))
    P=Panel(a,L,R,x0,x1,list(zip(rowsL,vals)),list(zip(rowsR,vals)),name=name)
    return P
def list_clusters(P, xs, ylim, maxt=45, ymin=None):
    for x in xs:
        cl=clusters2(P,x,ylim,halfw=2)
        print(x,[(round(c[0],2),c[1]) for c in cl if c[1]<maxt and (ymin is None or c[0]>=ymin)])
