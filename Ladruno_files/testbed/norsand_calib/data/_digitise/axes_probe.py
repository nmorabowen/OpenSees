import sys, time
from dtool import *
def probe(tag, n, pbox, name):
    a=page(tag,n); nb=prev2nat(tag,pbox); X0,Y0,X1,Y1=nb
    sub=a[Y0:Y1,X0:X1]
    # vertical lines: long runs in columns
    vl=lines(sub,0,300); hl=lines(sub,1,300)
    print(name,"box",nb)
    print(" vertical:",[(round(c+X0,1),s+Y0,e+Y0) for c,s,e,k in vl])
    print(" horizontal:",[(round(r+Y0,1),s+X0,e+X0) for r,s,e,k in hl])
if __name__=="__main__":
    t=time.time()
    probe("ft",7,(70,505,395,880),"b")
    probe("ft",7,(400,125,800,505),"c")
    probe("ft",7,(400,570,800,880),"d")
    print(time.time()-t)
