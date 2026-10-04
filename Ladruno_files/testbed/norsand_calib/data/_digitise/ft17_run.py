import math, json
from famrun import *
a=page("ft",16)
X0,Y0=297,436
L=fit_vline(a,(X0+150,X0+185),(Y0+40,Y0+440),40); R=fit_vline(a,(X0+1020,X0+1060),(Y0+40,Y0+440),40)
P=Panel(a,L,R,-2.0,math.log10(5.0),[(491.4,45),(871.4,30)],[(484.6,45),(864.6,30)],name="ft17")
# check tilt: right frame top ~ T fit at x=1330
xs=[-1.55+0.05*i for i in range(0,50)]
out={}
for nm,y0 in (("e0.70",42.5),("e0.85",36.2)):
    pts=follow_isolated(P,xs,y0,1.2,(33,45),maxthick=14,halfw=1)
    out[nm]=[(10**x,y) for x,y in pts]
    print(nm,len(pts),[(round(10**x,3),round(y,2)) for x,y in pts][::4])
json.dump(out,open("ft17_out.json","w"))
pts={k:((255,0,0) if k=="e0.70" else (0,0,255),[(math.log10(x),y) for x,y in v]) for k,v in out.items()}
zoom_view(P,"ov_ft17.png",(-2.0,0.7),(30,45),0.3,5,scale=(0.5,1.0),pts=pts,pad=5)
print("done")
