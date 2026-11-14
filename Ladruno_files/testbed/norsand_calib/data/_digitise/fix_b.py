import json
from famrun import *
P=make_panel("ft",7,356,1133,1600,2650,0,15,[1613.5,2427.5,2631.0],[1614.5,2428.5,2631.5],[8,0,-2],name="ft5b")
d=json.load(open("ft5b_out.json"))
xs=[1.0+0.25*i for i in range(0,17)]
left=family(P,{'0.02':(4,6.34,0.25),'0.05':(4,5.87,0.2),'0.1':(4,5.71,0.2)},{},xs,(4.5,7.0),{k:(1.0,5.0) for k in ('0.02','0.05','0.1')},base_tol=0.12,max_hold=3,maxthick=25)
for k,v in left.items():
    cur=d['sr'][k]; xmin=min(p[0] for p in cur) if cur else 99
    add=[p for p in v if p[0]<xmin-0.1]
    print(k,"adding",len(add),[(round(x,2),round(y,2)) for x,y,s in add])
    d['sr'][k]=sorted([list(p) for p in add]+cur)
json.dump(d,open("ft5b_out.json","w"))
pts={k:(PALETTE[i],[(x,y) for x,y,s in v]) for i,(k,v) in enumerate(d['sr'].items())}
zoom_view(P,"ov_ft5b_sr.png",(0,15),(2.5,7.0),1,0.5,scale=(1.3,1.5),pts=pts)
print("done")
