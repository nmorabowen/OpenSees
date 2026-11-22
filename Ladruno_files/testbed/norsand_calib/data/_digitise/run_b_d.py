import json, sys
from famrun import *
which=sys.argv[1]
if which=='d':
    P=make_panel("ft",7,1401,2419,1790,2780,0,20,[1825.5,2434.5,2636.0],[1835.5,2446.5,2649.5],[6,0,-2],name="ft5d")
    xs=[x*0.25 for x in range(4,64)]
    ylow=[2.0+0.1*i for i in range(0,20)]
    rowlow=row_track(P,ylow,1.0,0.8,0.5)
    print("rowlow",rowlow[:3],len(rowlow))
    tail=follow_isolated(P,[9+0.25*i for i in range(0,26)],4.12,0.3,(3.2,4.6))
    print("tail",len(tail))
    sr=family(P,{'0.02':(6,4.18,0.0),'low':(6,3.95,0.05)},{'0.02':(6,4.18)},xs,(2.0,4.6),{'0.02':(0.2,10.5),'low':(0.2,15.3)},order=['0.02'],maxthick=60,extra={'low':rowlow+tail})
    ev=family(P,{'0.02':(6,0.56,0.08),'low':(6,0.16,0.04)},{'0.02':(6,0.56)},xs,(-0.3,1.4),{'0.02':(0.5,10.5),'low':(0.5,15.2)},order=['0.02'],maxthick=60)
else:
    P=make_panel("ft",7,356,1133,1600,2650,0,15,[1613.5,2427.5,2631.0],[1614.5,2428.5,2631.5],[8,0,-2],name="ft5b")
    xs=[x*0.25 for x in range(4,60)]
    sr=family(P,{'0.02':(8,6.12,-0.12),'0.05':(8,5.52,-0.25),'0.1':(8,5.64,-0.08)},
              {'0.02':(5,6.38),'0.05':(5,5.91),'0.1':(5,5.75),'0.1#2':(11,5.2)},xs,(2.5,7.0),{'0.02':(0.3,9.6),'0.05':(0.3,9.3),'0.1':(0.3,13.6)},order=['0.1#2','0.1','0.05','0.02'],chain_area=(40,300))
    ev=family(P,{'pair':(8,2.75,0.33),'0.1':(8,2.45,0.27)},{'0.1':(8,2.45)},xs,(-0.3,3.8),{'pair':(0.3,10.0),'0.1':(0.3,13.6)},order=['0.1'],maxthick=40)
json.dump({'sr':sr,'ev':ev},open(f"ft5{which}_out.json","w"))
pts_sr={k:(PALETTE[i],[(x,y) for x,y,s in v]) for i,(k,v) in enumerate(sr.items())}
pts_ev={k:(PALETTE[i],[(x,y) for x,y,s in v]) for i,(k,v) in enumerate(ev.items())}
if which=='d':
    zoom_view(P,f"ov_ft5{which}_sr.png",(0,16),(2.0,4.6),2,0.5,scale=(1.0,1.6),pts=pts_sr)
    zoom_view(P,f"ov_ft5{which}_ev.png",(0,16),(-0.3,1.2),2,0.25,scale=(1.0,3.0),pts=pts_ev)
else:
    zoom_view(P,f"ov_ft5{which}_sr.png",(0,15),(2.5,7.0),1,0.5,scale=(1.3,1.5),pts=pts_sr)
    zoom_view(P,f"ov_ft5{which}_ev.png",(0,15),(-0.3,3.8),1,0.25,scale=(1.3,1.8),pts=pts_ev)
for k,v in list(sr.items())+list(ev.items()): print(k,len(v),round(v[0][0],2),round(v[-1][0],2))
print("done")
