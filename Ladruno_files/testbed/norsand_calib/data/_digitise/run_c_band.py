import json
from famrun import *
from lt9_setup import *
# ---- FT Fig 5(c): top two sr curves (loose, sigma_c' 0.1 and 0.2)
Pc=make_panel("ft",7,1396,2420,380,1500,0,20,[637.5,1250.3,1448.5],[646.5,1257.0,1460.5],[6,0,-2],name="ft5c")
xs=[0.5*i for i in range(2,40)]
src=family(Pc,{'0.1':(14,4.42,0.0),'0.2':(14,4.17,0.0)},{'0.1':(14,4.42),'0.2':(14,4.17)},xs,(2.0,4.8),{'0.1':(1.0,19.6),'0.2':(1.0,19.6)},base_tol=0.07,max_hold=4,maxthick=17,order=['0.2','0.1'],despike_thr=0.08)
pts={k:(PALETTE[i],[(x,y) for x,y,s in v]) for i,(k,v) in enumerate(src.items())}
zoom_view(Pc,"ov_ft5c_sr.png",(0,20),(2.5,4.8),2,0.25,scale=(1.0,3.0),pts=pts)
for k,v in src.items(): print("5c",k,len(v),round(v[0][0],2),round(v[-1][0],2))
json.dump({'sr':src},open("ft5c_out.json","w"))
# ---- LT Fig 9: eps_v band centreline D in [0.5,8.5]
P=lt9_panel()
rows=[]
for x in [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9]:
    cl=clusters2(P,x,(-0.9,1.6),halfw=2)
    rows.append((x,[(round(c[0],2),c[1],round(c[2],2),round(c[3],2)) for c in cl if c[1]<60]))
for r in rows: print("band",r[0],r[1])
print("done")
