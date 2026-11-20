import json, pickle
from ft5a_setup import *
from curate import despike
PF,P=filled_panel_ft5a()
cols={'0.1':(255,0,0),'0.2':(0,150,0),'0.5':(0,0,255),'1.0':(200,0,200),'2.0':(255,140,0),'4.0':(0,170,170)}
# ---- sr: marker chains (already computed) curated by x-range + despike
R=pickle.load(open("ft5a_chain.pkl","rb"))
rng={'0.1':(0.3,11.5),'0.2':(1.2,10.0),'0.5':(0.6,11.0),'1.0':(0.8,12.5),'2.0':(0.7,12.5),'4.0':(0.7,12.5)}
sr={}
for k in cols:
    pts=[(x,y) for x,y,a,f in R[k] if rng[k][0]<=x<=rng[k][1]]
    sr[k]=despike(pts)
# ---- eps_v: strip family tracker from seeds at 7 % with slopes
xs=[0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5]
seeds={'0.1':(7,2.55,0.40),'0.2':(7,2.41,0.40),'0.5':(7,2.27,0.40),'1.0':(7,1.95,0.38),'2.0':(7,1.53,0.32),'4.0':(7,1.17,0.26)}
E=track_family(PF,seeds,xs,(-0.3,4.4),base_tol=0.09,max_hold=6)
ev={k:[(x,-2*y) for x,(y,st,u) in r.items() if st in('ok','seed')] for k,r in E.items()}
json.dump({'sr':sr,'ev':ev},open("ft5a_out.json","w"))
pts_sr={k:(cols[k],v) for k,v in sr.items()}
pts_ev={k:(cols[k],[(x,-v/2) for x,v in vv]) for k,vv in ev.items()}
zoom_view(P,"ov2_ft5a_sr.png",(0,15),(3.0,6.4),1,0.5,scale=(1.35,1.6),pts=pts_sr)
zoom_view(P,"ov2_ft5a_ev.png",(0,15),(-0.3,4.6),1,0.5,scale=(1.35,1.3),pts=pts_ev)
for k in cols: print(k,"sr",len(sr[k]),"ev",len(ev[k]),[round(x,1) for x,_ in ev[k]][:3],'...',[round(x,1) for x,_ in ev[k]][-2:])
print("done")
