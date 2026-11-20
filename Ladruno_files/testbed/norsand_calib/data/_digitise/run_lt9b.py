import json, pickle
from lt9_setup import *
P=lt9_panel()
bl=pickle.load(open("lt9_blobs.pkl","rb"))
def xy(i):
    b=bl[i]; x,y=P.to_data(b['c'],b['r']); return (float(x),float(y))
IDS={
 'sr_tri':[77,60,40,31,30,55,68,67,59,54,57,70,76,74,72],
 'sr_ci':[79,66,45,34,29,36,46,64,75,73,69,65,62,61,56],
 's2_tri':[129,100,103,101],
 's2_ci':[124,112,109,107],
 'ev_tri':[146,137,134,132,131,127,125],
 'ev_ci':[145,135,130,126,120,118,117,116],
}
out={k:[xy(i) for i in v] for k,v in IDS.items()}
out={k:sorted(v) for k,v in out.items()}
json.dump(out,open("lt9_out.json","w"))
cols=dict(zip(IDS,[(0,150,0),(255,0,0),(0,100,0),(200,0,0),(0,170,170),(255,140,0)]))
pts={k:(cols[k],v) for k,v in out.items()}
zoom_view(P,"ov_lt9b.png",(0,30),(-0.8,7.9),5,1,scale=(1.0,1.5),pts=pts)
for k,v in out.items(): print(k,[(round(x,2),round(y,2)) for x,y in v])
print("done")
