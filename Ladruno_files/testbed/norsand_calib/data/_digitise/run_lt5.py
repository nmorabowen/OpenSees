import json, sys
from famrun import *
from lt5_setup import *
Pa,Pb=lt5_panels()
xs=[0.5*i for i in range(2,60)]
res={}
# (a) H/W = 1.0
sra=family(Pa,{'0':(10,5.44,-0.03),'30':(10,5.17,-0.08),'60':(10,4.58,-0.01),'90':(10,4.37,0.0)},{}, xs,(0.8,6.2),{k:(0.8,30) for k in ('0','30','60','90')},base_tol=0.09,max_hold=8,maxthick=25)
eva=family(Pa,{'0':(10,1.67,0.19),'30':(10,1.55,0.19),'60':(10,1.22,0.18),'90':(10,1.01,0.17)},{}, xs,(-0.2,5.4),{k:(1.0,30) for k in ('0','30','60','90')},base_tol=0.08,max_hold=8,maxthick=25)
# (b) H/W = 0.25
srb=family(Pb,{'0':(10,5.37,0.05),'30':(10,4.97,0.04),'6090':(10,4.27,0.04)},{}, xs,(0.8,6.2),{k:(0.8,30) for k in ('0','30','6090')},base_tol=0.09,max_hold=8,maxthick=25)
evb=family(Pb,{'0':(20,2.66,0.17),'30':(20,2.47,0.17),'60':(20,2.07,0.17),'90':(20,1.91,0.16)},{}, xs,(-0.2,5.4),{k:(1.0,30) for k in ('0','30','60','90')},base_tol=0.08,max_hold=8,maxthick=25)
json.dump({'a':{'sr':sra,'ev':eva},'b':{'sr':srb,'ev':evb}},open("lt5_out.json","w"))
def ov(P,sr,ev,name):
    pts={('s'+k):(PALETTE[i],[(x,y) for x,y,s in v]) for i,(k,v) in enumerate(sr.items())}
    pts.update({('e'+k):(PALETTE[i],[(x,y) for x,y,s in v]) for i,(k,v) in enumerate(ev.items())})
    zoom_view(P,name,(0,30),(-0.3,6.2),5,1,scale=(0.9,1.5),pts=pts)
ov(Pa,sra,eva,"ov_lt5a.png"); ov(Pb,srb,evb,"ov_lt5b.png")
for nm,d in (("a sr",sra),("a ev",eva),("b sr",srb),("b ev",evb)):
    for k,v in d.items(): print(nm,k,len(v),round(v[0][0],2) if v else None,round(v[-1][0],2) if v else None)
print("done")
