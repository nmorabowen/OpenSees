import sys; sys.path.insert(0,'.')
from dtool import *
def panel_ft5a():
    tag,n="ft",7
    a=page(tag,n)
    nb=prev2nat(tag,(70,125,395,500)); X0,Y0=nb[0],nb[1]
    L=fit_vline(a,(X0+120,X0+170),(Y0+60,Y0+1100),40); R=fit_vline(a,(X0+880,X0+960),(Y0+60,Y0+1100),40)
    yL=[(452.5,8),(1270.5,0),(1469.5,-2)]
    yR=[(75.5+Y0,8),(887.5+Y0,0),(1090+Y0,-2)]
    return Panel(a,L,R,0,15,yL,yR,name="ft5a"), (X0,Y0)
def filled_panel_ft5a():
    P,(X0,Y0)=panel_ft5a()
    F=fill_small_holes(P.a,(X0+100,Y0+60,X0+960,Y0+1120))
    return Panel(F,P.L,P.R,0,15,[(452.5,8),(1270.5,0),(1469.5,-2)],[(75.5+Y0,8),(887.5+Y0,0),(1090+Y0,-2)],name="ft5a_f"),P
