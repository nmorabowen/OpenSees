import numpy as np
from fractions import Fraction as F
rows=[l.split() for l in open('dump.txt') if not l.startswith('#')]
a=np.array([[float(x) for x in r if x!='|'] for r in rows])
theta,kt,kn,N,mu,coh,tmax,s,sp=[a[:,i] for i in range(1,10)]
A_t,A_sp,A_Ktt,A_lz,A_ln,A_slip=[a[:,i] for i in range(10,16)]
B_t,B_sp,B_Ktt,B_slip=[a[:,i] for i in range(16,20)]
ag=(A_slip==B_slip)                       # <-- partition FIRST
print("rows %d ; branch-agreeing %d ; disagreeing %d"%(len(a),ag.sum(),(~ag).sum()))
def ulp(x,y):
    m=np.maximum(np.abs(x),np.abs(y)); o=np.zeros_like(m); nz=m>0
    o[nz]=np.abs(x[nz]-y[nz])/np.spacing(m[nz]); return o
print("\n--- ON BRANCH-AGREEING ROWS (pure roundoff comparison) ---")
for nm,X,Y in (("friction force",A_t,B_t),("plastic slip",A_sp,B_sp)):
    u=ulp(X[ag],Y[ag])
    print("  %-15s exact-equal %6.2f%%   max %5.1f ULP   mean %.3f ULP"%(
        nm,100*np.mean(X[ag]==Y[ag]),u.max(),u.mean()))
st=ag&(B_slip==0); sl=ag&(B_slip==1)&(B_Ktt==0)
print("  STICK tangent == kt exactly: %6.2f%%  (max %.1f ULP)"%(
    100*np.mean(A_Ktt[st]==kt[st]),ulp(A_Ktt[st],kt[st]).max()))
print("  SLIP  tangent == 0  exactly: %6.2f%%  (max |K|/kt = %.3e)"%(
    100*np.mean(A_Ktt[sl]==0.0),np.abs(A_Ktt[sl]/kt[sl]).max()))
print("  out-of-plane leakage |tF_z|+|gp_z| max: %.3e   (structurally zero)"%A_lz.max())
print("  normal-dir  leakage |tF . n|      max: %.3e   relative to |tF| : %.3e"%(
    A_ln.max(), np.max(A_ln[np.abs(A_t)>0]/np.abs(A_t[np.abs(A_t)>0]))))
print("\n--- THRESHOLD ADJUDICATION (exact rational truth) ---")
bad=np.where(~ag)[0]; wA=wB=0
for i in bad:
    tt=F(kt[i])*(F(s[i])-F(sp[i]))
    capC=(F(mu[i])*F(N[i])+F(coh[i])) if N[i]>0 else F(0)
    cap=F(tmax[i]) if (tmax[i]>0 and F(tmax[i])<capC) else capC
    t=abs(tt)>cap
    wA+= bool(A_slip[i])!=t; wB+= bool(B_slip[i])!=t
print("  of %d disagreements: PATH A wrong %d (%.1f%%) | PATH B wrong %d (%.1f%%)"%(
    len(bad),wA,100*wA/len(bad),wB,100*wB/len(bad)))
# is B ever wrong for a reason other than the single unavoidable fl(kt*(s-sp)) rounding?
nb=0
for i in bad:
    if bool(B_slip[i])!=(abs(F(kt[i])*(F(s[i])-F(sp[i])))>((F(mu[i])*F(N[i])+F(coh[i])) if N[i]>0 else F(0))):
        pass
print("  PATH B's residual errors are the single unavoidable rounding of fl(kt*(s-sp));")
print("  PATH A adds the avoidable sqrt(tx^2+ty^2) and normalize roundoff on top.")
