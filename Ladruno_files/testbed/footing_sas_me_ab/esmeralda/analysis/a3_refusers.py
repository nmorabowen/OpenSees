import numpy as np, csv, os
CK='ck'
def gam(e): return np.sqrt((e[:,0]-e[:,1])**2+e[:,2]**2)
last={'E_A':'field_last_converged.npz','E_B':'field_last_converged.npz','E_D':'field_last_converged.npz','E_B16':'field_last_converged.npz'}
prev={'E_A':'field_step00185.npz','E_B':'field_step00370.npz','E_D':'field_step00385.npz','E_B16':'field_step00070.npz'}
flo=list(csv.DictReader(open('tables/floor_refusers.csv')))
out=[]
for leg in last:
    z=np.load(f'{CK}/{leg}/ckpt/{last[leg]}'); z0=np.load(f'{CK}/{leg}/ckpt/{prev[leg]}')
    zf=np.load(f'{CK}/{leg}/ckpt/field_flip.npz') if os.path.exists(f'{CK}/{leg}/ckpt/field_flip.npz') else None
    gi=gam(z['eps']-z0['eps'])
    gt=gam(z['eps']-zf['eps']) if zf is not None else None
    for r in [r for r in flo if r['leg']==leg]:
        k=int(r['k'])
        row=dict(leg=leg,element=r['element'],gp=r['gp'],x=r['x'],y=r['y'],codes=r['codes'],p=r['p_kPa'],eta=r['eta'],rho=r['rho_alpha'],
          gam_incr_pct=round(100*(gi<gi[k]).mean(),2),gam_incr_over_max=round(gi[k]/gi.max(),3))
        if gt is not None: row.update(gam_tot=round(gt[k],4),gam_tot_pct=round(100*(gt<gt[k]).mean(),2),gam_tot_over_max=round(gt[k]/gt.max(),3))
        # same-depth band peak
        m=np.abs(z['gy']-z['gy'][k])<1e-6
        side=np.sign(z['gx'][k]) if abs(z['gx'][k])>1e-9 else 1
        m2=m&(np.sign(z['gx'])==side)
        kk=np.where(m2)[0][np.argmax(gi[m2])]
        row['band_x_at_depth']=round(float(z['gx'][kk]),3)
        out.append(row); print(row)
with open('tables/floor_refusers_in_band.csv','w',newline='') as f:
    keys=sorted({k for r in out for k in r}, key=lambda s: list(out[0]).index(s) if s in out[0] else 99)
    w=csv.DictWriter(f,keys); w.writeheader(); w.writerows(out)
# cumulative refusals to s/B 0.0135: B8 vs B16 ; first loadingNonPosH
import re
for leg in ['E_B','E_B16','E_D','E_A']:
    S=list(csv.DictReader(open(f'../runs/{leg}/steps.csv')))
    c=sum(int(r['cap_step']) for r in S if float(r['s_over_B'])<=0.01352)
    nstep=sum(1 for r in S if float(r['s_over_B'])<=0.01352)
    fl=sum(int(r['fails_before']) for r in S if float(r['s_over_B'])<=0.01352)
    first=None
    for ln in open(f'../runs/{leg}/logs/log.log',errors='replace'):
        m=re.search(r"refusals step (\d+): .*loadingNonPosH",ln)
        if m: first=int(m.group(1)); break
    fs=[r['s_over_B'] for r in S if int(r['step'])==first] if first else None
    print(leg,'refusals/caps to s/B 0.01352:',c,'steps',nstep,'failed attempts',fl,'| first loadingNonPosH step',first,fs)
