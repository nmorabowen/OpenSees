"""Lam & Tatsuoka (1988) phi vs e0.3 point sets (air-pluviated Toyoura, sigma3' = 98 kPa).
 Fig. 4(a) p. 94 (PDF page 6): TC-3, H/W = 1.0, omega = 0/30/60/90 deg.   axis calibration (crop px of the 1:1 native crop): e 0.6 at 118, 0.7 at 285, 0.8 at 453, 0.9 at 620; phi 45 deg at 20, 40 at 227, 35 at 435, 30 at 640.
 Fig. 8(a) p. 96 (PDF page 8): PSC-3, omega = 0, xi = 90, H/W = 1.9 / 1.0 / 0.5 / 0.25.  e 0.6 at the left frame, 0.7 at 168 px, 0.1 per 168 px; phi 50 deg at the frame top, 30 deg at the frame bottom (40.6 px/deg).
Marker centres were read from the native-resolution crops (+-3 px => e +-0.002, phi +-0.08 deg; plus +-0.2 deg for marker overlap).
"""
import csv, os
HERE=os.path.dirname(os.path.abspath(__file__))
def tc(px,py): return (0.6+(px-118)/1675.0, 45-(py-20)/41.4)
TC={0:[(234,77),(243,140),(320,165),(485,322),(487,368)],
    30:[(172,134),(413,337),(513,424),(507,440)],
    60:[(193,212),(197,224),(492,447)],
    90:[(240,255),(257,270),(283,311),(521,457),(536,462)]}
with open(os.path.join(HERE,"fig4a_tc_phi_vs_e.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Lam & Tatsuoka (1988) Fig. 4(a): phi in TC-3, H/W = 1.0, sigma3' = 98 kPa vs e0.3.  kind TC.  tolerance e +-0.005, phi +-0.3 deg (reading).\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["sigma3_kPa","e","phi_peak_deg","eps_peak_pct","kind","omega_deg","H_over_W"])
    for om,pts in TC.items():
        for px,py in pts:
            e,ph=tc(px,py); w.writerow([98.0,round(e,4),round(ph,2),"","TC",om,1.0])
def ps(px,py): return (0.6+(px-123.5)/1680.0, 50-(py-42.5)/40.6)       # crop px -> native: +125, +433 ; native frame: left 248.5, top 475.5
PS={1.9:[(190,100),(219,158),(388,378),(458,432),(465,481),(494,495)],
    1.0:[(210,147),(490,473)],
    0.5:[(200,88),(484,359)],
    0.25:[(212,61),(493,278)]}
with open(os.path.join(HERE,"fig8a_psc_phi_vs_e.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Lam & Tatsuoka (1988) Fig. 8(a): phi in PSC-3, omega = 0, sigma3' = 98 kPa vs e0.3 for four H/W.  kind PS.  H/W = 1.9 is the one the paper uses for PSC.  tolerance e +-0.005, phi +-0.3 deg.\n")
    f.write("# b at failure is plotted on the same figure (0.20-0.34, nearly independent of e) but is not transcribed here; use fig9 curves for b.\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["sigma3_kPa","e","phi_peak_deg","eps_peak_pct","kind","omega_deg","H_over_W"])
    for hw,pts in PS.items():
        for px,py in pts:
            e,ph=ps(px,py); w.writerow([98.0,round(e,4),round(ph,2),"","PS",0,hw])
print(open(os.path.join(HERE,"fig8a_psc_phi_vs_e.csv")).read())
