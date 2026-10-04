"""Fukushima & Tatsuoka (1984) summary figures, digitised from native-resolution crops of the scan.
 Fig. 17 (p. 45, PDF page 16): corrected phi vs sigma3' (log axis) at e0.3 = 0.70 and 0.85 -- the two curves traced automatically (single isolated line, 1 px = 0.04 deg).
 Fig. 19 (p. 46, PDF page 17): volumetric strain v at eps_a = 5 and 10 % vs sigma3' corrected by Method III (log), e0.3 = 0.70 / 0.85 -- markers read by eye.
 Fig. 20 (p. 46, PDF page 17): axial strain at sigma1'/sigma3' = 3 and 4 vs sigma3' (log), e0.3 = 0.70 / 0.85 -- markers read by eye.
Axis calibration of the crops (px): Fig. 19: v -8 % at y 24, 0 at 427 (50.1 px per %), sigma3' 0.1 at x 355, 1 at 665 (310 px/decade).
 Fig. 20: eps_a 0 at y 622, 4 at 120 (125.5 px per %), sigma3' 0.1 at x 410, 1 at 732 (322 px/decade).   1 kgf/cm2 = 98 kPa.
Reading tolerance: v +-0.1 %, eps_a +-0.05 %, sigma3' +-4 % (marker centre +-4 px).
"""
import csv, json, math, os
HERE=os.path.dirname(os.path.abspath(__file__))
KG=98.0
# ---- Fig. 17
c=json.load(open(os.path.join(HERE,"_ft17_curves.json")))
with open(os.path.join(HERE,"fig17_phi_corrected_vs_sigma3.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Fukushima & Tatsuoka (1984) Fig. 17: phi corrected for membrane forces (average of Methods I and III), mid-height of sample, vs sigma3' (log axis); drained TC; saturated air-pluviated Toyoura.\n")
    f.write("# Traced automatically from the curve (1 px = 0.04 deg).  The first point of the e0.70 trace (sigma3' 0.028 kgf/cm2, 43.4 deg) is a tracing artefact and dropped.  The paper's own extrapolations to sigma3'=0: 42.4 deg (e0.3=0.70) and 36.1 deg (e0.3=0.85); the traces start at 42.2 and 36.1.\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["e0p3","sigma3_kgf_cm2","sigma3_kPa","phi_deg"])
    for k,e in (("e0.70",0.70),("e0.85",0.85)):
        for i,(s,p) in enumerate(c[k]):
            if k=="e0.70" and i==0: continue
            w.writerow([e,round(s,4),round(s*KG,2),round(p,2)])
# ---- Fig. 19
def f19(px,py): return (10**((px-665)/310.0), (py-427)/50.1)
F19={("0.70","10%"):[(233,135),(295,139),(383,150),(467,159),(582,174),(670,195),(766,233),(864,290)],
     ("0.70","5%"):[(227,277),(292,282),(371,296),(454,309),(582,326),(669,342),(762,362),(860,393)],
     ("0.85","10%"):[(233,338),(378,368),(459,372),(581,389),(669,402),(764,418),(860,460)],
     ("0.85","5%"):[(228,387),(291,396),(371,412),(454,416),(578,419),(667,438),(763,452),(858,485)]}
with open(os.path.join(HERE,"fig19_volumetric_strain_vs_sigma3.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Fukushima & Tatsuoka (1984) Fig. 19: volumetric strain v [%] (dilation NEGATIVE, as on the paper's axis) at eps_a = 5 % and 10 % vs sigma3' corrected by Method III (log axis).\n")
    f.write("# Markers read by eye (+-0.1 % in v, +-4 % in sigma3').  The paper notes the points below ~0.1 kgf/cm2 show slightly more dilation than the smooth trend (dashed) because of strain hardening before shear.\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["e0p3","eps_a_pct","sigma3_kgf_cm2","sigma3_kPa","eps_v_pct"])
    for (e,ea),pts in F19.items():
        for px,py in pts:
            s,v=f19(px,py); w.writerow([e,ea.rstrip('%'),round(s,4),round(s*KG,2),round(v,2)])
# ---- Fig. 20
def f20(px,py): return (10**((px-732)/322.0), (622-py)/125.5)
F20={("0.70","4"):[(283,521),(345,536),(430,508),(520,487),(639,445),(736,402),(832,355),(930,140)],
     ("0.70","3"):[(243,592),(336,601),(424,590),(517,580),(637,568),(731,558),(829,547),(924,513)],
     ("0.85","3"):[(245,473),(335,484),(423,486),(515,485),(640,430),(736,388),(832,321),(932,125)]}
with open(os.path.join(HERE,"fig20_axial_strain_at_sr_vs_sigma3.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Fukushima & Tatsuoka (1984) Fig. 20: axial strain eps_a [%] at which sigma1'/sigma3' reaches 3 or 4 (mid-height, corrected by Method III) vs sigma3' at that point (log axis).\n")
    f.write("# Markers read by eye (+-0.05 % in eps_a, +-4 % in sigma3').  e0.85 sr=4 is not plotted in the paper.\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["e0p3","sr_level","sigma3_kgf_cm2","sigma3_kPa","eps_a_pct"])
    for (e,sr),pts in F20.items():
        for px,py in pts:
            s,v=f20(px,py); w.writerow([e,sr,round(s,4),round(s*KG,2),round(v,3)])
print("ok")
