"""Lam & Tatsuoka (1988) Soils and Foundations 28(1):89-106, Fig. 21 (p. 101 = PDF page 13): phi of air-pluviated Toyoura sand,
e0.3 = 0.70, sigma3' = 1.0 kgf/cm2 = 98 kPa (TE corrected to sigma3' = 98 kPa and for local area reduction, text p. 98-99), vs
(a) the orientation of the bedding plane for TC (H/W 1.0), PSC (H/W 1.9) and TE (Fig. 16 values) and (b) b = (s2-s3)/(s1-s3) at failure.
Values were read from native-resolution crops with a data-coordinate grid overlay (phi axis 50 deg at the frame top, 30 deg at the frame bottom;
categorical x positions of panel (a) read at the tick marks).  Reading tolerance: phi +-0.25 deg, b +-0.02.
    python Ladruno_files/testbed/norsand_calib/data/lam_tatsuoka1988/fig21_phi_tables.py
"""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
# ---------------------------------------------------------------- panel (a): phi vs orientation
# categories k0..k9 and the (omega, xi) they denote (axis labels of Fig. 21(a): sector 1 'xi at omega=90', sector 2 'omega at xi=90', sector 3 'omega at xi=0'
# (the printed label of sector 3 reads 'xi = 90'; TC invariance (TC-3 values equal sector 2) and 'theta = 0 for TE' show it is xi = 0))
CAT = [(90,0),(90,30),(90,60),(90,90),(60,90),(30,90),(0,90),(30,0),(60,0),(90,0)]   # (omega, xi) for TC/PSC
CATLAB = ["xi=0@omega=90","xi=30@omega=90","xi=60@omega=90","xi=90@omega=90","omega=60@xi=90","omega=30@xi=90","omega=0@xi=90","omega=30@xi=0","omega=60@xi=0","omega=90@xi=0"]
A = {
 "TC-3 (H/W=1.0)": [38.52,38.31,38.34,38.34,38.34,39.76,41.74,39.85,38.47,38.39],
 "PSC-3 (H/W=1.9)":[43.70,43.98,41.28,40.00,39.70,42.78,45.32,44.49,43.69,43.76],
 "TE double (Fig.16)":[48.38,46.34,42.00,43.50,42.00,46.40,48.40,48.40,48.40,48.40],
 "TE single (Fig.16)":[46.50,43.70,38.97,None,38.97,43.77,46.53,46.50,46.50,46.50],
 "TE no (Fig.16)":[42.80,None,36.90,38.30,36.90,None,42.90,42.90,42.90,42.90],
}
with open(os.path.join(HERE,"fig21a_phi_vs_orientation.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Lam & Tatsuoka (1988) Fig. 21(a), p. 101.  phi [deg] of air-pluviated Toyoura, e0.3 = 0.70, sigma3' = 98 kPa.  tolerance +-0.25 deg (reading).\n")
    f.write("# TC/PSC: omega = angle of sigma1' from the deposition axis n; xi = angle of the n-axis projection (Fig. 1).  For TE the abscissa is theta (sin theta = sin omega sin xi): sectors 1-2 theta = xi / omega, sector 3 theta = 0 (flat lines).\n")
    f.write("# Empty = marker hidden/not legible.  PSC-1 (earlier Tatsuoka et al. 1986 data, omega 56-67 deg) not transcribed.\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["series","category_k","category_label","omega_deg","xi_deg","phi_deg"])
    for s,vals in A.items():
        for k,v in enumerate(vals):
            if v is None: continue
            w.writerow([s,k,CATLAB[k],CAT[k][0],CAT[k][1],v])
# ---------------------------------------------------------------- panel (b): phi vs b
# (series, omega, xi, mode, b, phi)
B=[
 ("A: xi=90, omega varies",0,90,"TC",0.00,41.86),("A: xi=90, omega varies",0,90,"PSC",0.25,45.36),("A: xi=90, omega varies",0,90,"TE single (theta=0)",1.00,46.62),
 ("A: xi=90, omega varies",30,90,"TC",0.00,39.72),("A: xi=90, omega varies",30,90,"PSC",0.25,42.88),("A: xi=90, omega varies",30,90,"TE single (theta=30)",1.00,43.94),
 ("A: xi=90, omega varies",60,90,"TC",0.00,38.31),("A: xi=90, omega varies",60,90,"PSC",0.24,39.72),("A: xi=90, omega varies",60,90,"TE single (theta=60)",1.00,38.90),
 ("A: xi=90, omega varies",90,90,"TC",0.01,38.34),("A: xi=90, omega varies",90,90,"PSC",0.29,40.20),("A: xi=90, omega varies",90,90,"TE single (theta=90)",1.00,41.24),
 ("B: xi=0, omega varies",30,0,"PSC",0.28,44.31),("B: xi=0, omega varies",60,0,"PSC",0.32,43.60),("B: xi=0, omega varies",90,0,"PSC",0.34,43.65),
 ("C: omega=90, xi varies",90,30,"PSC",0.36,44.10),("C: omega=90, xi varies",90,60,"PSC",0.29,41.30),
]
with open(os.path.join(HERE,"fig21b_phi_vs_b.csv"),"w",newline="",encoding="utf-8") as f:
    f.write("# Lam & Tatsuoka (1988) Fig. 21(b), p. 101: phi at failure vs b = (s2'-s3')/(s1'-s3'), SINGLE-intersection failure mode for TE, e0.3 = 0.70, sigma3' = 98 kPa.\n")
    f.write("# Group A = legend column 1 (omega varies, xi = 90), B = column 2 (xi = 0), C = column 3 (omega = 90, xi varies).  b = 0 is TC, b = 1 is TE (b at PSC failure 0.24-0.36, not 0.3 exactly).\n")
    f.write("# Points at b = 0 and b = 1 of groups B and C coincide with group A or are the same stress state and are not repeated.  Reading tolerance: phi +-0.3 deg, b +-0.02 (PSC), +-0.01 (TC/TE).\n")
    w=csv.writer(f,lineterminator="\n"); w.writerow(["series","omega_deg","xi_deg","mode","b","phi_deg"])
    w.writerows(B)
print("ok")
