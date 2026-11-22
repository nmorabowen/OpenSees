"""Chakraborty & Salgado (2009) 'Dilatancy and shear strength behavior of sand at low confining pressures', Proc. 17th ICSMGE (Alexandria),
IOS Press, pp. 652-655, doi 10.3233/978-1-60750-031-5-652.  PDF: Dropbox/SOILS_rev/WP144_calibration/Chakraborty_Salgado_2009_ISSMGE_STAL0652.pdf

Tables 1 and 2 (journal p. 654 = PDF page 4) transcribed from a 250-dpi render (the text layer of the PDF scrambles Table 2, so the numbers were read from the
image and compared with the text layer for Table 1).  Equations and the regression lines of Fig. 3 are from pp. 652-655.

    python Ladruno_files/testbed/norsand_calib/data/chakraborty2009/build_tables.py
"""
import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# sigma3p [kPa], sigma_mp [kPa], data points, Q_best, R_best, r2_best, Q(R=1), R, r2(R=1)
TX = [
    (4.0, 9.3, 11, 6.9, 0.47, 0.928, 7.7, 1, 0.914),
    (6.2, 14.3, 10, 6.2, -0.23, 0.943, 8.1, 1, 0.839),
    (11.2, 25.8, 12, 7.4, 0.13, 0.99, 8.7, 1, 0.954),
    (20.8, 47.2, 11, 7.5, 0.03, 0.987, 9.0, 1, 0.945),
    (50.3, 108.4, 7, 8.9, 0.79, 0.999, 9.3, 1, 0.997),
    (99.3, 207.5, 13, 9.3, 0.8, 0.997, 9.7, 1, 0.996),
    (197.2, 412.4, 8, 9.6, 0.73, 0.999, 10.0, 1, 0.997),
]
PSC = [  # + b = 0.25 assumed
    (4.9, 15.7, 10, 0.25, 9.2, 1.57, 0.971, 8.4, 1, 0.963),
    (9.8, 29.7, 5, 0.25, 10.4, 2.1, 0.985, 8.8, 1, 0.961),
    (49.0, 151.2, 4, 0.25, 10.2, 0.88, 0.998, 10.3, 1, 0.998),
    (68.6, 196.7, 6, 0.25, 10.1, 0.9, 0.993, 10.2, 1, 0.993),
    (98.0, 297.1, 6, 0.25, 11.0, 1.2, 0.992, 10.7, 1, 0.99),
]
HDR_COMMON = ("# Chakraborty & Salgado (2009) ICSMGE 17, pp. 652-655, {tab}, Toyoura sand ({kind}). I_R = I_D (Q - ln(100 sigma'_mp/p_A)) - R (eq. 3), p_A = 100 kPa, I_D = D_R/100;\n"
              "# fitted by eq. 9 to the peak-strength data of each sigma'_3p level (relative densities 30-90 %).  Transcribed by eye from the PDF page image; no tolerance expected on the digits.\n")
with open(os.path.join(HERE, "table1_tx_Q_R.csv"), "w", newline="", encoding="utf-8") as f:
    f.write(HDR_COMMON.format(tab="Table 1", kind="triaxial compression, Fukushima & Tatsuoka 1984 and Tatsuoka 1987 data; phi_c = 32.8 deg"))
    f.write("# sigma3p = minor principal effective stress at peak strength; sigma_mp = mean effective stress at peak (eq. 10: (s1+2 s3)/3).  'R=1' = trend line set with R = 1.\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["sigma3p_kPa", "sigma_mp_kPa", "data_points", "Q_best", "R_best", "r2_best", "Q_R1", "R_R1", "r2_R1"])
    w.writerows(TX)
with open(os.path.join(HERE, "table2_psc_Q_R.csv"), "w", newline="", encoding="utf-8") as f:
    f.write(HDR_COMMON.format(tab="Table 2", kind="plane-strain compression, Tatsuoka et al. 1986 and Tatsuoka 1987 data; phi_c = 36 deg"))
    f.write("# sigma_mp from eq. 11 ((s1+s2+s3)/3) with b = (s2-s3)/(s1-s3) assumed 0.25 (paper: b = 0.2-0.3 at peak).\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["sigma3p_kPa", "sigma_mp_kPa", "data_points", "b", "Q_best", "R_best", "r2_best", "Q_R1", "R_R1", "r2_R1"])
    w.writerows(PSC)

PARAMS = [
    # name, TX, PSC, where, note
    ("phi_c_deg", 32.8, 36.0, "p. 653 sect. 3", "critical-state friction angle chosen by the authors: TX 32.8 (average of 31.2-34.4 from Fukushima & Tatsuoka 1984 and Tatsuoka 1987; others quoted 31.6 Verdugo & Ishihara, 31.1 Wang et al.), PSC 36 (quoted range 34.5-38)"),
    ("A_psi", 3.8, 3.8, "Fig. 1(b), p. 653, eq. 8", "phi_p = phi_c + A_psi I_R fitted to ALL TX and PSC Toyoura data (the same 3.8 for both; Bolton's 3 TX / 5 PSC), i.e. I_R = (-d eps_v/d eps_1)/0.3 = 0.26 (phi_p - phi_c)"),
    ("dilatancy_relation", "phi_p = phi_c + 0.6 psi_p", "phi_p = phi_c + 0.6 psi_p", "eq. 7", "combining A_psi = 3.8 with eq. 2 for both TX and PSC"),
    ("I_R", "I_D (Q - ln(100 sigma_mp'/p_A)) - R", "same", "eq. 3", "p_A = 100 kPa; Bolton: Q = 10, R = 1; I_R capped at 4 by Bolton (no low-stress data)"),
    ("Q_trend_R1", "0.60 ln(sigma_c') + 7.4", "0.75 ln(sigma_c') + 7.1", "Fig. 3, p. 655", "sigma_c' = INITIAL confining stress in kPa (not sigma3p of the tables); R = 1; the tables' Q(R=1) rises 7.7-10.0 (TX) and 8.4-10.7 (PSC)"),
    ("Q_range_R1", "7.7-10.0", "8.4-10.7", "Tables 1, 2", "r2 0.839-0.997 (TX), 0.961-0.998 (PSC)"),
    ("confining_stress_range_kPa", "4-197 (sigma3p)", "4.9-98 (sigma3p)", "Tables 1, 2", "relative densities 30-90 %"),
    ("b_assumed_PSC", "n/a", 0.25, "p. 654", "only used for sigma_mp"),
]
with open(os.path.join(HERE, "fit_parameters.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Chakraborty & Salgado (2009): strength-dilatancy correlation parameters for Toyoura sand.  Source of the data fitted: Fukushima & Tatsuoka 1984 (TX), Tatsuoka et al. 1986 and Tatsuoka 1987 (PSC).\n")
    w = csv.writer(f, lineterminator="\n")
    w.writerow(["parameter", "TX", "PSC", "location", "note"])
    w.writerows(PARAMS)
print("ok")
