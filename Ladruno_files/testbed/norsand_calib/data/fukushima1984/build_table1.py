"""Fukushima & Tatsuoka (1984) Soils and Foundations 24(4):30-48, Table 1 (journal p. 32 = PDF page 3), transcribed by eye from a
native-resolution crop of the scan (every number checked twice against the page).  Columns of the source:
(2) sigma_c' [kgf/cm2] consolidation effective stress, (3) (sigma3')_f [kgf/cm2] at (s1'/s3')max, uncorrected for membrane forces,
(4) e_0.3 = void ratio at sigma_c' = 0.3 kgf/cm2, (5) t0 [mm] initial membrane thickness, (6) phi [deg] UNCORRECTED for membrane forces
(mid-height stress ratio).  1 kgf/cm2 = 98 kN/m2 (paper's own conversion).
    python Ladruno_files/testbed/norsand_calib/data/fukushima1984/build_table1.py
"""
import csv, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
# test, sigma_c', (sigma3')_f (None = not measured), e0.3, t0, phi
T = """
1 0.02 0.060 0.680 0.1 44.3
2 0.02 0.038 0.689 0.1 43.9
3 0.02 0.057 0.714 0.1 42.0
4 0.02 0.048 0.722 0.1 41.5
5 0.02 0.033 0.754 0.1 40.3
6 0.02 0.034 0.769 0.1 39.8
7 0.02 0.030 0.833 0.1 37.5
8 0.02 0.026 0.900 0.1 35.6
9 0.02 0.038 0.679 0.3 47.2
10 0.02 0.033 0.755 0.3 43.4
11 0.02 0.051 0.853 0.3 38.4
12 0.05 0.070 0.687 0.1 43.7
13 0.05 0.061 0.722 0.1 42.4
14 0.05 0.064 0.782 0.1 39.0
15 0.05 0.063 0.833 0.1 37.4
16 0.05 0.058 0.908 0.1 35.4
17 0.05 0.066 0.687 0.12 44.6
18 0.05 0.065 0.761 0.12 41.1
19 0.05 0.067 0.681 0.3 45.4
20 0.05 0.065 0.739 0.3 42.4
21 0.05 0.056 0.898 0.3 37.9
22 0.10 0.117 0.686 0.1 43.7
23 0.10 0.115 0.758 0.1 40.4
24 0.10 0.103 0.902 0.1 35.4
25 0.10 0.118 0.687 0.12 43.9
26 0.10 0.117 0.740 0.12 41.2
27 0.10 0.107 0.902 0.12 35.6
28 0.10 0.120 0.650 0.3 46.0
29 0.10 0.122 0.687 0.3 44.9
30 0.10 0.108 0.696 0.3 44.3
31 0.10 0.110 0.773 0.3 41.1
32 0.10 0.110 0.824 0.3 39.2
33 0.10 0.118 0.866 0.3 38.0
34 0.20 0.213 0.658 0.3 43.9
35 0.20 0.215 0.664 0.3 45.0
36 0.20 0.218 0.669 0.3 44.5
37 0.20 0.215 0.718 0.3 42.3
38 0.20 0.213 0.734 0.3 40.9
39 0.20 0.213 0.746 0.3 40.5
40 0.20 0.212 0.771 0.3 39.8
41 0.20 0.210 0.775 0.3 39.7
42 0.20 0.211 0.827 0.3 38.1
43 0.20 0.202 0.885 0.3 37.1
44 0.20 0.209 0.888 0.3 36.6
45 0.50 0.514 0.658 0.3 43.3
46 0.50 0.514 0.693 0.3 42.2
47 0.50 0.513 0.727 0.3 40.3
48 0.50 0.513 0.750 0.3 39.5
49 0.50 0.512 0.788 0.3 38.0
50 0.50 0.512 0.831 0.3 36.1
51 0.50 0.511 0.881 0.3 34.4
52 1.0 1.014 0.653 0.3 42.4
53 1.0 1.015 0.671 0.3 42.3
54 1.0 1.013 0.723 0.3 40.1
55 1.0 1.013 0.730 0.3 40.1
56 1.0 1.013 0.731 0.3 39.7
57 1.0 1.015 0.760 0.3 38.4
58 1.0 1.013 0.764 0.3 38.0
59 1.0 1.013 0.803 0.3 36.4
60 1.0 1.012 0.806 0.3 36.9
61 1.0 1.012 0.811 0.3 36.8
62 1.0 1.012 0.829 0.3 35.5
63 1.0 1.012 0.853 0.3 35.2
64 1.0 1.012 0.878 0.3 34.4
65 2.0 2.013 0.655 0.3 41.5
66 2.0 2.012 0.677 0.3 40.5
67 2.0 2.014 0.713 0.3 39.4
68 2.0 2.014 0.720 0.3 39.7
69 2.0 2.013 0.749 0.3 38.1
70 2.0 2.013 0.781 0.3 37.0
71 2.0 2.012 0.829 0.3 35.3
72 2.0 2.010 0.876 0.3 34.1
73 4.0 NA 0.679 0.3 39.4
74 4.0 NA 0.732 0.3 37.9
75 4.0 NA 0.751 0.3 36.5
76 4.0 NA 0.782 0.3 35.6
77 4.0 NA 0.835 0.3 34.7
78 4.0 NA 0.901 0.3 32.9
"""
KG = 98.0
rows = []
for ln in T.strip().splitlines():
    t, sc, s3f, e, t0, phi = ln.split()
    sc = float(sc); s3 = None if s3f == "NA" else float(s3f)
    phi = float(phi)
    s = math.sin(math.radians(phi))
    rows.append(dict(test=int(t), sigma_c_kPa=round(sc*KG, 2), sigma3_f_kPa=("" if s3 is None else round(s3*KG, 2)),
                     sigma3_used_kPa=round((s3 if s3 is not None else sc)*KG, 2), sigma3_used_note=("measured (sigma3')_f" if s3 is not None else "sigma_c' (not measured with the HC-DPT)"),
                     e_0p3=float(e), t0_mm=float(t0), phi_deg_uncorrected=phi, sr_derived=round((1+s)/(1-s), 4)))
with open(os.path.join(HERE, "table1_tests.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Fukushima & Tatsuoka (1984) S&F 24(4):30-48, Table 1 (p. 32), 78 drained TC tests on air-pluviated saturated Toyoura sand (e_max 0.977, e_min 0.605, Gs 2.64, D50 0.16 mm, Uc 1.46).\n")
    f.write("# phi_deg_uncorrected = from the MID-HEIGHT stress ratio, NOT corrected for membrane forces (the paper's Figs 14-17 give corrected values only as plotted points / two averaged curves).\n")
    f.write("# e_0p3 = void ratio at sigma_c' = 0.3 kgf/cm2 = 29.4 kPa (for tests sheared below 29.4 kPa it is ESTIMATED from the consolidation curve, p. 35); the void ratio at the shear stress is NOT reported.\n")
    f.write("# sigma_c_kPa, sigma3_f_kPa = kgf/cm2 x 98.  sr_derived = (1+sin phi)/(1-sin phi) (derived here).  Transcribed by eye from the scan; tolerance: none expected on the digits (printed numbers).\n")
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
# harness point-test schema (sigma3_kPa, e, phi_peak_deg, eps_peak_pct, kind) + provenance columns
with open(os.path.join(HERE, "point_tests_tc.csv"), "w", newline="", encoding="utf-8") as f:
    f.write("# Harness point-test schema (harness/data.py load_points).  sigma3_kPa = (sigma3')_f x 98 (sigma_c' for sigma_c' = 4.0 kgf/cm2, not measured).  eps_peak_pct empty (not tabulated).\n")
    f.write("# phi_peak_deg is UNCORRECTED for membrane forces: biased LOW at small sigma3' and large t0 (paper Figs 7, 14-17).  Weight accordingly; t0_mm and test carry that information.\n")
    w = csv.writer(f, lineterminator="\n"); w.writerow(["sigma3_kPa", "e", "phi_peak_deg", "eps_peak_pct", "kind", "test", "t0_mm", "sigma_c_kPa"])
    for r in rows:
        w.writerow([r["sigma3_used_kPa"], r["e_0p3"], r["phi_deg_uncorrected"], "", "TC", r["test"], r["t0_mm"], r["sigma_c_kPa"]])
print(len(rows), "tests")
