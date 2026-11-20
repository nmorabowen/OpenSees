"""Digitized Tatsuoka et al. (1986) S&F 26(1):65-84 data (delta = 90 deg, isotropically consolidated, air-pluviated,
saturated Toyoura; stresses by method C-2-T; axial strain from boundary (external) displacements).
Axis calibration: from frame lines/ticks at 300 dpi (see digit notes in the report). ratio = sigma1'/sigma3';
ev_pct: volumetric strain in %, dilation NEGATIVE (the figures' convention).
Resolution: 1 px = 0.019-0.022 % axial strain, 0.011 in ratio.
"""
import json
import numpy as np

f4a = json.load(open("data_f4a.json"))
# PSL-24: clean the automatic circle detections (drop the spurious 'eps_v' label hit and the overlapping-circle duplicates)
sr = [p for p in f4a["ratio"] if not (abs(p[0] - 0.037) < 1e-3 and abs(p[1] - 2.259) < 1e-3)
      and not (abs(p[0] - 0.093) < 1e-3) and not (abs(p[0] - 2.796) < 1e-3) and not (abs(p[0] - 0.815) < 1e-3 and p[1] > 5.845)
      and not (abs(p[0] - 2.185) < 1e-3 and p[1] > 5.68)]
pre = sorted([p for p in sr if p[0] <= 1.45], key=lambda p: p[1])
post = sorted([p for p in sr if p[0] > 1.45])
# pre-peak: monotone in ratio; strain made non-decreasing (overlapping circles give +-1 px jitter)
m = 0.0
pre2 = []
for x, y in pre:
    m = max(m, x)
    pre2.append((m, y))
ev = [p for p in f4a["ev_pct"] if not (abs(p[0] - 1.833) < 1e-3)]
ev = [(0.0, 0.0)] + [p for p in ev if p[0] > 0.05]

TESTS = {
    "PSL-24 (Fig. 4a)": dict(fig="4(a)", page=71, s3=4.9, e=0.714, marker="circles, auto-detected",
                            ratio=[(0.0, 1.0)] + pre2[1:] + post, ev=ev),
    "Fig. 6(a) delta90": dict(fig="6(a)", page=73, s3=4.9, e=0.700, marker="solid line, read on a 0.05-0.1 % grid",
                              ratio=[(0, 1.0), (0.03, 3.0), (0.11, 3.5), (0.32, 4.0), (0.50, 4.5), (0.68, 5.0), (0.86, 5.5),
                                     (1.15, 6.0), (1.47, 6.5), (1.74, 6.52), (2.45, 5.34), (2.77, 4.78), (3.10, 4.51),
                                     (3.64, 4.34), (4.18, 4.29), (4.93, 4.26), (5.90, 4.34), (6.90, 4.48), (7.97, 4.62),
                                     (8.83, 4.68)],
                              ev=[(0, 0), (0.5, 0.05), (1.0, -0.3), (1.47, -0.72), (2.34, -1.54), (3.10, -1.92), (4.18, -2.02),
                                  (5.26, -2.14), (6.34, -2.30), (7.4, -2.38), (8.9, -2.40)]),
    "Fig. 16(a) iso": dict(fig="16(a)", page=78, s3=49.0, e=0.716, marker="solid line, read on a 0.05-0.1 % grid",
                           ratio=[(0, 1.0), (0.17, 2.72), (0.26, 3.0), (0.60, 4.0), (1.01, 5.0), (1.30, 5.5), (1.71, 6.0),
                                  (2.05, 6.32), (2.62, 6.0), (3.02, 5.5), (3.22, 5.0), (3.48, 4.5), (3.95, 4.25), (5.0, 4.18),
                                  (7.0, 4.15), (9.0, 4.2), (11.9, 4.22)],
                           ev=[(0, 0), (1.0, 0.0), (1.98, -0.55), (3.06, -1.28), (3.9, -1.54), (6.3, -1.76), (8.4, -1.87),
                               (11.9, -2.04)]),
    "Fig. 8 s3=4.9": dict(fig="8", page=73, s3=4.9, e=0.755, marker="solid line, read on a 0.05-0.1 % grid",
                          ratio=[(0, 1.0), (0.066, 2.0), (0.22, 3.0), (0.50, 4.0), (0.96, 5.0), (1.33, 5.5), (1.99, 5.76),
                                 (3.5, 4.8), (7.0, 4.55), (9.0, 4.5), (12.5, 4.3)],
                          ev=[]),
}

# Fig. 22 (p. 82): axial strain at (s1/s3)max, delta 90, by sigma_c' [kgf/cm2] -> (e0.05, eps_peak %)
FIG22 = {0.05: [(0.674, 1.14), (0.677, 1.79), (0.700, 2.00), (0.715, 1.46), (0.755, 2.04), (0.800, 2.66), (0.861, 3.79)],
         0.1: [(0.670, 1.38), (0.709, 2.05), (0.746, 2.30), (0.805, 2.94)],
         0.5: [(0.669, 1.81), (0.717, 2.14), (0.741, 2.16), (0.810, 2.50)],
         1.0: [(0.653, 2.26), (0.691, 2.59), (0.721, 2.39), (0.744, 2.96), (0.775, 3.34), (0.821, 3.59)],
         4.0: [(0.674, 3.61), (0.713, 3.92), (0.752, 4.44), (0.821, 5.43)]}
# Fig. 9 (p. 75): phi_peak (method C-2-T), delta 90 -> (e0.05, phi deg)
FIG9 = {0.05: [(0.672, 49.9), (0.674, 48.6), (0.700, 47.4), (0.755, 44.7), (0.800, 41.3), (0.862, 37.5)],
        0.1: [(0.672, 49.9), (0.708, 47.8)],
        0.5: [(0.670, 47.7), (0.717, 46.7), (0.740, 45.3), (0.808, 41.6)],
        1.0: [(0.653, 49.0), (0.692, 47.3), (0.721, 45.6), (0.747, 44.1), (0.774, 42.5), (0.821, 40.7)],
        4.0: [(0.674, 46.1), (0.713, 44.3), (0.752, 42.3), (0.820, 39.6)]}

if __name__ == "__main__":
    json.dump(dict(TESTS=TESTS, FIG22={str(k): v for k, v in FIG22.items()}, FIG9={str(k): v for k, v in FIG9.items()}),
              open("tatsuoka_digitized.json", "w"), indent=1)
    for k, v in TESTS.items():
        print(k, len(v["ratio"]), len(v["ev"]))
