"""Round 3b: single-step pattern map around the K1.14b state (HAR TIMs, p = -0.6 kPa, eta = 1.2 M, dry side):
which (expansion, shear) increments give FPf (post floor on the dry side) and which give FP- (the dilative
volumetric term wins). Run: python -u fpf_region_har.py"""
import math
import numpy as np
import har_patch as HP
from har_patch import K, ONES
from conftest import make_params

H = HP.HarParams(**HP.TIMS_HAR); HP.install(H); PMIN = 5e-3 * H.p_a
kw = dict(p0=-H.p_a, M=1.3309, N=0.4, N_bar=0.2, chi=-3.5, h=280.0, rho=0.71, rho_bar=0.71, zeta="WW",
          csl_mode="fork", e0=0.83, lam_c=0.027, xi=0.45, p_a=H.p_a, cap="none")
P = make_params("O2", **kw)
sig, th, nh = HP.off_corner_sig(P, -0.6, 1.2 * P.M); pi = K.pi_of_eta(P, -0.6, 1.2 * P.M)
eps = K.invert_elastic(P, sig); v = 1.0 + P.e0 - P.lam_c * (-pi / P.p_a) ** P.xi - 0.10
def one(d):
    vn = v * math.exp(float(d.sum()))
    etf, at, _, _ = HP.floor_op(H, eps + d, PMIN)
    res = K.return_map(P, etf, pi, vn, vn)
    if res.refused: return "REF", float("nan")
    ef, ap, _, _ = HP.floor_op(H, res.eps_e, PMIN)
    return ("F" if at else "-") + ("P" if res.plastic else "E") + ("f" if ap else "-"), K.elastic(P, res.eps_e).p
avs = (1e-5, 2e-5, 4e-5, 8e-5, 1.6e-4, 3.2e-4)
ass = (0.0, 5e-6, 1e-5, 2e-5, 4e-5, 8e-5, 1.6e-4, 3.2e-4)
print("single-step pattern (p_c) at the K1.14b state; rows tr deps (expansion), columns shear amplitude along n^")
print("tr deps | shear " + "".join(f"{a:>16.0e}" for a in ass))
for av in avs:
    row = []
    for a in ass:
        pat, pc = one((av / 3.0) * ONES + a * nh)
        row.append(f"{pat} ({pc:+.4f})")
    print(f"{av:>14.0e}  " + "".join(f"{r:>16}" for r in row))
HP.uninstall()
