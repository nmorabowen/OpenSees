"""WP-150 T6: classical capacity band for the WP-138 strip deck, q_u(phi) = 1/2 gamma' B N_gamma + q N_q.

N_gamma: EXACT rough-strip values by the method of characteristics (Martin 2005; reproduced by Han et al. 2016,
SpringerPlus 5:1482, Table 2): 14.8 / 34.5 / 85.6 / 234 at 30 / 35 / 40 / 45 deg (associated flow). Values between
the four sourced angles are LOG-INTERPOLATED and flagged. No value is extrapolated beyond 45 deg.
N_q: exact Prandtl-Reissner, exp(pi tan phi) tan^2(45 + phi/2).
Deck (WP-138): gamma' 9.81 kN/m3, B 1.5 m, surcharge 7.65 kPa outside the footprint. Superposition is the usual
approximation (conservative for rough footings).
Martin also shows the NON-associated N_gamma (psi < phi) is lower: an upper band for a psi = 0 DP.
"""
import math

MARTIN = {30: 14.8, 35: 34.5, 40: 85.6, 45: 234.0}
import sys
kw = dict(a.split("=") for a in sys.argv[1:])
GAMMA, B, Q = float(kw.get("gamma", 9.81)), float(kw.get("B", 1.5)), float(kw.get("q", 7.65))
PR = float(kw.get("pr", 0.0))      # -Presidual: a shift of the MC envelope by p_r -> extra capacity p_r (N_q - 1)


def n_gamma(phi):
    if phi in MARTIN:
        return MARTIN[phi], "Martin"
    ks = sorted(MARTIN)
    if not ks[0] < phi < ks[-1]:
        return float("nan"), "out of the sourced range"
    lo = max(k for k in ks if k < phi); hi = min(k for k in ks if k > phi)
    w = (phi - lo) / (hi - lo)
    return math.exp((1 - w) * math.log(MARTIN[lo]) + w * math.log(MARTIN[hi])), "log-interp"


def n_q(phi):
    t = math.tan(math.radians(phi))
    return math.exp(math.pi * t) * math.tan(math.radians(45 + phi / 2)) ** 2


print(f"deck: gamma {GAMMA} kN/m3, B {B} m, surcharge {Q} kPa, p_r {PR} kPa")
print()
print("| φ′ ° | N_γ (rough, exact) | N_q | q_u = ½γ′BN_γ + qN_q kPa | + p_r (N_q − 1) kPa | note |")
print("|---|---|---|---|---|---|")
for phi in (30, 33, 35, 38, 40, 42, 43, 44, 45):
    ng, src = n_gamma(phi)
    nq = n_q(phi)
    qu = 0.5 * GAMMA * B * ng + Q * nq
    print(f"| {phi} | {ng:.1f} | {nq:.1f} | {qu:.0f} | {PR * (nq - 1):.0f} | {src} |")
