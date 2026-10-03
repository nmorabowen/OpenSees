"""WP-130: BackwardEuler_CPPM on WP-128's smallest reproducer.

WP-128 (128_sanisand_ring_trace.md §2.3): sigma = p_s*I, alpha = alpha_in =
z = 0, one plane-strain d_eps = (0, delta, 0, 0, 0, 0) (compression-positive).
Today's ModifiedEuler (IntScheme 1, the campaign set) returns alpha/alpha^b =
5.14 at p_s = 0.0101, delta = 1e-4 in ONE accepted substep, rc 0.  The
alpha-aware port / RK45 give ~0.25-0.27.

Here the same replay (`ladrunoSANISANDReplay`, WP-127) is run with the
campaign material under IntScheme 1 and under IntScheme 2 (vanilla CPPM and
the WP-130 variants), reporting rc, eta, alpha/alpha^b, the CPPM census.

alpha/alpha^b = sqrt(3/2)*||alpha|| / alpha^b_theta with alpha^b_theta =
g(theta)*Mc*exp(-nb*psi) - m evaluated at the RETURNED state (the model's own
GetStateDependent), g from the Lode angle of n = (s - p alpha)/|s - p alpha|.

usage: python -S <bootstrap> q128_reproducer.py   (WP-130 dist on the path)
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "Ladruno_scripts"))
import sanisand_replay as sr  # noqa: E402
import opensees as ops  # noqa: E402

P = dict(zip(["G0", "nu", "e_init", "Mc", "c", "lambda_c", "e0", "ksi", "P_atm", "m",
              "h0", "ch", "nb", "A0", "nd", "z_max", "cz", "rho"], sr.CAMPAIGN_PARAMS))


def _t(v):   # Voigt (xx, yy, zz, xy, yz, zx) -> 3x3
    return [[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]]


def ratio(res):
    s, a = res["sigma"], res["alpha"]
    p = (s[0] + s[1] + s[2]) / 3.0
    dev = [s[i] - (p if i < 3 else 0.0) for i in range(6)]
    r = [dev[i] - p * a[i] for i in range(6)]
    nr = math.sqrt(sum(r[i] ** 2 for i in range(3)) + 2 * sum(r[i] ** 2 for i in range(3, 6)))
    n = [x / nr for x in r] if nr > 0 else [0.0] * 6
    N = _t(n)
    n3 = sum(N[i][j] * N[j][k] * N[k][i] for i in range(3) for j in range(3) for k in range(3))
    c3 = max(-1.0, min(1.0, math.sqrt(6.0) * n3))
    c = P["c"]
    g = 2 * c / ((1 + c) - (1 - c) * c3)
    pp = max(p, 1e-10)
    psi = res["e"] - (P["e0"] - P["lambda_c"] * (pp / P["P_atm"]) ** P["ksi"])
    ab = g * P["Mc"] * math.exp(-P["nb"] * psi) - P["m"]
    na = math.sqrt(sum(a[i] ** 2 for i in range(3)) + 2 * sum(a[i] ** 2 for i in range(3, 6)))
    eta = math.sqrt(1.5) * (math.sqrt(sum(dev[i] ** 2 for i in range(3)) + 2 * sum(dev[i] ** 2 for i in range(3, 6)))) / pp
    return eta, math.sqrt(1.5) * na / ab


def material(tag, scheme, extra):
    opts = list(sr.CAMPAIGN_OPTS)
    opts[0] = scheme
    ops.nDMaterial("LadrunoSANISAND", tag, *sr.CAMPAIGN_PARAMS, *opts, *extra)


ARMS = [("ME (IntScheme 1, campaign)", 1, ()),
        ("CPPM vanilla (IntScheme 2)", 2, ()),
        ("CPPM refuse h0", 2, ("-cppmOnFail", "refuse", "-cppmHalvings", 0)),
        ("CPPM refuse h0 + start + LS", 2, ("-cppmOnFail", "refuse", "-cppmHalvings", 0,
                                           "-cppmStart", "explicit", "-cppmLineSearch", "on")),
        ("CPPM refuse h9 + start + LS", 2, ("-cppmOnFail", "refuse",
                                           "-cppmStart", "explicit", "-cppmLineSearch", "on")),
        ("ME cap 20 + -meFallback cppm", 1, ("-maxSubsteps", 20, "-meFallback", "cppm"))]

CASES = [(0.0101, 1e-5), (0.0101, 3e-5), (0.0101, 1e-4), (0.0101, 3e-4), (0.1, 1e-4),
         (1.0, 3e-4), (5.0, 3e-4)]

print(f"CAMPAIGN_OPTS = {sr.CAMPAIGN_OPTS}")
print("| p_s | delta | arm | rc | eta | alpha/alpha^b | ME substeps | CPPM newtonFail/halv/explFail/lowP/refused | guess ok/tries |")
print("|---|---|---|---|---|---|---|---|---|")
for ps, d in CASES:
    for k, (name, scheme, extra) in enumerate(ARMS):
        ops.wipe()
        ops.model("basic", "-ndm", 3, "-ndf", 3)
        material(1, scheme, extra)
        e = P["e_init"]
        try:
            res = sr.replay(ops, 1, [ps, ps, ps, 0, 0, 0], [0] * 6, [0] * 6, [0] * 6, e,
                            [0, d, 0, 0, 0, 0], "compressionPositive", mat_type="PlaneStrain")
        except Exception as exc:
            print(f"| {ps} | {d:g} | {name} | ERR {exc} |")
            continue
        eta, ar = ratio(res)
        s = res["stats"]
        print(f"| {ps} | {d:g} | {name} | {res['rc']} | {eta:.3f} | {ar:.3f} | {int(s['substeps'])} | "
              f"{int(s['cppmNewtonFail'])}/{int(s['cppmHalvings'])}/{int(s['cppmExplicitFail'])}/"
              f"{int(s['cppmExplicitLowP'])}/{int(s['cppmRefusals'])} | "
              f"{int(s['cppmGuessOk'])}/{int(s['cppmGuessTries'])} |", flush=True)
