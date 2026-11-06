"""WP-144 round 3b, the p' floor SENSITIVITY mini-BVP ("F vs F/2", sheet 9.7 "Determinism, counting, interfaces", plan 2.7, deck report K5).

Question.  The floor p_min is a regularisation (an energy source of at most W_f = p_min eps^f_v per Gauss point, sheet 9.7 item 1); the plan accepts
it for a deck if the limit load moves by < 2 % when p_min is halved.  This module runs that experiment on a SMALL strip so that the report exists
before the big decks do (about a minute per run): a rigid smooth strip footing on a dry sand half-space, plane strain, the TIMs HAR material
(sheet 2.3 / 15) in K0 geostatic state, displacement-controlled.  The first rows of Gauss points sit at p' = 0.1 ... 3 kPa, i.e. BELOW and around the
default floor of 0.505 kPa, which is exactly the regime the floor exists for (the free-surface ring of the Kimura / TIMs decks).

    python floor_sensitivity.py            # prints the report for p_min scales 1, 1/2 and 2 (default floor 0.505 kPa)

Model (all numbers are inputs, none is an expected value of any test):
  * half strip by symmetry: x in [0, 6] m, depth 4.5 m, footing half-width B = 1 m (nodes with x <= B prescribed in y; smooth: x free), roller sides,
    fixed base; 108 `quad ... PlaneStrain` elements on a graded mesh (rows 0.02 m near the surface -> 1.1 m at depth; columns 0.25 -> 1.5 m);
  * dry sand gamma = 15 kN/m^3, K0 = 0.4554 (sheet 9.7: geostatic p ~ 0.64 sigma_v), one material per element row with its mid-depth K0 stress as -sigma0
    (a stepwise-uniform row field is in exact discrete equilibrium with the tributary nodal weights, so the first equilibrium iteration only
    absorbs the floor's own O(p_min) projection); -pi0_auto (the unified rule (S.53));
  * the TIMs elastic set (HAR n = 1/2, g = 807.80, k = 1889.48, p_a = 101), M 1.3309, rho = rho_bar 0.71, fork CSL (e0 0.83, lambda_c 0.027, xi 0.45),
    e_init = 0.6944 (v0 = 1.6944), the K2 plastic constants N 0.4, N_bar 0.2, chi -3.5, h 280 (TIMs' own are a P3 refit);
  * the footing is displaced by delta = 0.1 m (10 % of B) in 40 equal steps with step halving on a failed step (up to 6 halvings).

Reported per run: the load-displacement curve (reaction of the footing nodes per metre of strip, kN/m, half strip), the reaction at the final
displacement and its peak, the number of Gauss points at the floor at the end (`floor` response at_floor), the sum of W_f over the mesh and the
number of floor events, the number of refused steps (must be zero), the wall time.  `compare()` gives the relative change of the limit load.
"""
from __future__ import annotations

import time

import g2_common as G

ops = G.ops

B_HALF = 1.0
XS = [0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0]
DEPTHS = [0.0, 0.02, 0.05, 0.10, 0.20, 0.35, 0.55, 0.80, 1.20, 1.80, 2.50, 3.40, 4.50]
GAMMA, K0 = 15.0, 0.4554
DELTA, NSTEPS, MAX_HALVINGS = 0.10, 40, 6
PA, K_HAR, G_HAR = 101.0, 1889.48104361, 807.80387674
PMIN_DEFAULT = 5.0e-3 * PA                       # 0.505 kPa


def material_args(sigma_v, pmin):
    sig = [-K0 * sigma_v, -sigma_v, -K0 * sigma_v, 0.0, 0.0, 0.0]          # 11 = x, 22 = y (vertical), 33 = out of plane
    return ["-energy", "HAR", "-k", K_HAR, "-g", G_HAR, "-n", 0.5, "-p_a", PA,
            "-M", 1.3309, "-N", 0.4, "-N_bar", 0.2, "-rho", 0.71, "-rho_bar", 0.71, "-chi", -3.5, "-h", 280.0,
            "-csl", "fork", "-e0", 0.83, "-lambda_c", 0.027, "-xi", 0.45, "-pmin", float(pmin),
            "-v0", 1.6944, "-pi0_auto", "-sigma0", *sig]


def node_id(i, j):
    return j * len(XS) + i + 1


def build(pmin):
    """Model, materials, elements, boundary conditions, gravity and footing patterns; returns (footing nodes, element tags)."""
    nx, ny = len(XS), len(DEPTHS)
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    for j in range(ny):
        for i in range(nx):
            ops.node(node_id(i, j), XS[i], -DEPTHS[j])
    eles, weight = [], {}
    tag = 0
    for j in range(ny - 1):
        zc = 0.5 * (DEPTHS[j] + DEPTHS[j + 1])
        ops.nDMaterial("LadrunoNorSand", 100 + j, *material_args(GAMMA * zc, pmin))
        for i in range(nx - 1):
            tag += 1
            n1, n2, n3, n4 = node_id(i, j + 1), node_id(i + 1, j + 1), node_id(i + 1, j), node_id(i, j)
            ops.element("quad", tag, n1, n2, n3, n4, 1.0, "PlaneStrain", 100 + j)
            eles.append(tag)
            area = (XS[i + 1] - XS[i]) * (DEPTHS[j + 1] - DEPTHS[j])
            for n in (n1, n2, n3, n4):
                weight[n] = weight.get(n, 0.0) + GAMMA * area / 4.0
    for i in range(nx):
        ops.fix(node_id(i, ny - 1), 1, 1)                 # base (the corners included)
    for j in range(ny - 1):
        ops.fix(node_id(0, j), 1, 0)                      # symmetry axis
        ops.fix(node_id(nx - 1, j), 1, 0)                 # roller side
    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n, w in weight.items():
        ops.load(n, 0.0, -w)
    ops.timeSeries("Linear", 2)
    ops.pattern("Plain", 2, 2)
    foot = [node_id(i, 0) for i in range(nx) if XS[i] <= B_HALF + 1e-12]
    for n in foot:
        ops.sp(n, 2, -DELTA)
    return foot, eles


def solver(dt):
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("FullGeneral")                     # the tangent is NON-symmetric
    ops.test("NormDispIncr", 1.0e-8, 40, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", dt)
    ops.analysis("Static")


def reaction(foot):
    ops.reactions()
    return -sum(ops.nodeReaction(n, 2) for n in foot)          # downward load carried by the half strip (kN/m)


def run_strip(pmin_scale=1.0, nsteps=NSTEPS, verbose=False):
    """One displacement-controlled run at p_min = pmin_scale x 0.505 kPa; returns the report dict."""
    t0 = time.time()
    pmin = pmin_scale * PMIN_DEFAULT
    foot, eles = build(pmin)
    dt0 = 1.0 / nsteps
    solver(0.0)
    assert ops.analyze(1) == 0, "the geostatic state did not equilibrate"      # absorbs the floor's O(p_min) projection
    r0 = reaction(foot)                                                      # the footing nodes' share of the dead weight: subtracted
    curve, t, dt, refused, failed = [(0.0, 0.0)], 0.0, dt0, 0, 0
    while t < 1.0 - 1e-12:
        dt = min(dt, 1.0 - t)
        ops.integrator("LoadControl", dt)
        halv = 0
        ok = ops.analyze(1) == 0
        while not ok and halv < MAX_HALVINGS:
            failed += 1
            halv += 1
            dt *= 0.5
            ops.integrator("LoadControl", dt)
            ok = ops.analyze(1) == 0
        if not ok:
            refused += 1
            break
        t += dt
        dt = min(dt * 1.5, dt0)
        curve.append((t * DELTA, reaction(foot) - r0))
        if verbose:
            print(f"  t = {t:.4f}  R = {curve[-1][1]:.3f}", flush=True)
    n_at, sum_wf, events = 0, 0.0, 0
    for e in eles:
        for g in range(1, 5):
            fl = list(ops.eleResponse(e, "material", g, "floor"))
            n_at += int(fl[0])
            sum_wf += fl[4]
            events += int(fl[1] + fl[2])
    return dict(pmin=pmin, curve=curve, R_final=curve[-1][1], R_peak=max(r for _, r in curve), reached=curve[-1][0],
                n_gp=4 * len(eles), n_at_floor=n_at, sum_Wf=sum_wf, n_events=events, failed_substeps=failed,
                aborted=refused > 0, wall=time.time() - t0)


def compare(a, b):
    """relative change of the final / peak load between two runs (|b - a| / a)."""
    return dict(final=abs(b["R_final"] - a["R_final"]) / abs(a["R_final"]), peak=abs(b["R_peak"] - a["R_peak"]) / abs(a["R_peak"]))


def report(runs):
    base = runs[1.0]
    lines = ["floor sensitivity, rigid strip footing on dry HAR sand (TIMs), displacement 0.1 m = 10 % B, half strip",
             f"{'p_min/0.505':>12s} {'p_min':>8s} {'R_final':>10s} {'R_peak':>10s} {'d_final%':>9s} {'d_peak%':>8s} {'GP at floor':>12s} {'sum W_f':>11s} "
             f"{'events':>7s} {'aborted':>8s} {'wall s':>7s}"]
    for s in sorted(runs, reverse=True):
        r = runs[s]
        c = compare(base, r)
        lines.append(f"{s:12.3g} {r['pmin']:8.4f} {r['R_final']:10.3f} {r['R_peak']:10.3f} {100 * c['final']:9.3f} {100 * c['peak']:8.3f} "
                     f"{r['n_at_floor']:5d}/{r['n_gp']:<5d} {r['sum_Wf']:11.3e} {r['n_events']:7d} {str(r['aborted']):>8s} {r['wall']:7.1f}")
    lines.append("acceptance (plan 2.7): the limit load moves by < 2 % when p_min is halved")
    return "\n".join(lines)


if __name__ == "__main__":
    if ops is None:
        raise SystemExit("opensees.pyd not found: build first (Ladruno_scripts\\build.bat) or set LADRUNO_OPENSEES_BIN")
    res = {s: run_strip(s) for s in (1.0, 0.5, 2.0)}
    print(report(res))
