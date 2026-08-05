"""The Chen-Han route, tested where it actually applies: PLANE-STRAIN strip
footings on the cone PDMY has, including the one problem with an EXACT answer.

`cone_probe.py`'s 52x estimate is a two-step argument: (a) the measured cone,
matched to Mohr-Coulomb IN PLANE STRAIN, is phi_ps = 53.7 deg; (b) Vesic's
square-footing formula at that angle gives 10.8 MPa of gamma-term. A
disagreement measured on the square footing cannot say which step failed,
because step (b) drags in a shape factor and an N_gamma extrapolated far past
its calibration range.

Plane strain separates them, and one plane-strain problem has an exact answer:

  * WEIGHTLESS soil (gamma = 0) under a surcharge q0 is the Prandtl-Reissner
    problem, q_u = q0*N_q, for which the limit-analysis upper and lower bounds
    COINCIDE. There is no shape factor, no N_gamma, no empirical fit — so
    numeric/q0*N_q(phi_ps) is a clean test of the plane-strain MATCHING itself.
  * The full strip (gamma' and q0 both on) then adds only N_gamma, so the two
    legs together separate "the matching is wrong" from "N_gamma is wrong".

WHY nu IS RAISED IN THESE LEGS. Under 1-D gravity/surcharge loading the elastic
K0 = nu/(1-nu) fixes the mobilisation of the initial state, m = sqrt(J2)/(3*a*p).
At PDMY's moduli nu = 0.337, K0 = 0.507 and m = 0.579 on the real cone — but
0.950 on a phi_txc = 20 deg cone, i.e. an easy-regime check at PDMY's own
moduli would start with the WHOLE DOMAIN at 95 % of yield and measure nothing
but that. (That is exactly what the first `verify_phi20` leg did, and why it is
reported as void rather than as a result.) A collapse load of an
elastic-perfectly-plastic body does not depend on its elastic constants, so nu
is a free numerical parameter here; `*_nuPDMY` legs measure that independence
instead of assuming it.

The mesh is generated here as a structured graded tensor grid (no gmsh): the
same x and z grading as `build_mesh_big.py` — 14.5 B of clearance, 10 B of
depth, 0.5 m elements under the footing — and ONE 0.5 m element in y with u_y
fixed on every node, which is plane strain exactly. 336 hexes, so a leg costs
minutes.

Run:  py -3.12 dp_strip.py [all|nq|full|<leg>]
"""
import os
import sys
import csv
import math
import time

import numpy as np

import dp_collapse as dc          # engine guard, cone conventions, Vesic
from dp_collapse import ops       # the SAME loaded module, never a second one

HERE = os.path.dirname(os.path.abspath(__file__))
B_FOOT = dc.B_FOOT
THICK = 0.5                       # out-of-plane slab thickness, m
SFRAC = 0.25                      # strips reach collapse deeper than squares
DS_BASE, DS_MIN, DS_MAX = 2.0e-4, 2.0e-6, 1.0e-3
GROW_AFTER = 6
K_EL = dc.K_EL
NU_SAFE = 0.45                    # see the module docstring


def g_from_nu(nu, k=K_EL):
    return 3.0 * k * (1.0 - 2.0 * nu) / (2.0 * (1.0 + nu))


def k0_mobilisation(alpha, nu):
    """m = sqrt(J2)/(3*alpha*p) of the 1-D elastic state — how close to yield
    the model already is before the footing has moved at all."""
    k0 = nu / (1.0 - nu)
    return (1.0 - k0) / (math.sqrt(3.0) * alpha * (1.0 + 2.0 * k0))


def graded(lo, hi, n, h0, fine_at_lo):
    """n intervals over [lo,hi] whose element at the FINE end is h0 and which
    grow geometrically outward — build_mesh.py's `growth`, inlined."""
    if abs(h0 * n - (hi - lo)) < 1e-12:
        return np.linspace(lo, hi, n + 1)
    a, b = 1.0, 4.0
    for _ in range(200):
        r = 0.5 * (a + b)
        tot = h0 * n if abs(r - 1) < 1e-12 else h0 * (r ** n - 1) / (r - 1)
        if tot < hi - lo:
            a = r
        else:
            b = r
    r = 0.5 * (a + b)
    e = np.concatenate([[0.0], np.cumsum(h0 * r ** np.arange(n))])
    e *= (hi - lo) / e[-1]
    return lo + e if fine_at_lo else hi - e[::-1]


XLIM, ZBOT = 30.0, -20.0
R_GRADE = 1.35          # target growth ratio of the far-field blocks


def n_graded(length, h0, r=R_GRADE):
    """Elements needed to span `length` starting at h0 and growing at ~r, so
    that refining h0 does NOT also shrink the domain or steepen the grading."""
    return max(1, int(np.ceil(np.log(1.0 + length * (r - 1.0) / h0)
                              / np.log(r))))


def strip_mesh(h0=0.5):
    nf = int(round(2.0 / h0))                          # elements across B
    no = n_graded(XLIM - 1.0, h0)
    nz = n_graded(-ZBOT - 3.0, h0)
    x = np.unique(np.concatenate([
        graded(-XLIM, -1.0, no, h0, False),            # fine at -1
        np.linspace(-1.0, 1.0, nf + 1),
        graded(1.0, XLIM, no, h0, True)]))             # fine at +1
    y = np.array([0.0, THICK])
    z = np.unique(np.concatenate([
        graded(ZBOT, -3.0, nz, h0, False),             # fine at -3
        np.linspace(-3.0, 0.0, int(round(3.0 / h0)) + 1)]))
    nx, ny, nz = len(x), len(y), len(z)
    idx = np.arange(nx * ny * nz).reshape(nx, ny, nz)
    nodes = np.stack(np.meshgrid(x, y, z, indexing="ij"), axis=-1).reshape(-1, 3)
    hexes = np.array([[idx[i, j, k], idx[i + 1, j, k],
                       idx[i + 1, j + 1, k], idx[i, j + 1, k],
                       idx[i, j, k + 1], idx[i + 1, j, k + 1],
                       idx[i + 1, j + 1, k + 1], idx[i, j + 1, k + 1]]
                      for i in range(nx - 1) for j in range(ny - 1)
                      for k in range(nz - 1)], dtype=np.int32)
    p = nodes[hexes]
    vol = np.einsum("ij,ij->i", p[:, 6] - p[:, 0],
                    np.cross(p[:, 1] - p[:, 0], p[:, 3] - p[:, 0]))
    assert (vol > 0).all(), "inverted hexes in the strip mesh"

    tol = 1e-9
    top = np.abs(nodes[:, 2]) < tol
    sets = dict(top=np.where(top)[0],
                bottom=np.where(np.abs(nodes[:, 2] - ZBOT) < tol)[0],
                xface=np.where(np.abs(np.abs(nodes[:, 0]) - XLIM) < tol)[0],
                footing=np.where(top & (np.abs(nodes[:, 0]) <= 1.0 + tol))[0])
    assert len(sets["footing"]) == 2 * (nf + 1), len(sets["footing"])

    trib = np.zeros(len(nodes))
    for h in hexes:
        face = h[[4, 5, 6, 7]]
        q = nodes[face]
        if np.any(np.abs(q[:, 2]) > tol):
            continue
        # shoelace on the diagonals; np.cross lost its 2-D form in numpy 2 and
        # this interpreter has numpy 2 (build_mesh.py runs under the apeGmsh
        # env, which still has numpy 1.x).
        d1, d2 = q[2, :2] - q[0, :2], q[3, :2] - q[1, :2]
        trib[face] += 0.5 * abs(d1[0] * d2[1] - d1[1] * d2[0]) / 4.0
    assert abs(trib.sum() - 2 * XLIM * THICK) < 1e-6, trib.sum()
    return nodes, hexes, trib, sets


def davis_phi(phi_deg, psi_deg=0.0):
    """Davis (1968) reduced friction angle for NON-ASSOCIATED flow:
    tan(phi*) = cos(psi) sin(phi) / (1 - sin(phi) sin(psi)); at psi = 0 that is
    simply tan(phi*) = sin(phi). It is the standard closed-form correction for
    exactly the assumption Chen & Han's matching makes, so it belongs beside
    the associated oracle rather than in a footnote."""
    p, s = math.radians(phi_deg), math.radians(psi_deg)
    return math.degrees(math.atan(math.cos(s) * math.sin(p)
                                  / (1.0 - math.sin(p) * math.sin(s))))


def vesic_strip(phi_deg, q0, gamma):
    """PLANE STRAIN — no shape factors. q0*N_q is Prandtl-Reissner and is EXACT
    for associated Mohr-Coulomb; the N_gamma term is the empirical part."""
    ph = math.radians(phi_deg)
    nq = math.exp(math.pi * math.tan(ph)) * math.tan(math.pi / 4 + ph / 2) ** 2
    ng = 2.0 * (nq + 1.0) * math.tan(ph)
    return dict(Nq=nq, Ngamma=ng, q_term=q0 * nq,
                g_term=0.5 * gamma * B_FOOT * ng,
                q_u=q0 * nq + 0.5 * gamma * B_FOOT * ng)


# --- legs -------------------------------------------------------------------
def _sleg(alpha, assoc, gamma, nu=NU_SAFE, q0=10.0, sy=dc.SY_BASE, h0=0.5):
    return dict(alpha=alpha, assoc=assoc, gamma=gamma, nu=nu, q0=q0, sy=sy,
                h0=h0)


A20 = dc.alpha_from_phi_txc(20.0)
NU_PDMY = 0.3366
SLEGS = {
    # Prandtl-Reissner: weightless, so q_u = q0*N_q EXACTLY (associated MC).
    # This is the test of the plane-strain MATCHING, with no empirical content.
    "nq_assoc":        _sleg(dc.ALPHA_PDMY, True, 0.0),
    "nq_nonassoc":     _sleg(dc.ALPHA_PDMY, False, 0.0),
    # same leg at PDMY's own nu: a collapse load must not depend on the elastic
    # constants, and this MEASURES that rather than asserting it.
    "nq_assoc_nuPDMY": _sleg(dc.ALPHA_PDMY, True, 0.0, nu=NU_PDMY),
    # machinery check in a regime Vesic/Prandtl coefficients were fitted near.
    # Only valid at the raised nu — at PDMY's nu this cone is 95 % mobilised by
    # the initial state alone.
    "nq20_assoc":      _sleg(A20, True, 0.0),
    "nq20_nonassoc":   _sleg(A20, False, 0.0),
    # the full strip: N_q and N_gamma together, against Vesic plane strain
    "full_assoc":      _sleg(dc.ALPHA_PDMY, True, dc.GAMMA_B),
    "full_nonassoc":   _sleg(dc.ALPHA_PDMY, False, dc.GAMMA_B),
}

# --- the mesh-refinement series --------------------------------------------
# Not one leg among the h0 = 0.5 legs above reaches a plateau, and the easy
# regime — where the answer is EXACT — walls out 12 % ABOVE it. That is the
# signature of a discretisation that cannot form the mechanism: 4 elements
# across B, and a phi_ps = 54 deg failure surface is a thin fan. Refinement is
# the experiment that separates "this soil punches" from "this mesh cannot
# make a band", and in plane strain it is affordable (336 -> ~1700 hexes).
# `nq20` is the control: its oracle is exact, so its convergence with h0 says
# how much of any residual gap is discretisation. The cone legs run at PDMY's
# own nu (m0 = 0.579, still safe) because that E is 3x stiffer and so reaches
# collapse in a third of the settlement.
for _h in (0.25, 0.125):
    _t = f"_h{str(_h).replace('0.', '')}"
    SLEGS[f"nq20_assoc{_t}"] = _sleg(A20, True, 0.0, h0=_h)
    SLEGS[f"nq_assoc{_t}"] = _sleg(dc.ALPHA_PDMY, True, 0.0, nu=NU_PDMY, h0=_h)
    SLEGS[f"full_nonassoc{_t}"] = _sleg(dc.ALPHA_PDMY, False, dc.GAMMA_B,
                                        nu=NU_PDMY, h0=_h)
    # NON-associated is both the physical match to PDMY (dil1 = dil2 = 0) and
    # the numerically well-behaved rung: nq20_nonassoc holds the exact Prandtl
    # answer to 0.2 % from s/B = 0.09 to 0.30, while every ASSOCIATED leg
    # hardens without ever plateauing — dilatancy at psi = phi has to push the
    # surrounding soil out of the way, and a bounded mesh resists that. So the
    # cone's Prandtl problem gets a non-associated series too.
    SLEGS[f"nq_nonassoc{_t}"] = _sleg(dc.ALPHA_PDMY, False, 0.0, nu=NU_PDMY,
                                      h0=_h)
    SLEGS[f"nq20_nonassoc{_t}"] = _sleg(A20, False, 0.0, h0=_h)
SLEGS["nq_nonassoc_h5"] = _sleg(dc.ALPHA_PDMY, False, 0.0, nu=NU_PDMY)
SLEGS["nq20_assoc_h5"] = _sleg(A20, True, 0.0)
SLEGS["nq_assoc_h5"] = _sleg(dc.ALPHA_PDMY, True, 0.0, nu=NU_PDMY)
SLEGS["full_nonassoc_h5"] = _sleg(dc.ALPHA_PDMY, False, dc.GAMMA_B, nu=NU_PDMY)

SGROUPS = {"all": list(SLEGS),
           "nq": ["nq_assoc", "nq_nonassoc", "nq_assoc_nuPDMY"],
           "easy": ["nq20_assoc", "nq20_nonassoc"],
           "full": ["full_assoc", "full_nonassoc"],
           # the refinement series, one group per family so they parallelise
           "ref20": ["nq20_assoc_h5", "nq20_assoc_h25", "nq20_assoc_h125"],
           "refnq": ["nq_assoc_h5", "nq_assoc_h25", "nq_assoc_h125"],
           "reffull": ["full_nonassoc_h5", "full_nonassoc_h25",
                       "full_nonassoc_h125"],
           "refnqna": ["nq_nonassoc_h5", "nq_nonassoc_h25",
                       "nq_nonassoc_h125"],
           "ref20na": ["nq20_nonassoc_h25", "nq20_nonassoc_h125"]}


# --- an algorithm ladder ----------------------------------------------------
# A perfectly plastic collapse is where Newton is worst: the tangent goes
# singular in the mechanism mode exactly when the answer arrives. Halving the
# step alone is not enough (measured: a strip leg walled out at 58 % of its
# initial tangent, nowhere near a plateau), so each increment gets a ladder of
# attempts before the step is shrunk, and any step that needed the RELAXED
# tolerance is flagged in the CSV so it cannot pass unnoticed.
def attempts(tol):
    return [("KrylovNewton", tol, 25, 0),
            ("NewtonLineSearch", tol, 40, 0),
            ("KrylovNewton", 10.0 * tol, 60, 1)]


def set_attempt(a):
    algo, tol, it, _ = a
    ops.test("NormUnbalance", tol, it, 0)
    if algo == "NewtonLineSearch":
        ops.algorithm("NewtonLineSearch")
    else:
        ops.algorithm("KrylovNewton")


def run(name):
    cfg = SLEGS[name]
    alpha, assoc, gamma, nu, q0, sy = (cfg["alpha"], cfg["assoc"], cfg["gamma"],
                                       cfg["nu"], cfg["q0"], cfg["sy"])
    nodes, hexes, trib, sets = strip_mesh(cfg["h0"])
    rho = math.sqrt(2.0) * alpha
    rho_bar = rho if assoc else 0.0
    g_el = g_from_nu(nu)
    mc = dc.mc_from_cone(alpha)
    vs = vesic_strip(mc["ps"], q0, gamma)
    m0 = k0_mobilisation(alpha, nu)
    dc.log(name, f"strip h0={cfg['h0']} m ({int(round(2.0 / cfg['h0']))} "
           f"elements across B), {len(nodes)} nodes / {len(hexes)} hex; "
           f"alpha={alpha:.5f} "
           f"(phi_txc={mc['txc']:.2f}, phi_ps={mc['ps']:.2f}), "
           f"{'associated' if assoc else 'NON-associated'}, "
           f"gamma'={gamma:.2f}, q0={q0}, nu={nu} (G={g_el:.4g}), sigma_y={sy}")
    dc.log(name, f"initial 1-D state is at m = {m0:.3f} of yield "
           + ("(SAFE)" if m0 < 0.8 else "(TOO CLOSE TO YIELD — leg is void)"))
    dc.log(name, f"oracle: Vesic/Prandtl plane strain q_u = {vs['q_u']:.1f} kPa "
           f"= q0*Nq {vs['q_term']:.1f} + gamma-term {vs['g_term']:.1f} "
           f"(N_q={vs['Nq']:.1f}, N_gamma={vs['Ngamma']:.1f})")
    if not assoc:
        pd = davis_phi(mc["ps"])
        vd = vesic_strip(pd, q0, gamma)
        dc.log(name, f"  Davis non-associated reduction: phi* = {pd:.2f} deg "
               f"=> q_u = {vd['q_u']:.1f} kPa (N_q={vd['Nq']:.1f}), "
               f"{vs['q_u'] / vd['q_u']:.1f}x below the associated oracle")

    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 3)
    for i, (x, y, z) in enumerate(nodes, start=1):
        ops.node(i, float(x), float(y), float(z))
    ops.nDMaterial("DruckerPrager", 1, K_EL, g_el, sy, rho, rho_bar,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for e, conn in enumerate(hexes, start=1):
        ops.element("LadrunoBrick", e, *[int(c) + 1 for c in conn], 1,
                    "-geom", "linear", "-b", 0.0, 0.0, -gamma,
                    "-formulation", "bbar")
    # one fix() per node — a second fix() on the same node adds a DUPLICATE
    # SP_Constraint. u_y = 0 on every node is what makes this plane strain.
    xf, bt = set(sets["xface"].tolist()), set(sets["bottom"].tolist())
    for n in range(len(nodes)):
        if n in bt:
            ops.fix(n + 1, 1, 1, 1)
        else:
            ops.fix(n + 1, 1 if n in xf else 0, 1, 0)

    ops.timeSeries("Constant", 1)
    ops.pattern("Plain", 1, 1)
    for n in sets["top"]:
        if trib[n] > 0:
            ops.load(int(n) + 1, 0.0, 0.0, -q0 * float(trib[n]))

    vol = 2 * XLIM * THICK * (-ZBOT)
    tol = 1.0e-5 * max(gamma * vol + q0 * float(trib.sum()), 1.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    dc.pick_system()
    ops.test("NormUnbalance", tol, 25, 0)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    assert ops.analyze(1) == 0, f"{name} gravity failed"
    ops.reactions()
    rz = sum(ops.nodeReaction(int(n) + 1, 3) for n in sets["bottom"])
    want = gamma * vol + q0 * float(trib.sum())
    dc.log(name, f"initial equilibrium: {rz:.4f} vs {want:.4f} kN "
           f"({100 * (rz / want - 1):+.5f} %)")
    assert abs(rz / want - 1) < 1e-6, "body force / tributary convention"

    foot = [int(n) + 1 for n in sets["footing"]]
    r_corr = q0 * float(trib[sets["footing"]].sum())
    area = B_FOOT * THICK
    uz0 = ops.nodeDisp(foot[0], 3)
    ops.loadConst("-time", 0.0)
    smax = SFRAC * B_FOOT
    ops.timeSeries("Path", 2, "-time", 0.0, 1.0,
                   "-values", uz0, uz0 - smax, "-useLast")
    ops.pattern("Plain", 2, 2)
    for t in foot:
        ops.sp(t, 3, 1.0)
    ops.wipeAnalysis()
    ops.constraints("Transformation")
    ops.numberer("RCM")
    dc.pick_system()
    ops.analysis("Static")

    path = os.path.join(HERE, f"dpstrip_{name}.csv")
    f = open(path, "w", newline="")
    w = csv.writer(f)
    w.writerow(["s_m", "R_kN", "q_kPa", "ds_mm", "relaxed", "wall_s"])
    rows, ds, good, nfail, nrelax, t0 = [], DS_BASE, 0, 0, 0, time.time()
    lad = attempts(tol)
    while True:
        s_now = ops.getTime() * smax
        if s_now >= smax - 1e-12:
            break
        ds = min(ds, smax - s_now)
        ops.integrator("LoadControl", ds / smax)
        ok, relaxed = False, 0
        for a in lad:
            set_attempt(a)
            if ops.analyze(1) == 0:
                ok, relaxed = True, a[3]
                break
            nfail += 1
        if not ok:
            good = 0
            ds *= 0.5
            if ds < DS_MIN:
                dc.log(name, f"no convergence at s/B={s_now / B_FOOT:.4f} "
                       f"(every ladder rung failed at ds={ds * 1000:.4f} mm)")
                break
            continue
        nrelax += relaxed
        good += 1
        if good >= GROW_AFTER and ds < DS_MAX:
            ds, good = min(2 * ds, DS_MAX), 0
        ops.reactions()
        R = -sum(ops.nodeReaction(t, 3) for t in foot) + r_corr
        s = ops.getTime() * smax
        rows.append((s, R, R / area, ds * 1000, relaxed, time.time() - t0))
        w.writerow([f"{v:.9g}" for v in rows[-1]])
        f.flush()
    f.close()

    a = np.array(rows)
    s, q = a[:, 0], a[:, 2]
    m = s >= 0.9 * s[-1]
    t_last = np.polyfit(s[m], q[m], 1)[0]
    n0 = max(4, len(s) // 50)
    t_init = np.polyfit(s[:n0], q[:n0], 1)[0]
    plateau = abs(t_last) < 0.02 * abs(t_init)
    dc.log(name, f"q_max = {q.max():.1f} kPa at s/B = "
           f"{s[int(q.argmax())] / B_FOOT:.4f}; end s/B = {s[-1] / B_FOOT:.4f}; "
           f"{len(rows)} steps, {nfail} failed attempts, {nrelax} relaxed")
    dc.log(name, f"tail dq/ds = {t_last:.0f} kPa/m = "
           f"{100 * t_last / t_init:.2f} % of initial => "
           f"{'PLATEAU (collapse reached)' if plateau else 'STILL HARDENING'}")
    dc.log(name, f"numeric / oracle = {q.max() / vs['q_u']:.4f}   "
           f"[{q.max():.1f} vs {vs['q_u']:.1f} kPa]")
    return dict(name=name, qmax=float(q.max()), oracle=vs, plateau=plateau,
                m0=m0, mc=mc, cfg=cfg)


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    names = SGROUPS.get(which, [which])
    for n in names:
        if n not in SLEGS:
            raise SystemExit(f"unknown strip leg {n!r}; pick from "
                             f"{list(SLEGS)} or a group {list(SGROUPS)}")
    out = []
    for n in names:
        dc.log("=== strip leg:", n)
        out.append(run(n))
    print()
    hdr = (f"{'leg':>18} {'flow':>6} {'m0':>6} {'q_num':>10} {'oracle':>10} "
           f"{'num/oracle':>11} {'plateau':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in out:
        print(f"{r['name']:>18} "
              f"{'assoc' if r['cfg']['assoc'] else 'non':>6} "
              f"{r['m0']:6.3f} {r['qmax']:10.1f} {r['oracle']['q_u']:10.1f} "
              f"{r['qmax'] / r['oracle']['q_u']:11.4f} "
              f"{'yes' if r['plateau'] else 'NO':>8}")


if __name__ == "__main__":
    main()
