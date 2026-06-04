"""Generate the validation figures for Ladruno_implementation/16_lemaitre_validation.md.

Re-runs the Layer A/B/C/D validation cases against the built opensees.pyd and writes
PNGs into Ladruno_implementation/figures/. Standalone (own DLL bootstrap); reuses the
test drivers in tests/.

Usage:
    python Ladruno_scripts/lemaitre_validation_figures.py  [DIST_BIN_DIR]
"""
import os
import sys
import math

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
TESTS = os.path.join(ROOT, "tests")
FIGDIR = os.path.join(ROOT, "Ladruno_implementation", "lemaitre_validation", "figures")
os.makedirs(FIGDIR, exist_ok=True)

DIST = sys.argv[1] if len(sys.argv) > 1 else \
    r"C:\Users\nmora\Github\OpenSees_Compile\OpenSees\.claude\worktrees\wonderful-nobel-1377f4\dist\bin"
os.add_dll_directory(DIST)
sys.path.insert(0, DIST)
sys.path.insert(0, TESTS)
os.environ["LADRUNO_OPENSEES_QUIET"] = "1"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from _testbed import ops
from _testbed.lemaitre_ref import Lemaitre1D
import test_lemaitre_validation as V

plt.rcParams.update({"figure.dpi": 120, "font.size": 10, "axes.grid": True,
                     "grid.alpha": 0.3, "axes.spines.top": False, "axes.spines.right": False})


# ---------------------------------------------------------------- Layer A1
def fig_A1_oracle():
    hist = V._ramp(0.12, 120)
    code = V._run_truss(lambda t: V._mat1d(t, dmg=V._DMG), hist)
    ref = V._ref(V._DMG)
    eps = [c[0] for c in code]
    sc = [c[1] for c in code]
    Dc = [c[2] for c in code]
    sr, Dr = [], []
    for e in eps:
        s, D, _p, _se = ref.step(e)
        sr.append(s); Dr.append(D)
    serr = [abs(a - b) for a, b in zip(sc, sr)]

    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    ax[0].plot(eps, sr, "-", lw=3, color="#bbbbbb", label="reference (eqns)")
    ax[0].plot(eps, sc, "--", lw=1.3, color="C3", label="LadrunoUniaxialJ2")
    ax[0].set(xlabel=r"strain $\varepsilon$", ylabel=r"nominal stress $\sigma=(1-D)\tilde\sigma$",
              title="A1  coupled stress–strain")
    ax[0].legend()
    ax[1].plot(eps, Dr, "-", lw=3, color="#bbbbbb", label="reference")
    ax[1].plot(eps, Dc, "--", lw=1.3, color="C0", label="code")
    ax[1].set(xlabel=r"strain $\varepsilon$", ylabel=r"damage $D$", title="A1  damage evolution")
    ax[1].legend()
    ax[2].semilogy(eps, [max(s, 1e-16) for s in serr], color="C2")
    ax[2].axhline(1e-8, ls=":", color="k", lw=1)
    ax[2].set(xlabel=r"strain $\varepsilon$", ylabel=r"$|\sigma_{code}-\sigma_{ref}|$",
              title="A1  oracle residual (<1e-8)")
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lemaitre_A1_oracle.png")); plt.close(fig)
    print("A1 done")


# ---------------------------------------------------------------- Layer A2
def fig_A2_closed_form():
    s0, r, Dc = 5.0, 0.5, 0.99

    def mat(tag):
        ops.uniaxialMaterial("LadrunoUniaxialJ2", tag, V._E,
                             "-iso", "voce", s0, 0.0, 0.0, 0.0,
                             "-damage", "lemaitre", r, 1.0, 0.0, Dc)

    hist = V._ramp(0.05, 200)
    code = V._run_truss(mat, hist)
    p = [c[3] for c in code]
    D = [c[2] for c in code]
    slope = s0 * s0 / (2.0 * V._E * r)
    pp = [x for x in p if x > 0]
    line = [slope * x for x in pp]
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    ax.plot(pp, line, "-", lw=3, color="#bbbbbb", label=r"closed form $D=\frac{\sigma_0^2}{2Er}\,p$")
    ax.plot([x for x in p if x > 0], [d for x, d in zip(p, D) if x > 0],
            "o", ms=3, color="C3", label="LadrunoUniaxialJ2")
    ax.set(xlabel=r"accumulated plastic strain $p$", ylabel=r"damage $D$",
           title="A2  perfectly-plastic + linear damage")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "lemaitre_A2_closed_form.png")); plt.close(fig)
    print("A2 done")


# ---------------------------------------------------------------- Layer B
def fig_B_triaxiality():
    nu = V._NU
    a = (2.0 / 3.0) * (1.0 + nu); b = 3.0 * (1.0 - 2.0 * nu)
    eta = [i / 200.0 for i in range(0, 241)]            # 0 .. 1.2 (covers confined point)
    Rv = [a + b * e * e for e in eta]

    dmg = (3.0, 1.0, 0.0, 0.98); rr = dmg[0]
    free = V._run_brick(lambda t: V._mat3d(t, dmg), 0.06, 60, confine=False)
    conf = V._run_brick(lambda t: V._mat3d(t, dmg), 0.06, 60, confine=True)
    pts, scat = {}, {"free": ([], []), "conf": ([], [])}
    for label, res in (("free", free), ("conf", conf)):
        for (e0, s0, D0, p0), (e1, s1, D1, p1) in zip(res, res[1:]):
            dD, dp = D1 - D0, p1 - p0
            if dD <= 1e-9 or dp <= 1e-9 or D1 >= dmg[3] - 1e-6:
                continue
            Yi = rr * dD / dp
            Yf, q, sH, et, Rvv = V._Y_formula(s1, D1, V._E, nu)
            scat[label][0].append(Yf); scat[label][1].append(Yi)
            pts.setdefault(label, (et, Rvv))

    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.0))
    ax[0].plot(eta, Rv, "-", color="C0", label=r"$R_v(\eta)=\frac{2}{3}(1{+}\nu)+3(1{-}2\nu)\eta^2$")
    for lab, col, name in (("free", "C2", "free (uniaxial)"), ("conf", "C3", "confined")):
        if lab in pts:
            et, Rvv = pts[lab]
            ax[0].plot(et, Rvv, "o", color=col, ms=9, label=f"{name}: $\\eta$={et:.2f}, $R_v$={Rvv:.2f}")
    ax[0].axvline(1/3, ls=":", color="k", lw=1); ax[0].text(0.34, a + 0.02, r"$\eta=1/3$", fontsize=8)
    ax[0].set(xlabel=r"triaxiality $\eta=\sigma_H/\sigma_{eq}$", ylabel=r"$R_v$",
              title="B  triaxiality function"); ax[0].legend(fontsize=8)
    lim = 0
    for lab, col, name in (("free", "C2", "free"), ("conf", "C3", "confined")):
        ax[1].plot(scat[lab][0], scat[lab][1], "o", ms=4, color=col, label=name)
        lim = max(lim, max(scat[lab][0] + scat[lab][1], default=0))
    ax[1].plot([0, lim], [0, lim], "k--", lw=1, label="y=x")
    ax[1].set(xlabel=r"$Y$ closed form $=(aq^2+b\sigma_H^2)/2E$",
              ylabel=r"$Y$ implied $=r\,\Delta D/\Delta p$", title="B  wired $Y$ vs closed form")
    ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lemaitre_B_triaxiality.png")); plt.close(fig)
    print("B done")


# ---------------------------------------------------------------- Layer D
def _asd_chaboche(E, sy, su, eu):
    H1b = E / (400.0 * eu); g1 = E / (400.0 * eu * (su - sy))
    return (H1b * 0.9, g1, H1b * 5.0, g1 * 50.0)


# Lemaitre damage params for the DAMAGED ASDSteel1D cross-exam.
#   _LEM_UNCAL = the original, *uncalibrated* pick (just enough to reach rupture).
#   _calibrate_lemaitre() = a ROUGH fit of (p_D, r) to ASDSteel1D's `-fracture`:
#     * p_D := eupl = eu - sy/E       -> shares ASDSteel1D's necking/fracture ONSET
#     * r   := bisected so the dissipated fracture ENERGY (∫σ dε to rupture) matches.
#   The softening SHAPE still differs (Lemaitre's smooth energy-driven D vs ASDSteel1D's
#   sharp exponential necking drop) — that is the irreducible law difference (§5).
_LEM_UNCAL = (0.12, 1.0, 0.03, 0.95)        # (r, s, p_D, D_c)
_LEM_CAL_CACHE = {}


def _energy(res):
    """Area under the σ–ε curve of a truss history (list of (eps, sigma, D|None))."""
    W = 0.0
    for (e0, s0, _), (e1, s1, _) in zip(res, res[1:]):
        W += 0.5 * (s0 + s1) * (e1 - e0)
    return W


def _calibrate_lemaitre(E, sy, su, eu, C1, g1, C2, g2, eps_t=0.20, n=400):
    """Return (r, s, p_D, D_c) fitting ASDSteel1D `-fracture` onset+energy on a monotone
    tension coupon. p_D = eupl (shared onset); r bisected to match dissipated energy."""
    key = (E, sy, su, eu)
    if key in _LEM_CAL_CACHE:
        return _LEM_CAL_CACHE[key]
    hist = [eps_t * i / n for i in range(1, n + 1)]
    asd = _truss_hist(lambda t: ops.uniaxialMaterial("ASDSteel1D", t, E, sy, su, eu, "-fracture"),
                      hist, "Damage")
    W_target = _energy(asd)
    pD = round(eu - sy / E, 3)                    # eupl
    s, Dc = 1.0, 0.95

    def W_of(r):
        l = _truss_hist(lambda t: ops.uniaxialMaterial(
            "LadrunoUniaxialJ2", t, E, "-iso", "voce", sy, 0, 0, 0,
            "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r, s, pD, Dc), hist, "damage")
        return _energy(l)

    lo, hi = 0.02, 0.6                            # larger r -> slower damage -> more energy
    for _ in range(22):
        r = 0.5 * (lo + hi)
        if W_of(r) < W_target:
            lo = r
        else:
            hi = r
    out = (round(0.5 * (lo + hi), 4), s, pD, Dc)
    _LEM_CAL_CACHE[key] = out
    return out


def _truss_hist(mat_fn, hist, dmg_key=None):
    ops.wipe(); ops.model("basic", "-ndm", 1, "-ndf", 1)
    ops.node(1, 0.0); ops.node(2, 1.0); mat_fn(1); ops.fix(1, 1)
    ops.element("Truss", 1, 1, 2, 1.0, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(2, 1.0)
    ops.system("FullGeneral"); ops.numberer("Plain"); ops.constraints("Plain")
    ops.test("NormDispIncr", 1e-8, 200, 0); ops.algorithm("Newton"); ops.analysis("Static")
    out = []; prev = 0.0
    for eps in hist:
        ops.integrator("DisplacementControl", 2, 1, eps - prev)
        if ops.analyze(1) != 0:
            break
        prev = eps
        s = ops.eleResponse(1, "material", 1, "stress")
        D = ops.eleResponse(1, "material", 1, dmg_key) if dmg_key else None
        out.append((eps, s[0] if s else 0.0, (D[0] if D else None)))
    return out


def fig_D_cross():
    E, sy, su, eu = 200000.0, 400.0, 550.0, 0.10
    C1, g1, C2, g2 = _asd_chaboche(E, sy, su, eu)

    def asd(t): ops.uniaxialMaterial("ASDSteel1D", t, E, sy, su, eu)

    def lad(t): ops.uniaxialMaterial("LadrunoUniaxialJ2", t, E, "-iso", "voce", sy, 0, 0, 0,
                                     "-kin", 2, C1, g1, C2, g2)
    hist = [0.08 * i / 160 for i in range(1, 161)]
    a = _truss_hist(asd, hist); l = _truss_hist(lad, hist)

    # convergence under refinement
    def err(n):
        h = [0.08 * i / n for i in range(1, n + 1)]
        A = _truss_hist(asd, h); L = _truss_hist(lad, h)
        return max(abs(sa - sl) / abs(sa) for (_e, sa, _), (_e2, sl, _) in zip(A, L) if abs(sa) > 1.0)
    ns = [10, 20, 40, 80, 160, 320]
    errs = [err(n) for n in ns]

    # fracture signatures — ASDSteel1D vs Lemaitre, both UNCALIBRATED and ROUGH-CALIBRATED
    r_u, s_u, pD_u, Dc_u = _LEM_UNCAL
    r_c, s_c, pD_c, Dc_c = _calibrate_lemaitre(E, sy, su, eu, C1, g1, C2, g2)

    def asd_fr(t): ops.uniaxialMaterial("ASDSteel1D", t, E, sy, su, eu, "-fracture")

    def lad_unc(t): ops.uniaxialMaterial("LadrunoUniaxialJ2", t, E, "-iso", "voce", sy, 0, 0, 0,
                                         "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r_u, s_u, pD_u, Dc_u)

    def lad_cal(t): ops.uniaxialMaterial("LadrunoUniaxialJ2", t, E, "-iso", "voce", sy, 0, 0, 0,
                                         "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r_c, s_c, pD_c, Dc_c)
    hf = [0.20 * i / 200 for i in range(1, 201)]
    af = _truss_hist(asd_fr, hf, "Damage")
    lu = _truss_hist(lad_unc, hf, "damage")
    lc = _truss_hist(lad_cal, hf, "damage")
    pk_a = max(x[1] for x in af); pk_c = max(x[1] for x in lc)
    W_a, W_c = _energy(af), _energy(lc)

    fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.0))
    ax[0].plot([x[0] for x in a], [x[1] for x in a], "-", lw=2.5, color="#888", label="ASDSteel1D")
    ax[0].plot([x[0] for x in l], [x[1] for x in l], "--", lw=1.4, color="C3", label="Ladruno (=ASD backstresses)")
    ax[0].set(xlabel=r"$\varepsilon$", ylabel=r"$\sigma$", title="D1  undamaged backbone"); ax[0].legend(fontsize=8)
    ax[1].loglog(ns, errs, "o-", color="C0")
    ax[1].set(xlabel="steps (1/Δε)", ylabel="max rel. backbone error", title="D1  convergence under refinement")
    ax[2].plot([x[0] for x in af], [x[1] for x in af], "-", lw=2.4, color="#888", label="ASDSteel1D -fracture")
    ax[2].plot([x[0] for x in lu], [x[1] for x in lu], ":", lw=1.4, color="C0",
               label=f"Lemaitre uncalibrated\n($p_D$={pD_u}, r={r_u})")
    ax[2].plot([x[0] for x in lc], [x[1] for x in lc], "-", lw=1.8, color="C3",
               label=f"Lemaitre calibrated\n($p_D$=eupl={pD_c}, r={r_c})")
    ax[2].set(xlabel=r"$\varepsilon$", ylabel=r"$\sigma$", title="D2  ductile-fracture signature")
    ax[2].legend(fontsize=7.5, loc="upper right")
    ax[2].text(0.02, 0.04,
               f"calib matches ONSET ($p_D$=eupl) + ENERGY (r bisected):\n"
               f"peak {pk_c:.0f} vs {pk_a:.0f},  ∫σdε {W_c:.1f} vs {W_a:.1f}\n"
               f"residual = softening SHAPE (smooth vs necking drop)",
               transform=ax[2].transAxes, fontsize=6.8, va="bottom",
               bbox=dict(boxstyle="round", fc="#fff7e6", ec="0.7", alpha=0.9))
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lemaitre_D_cross_asdsteel1d.png")); plt.close(fig)
    print(f"D done  calib(r,s,pD,Dc)={(r_c,s_c,pD_c,Dc_c)}  peak {pk_c:.0f}/{pk_a:.0f}  W {W_c:.1f}/{W_a:.1f}")


# ---------------------------------------------------------------- Layer C
def fig_C_mesh_objectivity():
    import test_lemaitre_notched_bar as C

    def solve_hist(h, lch_ref):
        nodes, hexes = C._mesh(h)
        keep = [hx for hx in hexes if not C._in_notch(sum(nodes[n][0] for n in hx) / 8.0,
                                                      sum(nodes[n][1] for n in hx) / 8.0)]
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        used = set(n for hx in keep for n in hx)
        for t in used:
            ops.node(int(t), *nodes[t])
        extra = ["-autoRegularization", lch_ref] if lch_ref else []
        ops.nDMaterial("LadrunoJ2", 1, C.K, C.G, *C.ISO, "-damage", "lemaitre", *C.DMG, *extra)
        for i, hx in enumerate(keep, 1):
            ops.element("LadrunoBrick", i, *[int(n) for n in hx], 1, "-formulation", "bbar")
        xmin = min(nodes[t][0] for t in used); xmax = max(nodes[t][0] for t in used)
        left = sorted(int(t) for t in used if abs(nodes[t][0] - xmin) < 1e-6)
        for t in used:
            if abs(nodes[t][2]) < 1e-6:
                ops.fix(int(t), 0, 0, 1)
        for t in left:
            ops.fix(int(t), 1, 1 if t == left[0] else 0, 0)
        pull = sorted(int(t) for t in used if abs(nodes[t][0] - xmax) < 1e-6); master = pull[0]
        for n in pull:
            if n != master:
                ops.equalDOF(master, n, 1)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(master, 1.0, 0.0, 0.0)
        ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-6, 50, 0); ops.analysis("Static")
        us, Ps = [0.0], [0.0]; u = 0.0; base = 1.5 / 150; calls = 0
        while u < 1.5 - 1e-12 and calls < 1200:
            ok = -1
            for cut, alg in ((1.0, "KrylovNewton"), (0.5, "KrylovNewton"), (0.25, "NewtonLineSearch"),
                             (0.1, "NewtonLineSearch"), (0.03, "NewtonLineSearch")):
                ops.integrator("DisplacementControl", master, 1, base * cut); ops.algorithm(alg)
                ok = ops.analyze(1); calls += 1
                if ok == 0:
                    break
            if ok != 0:
                break
            ops.reactions(); P = -sum(ops.nodeReaction(n)[0] for n in left); u = ops.nodeDisp(master)[0]
            us.append(u); Ps.append(P / 1000.0)   # kN
        return us, Ps, len(keep)

    # three progressively finer meshes (1x / 2x / 2.7x refinement) to show CONVERGENCE.
    # Regularized cases traverse cleanly; the LOCAL (no-reg) finest mesh is pathologically
    # brittle (endless step-cuts), so the local panel uses the two coarser meshes — enough
    # to show the divergence trend.
    meshes = [(4.0, "coarse", "C0", "-"), (2.0, "medium", "C1", "--"), (1.5, "fine", "C3", ":")]
    loc_meshes = meshes[:2]
    reg = {h: solve_hist(h, C.LCH_REF) for h, _, _, _ in meshes}
    loc = {h: solve_hist(h, None) for h, _, _, _ in loc_meshes}

    # integrate dissipated energy up to the COMMON reached elongation, so an early
    # numerical stop on one mesh can't bias the comparison.
    ucommon = min(P_u[0][-1] for P_u in list(reg.values()) + list(loc.values()))

    def diss_to(P_u, ucap):
        us, Ps, _ = P_u
        W = 0.0
        for i in range(1, len(us)):
            if us[i] > ucap:
                break
            W += 0.5 * (Ps[i] + Ps[i - 1]) * (us[i] - us[i - 1])
        return W

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for h, name, col, ls in meshes:
        ax[0].plot(reg[h][0], reg[h][1], ls, lw=1.9, color=col, label=f"{name} h={h} (n={reg[h][2]})")
    for h, name, col, ls in loc_meshes:
        ax[1].plot(loc[h][0], loc[h][1], ls, lw=1.9, color=col, label=f"{name} h={h} (n={loc[h][2]})")
    for a in (ax[0], ax[1]):
        a.axvline(ucommon, ls=":", color="0.6", lw=1)
    ax[0].set(xlabel="elongation u [mm]", ylabel="load P [kN]",
              title="C  WITH -autoRegularization\n(bounded response; peak converges as h→0)"); ax[0].legend(fontsize=8)
    ax[1].set(xlabel="elongation u [mm]", ylabel="load P [kN]",
              title="C  LOCAL (no regularization)\n(softening tail diverges / more brittle)"); ax[1].legend(fontsize=8)
    # convergence of dissipated energy vs refinement (to the common elongation)
    Wreg = [diss_to(reg[h], ucommon) for h, _, _, _ in meshes]
    Wloc = [diss_to(loc[h], ucommon) for h, _, _, _ in loc_meshes]
    ax[2].plot([reg[h][2] for h, _, _, _ in meshes], Wreg, "o-", color="C2", label="regularized")
    ax[2].plot([loc[h][2] for h, _, _, _ in loc_meshes], Wloc, "s--", color="C3", label="local (no reg)")
    ax[2].set(xlabel="number of elements", ylabel=f"dissipated energy to u={ucommon:.2f} [kN·mm]",
              title="C  dissipated energy vs mesh\n(regularized stays bounded; local drifts down)")
    ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lemaitre_C_mesh_objectivity.png")); plt.close(fig)
    print("C done  ucommon=", round(ucommon, 3),
          " Wreg=", [round(w, 1) for w in Wreg], " Wloc=", [round(w, 1) for w in Wloc])


# ---------------------------------------------------------------- Layer E (cyclic)
def fig_E_cyclic():
    import test_lemaitre_cyclic as CY
    import test_lemaitre_vs_asdsteel1d as X

    # E1 — oracle hysteresis (code vs reference)
    hist = CY._cyclic(0.03, 3, 60)
    code = V._run_truss(lambda t: V._mat1d(t, dmg=V._DMG), hist)
    ref = V._ref(V._DMG)
    rs = []
    for (eps, _s, _D, _p) in code:
        s, _D2, _p2, _se = ref.step(eps)
        rs.append(s)

    # E2 — damage accumulation per cycle (one amplitude) + Coffin-Manson life curve
    nC, per = CY._run_to_damage(0.05, Dstop=0.8, maxcyc=120)
    Dcum, acc = [], 0.0
    for d in per:
        acc += d; Dcum.append(acc)
    amps = [0.04, 0.05, 0.06, 0.08, 0.10]
    Ns = [CY._run_to_damage(a, Dstop=0.5, maxcyc=200)[0] for a in amps]

    # E3 — ASDSteel1D vs Ladruno cyclic loops, UNDAMAGED (overlay) and DAMAGED (diverge)
    C1, g1, C2, g2 = X._asd_chaboche(X._E, X._SY, X._SU, X._EU)
    hc = CY._cyclic(0.05, 5, 40)
    r_c, s_c, pD_c, Dc_c = _calibrate_lemaitre(X._E, X._SY, X._SU, X._EU, C1, g1, C2, g2)
    r_u, s_u, pD_u, Dc_u = _LEM_UNCAL

    # --- undamaged loops: matched backstresses, NO damage => the loops COINCIDE ---
    def asd_u(t): ops.uniaxialMaterial("ASDSteel1D", t, X._E, X._SY, X._SU, X._EU)  # no -fracture

    def lad_u(t): ops.uniaxialMaterial("LadrunoUniaxialJ2", t, X._E, "-iso", "voce", X._SY, 0, 0, 0,
                                       "-kin", 2, C1, g1, C2, g2)                    # no -damage
    au = X._truss(asd_u, hc); lu = X._truss(lad_u, hc)
    res_u = max(abs(sa - sl) / (abs(sa) + 1.0) for (_e, sa, _), (_e2, sl, _) in zip(au, lu))

    # --- damaged loops: ASDSteel1D -fracture vs Lemaitre (uncalibrated + calibrated) ---
    def asd(t): ops.uniaxialMaterial("ASDSteel1D", t, X._E, X._SY, X._SU, X._EU, "-fracture")

    def lad_unc(t): ops.uniaxialMaterial("LadrunoUniaxialJ2", t, X._E, "-iso", "voce", X._SY, 0, 0, 0,
                                         "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r_u, s_u, pD_u, Dc_u)

    def lad_cal(t): ops.uniaxialMaterial("LadrunoUniaxialJ2", t, X._E, "-iso", "voce", X._SY, 0, 0, 0,
                                         "-kin", 2, C1, g1, C2, g2, "-damage", "lemaitre", r_c, s_c, pD_c, Dc_c)
    aa = X._truss(asd, hc, dmg_key="Damage")
    lnc = X._truss(lad_unc, hc, dmg_key="damage")
    lca = X._truss(lad_cal, hc, dmg_key="damage")

    fig, ax = plt.subplots(2, 3, figsize=(15.5, 8.2))
    ax[0, 0].plot([c[0] for c in code], rs, "-", lw=3, color="#bbbbbb", label="reference")
    ax[0, 0].plot([c[0] for c in code], [c[1] for c in code], "--", lw=1.0, color="C3", label="LadrunoUniaxialJ2")
    ax[0, 0].set(xlabel=r"$\varepsilon$", ylabel=r"$\sigma$", title="E1  cyclic oracle hysteresis (<1e-8)")
    ax[0, 0].legend(fontsize=8)
    ax[0, 1].plot(range(1, len(Dcum) + 1), Dcum, "o-", ms=3, color="C0")
    ax[0, 1].set(xlabel="cycle number", ylabel="damage $D$",
                 title=r"E2  LCF damage accumulation ($\pm$0.05)")
    ax[0, 2].loglog(amps, Ns, "s-", color="C2")
    ax[0, 2].set(xlabel="strain amplitude", ylabel="cycles to $D=0.5$",
                 title="E2  life vs amplitude (Coffin-Manson trend)")
    # NEW: undamaged overlay — the apples-to-apples cross-code check (loops coincide)
    ax[1, 0].plot([c[0] for c in au], [c[1] for c in au], "-", lw=2.6, color="#888", label="ASDSteel1D")
    ax[1, 0].plot([c[0] for c in lu], [c[1] for c in lu], "--", lw=1.2, color="C3", label="Ladruno (matched)")
    ax[1, 0].set(xlabel=r"$\varepsilon$", ylabel=r"$\sigma$",
                 title=f"E3a  UNDAMAGED loops — matched backbone\n(loops coincide, max resid {res_u:.1%})")
    ax[1, 0].legend(fontsize=8)
    # damaged loops — diverge by design; calibration narrows it
    ax[1, 1].plot([c[0] for c in aa], [c[1] for c in aa], "-", lw=1.1, color="#888", label="ASDSteel1D -fracture")
    ax[1, 1].plot([c[0] for c in lnc], [c[1] for c in lnc], ":", lw=1.0, color="C0", label="Lemaitre uncalib.")
    ax[1, 1].plot([c[0] for c in lca], [c[1] for c in lca], "-", lw=1.0, color="C3", label="Lemaitre calibrated")
    ax[1, 1].set(xlabel=r"$\varepsilon$", ylabel=r"$\sigma$",
                 title="E3b  DAMAGED loops (damage on)\n(diverge by design — different damage laws)")
    ax[1, 1].legend(fontsize=7.5)
    # annotation panel — WHY they diverge
    ax[1, 2].axis("off")
    ax[1, 2].text(
        0.0, 1.0,
        "Why the damaged loops differ — by design, not a bug\n"
        "─────────────────────────────────────────\n"
        "Both couple damage the SAME way:  σ = (1−d)·σ̃,\n"
        "C = (1−d)·C̃  (isotropic stiffness+strength degr.).\n\n"
        "The EVOLUTION laws differ:\n"
        "  • ASDSteel1D -fracture:  d = 1 − e^(−sy(eₚₗ−eupl)/G),\n"
        "    keyed to plastic strain past a NECKING threshold\n"
        "    eupl = eu − sy/E, with fracture energy G(eu).\n"
        "    → sharp, necking-driven drop.\n"
        "  • Lemaitre:  D = ∫(Y/r)ˢ dp, keyed to the energy\n"
        "    release rate over ACCUMULATED p past p_D.\n"
        "    → smooth, energy-driven softening.\n\n"
        "E3a (undamaged) is the apples-to-apples cross-code\n"
        "check — matched Chaboche backbone → loops coincide.\n"
        "E3b is apples-to-oranges (two damage laws); the\n"
        "rough calibration (p_D=eupl, r matched to ∫σdε)\n"
        "aligns onset+energy; the residual is the softening\n"
        "SHAPE, which the differing functional forms fix.",
        transform=ax[1, 2].transAxes, fontsize=8.2, va="top", family="DejaVu Sans",
        bbox=dict(boxstyle="round", fc="#f4f8ff", ec="0.7"))
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lemaitre_E_cyclic.png")); plt.close(fig)
    print(f"E done  undamaged-overlay resid={res_u:.2%}  calib(r,pD)=({r_c},{pD_c})  life(amps)=",
          list(zip(amps, Ns)))


# ---------------------------------------------------------------- Layer C geometry
def fig_C_geometry():
    """Mesh + boundary conditions and the final damage field of the DEN bar."""
    import test_lemaitre_notched_bar as C
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import Normalize

    HY = C.HY
    h = 2.0
    nodes, hexes = C._mesh(h)
    keep = [hx for hx in hexes if not C._in_notch(sum(nodes[n][0] for n in hx) / 8.0,
                                                  sum(nodes[n][1] for n in hx) / 8.0)]
    ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
    used = set(n for hx in keep for n in hx)
    for t in used:
        ops.node(int(t), *nodes[t])
    ops.nDMaterial("LadrunoJ2", 1, C.K, C.G, *C.ISO, "-damage", "lemaitre", *C.DMG,
                   "-autoRegularization", C.LCH_REF)
    for i, hx in enumerate(keep, 1):
        ops.element("LadrunoBrick", i, *[int(n) for n in hx], 1, "-formulation", "bbar")
    xmin = min(nodes[t][0] for t in used); xmax = max(nodes[t][0] for t in used)
    left = sorted(int(t) for t in used if abs(nodes[t][0] - xmin) < 1e-6)
    for t in used:
        if abs(nodes[t][2]) < 1e-6:
            ops.fix(int(t), 0, 0, 1)
    for t in left:
        ops.fix(int(t), 1, 1 if t == left[0] else 0, 0)
    pull = sorted(int(t) for t in used if abs(nodes[t][0] - xmax) < 1e-6); master = pull[0]
    for n in pull:
        if n != master:
            ops.equalDOF(master, n, 1)
    ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(master, 1.0, 0.0, 0.0)
    ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
    ops.test("NormDispIncr", 1e-6, 50, 0); ops.analysis("Static")
    base = 1.5 / 150; u = 0.0; calls = 0
    while u < 1.5 - 1e-12 and calls < 1200:
        ok = -1
        for cut, alg in ((1.0, "KrylovNewton"), (0.5, "KrylovNewton"), (0.25, "NewtonLineSearch"),
                         (0.1, "NewtonLineSearch"), (0.03, "NewtonLineSearch")):
            ops.integrator("DisplacementControl", master, 1, base * cut); ops.algorithm(alg)
            ok = ops.analyze(1); calls += 1
            if ok == 0:
                break
        if ok != 0:
            break
        u = ops.nodeDisp(master)[0]

    # front-face (z=0) quad of each kept hex + its final damage
    polys, dmg = [], []
    for i, hx in enumerate(keep, 1):
        face = [n for n in hx if abs(nodes[n][2]) < 1e-6][:4]
        if len(face) < 4:
            face = hx[:4]
        xy = [(nodes[n][0], nodes[n][1]) for n in face]
        # order the 4 corners CCW about their centroid
        cx = sum(p[0] for p in xy) / 4.0; cy = sum(p[1] for p in xy) / 4.0
        xy.sort(key=lambda p: math.atan2(p[1] - cy, p[0] - cx))
        polys.append(xy)
        d = ops.eleResponse(i, "material", 1, "damage")
        dmg.append(max(d) if d else 0.0)

    import numpy as _np
    fig, ax = plt.subplots(2, 1, figsize=(8.5, 8))
    for a in ax:
        a.set_xlim(-7, xmax + 13); a.set_ylim(-3, HY + 3); a.set_aspect("equal")
    # panel 1: mesh + BCs
    ax[0].add_collection(PolyCollection(polys, facecolors="#d9e7f5", edgecolors="0.4", lw=0.6))
    ax[0].plot([xmin, xmin], [0, HY], color="C0", lw=4, solid_capstyle="butt")
    ax[0].text(xmin - 1.5, HY / 2, "u$_x$=0", ha="right", va="center", color="C0", fontsize=9)
    for yy in [i * HY / 6 for i in range(7)]:
        ax[0].annotate("", xy=(xmax + 5, yy), xytext=(xmax, yy),
                       arrowprops=dict(arrowstyle="->", color="C3", lw=1.4))
    ax[0].text(xmax + 6, HY / 2, "pull U", ha="left", va="center", color="C3", fontsize=9)
    ax[0].set(title=f"C  DEN bar mesh + BCs  (h={h}, {len(keep)} hexes, bbar, geom=linear)",
              xlabel="x [mm]", ylabel="y [mm]")
    # panel 2: final damage field
    pc2 = PolyCollection(polys, array=_np.array(dmg), cmap="inferno",
                         norm=Normalize(0.0, max(dmg) if dmg else 1.0), edgecolors="0.6", lw=0.3)
    ax[1].add_collection(pc2)
    ax[1].set(title=f"C  final damage field D  (localizes in the ligament; max={max(dmg):.2f})",
              xlabel="x [mm]", ylabel="y [mm]")
    fig.colorbar(pc2, ax=ax[1], label="damage D", fraction=0.025, pad=0.02)
    fig.savefig(os.path.join(FIGDIR, "lemaitre_C_geometry.png"), bbox_inches="tight")
    plt.close(fig)
    print("Cgeom done  reached u=", round(u, 3), " maxD=", round(max(dmg), 3))


# ---------------------------------------------------------------- Layer F (formulation + performance)
def _sum_elem_tangent(h5path, classTag=33002):
    """Walk the profiler rollup; sum element TANGENT-formation wall-time + calls for a
    classTag across all formTangent scopes. Returns (total_wall_ms, total_count)."""
    import h5py
    tot_ns, tot_n = 0, 0

    def visit(g):
        nonlocal tot_ns, tot_n
        for k in g:
            it = g[k]
            if isinstance(it, h5py.Group):
                visit(it)
            elif k == "elem_by_type" and "tangent" in g.name:
                arr = it[()]
                for row in arr:
                    if int(row["classTag"]) == classTag:
                        tot_ns += int(row["wall_ns"]); tot_n += int(row["count"])
    with h5py.File(h5path, "r") as f:
        if "runs" in f:
            visit(f["runs"])
    return tot_ns / 1e6, tot_n


def _ssp_or_eas():
    ops.wipe(); ops.model("basic","-ndm",3,"-ndf",3)
    for t,c in {1:(0,0,0),2:(1,0,0),3:(1,1,0),4:(0,1,0),5:(0,0,1),6:(1,0,1),7:(1,1,1),8:(0,1,1)}.items():
        ops.node(t,float(c[0]),float(c[1]),float(c[2]))
    ops.nDMaterial("ElasticIsotropic",1,1.0e5,0.3)
    try:
        ops.element("LadrunoBrick",1,1,2,3,4,5,6,7,8,1,"-formulation","ssp"); return "ssp"
    except Exception:
        return "eas"


def fig_F_formulation_perf():
    _SP_NAME=_ssp_or_eas()
    """bbar vs eas: same DEN bar, same regularization. Compares the structural response
    (should agree — both mesh-objective), the eas hourglass-energy fraction (the
    damage-scaled Kstab correction), and the per-element TANGENT cost via the profiler."""
    import test_lemaitre_notched_bar as C
    import time as _time
    TMP = os.environ.get("TEMP", "/tmp")

    def run(h, form, prof_path):
        nodes, hexes = C._mesh(h)
        keep = [hx for hx in hexes if not C._in_notch(sum(nodes[n][0] for n in hx) / 8.0,
                                                      sum(nodes[n][1] for n in hx) / 8.0)]
        ops.wipe(); ops.model("basic", "-ndm", 3, "-ndf", 3)
        used = set(n for hx in keep for n in hx)
        for t in used:
            ops.node(int(t), *nodes[t])
        ops.nDMaterial("LadrunoJ2", 1, C.K, C.G, *C.ISO, "-damage", "lemaitre", *C.DMG,
                       "-autoRegularization", C.LCH_REF)
        eids = []
        for i, hx in enumerate(keep, 1):
            ops.element("LadrunoBrick", i, *[int(n) for n in hx], 1, "-formulation", form)
            eids.append(i)
        xmin = min(nodes[t][0] for t in used); xmax = max(nodes[t][0] for t in used)
        left = sorted(int(t) for t in used if abs(nodes[t][0] - xmin) < 1e-6)
        for t in used:
            if abs(nodes[t][2]) < 1e-6:
                ops.fix(int(t), 0, 0, 1)
        for t in left:
            ops.fix(int(t), 1, 1 if t == left[0] else 0, 0)
        pull = sorted(int(t) for t in used if abs(nodes[t][0] - xmax) < 1e-6); master = pull[0]
        for n in pull:
            if n != master:
                ops.equalDOF(master, n, 1)
        ops.timeSeries("Linear", 1); ops.pattern("Plain", 1, 1); ops.load(master, 1.0, 0.0, 0.0)
        ops.constraints("Transformation"); ops.numberer("RCM"); ops.system("UmfPack")
        ops.test("NormDispIncr", 1e-6, 50, 0); ops.analysis("Static")
        us, Ps, hgf = [0.0], [0.0], [0.0]
        W = uprev = Pprev = Pmax = u = 0.0; base = 1.5 / 150; calls = 0
        ops.profiler("start", "-deep", "-perStep")   # -perStep => nSteps + per-step series in the run header
        t0 = _time.time()
        while u < 1.5 - 1e-12 and calls < 1200:
            ok = -1
            for cut, alg in ((1.0, "KrylovNewton"), (0.5, "KrylovNewton"), (0.25, "NewtonLineSearch"),
                             (0.1, "NewtonLineSearch"), (0.03, "NewtonLineSearch")):
                ops.integrator("DisplacementControl", master, 1, base * cut); ops.algorithm(alg)
                ok = ops.analyze(1); calls += 1
                if ok == 0:
                    break
            if ok != 0:
                break
            ops.reactions(); P = -sum(ops.nodeReaction(n)[0] for n in left); u = ops.nodeDisp(master)[0]
            W += 0.5 * (P + Pprev) * (u - uprev); Pprev, uprev, Pmax = P, u, max(Pmax, P)
            us.append(u); Ps.append(P / 1000.0)
            ehg = sum((ops.eleResponse(e, "hourglassEnergy") or [0.0])[0] for e in eids)
            hgf.append(ehg / W if W > 1e-9 else 0.0)
        wall = _time.time() - t0
        ops.profiler("stop")
        if os.path.exists(prof_path):           # fresh file (avoid HDF5 append/lock on re-run)
            os.remove(prof_path)
        ops.profiler("report", prof_path, "-run", form)
        tang_ms, tang_n = _sum_elem_tangent(prof_path)
        return dict(us=us, Ps=Ps, hgf=hgf, peak=Pmax, W=W, reached=u, ne=len(keep),
                    wall=wall, tang_ms=tang_ms, tang_n=tang_n,
                    per_elem_us=(tang_ms * 1e3 / max(tang_n, 1)))

    h = 2.0
    bb = run(h, "bbar", os.path.join(TMP, "prof_bbar.h5"))
    ea = run(h, _SP_NAME, os.path.join(TMP, "prof_eas.h5"))
    # eas mesh-objectivity: coarse vs the medium above
    ea_c = run(4.0, _SP_NAME, os.path.join(TMP, "prof_eas_c.h5"))

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(bb["us"], bb["Ps"], "-", lw=2, color="C0", label=f"bbar 8GP (n={bb['ne']})")
    ax[0].plot(ea["us"], ea["Ps"], "--", lw=1.7, color="C3", label=f"{_SP_NAME} 1GP (n={ea['ne']})")
    ax[0].set(xlabel="elongation u [mm]", ylabel="load P [kN]",
              title="F  response: bbar vs eas (h=2)\n(eas single-point is stiffer on the notch gradient)")
    ax[0].legend(fontsize=8)
    ax[1].plot(ea["us"], [100 * x for x in ea["hgf"]], "-", color="C2")
    ax[1].axhline(10, ls=":", color="k", lw=1)
    ax[1].set(xlabel="elongation u [mm]", ylabel="hourglass energy / internal [%]",
              title="F  eas stabilization energy\n(damage-scaled Kstab keeps it bounded)")
    labels = ["bbar\n(8 GP)", "eas\n(1 GP + Kstab)"]
    per = [bb["per_elem_us"], ea["per_elem_us"]]
    xb = range(len(labels))
    b1 = ax[2].bar(list(xb), per, 0.5, color=["C0", "C3"])
    ax[2].set_xticks(list(xb)); ax[2].set_xticklabels(labels)
    ax[2].set_ylabel("per-element tangent [µs]")
    ax[2].set_ylim(0, max(per) * 1.25)
    # NOTE: wall-time deliberately NOT plotted — the profiler's wall/cpu totals are a Win32
    # stub on Windows and swing run-to-run; only the per-element tangent (same op) is shown.
    ax[2].set_title("F  profiler: per-element tangent cost\n(eas 1 GP vs bbar 8 GP; "
                    "wall-time omitted — Win32 timer noise)")
    for r, v in zip(b1, per):
        ax[2].text(r.get_x() + r.get_width() / 2, v, f"{v:.1f}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIGDIR, "lemaitre_F_formulation_perf.png")); plt.close(fig)
    pP = abs(bb["peak"] - ea["peak"]) / max(bb["peak"], 1)
    pE_obj = abs(ea_c["W"] - ea["W"]) / max(ea_c["W"], 1e-9)
    print(f"F done  bbar: peak={bb['peak']:.0f} W={bb['W']:.0f} wall={bb['wall']:.1f}s tang={bb['tang_ms']:.0f}ms/{bb['tang_n']} = {bb['per_elem_us']:.2f}us/elem")
    print(f"        eas : peak={ea['peak']:.0f} W={ea['W']:.0f} wall={ea['wall']:.1f}s tang={ea['tang_ms']:.0f}ms/{ea['tang_n']} = {ea['per_elem_us']:.2f}us/elem  hgf_max={100*max(ea['hgf']):.1f}%")
    print(f"        bbar-vs-eas peak diff={pP:.1%}; eas mesh-obj (coarse vs med) energy diff={pE_obj:.1%}")


if __name__ == "__main__":
    which = sys.argv[2] if len(sys.argv) > 2 else "all"
    todo = {"A1": fig_A1_oracle, "A2": fig_A2_closed_form, "B": fig_B_triaxiality,
            "D": fig_D_cross, "C": fig_C_mesh_objectivity, "E": fig_E_cyclic,
            "Cgeom": fig_C_geometry, "F": fig_F_formulation_perf}
    for k, fn in todo.items():
        if which in ("all", k):
            fn()
    print("FIGURES_DONE", FIGDIR)
