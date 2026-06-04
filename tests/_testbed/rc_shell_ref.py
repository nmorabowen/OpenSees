"""
rc_shell_ref.py — numpy reference oracle for the Ladruno RC plane-stress / plate-fiber
constitutive kernel (Phase 1).

This is a faithful, OpenSees-free clone of ASDConcrete3DMaterial::compute()
(SRC/material/nD/ASDConcrete3DMaterial.cpp:2322-2523), plus the Phase-1 MCFT
compression-softening addition derived + adversarially verified in the Phase-1 blueprint
(ADR 19_ladruno_rc_shell_adr.md). It exists to PROVE the blocking correctness gate (A1)
before any C++ is written, and later to serve as the trajectory oracle for the C++ tests.

The single Phase-1 physics change (Derivation B, verified gate_holds=true):
    beta multiplies the ASSEMBLED compressive cone at final nominal-stress assembly only:
        sigma = (1 - dt_bar)*ST  +  beta * (1 - dc_bar)*SC
    beta = clamp(1/(0.8 + 170*eps_1), betaFloor, 1),  eps_1 = in-plane principal TENSILE strain.
    It is NOT applied to the strain abscissa, the equivalent-strain measure, the effective
    stress reassembly, dc_plastic, or dc_bar. `betaMode='abscissa'` reproduces the FORBIDDEN
    insertion (scale the compressive abscissa) so the test can show it misses the peak.

Voigt order matches the C++ spine exactly: {00,11,22,01,12,02} (indices 0..5).
Run standalone (any numpy interpreter, no OpenSees):  python tests/_testbed/rc_shell_ref.py
"""

import numpy as np

YTOL = 1.0e-12
XTOL = 1.0e-12
STRESS_TOL = 1.0e-6   # ht.stressTolerance() analog for the tensile-measure gate


# ---------------------------------------------------------------------------
# Hardening backbone  (clone of HardeningLaw::evaluateAt, cpp:1010-1054)
# A backbone is a list of (x, y, q) triplets:  x=total equiv strain,
# y=nominal stress, q=effective(plasticity) stress.  d = 1 - y/q.
# ---------------------------------------------------------------------------
class Backbone:
    def __init__(self, pts):
        self.x = np.array([p[0] for p in pts], float)
        self.y = np.array([p[1] for p in pts], float)
        self.q = np.array([p[2] for p in pts], float)
        assert len(self.x) >= 2

    def max_stress(self):
        return float(np.max(self.y))

    def evaluate_at(self, x):
        """Return (x, y, d, q) at abscissa x — piecewise linear, with the
        positive-tangent extrapolation rule beyond the last point."""
        xs, ys, qs = self.x, self.y, self.q
        n = len(xs)
        found = False
        for i in range(1, n):
            if x <= xs[i] + XTOL:
                x1, x2 = xs[i - 1], xs[i]
                y1, y2 = ys[i - 1], ys[i]
                q1, q2 = qs[i - 1], qs[i]
                found = True
                break
        if not found:
            x1 = xs[-1]
            x2 = x
            span = x - x1
            y1 = ys[-1]
            tan_y = (ys[-1] - ys[-2]) / (xs[-1] - xs[-2])
            y2 = y1 + span * tan_y if tan_y > 0.0 else y1
            q1 = qs[-1]
            tan_q = (qs[-1] - qs[-2]) / (xs[-1] - xs[-2])
            q2 = q1 + span * tan_q if tan_q > 0.0 else q1
        xspan = x2 - x1
        xratio = (x - x1) / xspan if xspan > 0.0 else 0.0
        y = max(YTOL, y1 + (y2 - y1) * xratio)
        q = max(YTOL, q1 + (q2 - q1) * xratio)
        d = 1.0 - y / q
        return (x, y, d, q)


def plastic_strain(pt, E):
    # HardeningLawPoint::plasticStrain(E) = x - q/E
    return pt[0] - pt[3] / E


# ---------------------------------------------------------------------------
# isotropic elastic 6x6  (clone of getInitialTangent, cpp:1733-1743)
# Voigt {00,11,22,01,12,02}; shear rows carry mu (engineering shear strain).
# ---------------------------------------------------------------------------
def elastic_C0(E, nu):
    C = np.zeros((6, 6))
    mu2 = E / (1.0 + nu)
    lam = nu * mu2 / (1.0 - 2.0 * nu)
    mu = 0.5 * mu2
    mu2 += lam
    C[0, 0] = C[1, 1] = C[2, 2] = mu2
    C[0, 1] = C[1, 0] = C[0, 2] = C[2, 0] = C[1, 2] = C[2, 1] = lam
    C[3, 3] = C[4, 4] = C[5, 5] = mu
    return C


def macauley(x):
    return x if x > 0.0 else 0.0


def lubliner(s1, s2, s3, ft, fc, k1, scale, Kc):
    """clone of lublinerCriterion (cpp:2483-2497). s1>=s2>=s3 (sorted desc).
    Local Lubliner shape parameter renamed beta_lub to avoid colliding with MCFT beta."""
    fb = 1.16 * fc
    gamma = 3.0 * (1.0 - Kc) / (2.0 * Kc - 1.0)
    alpha = (fb - fc) / (2.0 * fb - fc)
    I1 = s1 + s2 + s3
    p = I1 / 3.0
    sdev = np.array([s1 - p, s2 - p, s3 - p])
    J2 = 0.5 * float(sdev @ sdev)
    beta_lub = fc / ft * (1.0 - alpha) - (1.0 + alpha)
    smax = macauley(s1)
    smin = macauley(-s1)
    return (1.0 / (1.0 - alpha) * (alpha * I1 + np.sqrt(3.0 * J2)
            + k1 * beta_lub * smax - gamma * smin)) * scale


# ---------------------------------------------------------------------------
# spectral split  (clone of StressDecomposition::compute, cpp:781-851)
# returns Si (desc), eigvecs V (cols), ST, SC (6-Voigt), R
# ---------------------------------------------------------------------------
def stress_decomposition(S6, cdf):
    M = np.array([[S6[0], S6[3], S6[5]],
                  [S6[3], S6[1], S6[4]],
                  [S6[5], S6[4], S6[2]]])
    w, V = np.linalg.eigh(M)          # ascending
    order = np.argsort(w)[::-1]       # -> descending, so Si[0] is the largest principal
    Si = w[order]
    V = V[:, order]

    def heavyside(x):
        return 1.0 if x > 0.0 else (0.0 if x < 0.0 else 0.5)

    # split each principal value into +/- parts (method 2: PT, PC, with PO split 50/50)
    st_pr = np.zeros(3)
    sc_pr = np.zeros(3)
    for j in range(3):
        hp = heavyside(Si[j])
        hn = heavyside(-Si[j])
        po = 1.0 - hp - hn            # 1 only if Si[j]==0
        st_pr[j] = (hp + 0.5 * po) * Si[j]
        sc_pr[j] = (hn + 0.5 * po) * Si[j]

    def recompose(spr):
        # clone of recompose(): build 6-Voigt stress from principal values spr & V
        Sv = np.zeros(6)
        Sv[0] = spr[0]*V[0,0]**2 + spr[1]*V[0,1]**2 + spr[2]*V[0,2]**2
        Sv[1] = spr[0]*V[1,0]**2 + spr[1]*V[1,1]**2 + spr[2]*V[1,2]**2
        Sv[2] = spr[0]*V[2,0]**2 + spr[1]*V[2,1]**2 + spr[2]*V[2,2]**2
        Sv[3] = spr[0]*V[0,0]*V[1,0] + spr[1]*V[0,1]*V[1,1] + spr[2]*V[0,2]*V[1,2]
        Sv[4] = spr[0]*V[1,0]*V[2,0] + spr[1]*V[1,1]*V[2,1] + spr[2]*V[1,2]*V[2,2]
        Sv[5] = spr[0]*V[0,0]*V[2,0] + spr[1]*V[0,1]*V[2,1] + spr[2]*V[0,2]*V[2,2]
        return Sv

    ST = recompose(st_pr)
    SC = recompose(sc_pr)

    Rnum = Rden = 0.0
    for j in range(3):
        Sj = Si[j]
        if Sj > 0.0:
            Sj = Sj * cdf
        if Sj > 0.0:
            Rnum += Sj
        Rden += abs(Sj)
    R = Rnum / Rden if Rden > 0.0 else 0.0
    return Si, V, ST, SC, R


def beta_compr(e1, floorv):
    b = 1.0 / (0.8 + 170.0 * e1)
    if b > 1.0:
        b = 1.0
    if b < floorv:
        b = floorv
    return b


def eps1_from_membrane(exx, eyy, gxy):
    """Largest eigenvalue of the in-plane 2x2 TOTAL strain tensor, clamped >=0.
    gxy is engineering shear; tensor shear = gxy/2."""
    c = 0.5 * (exx + eyy)
    r = np.sqrt((0.5 * (exx - eyy))**2 + (0.5 * gxy)**2)
    return max(0.0, c + r)


# ---------------------------------------------------------------------------
# Material params + state
# ---------------------------------------------------------------------------
class Params:
    def __init__(self, E, nu, Kc, ht, hc,
                 beta_on=False, lubliner_tc_reduced=False, beta_floor=0.1,
                 cdf=0.0, eta=0.0):
        self.E = E
        self.nu = nu
        self.Kc = Kc
        self.ht = ht
        self.hc = hc
        self.beta_on = beta_on
        self.lubliner_tc_reduced = lubliner_tc_reduced
        self.beta_floor = beta_floor
        self.cdf = cdf
        self.eta = eta
        # fcft_ratio = max(5, fcmax/ftmax)   (cpp:1585-1590)
        fcmax = hc.max_stress()
        ftmax = ht.max_stress()
        self.fcft_ratio = max(5.0, fcmax / ftmax) if ftmax > 0.0 else 1.0


class State:
    """committed state, mirrors the RCHist fields used in Phase 1"""
    def __init__(self):
        self.stress_eff = np.zeros(6)
        self.strain = np.zeros(6)
        self.xt = 0.0   # committed tensile equiv strain (svt)
        self.xc = 0.0   # committed compressive equiv strain (svc)

    def copy(self):
        s = State()
        s.stress_eff = self.stress_eff.copy()
        s.strain = self.strain.copy()
        s.xt = self.xt
        s.xc = self.xc
        return s


def equiv_tensile(Si, fcft_ratio, Kc):
    if Si[0] < STRESS_TOL:
        return 0.0
    return lubliner(Si[0], Si[1], Si[2], 1.0 / fcft_ratio, 1.0, 1.0, 1.0 / fcft_ratio, Kc)


def equiv_compressive(Si, fcft_ratio, Kc, k1):
    s = [min(v, 0.0) for v in Si]
    return lubliner(s[0], s[1], s[2], 1.0 / fcft_ratio, 1.0, k1, 1.0, Kc)


def compute(P, st_committed, strain6, betaMode='strength'):
    """One trial evaluation at total strain `strain6` from committed state `st_committed`.
    Returns (sigma6, info_dict, new_state).  Phase-1: implex off, eta-rate off by default.
    betaMode: 'strength' (verified), 'abscissa' (forbidden, for the discrimination test),
              'off' (beta=1 -> byte-faithful spine)."""
    E, nu = P.E, P.nu
    C0 = elastic_C0(E, nu)

    st = st_committed.copy()
    # elastic predictor: SEFFn = SEFF_commit + C0:(En - En-1)
    dstrain = strain6 - st_committed.strain
    stress_eff = st_committed.stress_eff + C0 @ dstrain

    Si, V, ST, SC, R = stress_decomposition(stress_eff, P.cdf)

    # MCFT beta from the in-plane membrane tensile principal strain (independent of eps_33)
    e1 = eps1_from_membrane(strain6[0], strain6[1], 2.0 * strain6[3])  # strain6[3] is tensor; *2 -> engineering
    if betaMode == 'off' or not P.beta_on:
        beta = 1.0
    else:
        beta = beta_compr(e1, P.beta_floor)

    k1_c = 0.0 if (P.beta_on and P.lubliner_tc_reduced and betaMode != 'off') else 1.0

    # committed hardening points
    pt = P.ht.evaluate_at(st_committed.xt)
    pc = P.hc.evaluate_at(st_committed.xc)
    xt_pl = plastic_strain(pt, E)
    xc_pl = plastic_strain(pc, E)

    # trial measures (/E converts the Lubliner STRESS measure to the strain abscissa, cpp:2509/2522)
    xt_meas = equiv_tensile(Si, P.fcft_ratio, P.Kc) / E
    xc_meas = equiv_compressive(Si, P.fcft_ratio, P.Kc, k1_c) / E
    # FORBIDDEN insertion (for discrimination only): scale the compressive abscissa
    if betaMode == 'abscissa' and P.beta_on:
        xc_meas = xc_meas * beta

    xt_trial = xt_meas + xt_pl
    xc_trial = xc_meas + xc_pl

    rc1 = P.eta / (P.eta + 1.0) if P.eta > 0.0 else 0.0   # dtime=1 convention; eta=0 -> rc1=0
    rc2 = 1.0 - rc1

    new_xt = st_committed.xt
    new_xc = st_committed.xc
    if xt_trial > pt[0]:
        xt_new = rc1 * pt[0] + rc2 * xt_trial
        pt = P.ht.evaluate_at(xt_new)
        new_xt = pt[0]
    if xc_trial > pc[0]:
        xc_new = rc1 * pc[0] + rc2 * xc_trial
        pc = P.hc.evaluate_at(xc_new)
        new_xc = pc[0]

    # plastic damage (uses committed xt_pl/xc_pl, advanced pt/pc) — cpp:2430-2433
    seff_eq_t = (pt[0] - xt_pl) * E
    dt_plastic = 1.0 - pt[3] / seff_eq_t if seff_eq_t > 0.0 else 0.0
    seff_eq_c = (pc[0] - xc_pl) * E
    dc_plastic = 1.0 - pc[3] / seff_eq_c if seff_eq_c > 0.0 else 0.0

    # mix damages (cpp:2436-2449); with cdf=0 -> R=0 -> identity
    def mix_dam(dt, dc):
        return 1.0 - (1.0 - dc) * (1.0 - R * dt)
    auxp = P.hc.evaluate_at(pt[0] * P.fcft_ratio)
    aux_pl = plastic_strain(auxp, E)
    aux_seff = (auxp[0] - aux_pl) * E
    aux_d = 1.0 - auxp[3] / aux_seff if aux_seff > 0.0 else 0.0
    dc_plastic = mix_dam(aux_d, dc_plastic)
    pc_d = mix_dam(auxp[2], pc[2])

    dt_bar = pt[2] + dt_plastic - pt[2] * dt_plastic
    dc_bar = pc_d + dc_plastic - pc_d * dc_plastic

    # nominal stress — THE Phase-1 insertion: beta multiplies the assembled compressive cone
    sigma = (1.0 - dt_bar) * ST + beta * (1.0 - dc_bar) * SC

    # effective stress reassembly for next-step commit (beta NOT applied — cpp:2451-2453)
    stress_eff = (1.0 - dt_plastic) * ST + (1.0 - dc_plastic) * SC

    st.stress_eff = stress_eff
    st.strain = strain6.copy()
    st.xt = new_xt
    st.xc = new_xc

    info = dict(Si=Si, beta=beta, dt_bar=dt_bar, dc_bar=dc_bar, e1=e1,
                sigma=sigma, dc_plastic=dc_plastic)
    return sigma, info, st


# ---------------------------------------------------------------------------
# reference fixture (the curved backbone from the blueprint: fc'=30 @ x*=0.002,
# q*=40 -> d*=0.25, softening to y=5 @ x=0.01)
# ---------------------------------------------------------------------------
def reference_params(beta_on=False, lubliner_tc_reduced=False):
    hc = Backbone([
        (0.0,    0.0,   1.0e-12),
        (0.0007, 24.0,  24.0),    # ~elastic, d=0
        (0.0020, 30.0,  40.0),    # peak nominal fc'=30, d=0.25
        (0.0100,  5.0,  45.0),    # softened, d=1-5/45
    ])
    ht = Backbone([
        (0.0,     0.0,  1.0e-12),
        (0.0001,  3.0,  3.0),     # peak tensile ft=3
        (0.0010,  0.5,  5.0),
    ])
    return Params(E=30000.0, nu=0.2, Kc=0.667, ht=ht, hc=hc,
                  beta_on=beta_on, lubliner_tc_reduced=lubliner_tc_reduced)


def peak_compressive_under_transverse_tension(eps1_fixed, betaMode, nsteps=400, eyy_max=-0.006):
    """A1 driver: hold transverse tensile membrane strain eps_xx = eps1_fixed (axis 0),
    ramp compression eps_yy 0 -> eyy_max (axis 1, the in-plane strut), eps_zz=0, no shear.
    Return the peak magnitude of the most-compressive principal of the NOMINAL stress."""
    beta_on = (betaMode != 'off')
    P = reference_params(beta_on=beta_on, lubliner_tc_reduced=beta_on)
    st = State()
    peak = 0.0
    for k in range(1, nsteps + 1):
        eyy = eyy_max * k / nsteps
        strain6 = np.array([eps1_fixed, eyy, 0.0, 0.0, 0.0, 0.0])
        sigma, info, st = compute(P, st, strain6, betaMode=betaMode)
        M = np.array([[sigma[0], sigma[3], sigma[5]],
                      [sigma[3], sigma[1], sigma[4]],
                      [sigma[5], sigma[4], sigma[2]]])
        smin = float(np.min(np.linalg.eigvalsh(M)))   # most compressive principal
        peak = max(peak, -smin)                         # magnitude
    return peak


def run_A1(verbose=True):
    """THE blocking gate (telescoping identity). For each eps_1 the baseline is the SAME
    eps_1 with beta OFF (identical spine/confinement state); beta is then the ONLY
    difference, so:
        peak(strength, eps_1) / peak(off, eps_1)  ==  beta(eps_1)   EXACTLY.
    The abscissa-mode insertion shifts the hardening lookup instead of scaling the cone,
    so its ratio must NOT equal beta wherever beta<1.  (beta clamps to 1 for
    eps_1 < ~1.2e-3, where all three modes coincide trivially.)"""
    sweep = [0.0, 1.0e-4, 5.0e-4, 1.0e-3, 3.0e-3, 1.0e-2]
    rows = []
    ok = True
    for e1 in sweep:
        b = beta_compr(e1, 0.1)
        base_e1 = peak_compressive_under_transverse_tension(e1, 'off')
        p_str = peak_compressive_under_transverse_tension(e1, 'strength')
        p_abs = peak_compressive_under_transverse_tension(e1, 'abscissa')
        ratio_str = p_str / base_e1
        ratio_abs = p_abs / base_e1
        err_str = abs(ratio_str - b)
        if b < 1.0 - 1.0e-9:
            # discriminating point: strength==beta exactly; abscissa must miss it
            gate = (err_str < 1.0e-6) and (abs(ratio_abs - b) > 1.0e-3)
        else:
            # beta clamped to 1: strength / abscissa / off must all coincide
            gate = (err_str < 1.0e-6) and (abs(ratio_abs - 1.0) < 1.0e-9)
        ok = ok and gate
        rows.append((e1, b, base_e1, ratio_str, ratio_abs, err_str, gate))
    if verbose:
        print(f"  {'eps_1':>8} {'beta':>8} {'base(off)':>11} {'str/base':>10} {'abs/base':>10} {'|str-beta|':>11}  gate")
        for e1, b, base_e1, rs, ra, es, g in rows:
            print(f"  {e1:8.1e} {b:8.4f} {base_e1:11.4f} {rs:10.6f} {ra:10.6f} {es:11.2e}  {'PASS' if g else 'FAIL'}")
    return ok, rows


if __name__ == '__main__':
    print("=== A1: beta-on-strength closed-form gate (numpy oracle) ===")
    ok, _ = run_A1()
    print(f"\nA1 GATE: {'PASS — beta scales the strength axis exactly' if ok else 'FAIL'}")
    raise SystemExit(0 if ok else 1)
