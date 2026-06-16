"""
concrete3d_ref.py — numpy reference oracle for LadrunoConcrete3D (CDPM2-grade solid concrete).
ADR: Ladruno_implementation/31_ladruno_concrete3d_adr.md.

P0 SCOPE (this file, gateable BEFORE any return-map / C++ code — red-team blocking fix):
    Implement the fc-NORMALIZED Menetrey-Willam three-invariant FAILURE surface and PROVE it is
    correctly normalized:
      (G1) uniaxial compression  sigma=(-fc,0,0)   lands on f=0 EXACTLY (independent of m0, e),
      (G2) uniaxial tension      sigma=(+ft,0,0)   lands on f=0 EXACTLY  => fixes m0,
      (G3) deviatoric meridian ratio  rho_t/rho_c == e   (the eccentricity IS that ratio),
      (G4) KUPFER equibiaxial: solve e from the target equibiaxial strength fcc/fc (default 1.16),
           recovering the canonical concrete eccentricity e ~ 0.52  (ADR 4.1 option (b)).
    Also: the lateral-sigma3 MIXED-CONTROL residual contract (a stub here — the constant-/passive-
    confinement boundary condition that P1's return map and §4.6's confined-fiber view consume).

Surface (ADR 4.1; PIN against Grassl et al. 2013 IJSS 50:3805-3816 + Menetrey-Willam 1995 by Eq.
number at code-review). With xi = I1/sqrt(3), rho = sqrt(2*J2), Lode angle theta in [0, pi/3]:

    f(xi,rho,theta) = (sqrt(1.5)*rho/fc)^2
                    + m0 * ( rho*r(theta,e)/(sqrt(6)*fc) + xi/(sqrt(3)*fc) )
                    - 1                                   (un-hardened: qh1=qh2=1)
    m0 = 3*(fc^2 - ft^2)/(fc*ft) * e/(e+1)
    r(theta,e) = Willam-Warnke elliptic Lode function, e in [0.5, 1]  (0.5 = convexity limit)

Sign convention: COMPRESSION NEGATIVE (continuum-mechanics / OpenSees). fc, ft are entered as
POSITIVE magnitudes; compressive stresses are negative.
Voigt order matches the C++ spine: {00,11,22,01,12,02}.
Run standalone (any numpy interpreter, no OpenSees):  python tests/_testbed/concrete3d_ref.py
"""

import numpy as np

SQRT3 = np.sqrt(3.0)
SQRT6 = np.sqrt(6.0)
SQRT1_5 = np.sqrt(1.5)
TOL_ONSURF = 1.0e-12   # algebraic on-surface identities (G1/G2) must be machine-exact
TOL_RATIO = 1.0e-10    # meridian-ratio identity (G3)


# ---------------------------------------------------------------------------
# Stress invariants  (xi, rho, Lode angle theta)
# ---------------------------------------------------------------------------
def invariants(sig):
    """sig = Voigt {s00,s11,s22,s01,s12,s02}. Returns (xi, rho, theta, I1, J2, J3)."""
    s = np.asarray(sig, float)
    s00, s11, s22, s01, s12, s02 = s
    I1 = s00 + s11 + s22
    p = I1 / 3.0
    # deviator
    d00, d11, d22 = s00 - p, s11 - p, s22 - p
    J2 = 0.5 * (d00 * d00 + d11 * d11 + d22 * d22) + s01 * s01 + s12 * s12 + s02 * s02
    # J3 = det of the deviatoric tensor
    J3 = (d00 * (d11 * d22 - s12 * s12)
          - s01 * (s01 * d22 - s12 * s02)
          + s02 * (s01 * s12 - d11 * s02))
    xi = I1 / SQRT3
    rho = np.sqrt(max(2.0 * J2, 0.0))
    # Lode angle: cos(3 theta) = (3*sqrt(3)/2) * J3 / J2^(3/2),  theta in [0, pi/3]
    if J2 <= 1.0e-300:
        theta = 0.0
    else:
        c3 = (3.0 * SQRT3 / 2.0) * J3 / (J2 ** 1.5)
        c3 = min(1.0, max(-1.0, c3))
        theta = np.arccos(c3) / 3.0
    return xi, rho, theta, I1, J2, J3


def principal_to_voigt(s1, s2, s3):
    return np.array([s1, s2, s3, 0.0, 0.0, 0.0], float)


# ---------------------------------------------------------------------------
# Willam-Warnke elliptic Lode function r(theta, e)
#   verified by hand (ADR review):  r(theta=0)=1/e (tensile meridian),
#                                    r(theta=pi/3)=1 (compressive meridian)
#   => deviatoric meridian ratio  rho_t/rho_c = r_c/r_t = e
# ---------------------------------------------------------------------------
def lode_r(theta, e):
    # Convex for e in [0.5,1] (Willam-Warnke 1975). NOTE: a naive per-sextant 1/r second-derivative
    # ("g+g''") polar test reports false violations near theta~52deg for larger e and is NOT the
    # correct convexity criterion for the elliptic interpolation — do not "fix" a non-bug.
    ct = np.cos(theta)
    num = 4.0 * (1.0 - e * e) * ct * ct + (2.0 * e - 1.0) ** 2
    rad = 4.0 * (1.0 - e * e) * ct * ct + 5.0 * e * e - 4.0 * e
    rad = max(rad, 0.0)
    den = 2.0 * (1.0 - e * e) * ct + (2.0 * e - 1.0) * np.sqrt(rad)
    return num / den


def m0_of(fc, ft, e):
    return 3.0 * (fc * fc - ft * ft) / (fc * ft) * e / (e + 1.0)


def yield_f(sig, fc, ft, e, qh1=1.0, qh2=1.0):
    """fc-normalized Menetrey-Willam. qh1=qh2=1 => failure surface. COMPRESSION NEGATIVE."""
    xi, rho, theta, *_ = invariants(sig)
    m0 = m0_of(fc, ft, e)
    r = lode_r(theta, e)
    term_quad = (SQRT1_5 * rho / fc) ** 2
    term_lin = m0 * qh1 * qh2 * (rho * r / (SQRT6 * fc) + xi / (SQRT3 * fc))
    return term_quad + term_lin - (qh1 * qh1) * (qh2 * qh2)


# ---------------------------------------------------------------------------
# Meridian radii rho(xi) on the tensile (theta=0) and compressive (theta=pi/3)
# meridians, by root-finding rho on the failure surface at a fixed xi.
# ---------------------------------------------------------------------------
def _rho_on_meridian(xi, theta, fc, ft, e):
    """Solve f=0 for rho at fixed (xi, theta). f is quadratic in rho => closed form."""
    m0 = m0_of(fc, ft, e)
    r = lode_r(theta, e)
    A = 1.5 / (fc * fc)                  # coeff of rho^2
    B = m0 * r / (SQRT6 * fc)            # coeff of rho^1
    C = m0 * xi / (SQRT3 * fc) - 1.0     # constant
    disc = B * B - 4.0 * A * C
    if disc < 0.0:
        return np.nan
    return (-B + np.sqrt(disc)) / (2.0 * A)


def meridian_ratio(fc, ft, e, xi_over_fc=-1.0):
    """rho_t/rho_c at a representative compressive xi (default xi = -fc)."""
    xi = xi_over_fc * fc
    rho_t = _rho_on_meridian(xi, 0.0, fc, ft, e)         # tensile meridian
    rho_c = _rho_on_meridian(xi, np.pi / 3.0, fc, ft, e)  # compressive meridian
    return rho_t / rho_c


# ---------------------------------------------------------------------------
# Kupfer equibiaxial strength  fcc = |sigma| at sigma=(-fcc,-fcc,0) on the surface,
# and the inverse: solve eccentricity e so fcc/fc == target (ADR 4.1 option b).
# ---------------------------------------------------------------------------
def equibiaxial_strength(fc, ft, e):
    """Return fcc (positive magnitude) at biaxial compression on the failure surface."""
    def f_of_b(b):
        return yield_f(principal_to_voigt(-b, -b, 0.0), fc, ft, e)
    lo, hi = 0.5 * fc, 3.0 * fc
    flo, fhi = f_of_b(lo), f_of_b(hi)
    # f increases with b (more compression -> off the surface positive); bracket sign change
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = f_of_b(mid)
        if flo * fm <= 0.0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        if hi - lo < 1.0e-12 * fc:
            break
    return 0.5 * (lo + hi)


def eccentricity_from_kupfer(fc, ft, target_fcc_ratio=1.16):
    """Solve e in [0.5,1] so equibiaxial fcc/fc == target. ADR 4.1(b)."""
    def g(e):
        return equibiaxial_strength(fc, ft, e) / fc - target_fcc_ratio
    lo, hi = 0.5 + 1.0e-6, 1.0 - 1.0e-9
    glo, ghi = g(lo), g(hi)
    if glo * ghi > 0.0:
        raise RuntimeError(f"no e in [0.5,1] hits fcc/fc={target_fcc_ratio} for ft/fc={ft/fc:.3f}"
                           f" (g(lo)={glo:+.4f}, g(hi)={ghi:+.4f})")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        gm = g(mid)
        if glo * gm <= 0.0:
            hi, ghi = mid, gm
        else:
            lo, glo = mid, gm
        if hi - lo < 1.0e-12:
            break
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Lateral-sigma3 MIXED-CONTROL residual contract (stub for P1 + §4.6).
# The confined uniaxial / triaxial driver holds the axial strain and Newton-solves the
# lateral strains so a chosen lateral residual vanishes:
#     free surface (plain reduction):   R = sigma_lat                  = 0
#     ACTIVE confinement (prescribed):  R = sigma_lat + p              = 0
#     PASSIVE confinement (hoop spring): R = sigma_lat + sigma_hoop(eps_lat) = 0
# At P0 there is no stress update yet; this documents the contract the return map will satisfy.
# ---------------------------------------------------------------------------
def lateral_residual(sigma_lat, eps_lat, mode="free", p=0.0, hoop_law=None):
    if mode == "free":
        return sigma_lat
    if mode == "active":
        return sigma_lat + p
    if mode == "passive":
        if hoop_law is None:
            raise ValueError("passive confinement needs hoop_law(eps_lat) -> sigma_hoop")
        return sigma_lat + hoop_law(eps_lat)
    raise ValueError(f"unknown confinement mode {mode!r}")


# ---------------------------------------------------------------------------
# P0 GATE
# ---------------------------------------------------------------------------
def run_p0_gate(fc=30.0, ft=3.0, target_fcc_ratio=1.16, verbose=True):
    results = {}

    # derive e from Kupfer (ADR 4.1 option b), then m0 follows
    e = eccentricity_from_kupfer(fc, ft, target_fcc_ratio)
    m0 = m0_of(fc, ft, e)
    results["e"] = e
    results["m0"] = m0

    # G1: uniaxial compression on surface (must be machine-exact, independent of m0/e)
    f_uc = yield_f(principal_to_voigt(-fc, 0.0, 0.0), fc, ft, e)
    results["G1_uniaxial_compression_f"] = f_uc

    # G2: uniaxial tension on surface (must be machine-exact; this is what fixes m0)
    f_ut = yield_f(principal_to_voigt(ft, 0.0, 0.0), fc, ft, e)
    results["G2_uniaxial_tension_f"] = f_ut

    # G3: ECCENTRICITY identity of the Lode function (exact, machine-precision).
    # e is the DEVIATORIC out-of-roundness: r(theta=0)=1/e, r(theta=pi/3)=1 => r_c/r_t = e.
    # (NOT the full-surface meridian ratio rho_t/rho_c, which the quadratic rho^2 term shifts
    #  away from e and which is reported below as physical output.)
    r_t = lode_r(0.0, e)
    r_c = lode_r(np.pi / 3.0, e)
    g3_err = max(abs(r_t * e - 1.0), abs(r_c - 1.0))
    results["G3_r_tensile"] = r_t
    results["G3_r_compressive"] = r_c
    results["G3_ecc_identity_err"] = g3_err
    # physical surface meridian ratios (informational): should rise toward 1 with confinement
    surf_mr = {x: meridian_ratio(fc, ft, e, xi_over_fc=x) for x in (-1.0, -3.0, -5.0)}
    results["surface_meridian_ratio"] = surf_mr

    # G4: Kupfer equibiaxial fcc/fc reproduces target
    fcc = equibiaxial_strength(fc, ft, e)
    results["G4_fcc_over_fc"] = fcc / fc

    ok = (abs(f_uc) < TOL_ONSURF and abs(f_ut) < TOL_ONSURF
          and g3_err < TOL_RATIO
          and abs(fcc / fc - target_fcc_ratio) < 1.0e-10   # bisection-tight, not the loose 1e-6
          and 0.50 < e < 0.55)                             # canonical concrete band (independent fact)
    results["PASS"] = bool(ok)

    if verbose:
        smr = "  ".join(f"xi={x:+.0f}fc:{v:.3f}" for x, v in surf_mr.items())
        print(f"  fc={fc}  ft={ft}  ft/fc={ft/fc:.3f}  target fcc/fc={target_fcc_ratio}")
        print(f"  e (from Kupfer)              = {e:.6f}   (canonical concrete ~0.52, convex>0.5)")
        print(f"  m0                           = {m0:.6f}")
        print(f"  G1 |f| uniaxial compression  = {abs(f_uc):.3e}   (tol {TOL_ONSURF:.0e})")
        print(f"  G2 |f| uniaxial tension      = {abs(f_ut):.3e}   (tol {TOL_ONSURF:.0e})")
        print(f"  G3 ecc identity r(0)e=1,r(60)=1 err = {g3_err:.3e}   (tol {TOL_RATIO:.0e})")
        print(f"     surface meridian rho_t/rho_c: {smr}   (rises to 1 with confinement)")
        print(f"  G4 equibiaxial fcc/fc        = {fcc/fc:.6f}  (target {target_fcc_ratio})")
        print(f"  => P0 GATE {'PASS' if ok else 'FAIL'}")
    return results


# ===========================================================================
# P1 — semi-implicit Menetrey-Willam return map (perfect plasticity, qh=1).
#
# Return in INVARIANT (xi, rho) space with the Lode angle theta FROZEN at the trial value
# (the ADR semi-implicit scheme): the deviatoric return is radial, so it sidesteps the
# d(theta)/d(sigma) Lode-corner singularity. Unknowns (xi, rho, dlam); 3x3 Newton on:
#     xi  = xi_tr  - 3K * dlam * m_v        (m_v = dg/dxi  = Df * m0/(sqrt3 fc), dilatancy-scaled)
#     rho = rho_tr - 2G * dlam * m_s        (m_s = dg/drho = 3 rho/fc^2 + m0 r/(sqrt6 fc))
#     f(xi, rho; theta_tr) = 0
# NON-ASSOCIATED flow: Df in (0,1] scales the volumetric (dilatant) flow; Df=1 => associated.
# Apex: if the deviatoric return would drive rho<=0, return to the hydrostatic-tension apex
#       xi_apex = sqrt3 fc / m0 (where f(xi,0)=0), rho=0.
# Hardening qh1/qh2 + ductility measure x(sigma) = the NEXT P1 increment (this slice gives the
# PEAK-STRENGTH envelope, which is the failure surface and the headline confined-triaxial gate).
# ===========================================================================
def make_material(E, nu, fc, ft, Df=1.0, target_fcc_ratio=1.16, e=None):
    if e is None:
        e = eccentricity_from_kupfer(fc, ft, target_fcc_ratio)
    K = E / (3.0 * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    return dict(E=E, nu=nu, fc=fc, ft=ft, e=e, m0=m0_of(fc, ft, e), Df=Df, K=K, G=G)


def _yf_inv(xi, rho, r, mp):
    """MW yield value from invariants directly (theta enters via the frozen r)."""
    return (SQRT1_5 * rho / mp["fc"]) ** 2 \
        + mp["m0"] * (rho * r / (SQRT6 * mp["fc"]) + xi / (SQRT3 * mp["fc"])) - 1.0


def return_map_principal(sig_tr, mp, tol=1.0e-12):
    """sig_tr = 3 trial principal stresses.
    Returns (sig_new[3], plastic, f_after, dlam, converged)."""
    fc, m0, Df, K, G = mp["fc"], mp["m0"], mp["Df"], mp["K"], mp["G"]
    sv = np.array([sig_tr[0], sig_tr[1], sig_tr[2], 0.0, 0.0, 0.0])
    xi_tr, rho_tr, th_tr, *_ = invariants(sv)
    f_tr = yield_f(sv, fc, mp["ft"], mp["e"])
    if f_tr <= tol * fc:
        return np.array(sig_tr, float), False, f_tr, 0.0, True
    r = lode_r(th_tr, mp["e"])
    m_v = Df * m0 / (SQRT3 * fc)            # dg/dxi (constant, dilatancy-scaled)
    p_tr = xi_tr / SQRT3
    s_tr = np.array(sig_tr, float) - p_tr   # trial deviator principals

    # smooth semi-implicit Newton in (xi, rho, dlam); theta frozen at th_tr
    xi, rho, dlam = xi_tr, rho_tr, 0.0
    apex = False
    converged = False
    for _ in range(100):
        m_s = 3.0 * rho / (fc * fc) + m0 * r / (SQRT6 * fc)
        R1 = xi - xi_tr + 3.0 * K * dlam * m_v
        R2 = rho - rho_tr + 2.0 * G * dlam * m_s
        R3 = _yf_inv(xi, rho, r, mp)
        if abs(R1) + abs(R2) + abs(R3) < tol * (fc + 1.0):
            converged = True
            break
        J = np.array([
            [1.0, 0.0, 3.0 * K * m_v],
            [0.0, 1.0 + 2.0 * G * dlam * (3.0 / (fc * fc)), 2.0 * G * m_s],
            [m0 / (SQRT3 * fc), 3.0 * rho / (fc * fc) + m0 * r / (SQRT6 * fc), 0.0],
        ])
        dx = np.linalg.solve(J, -np.array([R1, R2, R3]))
        xi, rho, dlam = xi + dx[0], rho + dx[1], dlam + dx[2]
        if rho < 0.0:
            apex = True
            break

    if apex:
        # Hydrostatic-tension APEX return: the deviatoric flow would drive rho<0, so project to
        # the cone vertex (rho=0, xi=sqrt3 fc/m0 where f(xi,0)=0). dlam is reconciled to the
        # VOLUMETRIC-only flow that reaches the vertex (NOT the stale deviatoric iterate):
        # dlam = (xi_tr - xi_apex)/(3K m_v). The rigorous non-unique vertex multiplier
        # (Koiter multi-surface, ADR 4.2) is deferred to the dedicated vertex sub-algorithm; the
        # trigger here is a transient-overshoot heuristic (a-priori cone test = P2+). See LEDGER_quirks.
        xi = SQRT3 * fc / m0
        p_new = xi / SQRT3
        sig_new = np.array([p_new, p_new, p_new])
        f_apex = _yf_inv(xi, 0.0, r, mp)
        dlam_apex = (xi_tr - xi) / (3.0 * K * m_v) if m_v > 0.0 else 0.0
        return sig_new, True, f_apex, dlam_apex, abs(f_apex) < tol * (fc + 1.0)

    p_new = xi / SQRT3
    dev_scale = rho / rho_tr if rho_tr > 0.0 else 0.0
    sig_new = s_tr * dev_scale + p_new
    f_after = yield_f(np.array([sig_new[0], sig_new[1], sig_new[2], 0, 0, 0]),
                      fc, mp["ft"], mp["e"])
    return sig_new, True, f_after, dlam, converged


def _elastic_pred(sig_n, deps, mp):
    lam = mp["K"] - 2.0 * mp["G"] / 3.0
    tr = deps.sum()
    return sig_n + lam * tr + 2.0 * mp["G"] * deps


def drive_axial(mp, eps11_path, confine="free", sigma3=0.0, hoop_k=0.0):
    """Strain-driven principal-stress integrator with lateral mixed control (sym eps22=eps33).
    confine: 'free' (sigma_lat=0), 'active' (sigma_lat=-sigma3), 'passive' (sigma_lat=-hoop_k*eps_lat).
    Returns dict of arrays."""
    def _hoop():
        return (lambda x: hoop_k * x) if confine == "passive" else None
    eps = np.zeros(3)
    sig = np.zeros(3)
    el = 0.0
    out = {k: [] for k in ("eps11", "sig11", "eps_lat", "sig_lat", "f", "epl_lat", "plastic", "conv")}
    for e11 in eps11_path:
        # inner Newton on lateral strain el so the lateral residual vanishes (FD Jacobian, |J| floor)
        inner_conv = False
        for _ in range(80):
            deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
            snew, _, _, _, _ = return_map_principal(_elastic_pred(sig, deps, mp), mp)
            slat = 0.5 * (snew[1] + snew[2])
            res = lateral_residual(slat, el, confine, p=sigma3, hoop_law=_hoop())
            if abs(res) < 1.0e-10 * (mp["fc"] + 1.0):
                inner_conv = True
                break
            d = 1.0e-8 * (abs(el) + 1.0e-6)
            deps2 = np.array([e11 - eps[0], (el + d) - eps[1], (el + d) - eps[2]])
            snew2, _, _, _, _ = return_map_principal(_elastic_pred(sig, deps2, mp), mp)
            slat2 = 0.5 * (snew2[1] + snew2[2])
            res2 = lateral_residual(slat2, el + d, confine, p=sigma3, hoop_law=_hoop())
            Jd = (res2 - res) / d
            if abs(Jd) < 1.0e-12:                       # |J| floor through the elastic/plastic kink
                Jd = 1.0e-12 if Jd >= 0.0 else -1.0e-12
            el -= res / Jd
        deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
        snew, plastic, f_after, _, map_conv = return_map_principal(_elastic_pred(sig, deps, mp), mp)
        sig = snew
        eps = np.array([e11, el, el])
        # lateral plastic strain PROXY (ordering-only: conflates volumetric + deviatoric parts)
        epl_lat = el - (sig[1] - mp["nu"] * (sig[0] + sig[2])) / mp["E"]
        for k, v in (("eps11", e11), ("sig11", sig[0]), ("eps_lat", el),
                     ("sig_lat", 0.5 * (sig[1] + sig[2])), ("f", f_after),
                     ("epl_lat", epl_lat), ("plastic", plastic),
                     ("conv", bool(inner_conv and map_conv))):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def confined_strength_analytic(mp, sigma3):
    """fcc: axial magnitude where sigma=(-s,-sigma3,-sigma3) crosses f=0 (oracle-of-oracle)."""
    def f_of_s(s):
        return yield_f(principal_to_voigt(-s, -sigma3, -sigma3), mp["fc"], mp["ft"], mp["e"])
    lo, hi = sigma3, sigma3 + 5.0 * mp["fc"]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f_of_s(lo) * f_of_s(mid) <= 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < 1.0e-12 * mp["fc"]:
            break
    return 0.5 * (lo + hi)


def run_p1_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, verbose=True):
    mp = make_material(E, nu, fc, ft, Df=1.0)
    res = {}

    # H1: return-to-surface — f recomputed INDEPENDENTLY from the recorded principal stresses
    # (NOT the map's self-reported f_after, which is tautological for the reconstruction step).
    # Never OUTSIDE the surface (signed max f<=tol; elastic states sit INSIDE with f<0); plastic
    # steps land ON it (|f|~0); every step's Newton CONVERGED.
    d = drive_axial(mp, np.linspace(0, -0.01, 400), confine="active", sigma3=0.1 * fc)
    pmask = d["plastic"].astype(bool)
    f_indep = np.array([yield_f(principal_to_voigt(d["sig11"][i], d["sig_lat"][i], d["sig_lat"][i]),
                                fc, ft, mp["e"]) for i in range(len(d["sig11"]))])
    res["H1_max_f_signed"] = float(np.max(f_indep))                  # must be <= tol (no overshoot)
    res["H1_max_f_on_plastic"] = float(np.max(np.abs(f_indep[pmask]))) if pmask.any() else 0.0
    res["H1_all_converged"] = bool(d["conv"].all())
    # also exercise a COARSE-step confined case to stress the elastic/plastic-kink inner Newton
    dc = drive_axial(mp, np.linspace(0, -0.012, 8), confine="active", sigma3=0.1 * fc)
    res["H1_coarse_converged"] = bool(dc["conv"].all())

    # H2: uniaxial compression strength -> fc
    duc = drive_axial(mp, np.linspace(0, -0.006, 300), confine="free")
    fc_sim = -np.min(duc["sig11"])
    res["H2_uniaxial_fc"] = fc_sim
    res["H2_err"] = abs(fc_sim / fc - 1.0)

    # H3: uniaxial tension strength -> ft
    dut = drive_axial(mp, np.linspace(0, 0.0008, 300), confine="free")
    ft_sim = np.max(dut["sig11"])
    res["H3_uniaxial_ft"] = ft_sim
    res["H3_err"] = abs(ft_sim / ft - 1.0)

    # H4: confined triaxial fcc(sigma3) — CONSISTENCY gate: the driver return map lands on the
    # SAME surface the analytic root-finder uses (both share yield_f/invariants; this is NOT an
    # independent oracle — see the Richart band test for the independent physical check), monotone.
    fcc_sim, fcc_ana = {}, {}
    for p in (0.0, 0.05 * fc, 0.10 * fc, 0.20 * fc):
        dd = drive_axial(mp, np.linspace(0, -0.012, 500), confine="active", sigma3=p)
        fcc_sim[p] = -np.min(dd["sig11"])
        fcc_ana[p] = confined_strength_analytic(mp, p)
    res["H4_fcc_sim"] = fcc_sim
    res["H4_fcc_ana"] = fcc_ana
    res["H4_max_rel_err"] = max(abs(fcc_sim[p] / fcc_ana[p] - 1.0) for p in fcc_sim)
    ps = sorted(fcc_sim)
    res["H4_monotone"] = all(fcc_sim[ps[i]] < fcc_sim[ps[i + 1]] for i in range(len(ps) - 1))

    # H6: APEX / hydrostatic-tension return (exercises the apex branch the drivers never reach)
    xi_apex = SQRT3 * fc / mp["m0"]
    p_apex = xi_apex / SQRT3
    sa, _, fa, _, conv_a = return_map_principal([10.0 * ft, 10.0 * ft, 10.0 * ft], mp)
    res["H6_apex_p"] = float(sa[0])
    res["H6_apex_p_target"] = float(p_apex)
    res["H6_apex_f"] = float(abs(fa))
    res["H6_apex_ok"] = bool(abs(sa[0] - p_apex) < 1.0e-9 * fc and abs(fa) < 1.0e-9 * fc
                             and np.allclose(sa, sa[0]) and conv_a)

    # H5: dilatancy knob — Df<1 reduces lateral plastic dilation vs Df=1 (associated)
    mp_assoc = make_material(E, nu, fc, ft, Df=1.0)
    mp_low = make_material(E, nu, fc, ft, Df=0.3)
    da = drive_axial(mp_assoc, np.linspace(0, -0.006, 300), confine="free")
    dl = drive_axial(mp_low, np.linspace(0, -0.006, 300), confine="free")
    res["H5_lat_plastic_assoc"] = float(da["epl_lat"][-1])
    res["H5_lat_plastic_low"] = float(dl["epl_lat"][-1])
    res["H5_dilatancy_reduced"] = abs(dl["epl_lat"][-1]) < abs(da["epl_lat"][-1])

    ok = (res["H1_max_f_signed"] < 1.0e-7 and res["H1_max_f_on_plastic"] < 1.0e-7
          and res["H1_all_converged"] and res["H1_coarse_converged"]
          and res["H2_err"] < 1.0e-8 and res["H3_err"] < 1.0e-8
          and res["H4_max_rel_err"] < 1.0e-8 and res["H4_monotone"]
          and res["H6_apex_ok"]
          and res["H5_dilatancy_reduced"])
    res["PASS"] = bool(ok)

    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft}  e={mp['e']:.4f} m0={mp['m0']:.4f}")
        print(f"  H1 return-to-surface (INDEPENDENT f): max f signed={res['H1_max_f_signed']:.3e},"
              f"  max|f| plastic={res['H1_max_f_on_plastic']:.3e},  converged={res['H1_all_converged']}"
              f"/{res['H1_coarse_converged']}(coarse)")
        print(f"  H2 uniaxial fc_sim = {fc_sim:.4f}  (target {fc})  err={res['H2_err']:.2e}")
        print(f"  H3 uniaxial ft_sim = {ft_sim:.4f}  (target {ft})  err={res['H3_err']:.2e}")
        print("  H4 confined fcc(sigma3): driver vs analytic surface crossing —")
        for p in sorted(fcc_sim):
            print(f"       sigma3={p:6.2f} (p/fc={p/fc:.2f}):  sim={fcc_sim[p]:7.3f}"
                  f"  ana={fcc_ana[p]:7.3f}  fcc/fc={fcc_sim[p]/fc:.3f}")
        print(f"       max rel err={res['H4_max_rel_err']:.2e} (consistency)  monotone={res['H4_monotone']}")
        print(f"  H6 apex (hydrostatic tension): p={res['H6_apex_p']:.4f} target={res['H6_apex_p_target']:.4f}"
              f"  |f|={res['H6_apex_f']:.2e}  ok={res['H6_apex_ok']}")
        print(f"  H5 dilatancy (ordering-only): lat plastic Df=1.0 -> {res['H5_lat_plastic_assoc']:+.3e},"
              f"  Df=0.3 -> {res['H5_lat_plastic_low']:+.3e}  reduced={res['H5_dilatancy_reduced']}")
        print(f"  => P1 GATE {'PASS' if ok else 'FAIL'}")
    return res


if __name__ == "__main__":
    print("=" * 74)
    print("LadrunoConcrete3D P0 oracle gate — Menetrey-Willam surface normalization")
    print("=" * 74)
    all_ok = True
    for (fc, ft) in [(30.0, 3.0), (40.0, 3.5), (50.0, 4.0), (25.0, 2.0)]:
        r = run_p0_gate(fc, ft, verbose=True)
        all_ok = all_ok and r["PASS"]
        print("-" * 74)
    print(f"P0 ALL CASES: {'PASS' if all_ok else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P1 oracle gate — semi-implicit MW return map (perfect plasticity)")
    print("=" * 74)
    p1 = run_p1_gate(verbose=True)
    print("-" * 74)
    print(f"P1: {'PASS' if p1['PASS'] else 'FAIL'}")
