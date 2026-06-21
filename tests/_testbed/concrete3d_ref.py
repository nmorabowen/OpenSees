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


def yield_f(sig, fc, ft, e, qh1v=1.0, qh2v=1.0):
    """CDPM2 yield function f_p — Grassl et al. 2013 IJSS Eq. (18). qh1v=qh2v=1 => failure surface
    Eq. (21) = Menetrey-Willam. COMPRESSION NEGATIVE. Note xi/(sqrt3 fc) = (I1/3)/fc = sigV/fc."""
    xi, rho, theta, *_ = invariants(sig)
    m0 = m0_of(fc, ft, e)
    r = lode_r(theta, e)
    sigV_fc = xi / (SQRT3 * fc)                       # = sigma_V/fc  (Eq.12: sigV=I1/3)
    AV = rho / (SQRT6 * fc) + sigV_fc                 # the ellipsoidal hardening-cap base
    RR = rho * r / (SQRT6 * fc) + sigV_fc             # the m0-friction bracket (carries Lode r)
    quad = SQRT1_5 * rho / fc                         # sqrt(3/2) rho/fc
    # Eq. (18): {[1-qh1]AV^2 + quad}^2 + m0 qh1^2 qh2 RR - qh1^2 qh2^2
    return ((1.0 - qh1v) * AV * AV + quad) ** 2 \
        + m0 * qh1v * qh1v * qh2v * RR - (qh1v * qh1v) * (qh2v * qh2v)


# ---------------------------------------------------------------------------
# CDPM2 plasticity HARDENING building blocks (Grassl et al. 2013, Sec. 2.2.3-2.2.4).
# Unit-verified here; wired into the hardening return map in the next increment.
#   qh1(kp)  Eq. (30)   pre-peak size: qh0 -> 1 over kp in [0,1], slope Hp at kp=1^-
#   qh2(kp)  Eq. (31)   post-peak: 1 (kp<1), then 1+Hp(kp-1)
#   xh(sigV) Eqs. (33-36)  hardening DUCTILITY measure; Rh=-sigV/fc-1/3 (Eq.34) => more ductile
#            under compression (sigV<0 => Rh>0 => larger xh => slower kp => more plastic strain to peak)
# ---------------------------------------------------------------------------
def qh1(kp, qh0, Hp):
    if kp < 1.0:
        return qh0 + (1.0 - qh0) * (kp**3 - 3.0 * kp**2 + 3.0 * kp) \
            - Hp * (kp**3 - 3.0 * kp**2 + 2.0 * kp)
    return 1.0


def qh2(kp, Hp):
    return 1.0 if kp < 1.0 else 1.0 + Hp * (kp - 1.0)


def ductility_xh(sigV, fc, Ah, Bh, Ch, Dh):
    Rh = -sigV / fc - 1.0 / 3.0                       # Eq. (34)
    if Rh >= 0.0:
        return Ah - (Ah - Bh) * np.exp(-Rh / Ch)      # Eq. (33) upper branch
    Eh = Bh - Dh                                      # Eq. (35)
    Fh = (Bh - Dh) * Ch / (Ah - Bh)                   # Eq. (36)
    return Eh * np.exp(Rh / Fh) + Dh                  # Eq. (33) lower branch


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
def make_material(E, nu, fc, ft, Df=1.0, target_fcc_ratio=1.16, e=None,
                  qh0=0.3, Hp=0.5, Ah=0.08, Bh=0.003, Ch=2.0, Dh=1.0e-6, eta=0.0,
                  ct_temper="none"):
    # qh0,Hp: hardening laws Eq.30-31.  Ah,Bh,Ch,Dh: ductility measure Eq.33 (literature defaults;
    # calibrated per-concrete from peak strains — flagged in ADR 6 as recalibrate-for-fork-data).
    # eta: Duvaut-Lions viscoplastic relaxation time (ADR 4.4). eta=0 => inviscid, BYTE-identical to
    #   the Tier-1 path. eta>0 relaxes the plastic return toward the inviscid solution with the factor
    #   beta = dt/(eta+dt); see `damaged_step_tensor`. Units of TIME (consistent with the analysis dt).
    if e is None:
        e = eccentricity_from_kupfer(fc, ft, target_fcc_ratio)
    K = E / (3.0 * (1.0 - 2.0 * nu))
    G = E / (2.0 * (1.0 + nu))
    return dict(E=E, nu=nu, fc=fc, ft=ft, e=e, m0=m0_of(fc, ft, e), Df=Df, K=K, G=G,
                qh0=qh0, Hp=Hp, Ah=Ah, Bh=Bh, Ch=Ch, Dh=Dh, eta=eta, ct_temper=ct_temper)


# ---------------------------------------------------------------------------
# P2h — compression->tension damage-coupling TEMPER (the `-ctTemper` modes). Per literal CDPM2 (Eq.43,
# kappa_dt-dot = eps_tilde-dot, NO (1-alpha_c) factor) the tensile damage history accumulates during
# compression too, so a compression excursion PRE-DAMAGES a subsequent tension reload to ~0 (the DT5
# diagnostic). `ct_temper` scales the tensile-history (kdt1/kdt2) accumulation by a tensile weight w_t:
#   'none'   : w_t = 1            -> literal CDPM2 (default; BYTE-identical to the un-tempered drivers)
#   'alphat' : w_t = 1 - alpha_c  -> stress-state weight; compression (alpha_c~1) => w_t~0, no tension
#              pre-damage. Leaves BOTH monotonic backbones EXACT (alpha_c=0 in tension => w_t=1; the
#              compression backbone rides the kdc channel) and removes ONLY the cross-coupling.
#   'proj'   : w_t = ||proj of d eps_p onto TENSILE-stress directions|| / ||d eps_p||  -> the fraction of
#              the plastic-strain increment that acts along POSITIVE effective-stress principal directions.
#              In compression ALL principals are compressive => the projection is empty => w_t=0 (full
#              shield). NB the plastic strain's OWN positive part is NOT a valid shield: compression's
#              dilatant flow makes the lateral plastic strains positive, so ||<deps_p>+||/||deps_p|| stays
#              large in compression. Projecting onto the tensile-STRESS frame is the correct mechanistic
#              measure. 'proj' lightly softens the monotonic TENSION backbone (the loaded axial direction
#              carries < 100% of ||deps_p|| because of the lateral plastic flow).
# Both temper modes -> 0 in pure compression and -> 1 (alphat) / <1 (proj) in pure tension.
# ---------------------------------------------------------------------------
def tensile_damage_weight(mp, ac, depl6, w_stress, V_stress):
    """Tensile-damage-history weight w_t for the -ctTemper modes. `ac`=alpha_c (Eq.46); `depl6`=the
    plastic-strain increment (Voigt tensor); `w_stress`,`V_stress`=the EFFECTIVE-stress eigenvalues and
    eigenvectors (proj projects depl onto the directions where w_stress>0)."""
    mode = mp.get("ct_temper", "none")
    if mode == "alphat":
        w = 1.0 - ac
        return w if w > 0.0 else 0.0
    if mode == "proj":
        M = voigt_to_mat(depl6)
        nrm = float(np.sqrt(np.sum(M * M)))
        if nrm <= 1.0e-300:
            return 1.0
        Deig = V_stress.T @ M @ V_stress                 # plastic-strain increment in the stress eigenframe
        floor = 1.0e-6 * mp["ft"]
        tens = sum(Deig[a, a] ** 2 for a in range(3) if w_stress[a] > floor)   # tensile-stress directions
        return float(np.sqrt(tens) / nrm)
    return 1.0   # 'none' (literal CDPM2)


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

    # smooth semi-implicit Newton in (xi, rho, dlam); theta frozen at th_tr.
    # FLOW gradient m_s = dg/drho is Lode-INDEPENDENT (CDPM2 potential has no r; Eq.22-25) — NOT
    # df/drho. The YIELD gradient (Jacobian row 3) keeps r. Volumetric flow m_v carries dilatancy Df.
    xi, rho, dlam = xi_tr, rho_tr, 0.0
    apex = False
    converged = False
    for _ in range(100):
        m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc)          # dg/drho (flow, NO Lode r)
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


# ===========================================================================
# P1-hardening — CDPM2 plasticity with hardening (qh1/qh2/kappa_p/ductility).
#
# Semi-implicit (theta frozen) return map with FOUR unknowns (xi, rho, dlam, kappa_p):
#     xi  = xi_tr  - 3K dlam m_v                          (m_v = Df m0/(sqrt3 fc), dilatancy)
#     rho = rho_tr - 2G dlam m_s                          (m_s = 3 rho/fc^2 + m0/(sqrt6 fc), NO Lode r)
#     f_p(xi, rho, theta_tr, kappa_p) = 0                 (Eq.18 with qh1(kp), qh2(kp))
#     kappa_p = kappa_p,n + dlam (||m||/xh(sigV)) (2 cos theta)^2   (Eq.32)
# with ||m|| = sqrt(m_v^2 + m_s^2). Jacobian = NUMERICAL (oracle is a correctness reference; the
# analytic 4x4 tangent is the C++ kernel's job, FD-checked against this). FLOW kept simplified
# (constant-Df dilatancy, not the full mg potential Eq.22-29) — documented; mg = follow-on.
# Reduce-to-P1 (qh0=1, Hp=0 => qh1=qh2=1 => Eq.21) is the safety gate.
# ===========================================================================
def _yf_inv_hard(xi, rho, r, kp, mp):
    """CDPM2 Eq.18 in invariant space with hardening qh1(kp), qh2(kp). _yf_inv == this at kp s.t.
    qh1=qh2=1 (the failure surface Eq.21)."""
    fc, m0 = mp["fc"], mp["m0"]
    q1 = qh1(kp, mp["qh0"], mp["Hp"])
    q2 = qh2(kp, mp["Hp"])
    sigV_fc = xi / (SQRT3 * fc)
    AV = rho / (SQRT6 * fc) + sigV_fc
    RR = rho * r / (SQRT6 * fc) + sigV_fc
    quad = SQRT1_5 * rho / fc
    return ((1.0 - q1) * AV * AV + quad) ** 2 + m0 * q1 * q1 * q2 * RR - (q1 * q1) * (q2 * q2)


def return_map_hardening(sig_tr, mp, kp_n, tol=1.0e-11):
    """sig_tr = 3 trial principal stresses, kp_n = committed kappa_p.
    Returns (sig_new[3], kp_new, plastic, f_after, converged)."""
    fc, m0, Df, K, G = mp["fc"], mp["m0"], mp["Df"], mp["K"], mp["G"]
    sv = np.array([sig_tr[0], sig_tr[1], sig_tr[2], 0.0, 0.0, 0.0])
    xi_tr, rho_tr, th_tr, *_ = invariants(sv)
    r = lode_r(th_tr, mp["e"])
    f_tr = _yf_inv_hard(xi_tr, rho_tr, r, kp_n, mp)
    if f_tr <= tol * fc:
        return np.array(sig_tr, float), kp_n, False, f_tr, True
    p_tr = xi_tr / SQRT3
    s_tr = np.array(sig_tr, float) - p_tr
    cos2 = (2.0 * np.cos(th_tr)) ** 2
    m_v = Df * m0 / (SQRT3 * fc)

    def resid(u):
        xi, rho, dlam, kp = u
        m_s = 3.0 * rho / (fc * fc) + m0 / (SQRT6 * fc)
        mnorm = np.sqrt(m_v * m_v + m_s * m_s)
        xh = ductility_xh(xi / SQRT3, fc, mp["Ah"], mp["Bh"], mp["Ch"], mp["Dh"])
        return np.array([
            xi - xi_tr + 3.0 * K * dlam * m_v,
            rho - rho_tr + 2.0 * G * dlam * m_s,
            _yf_inv_hard(xi, rho, r, kp, mp),
            kp - kp_n - dlam * mnorm / xh * cos2,
        ])

    u = np.array([xi_tr, rho_tr, 0.0, kp_n])
    converged = False
    apex = False
    for _ in range(100):
        R = resid(u)
        if (abs(R[0]) < tol * fc and abs(R[1]) < tol * fc
                and abs(R[2]) < tol and abs(R[3]) < tol):
            converged = True
            break
        J = np.zeros((4, 4))
        for j in range(4):
            du = 1.0e-8 * (abs(u[j]) + 1.0e-6)
            up = u.copy()
            up[j] += du
            J[:, j] = (resid(up) - R) / du
        u = u + np.linalg.solve(J, -R)
        if u[1] < 0.0:                       # deviatoric return overshoots => hydrostatic apex
            apex = True
            break

    xi, rho, dlam, kp = u
    if apex:
        # Hydrostatic-tension APEX return (mirror return_map_principal): project to the cone vertex
        # rho=0, f(xi,0;kp)=0. For the hardened surface qh2=1 (kp<1) this is xi_apex=sqrt3 fc/m0;
        # solve the 1-D f(xi,0)=0 generally for robustness. (Rigorous non-unique vertex multiplier =
        # Koiter, ADR 4.2 — deferred.)
        lo, hi = 0.0, SQRT3 * fc / m0 * 4.0
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if _yf_inv_hard(lo, 0.0, r, kp, mp) * _yf_inv_hard(mid, 0.0, r, kp, mp) <= 0.0:
                hi = mid
            else:
                lo = mid
        xi = 0.5 * (lo + hi)
        p_new = xi / SQRT3
        sig_new = np.array([p_new, p_new, p_new])
    else:
        p_new = xi / SQRT3
        dev_scale = rho / rho_tr if rho_tr > 0.0 else 0.0
        sig_new = s_tr * dev_scale + p_new

    # HONEST convergence: recompute f INDEPENDENTLY at the returned stress with its OWN (updated)
    # Lode angle — the frozen-theta residual R[2] can read ~0 while the point is off-surface near
    # the apex (the frozen theta != the returned stress's true theta). Never report converged for
    # an off-surface point.
    f_indep = yield_f(np.array([sig_new[0], sig_new[1], sig_new[2], 0.0, 0.0, 0.0]),
                      fc, mp["ft"], mp["e"], qh1(kp, mp["qh0"], mp["Hp"]), qh2(kp, mp["Hp"]))
    converged = bool((converged or apex) and abs(f_indep) < 1.0e-7 * (fc + 1.0))
    return sig_new, kp, True, f_indep, converged


def drive_hardening(mp, eps11_path, confine="free", sigma3=0.0):
    """Hardening counterpart of drive_axial (tracks kappa_p). Principal-stress, lateral mixed control."""
    eps = np.zeros(3)
    sig = np.zeros(3)
    el = 0.0
    kp = 0.0
    out = {k: [] for k in ("eps11", "sig11", "kp", "plastic", "conv")}
    for e11 in eps11_path:
        inner_conv = False
        for _ in range(80):
            deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
            snew, _, _, _, _ = return_map_hardening(_elastic_pred(sig, deps, mp), mp, kp)
            res = lateral_residual(0.5 * (snew[1] + snew[2]), el, confine, p=sigma3)
            if abs(res) < 1.0e-10 * (mp["fc"] + 1.0):
                inner_conv = True
                break
            d = 1.0e-8 * (abs(el) + 1.0e-6)
            deps2 = np.array([e11 - eps[0], (el + d) - eps[1], (el + d) - eps[2]])
            snew2, _, _, _, _ = return_map_hardening(_elastic_pred(sig, deps2, mp), mp, kp)
            res2 = lateral_residual(0.5 * (snew2[1] + snew2[2]), el + d, confine, p=sigma3)
            Jd = (res2 - res) / d
            if abs(Jd) < 1.0e-12:
                Jd = 1.0e-12 if Jd >= 0.0 else -1.0e-12
            el -= res / Jd
        deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
        snew, kp_new, plastic, _, map_conv = return_map_hardening(_elastic_pred(sig, deps, mp), mp, kp)
        sig, kp = snew, kp_new
        eps = np.array([e11, el, el])
        for k, v in (("eps11", e11), ("sig11", sig[0]), ("kp", kp),
                     ("plastic", plastic), ("conv", bool(inner_conv and map_conv))):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def drive_triaxial_strain(mp, eps_end, nsteps):
    """Fully strain-controlled OFF-MERIDIAN hardening path (3 DISTINCT principal strains, no lateral
    solve) — the axisymmetric drive_hardening cannot generate off-meridian (theta != 60deg) states.
    Returns (max independent-f drift over the path, final kappa_p, max |theta-theta_frozen| proxy)."""
    eps = np.zeros(3)
    sig = np.zeros(3)
    kp = 0.0
    max_f = 0.0
    for i in range(1, nsteps + 1):
        et = np.array(eps_end, float) * (i / nsteps)
        sig, kp, _, f_indep, _ = return_map_hardening(_elastic_pred(sig, et - eps, mp), mp, kp)
        eps = et
        max_f = max(max_f, abs(f_indep))          # f_indep is the INDEPENDENT-theta drift (post-M1)
    return max_f, kp


def run_hardening_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, verbose=True):
    res = {}
    # U1: hardening building-block unit identities (Eqs.30-36)
    q1_0 = qh1(0.0, 0.3, 0.5)
    q1_1 = qh1(1.0, 0.3, 0.5)
    ks = np.linspace(0.0, 1.0, 60)
    q1_mono = all(qh1(ks[i], 0.3, 0.5) <= qh1(ks[i + 1], 0.3, 0.5) + 1e-12 for i in range(len(ks) - 1))
    xh_comp = ductility_xh(-fc, fc, 0.08, 0.003, 2.0, 1.0e-6)        # compression (Rh>0)
    xh_tens = ductility_xh(+0.2 * fc, fc, 0.08, 0.003, 2.0, 1.0e-6)  # tension (Rh<0)
    res["U1_ok"] = bool(abs(q1_0 - 0.3) < 1e-14 and abs(q1_1 - 1.0) < 1e-14 and q1_mono
                        and qh2(0.5, 0.5) == 1.0 and qh2(1.5, 0.5) > 1.0 and xh_comp > xh_tens)

    # HA: reduce-to-P1 — qh0=1, Hp=0 => qh1=qh2=1 => hardening map == perfect-plastic map.
    # (compression/deviatoric states only; hydrostatic-tension apex+hardening is deferred to P2.)
    mp_rigid = make_material(E, nu, fc, ft, qh0=1.0, Hp=0.0)
    maxd = 0.0
    for st in ([-40, -3, -3], [-50, 0, 0], [-30, -10, -8], [-35, -5, -2], [3.5, 0, 0]):
        sh, _, _, _, _ = return_map_hardening(st, mp_rigid, 0.0)
        sp, _, _, _, _ = return_map_principal(st, mp_rigid)
        maxd = max(maxd, float(np.max(np.abs(sh - sp))))
    res["HA_reduce_to_P1"] = maxd

    # NB: effective-stress plasticity is MONOTONIC hardening (no peak — the peak/softening comes
    # from DAMAGE, P2). The failure surface is reached EXACTLY at kappa_p=1, where uniaxial
    # compression sits at sig11 = -fc. Gate on that, not on a (nonexistent) plasticity peak.
    mp = make_material(E, nu, fc, ft)
    d0 = drive_hardening(mp, np.linspace(0, -0.006, 800), confine="free")
    pl = d0["plastic"].astype(bool)
    first_yield = -d0["sig11"][pl][0] if pl.any() else 0.0
    # HB: sig11 at the step where kappa_p first crosses 1  ==  fc (failure surface reached)
    ik = int(np.argmax(d0["kp"] >= 1.0)) if (d0["kp"] >= 1.0).any() else -1
    sig_at_kp1 = -d0["sig11"][ik] if ik >= 0 else np.nan
    res["HB_sig_at_kp1"] = float(sig_at_kp1)
    res["HB_err"] = abs(sig_at_kp1 / fc - 1.0) if ik >= 0 else 1.0
    res["HB_all_conv"] = bool(d0["conv"].all())
    # HC: pre-peak nonlinearity — first yield near qh0*fc (<< fc), monotone-rising hardening to fc
    res["HC_first_yield_over_fc"] = first_yield / fc
    mono = all(-d0["sig11"][i] <= -d0["sig11"][i + 1] + 1e-9 for i in range(np.argmax(d0["kp"] >= 1.0)))
    res["HC_nonlinear"] = bool(first_yield < 0.6 * fc and mono and ik >= 0)

    # HD: confinement ductility — axial strain at kappa_p=1 grows with sigma3 (Eq.32-34 ductility)
    eps_at_kp1 = {}
    for p in (0.0, 0.10 * fc):
        dd = drive_hardening(mp, np.linspace(0, -0.02, 1000), confine="active", sigma3=p)
        j = int(np.argmax(dd["kp"] >= 1.0)) if (dd["kp"] >= 1.0).any() else -1
        eps_at_kp1[p] = abs(dd["eps11"][j]) if j >= 0 else np.nan
    res["HD_eps_at_kp1"] = eps_at_kp1
    res["HD_confinement_more_ductile"] = bool(eps_at_kp1[0.10 * fc] > eps_at_kp1[0.0])

    # HE (DIAGNOSTIC, not a pass gate): off-meridian hardening drift. The axisymmetric HB/HC/HD are
    # stuck on the compressive meridian (theta=60, frozen=true, no drift). A 3-distinct-principal
    # strain path exposes a KNOWN LIMITATION: near FIRST YIELD off-meridian (kappa_p~0, the surface
    # is small + qh1 ramps steeply) the semi-implicit return leaves an off-surface drift that the
    # honest independent-f flag (M1) correctly reports as converged=False. The build PR must address
    # it (sub-stepping near first yield / off-meridian, and the apex sub-algorithm). Recorded so the
    # handoff states it honestly; NOT gated because the fix is build-PR scope.
    eps_end = [-3.0e-3, 8.0e-4, 1.0e-4]
    res["HE_drift_offmeridian"] = drive_triaxial_strain(mp, eps_end, 256)[0]

    ok = (res["U1_ok"] and res["HA_reduce_to_P1"] < 1.0e-8
          and res["HB_err"] < 0.02 and res["HB_all_conv"]
          and res["HC_nonlinear"] and res["HD_confinement_more_ductile"])
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft}  qh0={mp['qh0']} Hp={mp['Hp']}")
        print(f"  U1 hardening unit identities (qh1/qh2/ductility) ok={res['U1_ok']}")
        print(f"  HA reduce-to-P1 (qh0=1,Hp=0): max|sig_hard-sig_pp|={res['HA_reduce_to_P1']:.2e}")
        print(f"  HB sig11 at kappa_p=1 = {res['HB_sig_at_kp1']:.4f} (target fc={fc}) err={res['HB_err']:.2e}"
              f"  conv={res['HB_all_conv']}")
        print(f"  HC pre-peak: first yield/fc={res['HC_first_yield_over_fc']:.3f} (~qh0), monotone hardening"
              f"  ok={res['HC_nonlinear']}")
        print(f"  HE off-meridian drift (DIAGNOSTIC, known limitation, not gated) = "
              f"{res['HE_drift_offmeridian']:.2e}  (honest-f flags these converged=False; build-PR fix)")
        print(f"  HD confinement ductility: eps(kp=1) unconf={eps_at_kp1[0.0]:.5f}"
              f"  conf(0.1fc)={eps_at_kp1[0.10*fc]:.5f}  more_ductile={res['HD_confinement_more_ductile']}")
        print(f"  => HARDENING GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P1-tangent — general-tensor return map (spectral recompose) + CONSISTENT (algorithmic) tangent.
#
# The principal-stress return maps above are lifted to full 6-tensor strain control via the
# spectral split: eigendecompose the trial stress, return the principal stresses (radial
# deviatoric return preserves eigenvectors), recompose with the trial eigenvectors. The consistent
# tangent dsigma/depsilon is NUMERICAL here (the algorithmic tangent of the discrete map — the
# oracle reference). The C++ kernel's ANALYTIC tangent is FD-checked against this in the build PR.
# Tensor (NOT engineering) shear convention throughout: voigt {00,11,22,01,12,02}, off-diag = eps_ij.
# ===========================================================================
def voigt_to_mat(v):
    return np.array([[v[0], v[3], v[5]], [v[3], v[1], v[4]], [v[5], v[4], v[2]]], float)


def mat_to_voigt(M):
    return np.array([M[0, 0], M[1, 1], M[2, 2], M[0, 1], M[1, 2], M[0, 2]], float)


def elastic_C(mp):
    lam = mp["K"] - 2.0 * mp["G"] / 3.0
    G = mp["G"]
    C = np.zeros((6, 6))
    for i in range(3):
        for j in range(3):
            C[i, j] = lam
        C[i, i] += 2.0 * G
    for i in range(3, 6):
        C[i, i] = 2.0 * G                      # tensor: dsig_ij = 2G deps_ij
    return C


def elastic_pred_tensor(sig_n, deps, mp):
    lam = mp["K"] - 2.0 * mp["G"] / 3.0
    G = mp["G"]
    tr = deps[0] + deps[1] + deps[2]
    s = np.array(sig_n, float).copy()
    for i in range(3):
        s[i] += lam * tr + 2.0 * G * deps[i]
    for i in range(3, 6):
        s[i] += 2.0 * G * deps[i]
    return s


def return_map_tensor(sig_n, deps, mp, kp_n, hardening=True):
    """6-tensor return: elastic predict -> eigendecompose -> principal return -> recompose.
    Returns (sig_new[6], kp_new, plastic, converged)."""
    sig_tr = elastic_pred_tensor(sig_n, deps, mp)
    w, V = np.linalg.eigh(voigt_to_mat(sig_tr))           # w ascending, V columns = eigenvectors
    if hardening:
        sp, kp_new, plastic, _, conv = return_map_hardening(w, mp, kp_n)
    else:
        sp, plastic, _, _, conv = return_map_principal(w, mp)
        kp_new = kp_n
    sig_new = mat_to_voigt(V @ np.diag(sp) @ V.T)          # radial return preserves V
    return sig_new, kp_new, plastic, conv


def consistent_tangent(sig_n, deps, mp, kp_n, hardening=True, rel_step=1.0e-6):
    """Numerical algorithmic tangent dsigma/depsilon (6x6) by central difference of the map."""
    base = mp["fc"] / mp["E"]
    C = np.zeros((6, 6))
    for j in range(6):
        d = rel_step * (abs(deps[j]) + base)
        dp = np.array(deps, float).copy(); dp[j] += d
        dm = np.array(deps, float).copy(); dm[j] -= d
        sp, _, _, _ = return_map_tensor(sig_n, dp, mp, kp_n, hardening)
        sm, _, _, _ = return_map_tensor(sig_n, dm, mp, kp_n, hardening)
        C[:, j] = (sp - sm) / (2.0 * d)
    return C


def run_tangent_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, verbose=True):
    res = {}
    mp = make_material(E, nu, fc, ft)
    Cel = elastic_C(mp)
    nrm = lambda A: float(np.sqrt(np.sum(A * A)))

    # T1: tensor return reduces to the principal map for a DIAGONAL plastic trial
    deps_diag = np.array([-2.0e-3, 5.0e-4, 7.0e-4, 0.0, 0.0, 0.0])
    st, _, _, _ = return_map_tensor(np.zeros(6), deps_diag, mp, 0.0)
    sig_tr = elastic_pred_tensor(np.zeros(6), deps_diag, mp)
    w = np.sort(sig_tr[:3])
    sp_ref, _, _, _, _ = return_map_hardening(w, mp, 0.0)
    res["T1_reduce_diag"] = float(np.max(np.abs(np.sort(st[:3]) - np.sort(sp_ref))))

    # Tangent gate runs on the PERFECT-PLASTIC map because it has an ANALYTIC inner Jacobian
    # (machine-clean inner solve). NOT because the hardening FD tangent is noisy — the outer
    # central-difference tangent of the hardening map is actually just as step-stable (~2.8e-10).
    # The hardening map's analytic 4x4 tangent is the C++ build-PR deliverable, FD-checked against
    # this machinery. hardening=False.
    HP = dict(hardening=False)

    # T2: elastic step -> tangent == elastic C
    C_el = consistent_tangent(np.zeros(6), np.array([-1.0e-5, 0, 0, 0, 0, 0]), mp, 0.0, **HP)
    res["T2_elastic_err"] = nrm(C_el - Cel) / nrm(Cel)

    # PLASTIC, off-meridian (Lode ~52deg), low-confinement (lateral expansion) so it actually
    # yields, 3 distinct principals + shear (clean eigvecs). Df=0.3 = non-associated DILATANCY.
    # Same state at (e=1,Df=1) is fully associated (T3b) => the symmetric/asymmetric pair isolates
    # exactly the non-associativity that forces an unsymmetric solver (ADR 4.4/4.5).
    deps_pl = np.array([-2.0e-3, 6.0e-4, 2.0e-4, 1.0e-4, 0.0, 0.0])
    mp_na = make_material(E, nu, fc, ft, Df=0.3)
    C_pl = consistent_tangent(np.zeros(6), deps_pl, mp_na, 0.0, **HP)
    res["T3_asym_nonassoc"] = nrm(C_pl - C_pl.T) / nrm(C_pl)        # non-associated => NON-symmetric

    # T3b: associated limit e=1 (no Lode r), Df=1 => flow == yield gradient. The asymmetry drops
    # ~20x vs the non-associated model, confirming Df+Lode drive the BULK asymmetry. The residual
    # (~2.4%) is the FROZEN-EIGENVECTOR SPECTRAL RECOMPOSE (V diag(sp) V^T with V held at the trial
    # drops the eigenprojection/spin terms dV/deps), NOT the Lode theta-freeze. FALSIFIED that it is
    # the theta-freeze (see T3c probe below): in principal space (no recompose) the same off-meridian
    # associated state is machine-symmetric (~2e-10); the 2.4% appears ONLY with shear and scales
    # LINEARLY with shear magnitude, FD-step-independent. Conclusion is unchanged: Tier-1 is non-
    # symmetric even when associated => unsymmetric solver required UNCONDITIONALLY (ADR 4.4).
    mp_as = make_material(E, nu, fc, ft, Df=1.0, e=1.0)
    C_as = consistent_tangent(np.zeros(6), deps_pl, mp_as, 0.0, **HP)
    res["T3b_sym_assoc"] = nrm(C_as - C_as.T) / nrm(C_as)

    # T3c: FALSIFICATION probe — the associated-limit asymmetry is the spectral recompose, not the
    # theta-freeze. (i) principal-space associated off-meridian state is machine-symmetric; (ii) the
    # full-tensor asymmetry scales linearly with shear and ->0 as shear->0.
    deps_noshear = np.array([-2.0e-3, 6.0e-4, 2.0e-4, 0.0, 0.0, 0.0])
    C_ns = consistent_tangent(np.zeros(6), deps_noshear, mp_as, 0.0, **HP)
    res["T3c_assoc_noshear_sym"] = nrm(C_ns - C_ns.T) / nrm(C_ns)            # ~0 (no recompose spin)
    asy = []
    for g in (1.0e-5, 5.0e-5):
        dps = np.array([-2.0e-3, 6.0e-4, 2.0e-4, g, 0.0, 0.0])
        Cg = consistent_tangent(np.zeros(6), dps, mp_as, 0.0, **HP)
        asy.append(nrm(Cg - Cg.T) / nrm(Cg))
    res["T3c_shear_linear"] = asy[1] / asy[0] if asy[0] > 0 else 0.0        # ~5 (linear in shear)

    # T4: quadratic-Taylor consistency — sigma(eps+d*v) - sigma(eps) - C:(d*v) = O(d^2), so the
    # residual ratio -> ~4 as d halves. Strictly stronger than step-stability: it proves C is the
    # ACTUAL derivative, not merely a step-converged FD operator. (Base state is deep plastic so
    # perturbations stay C1.) This is the gate the C++ analytic tangent will be FD-checked against.
    sig0t, _, _, _ = return_map_tensor(np.zeros(6), deps_pl, mp_na, 0.0, **HP)
    vdir = np.array([1.0, -0.3, 0.2, 0.15, -0.1, 0.05])
    tay = []
    for dd in (1.0e-4, 5.0e-5):
        s1, _, _, _ = return_map_tensor(np.zeros(6), deps_pl + dd * vdir, mp_na, 0.0, **HP)
        tay.append(nrm(s1 - (sig0t + C_pl @ (dd * vdir))))
    res["T4_taylor_ratio"] = tay[0] / tay[1] if tay[1] > 0 else 0.0

    # T5: frame objectivity — rotate strain by Q => stress rotates by Q (isotropic map)
    sig0, _, _, _ = return_map_tensor(np.zeros(6), deps_pl, mp_na, 0.0, **HP)
    th = 0.5
    Q = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
    Emat = voigt_to_mat(deps_pl)
    s_rot, _, _, _ = return_map_tensor(np.zeros(6), mat_to_voigt(Q @ Emat @ Q.T), mp_na, 0.0, **HP)
    s_base = voigt_to_mat(sig0)
    res["T5_objectivity"] = nrm(voigt_to_mat(s_rot) - Q @ s_base @ Q.T) / nrm(s_base)

    ok = (res["T1_reduce_diag"] < 1.0e-9 and res["T2_elastic_err"] < 1.0e-6
          and res["T3_asym_nonassoc"] > 1.0e-1                       # strongly non-symmetric
          and res["T3b_sym_assoc"] < 0.05                            # ~20x smaller (recompose-spin residual)
          and res["T3_asym_nonassoc"] > 5.0 * res["T3b_sym_assoc"]   # non-association dominates
          and res["T3c_assoc_noshear_sym"] < 1.0e-6                  # associated + NO shear => symmetric
          and 4.0 < res["T3c_shear_linear"] < 6.0                    # asymmetry LINEAR in shear (=> recompose)
          and res["T4_taylor_ratio"] > 3.5                           # quadratic Taylor convergence (~4)
          and res["T5_objectivity"] < 1.0e-9)
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft}")
        print(f"  T1 tensor-reduce-to-principal (diag)   = {res['T1_reduce_diag']:.2e}  (<1e-9)")
        print(f"  T2 elastic tangent == C_elastic        = {res['T2_elastic_err']:.2e}  (<1e-6)")
        print(f"  T3 asymmetry ||C-C^T||/||C|| (Df=0.3)  = {res['T3_asym_nonassoc']:.3e}  (non-associated)")
        print(f"  T3b associated limit (e=1,Df=1)        = {res['T3b_sym_assoc']:.2e}  (~20x smaller)")
        print(f"  T3c assoc+NO-shear symmetric           = {res['T3c_assoc_noshear_sym']:.2e}  (<1e-6) ;"
              f" shear-linearity = {res['T3c_shear_linear']:.2f} (~5 => SPECTRAL RECOMPOSE, not theta-freeze)")
        print(f"  T4 Taylor quadratic ratio (~4)         = {res['T4_taylor_ratio']:.2f}  (>3.5)")
        print(f"  T5 frame objectivity                   = {res['T5_objectivity']:.2e}  (<1e-9)")
        print(f"  => TANGENT GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P2 — dual scalar DAMAGE (CDPM2 §2.3, Grassl et al. 2013), crack-band regularized.
#
# This is where the stress PEAK + softening come from. P1 effective-stress plasticity is
# MONOTONIC (no peak — the failure surface is reached at kappa_p=1); damage then knocks the
# NOMINAL stress down below the effective stress. PINNED to Grassl 2013 by Eq. number:
#   Eq.1  sigma = (1-wt) sig_t_eff + (1-wc) sig_c_eff   (Macauley split of the EFFECTIVE stress)
#   Eq.38 uniaxial-tension equivalent strain  eps_eq = sig_t_eff/E ; onset eps0 = ft/E
#         (<=> qh2 = eps_eq/eps0 = 1 <=> kappa_p = 1 — damage starts exactly at the P1 failure
#          surface, so pre-peak is pure plasticity and post-peak is damage; no double-count)
#   Eq.6-7 history kappa_dt = max(eps_eq) ; Eq.51-59 softening, crack-band eps_f = wf/lch
#          with wf = Gf/ft (exponential)  =>  damage dissipation over the INELASTIC driver
#          (kappa_dt-eps0) times lch  ==  Gf  (Bazant crack-band size-objectivity)
#
# P2a SCOPE (this slice): TENSILE damage wt on the uniaxial path + the ADR §4.3 BLOCKING Gf
# energy-objectivity gate. Compression wc + the alpha_c spectral split (Eq.46) + unilateral
# crack-closure recovery + the dual-projector damaged tangent are the P2b+ increments.
# ===========================================================================
def _solve_omega_bracketed(kd1, kd2, sig_eff, f, eps_f):
    """Solve the implicit damage residual F(w) = (1-w)*sig_eff - f*exp(-(kd1+w*kd2)/eps_f) = 0 for
    w in [0,1]. The exponential softening makes eps_i = kd1 + w*kd2 carry w, so F is NON-MONOTONE in
    w when eps_f is small (steep softening) — a raw clamped Newton can step the wrong way and STALL at
    w=0 (the cracked material spuriously HEALS: PR #261 adversarial-review CRITICAL). During damage
    the root is ALWAYS bracketed: F(0)=sig_eff - f*exp(-kd1/eps_f) > 0, F(1) = -f*exp(...) < 0 — so use
    a SAFEGUARDED Newton-with-bisection-fallback that never leaves the bracket (Brent-lite)."""
    def Fof(w):
        return (1.0 - w) * sig_eff - f * np.exp(-(kd1 + w * kd2) / eps_f)
    lo, hi = 0.0, 1.0
    if Fof(lo) <= 0.0:       # not damage-loading (or already relaxed) => no damage
        return 0.0
    if Fof(hi) >= 0.0:       # fully damaged
        return 1.0
    w = 0.5
    for _ in range(100):
        Fw = Fof(w)
        if abs(Fw) < 1.0e-13 * (f + 1.0):
            break
        if Fw > 0.0:
            lo = w
        else:
            hi = w
        dF = -sig_eff + f * np.exp(-(kd1 + w * kd2) / eps_f) * kd2 / eps_f   # safeguarded Newton...
        wn = w - Fw / dF if dF != 0.0 else w
        w = wn if (lo < wn < hi) else 0.5 * (lo + hi)                        # ...fall back to bisection
    return min(1.0, max(0.0, w))


def _solve_omega_t_exp(kappa_dt, kdt1, kdt2, sig_t_eff, E, ft, eps_f):
    """CDPM2 tensile damage, EXPONENTIAL softening (Eq.55) with the inelastic-strain split Eq.52.
    eps_i = kdt1 + wt*kdt2 (Eq.52); nominal sig_t_nom = ft*exp(-eps_i/eps_f) (Eq.55); and
    sig_t_nom = (1-wt)*sig_t_eff (Eq.1). wt IMPLICIT (eps_i carries wt) => bracketed root solve.
    eps_f = wf/lch = Gf/(ft*lch) so int sig_t_nom d eps_i == ft*eps_f == Gf/lch (size-objective)."""
    if sig_t_eff <= 0.0 or kappa_dt <= ft / E:
        return 0.0
    return _solve_omega_bracketed(kdt1, kdt2, sig_t_eff, ft, eps_f)


def drive_uniaxial_tension_damaged(mp, eps11_path, Gf, lch):
    """Uniaxial-STRESS tension: P1 effective-stress return + CDPM2 tensile damage (Eq.1,38,44-45,
    52,55). Damage driver = the equivalent strain kappa_dt = max(eps_eq), eps_eq = sig_t_eff/E
    (Eq.38). The INELASTIC strain eps_i = kdt1 + wt*kdt2 (Eq.52): kdt1 = accumulated plastic-strain
    norm /x_s (Eq.44), kdt2 = (kappa_dt-eps0)/x_s (Eq.45); x_s = softening ductility (Eq.56, =1 in
    uniaxial tension since Rs=0). wt solved implicitly per step. Tracks eps_i for the energy gate."""
    E, ft, nu = mp["E"], mp["ft"], mp["nu"]
    eps0 = ft / E
    eps_f = Gf / (ft * lch)
    eps = np.zeros(3); sig_eff = np.zeros(3); kp = 0.0; el = 0.0
    kappa_dt = 0.0; kdt1 = 0.0; sigt_max = 0.0             # P2g: monotone running-max drive (no heal)
    epl_prev = np.zeros(3)
    out = {k: [] for k in ("eps11", "sig11", "wt", "kp", "sig_eff", "epsi", "kappa_dt")}
    for e11 in eps11_path:
        for _ in range(80):                       # lateral Newton: EFFECTIVE lateral stress -> 0
            deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
            snew, _, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps, mp), mp, kp)
            res = 0.5 * (snew[1] + snew[2])
            if abs(res) < 1.0e-10 * (mp["fc"] + 1.0):
                break
            d = 1.0e-8 * (abs(el) + 1.0e-6)
            deps2 = np.array([e11 - eps[0], (el + d) - eps[1], (el + d) - eps[2]])
            snew2, _, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps2, mp), mp, kp)
            Jd = (0.5 * (snew2[1] + snew2[2]) - res) / d
            if abs(Jd) < 1.0e-12:
                Jd = 1.0e-12 if Jd >= 0 else -1.0e-12
            el -= res / Jd
        deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
        sig_eff, kp, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps, mp), mp, kp)
        eps = np.array([e11, el, el])
        s_eff = sig_eff[0]
        # equivalent strain (Eq.38 uniaxial) + history (Eq.6-7)
        eps_eq = max(s_eff, 0.0) / E
        kappa_dt_new = max(kappa_dt, eps_eq)
        loading = kappa_dt_new > kappa_dt + 1.0e-300
        # plastic strain (principal): eps_p = eps_total - effective-elastic strain
        epl = eps - np.array([(sig_eff[0] - nu * (sig_eff[1] + sig_eff[2])) / E,
                              (sig_eff[1] - nu * (sig_eff[0] + sig_eff[2])) / E,
                              (sig_eff[2] - nu * (sig_eff[0] + sig_eff[1])) / E])
        if kappa_dt_new > eps0 and loading:                 # Eq.44: accumulate ||d eps_p||/x_s (x_s=1 tension)
            kdt1 += float(np.linalg.norm(epl - epl_prev))
        epl_prev = epl
        kappa_dt = kappa_dt_new
        kdt2 = max(kappa_dt - eps0, 0.0)                    # Eq.45 (x_s=1)
        sigt_max = max(sigt_max, max(s_eff, 0.0))           # P2g: monotone drive => omega_t never heals
        wt = _solve_omega_t_exp(kappa_dt, kdt1, kdt2, sigt_max, E, ft, eps_f)
        epsi = kdt1 + wt * kdt2                             # Eq.52
        sig11_nom = (1.0 - wt) * s_eff                      # Eq.1
        for k, v in (("eps11", e11), ("sig11", sig11_nom), ("wt", wt),
                     ("kp", kp), ("sig_eff", s_eff), ("epsi", epsi), ("kappa_dt", kappa_dt)):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


# ---------------------------------------------------------------------------
# P2b — COMPRESSIVE damage wc + the alpha_c tension/compression split (CDPM2 Eq.37,46-57).
# ---------------------------------------------------------------------------
def equiv_strain_general(sig_pr, mp):
    """CDPM2 equivalent strain eps_tilde (Eq.37): the closed form from setting f=0 with qh1=1,
    qh2=eps_tilde/eps0. Reduces to eps0 on the failure surface (ANY state) and to sig_t/E in
    uniaxial tension (Eq.38). sig_pr = 3 EFFECTIVE principal stresses."""
    fc, ft, E, m0, e = mp["fc"], mp["ft"], mp["E"], mp["m0"], mp["e"]
    eps0 = ft / E
    sv = np.array([sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0])
    xi, rho, theta, *_ = invariants(sv)
    sigV = xi / SQRT3
    A = rho * lode_r(theta, e) / (SQRT6 * fc) + sigV / fc
    rad = (eps0 * eps0 * m0 * m0 / 4.0) * A * A + 3.0 * eps0 * eps0 * rho * rho / (2.0 * fc * fc)
    return (eps0 * m0 / 2.0) * A + np.sqrt(max(rad, 0.0))


def alpha_compression(sig_pr):
    """alpha_c (Eq.46): 0 for pure tension, 1 for pure compression. sig_pc=negative principal parts,
    sig_pt=positive parts; alpha_c = sum_i sig_pc,i*(sig_pt,i+sig_pc,i)/||sig_p||^2 = sum sig_pc,i*sig_i/||sig||^2."""
    s = np.asarray(sig_pr, float)
    nrm2 = float(np.dot(s, s))
    if nrm2 <= 1.0e-300:
        return 0.0
    spc = np.minimum(s, 0.0)                       # negative parts
    return float(np.sum(spc * s) / nrm2)


def beta_c(sig_pr, kp, mp):
    """CDPM2 beta_c (Eq.50): the factor scaling the PLASTIC-strain part of the compressive damage
    driver kappa_dc1 (Eq.48). beta_c = ft * qh2(kp) * sqrt(2/3) / (rho_bar * sqrt(1 + 2 Df^2)),
    rho_bar = sqrt(2 J2) of the EFFECTIVE stress. It provides "a smooth transition from pure damage to
    damage-plasticity softening during cyclic loading" (Grassl 2013 2.3.5): the more deviatoric the
    effective stress (large rho_bar), the LESS the plastic strain feeds compressive damage. In monotonic
    uniaxial compression on the failure surface (qh2=1, rho_bar=fc*sqrt(2/3)) it is ~ft/(fc*sqrt(1+2Df^2))
    (~0.058 for ft/fc=0.1, Df=1) — so it makes compression markedly MORE DUCTILE than the beta_c=1
    simplification the monotonic slice used (P2b/P2c). Guarded for rho_bar->0 (hydrostatic); clamped to
    [0,1] (a plastic-contribution fraction is physically <=1 — inactive in the damaging regime where
    rho_bar is large, so the clamp never binds there and the analytic tangent stays smooth).
    sig_pr = 3 EFFECTIVE principal stresses (or the deviatoric invariant can be passed via rho)."""
    ft, Df, Hp = mp["ft"], mp["Df"], mp["Hp"]
    sv = np.array([sig_pr[0], sig_pr[1], sig_pr[2], 0.0, 0.0, 0.0])
    _, rho, *_ = invariants(sv)
    if rho <= 1.0e-12:
        return 1.0                                    # hydrostatic: no deviatoric drive -> fall back (clamp)
    bc = ft * qh2(kp, Hp) * np.sqrt(2.0 / 3.0) / (rho * np.sqrt(1.0 + 2.0 * Df * Df))
    return min(max(bc, 0.0), 1.0)


def _solve_omega_c_exp(kdc1, kdc2, sig_c_mag, fc, eps_fc):
    """Compressive damage, exponential softening (Eq.55 analog with fc, Gc). eps_i = kdc1 + wc*kdc2
    (Eq.52 analog); |sig_c_nom| = fc*exp(-eps_i/eps_fc) = (1-wc)*|sig_c_eff|; wc implicit -> Newton.
    eps_fc = wfc/lch = Gc/(fc*lch) => int |sig_c_nom| d eps_i == fc*eps_fc == Gc/lch (objective).
    wc IMPLICIT (eps_i carries wc) => bracketed root solve (same non-monotone-F healing guard as wt)."""
    if sig_c_mag <= 0.0:
        return 0.0
    return _solve_omega_bracketed(kdc1, kdc2, sig_c_mag, fc, eps_fc)


def drive_uniaxial_compression_damaged(mp, eps11_path, Gc, lch, As=5.0, beta_c_on=True):
    """Uniaxial-STRESS compression: P1 effective return + CDPM2 compressive damage (Eq.37,46-57).
    Tracks the compressive history kdc1 (Eq.48, alpha_c*beta_c-scaled plastic strain) and
    kdc2 (Eq.49, alpha_c-scaled equivalent strain /x_s). x_s = 1+(As-1)*Rs, Rs=-sqrt6 sigV/rho for
    sigV<=0 (Eq.56-57, =As in uniaxial compression -> ductility). Returns nominal axial stress.
    `beta_c_on` (default True) applies the full CDPM2 beta_c (Eq.50); False forces beta_c=1 (the old
    monotonic simplification) for the P2f non-tautology comparison (real beta_c is markedly more ductile)."""
    E, fc, nu, Df = mp["E"], mp["fc"], mp["nu"], mp["Df"]
    ft = mp["ft"]; eps0 = ft / E
    eps_fc = Gc / (fc * lch)
    eps = np.zeros(3); sig_eff = np.zeros(3); kp = 0.0; el = 0.0
    et_prev = 0.0; kdc = 0.0; kdc1 = 0.0; sigc_max = 0.0   # P2g: monotone running-max drive (no heal)
    epl_prev = np.zeros(3)
    out = {k: [] for k in ("eps11", "sig11", "wc", "kp", "sig_eff", "epsi", "alpha_c")}
    for e11 in eps11_path:
        for _ in range(80):                       # lateral Newton: EFFECTIVE lateral stress -> 0
            deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
            snew, _, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps, mp), mp, kp)
            res = 0.5 * (snew[1] + snew[2])
            if abs(res) < 1.0e-10 * (mp["fc"] + 1.0):
                break
            d = 1.0e-8 * (abs(el) + 1.0e-6)
            deps2 = np.array([e11 - eps[0], (el + d) - eps[1], (el + d) - eps[2]])
            snew2, _, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps2, mp), mp, kp)
            Jd = (0.5 * (snew2[1] + snew2[2]) - res) / d
            if abs(Jd) < 1.0e-12:
                Jd = 1.0e-12 if Jd >= 0 else -1.0e-12
            el -= res / Jd
        deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
        sig_eff, kp, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps, mp), mp, kp)
        eps = np.array([e11, el, el])
        et = equiv_strain_general(sig_eff, mp)               # Eq.37
        det = max(et - et_prev, 0.0); et_prev = et
        ac = alpha_compression(sig_eff)                      # Eq.46
        sv6 = np.array([sig_eff[0], sig_eff[1], sig_eff[2], 0, 0, 0])
        xi, rho, theta, *_ = invariants(sv6); sigV = xi / SQRT3
        Rs = (-SQRT6 * sigV / rho) if (sigV <= 0.0 and rho > 1.0e-12) else 0.0     # Eq.57
        xs = 1.0 + (As - 1.0) * Rs                                                 # Eq.56
        epl = eps - np.array([(sig_eff[0] - nu * (sig_eff[1] + sig_eff[2])) / E,
                              (sig_eff[1] - nu * (sig_eff[0] + sig_eff[2])) / E,
                              (sig_eff[2] - nu * (sig_eff[0] + sig_eff[1])) / E])
        if et > eps0 and det > 0.0:                          # damage active (kappa_p>1)
            kdc += ac * det / xs                             # Eq.47+49: kdc2 driver (alpha_c eps_tilde /x_s)
            bc = beta_c(sig_eff, kp, mp) if beta_c_on else 1.0   # Eq.50 (P2f) cyclic damage<->plasticity factor
            kdc1 += ac * bc / xs * float(np.linalg.norm(epl - epl_prev))       # Eq.48 (full CDPM2 beta_c)
        epl_prev = epl
        kdc2 = kdc
        sig_c_mag = max(-sig_eff[0], 0.0)                    # compressive axial magnitude (uniaxial)
        sigc_max = max(sigc_max, sig_c_mag)                  # P2g: monotone drive => omega_c never heals
        wc = _solve_omega_c_exp(kdc1, kdc2, sigc_max, fc, eps_fc)
        epsi = kdc1 + wc * kdc2
        sig11_nom = (1.0 - wc) * sig_eff[0]                  # Eq.1 (compression branch; sig_eff[0]<0)
        for k, v in (("eps11", e11), ("sig11", sig11_nom), ("wc", wc),
                     ("kp", kp), ("sig_eff", sig_eff[0]), ("epsi", epsi), ("alpha_c", ac)):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def run_p2b_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gc=5.0, As=2.0, verbose=True):
    mp = make_material(E, nu, fc, ft)
    res = {}
    # C0: the alpha_c split — ~0 in uniaxial tension, ~1 in uniaxial compression (Eq.46)
    res["C0_alpha_tens"] = alpha_compression([ft, 0.0, 0.0])
    res["C0_alpha_comp"] = alpha_compression([-fc, 0.0, 0.0])
    res["C0_ok"] = bool(res["C0_alpha_tens"] < 1e-9 and abs(res["C0_alpha_comp"] - 1.0) < 1e-9)
    # equivalent strain sanity: eps_tilde == eps0 on the failure surface (uniaxial comp -fc, tens +ft)
    res["C0_eqstrain_comp"] = equiv_strain_general([-fc, 0, 0], mp) / (ft / E)   # -> ~1
    res["C0_eqstrain_tens"] = equiv_strain_general([ft, 0, 0], mp) / (ft / E)    # -> ~1
    res["C0_eqstrain_ok"] = bool(abs(res["C0_eqstrain_comp"] - 1.0) < 1e-6 and abs(res["C0_eqstrain_tens"] - 1.0) < 1e-6)

    # C1: nominal uniaxial-compression PEAKS at fc then softens (the peak is damage); eff monotone.
    # (P2f: range extended 0.15->0.4 — the full beta_c makes compression more ductile, so it needs more
    # strain to soften to <10% fc.)
    d = drive_uniaxial_compression_damaged(mp, np.linspace(0, -0.4, 4000), Gc, lch=50.0, As=As)
    res["C1_peak"] = float(-np.min(d["sig11"]))
    res["C1_peak_err"] = abs(res["C1_peak"] / fc - 1.0)
    res["C1_softens"] = bool(-d["sig11"][-1] < 0.1 * fc and d["wc"][-1] > 0.85)
    res["C1_eff_monotone"] = bool(-d["sig_eff"][-1] > -d["sig_eff"][np.argmin(d["sig11"])])
    # C0b onset coincidence (PR #261 review): compression damage initiates at kappa_p=1 / sig_eff=-fc
    ci = np.argmax(d["wc"] > 1.0e-9)
    res["C0b_kp_at_onset"] = float(d["kp"][ci]) if d["wc"][ci] > 1.0e-9 else -1.0
    res["C0b_ok"] = bool(abs(res["C0b_kp_at_onset"] - 1.0) < 0.05
                         and abs(-d["sig_eff"][ci] / fc - 1.0) < 0.05)

    # C2 — crack-band SOFTENING-LAW lch-scaling (the eps_fc WIRING; Gc/lch BY CONSTRUCTION), NOT an
    # independent objectivity proof (PR #261 review: int over eps_i is tautological). C3 reports the
    # honest FE-visible total. NOTE (review MEDIUM): uniaxial compression has Rs=1 exactly, so x_s=As
    # is constant across all lch here — this leg does NOT exercise the confinement-ductility (x_s) effect
    # on dissipation; a multi-confinement Gc leg (sigma3!=0, ADR §4.3 [MAJOR]) is a documented P2c gap.
    # P2f re-gate: the full CDPM2 beta_c (Eq.50) makes compression markedly more ductile, so over a
    # FIXED compression strain the softening no longer fully completes at small lch (large eps_fc). The
    # omega_c solve enforces |sig_c_nom| = fc*exp(-eps_i/eps_fc) EXACTLY, so the integral tail from the
    # path's last eps_i to infinity is the analytic exp tail eps_fc*|sig_c_last|; adding it removes the
    # truncation (the dominant error) and keeps the path short/fast — the residual is the trapezoidal
    # discretization of the post-peak exponential.
    gc_lch = {}
    for lch in (50.0, 100.0, 200.0):
        dd = drive_uniaxial_compression_damaged(mp, np.linspace(0, -0.4, 4000), Gc, lch, As=As)
        epsi = dd["epsi"]
        area = float(np.sum(0.5 * (-dd["sig11"][1:] - dd["sig11"][:-1]) * np.diff(epsi)))  # |sig| d epsi
        tail = (Gc / (fc * lch)) * float(-dd["sig11"][-1])     # eps_fc * |sig_c_last| (analytic exp tail)
        gc_lch[lch] = (area + tail) * lch
    res["C2_gc_lch"] = gc_lch
    res["C2_max_rel_err"] = max(abs(gc_lch[l] / Gc - 1.0) for l in gc_lch)

    # C3 — HONEST FE-visible compression energy (REPORTED, not gated): same CDPM2 damage-only-
    # regularization caveat as tension D3 (un-regularized effective-plasticity dissipation makes the
    # FE-visible total lch-dependent).
    wtot = {}
    for lch in (50.0, 100.0, 200.0):
        dd = drive_uniaxial_compression_damaged(mp, np.linspace(0, -0.15, 4000), Gc, lch, As=As)
        s = -dd["sig11"]; et = -dd["eps11"]
        Wel = 0.5 * float(s[-1]) ** 2 / E
        wtot[lch] = (float(np.sum(0.5 * (s[1:] + s[:-1]) * np.diff(et))) - Wel) * lch
    res["C3_wtot_lch"] = wtot
    res["C3_total_spread"] = (max(wtot.values()) - min(wtot.values())) / Gc

    ok = (res["C0_ok"] and res["C0_eqstrain_ok"] and res["C0b_ok"] and res["C1_peak_err"] < 0.03
          and res["C1_softens"] and res["C1_eff_monotone"] and res["C2_max_rel_err"] < 0.02)
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gc={Gc} As={As}")
        print(f"  C0 alpha_c: tension={res['C0_alpha_tens']:.3e} compression={res['C0_alpha_comp']:.4f}"
              f"  eps_tilde/eps0 on surface: comp={res['C0_eqstrain_comp']:.5f} tens={res['C0_eqstrain_tens']:.5f}")
        print(f"  C0b onset: kappa_p={res['C0b_kp_at_onset']:.4f} (=1)  ok={res['C0b_ok']}")
        print(f"  C1 nominal peak = {res['C1_peak']:.4f} (target fc={fc}) err={res['C1_peak_err']:.2e}"
              f"  eff-monotone={res['C1_eff_monotone']}  softens={res['C1_softens']}")
        print(f"  C2 crack-band softening-law wiring (Gc/lch BY CONSTRUCTION): "
              + "  ".join(f"lch={int(l)}:{gc_lch[l]:.3f}" for l in sorted(gc_lch)) + f"  err={res['C2_max_rel_err']:.2e}")
        print(f"  C3 [REPORTED, not gated] FE-visible TOTAL/crack-area: "
              + "  ".join(f"lch={int(l)}:{wtot[l]:.3f}" for l in sorted(wtot))
              + f"  spread={res['C3_total_spread']:.2f}*Gc  [CDPM2 damage-only; confined x_s leg = P2c gap]")
        print(f"  => P2b GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P2c — UNIFIED TENSOR dual damage: the spectral T/C split on a GENERAL effective
# stress + automatic per-step unilateral crack-closure (CDPM2 Eq.1, Grassl 2013).
#
# P2a/P2b validated the two scalar softening laws on separate uniaxial-STRESS drivers.
# P2c fuses them into ONE constitutive update (the shape the C++ nDMaterial wrapper calls):
#   1. effective-stress return map (P1 kernel) -> sig_bar (6-tensor / principal)
#   2. spectral split sig_bar = sig_bar_t + sig_bar_c  (Macaulay on eigenvalues, Eq.1)
#   3. ONE equivalent strain eps_tilde (Eq.37) + alpha_c (Eq.46) drive BOTH histories:
#        tension  kappa_dt  <- FULL eps_tilde         (Eq.43:  d kappa_dt = d eps_tilde)
#        compr.   kappa_dc  <- alpha_c * eps_tilde    (Eq.47:  d kappa_dc = alpha_c d eps_tilde)
#      with kdt1/kdc1 the plastic-strain parts (Eq.44/48, beta_c=1 monotonic) and
#      kdt2/kdc2 the damage-scaled parts (Eq.45/49), all /x_s ductility (Eq.56-57).
#   4. omega_t, omega_c solved implicitly per channel (Eq.52,55; bracketed root, no healing)
#   5. NOMINAL sigma = (1-omega_t) sig_bar_t + (1-omega_c) sig_bar_c   (Eq.1)
#
# UNILATERAL recovery is AUTOMATIC and tier-independent (ADR 4.3 BLOCKING): because the
# split is recomputed from the CONVERGED effective stress EVERY step, a principal that flips
# negative (crack closing) is routed into the (1-omega_c) channel — it is no longer multiplied
# by (1-omega_t), so the cracked stiffness is recovered with NO extra state. The partial-recovery
# knob s_rec * g_close (ADR 4.3) and the dual-projector analytic tangent + beta_c cyclic are P2d.
#
# DRIVE CHOICE (reduces EXACTLY to the gated P2a/P2b): each channel's scalar softening law is
# driven by the EXTREME effective principal of its sign — sig_t_drive = max_i<sig_bar_i>+ ,
# sig_c_drive = max_i<-sig_bar_i>+ . In uniaxial tension/compression these collapse to the single
# axial effective stress P2a/P2b used. The genuinely MULTIAXIAL apportioning (biaxial/triaxial
# peak: extreme-principal vs ||sig_bar_t|| norm) is gated coaxially here and pinned in P2d.
# ===========================================================================
def spectral_split_principal(sig_pr):
    """Macaulay split of principal stresses into tension/compression parts (Eq.1). Returns
    (st, sc) with st = <sig>+ (>=0), sc = <sig>- (<=0), st + sc == sig elementwise."""
    s = np.asarray(sig_pr, float)
    return np.maximum(s, 0.0), np.minimum(s, 0.0)


def apply_damage_principal(sig_pr, wt, wc):
    """Nominal principal stress from the dual-scalar spectral split (Eq.1):
    sig = (1-wt) <sig_bar>+ + (1-wc) <sig_bar>- . Unilateral recovery is automatic — a
    principal < 0 sits in the (1-wc) channel regardless of wt."""
    st, sc = spectral_split_principal(sig_pr)
    return (1.0 - wt) * st + (1.0 - wc) * sc


def damaged_stress_tensor(sig_eff6, wt, wc):
    """Full 6-tensor nominal stress: eigendecompose the effective stress, split its eigenvalues,
    apply (1-wt)/(1-wc), recompose on the SAME eigenvectors. The frame-objective form of
    apply_damage_principal — the eigenvectors carry the rotation."""
    w, V = np.linalg.eigh(voigt_to_mat(sig_eff6))
    sp = apply_damage_principal(w, wt, wc)
    return mat_to_voigt(V @ np.diag(sp) @ V.T)


def _damage_drivers(sig_eff_pr, mp, As):
    """Per-step damage kinematics from the EFFECTIVE principal stress: equivalent strain
    eps_tilde (Eq.37), compressive weight alpha_c (Eq.46), softening ductility x_s (Eq.56-57)."""
    et = equiv_strain_general(sig_eff_pr, mp)                       # Eq.37
    ac = alpha_compression(sig_eff_pr)                             # Eq.46
    sv6 = np.array([sig_eff_pr[0], sig_eff_pr[1], sig_eff_pr[2], 0, 0, 0])
    xi, rho, _, *_ = invariants(sv6)
    sigV = xi / SQRT3
    Rs = (-SQRT6 * sigV / rho) if (sigV <= 0.0 and rho > 1.0e-12) else 0.0   # Eq.57
    xs = 1.0 + (As - 1.0) * Rs                                               # Eq.56
    return et, ac, xs


def drive_damaged_unified(mp, eps11_path, Gf, Gc, lch, As=2.0, sigma3=0.0):
    """UNIFIED dual-damage uniaxial-STRESS driver (lateral mixed control). Default `sigma3=0` holds the
    EFFECTIVE lateral stress at 0 (uniaxial stress); `sigma3>0` holds it at -sigma3 (ACTIVE confinement
    — a small confining pressure keeps the lateral principals clearly compressive, off the sigma_lat=0
    Macaulay kink, and is the STRESS-controlled = cross-platform-robust way to a confined-compression-
    damaged state). Handles tension, compression, AND tension<->compression reversal in one update —
    the spectral split routes each principal to its channel every step (automatic unilateral).
    Returns nominal axial stress + both damage variables + effective stress (for ratio checks)."""
    E, fc, ft, nu = mp["E"], mp["fc"], mp["ft"], mp["nu"]
    eps0 = ft / E
    eps_f = Gf / (ft * lch)
    eps_fc = Gc / (fc * lch)
    eps = np.zeros(3); sig_eff = np.zeros(3); kp = 0.0; el = 0.0
    et_max = 0.0                                          # running max of eps_tilde (Eq.43 history)
    sigt_max = sigc_max = 0.0                             # P2g: monotone running-max effective drive (no heal)
    kdt1 = kdt2 = kdc = kdc1 = kdc2 = 0.0
    epl_prev = np.zeros(3)
    out = {k: [] for k in ("eps11", "eps_lat", "sig11", "wt", "wc", "sig_eff", "kp", "epsi_t", "epsi_c")}
    for e11 in eps11_path:
        for _ in range(80):                              # lateral Newton: EFFECTIVE lateral -> -sigma3
            deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
            snew, _, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps, mp), mp, kp)
            res = 0.5 * (snew[1] + snew[2]) + sigma3        # target EFFECTIVE lateral stress = -sigma3
            if abs(res) < 1.0e-10 * (mp["fc"] + 1.0):
                break
            d = 1.0e-8 * (abs(el) + 1.0e-6)
            deps2 = np.array([e11 - eps[0], (el + d) - eps[1], (el + d) - eps[2]])
            snew2, _, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps2, mp), mp, kp)
            Jd = (0.5 * (snew2[1] + snew2[2]) + sigma3 - res) / d
            if abs(Jd) < 1.0e-12:
                Jd = 1.0e-12 if Jd >= 0 else -1.0e-12
            el -= res / Jd
        deps = np.array([e11 - eps[0], el - eps[1], el - eps[2]])
        sig_eff, kp, _, _, _ = return_map_hardening(_elastic_pred(sig_eff, deps, mp), mp, kp)
        eps = np.array([e11, el, el])

        # damage kinematics (frame-invariant; computed from the EFFECTIVE stress)
        et, ac, xs = _damage_drivers(sig_eff, mp, As)
        det_raw = max(et - et_max, 0.0)                   # Eq.43 raw increment d kappa_dt
        above = max(et - max(et_max, eps0), 0.0)          # increment ABOVE the onset eps0 (kappa_p>1)
        loading = det_raw > 0.0 and et > eps0             # damage active <=> past the failure surface
        # plastic principal strain (effective elastic strain removed)
        epl = eps - np.array([(sig_eff[0] - nu * (sig_eff[1] + sig_eff[2])) / E,
                              (sig_eff[1] - nu * (sig_eff[0] + sig_eff[2])) / E,
                              (sig_eff[2] - nu * (sig_eff[0] + sig_eff[1])) / E])
        dnorm_epl = float(np.linalg.norm(epl - epl_prev))
        if loading:
            # kdt2/kdc2 integrate d kappa_d/x_s from the onset eps0 (Eq.45/49); clamping the lower
            # limit at eps0 telescopes to max(kappa_dt-eps0,0) when x_s=1, so the tension channel is
            # byte-identical to the P2a driver. kdt1/kdc1 take the FULL plastic-strain increment
            # (Eq.44/48), matching P2a/P2b. (Variable-x_s onset harmonization across channels = P2d.)
            # coaxial uniaxial-stress path: stress frame = principal axes (V=I), depl principal = epl-epl_prev
            wt_w = tensile_damage_weight(mp, ac, np.array([(epl - epl_prev)[0], (epl - epl_prev)[1],
                                                           (epl - epl_prev)[2], 0.0, 0.0, 0.0]),
                                         sig_eff, np.eye(3))    # P2h ctTemper weight (1 if 'none')
            kdt2 += wt_w * above / xs                     # Eq.45  d kappa_dt2 = d kappa_dt / x_s
            kdt1 += wt_w * dnorm_epl / xs                 # Eq.44  (w_t tempers compression->tension coupling)
            kdc += ac * above                             # Eq.47  d kappa_dc = alpha_c d eps_tilde
            kdc2 += ac * above / xs                       # Eq.49
            bc = beta_c(sig_eff, kp, mp)                  # Eq.50  (P2f) cyclic damage<->plasticity factor
            kdc1 += ac * bc * dnorm_epl / xs              # Eq.48  (full CDPM2 beta_c)
        et_max = max(et_max, et)                          # ALWAYS track the eps_tilde history (Eq.43)
        epl_prev = epl

        # P2i — multiaxial-consistent TENSILE drive E*eps_tilde (Eq.37), gated by a real tensile
        # principal; compressive drive stays the extreme principal (ft-scaled eps_tilde can't onset wc).
        # Reduces to the extreme principal uniaxially (E*eps_tilde == sig_bar_t). See damaged_step_tensor.
        sig_t_drive = E * et if float(np.max(sig_eff)) > 1.0e-6 * ft else 0.0
        sig_c_drive = max(-float(np.min(sig_eff)), 0.0)
        # P2g — MONOTONE (no-heal) cyclic damage: drive each omega with the running MAX of its drive
        # stress (see damaged_step_tensor). On an elastic unload the live drive drops but the max (and the
        # frozen histories) keep omega fixed => no healing; the nominal stress unloads along the secant
        # (1-omega)*sig_bar. Identical to the live drive on monotonic loading (max == live).
        sigt_max = max(sigt_max, sig_t_drive)
        sigc_max = max(sigc_max, sig_c_drive)
        # PHYSICAL FLOOR on the softening drive: only solve omega when the extreme principal of that
        # sign is a REAL stress (> 1e-6 * strength), never on the ~1e-10 MPa lateral-Newton residual.
        # Without it the residual's SIGN spuriously flips wt 0<->1 in pure compression (review-fix).
        wt = _solve_omega_bracketed(kdt1, kdt2, sigt_max, ft, eps_f) if (et_max > eps0 and sigt_max > 1.0e-6 * ft) else 0.0
        wc = _solve_omega_bracketed(kdc1, kdc2, sigc_max, fc, eps_fc) if (kdc > 0.0 and sigc_max > 1.0e-6 * fc) else 0.0

        sig_nom = apply_damage_principal(sig_eff, wt, wc)            # Eq.1
        for k, v in (("eps11", e11), ("eps_lat", el), ("sig11", sig_nom[0]), ("wt", wt), ("wc", wc),
                     ("sig_eff", sig_eff[0]), ("kp", kp),
                     ("epsi_t", kdt1 + wt * kdt2), ("epsi_c", kdc1 + wc * kdc2)):
            out[k].append(v)
    return {k: np.array(v) for k, v in out.items()}


def run_p2c_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    mp = make_material(E, nu, fc, ft)
    res = {}

    # DT0: the pure-kinematics identities the unified update rests on.
    #   - spectral split is a partition: st + sc == sig, st>=0>=sc
    #   - apply_damage with wt=wc=0 is the identity (no damage => effective stress unchanged)
    sp = np.array([5.0, -2.0, -7.0])
    st, sc = spectral_split_principal(sp)
    res["DT0_split_partition"] = float(np.max(np.abs((st + sc) - sp)))
    res["DT0_identity"] = float(np.max(np.abs(apply_damage_principal(sp, 0.0, 0.0) - sp)))
    # DT0c: the UNILATERAL routing mechanism, tested DIRECTLY (the sharp check the path-level DT3
    # cannot make — under reversal the driver floors wt to 0 once all principals are compressive).
    # A high tensile damage wt=0.95 must leave the COMPRESSIVE principals scaled by (1-wc) ONLY
    # (wt-invariant = the crack closes), and the tensile principal scaled by (1-wt).
    mixed = np.array([-5.0, -1.0, 2.0])
    expect = np.array([0.8 * -5.0, 0.8 * -1.0, 0.05 * 2.0])           # (1-0.2)*compr, (1-0.95)*tens
    res["DT0_unilateral"] = float(np.max(np.abs(apply_damage_principal(mixed, 0.95, 0.2) - expect)))
    wt_invar = apply_damage_principal(mixed, 0.95, 0.2)[:2]
    wt_invar2 = apply_damage_principal(mixed, 0.10, 0.2)[:2]          # change wt only
    res["DT0_compr_wt_invariant"] = float(np.max(np.abs(wt_invar - wt_invar2)))   # compressive entries unchanged
    # DT0d: the physical FLOOR — pure monotonic compression must NOT spuriously activate tensile
    # damage off the ~1e-10 MPa lateral-Newton residual (review-fix; was wt->1, mask-hidden behind
    # the compressive channel). Mirror: pure tension must not activate wc.
    dpc = drive_damaged_unified(mp, np.linspace(0, -0.05, 1500), Gf, Gc, 50.0, As)
    dpt = drive_damaged_unified(mp, np.linspace(0, 0.004, 1500), Gf, Gc, 50.0, As)
    res["DT0_pure_compression_wt"] = float(np.max(dpc["wt"]))         # ~0
    res["DT0_pure_tension_wc"] = float(np.max(dpt["wc"]))             # ~0

    # DT1: REDUCE-TO-P2a — pure uniaxial tension. The unified driver's nominal stress + wt must be
    # byte-identical to the validated tensile-only driver (same crack-band eps_f, same softening).
    lch = 50.0
    du = drive_damaged_unified(mp, np.linspace(0, 0.008, 3000), Gf, Gc, lch, As)
    da = drive_uniaxial_tension_damaged(mp, np.linspace(0, 0.008, 3000), Gf, lch)
    res["DT1_sig_maxdiff"] = float(np.max(np.abs(du["sig11"] - da["sig11"])))
    res["DT1_wt_maxdiff"] = float(np.max(np.abs(du["wt"] - da["wt"])))
    res["DT1_peak"] = float(np.max(du["sig11"]))
    res["DT1_peak_err"] = abs(res["DT1_peak"] / ft - 1.0)

    # DT2: REDUCE-TO-P2b — pure uniaxial compression vs the validated compressive-only driver.
    duc = drive_damaged_unified(mp, np.linspace(0, -0.15, 4000), Gf, Gc, lch, As)
    db = drive_uniaxial_compression_damaged(mp, np.linspace(0, -0.15, 4000), Gc, lch, As=As)
    res["DT2_sig_maxdiff"] = float(np.max(np.abs(duc["sig11"] - db["sig11"])))
    res["DT2_wc_maxdiff"] = float(np.max(np.abs(duc["wc"] - db["wc"])))
    res["DT2_peak"] = float(-np.min(duc["sig11"]))
    res["DT2_peak_err"] = abs(res["DT2_peak"] / fc - 1.0)

    # DT3: UNILATERAL crack-closure at the PATH level — load tension into softening (wt->1), reverse
    # through 0 into compression, confirm the nominal compressive branch recovers full stiffness
    # (nominal/effective -> 1). NOTE (review-fix, honest): once the reversal makes ALL principals
    # compressive the driver floors wt to 0, so in this window the recovery is the per-step RE-SPLIT
    # routing the (now-compressive) axial principal into the (1-wc~0) channel. The wt-INVARIANT
    # closed-crack mechanism itself (a compressive principal carried by (1-wc) regardless of a LIVE
    # wt) is tested directly + unconditionally by DT0_unilateral above; DT3 is the end-to-end check.
    path = np.concatenate([np.linspace(0, 0.006, 600),       # tension to deep softening
                           np.linspace(0.006, -0.012, 1200)])  # reverse into compression
    dr = drive_damaged_unified(mp, path, Gf, Gc, lch, As)
    safe = np.where(dr["sig_eff"] == 0.0, 1.0, dr["sig_eff"])
    ratio = dr["sig11"] / safe                              # nominal/effective = (1-omega) on the live channel
    ten = dr["sig_eff"] > 0.01 * ft
    comp_early = (dr["sig_eff"] < -0.02 * fc) & (dr["wc"] < 1.0e-3)   # real compression, pre-compression-damage
    res["DT3_wt_peak"] = float(np.max(dr["wt"]))
    res["DT3_min_ratio_tension"] = float(np.min(ratio[ten])) if ten.any() else 1.0       # << 1 (damaged)
    res["DT3_recovery_ratio"] = float(np.max(ratio[comp_early])) if comp_early.any() else 0.0  # ~1 (recovered)
    res["DT3_recovered"] = bool(res["DT3_min_ratio_tension"] < 0.5 and res["DT3_recovery_ratio"] > 0.98
                                and res["DT3_wt_peak"] > 0.8)

    # DT4: FRAME OBJECTIVITY of the damaged stress — rotate a damaged tensor state by Q, the nominal
    # stress rotates by Q (the spectral recompose carries the eigenvectors; wt/wc are invariant).
    sig_eff6 = np.array([4.0, -1.5, -6.0, 0.8, -0.3, 0.5])    # an arbitrary off-axis effective state
    wt0, wc0 = 0.6, 0.25
    th = 0.7
    Q = np.array([[np.cos(th), -np.sin(th), 0.0], [np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
    base = voigt_to_mat(damaged_stress_tensor(sig_eff6, wt0, wc0))
    rot_in = mat_to_voigt(Q @ voigt_to_mat(sig_eff6) @ Q.T)
    got = voigt_to_mat(damaged_stress_tensor(rot_in, wt0, wc0))
    nrm = lambda A: float(np.sqrt(np.sum(A * A)))
    res["DT4_objectivity"] = nrm(got - Q @ base @ Q.T) / nrm(base)

    # DT5 [DIAGNOSTIC, REPORTED not gated] — compression->tension damage COUPLING. Per CDPM2 Eq.43 the
    # tensile history kappa_dt tracks the FULL equivalent strain (no (1-alpha_c) factor), so a
    # compression excursion accumulates kappa_dt1/kdt2 and PRE-DAMAGES a subsequent tension reload. In
    # this MONOTONIC slice (beta_c=1, no tensile/compressive plastic-strain projection) that coupling
    # is un-tempered, so compression-then-tension currently loses ALL tensile strength. The honest
    # tempering (the cyclic beta_c transition Eq.50 + the alpha_t-weighting question — literal CDPM2 vs
    # a tensile-plastic-strain projection) is the P2f cyclic increment. Reported so the limitation is
    # tracked, NOT gated (the monotonic tension/compression responses DT1/DT2 are correct).
    cpath = np.concatenate([np.linspace(0, -0.01, 800), np.linspace(-0.01, 0.004, 800)])
    dcyc = drive_damaged_unified(mp, cpath, Gf, Gc, lch, As)
    res["DT5_tension_after_compression_peak"] = float(np.max(dcyc["sig11"][800:]))   # ~0 today (vs ft fresh)

    # DT1 (tension) telescopes EXACTLY to the P2a driver (x_s=1 => clamped kdt2 == max(kdt-eps0,0)).
    # DT2 (compression) matches P2b up to ONE onset-crossing step (the kdc2 sliver between the last
    # sub-eps0 eps_tilde and eps0; P2b doesn't clamp it) => a tight non-zero floor, vanishing under
    # step refinement. Stress diffs are in MPa (fc=30), so 1e-3 is ~3e-5 relative.
    ok = (res["DT0_split_partition"] < 1.0e-14 and res["DT0_identity"] < 1.0e-14
          and res["DT0_unilateral"] < 1.0e-14 and res["DT0_compr_wt_invariant"] < 1.0e-14
          and res["DT0_pure_compression_wt"] < 1.0e-6 and res["DT0_pure_tension_wc"] < 1.0e-6
          and res["DT1_sig_maxdiff"] < 1.0e-7 and res["DT1_wt_maxdiff"] < 1.0e-7
          and res["DT1_peak_err"] < 0.02
          and res["DT2_sig_maxdiff"] < 1.0e-2 and res["DT2_wc_maxdiff"] < 1.0e-2
          and res["DT2_peak_err"] < 0.03
          and res["DT3_recovered"]
          and res["DT4_objectivity"] < 1.0e-9)
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As}")
        print(f"  DT0 split partition={res['DT0_split_partition']:.1e}  wt=wc=0 identity={res['DT0_identity']:.1e}"
              f"  unilateral(compr wt-invariant)={res['DT0_unilateral']:.1e}/{res['DT0_compr_wt_invariant']:.1e}"
              f"  pure-compr wt={res['DT0_pure_compression_wt']:.2e} pure-tens wc={res['DT0_pure_tension_wc']:.2e}")
        print(f"  DT1 reduce->P2a (tension): sig maxdiff={res['DT1_sig_maxdiff']:.2e}"
              f"  wt maxdiff={res['DT1_wt_maxdiff']:.2e}  peak={res['DT1_peak']:.3f} (ft={ft})")
        print(f"  DT2 reduce->P2b (compression): sig maxdiff={res['DT2_sig_maxdiff']:.2e}"
              f"  wc maxdiff={res['DT2_wc_maxdiff']:.2e}  peak={res['DT2_peak']:.3f} (fc={fc})")
        print(f"  DT3 unilateral: wt_peak={res['DT3_wt_peak']:.3f}  tension nom/eff(min)={res['DT3_min_ratio_tension']:.3f}"
              f"  early-compression nom/eff(recovery)={res['DT3_recovery_ratio']:.3f}  recovered={res['DT3_recovered']}")
        print(f"  DT4 frame objectivity (rotated damaged stress) = {res['DT4_objectivity']:.2e}")
        print(f"  DT5 [DIAGNOSTIC] tension-after-compression peak = {res['DT5_tension_after_compression_peak']:.3f}"
              f" (fresh ft={ft}; ~0 = un-tempered CDPM2 T/C coupling, P2f beta_c/alpha_t task)")
        print(f"  => P2c GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P2d — the SINGLE-STEP tensor damaged update (what the C++ setTrialStrain mirrors) + its
# numerical DAMAGED CONSISTENT TANGENT (the reference the C++ analytic dual-projector tangent
# is FD-checked against — exactly the P1-tangent slice pattern: oracle = numerical, C++ = analytic).
#
# The path drivers above carry their state in locals; here it is an explicit committed `state` dict
# so a single update can be perturbed for the tangent (fixed committed state, vary the strain
# increment) and chained step-by-step (== the path driver, gate TD0). The damaged consistent tangent
#   dsigma/deps = (1-w_t) dsig_t/deps + (1-w_c) dsig_c/deps  -  sig_t (x) dw_t/deps  -  sig_c (x) dw_c/deps
# is the full dual-projector operator (ADR 4.3 [MAJOR]); its ANALYTIC form (spectral projector
# derivatives dP_T/dsig + the implicit-omega chain rule) is the C++ deliverable. Here it is the
# central difference of `damaged_step_tensor` — the honest reference + the gates that characterize it
# (degraded, INDEFINITE on softening => the Tier-2 IMPL-EX motivation, non-symmetric, frame-objective,
# finite across a reversal AND an eigenvalue-crossing).
# ===========================================================================
def make_damage_state(mp):
    """Committed state for the single-step tensor damaged update (all tensors Voigt {00,11,22,01,12,02},
    tensor shear convention). sigt_max/sigc_max (P2g) are the MONOTONE running maxima of the per-channel
    effective drive stress — they make omega_t/omega_c monotone (no healing on elastic unload)."""
    return dict(sig_bar=np.zeros(6), kp=0.0, eps=np.zeros(6),
                et_max=0.0, kdt1=0.0, kdt2=0.0, kdc=0.0, kdc1=0.0, kdc2=0.0,
                sigt_max=0.0, sigc_max=0.0)


def _plastic_strain6(sig_bar6, eps6, mp):
    """Plastic strain tensor eps_p = eps - C^-1 : sig_bar (tensor Voigt). Generalizes the principal
    compliance the path drivers use; off-diagonal handled via the tensor elastic C."""
    return np.asarray(eps6, float) - np.linalg.solve(elastic_C(mp), np.asarray(sig_bar6, float))


def _dl_beta(mp, dt):
    """Duvaut-Lions relaxation factor beta = dt/(eta+dt) (ADR 4.4). Returns 1.0 (=> inviscid, the
    Tier-1 path BYTE-for-byte) whenever eta<=0 OR dt<=0 — so eta=0 is byte-identical AND a missing
    time increment (static / no pseudo-time) safely falls back to the inviscid return rather than the
    elastic limit beta->0. The viscous regime needs BOTH a positive viscosity and a positive dt."""
    eta = mp.get("eta", 0.0)
    if eta > 0.0 and dt > 0.0:
        return dt / (eta + dt)
    return 1.0


def damaged_step_tensor(state, deps6, mp, Gf, Gc, lch, As=2.0, dt=0.0):
    """ONE constitutive update: committed `state` + strain increment `deps6` -> (nominal sigma[6],
    NEW state, diagnostics). Pure (does not mutate `state`). Identical kinematics to
    drive_damaged_unified (gate TD0), lifted to a general 6-strain via the spectral split.

    Duvaut-Lions viscoplasticity (ADR 4.4, applied at the PLASTIC/effective-stress level): with
    `mp["eta"]>0` and `dt>0` the inviscid effective return `sig_bar` and its hardening variable `kp`
    are RELAXED toward the elastic trial by `beta=dt/(eta+dt)`:
        sig_bar <- (1-beta) sig_tr + beta sig_bar_inviscid ,  kp <- (1-beta) kp_n + beta kp_inviscid
    (the Simo-Hughes closed form: beta->1 as eta->0 recovers the inviscid return EXACTLY; beta->0 as
    eta->inf freezes the elastic trial). Damage is then computed from the RELAXED effective stress, so
    the implied plastic strain (and thus the damage drivers) shrink consistently with the relaxation."""
    E, fc, ft = mp["E"], mp["fc"], mp["ft"]
    eps0 = ft / E
    eps_f = Gf / (ft * lch)
    eps_fc = Gc / (fc * lch)
    deps6 = np.asarray(deps6, float)
    sig_bar, kp_new, plastic, conv = return_map_tensor(state["sig_bar"], deps6, mp, state["kp"])
    beta = _dl_beta(mp, dt)
    if beta < 1.0:                                        # Duvaut-Lions: relax toward the inviscid return
        sig_tr = elastic_pred_tensor(state["sig_bar"], deps6, mp)
        sig_bar = (1.0 - beta) * sig_tr + beta * sig_bar
        kp_new = (1.0 - beta) * state["kp"] + beta * kp_new
    eps_new = state["eps"] + deps6
    w, V = np.linalg.eigh(voigt_to_mat(sig_bar))          # effective principal stresses (ascending)
    et, ac, xs = _damage_drivers(w, mp, As)
    # plastic-strain increment (TENSOR Frobenius norm; reduces to the principal norm when coaxial)
    depl = _plastic_strain6(sig_bar, eps_new, mp) - _plastic_strain6(state["sig_bar"], state["eps"], mp)
    dnorm_epl = float(np.sqrt(np.sum(voigt_to_mat(depl) ** 2)))
    et_max = state["et_max"]
    det_raw = max(et - et_max, 0.0)
    above = max(et - max(et_max, eps0), 0.0)
    loading = det_raw > 0.0 and et > eps0
    kdt1, kdt2 = state["kdt1"], state["kdt2"]
    kdc, kdc1, kdc2 = state["kdc"], state["kdc1"], state["kdc2"]
    if loading:                                           # same accumulation as drive_damaged_unified
        wt_w = tensile_damage_weight(mp, ac, depl, w, V)     # P2h ctTemper weight (1 if 'none')
        kdt2 += wt_w * above / xs
        kdt1 += wt_w * dnorm_epl / xs
        kdc += ac * above
        kdc2 += ac * above / xs
        kdc1 += ac * beta_c(w, kp_new, mp) * dnorm_epl / xs   # Eq.48 with the full CDPM2 beta_c (Eq.50, P2f)
    et_max = max(et_max, et)
    # P2i — multiaxial-consistent TENSILE drive: E*eps_tilde (Eq.37 equivalent strain) instead of the
    # extreme tensile principal, GATED by the presence of a real tensile principal. In uniaxial tension
    # E*eps_tilde == sig_bar_t so it reduces to the extreme principal; in biaxial/triaxial tension
    # eps_tilde folds in ALL tensile principals (E*eps_tilde > the extreme principal) => damage onsets at
    # a lower per-principal stress (the CDPM2-consistent tensile envelope, Eq.37). The COMPRESSIVE drive
    # stays the extreme principal: eps_tilde is ft-scaled (== eps0 on ANY failure surface), so E*eps_tilde
    # can never reach fc and would NEVER onset wc. The tension gate keeps the tensile history clean in
    # pure/dominant compression (wt affects the stress only through the positive-part split anyway).
    # USER decision 2026-06-21 ([[LadrunoConcrete3D_dev_handoff]] §0).
    sig_t_drive = E * et if float(np.max(w)) > 1.0e-6 * ft else 0.0
    sig_c_drive = max(-float(np.min(w)), 0.0)
    # P2g — MONOTONE (no-heal) cyclic damage. Drive each omega with the running MAXIMUM of its effective
    # drive stress, not the live value. The histories kd*1/kd*2 are already monotone (they accumulate only
    # when loading), so the live drive is the ONLY non-monotone input to the bracketed omega-solve: on an
    # elastic UNLOAD the drive drops while the histories are frozen, and F(0)=drive-f*exp(-kd1/eps_f) goes
    # <=0 => the solve relaxes omega back ("heals", the F4 diagnostic). The running max freezes omega on
    # unload (no heal; the nominal stress then unloads along the damage secant (1-omega)*sig_bar) and is
    # BYTE-IDENTICAL on any monotonic path (max == live), so DT1/DT2 reduce-to-P2a/P2b still hold. This is
    # the CDPM2 statement that omega = omega(kappa_d) is a function of the monotone history only.
    sigt_max = max(state.get("sigt_max", 0.0), sig_t_drive)
    sigc_max = max(state.get("sigc_max", 0.0), sig_c_drive)
    # physical floor (see drive_damaged_unified): no omega-solve on a numerical-residual stress
    wt = _solve_omega_bracketed(kdt1, kdt2, sigt_max, ft, eps_f) if (et_max > eps0 and sigt_max > 1.0e-6 * ft) else 0.0
    wc = _solve_omega_bracketed(kdc1, kdc2, sigc_max, fc, eps_fc) if (kdc > 0.0 and sigc_max > 1.0e-6 * fc) else 0.0
    sig_nom = mat_to_voigt(V @ np.diag(apply_damage_principal(w, wt, wc)) @ V.T)   # Eq.1, recompose
    new_state = dict(sig_bar=sig_bar, kp=kp_new, eps=eps_new, et_max=et_max,
                     kdt1=kdt1, kdt2=kdt2, kdc=kdc, kdc1=kdc1, kdc2=kdc2,
                     sigt_max=sigt_max, sigc_max=sigc_max)
    return sig_nom, new_state, dict(wt=wt, wc=wc, plastic=plastic, conv=conv)


def damaged_consistent_tangent(state, deps6, mp, Gf, Gc, lch, As=2.0, rel_step=1.0e-6, dt=0.0):
    """Numerical algorithmic tangent dsigma/deps (6x6) of the damaged update by central difference,
    at fixed COMMITTED state (the algorithmic operator the global Newton consumes). `dt` forwards the
    Duvaut-Lions relaxation (beta is constant in deps, so the FD picks up the blended effective tangent)."""
    base = mp["fc"] / mp["E"]
    C = np.zeros((6, 6))
    for j in range(6):
        d = rel_step * (abs(deps6[j]) + base)
        dp = np.array(deps6, float); dp[j] += d
        dm = np.array(deps6, float); dm[j] -= d
        sp, _, _ = damaged_step_tensor(state, dp, mp, Gf, Gc, lch, As, dt=dt)
        sm, _, _ = damaged_step_tensor(state, dm, mp, Gf, Gc, lch, As, dt=dt)
        C[:, j] = (sp - sm) / (2.0 * d)
    return C


def _advance_damaged(state, eps_path6, mp, Gf, Gc, lch, As=2.0):
    """Chain damaged_step_tensor along a prescribed list of TOTAL strain tensors; return the state +
    the per-step (sigma, wt, wc). Used to build committed softening/reversal states for the tangent."""
    sig, wt, wc = [], [], []
    for eps_t in eps_path6:
        deps = np.asarray(eps_t, float) - state["eps"]
        s, state, info = damaged_step_tensor(state, deps, mp, Gf, Gc, lch, As)
        sig.append(s); wt.append(info["wt"]); wc.append(info["wc"])
    return state, np.array(sig), np.array(wt), np.array(wc)


def run_p2d_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    mp = make_material(E, nu, fc, ft)
    nrm = lambda A: float(np.sqrt(np.sum(A * A)))
    res = {}
    lch = 50.0

    # TD0: the single-step tensor update == the validated P2c path driver. Feed the unified driver's
    # OWN (axial, lateral) strain history through damaged_step_tensor step-by-step and match sig11.
    du = drive_damaged_unified(mp, np.linspace(0, 0.008, 1500), Gf, Gc, lch, As)
    st = make_damage_state(mp)
    sig_step = []
    for i in range(len(du["eps11"])):
        eps_i = np.array([du["eps11"][i], du["eps_lat"][i], du["eps_lat"][i], 0.0, 0.0, 0.0])
        s, st, _ = damaged_step_tensor(st, eps_i - st["eps"], mp, Gf, Gc, lch, As)
        sig_step.append(s[0])
    res["TD0_tension_maxdiff"] = float(np.max(np.abs(np.array(sig_step) - du["sig11"])))
    duc = drive_damaged_unified(mp, np.linspace(0, -0.12, 2000), Gf, Gc, lch, As)
    st = make_damage_state(mp); sigc_step = []
    for i in range(len(duc["eps11"])):
        eps_i = np.array([duc["eps11"][i], duc["eps_lat"][i], duc["eps_lat"][i], 0.0, 0.0, 0.0])
        s, st, _ = damaged_step_tensor(st, eps_i - st["eps"], mp, Gf, Gc, lch, As)
        sigc_step.append(s[0])
    res["TD0_compression_maxdiff"] = float(np.max(np.abs(np.array(sigc_step) - duc["sig11"])))

    # TD1: with NO damage (w=0) the damaged tangent is the P1 EFFECTIVE consistent tangent.
    #   (a) tiny elastic step -> elastic C; (b) a plastic-but-pre-peak step (kappa_p<1, w=0) -> the
    #   effective P1 tangent (consistent_tangent on the same trial). Confirms the (1-w) factor and the
    #   -sig (x) dw rank-update both vanish before onset.
    st0 = make_damage_state(mp)
    Cel = elastic_C(mp)
    C_el = damaged_consistent_tangent(st0, np.array([-1.0e-6, 0, 0, 0, 0, 0]), mp, Gf, Gc, lch, As)
    res["TD1a_elastic_err"] = nrm(C_el - Cel) / nrm(Cel)
    deps_pp = np.array([-9.0e-4, 2.0e-4, 2.0e-4, 0.0, 0.0, 0.0])        # plastic, pre-peak (kp<1)
    _, _, info_pp = damaged_step_tensor(st0, deps_pp, mp, Gf, Gc, lch, As)
    res["TD1b_predamage"] = bool(info_pp["wt"] == 0.0 and info_pp["wc"] == 0.0 and info_pp["plastic"])
    C_dmg_pp = damaged_consistent_tangent(st0, deps_pp, mp, Gf, Gc, lch, As)
    C_eff_pp = consistent_tangent(np.zeros(6), deps_pp, mp, 0.0, hardening=True)
    res["TD1b_match_effective"] = nrm(C_dmg_pp - C_eff_pp) / nrm(C_eff_pp)

    # TD2: in the SOFTENING regime the damaged tangent is DEGRADED and INDEFINITE. Drive a prescribed
    # uniaxial-STRAIN tension path just past the nominal peak, freeze that committed state, take the
    # axial tangent for a further tension increment: C[0,0] < 0 (snap-back) and lambda_min(sym C) < 0.
    path = [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(0, 6.0e-4, 400)]
    sstate, sig_path, wt_path, _ = _advance_damaged(make_damage_state(mp), path, mp, Gf, Gc, lch, As)
    ipk = int(np.argmax(sig_path[:, 0]))
    res["TD2_softening_reached"] = bool(wt_path[-1] > 0.5 and sig_path[-1, 0] < sig_path[ipk, 0])
    # rebuild the committed state at a post-peak point (2/3 into the softening tail)
    isoft = (ipk + len(path)) // 2 + (len(path) - ipk) // 6
    sstate2, _, wt2, _ = _advance_damaged(make_damage_state(mp), path[:isoft], mp, Gf, Gc, lch, As)
    C_soft = damaged_consistent_tangent(sstate2, np.array([2.0e-7, 0, 0, 0, 0, 0]), mp, Gf, Gc, lch, As)
    res["TD2_wt_at_state"] = float(wt2[-1])
    res["TD2_axial_tangent"] = float(C_soft[0, 0])
    sym = 0.5 * (C_soft + C_soft.T)
    res["TD2_min_eig_sym"] = float(np.min(np.linalg.eigvalsh(sym)))
    res["TD2_degraded_indefinite"] = bool(C_soft[0, 0] < 0.0 and res["TD2_min_eig_sym"] < 0.0)

    # TD3: NON-symmetric — a damaged, sheared, non-associated state. Both the non-associated plastic
    # flow AND the -sig (x) dw damage rank-update break symmetry => unsymmetric solver (ADR 4.4).
    mp_na = make_material(E, nu, fc, ft, Df=0.3)
    st_sh = make_damage_state(mp_na)
    pre = [np.array([e, -0.2 * e, -0.2 * e, 0.3 * e, 0, 0]) for e in np.linspace(0, 5.0e-4, 300)]
    st_sh, _, wts, _ = _advance_damaged(st_sh, pre, mp_na, Gf, Gc, lch, As)
    C_sh = damaged_consistent_tangent(st_sh, np.array([3.0e-6, -0.6e-6, -0.6e-6, 0.9e-6, 0, 0]),
                                      mp_na, Gf, Gc, lch, As)
    res["TD3_wt"] = float(wts[-1])
    res["TD3_asym"] = nrm(C_sh - C_sh.T) / nrm(C_sh)

    # TD4: frame objectivity of the UPDATE (state + increment rotated by Q -> stress rotates by Q).
    # Exercises the histories (invariant) + the spectral recompose under a non-trivial committed state.
    base_state = make_damage_state(mp)
    bp = [np.array([e, -0.1 * e, 0.05 * e, 0.2 * e, -0.1 * e, 0.05 * e]) for e in np.linspace(0, 4.0e-4, 200)]
    base_state, _, _, _ = _advance_damaged(base_state, bp, mp, Gf, Gc, lch, As)
    deps_o = np.array([2.0e-6, -0.5e-6, 0.3e-6, 0.4e-6, -0.2e-6, 0.1e-6])
    s_base, _, _ = damaged_step_tensor(base_state, deps_o, mp, Gf, Gc, lch, As)
    th = 0.6
    Q = np.array([[np.cos(th), -np.sin(th), 0.0], [np.sin(th), np.cos(th), 0.0], [0.0, 0.0, 1.0]])
    rot_state = dict(base_state)
    rot_state["sig_bar"] = mat_to_voigt(Q @ voigt_to_mat(base_state["sig_bar"]) @ Q.T)
    rot_state["eps"] = mat_to_voigt(Q @ voigt_to_mat(base_state["eps"]) @ Q.T)
    deps_rot = mat_to_voigt(Q @ voigt_to_mat(deps_o) @ Q.T)
    s_rot, _, _ = damaged_step_tensor(rot_state, deps_rot, mp, Gf, Gc, lch, As)
    res["TD4_objectivity"] = nrm(voigt_to_mat(s_rot) - Q @ voigt_to_mat(s_base) @ Q.T) / nrm(voigt_to_mat(s_base))

    # TD5: the HARD cases stay finite (no NaN / blow-up). (i) load REVERSAL — committed state is
    # tension-damaged (wt high), tangent for a COMPRESSION increment (the (1-wt) channel switches off,
    # the (1-wc) channel switches on). (ii) near-DEGENERATE eigenvalues (sig_bar ~ hydrostatic) where
    # the analytic dP_T/dsig term is singular — the numerical tangent must still be finite (the C++
    # analytic tangent must regularize the eigenvalue-crossing or accept Tier-2 drops; ADR 4.3 [MAJOR]).
    rev_state, _, wt_rev, _ = _advance_damaged(make_damage_state(mp),
                                               [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(0, 5.0e-4, 300)],
                                               mp, Gf, Gc, lch, As)
    C_rev = damaged_consistent_tangent(rev_state, np.array([-2.0e-6, 0, 0, 0, 0, 0]), mp, Gf, Gc, lch, As)
    res["TD5_reversal_wt"] = float(wt_rev[-1])
    res["TD5_reversal_finite"] = bool(np.all(np.isfinite(C_rev)) and nrm(C_rev) < 1.0e3 * nrm(Cel))
    deg_state = make_damage_state(mp)
    deg_state["sig_bar"] = np.array([2.0, 2.0, 2.0 + 1.0e-7, 0.0, 0.0, 0.0])     # near-triple eigenvalue
    C_deg = damaged_consistent_tangent(deg_state, np.array([1.0e-6, 1.0e-6, -2.0e-6, 0, 0, 0]), mp, Gf, Gc, lch, As)
    res["TD5_degenerate_finite"] = bool(np.all(np.isfinite(C_deg)))

    # TD0 tension is exact (1e-14); the compression residual (~5e-9 MPa, ~2e-10 relative to fc) is the
    # eigendecompose route (return_map_tensor) vs the direct principal path driver accumulated over the
    # long deep-compression path — numerically "the same stress", with a platform-robust floor for the
    # LAPACK eigensolver. The genuine physics-equivalence is the tension leg + TD1's tangent reductions.
    ok = (res["TD0_tension_maxdiff"] < 1.0e-9 and res["TD0_compression_maxdiff"] < 1.0e-6
          and res["TD1a_elastic_err"] < 1.0e-6
          and res["TD1b_predamage"] and res["TD1b_match_effective"] < 1.0e-6
          and res["TD2_softening_reached"] and res["TD2_degraded_indefinite"]
          and res["TD3_wt"] > 0.05 and res["TD3_asym"] > 1.0e-2
          and res["TD4_objectivity"] < 1.0e-9
          and res["TD5_reversal_finite"] and res["TD5_degenerate_finite"])
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As}")
        print(f"  TD0 single-step==path-driver: tension={res['TD0_tension_maxdiff']:.2e}"
              f"  compression={res['TD0_compression_maxdiff']:.2e}")
        print(f"  TD1a elastic tangent == C_el = {res['TD1a_elastic_err']:.2e}")
        print(f"  TD1b pre-damage (w=0,plastic)={res['TD1b_predamage']}  tangent==P1 effective:"
              f" {res['TD1b_match_effective']:.2e}")
        print(f"  TD2 softening: wt@state={res['TD2_wt_at_state']:.3f}  C[0,0]={res['TD2_axial_tangent']:.1f}"
              f"  lambda_min(sym)={res['TD2_min_eig_sym']:.1f}  degraded+indefinite={res['TD2_degraded_indefinite']}")
        print(f"  TD3 non-symmetric (Df=0.3, damaged, shear): wt={res['TD3_wt']:.3f}"
              f"  ||C-C^T||/||C||={res['TD3_asym']:.3e}")
        print(f"  TD4 frame objectivity of the update = {res['TD4_objectivity']:.2e}")
        print(f"  TD5 reversal finite (wt={res['TD5_reversal_wt']:.2f})={res['TD5_reversal_finite']}"
              f"  near-degenerate-eig finite={res['TD5_degenerate_finite']}")
        print(f"  => P2d GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P2e — the ANALYTIC dual-projector DAMAGED CONSISTENT TANGENT (the C++ deliverable), FD-verified
# against the P2d numerical reference (`damaged_consistent_tangent`). Structure (ADR 4.3 [MAJOR]):
#     dsigma/deps = D_dam : C_eff  -  sig_t (x) dw_t/deps  -  sig_c (x) dw_c/deps
#   * C_eff = the P1 EFFECTIVE consistent tangent dsig_bar/deps (numerical here; ANALYTIC in the C++
#     kernel #249). P2e adds the analytic DAMAGE linearization on top.
#   * D_dam = the spectral derivative of the per-principal damaged stress with omega FROZEN
#     (de Souza Neto Box A.6 isotropic tensor-function derivative; the same machinery as the P1
#     spectral tangent). This IS the (1-w_t)dsig_t/deps + (1-w_c)dsig_c/deps dual-projector secant.
#   * dw_t/deps, dw_c/deps via the implicit-function theorem on the bracketed omega-solve, chained
#     through the histories. The scalar sub-gradients d(eps_tilde)/dsig, d(x_s)/dsig, d(alpha_c)/dsig
#     are ISOLATED scalar micro-FDs (the LadrunoJ2/P1 "Lode directional gradient by micro-FD" pattern);
#     d(lambda_extreme)/dsig is the analytic eigenprojection; d||deps_p||/deps is closed form.
#
# KEY non-smooth point (NOT a bug): at sigma_lat = 0 (uniaxial-STRESS compression) the lateral
# eigenvalues sit on the Macaulay kink, so the central-difference tangent crosses into tension and
# disagrees with the analytic subgradient by O(w_c-w_t) on the ~zero-stress lateral directions. The
# analytic tangent returns a VALID subgradient (the committed-sign side); the gate verifies it at
# smooth states (eigenvalues bounded from 0) and characterizes the kink on the loaded component only.
# ===========================================================================
_TENSOR_W6 = np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0])   # double-contraction weights, tensor-shear Voigt


def ddot6(a, b):
    """Tensor double contraction A:B for symmetric tensors in Voigt {00,11,22,01,12,02} (tensor shear)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    return float(a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + 2.0*(a[3]*b[3] + a[4]*b[4] + a[5]*b[5]))


def isotropic_tangent(lam, V, yv, ypv, tol=1.0e-8):
    """6x6 derivative dY/dX of an isotropic symmetric-tensor function Y = sum_a y(lam_a) E_a, given the
    eigenvalues `lam`, eigenvectors `V` (columns), and the per-eigenvalue values `yv`=y(lam_a) and
    `ypv`=y'(lam_a). de Souza Neto Box A.6 operational form:
        (dY/dX : S) = sum_a ypv_a (E_a:S) E_a + sum_{a!=b} G_ab E_a S E_b,
        G_ab = (yv_a - yv_b)/(lam_a - lam_b)   (-> ypv_a as lam_b -> lam_a, l'Hopital)."""
    E = [np.outer(V[:, a], V[:, a]) for a in range(3)]

    def apply(S):
        out = np.zeros((3, 3))
        for a in range(3):
            out += ypv[a] * np.sum(E[a] * S) * E[a]
        for a in range(3):
            for b in range(3):
                if a == b:
                    continue
                dl = lam[a] - lam[b]
                G = ypv[a] if abs(dl) < tol else (yv[a] - yv[b]) / dl
                out += G * (E[a] @ S @ E[b])
        return out

    D = np.zeros((6, 6))
    for j in range(6):
        Sj = np.zeros(6); Sj[j] = 1.0
        D[:, j] = mat_to_voigt(apply(voigt_to_mat(Sj)))
    return D


def _dscalar_dsig(fn, sig6, h=1.0e-6):
    """Isolated micro-FD gradient d(fn)/d(sig) of a scalar-of-stress function (per Voigt component)."""
    g = np.zeros(6)
    for k in range(6):
        dp = np.array(sig6, float); dp[k] += h
        dm = np.array(sig6, float); dm[k] -= h
        g[k] = (fn(dp) - fn(dm)) / (2.0 * h)
    return g


def damaged_tangent_analytic(state, deps6, mp, Gf, Gc, lch, As=2.0, dt=0.0):
    """ANALYTIC dual-projector damaged consistent tangent (6x6). Recomputes the same step as
    `damaged_step_tensor` (so it is self-contained) and assembles D_dam:C_eff - sig_t(x)dw_t -
    sig_c(x)dw_c. FD-verified against `damaged_consistent_tangent` at smooth states (run_p2e_gate).

    Duvaut-Lions (ADR 4.4): with `mp["eta"]>0` and `dt>0` the effective stress is relaxed exactly as in
    `damaged_step_tensor`, and the effective consistent tangent blends `C_eff <- (1-beta)C_elastic +
    beta C_eff_inviscid` (since beta is constant in deps and sig_tr is linear in deps). The damage
    linearization then chains through the RELAXED effective stress and its blended tangent automatically."""
    E, fc, ft = mp["E"], mp["fc"], mp["ft"]
    eps0 = ft / E
    eps_f = Gf / (ft * lch)
    eps_fc = Gc / (fc * lch)
    deps6 = np.asarray(deps6, float)
    sig_bar, kp_new, plastic, conv = return_map_tensor(state["sig_bar"], deps6, mp, state["kp"])
    beta = _dl_beta(mp, dt)
    if beta < 1.0:                                        # relax the effective stress + kp (mirror the update)
        sig_tr = elastic_pred_tensor(state["sig_bar"], deps6, mp)
        sig_bar = (1.0 - beta) * sig_tr + beta * sig_bar
        kp_new = (1.0 - beta) * state["kp"] + beta * kp_new
    eps_new = state["eps"] + deps6
    lam, V = np.linalg.eigh(voigt_to_mat(sig_bar))
    et, ac, xs = _damage_drivers(lam, mp, As)
    depl = _plastic_strain6(sig_bar, eps_new, mp) - _plastic_strain6(state["sig_bar"], state["eps"], mp)
    dnorm = float(np.sqrt(ddot6(depl, depl)))
    et_max = state["et_max"]
    det_raw = max(et - et_max, 0.0)
    above = max(et - max(et_max, eps0), 0.0)
    loading = det_raw > 0.0 and et > eps0
    kdt1, kdt2 = state["kdt1"], state["kdt2"]
    kdc, kdc1, kdc2 = state["kdc"], state["kdc1"], state["kdc2"]
    bc = beta_c(lam, kp_new, mp)                          # Eq.50 (P2f): scales the kdc1 plastic part
    wt_w = tensile_damage_weight(mp, ac, depl, lam, V)    # P2h ctTemper weight
    if loading:
        kdt2 += wt_w * above / xs; kdt1 += wt_w * dnorm / xs
        kdc += ac * above; kdc2 += ac * above / xs; kdc1 += ac * bc * dnorm / xs
    et_max2 = max(et_max, et)
    Dt = E * et if float(np.max(lam)) > 1.0e-6 * ft else 0.0   # P2i: E*eps_tilde tensile drive (Eq.37)
    Dc = max(-float(np.min(lam)), 0.0)                          # compressive drive: extreme principal
    # P2g — MONOTONE drive (mirror damaged_step_tensor). Solve omega against the running max, and mark
    # whether each channel is ADVANCING its max (== loading). On UNLOAD (live drive < committed max) the
    # drive is frozen AND the histories are frozen, so d(omega)/d(eps)=0 and the tangent collapses to the
    # damage SECANT D_dam:C_eff (SPD) — the well-conditioned unloading branch. During loading max == live
    # and the channel advances, so the tangent is byte-identical to the pre-P2g analytic tangent (the P2e
    # gates are unaffected).
    sigt_max = max(state.get("sigt_max", 0.0), Dt)
    sigc_max = max(state.get("sigc_max", 0.0), Dc)
    t_loading = Dt >= state.get("sigt_max", 0.0)          # live tensile drive at/above the committed max
    c_loading = Dc >= state.get("sigc_max", 0.0)
    # SAME physical floor as damaged_step_tensor (keeps the analytic tangent's omega == the update's)
    wt = _solve_omega_bracketed(kdt1, kdt2, sigt_max, ft, eps_f) if (et_max2 > eps0 and sigt_max > 1.0e-6 * ft) else 0.0
    wc = _solve_omega_bracketed(kdc1, kdc2, sigc_max, fc, eps_fc) if (kdc > 0.0 and sigc_max > 1.0e-6 * fc) else 0.0

    Ceff = consistent_tangent(state["sig_bar"], deps6, mp, state["kp"], hardening=True)
    if beta < 1.0:                                        # Duvaut-Lions blended effective tangent
        Ceff = (1.0 - beta) * elastic_C(mp) + beta * Ceff

    # D_dam: spectral derivative of the per-principal damage with omega FROZEN
    yv = [(1.0 - wt) * max(lam[a], 0.0) + (1.0 - wc) * min(lam[a], 0.0) for a in range(3)]
    ypv = [(1.0 - wt) if lam[a] > 0.0 else (1.0 - wc) for a in range(3)]
    Ddam = isotropic_tangent(lam, V, yv, ypv)
    sig_t = mat_to_voigt(V @ np.diag(np.maximum(lam, 0.0)) @ V.T)
    sig_c = mat_to_voigt(V @ np.diag(np.minimum(lam, 0.0)) @ V.T)
    C = Ddam @ Ceff

    # chain-rule pieces (each a 6-vector d(.)/deps)
    Cinv = np.linalg.inv(elastic_C(mp))
    depl_deps = np.eye(6) - Cinv @ Ceff                                  # d(eps_p)/d(eps)
    det_deps = Ceff.T @ _dscalar_dsig(lambda s: equiv_strain_general(np.linalg.eigvalsh(voigt_to_mat(s)), mp), sig_bar)
    dxs_deps = Ceff.T @ _dscalar_dsig(lambda s: _damage_drivers(np.linalg.eigvalsh(voigt_to_mat(s)), mp, As)[2], sig_bar)
    dac_deps = Ceff.T @ _dscalar_dsig(lambda s: alpha_compression(np.linalg.eigvalsh(voigt_to_mat(s))), sig_bar)
    Emax = np.outer(V[:, 2], V[:, 2]); Emin = np.outer(V[:, 0], V[:, 0])  # eigh ascending
    # d(drive_max)/d(eps): the live eigenprojection ONLY while the channel advances its running max
    # (P2g); frozen (zero) on unload so the -sig(x)d(omega) rank-update vanishes => secant tangent.
    # P2i: d(E*eps_tilde)/d(eps) = E * d(eps_tilde)/d(eps) = E * det_deps (the equiv-strain gradient
    # already assembled below), REPLACING the extreme-principal eigenprojection (Emax) for the tensile
    # channel. The compressive channel keeps the extreme-principal eigenprojection (Emin).
    dDt_deps = E * det_deps if (Dt > 0.0 and t_loading) else np.zeros(6)
    dDc_deps = -(Ceff.T @ (_TENSOR_W6 * mat_to_voigt(Emin))) if (Dc > 0.0 and c_loading) else np.zeros(6)
    dnorm_deps = (depl_deps.T @ (_TENSOR_W6 * depl)) / dnorm if dnorm > 1.0e-14 else np.zeros(6)
    # d(beta_c)/d(eps): beta_c = beta_c(sig_bar(eps), kp(eps)) depends on BOTH rho_bar(sig_bar) AND
    # qh2(kp), so a single composite micro-FD THROUGH the return map captures the full gradient (the
    # rho_bar part is smooth; the qh2(kp) part needs d kp/d eps, which the C++ kernel gets analytically
    # from the 4-unknown return-map IFT — here it is the legitimate oracle micro-FD, same step as the
    # numerical reference so the FD truncation correlates). Only needed under compressive loading.
    dbc_deps = np.zeros(6)
    if loading and bc > 0.0:
        base = fc / E

        def _bc_of(d6):
            sb, kpv, _, _ = return_map_tensor(state["sig_bar"], d6, mp, state["kp"])
            if beta < 1.0:
                sb = (1.0 - beta) * elastic_pred_tensor(state["sig_bar"], d6, mp) + beta * sb
                kpv = (1.0 - beta) * state["kp"] + beta * kpv
            return beta_c(np.linalg.eigvalsh(voigt_to_mat(sb)), kpv, mp)

        for j in range(6):
            hh = 1.0e-6 * (abs(deps6[j]) + base)
            dp = np.array(deps6, float); dp[j] += hh
            dm = np.array(deps6, float); dm[j] -= hh
            dbc_deps[j] = (_bc_of(dp) - _bc_of(dm)) / (2.0 * hh)

    # d(w_t)/d(eps) for the ctTemper modes (P2h): 'none' -> 0 (unchanged tangent); 'alphat' ->
    # -d(alpha_c)/deps (analytic, reuses dac_deps); 'proj' -> composite micro-FD through the return map
    # (w_t is the tensile fraction of the plastic-strain increment). Only under loading.
    mode = mp.get("ct_temper", "none")
    dwt_w_deps = np.zeros(6)
    if loading and mode == "alphat":
        dwt_w_deps = -dac_deps
    elif loading and mode == "proj":
        base = fc / E

        def _wtw_of(d6):
            sb, _, _, _ = return_map_tensor(state["sig_bar"], d6, mp, state["kp"])
            if beta < 1.0:
                sb = (1.0 - beta) * elastic_pred_tensor(state["sig_bar"], d6, mp) + beta * sb
            dpl = _plastic_strain6(sb, state["eps"] + d6, mp) - _plastic_strain6(state["sig_bar"], state["eps"], mp)
            wv, Vv = np.linalg.eigh(voigt_to_mat(sb))
            return tensile_damage_weight(mp, ac, dpl, wv, Vv)

        for j in range(6):
            hh = 1.0e-6 * (abs(deps6[j]) + base)
            dp = np.array(deps6, float); dp[j] += hh
            dm = np.array(deps6, float); dm[j] -= hh
            dwt_w_deps[j] = (_wtw_of(dp) - _wtw_of(dm)) / (2.0 * hh)

    if loading:
        # kdt2 = w_t * above / xs ; kdt1 = w_t * dnorm / xs  (Eq.45/44 with the ctTemper weight) => product rule
        dkdt2 = wt_w * (det_deps / xs - above * dxs_deps / xs**2) + (above / xs) * dwt_w_deps
        dkdt1 = wt_w * (dnorm_deps / xs - dnorm * dxs_deps / xs**2) + (dnorm / xs) * dwt_w_deps
        dkdc2 = (dac_deps * above + ac * det_deps) / xs - ac * above * dxs_deps / xs**2
        # kdc1 = ac * bc * dnorm / xs  (Eq.48 with beta_c) => product rule over ac, bc, dnorm, xs
        dkdc1 = (dac_deps * bc * dnorm + ac * dbc_deps * dnorm + ac * bc * dnorm_deps) / xs \
            - ac * bc * dnorm * dxs_deps / xs**2
    else:
        dkdt2 = dkdt1 = dkdc2 = dkdc1 = np.zeros(6)

    # omega via IFT on F(w) = (1-w)D - f exp(-(kd1+w kd2)/eps_f) = 0 ; H = dF/dw = D[(1-w)kd2/eps_f - 1]
    # (P2g) D = the MONOTONE drive sigt_max/sigc_max; dDt_deps/dDc_deps are already zeroed on unload, so
    # an unloading channel contributes d(omega)=0 (secant). On loading D == live drive => unchanged.
    if 0.0 < wt < 1.0:
        Ht = sigt_max * ((1.0 - wt) * kdt2 / eps_f - 1.0)
        dwt = (-(1.0 - wt) / Ht) * dDt_deps + (-(1.0 - wt) * sigt_max / (eps_f * Ht)) * dkdt1 \
            + (-(1.0 - wt) * sigt_max * wt / (eps_f * Ht)) * dkdt2
    else:
        dwt = np.zeros(6)
    if 0.0 < wc < 1.0:
        Hc = sigc_max * ((1.0 - wc) * kdc2 / eps_fc - 1.0)
        dwc = (-(1.0 - wc) / Hc) * dDc_deps + (-(1.0 - wc) * sigc_max / (eps_fc * Hc)) * dkdc1 \
            + (-(1.0 - wc) * sigc_max * wc / (eps_fc * Hc)) * dkdc2
    else:
        dwc = np.zeros(6)

    return C - np.outer(sig_t, dwt) - np.outer(sig_c, dwc)


def run_p2e_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    mp = make_material(E, nu, fc, ft)
    nrm = lambda A: float(np.sqrt(np.sum(A * A)))
    rel = lambda Ca, Cn: nrm(Ca - Cn) / nrm(Cn)
    lch = 50.0
    res = {}

    def tangents(state, deps, m=mp):
        return (damaged_tangent_analytic(state, deps, m, Gf, Gc, lch, As),
                damaged_consistent_tangent(state, deps, m, Gf, Gc, lch, As))

    # PE0: the spectral isotropic-function derivative is correct vs a numerical dY/dX — for the known
    # Y=X^2 (distinct eigvals) AND the damage function (incl. a near-degenerate eigenvalue, l'Hopital).
    def Yof(X6, y):
        l, Vv = np.linalg.eigh(voigt_to_mat(X6)); return mat_to_voigt(Vv @ np.diag([y(z) for z in l]) @ Vv.T)
    def numD(X6, y, h=1.0e-6):
        D = np.zeros((6, 6))
        for j in range(6):
            dp = X6.copy(); dp[j] += h; dm = X6.copy(); dm[j] -= h
            D[:, j] = (Yof(dp, y) - Yof(dm, y)) / (2.0 * h)
        return D
    def specD(X6, y, yp):
        l, Vv = np.linalg.eigh(voigt_to_mat(X6))
        return isotropic_tangent(l, Vv, [y(z) for z in l], [yp(z) for z in l])
    Xd = np.array([3.0, -1.0, 2.0, 0.4, -0.2, 0.3])
    res["PE0_sq"] = nrm(specD(Xd, lambda z: z*z, lambda z: 2*z) - numD(Xd, lambda z: z*z))
    yd = lambda z: 0.4*max(z, 0.0) + 0.8*min(z, 0.0); ypd = lambda z: 0.4 if z > 0 else 0.8
    res["PE0_damage"] = nrm(specD(Xd, yd, ypd) - numD(Xd, yd))
    Xdeg = np.array([2.0, 2.0 + 1.0e-7, -1.0, 0.0, 0.0, 0.0])
    res["PE0_degenerate"] = nrm(specD(Xdeg, yd, ypd) - numD(Xdeg, yd, h=1.0e-7))

    # PE1: TENSION-damaged — analytic == numerical reference.
    st_t, _, wt_t, _ = _advance_damaged(make_damage_state(mp),
                                        [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(0, 4.5e-4, 300)],
                                        mp, Gf, Gc, lch, As)
    Ca, Cn = tangents(st_t, np.array([1.0e-6, 0, 0, 0, 0, 0]))
    res["PE1_tension_rel"] = rel(Ca, Cn); res["PE1_wt"] = float(wt_t[-1])

    # PE2: CONFINED COMPRESSION-damaged (smooth — lateral effective stress held at -sigma3 so the
    # principals stay clearly compressive, OFF the sigma_lat=0 Macaulay kink, AND alpha_c<1 so the
    # d(alpha_c)/dsig term is exercised). Built STRESS-controlled (cross-platform-robust, unlike the
    # confined-STRAIN return-map regime which is chaotic near the apex); the committed state is taken
    # at the first step past wc>0.5 so it is well inside softening with margin.
    dconf = drive_damaged_unified(mp, np.linspace(0, -0.05, 2000), Gf, Gc, lch, As, sigma3=0.05 * fc)
    ic = int(np.argmax(dconf["wc"] > 0.5)) if (dconf["wc"] > 0.5).any() else -1
    st_c, _, _, _ = _advance_damaged(make_damage_state(mp),
                                     [np.array([dconf["eps11"][i], dconf["eps_lat"][i], dconf["eps_lat"][i], 0, 0, 0]) for i in range(ic)],
                                     mp, Gf, Gc, lch, As)
    lam_c = np.linalg.eigvalsh(voigt_to_mat(st_c["sig_bar"]))
    dc = np.array([dconf["eps11"][ic], dconf["eps_lat"][ic], dconf["eps_lat"][ic], 0, 0, 0]) - st_c["eps"]
    Ca, Cn = tangents(st_c, dc)
    res["PE2_compression_rel"] = rel(Ca, Cn); res["PE2_wc"] = float(dconf["wc"][ic]) if ic > 0 else 0.0
    res["PE2_lam_max"] = float(lam_c.max())

    # PE3: SHEAR, NON-ASSOCIATED (Df=0.3), damaged — exercises the off-diagonal spectral terms +
    # the non-symmetric coupling.
    mp_na = make_material(E, nu, fc, ft, Df=0.3)
    st_s, _, wt_s, _ = _advance_damaged(make_damage_state(mp_na),
                                        [np.array([e, -0.2*e, -0.2*e, 0.3*e, 0, 0]) for e in np.linspace(0, 5.0e-4, 300)],
                                        mp_na, Gf, Gc, lch, As)
    Ca, Cn = tangents(st_s, np.array([3.0e-6, -0.6e-6, -0.6e-6, 0.9e-6, 0, 0]), m=mp_na)
    res["PE3_shear_rel"] = rel(Ca, Cn); res["PE3_wt"] = float(wt_s[-1])

    # PE4: load REVERSAL — committed tension-damaged, tangent for a compression increment.
    Ca, Cn = tangents(st_t, np.array([-2.0e-6, 0, 0, 0, 0, 0]))
    res["PE4_reversal_rel"] = rel(Ca, Cn)

    # PE5: reduce-to — no damage => analytic tangent is the elastic C / the P1 effective tangent.
    res["PE5_elastic_rel"] = nrm(damaged_tangent_analytic(make_damage_state(mp), np.array([-1.0e-6, 0, 0, 0, 0, 0]),
                                                          mp, Gf, Gc, lch, As) - elastic_C(mp)) / nrm(elastic_C(mp))
    deps_pp = np.array([-9.0e-4, 2.0e-4, 2.0e-4, 0, 0, 0])
    _, _, info_pp = damaged_step_tensor(make_damage_state(mp), deps_pp, mp, Gf, Gc, lch, As)
    Ca = damaged_tangent_analytic(make_damage_state(mp), deps_pp, mp, Gf, Gc, lch, As)
    Ceff_pp = consistent_tangent(np.zeros(6), deps_pp, mp, 0.0, hardening=True)
    res["PE5_predamage_rel"] = rel(Ca, Ceff_pp); res["PE5_predamage_w0"] = bool(info_pp["wt"] == 0.0 and info_pp["wc"] == 0.0)

    # PE6: the Macaulay-KINK at sigma_lat=0 (uniaxial-STRESS compression, the lateral eigenvalues
    # degenerate near zero) is a valid-subgradient point, NOT a bug. The analytic tangent matches the
    # numerical on the LOADED axial column (gated below); the disagreement (~33%) is in the columns
    # whose perturbation crosses the kink — the near-zero-stress LATERAL-normal directions AND the
    # coupled in-plane SHEAR (the central diff straddles the t/c boundary there). PE1-PE5 gate the FULL
    # 6x6 at smooth states (eigenvalues bounded from 0, incl. confined compression PE2); PE6 only
    # asserts the axial column + REPORTS the kink spread.
    duc = drive_damaged_unified(mp, np.linspace(0, -0.12, 2000), Gf, Gc, lch, As)
    ik = int(np.argmax(duc["wc"] > 0.5))
    st_k, _, _, _ = _advance_damaged(make_damage_state(mp),
                                     [np.array([duc["eps11"][i], duc["eps_lat"][i], duc["eps_lat"][i], 0, 0, 0]) for i in range(ik)],
                                     mp, Gf, Gc, lch, As)
    # the PHYSICAL next increment (lateral expanding to hold sigma_lat=0) — perturbing the near-zero
    # lateral directions is what crosses the Macaulay kink under the central difference.
    dk = np.array([duc["eps11"][ik], duc["eps_lat"][ik], duc["eps_lat"][ik], 0, 0, 0]) - st_k["eps"]
    Ca = damaged_tangent_analytic(st_k, dk, mp, Gf, Gc, lch, As)
    Cn = damaged_consistent_tangent(st_k, dk, mp, Gf, Gc, lch, As)
    res["PE6_kink_axial_rel"] = abs(Ca[0, 0] - Cn[0, 0]) / abs(Cn[0, 0])      # loaded axial component agrees
    res["PE6_kink_full_rel"] = rel(Ca, Cn)                                     # full tangent: the kink shows here

    TOL = 1.0e-6
    ok = (res["PE0_sq"] < 1.0e-5 and res["PE0_damage"] < 1.0e-7 and res["PE0_degenerate"] < 1.0e-7
          and res["PE1_tension_rel"] < TOL and res["PE1_wt"] > 0.5
          and res["PE2_compression_rel"] < TOL and res["PE2_wc"] > 0.3 and res["PE2_lam_max"] < -0.5
          and res["PE3_shear_rel"] < TOL and res["PE3_wt"] > 0.5
          and res["PE4_reversal_rel"] < TOL
          and res["PE5_elastic_rel"] < 1.0e-9 and res["PE5_predamage_w0"] and res["PE5_predamage_rel"] < TOL
          and res["PE6_kink_axial_rel"] < 1.0e-3)
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As}")
        print(f"  PE0 spectral dY/dX: X^2={res['PE0_sq']:.2e}  damage={res['PE0_damage']:.2e}"
              f"  near-degenerate(l'Hopital)={res['PE0_degenerate']:.2e}")
        print(f"  PE1 TENSION (wt={res['PE1_wt']:.3f}) analytic==numerical: rel={res['PE1_tension_rel']:.2e}")
        print(f"  PE2 CONFINED-COMPRESSION (wc={res['PE2_wc']:.3f}, lam_max={res['PE2_lam_max']:.2f}): rel={res['PE2_compression_rel']:.2e}")
        print(f"  PE3 SHEAR non-assoc (wt={res['PE3_wt']:.3f}, Df=0.3): rel={res['PE3_shear_rel']:.2e}")
        print(f"  PE4 REVERSAL (tension-damaged + compression incr): rel={res['PE4_reversal_rel']:.2e}")
        print(f"  PE5 reduce: elastic={res['PE5_elastic_rel']:.2e}  pre-damage(w=0)={res['PE5_predamage_w0']}"
              f" tangent==P1-effective rel={res['PE5_predamage_rel']:.2e}")
        print(f"  PE6 Macaulay-kink (uniaxial-stress compression sigma_lat=0): axial rel={res['PE6_kink_axial_rel']:.2e}"
              f" (matches); full rel={res['PE6_kink_full_rel']:.2e} (kink on ~zero-stress lateral — valid subgradient)")
        print(f"  => P2e GATE {'PASS' if ok else 'FAIL'}")
    return res


def run_p2f_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    """P2f — the CDPM2 cyclic beta_c (Eq.50), restored into the compressive-damage plastic driver
    kappa_dc1 (Eq.48). beta_c = ft*qh2(kp)*sqrt(2/3) / (rho_bar*sqrt(1+2*Df^2)) gives "a smooth
    transition from pure damage to damage-plasticity during cyclic loading" (Grassl 2013 2.3.5); in
    MONOTONIC compression it is ~ft/(fc*sqrt(1+2*Df^2)) << 1, so it makes compression markedly MORE
    DUCTILE than the beta_c=1 simplification the P2b/P2c monotonic slice used (the chosen 'faithful'
    direction). Gates: F1 the closed-form beta_c at the uniaxial-compression peak + bounds; F2 the
    monotonic backbone STILL valid (peak=fc, softens, effective-stress monotone); F3 NON-tautology — the
    real beta_c suppresses damage vs beta_c=1 (more ductile) AND is a genuine state-dependent factor in
    (0,1); F4 cyclic compression load->unload->reload is consistent (no damage healing, reload follows
    the degraded secant, reload stays at/under the monotonic envelope)."""
    mp = make_material(E, nu, fc, ft); Df = mp["Df"]; lch = 50.0; res = {}

    # F1: beta_c at the uniaxial-compression peak (kp=1 => qh2=1, rho_bar = fc*sqrt(2/3))
    bc_peak = beta_c([-fc, 0.0, 0.0], 1.0, mp)
    bc_expect = ft / (fc * np.sqrt(1.0 + 2.0 * Df * Df))     # ft*qh2*sqrt(2/3)/(fc*sqrt(2/3)*sqrt(1+2Df^2))
    res["F1_bc_peak"] = bc_peak; res["F1_bc_expect"] = bc_expect
    res["F1_ok"] = bool(abs(bc_peak - bc_expect) < 1.0e-12 and 0.0 < bc_peak < 1.0)

    # F2 + F3: monotonic compression with the real beta_c vs the beta_c=1 simplification (fine resolution
    # — the damage onset region is step-sensitive, so use a resolved path)
    path = np.linspace(0, -0.3, 5000)
    dr = drive_uniaxial_compression_damaged(mp, path, Gc, lch, As=As, beta_c_on=True)
    d1 = drive_uniaxial_compression_damaged(mp, path, Gc, lch, As=As, beta_c_on=False)
    res["F2_peak"] = float(-np.min(dr["sig11"])); res["F2_peak_err"] = abs(res["F2_peak"] / fc - 1.0)
    res["F2_softens"] = bool(-dr["sig11"][-1] < 0.1 * fc and dr["wc"][-1] > 0.85)
    res["F2_eff_monotone"] = bool(-dr["sig_eff"][-1] > -dr["sig_eff"][int(np.argmin(dr["sig11"]))])
    res["F2_ok"] = bool(res["F2_peak_err"] < 0.03 and res["F2_softens"] and res["F2_eff_monotone"])
    # F3 non-tautology: the real beta_c is markedly MORE DUCTILE than beta_c=1 — the post-peak STRESS
    # differs by tens of MPa (omega_c eventually saturates to ~1 for both at deep strain, so the wc gap
    # is small there; the STRESS gap is the unambiguous metric). beta_c is a genuine state-dependent
    # factor in (0,1), not a constant.
    res["F3_stress_gap"] = float(np.max(np.abs(dr["sig11"] - d1["sig11"])))
    bcs = [beta_c([dr["sig_eff"][k], 0.0, 0.0], dr["kp"][k], mp) for k in range(len(path)) if dr["wc"][k] > 1.0e-9]
    res["F3_bc_min"] = float(min(bcs)); res["F3_bc_max"] = float(max(bcs))
    res["F3_ok"] = bool(res["F3_stress_gap"] > 0.1 * fc                       # real beta_c clearly more ductile
                        and 0.0 < res["F3_bc_min"] <= res["F3_bc_max"] < 1.0)  # beta_c in (0,1) (~const in uniaxial comp)

    # F4 — full cyclic compression load->unload->reload. beta_c (Eq.50) supplies the correct compressive-
    # damage RATE but does NOT by itself make the damage variable monotone. P2g drives omega_c with the
    # MONOTONE running-max of the effective drive stress, so an elastic UNLOAD (drive drops, histories
    # frozen) keeps omega_c FIXED (no healing) and the nominal stress unloads along the degraded secant.
    # Now a real gate (was a reported diagnostic in #321).
    seg1 = np.linspace(0.0, -2.0e-3, 700)
    seg2 = np.linspace(-2.0e-3, -0.5e-3, 250)        # elastic unload
    seg3 = np.linspace(-0.5e-3, -4.0e-3, 1200)       # reload past the previous max
    cyc = np.concatenate([seg1, seg2[1:], seg3[1:]])
    dc = drive_uniaxial_compression_damaged(mp, cyc, Gc, lch, As=As, beta_c_on=True)
    res["F4_wc_heals_on_unload"] = bool(np.min(np.diff(dc["wc"])) < -1.0e-9)   # False now (monotone omega_c)
    res["F4_wc_peak"] = float(np.max(dc["wc"]))
    res["F4_ok"] = bool(not res["F4_wc_heals_on_unload"] and res["F4_wc_peak"] > 0.1)

    ok = res["F1_ok"] and res["F2_ok"] and res["F3_ok"] and res["F4_ok"]
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gc={Gc} As={As} Df={Df} lch={lch}")
        print(f"  F1 beta_c peak={bc_peak:.5f} (closed form {bc_expect:.5f}) in (0,1) ({res['F1_ok']})")
        print(f"  F2 monotonic backbone: peak={res['F2_peak']:.3f} (fc={fc}) softens={res['F2_softens']} "
              f"eff-monotone={res['F2_eff_monotone']} ({res['F2_ok']})")
        print(f"  F3 non-tautology: max stress gap (real beta_c vs beta_c=1)={res['F3_stress_gap']:.2f} MPa (>{0.1*fc:.1f}); "
              f"beta_c in [{res['F3_bc_min']:.3f},{res['F3_bc_max']:.3f}] ({res['F3_ok']})")
        print(f"  F4 cyclic omega_c heals on unload = {res['F4_wc_heals_on_unload']} (=False; monotone-omega_c, P2g) "
              f"wc_peak={res['F4_wc_peak']:.3f} ({res['F4_ok']})")
        print(f"  => P2f GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P2g — MONOTONE (no-heal) cyclic damage. omega_t/omega_c were re-solved every step against the LIVE
# effective drive stress (sig_t_drive/sig_c_drive). The kappa-histories are already monotone (they
# accumulate only when loading), so the live drive is the ONLY non-monotone input: on an elastic UNLOAD
# the drive drops, F(0)=drive - f*exp(-kd1/eps_f) goes <=0, and the bracketed solve relaxes omega back
# (the material spuriously HEALS — the #321 F4 diagnostic). CDPM2 states omega = omega(kappa_d): a
# function of the MONOTONE history only. The fix tracks the running max of each channel's drive stress
# (sigt_max/sigc_max) and solves omega against THAT, so omega is monotone-nondecreasing and the nominal
# stress unloads along the degraded damage secant (1-omega)*sig_bar. On any MONOTONIC path max == live,
# so the change is byte-identical to the pre-P2g drivers (DT1/DT2/P2e/P2f backbones are unaffected).
# The analytic damaged tangent drops the -sig(x)d(omega) rank-update on an unloading channel (frozen
# drive => d(omega)=0), so the unload tangent is the SPD damage secant D_dam:C_eff (well-conditioned,
# contrast the INDEFINITE Tier-1 loading tangent of gate TD2).
# ===========================================================================
def run_p2g_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    """P2g — monotone (no-heal) cyclic damage. Gates:
      G1 tension load->unload: omega_t monotone (no heal) + secant unload (sig_nom == (1-wt)*sig_bar)
      G2 tension reload BELOW the previous peak retraces the secant (omega_t frozen, no strength regain)
      G3 compression load->unload: omega_c monotone (no heal) + secant unload (the F4 fix)
      G4 omega never decreases across a tension->compression reversal (persistent crack, no cross-heal)
      G5 reduce-to-monotonic: the single-step tensor update == the path driver on a monotonic path
      G6 the UNLOAD damaged tangent is the SPD secant (analytic == numerical; -sig(x)d(omega) vanished)"""
    mp = make_material(E, nu, fc, ft); lch = 50.0; res = {}

    # ---- G1/G2 tension: load deep into softening, elastic-unload to half, reload past the previous peak
    pk = 0.0012
    nL, nU, nR = 600, 300, 600
    tpath = np.concatenate([np.linspace(0.0, pk, nL),
                            np.linspace(pk, 0.5 * pk, nU)[1:],            # elastic unload
                            np.linspace(0.5 * pk, 1.5 * pk, nR)[1:]])     # reload past previous max
    d = drive_damaged_unified(mp, tpath, Gf, Gc, lch, As)
    wt, se, sn = d["wt"], d["sig_eff"], d["sig11"]
    wt_peak = float(wt[nL - 1])
    unl = slice(nL, nL + nU - 1)                                          # the unload window
    res["G1_wt_min_diff"] = float(np.min(np.diff(wt)))                    # >= ~0 over the WHOLE path
    res["G1_wt_peak"] = wt_peak
    res["G1_secant_err"] = float(np.max(np.abs(sn[unl] - (1.0 - wt_peak) * se[unl])) / ft)
    res["G1_ok"] = bool(res["G1_wt_min_diff"] > -1.0e-9 and 0.3 < wt_peak < 0.999
                        and res["G1_secant_err"] < 1.0e-9)
    rel = slice(nL + nU - 1, nL + nU - 1 + 200)                          # reload, safely below the peak
    res["G2_wt_frozen"] = float(np.max(wt[rel]) - np.min(wt[rel]))        # ~0 (omega frozen on reload<peak)
    res["G2_ok"] = bool(res["G2_wt_frozen"] < 1.0e-9)

    # ---- G3 compression: load past the damage onset into softening, elastic-unload (the F4 fix). The
    #      unload stays COMPRESSIVE (a small dEps; a large unload would elastically overshoot the big
    #      residual plastic strain into tension, where the stress routes through the (1-wc) -> (1-wt)
    #      channel and the single-channel secant identity no longer applies).
    pc = -0.06
    cpath = np.concatenate([np.linspace(0.0, pc, 1500), np.linspace(pc, pc + 0.002, 400)[1:]])
    dc = drive_damaged_unified(mp, cpath, Gf, Gc, lch, As)
    wc, sec, snc = dc["wc"], dc["sig_eff"], dc["sig11"]
    wc_peak = float(wc[1500 - 1])
    unlc = slice(1500, 1500 + 399)
    res["G3_wc_min_diff"] = float(np.min(np.diff(wc)))
    res["G3_wc_peak"] = wc_peak
    res["G3_secant_err"] = float(np.max(np.abs(snc[unlc] - (1.0 - wc_peak) * sec[unlc])) / fc)
    res["G3_ok"] = bool(res["G3_wc_min_diff"] > -1.0e-9 and 0.1 < wc_peak < 0.999
                        and res["G3_secant_err"] < 1.0e-9)

    # ---- G4 monotone omega across a tension->compression reversal (the crack persists, no cross-heal)
    rpath = np.concatenate([np.linspace(0.0, 0.0009, 400), np.linspace(0.0009, -0.02, 800)[1:]])
    dr4 = drive_damaged_unified(mp, rpath, Gf, Gc, lch, As)
    res["G4_wt_min_diff"] = float(np.min(np.diff(dr4["wt"])))
    res["G4_wc_min_diff"] = float(np.min(np.diff(dr4["wc"])))
    res["G4_ok"] = bool(res["G4_wt_min_diff"] > -1.0e-9 and res["G4_wc_min_diff"] > -1.0e-9)

    # ---- G5 reduce-to-monotonic: chain the single-step tensor update (the C++ setTrialStrain contract)
    #      along the recorded uniaxial-stress monotonic tension path; nominal axial == the path driver.
    mono = np.linspace(0.0, 0.0009, 500)
    dm = drive_damaged_unified(mp, mono, Gf, Gc, lch, As)
    st = make_damage_state(mp); chain = []
    for k in range(len(mono)):
        eps_t = np.array([mono[k], dm["eps_lat"][k], dm["eps_lat"][k], 0.0, 0.0, 0.0])
        s, st, _ = damaged_step_tensor(st, eps_t - st["eps"], mp, Gf, Gc, lch, As)
        chain.append(s[0])
    res["G5_maxdiff"] = float(np.max(np.abs(np.array(chain) - dm["sig11"])))
    res["G5_ok"] = bool(res["G5_maxdiff"] < 1.0e-9)

    # ---- G6 the UNLOAD tangent is the SPD damage secant. Build a tensile-damaged committed state on a
    #      uniaxial-STRAIN path (no apex fragility), take one elastic-UNLOAD step, and check the analytic
    #      damaged tangent == the numerical central difference AND its symmetric part is SPD (lambda_min>0)
    #      — i.e. the -sig(x)d(omega) rank-update vanished (frozen omega) leaving D_dam:C_eff.
    st = make_damage_state(mp)
    for e in np.linspace(0.0, 0.0009, 400):
        st = damaged_step_tensor(st, np.array([e, 0, 0, 0, 0, 0]) - st["eps"], mp, Gf, Gc, lch, As)[1]
    deps_unl = np.array([-2.0e-5, 0.0, 0.0, 0.0, 0.0, 0.0])
    _, _, info = damaged_step_tensor(st, deps_unl, mp, Gf, Gc, lch, As)
    Ca = damaged_tangent_analytic(st, deps_unl, mp, Gf, Gc, lch, As)
    Cn = damaged_consistent_tangent(st, deps_unl, mp, Gf, Gc, lch, As)
    res["G6_wt_at_unload"] = float(info["wt"])
    res["G6_tangent_relerr"] = float(np.max(np.abs(Ca - Cn)) / (np.max(np.abs(Cn)) + 1.0e-30))
    res["G6_lambda_min"] = float(np.min(np.linalg.eigvalsh(0.5 * (Ca + Ca.T))))
    res["G6_ok"] = bool(info["wt"] > 0.3 and res["G6_tangent_relerr"] < 1.0e-6 and res["G6_lambda_min"] > 0.0)

    ok = (res["G1_ok"] and res["G2_ok"] and res["G3_ok"] and res["G4_ok"]
          and res["G5_ok"] and res["G6_ok"])
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As} lch={lch}")
        print(f"  G1 tension no-heal: wt_peak={wt_peak:.3f} min d(wt)={res['G1_wt_min_diff']:.1e} "
              f"secant_err={res['G1_secant_err']:.1e} ({res['G1_ok']})")
        print(f"  G2 reload<peak omega frozen: span={res['G2_wt_frozen']:.1e} ({res['G2_ok']})")
        print(f"  G3 compression no-heal (F4 fix): wc_peak={wc_peak:.3f} min d(wc)={res['G3_wc_min_diff']:.1e} "
              f"secant_err={res['G3_secant_err']:.1e} ({res['G3_ok']})")
        print(f"  G4 reversal monotone: min d(wt)={res['G4_wt_min_diff']:.1e} min d(wc)={res['G4_wc_min_diff']:.1e} "
              f"({res['G4_ok']})")
        print(f"  G5 single-step == path driver (monotonic): maxdiff={res['G5_maxdiff']:.1e} ({res['G5_ok']})")
        print(f"  G6 unload tangent SPD secant: wt={res['G6_wt_at_unload']:.3f} "
              f"analytic-vs-FD relerr={res['G6_tangent_relerr']:.1e} lambda_min={res['G6_lambda_min']:.3e} "
              f"({res['G6_ok']})")
        print(f"  => P2g GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P2h — compression->tension damage-coupling TEMPER (the `-ctTemper` modes). Per literal CDPM2 (Eq.43,
# kappa_dt-dot = eps_tilde-dot, no (1-alpha_c)) a compression excursion accumulates the TENSILE damage
# history, so a subsequent tension reload comes back PRE-DAMAGED to ~0 (the DT5 diagnostic). The temper
# scales the tensile-history accumulation by a weight w_t (see tensile_damage_weight):
#   'none'   (default) w_t=1            -> literal CDPM2 (byte-identical to the shipped un-tempered material)
#   'alphat'           w_t=1-alpha_c    -> compression (alpha_c~1)=>w_t~0; restores tension after compression
#                                          AND leaves both monotonic backbones exact (alpha_c=0 in tension)
#   'proj'             w_t=||<deps_p>+||/||deps_p|| -> kinematic tensile fraction; restores tension but also
#                                          lightly softens the monotonic tension backbone.
# ===========================================================================
def run_p2h_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    """P2h — the ctTemper compression->tension coupling modes. H0 'none' keeps faithful CDPM2 (tension
    dies after compression, byte-identical monotonic tension to the default); H1 'alphat' restores
    tension-after-compression to ~ft AND keeps the monotonic tension backbone byte-identical; H2 'proj'
    restores tension (looser) with a lightly-softened monotonic tension backbone; H3 both modes keep
    omega monotone (P2g no-heal preserved); H4 analytic == numerical damaged tangent for both modes."""
    lch = 50.0; res = {}
    cpath = np.concatenate([np.linspace(0, -0.01, 800), np.linspace(-0.01, 0.004, 800)[1:]])
    mono_t = np.linspace(0, 0.004, 800)

    def _tac(mode):                                          # tension-after-compression peak + the path
        mp = make_material(E, nu, fc, ft, ct_temper=mode)
        d = drive_damaged_unified(mp, cpath, Gf, Gc, lch, As)
        return float(np.max(d["sig11"][800:])), d

    def _mono(mode):
        mp = make_material(E, nu, fc, ft, ct_temper=mode)
        return drive_damaged_unified(mp, mono_t, Gf, Gc, lch, As)["sig11"]

    p_none, _ = _tac("none")
    p_at, d_at = _tac("alphat")
    p_pr, d_pr = _tac("proj")
    res["H0_none_tac"] = p_none; res["H1_alphat_tac"] = p_at; res["H2_proj_tac"] = p_pr
    res["H0_ok"] = bool(p_none < 0.2 * ft)                   # faithful: tension dies after compression

    mn = _mono("none")
    res["H1_mono_maxdiff"] = float(np.max(np.abs(mn - _mono("alphat"))))
    res["H1_restored"] = bool(p_at > 0.7 * ft)
    res["H1_ok"] = bool(res["H1_restored"] and res["H1_mono_maxdiff"] < 1.0e-7)

    res["H2_mono_reldiff"] = float(np.max(np.abs(_mono("proj") - mn)) / ft)
    res["H2_ok"] = bool(p_pr > 0.5 * ft and res["H2_mono_reldiff"] < 0.2)

    res["H3_alphat_wt_mono"] = float(np.min(np.diff(d_at["wt"])))
    res["H3_proj_wt_mono"] = float(np.min(np.diff(d_pr["wt"])))
    res["H3_ok"] = bool(res["H3_alphat_wt_mono"] > -1.0e-9 and res["H3_proj_wt_mono"] > -1.0e-9)

    def _tan_rel(mode):                                      # analytic vs numerical damaged tangent
        mp = make_material(E, nu, fc, ft, ct_temper=mode)
        st = make_damage_state(mp)
        for e in np.linspace(0, 0.0009, 400):
            st = damaged_step_tensor(st, np.array([e, 0, 0, 0, 0, 0]) - st["eps"], mp, Gf, Gc, lch, As)[1]
        deps = np.array([3.0e-5, 0, 0, 0, 0, 0])
        Ca = damaged_tangent_analytic(st, deps, mp, Gf, Gc, lch, As)
        Cn = damaged_consistent_tangent(st, deps, mp, Gf, Gc, lch, As)
        return float(np.max(np.abs(Ca - Cn)) / (np.max(np.abs(Cn)) + 1.0e-30))
    res["H4_alphat_tan"] = _tan_rel("alphat")
    res["H4_proj_tan"] = _tan_rel("proj")
    res["H4_ok"] = bool(res["H4_alphat_tan"] < 1.0e-5 and res["H4_proj_tan"] < 1.0e-4)

    ok = res["H0_ok"] and res["H1_ok"] and res["H2_ok"] and res["H3_ok"] and res["H4_ok"]
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As} lch={lch}")
        print(f"  H0 none (faithful): tension-after-compression peak={p_none:.4f} (<{0.2*ft:.2f}; dies) ({res['H0_ok']})")
        print(f"  H1 alphat: tac peak={p_at:.4f} (>{0.7*ft:.2f}=restored)  monotonic-tension maxdiff vs none="
              f"{res['H1_mono_maxdiff']:.1e} (byte) ({res['H1_ok']})")
        print(f"  H2 proj: tac peak={p_pr:.4f} (>{0.5*ft:.2f}=restored)  monotonic-tension reldiff vs none="
              f"{res['H2_mono_reldiff']:.3f} (<0.2) ({res['H2_ok']})")
        print(f"  H3 omega monotone (no heal): alphat min d(wt)={res['H3_alphat_wt_mono']:.1e} "
              f"proj min d(wt)={res['H3_proj_wt_mono']:.1e} ({res['H3_ok']})")
        print(f"  H4 analytic==numerical tangent: alphat rel={res['H4_alphat_tan']:.1e} "
              f"proj rel={res['H4_proj_tan']:.1e} ({res['H4_ok']})")
        print(f"  => P2h GATE {'PASS' if ok else 'FAIL'}")
    return res


def run_p2i_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    """P2i — MULTIAXIAL-DAMAGE APPORTIONING. The tensile omega-solve is now driven by E*eps_tilde (the
    CDPM2 equivalent strain, Eq.37) instead of the extreme tensile principal max<sig_bar_i>+; the
    COMPRESSIVE channel keeps the extreme principal (eps_tilde is ft-scaled => E*eps_tilde could never
    reach fc, so it would never onset wc). USER decision 2026-06-21.

    I1  reduce-to-uniaxial (byte): in uniaxial tension E*eps_tilde == sig_bar_t, so the unified driver
        still matches the uniaxial-tension reference driver byte-for-byte (the change is invisible there,
        hence DT1/D1/P2a are unaffected).
    I2  multiaxial escalation: BOTH equibiaxial and triaxial tension drive E*eps_tilde ABOVE the uniaxial
        ft, so they reach the damage onset eps0 at a LOWER per-principal stress than ft — the CDPM2-
        consistent tensile envelope the extreme-principal model can't reproduce (it keeps every principal
        at ft). Ordering is uni < TRI < BI (NOT uni<bi<tri): hydrostatic-triaxial tension is APEX-CAPPED
        by the MW hydrostatic-tension vertex, so it escalates LESS than the deviatoric equibiaxial state.
    I3  no spurious compression->tension damage: the tension gate keeps wt == 0 in pure compression
        (E*eps_tilde > 0 there, but the gate blocks it; wt only affects the stress via the positive-part
        split anyway).
    I4  the analytic damaged tangent, with its tensile drive-gradient now d(E*eps_tilde)/deps = E*det_deps
        (replacing the extreme-principal eigenprojection), == the numerical reference at a genuinely
        biaxial-tension DAMAGED state."""
    mp = make_material(E, nu, fc, ft)
    eps0 = ft / E
    res = {}

    # I1 — reduce to the uniaxial-tension reference driver byte-for-byte.
    path = np.linspace(0.0, 6.0 * eps0, 120)
    uni = drive_uniaxial_tension_damaged(mp, path, Gf, lch=50.0)
    unified = drive_damaged_unified(mp, path, Gf, Gc, lch=50.0, As=As)
    i1 = float(np.max(np.abs(uni["sig11"] - unified["sig11"])))
    res["I1_reduce_uniaxial"] = i1

    # I2 — E*eps_tilde escalates with added tensile principals (extreme-principal would stay == ft).
    Et_uni = E * equiv_strain_general(np.array([ft, 0.0, 0.0]), mp)       # == ft (Eq.38)
    Et_bi = E * equiv_strain_general(np.array([ft, ft, 0.0]), mp)         # > ft
    Et_tri = E * equiv_strain_general(np.array([ft, ft, ft]), mp)         # > equibiaxial
    res["I2_uni"], res["I2_bi"], res["I2_tri"] = Et_uni, Et_bi, Et_tri

    # I3 — pure (uniaxial-strain) compression: the tensile drive gate keeps wt == 0.
    st = make_damage_state(mp)
    wt_comp = 0.0
    for e in np.linspace(0.0, -0.01, 80):
        _, st, info = damaged_step_tensor(st, np.array([e, 0.0, 0.0, 0.0, 0.0, 0.0]) - st["eps"],
                                          mp, Gf, Gc, 50.0, As)
        wt_comp = max(wt_comp, info["wt"])
    res["I3_pure_compression_wt"] = wt_comp

    # I4 — analytic == numerical tangent at an UNEQUAL biaxial-tension DAMAGED state, built like the P2e
    # gates: advance to a damaged committed state, then evaluate the tangent at a SMALL loading PROBE (NOT
    # at the large onset-crossing step). Unequal biaxial (e, 0.6e, 0) keeps the tensile principals DISTINCT
    # (an equal-biaxial (e,e,0) state has an exact eigenvalue degeneracy sig_bar_11==sig_bar_22, where the
    # FD rotates the degenerate eigenvectors while the analytic holds them frozen — a known frozen-
    # eigenvector limitation, not a drive-gradient bug). Here E*eps_tilde > the extreme principal, so this
    # genuinely exercises the new tensile drive-gradient E*det_deps in a multiaxial regime.
    nrm = lambda A: float(np.sqrt(np.sum(A * A)))
    st_b, _, wt_b, _ = _advance_damaged(
        make_damage_state(mp),
        [np.array([e, 0.6 * e, 0.0, 0.0, 0.0, 0.0]) for e in np.linspace(0.0, 6.0e-4, 300)],
        mp, Gf, Gc, 50.0, As)
    probe = np.array([1.0e-6, 0.6e-6, 0.0, 0.0, 0.0, 0.0])
    Ca = damaged_tangent_analytic(st_b, probe, mp, Gf, Gc, 50.0, As)
    Cn = damaged_consistent_tangent(st_b, probe, mp, Gf, Gc, 50.0, As)
    i4 = nrm(Ca - Cn) / nrm(Cn)
    res["I4_tangent_rel"] = i4
    res["I4_wt"] = float(wt_b[-1])

    # I2 — both equibiaxial AND triaxial tension ESCALATE above uniaxial (damage onsets at a lower
    # per-principal stress than ft), the CDPM2-consistent multiaxial behavior the extreme-principal model
    # CANNOT reproduce (it keeps every principal at ft). NB the ordering is uni < TRI < BI, NOT uni<bi<tri:
    # equibiaxial (deviatoric) drives eps_tilde HIGHER than hydrostatic-triaxial tension, which is
    # APEX-CAPPED (the MW hydrostatic-tension vertex limits how far eps_tilde grows under (ft,ft,ft)). So
    # assert each escalates above uniaxial; do NOT assume tri>bi. (Plus the Eq.38 identity uni==ft.)
    ok = (i1 < 1.0e-9
          and abs(Et_uni - ft) < 1.0e-6 * ft
          and Et_bi > Et_uni * (1.0 + 1.0e-4) and Et_tri > Et_uni * (1.0 + 1.0e-4)
          and wt_comp < 1.0e-9
          and i4 < 1.0e-6)
    res["PASS"] = ok
    if verbose:
        print(f"  I1 reduce-to-uniaxial       max|dsig11| = {i1:.2e}   (< 1e-9)")
        print(f"  I2 E*eps_tilde escalation   uni={Et_uni:.4f} (=ft) bi={Et_bi:.4f} tri={Et_tri:.4f}")
        print(f"  I3 pure-compression wt      = {wt_comp:.2e}   (< 1e-9, no spurious tension damage)")
        print(f"  I4 analytic==numerical tan  rel = {i4:.2e}   (< 1e-6)")
        print(f"  => P2i GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P3 — Tier-2 IMPL-EX robustness (ADR 4.4). The Tier-1 implicit return map + dual-damage solve is
# accurate, but its consistent tangent is NON-SYMMETRIC and INDEFINITE on the softening branch
# (gate TD2: C[0,0]<0, lambda_min(sym)<0) => the global Newton stalls past the limit point (the
# snap-back stalls we hit on the single-element softening tests). IMPL-EX (Oliver/Huespe) reports an
# EXPLICIT stress assembled from EXTRAPOLATED internal variables (the plastic-strain increment AND
# the dual damage frozen at the committed-history rate) so the algorithmic tangent is a CONSTANT
# degraded-elastic secant D_dam(w~):C0 — symmetric-part SPD across the snap-back — while the
# COMMITTED internal variables are the exact IMPLICIT ones (so the explicit error -> 0 as dt -> 0 and
# the committed trajectory is byte-identical to Tier-1). ADR 4.4 [BLOCKING]: freeze the PLASTIC state
# too, not just damage — freezing only omega leaves Ceff non-associated/non-symmetric and the SPD
# promise is false. Here freezing the plastic-strain increment collapses the effective tangent to the
# elastic C0; freezing omega drops the -sig(x)dw rank-update; the product is the SPD secant.
#
# Extrapolation (ASDConcrete3D pattern): x~ = x_n + (dt/dt_n)*(x_n - x_{n-1}). We carry the committed
# IMPLICIT increment of each frozen variable (dwt, dwc, depl) + the committed dt, so x~ = x_n +
# (dt/dt_n)*dx_n. First step / no history (dt_n=0) => factor 0 => pure elastic-predictor explicit.
# ===========================================================================
# IMPL-EX extrapolation time-ratio cap r = dt/dt_n (adversarial-review ALG-2/NUM-2/NUM-3): an adaptive
# step-GROWTH (a small step then a much larger one, exactly what a step-cutting solver produces) makes
# r large, over-extrapolates the bounded damage past [0,1) (reported stress collapses, tangent goes
# near-singular/indefinite) and injects an unbounded spurious plastic strain via depl_x = r*depl.
# ASDConcrete3D likewise guards the IMPL-EX time ratio (there via error-control step reduction); we use
# a hard cap. r<0 (a backward/negative dt) is floored to 0 — extrapolation only advances forward.
_IMPLEX_RMAX = 2.0


def make_implex_state(mp):
    """Committed IMPL-EX state = the Tier-1 damage state + the frozen-variable increments and dt."""
    s = make_damage_state(mp)
    s.update(wt=0.0, wc=0.0, dwt=0.0, dwc=0.0, depl=np.zeros(6), dt_n=0.0)
    return s


def _implex_secant_tangent(sig_bar_x, wt_x, wc_x, mp):
    """The Tier-2 secant tangent: D_dam(omega FROZEN) : C_elastic. With the plastic flow frozen the
    effective tangent collapses to the elastic C0, and freezing omega drops the -sig(x)dw rank-update.
    SPD CAVEAT (adversarial-review NUM-1): D_dam is symmetric, but D_dam @ C0 does NOT commute, so the
    secant is symmetric-part SPD only in SINGLE-SIGN principal regimes (all-tensile or all-compressive
    sig_bar_x, where D_dam is a single positive scaling). On a MIXED-SIGN, high-omega direction-contrast
    state (a tensile-damaged principal beside an undamaged compressive one, omega_t > ~0.97) the two
    branch slopes (1-wt) != (1-wc) make the symmetric part lose definiteness — the intrinsic dual-damage
    IMPL-EX limitation (gate PI5 pins it). It is still far better-conditioned than the Tier-1 tangent and
    the COMMITTED physics is exact (PI2). This IS the consistent tangent of the reported explicit stress
    (gate PI5 FD-verifies), so the global Newton consumes exactly d(sig_rep)/d(deps)."""
    lam, V = np.linalg.eigh(voigt_to_mat(sig_bar_x))
    yv = [(1.0 - wt_x) * max(lam[a], 0.0) + (1.0 - wc_x) * min(lam[a], 0.0) for a in range(3)]
    ypv = [(1.0 - wt_x) if lam[a] > 0.0 else (1.0 - wc_x) for a in range(3)]
    return isotropic_tangent(lam, V, yv, ypv) @ elastic_C(mp)


def damaged_step_implex(state, deps6, dt, mp, Gf, Gc, lch, As=2.0):
    """ONE Tier-2 IMPL-EX update (pure). Returns (sig_reported[6], C_spd[6,6], new_state, diag).
      * EXPLICIT (reported): effective stress with the plastic-strain increment FROZEN (extrapolated)
        + dual damage FROZEN (extrapolated) => sig_bar_x is LINEAR in deps => tangent D_dam(w~):C0.
      * IMPLICIT (committed): the exact `damaged_step_tensor`; commit its internal variables + the
        per-variable increments (dwt, dwc, depl) and dt for the NEXT step's extrapolation.
    The extrapolation time-ratio is CLAMPED to [0, _IMPLEX_RMAX] (see the module note) so adaptive /
    step-cutting stepping cannot over-extrapolate the bounded damage or inject unbounded plastic strain."""
    deps6 = np.asarray(deps6, float)
    # extrapolation factor r = dt/dt_n, clamped to [0, R_MAX]; r<0 (negative dt) -> 0 (forward-only)
    r = min(dt / state["dt_n"], _IMPLEX_RMAX) if (state["dt_n"] > 0.0 and dt > 0.0) else 0.0
    wt_x = min(max(state["wt"] + r * state["dwt"], 0.0), 1.0 - 1.0e-12)
    wc_x = min(max(state["wc"] + r * state["dwc"], 0.0), 1.0 - 1.0e-12)
    depl_x = r * np.asarray(state["depl"], float)                     # frozen plastic-strain increment
    # explicit effective stress = sig_bar_n + C0:(deps - depl_x)  (LINEAR in deps => elastic tangent)
    sig_bar_x = elastic_pred_tensor(state["sig_bar"], deps6 - depl_x, mp)
    lam_x, V_x = np.linalg.eigh(voigt_to_mat(sig_bar_x))
    sig_rep = mat_to_voigt(V_x @ np.diag(apply_damage_principal(lam_x, wt_x, wc_x)) @ V_x.T)
    C_spd = _implex_secant_tangent(sig_bar_x, wt_x, wc_x, mp)
    # implicit solve -> the COMMITTED internal variables (accuracy + the next extrapolation source)
    sig_impl, st_impl, info = damaged_step_tensor(state, deps6, mp, Gf, Gc, lch, As)
    depl_new = _plastic_strain6(st_impl["sig_bar"], st_impl["eps"], mp) \
        - _plastic_strain6(state["sig_bar"], state["eps"], mp)
    new_state = dict(st_impl)
    new_state.update(wt=info["wt"], wc=info["wc"],
                     dwt=info["wt"] - state["wt"], dwc=info["wc"] - state["wc"],
                     depl=depl_new, dt_n=dt)
    diag = dict(wt_x=wt_x, wc_x=wc_x, wt_impl=info["wt"], wc_impl=info["wc"],
                implex_err=max(abs(wt_x - info["wt"]), abs(wc_x - info["wc"])),
                sig_impl=sig_impl, min_eig_sym=float(np.min(np.linalg.eigvalsh(0.5 * (C_spd + C_spd.T)))))
    return sig_rep, C_spd, new_state, diag


def drive_damaged_implex(mp, eps_path6, Gf, Gc, lch, As=2.0, dt=1.0):
    """Chain `damaged_step_implex` along a list of TOTAL strain tensors (uniform dt). Returns dict of
    per-step arrays: reported sigma, implicit sigma, wt/wc (implicit), implex_err, min_eig_sym; and
    the final committed state (for a Tier-1 cross-check)."""
    st = make_implex_state(mp)
    sig_rep, sig_impl, wt, wc, err, mineig = [], [], [], [], [], []
    for eps_t in eps_path6:
        deps = np.asarray(eps_t, float) - st["eps"]
        s, _C, st, d = damaged_step_implex(st, deps, dt, mp, Gf, Gc, lch, As)
        sig_rep.append(s); sig_impl.append(d["sig_impl"])
        wt.append(d["wt_impl"]); wc.append(d["wc_impl"])
        err.append(d["implex_err"]); mineig.append(d["min_eig_sym"])
    return dict(sig_rep=np.array(sig_rep), sig_impl=np.array(sig_impl),
                wt=np.array(wt), wc=np.array(wc), err=np.array(err),
                min_eig_sym=np.array(mineig), state=st)


def _advance_implex(mp, eps_path6, Gf, Gc, lch, As=2.0, dt=1.0):
    """Chain IMPL-EX steps and return the final COMMITTED implex state (for the FD/robustness gates)."""
    st = make_implex_state(mp)
    for eps_t in eps_path6:
        _, _, st, _ = damaged_step_implex(st, np.asarray(eps_t, float) - st["eps"], dt, mp, Gf, Gc, lch, As)
    return st


def run_p3_implex_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    """Tier-2 IMPL-EX falsification battery (ADR 4.4). Oracle-only (numpy); the C++ kernel port is a
    follow-up build PR. PI1: the symmetrized Tier-2 secant stays POSITIVE-DEFINITE across a SINGLE-SIGN
    softening snap-back where the Tier-1 tangent is INDEFINITE (gate TD2). PI2 committed==Tier-1; PI3
    smooth-region O(dt) order across >=3 levels; PI4 error monitor; PI5 secant==d(sig_rep)/d(deps) +
    the PINNED dual-damage SPD limitation (mixed-sign high-omega is NOT SPD — conditional, see NUM-1);
    PI6 robustness under non-uniform dt (the extrapolation-ratio clamp). Adversarial-review-hardened."""
    mp = make_material(E, nu, fc, ft)
    lch = 100.0
    nrm = lambda A: float(np.sqrt(np.sum(np.asarray(A) ** 2)))
    res = {}

    # ---- a uniaxial-STRAIN tension path that drives WELL past the nominal peak into deep softening
    # (same construction as gate TD2, where the Tier-1 tangent is degraded + indefinite). ----
    eps_end = 6.0e-4
    nN = 400
    path = [np.array([e, 0, 0, 0, 0, 0]) for e in np.linspace(eps_end / nN, eps_end, nN)]
    run = drive_damaged_implex(mp, path, Gf, Gc, lch, As)
    ipk = int(np.argmax(run["sig_rep"][:, 0]))
    res["softening_reached"] = bool(run["wt"][-1] > 0.5 and run["sig_rep"][-1, 0] < run["sig_rep"][ipk, 0])

    # PI1 [HEADLINE]: SPD across the snap-back. Every post-onset IMPL-EX tangent has lambda_min(sym)>0,
    # while the Tier-1 tangent at the SAME deep-softening committed state is INDEFINITE (lambda_min<0).
    post = run["min_eig_sym"][run["wt"] > 1.0e-6]
    res["PI1_implex_min_eig"] = float(np.min(post)) if post.size else 0.0
    res["PI1_implex_SPD"] = bool(post.size > 0 and np.all(post > 0.0))
    # the Tier-1 reference at a deep-softening state (rebuild the committed state ~2/3 down the tail)
    isoft = (ipk + nN) // 2 + (nN - ipk) // 6
    s_t1, _, wt1, _ = _advance_damaged(make_damage_state(mp), path[:isoft], mp, Gf, Gc, lch, As)
    C_t1 = damaged_consistent_tangent(s_t1, np.array([2.0e-7, 0, 0, 0, 0, 0]), mp, Gf, Gc, lch, As)
    res["PI1_tier1_min_eig"] = float(np.min(np.linalg.eigvalsh(0.5 * (C_t1 + C_t1.T))))
    res["PI1_tier1_indefinite"] = bool(res["PI1_tier1_min_eig"] < 0.0 and float(wt1[-1]) > 0.5)

    # PI2 consistency: the COMMITTED internal trajectory is the IMPLICIT one => byte-identical to a
    # pure Tier-1 run (IMPL-EX must not corrupt the implicit state; ADR 4.4 finite-strain contract).
    s_ref, sig_ref, wt_ref, wc_ref = _advance_damaged(make_damage_state(mp), path, mp, Gf, Gc, lch, As)
    res["PI2_commit_sig_bar_err"] = nrm(run["state"]["sig_bar"] - s_ref["sig_bar"]) / (nrm(s_ref["sig_bar"]) + 1e-30)
    res["PI2_commit_wt_err"] = abs(run["wt"][-1] - wt_ref[-1])
    res["PI2_commit_matches_tier1"] = bool(res["PI2_commit_sig_bar_err"] < 1e-12 and res["PI2_commit_wt_err"] < 1e-12)

    # PI3 explicit error -> 0 under step refinement. (adversarial-review ALG-1/GAT-2) The GLOBAL-max gap
    # is dominated by an IRREDUCIBLE one-step lag at the damage-ONSET C0 kink (dwt jumps 0->finite, an
    # O(1)-in-rate event that does NOT refine away and is non-monotone) — gating it at a hardcoded N was
    # brittle. Instead measure the gap at a FIXED strain in the SMOOTH softening tail and assert a clean
    # convergence ORDER across >=3 refinement levels.
    eps_eval = 3.5e-4                                                 # smooth softening tail (onset ~ ft/E ~ 1e-4)
    def _gap_at_strain(nsteps):
        grid = np.linspace(eps_end / nsteps, eps_end, nsteps)
        p = [np.array([e, 0, 0, 0, 0, 0]) for e in grid]
        rr = drive_damaged_implex(mp, p, Gf, Gc, lch, As, dt=eps_end / nsteps)
        gap = np.abs(rr["sig_rep"][:, 0] - rr["sig_impl"][:, 0])
        return float(np.interp(eps_eval, grid, gap)), float(np.max(gap))
    levels = [200, 400, 800]
    gaps, gmax = zip(*[_gap_at_strain(n) for n in levels])
    orders = [np.log(gaps[i] / gaps[i + 1]) / np.log(levels[i + 1] / levels[i]) for i in range(len(levels) - 1)]
    res["PI3_gaps_smooth"] = list(gaps)
    res["PI3_orders"] = [float(o) for o in orders]
    res["PI3_min_order"] = float(np.min(orders))
    res["PI3_smooth_monotone"] = bool(all(gaps[i] > gaps[i + 1] for i in range(len(gaps) - 1)))
    res["PI3_onset_lag_floors"] = bool(min(gmax) < max(gmax) * 1.5 and min(gmax) > 0.0)  # global-max does NOT refine (documented)
    res["PI3_converges"] = bool(res["PI3_smooth_monotone"] and res["PI3_min_order"] > 0.8 and max(gaps) < ft)

    # PI4 error monitor is meaningful: > 0 while damage evolves, ~0 in the pre-onset elastic regime.
    res["PI4_err_softening"] = float(np.max(run["err"]))
    res["PI4_err_preonset"] = float(np.max(run["err"][run["wt"] <= 0.0])) if np.any(run["wt"] <= 0.0) else 0.0
    res["PI4_monitor_ok"] = bool(res["PI4_err_softening"] > 1e-4 and res["PI4_err_preonset"] < 1e-12)

    # PI5 the secant IS the consistent tangent of the reported explicit stress (the IMPL-EX quadratic-
    # convergence contract) + the HONEST SPD scope (adversarial-review NUM-1). (a) FD-verify C_spd ==
    # d(sig_rep)/d(deps) at a deep-softening committed state; (b) PIN the dual-damage SPD limitation:
    # the secant is SPD in SINGLE-SIGN regimes but the symmetric part goes INDEFINITE on a MIXED-SIGN
    # high-omega state (tension crack carrying lateral compression), so PI1's SPD is conditional.
    st_soft = _advance_implex(mp, path[:isoft], Gf, Gc, lch, As)
    deps_fd = np.array([2.0e-7, -0.5e-7, 0.3e-7, 1.0e-7, 0.0, 0.0])
    s0, C0c, _, _ = damaged_step_implex(st_soft, deps_fd, 1.0, mp, Gf, Gc, lch, As)
    Cfd = np.zeros((6, 6)); h = 1.0e-9
    for j in range(6):
        dp = np.array(deps_fd, float); dp[j] += h
        dm = np.array(deps_fd, float); dm[j] -= h
        sp, _, _, _ = damaged_step_implex(st_soft, dp, 1.0, mp, Gf, Gc, lch, As)
        sm, _, _, _ = damaged_step_implex(st_soft, dm, 1.0, mp, Gf, Gc, lch, As)
        Cfd[:, j] = (sp - sm) / (2.0 * h)
    res["PI5_fd_consistency_rel"] = nrm(C0c - Cfd) / (nrm(Cfd) + 1e-30)
    # single-sign vs mixed-sign secant at high omega_t (the NUM-1 boundary, directly probed)
    single = _implex_secant_tangent(np.array([2.0, 1.0, 0.5, 0, 0, 0]), 0.99, 0.0, mp)   # all-tensile
    mixed = _implex_secant_tangent(np.array([1.0, -2.0, -2.0, 0, 0, 0]), 0.99, 0.0, mp)   # tension + lateral compression
    res["PI5_singlesign_min_eig"] = float(np.min(np.linalg.eigvalsh(0.5 * (single + single.T))))
    res["PI5_mixedsign_min_eig"] = float(np.min(np.linalg.eigvalsh(0.5 * (mixed + mixed.T))))
    res["PI5_consistent_and_scope_pinned"] = bool(
        res["PI5_fd_consistency_rel"] < 1e-5            # C_spd == d(sig_rep)/d(deps)
        and res["PI5_singlesign_min_eig"] > 0.0         # SPD in single-sign regimes (PI1 generalizes)
        and res["PI5_mixedsign_min_eig"] < 0.0)         # NOT SPD on mixed-sign high-omega (the pinned limitation)

    # PI6 robustness under NON-UNIFORM dt (the r-clamp; adversarial-review ALG-2/NUM-2/NUM-3). From a
    # softening committed state: (a) a large dt-GROWTH step keeps the reported stress finite + sane and
    # SPD (the clamp stops omega over-extrapolating past [0,1) / a plastic-strain blow-up); (b) a
    # NEGATIVE dt floors r to 0 => no backward damage (wt_x == committed wt).
    base = _peak_state = _advance_implex(mp, path[:isoft], Gf, Gc, lch, As)
    dt0 = 1.0
    deps_next = np.array([eps_end / nN, 0, 0, 0, 0, 0])
    sj, Cj, _, dj = damaged_step_implex(base, deps_next, 100.0 * dt0, mp, Gf, Gc, lch, As)  # r would be 100 unclamped
    res["PI6_jump_sig_finite"] = bool(np.all(np.isfinite(sj)) and nrm(sj) < 10.0 * ft)
    res["PI6_jump_wt_x_bounded"] = bool(dj["wt_x"] < 1.0 - 1e-9)                  # NOT saturated to 1
    res["PI6_jump_min_eig"] = float(np.min(np.linalg.eigvalsh(0.5 * (Cj + Cj.T))))
    res["PI6_jump_SPD"] = bool(res["PI6_jump_min_eig"] > 0.0)                     # single-sign path stays SPD
    _, _, _, dn = damaged_step_implex(base, deps_next, -dt0, mp, Gf, Gc, lch, As)  # negative dt
    res["PI6_negdt_no_backward"] = bool(abs(dn["wt_x"] - base["wt"]) < 1e-14)
    res["PI6_robust"] = bool(res["PI6_jump_sig_finite"] and res["PI6_jump_wt_x_bounded"]
                             and res["PI6_jump_SPD"] and res["PI6_negdt_no_backward"])

    ok = (res["softening_reached"] and res["PI1_implex_SPD"] and res["PI1_tier1_indefinite"]
          and res["PI2_commit_matches_tier1"] and res["PI3_converges"] and res["PI4_monitor_ok"]
          and res["PI5_consistent_and_scope_pinned"] and res["PI6_robust"])
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As} lch={lch}")
        print(f"  softening reached (wt_final={run['wt'][-1]:.3f}): {res['softening_reached']}")
        print(f"  PI1 [HEADLINE, single-sign] IMPL-EX SPD across snap-back: lambda_min(sym)={res['PI1_implex_min_eig']:.3e}>0 "
              f"({res['PI1_implex_SPD']}); Tier-1 at same state INDEFINITE lambda_min={res['PI1_tier1_min_eig']:.3e}<0 "
              f"({res['PI1_tier1_indefinite']})")
        print(f"  PI2 committed trajectory == Tier-1: sig_bar rel={res['PI2_commit_sig_bar_err']:.2e} "
              f"wt err={res['PI2_commit_wt_err']:.2e} ({res['PI2_commit_matches_tier1']})")
        print(f"  PI3 smooth-region order (>=3 levels): orders={[f'{o:.2f}' for o in res['PI3_orders']]} "
              f"min={res['PI3_min_order']:.2f} monotone={res['PI3_smooth_monotone']} "
              f"(onset global-max lag floors={res['PI3_onset_lag_floors']}) ({res['PI3_converges']})")
        print(f"  PI4 error monitor: softening max={res['PI4_err_softening']:.3e} pre-onset={res['PI4_err_preonset']:.2e} "
              f"({res['PI4_monitor_ok']})")
        print(f"  PI5 secant==d(sig_rep)/d(deps) rel={res['PI5_fd_consistency_rel']:.2e}; "
              f"SPD scope: single-sign min_eig={res['PI5_singlesign_min_eig']:.2e}>0, "
              f"MIXED-SIGN min_eig={res['PI5_mixedsign_min_eig']:.2e}<0 (pinned limitation) "
              f"({res['PI5_consistent_and_scope_pinned']})")
        print(f"  PI6 non-uniform dt (r-clamp): jump sig finite={res['PI6_jump_sig_finite']} "
              f"wt_x bounded={res['PI6_jump_wt_x_bounded']} SPD={res['PI6_jump_SPD']} "
              f"neg-dt no-backward={res['PI6_negdt_no_backward']} ({res['PI6_robust']})")
        print(f"  => P3 IMPL-EX GATE {'PASS' if ok else 'FAIL'}")
    return res


# ===========================================================================
# P3 Duvaut-Lions viscoplastic regularization (-eta, ADR 4.4)
# ---------------------------------------------------------------------------
# The Duvaut-Lions model relaxes the inviscid (rate-independent) plastic return toward the elastic
# trial over a characteristic time eta. Simo & Hughes (Computational Inelasticity, Box) backward-Euler
# closed form:  sigma_{n+1} = (sigma_trial + (dt/eta) sigma_bar_inviscid) / (1 + dt/eta)
#            =  (1-beta) sigma_trial + beta sigma_bar_inviscid ,   beta = dt/(eta+dt) in [0,1].
#   eta -> 0   (beta -> 1): the inviscid return, recovered BYTE-for-byte (the byte gate).
#   eta -> inf (beta -> 0): the frozen elastic trial (infinitely fast loading, no time to relax).
# We apply it at the PLASTIC (effective-stress) level (ADR: "Duvaut-Lions at the PLASTIC level"), so
# the EFFECTIVE stress and its hardening variable kappa_p relax; the damage then follows from the
# relaxed effective stress. The consistent tangent blends C_eff <- (1-beta)C_elastic + beta C_inviscid.
#
# Closed-form 1-D overstress oracle (the ADR's named validation): a scalar perfectly-plastic bar
# (yield sigY) under constant strain rate. The DISCRETE backward-Euler fixed point has overstress
#   sigma_steady - sigY = E * eps_rate * eta   EXACTLY (independent of dt) — see derivation in PV3.
# The CONTINUOUS solution is sigma(t) = sigY + E eps_rate eta (1 - exp(-(t-t_y)/eta)) for t > t_y.
# ===========================================================================
def duvaut_lions_1d_discrete(E, sigY, eps_rate, eta, dt, n):
    """Independent scalar Duvaut-Lions backward-Euler integrator, 1-D perfect plasticity from rest.
    Inviscid projection = clip(sigma_trial, -sigY, sigY); blend with beta=dt/(eta+dt). Returns (t, sigma)."""
    beta = dt / (eta + dt) if (eta > 0.0 and dt > 0.0) else 1.0
    sig = 0.0
    t = np.empty(n); s = np.empty(n)
    for k in range(n):
        sig_tr = sig + E * eps_rate * dt
        sig_bar = min(max(sig_tr, -sigY), sigY)          # 1-D inviscid return (perfect plasticity)
        sig = (1.0 - beta) * sig_tr + beta * sig_bar
        t[k] = (k + 1) * dt; s[k] = sig
    return t, s


def duvaut_lions_1d_analytic(E, sigY, eps_rate, eta, t):
    """Continuous Duvaut-Lions solution, 1-D perfect plasticity, constant strain rate from rest:
    elastic sigma=E eps_rate t until yield at t_y=sigY/(E eps_rate); then the overstress relaxes
    sigma(t) = sigY + E eps_rate eta (1 - exp(-(t-t_y)/eta))."""
    t = np.asarray(t, float)
    t_y = sigY / (E * eps_rate)
    return np.where(t <= t_y, E * eps_rate * t,
                    sigY + E * eps_rate * eta * (1.0 - np.exp(-(t - t_y) / eta)))


def _confined_compression_state(mp, Gf, Gc, lch, As, e11_end=-3.0e-3, n=120):
    """A CLEAN plastic, PRE-onset (omega=0) committed state + the next strain increment, for the eta
    tangent/overstress gates. Built from the uniaxial-STRESS mixed-control driver (which holds the
    effective lateral stress at ~0 and stays on the compressive meridian) — its diagonal total-strain
    path [eps11, eps_lat, eps_lat] is replayed through `damaged_step_tensor`, reproducing a genuinely
    plastic hardening state WITHOUT the oracle's uniaxial-STRAIN deep-compression apex chaos (kp<0,
    spuriously elastic). Returns (state, deps_probe, kp_probe, plastic_probe)."""
    d = drive_damaged_unified(mp, np.linspace(0.0, e11_end, n), Gf, Gc, lch, As)
    path = [np.array([d["eps11"][k], d["eps_lat"][k], d["eps_lat"][k], 0.0, 0.0, 0.0]) for k in range(n)]
    idx = next(k for k in range(1, n) if 0.5 < d["kp"][k] < 0.98 and d["wc"][k] < 1.0e-12)  # hardening, pre-damage
    st, _, _, _ = _advance_damaged(make_damage_state(mp), path[:idx], mp, Gf, Gc, lch, As)
    deps = path[idx] - path[idx - 1]
    _, kp_probe, plastic_probe, _ = return_map_tensor(st["sig_bar"], deps, mp, st["kp"])
    return st, deps, kp_probe, plastic_probe


def run_p3_eta_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, Gc=5.0, As=2.0, verbose=True):
    """Duvaut-Lions viscoplastic (-eta) falsification battery (ADR 4.4). Oracle-only (numpy); the C++
    kernel port is a follow-up build PR. PV1 eta=0 byte-identical to the inviscid Tier-1 path; PV2 the
    scalar backward-Euler integrator converges to the CONTINUOUS closed-form overstress as dt->0; PV3
    the EXACT discrete steady overstress = E eps_rate eta (the closed-form 1-D oracle, dt-independent);
    PV4 the tensor kernel's relaxed effective stress IS the (1-beta)trial + beta inviscid blend;
    PV5 the viscous damaged tangent matches its numerical FD (PV5a) and reduces to the blended effective
    tangent pre-onset (PV5b); PV6 the overstress grows monotonically with eta (regularization signature)."""
    mp = make_material(E, nu, fc, ft)
    mp_eta = lambda eta: make_material(E, nu, fc, ft, eta=eta)
    nrm = lambda A: float(np.sqrt(np.sum(np.asarray(A) ** 2)))
    lch = 50.0
    res = {}

    # --- PV1: eta=0 is BYTE-identical to the inviscid Tier-1 damaged update (over a damaging path) ---
    st0 = make_damage_state(mp); st_e = make_damage_state(mp_eta(0.0))
    max_byte = 0.0
    for k in range(200):
        deps = np.array([4.0e-5, 0, 0, 0, 0, 0])          # uniaxial-strain tension into softening
        s_inv, st0, _ = damaged_step_tensor(st0, deps, mp, Gf, Gc, lch, As)
        s_eta, st_e, _ = damaged_step_tensor(st_e, deps, mp_eta(0.0), Gf, Gc, lch, As, dt=1.0)
        max_byte = max(max_byte, nrm(s_eta - s_inv))
    res["PV1_eta0_byte"] = max_byte
    res["PV1_byte_exact"] = bool(max_byte == 0.0)

    # --- PV2 + PV3: the closed-form 1-D overstress oracle ---
    sigY, eps_rate, eta = 30.0, 1.0e-3, 0.5
    # PV3 (EXACT): the discrete backward-Euler fixed point. In the plastic regime sig_bar=sigY:
    #   sig* = (1-beta) sig* + (1-beta) E eps_rate dt + beta sigY  =>  sig* - sigY = ((1-beta)/beta) E eps_rate dt
    #   (1-beta)/beta = eta/dt  =>  sig* - sigY = (eta/dt) E eps_rate dt = E eps_rate eta  (dt CANCELS).
    over_exact = E * eps_rate * eta
    res["PV3_overstress_target"] = over_exact
    pv3_err = []
    for dt in (1.0, 0.25, 0.05):                           # many steps to reach the fixed point
        n = int(60.0 / (eps_rate * E) / dt) + 4000         # well past yield, into steady state
        _, s = duvaut_lions_1d_discrete(E, sigY, eps_rate, eta, dt, n)
        pv3_err.append(abs((s[-1] - sigY) - over_exact) / over_exact)
    res["PV3_steady_rel_err"] = max(pv3_err)
    res["PV3_exact"] = bool(res["PV3_steady_rel_err"] < 1.0e-10)
    # PV2 (transient): the discrete blend IS backward-Euler of the Duvaut-Lions ODE (1/(1+dt/eta) per
    # step vs the exact exp(-dt/eta)) so it carries an O(dt) GLOBAL error that must converge at order 1.
    # Evaluate MID-transient (t_y + 2*eta, where exp(-2)~0.14 of the overstress is still relaxing — NOT
    # at steady state, where both would be exact and the order is undefined). t_y = sigY/(E*eps_rate) =
    # 1.0; pick T = t_y + 2*eta = 2.0 (a clean multiple of every dt so yield lands on a step boundary,
    # isolating the ODE-integration error from the yield-crossing error).
    t_y = sigY / (E * eps_rate)
    T = t_y + 2.0 * eta
    pv2_err = []
    for dt in (0.5, 0.25, 0.125):
        n = int(round(T / dt))
        t, s = duvaut_lions_1d_discrete(E, sigY, eps_rate, eta, dt, n)
        pv2_err.append(abs(s[-1] - float(duvaut_lions_1d_analytic(E, sigY, eps_rate, eta, t[-1]))))
    res["PV2_transient_errs"] = pv2_err
    res["PV2_order"] = float(np.log(pv2_err[0] / pv2_err[-1]) / np.log(4.0))   # dt halved twice => /4
    res["PV2_converges"] = bool(pv2_err[-1] < pv2_err[0] and res["PV2_order"] > 0.8)

    # --- PV4: the TENSOR kernel relaxed effective stress IS the independent (1-beta)trial+beta inviscid blend ---
    eta4, dt4 = 0.3, 1.0
    beta4 = dt4 / (eta4 + dt4)
    m4 = mp_eta(eta4)
    pv4 = []
    # tension and compression committed states, plus a fresh state, each probed with a plastic increment
    stT, _, _, _ = _advance_damaged(make_damage_state(m4), [np.array([2.0e-4 * i, 0, 0, 0, 0, 0]) for i in range(1, 30)], m4, Gf, Gc, lch, As)
    stC, _, _, _ = _advance_damaged(make_damage_state(m4), [np.array([-3.0e-4 * i, 0, 0, 0, 0, 0]) for i in range(1, 30)], m4, Gf, Gc, lch, As)
    for st, deps in ((make_damage_state(m4), np.array([3.0e-4, 0, 0, 0, 0, 0])),
                     (stT, np.array([2.0e-4, 0, 0, 0, 0, 0])),
                     (stC, np.array([-3.0e-4, 0, 0, 0, 0, 0])),
                     (stC, np.array([1.0e-4, 1.0e-4, 0, 0.5e-4, 0, 0]))):   # mixed/shear increment
        sig_bar_inv, _, _, _ = return_map_tensor(st["sig_bar"], deps, m4, st["kp"])
        sig_tr = elastic_pred_tensor(st["sig_bar"], deps, m4)
        blend = (1.0 - beta4) * sig_tr + beta4 * sig_bar_inv
        _, new_st, _ = damaged_step_tensor(st, deps, m4, Gf, Gc, lch, As, dt=dt4)
        pv4.append(nrm(new_st["sig_bar"] - blend) / max(nrm(blend), 1.0))
    res["PV4_blend_rel"] = max(pv4)
    res["PV4_is_blend"] = bool(res["PV4_blend_rel"] < 1.0e-12)

    # --- PV5: the viscous damaged consistent tangent ---
    # PV5a: analytic == numerical FD at a smooth viscous-damaged state (the C++ deliverable's check).
    eta5, dt5 = 0.4, 1.0
    m5 = mp_eta(eta5)
    st5, _, _, _ = _advance_damaged(make_damage_state(m5), [np.array([5.0e-5 * i, 0, 0, 0, 0, 0]) for i in range(1, 60)], m5, Gf, Gc, lch, As)
    deps5 = np.array([5.0e-5, 0, 0, 0, 0, 0])
    Ca = damaged_tangent_analytic(st5, deps5, m5, Gf, Gc, lch, As, dt=dt5)
    Cn = damaged_consistent_tangent(st5, deps5, m5, Gf, Gc, lch, As, dt=dt5)
    res["PV5a_fd_rel"] = nrm(Ca - Cn) / nrm(Cn)
    res["PV5a_ok"] = bool(res["PV5a_fd_rel"] < 1.0e-5)
    # PV5b: at a genuinely PLASTIC but PRE-onset step (omega=0) the viscous damaged tangent reduces to
    # the blended effective tangent (1-beta)C_el + beta C_eff_inviscid (D_dam=I and the -sig(x)domega
    # rank-updates vanish at omega=0). The state is a CLEAN confined-compression hardening state (NOT a
    # tautological elastic step — PV5b_plastic asserts the probe yields).
    beta5 = dt5 / (eta5 + dt5)
    stc, deps5b, _, plastic5b = _confined_compression_state(m5, Gf, Gc, lch, As)
    _, _, infob = damaged_step_tensor(stc, deps5b, m5, Gf, Gc, lch, As, dt=dt5)
    Cb = damaged_tangent_analytic(stc, deps5b, m5, Gf, Gc, lch, As, dt=dt5)
    Ceff_inv = consistent_tangent(stc["sig_bar"], deps5b, m5, stc["kp"], hardening=True)
    Cblend = (1.0 - beta5) * elastic_C(m5) + beta5 * Ceff_inv
    res["PV5b_plastic"] = bool(plastic5b)
    res["PV5b_omega"] = max(infob["wt"], infob["wc"])
    res["PV5b_rel"] = nrm(Cb - Cblend) / nrm(Cblend)
    res["PV5b_ok"] = bool(plastic5b and res["PV5b_omega"] < 1.0e-9 and res["PV5b_rel"] < 1.0e-9)

    # --- PV6: the overstress NORM grows monotonically with eta (the viscous-regularization signature) ---
    # At the clean confined-compression PLASTIC state, ||sig_visc - sig_inv|| = (1-beta)||sig_tr-sig_inv||
    # (the effective stress lags the rate-independent backbone) grows as eta increases (beta decreases)
    # and is > 0 whenever the step is plastic. The committed state is inviscid (built with eta=0).
    st6, deps6v, _, _ = _confined_compression_state(mp, Gf, Gc, lch, As)
    s_inv6, _, _, _ = return_map_tensor(st6["sig_bar"], deps6v, mp, st6["kp"])   # eta-independent inviscid backbone
    over = []
    for eta in (0.1, 0.5, 2.0, 10.0):
        me = mp_eta(eta)
        _, nst, _ = damaged_step_tensor(st6, deps6v, me, Gf, Gc, lch, As, dt=1.0)
        over.append(nrm(nst["sig_bar"] - s_inv6))           # overstress tensor norm above the backbone
    res["PV6_overstress"] = over
    res["PV6_monotone"] = bool(all(over[i + 1] > over[i] for i in range(len(over) - 1)) and over[0] > 0.0)

    ok = (res["PV1_byte_exact"] and res["PV2_converges"] and res["PV3_exact"]
          and res["PV4_is_blend"] and res["PV5a_ok"] and res["PV5b_ok"] and res["PV6_monotone"])
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf} Gc={Gc} As={As} lch={lch}")
        print(f"  PV1 eta=0 BYTE-identical to inviscid: max|dsig|={res['PV1_eta0_byte']:.1e} ({res['PV1_byte_exact']})")
        print(f"  PV2 1-D transient -> continuous closed form: errs={[f'{e:.2e}' for e in pv2_err]} "
              f"order={res['PV2_order']:.2f} ({res['PV2_converges']})")
        print(f"  PV3 [HEADLINE] 1-D EXACT steady overstress = E*eps_rate*eta = {over_exact:.4f}: "
              f"rel err={res['PV3_steady_rel_err']:.2e} ({res['PV3_exact']})")
        print(f"  PV4 tensor relaxed sig_bar == (1-beta)trial+beta inviscid: rel={res['PV4_blend_rel']:.2e} ({res['PV4_is_blend']})")
        print(f"  PV5a viscous damaged tangent analytic==FD: rel={res['PV5a_fd_rel']:.2e} ({res['PV5a_ok']})")
        print(f"  PV5b plastic={res['PV5b_plastic']} pre-onset (omega={res['PV5b_omega']:.1e}) == blended effective tangent: "
              f"rel={res['PV5b_rel']:.2e} ({res['PV5b_ok']})")
        print(f"  PV6 overstress NORM monotone in eta: {[f'{o:.3f}' for o in over]} ({res['PV6_monotone']})")
        print(f"  => P3 ETA GATE {'PASS' if ok else 'FAIL'}")
    return res


def run_p2_gate(E=30000.0, nu=0.2, fc=30.0, ft=3.0, Gf=0.1, verbose=True):
    mp = make_material(E, nu, fc, ft)
    eps0 = ft / E
    res = {}

    # UNITS: ft,E in MPa, Gf in N/mm, lch in mm => eps_f = Gf/(ft*lch) is a small strain (mm-scale
    # lch). This keeps softening at small strain, well BELOW the tensile-apex regime (the oracle's
    # return still apex-teleports in deep tension; the C++ kernel bails safe — kept apart here).
    # D1: nominal uniaxial-tension response PEAKS at ft (P1 alone is monotonic — no peak) then
    # softens to ~0 with wt -> 1.  (lch=50 => eps_f=6.7e-4 => full softening ~ eps_tot 0.005)
    d = drive_uniaxial_tension_damaged(mp, np.linspace(0, 0.008, 3000), Gf, lch=50.0)
    res["D1_peak"] = float(np.max(d["sig11"]))
    res["D1_peak_err"] = abs(res["D1_peak"] / ft - 1.0)
    res["D1_softens"] = bool(d["sig11"][-1] < 0.05 * ft and d["wt"][-1] > 0.9)
    # plasticity-alone effective stress keeps RISING (monotonic, no peak) — the peak is damage
    res["D1_eff_monotone"] = bool(d["sig_eff"][-1] > d["sig_eff"][np.argmax(d["sig11"])])
    # D0 (onset coincidence, PR #261 review): damage must initiate at kappa_p=1 / sig_eff=ft, NOT a
    # by-product. Directly gate kp at the first damaging step (was unverified — onset eps0*2 used to pass).
    di = np.argmax(d["wt"] > 1.0e-9)
    res["D0_kp_at_onset"] = float(d["kp"][di]) if d["wt"][di] > 1.0e-9 else -1.0
    res["D0_seff_at_onset"] = float(d["sig_eff"][di]) if d["wt"][di] > 1.0e-9 else 0.0
    res["D0_ok"] = bool(abs(res["D0_kp_at_onset"] - 1.0) < 0.05 and abs(res["D0_seff_at_onset"] / ft - 1.0) < 0.05)

    # D2 — crack-band SOFTENING-LAW lch-scaling (NOT independent objectivity). The damage softening
    # sig_nom = ft*exp(-eps_i/eps_f) with eps_f = Gf/(ft*lch) gives int sig_nom d eps_i = ft*eps_f =
    # Gf/lch BY CONSTRUCTION — so this verifies the eps_f WIRING (the DAMAGE dissipation is correctly
    # regularized to Gf), it does NOT independently prove FE size-objectivity (PR #261 review: that
    # integral is tautological w.r.t. eps_i). The honest FE-objectivity check is D3.
    gf_lch = {}
    for lch in (50.0, 100.0, 200.0):
        dd = drive_uniaxial_tension_damaged(mp, np.linspace(0, 0.008, 3000), Gf, lch)
        epsi = dd["epsi"]
        area = float(np.sum(0.5 * (dd["sig11"][1:] + dd["sig11"][:-1]) * np.diff(epsi)))   # trapz
        gf_lch[lch] = area * lch
    res["D2_gf_lch"] = gf_lch
    res["D2_max_rel_err"] = max(abs(gf_lch[l] / Gf - 1.0) for l in gf_lch)

    # D3 — HONEST FE-visible energy diagnostic (PR #261 review; REPORTED, not gated). The FE-visible
    # fracture energy per crack area is the nominal stress integrated over the PHYSICAL total strain,
    # times lch:  Wtot*lch = int sig_nom d eps_tot * lch. CDPM2 regularizes the DAMAGE softening ONLY
    # (D2), so the effective-PLASTICITY dissipation — which is per-VOLUME, not localized to the crack
    # band — is NOT lch-regularized and makes Wtot*lch lch-DEPENDENT (~30% over a 4x lch range here).
    # This is a CDPM2 damage-only-regularization characteristic, NOT a bug: it is small in realistic
    # regimes (damage softens before large plastic strain accrues) and shrinks if the post-peak
    # effective hardening is bounded. Reported honestly rather than hidden by the D2 construction;
    # regularizing the plastic dissipation too is a documented follow-on (P2c / ADR §4.3 [MAJOR]).
    wtot = {}
    for lch in (50.0, 100.0, 200.0):
        dd = drive_uniaxial_tension_damaged(mp, np.linspace(0, 0.008, 3000), Gf, lch)
        s = dd["sig11"]; et = dd["eps11"]
        Wel = 0.5 * float(s[-1]) ** 2 / E
        wtot[lch] = (float(np.sum(0.5 * (s[1:] + s[:-1]) * np.diff(et))) - Wel) * lch
    res["D3_wtot_lch"] = wtot
    res["D3_total_spread"] = (max(wtot.values()) - min(wtot.values())) / Gf       # FE-visible spread (reported)

    # PASS gate (honest): D0 onset at kappa_p=1; D1 the damage peak mechanism; D2 the crack-band
    # softening-law wiring (Gf/lch by construction — catches eps_f errors, not an independent
    # objectivity proof). The FE-visible TOTAL non-objectivity (D3) is a documented CDPM2
    # damage-only-regularization characteristic and is REPORTED, not gated.
    ok = (res["D0_ok"] and res["D1_peak_err"] < 0.02 and res["D1_eff_monotone"] and res["D1_softens"]
          and res["D2_max_rel_err"] < 0.02)
    res["PASS"] = bool(ok)
    if verbose:
        print(f"  E={E} nu={nu} fc={fc} ft={ft} Gf={Gf}  eps0={eps0:.3e}")
        print(f"  D0 onset: kappa_p={res['D0_kp_at_onset']:.4f} (=1) sig_eff/ft={res['D0_seff_at_onset']/ft:.4f}  ok={res['D0_ok']}")
        print(f"  D1 nominal peak = {res['D1_peak']:.4f} (target ft={ft}) err={res['D1_peak_err']:.2e}"
              f"  eff-monotone={res['D1_eff_monotone']}  softens-to-0={res['D1_softens']}")
        print(f"  D2 crack-band softening-law wiring (int sig d eps_i *lch == Gf BY CONSTRUCTION): "
              + "  ".join(f"lch={int(l)}:{gf_lch[l]:.4f}" for l in sorted(gf_lch)) + f"  err={res['D2_max_rel_err']:.2e}")
        print(f"  D3 [REPORTED, not gated] FE-visible TOTAL dissipation/crack-area (incl. un-regularized plastic): "
              + "  ".join(f"lch={int(l)}:{wtot[l]:.4f}" for l in sorted(wtot))
              + f"  spread={res['D3_total_spread']:.2f}*Gf  [CDPM2 damage-only regularization]")
        print(f"  => P2 GATE {'PASS' if ok else 'FAIL'}")
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
    print("=" * 74)
    print("LadrunoConcrete3D P1-hardening gate — CDPM2 qh1/qh2/kappa_p/ductility (Eq.18,30-36)")
    print("=" * 74)
    h = run_hardening_gate(verbose=True)
    print("-" * 74)
    print(f"HARDENING: {'PASS' if h['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P1-tangent gate — consistent tangent (spectral tensor return)")
    print("=" * 74)
    t = run_tangent_gate(verbose=True)
    print("-" * 74)
    print(f"TANGENT: {'PASS' if t['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2 gate — dual-damage (P2a: tensile wt + crack-band Gf objectivity)")
    print("=" * 74)
    p2 = run_p2_gate(verbose=True)
    print("-" * 74)
    print(f"P2: {'PASS' if p2['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2b gate — compressive damage wc + alpha_c split + crack-band Gc")
    print("=" * 74)
    p2b = run_p2b_gate(verbose=True)
    print("-" * 74)
    print(f"P2b: {'PASS' if p2b['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2c gate — unified TENSOR dual-damage split + unilateral crack-closure")
    print("=" * 74)
    p2c = run_p2c_gate(verbose=True)
    print("-" * 74)
    print(f"P2c: {'PASS' if p2c['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2d gate — single-step tensor update + damaged consistent tangent")
    print("=" * 74)
    p2d = run_p2d_gate(verbose=True)
    print("-" * 74)
    print(f"P2d: {'PASS' if p2d['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2e gate — ANALYTIC dual-projector damaged tangent (FD == numerical ref)")
    print("=" * 74)
    p2e = run_p2e_gate(verbose=True)
    print("-" * 74)
    print(f"P2e: {'PASS' if p2e['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2f gate — cyclic beta_c (Eq.50): faithful CDPM2 compressive-damage ductility")
    print("=" * 74)
    p2f = run_p2f_gate(verbose=True)
    print("-" * 74)
    print(f"P2f: {'PASS' if p2f['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2g gate — monotone (no-heal) cyclic damage + SPD unload secant tangent")
    print("=" * 74)
    p2g = run_p2g_gate(verbose=True)
    print("-" * 74)
    print(f"P2g: {'PASS' if p2g['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2h gate — compression->tension damage-coupling temper (-ctTemper modes)")
    print("=" * 74)
    p2h = run_p2h_gate(verbose=True)
    print("-" * 74)
    print(f"P2h: {'PASS' if p2h['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P2i gate — multiaxial-damage apportioning (tensile drive = E*eps_tilde, Eq.37)")
    print("=" * 74)
    p2i = run_p2i_gate(verbose=True)
    print("-" * 74)
    print(f"P2i: {'PASS' if p2i['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P3 Tier-2 IMPL-EX gate — explicit extrapolated stress + SPD secant tangent")
    print("=" * 74)
    p3i = run_p3_implex_gate(verbose=True)
    print("-" * 74)
    print(f"P3 IMPL-EX: {'PASS' if p3i['PASS'] else 'FAIL'}")
    print("=" * 74)
    print("LadrunoConcrete3D P3 Duvaut-Lions gate — viscoplastic -eta + closed-form 1-D overstress")
    print("=" * 74)
    p3e = run_p3_eta_gate(verbose=True)
    print("-" * 74)
    print(f"P3 ETA: {'PASS' if p3e['PASS'] else 'FAIL'}")
