"""Tier 1 & Tier 2 verification suite for the LadrunoHS / HSsmall Python oracle.

Covers:
  * Gate T1-A: Closed-form (alpha, beta_cap) calibration & 1D virgin oedometer
    across Loose/Dense Hostun sand (Schanz et al. 1999), Berlin sand (Schweiger
    2002), Soft/Stiff Bangkok clay (Surarak et al. 2012), and Brinkgreve et al.
    (2010) RD = 30%, 50%, 80% empirical sand suite.
  * Gate T1-B: Drained triaxial compression (CIDC) exact hyperbolic curve q(eps_1),
    ultimate failure truncation q = q_f, and exact parity of -lodeShape mc and mn
    at triaxial compression (theta = +30 deg) and extension (theta = -30 deg).
  * Gate T1-C: Rowe (1962) stress-dilatancy ratio d(eps_v^p,s)/d(gamma^p) = -sin(psi_m)
    and void-ratio dilatancy cutoff at e >= e_max.
  * Gate T1-D: Benz (2007) HSsmall modified Hardin-Drnevich secant/tangent
    backbone G_s/G_0(gamma), step-size independence, and Masing cyclic loop.
  * Gate T1-E: Benz (2007) multiaxial 180-deg and 90-deg strain-path reversals.
  * Gate T2-A: pi-plane radial Lode sweep (-30 deg to +30 deg) verifying
    Panteghini & Lagioia (2014) exact I1*I2/I3 = m_MN identity and strict
    circumscribing q_MN(theta) > q_MC(theta) on (-30 deg, +30 deg).
  * Gate T2-B: Plane-Strain Compression (PSC, eps_33 = 0) bounding comparison
    between -lodeShape mc (lower bound) and -lodeShape mn (upper bound).
  * Gate T2-C: Central finite-difference check of the 6x6 consistent algorithmic
    tangent C_alg across all active-set regimes for both mc and mn.
"""
from __future__ import annotations

import math
import numpy as np
import pytest

import ladruno_hs_reference as hs

pytestmark = [pytest.mark.zone_a, pytest.mark.t1]


# --------------------------------------------------------------------------- #
#  Gate T1-A: Auto-Calibration of (alpha, beta_cap) & 1D Virgin Oedometer     #
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "name,E50,Eoed,Eur,m,c,phi,psi,K0nc",
    [
        ("Loose_Hostun_Schanz1999", 20000.0, 16000.0, 60000.0, 0.65, 0.0, 34.0, 0.0, 0.44),
        ("Dense_Hostun_Schanz1999", 37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, 0.34),
        ("Berlin_Sand_Schweiger2002", 45000.0, 45000.0, 180000.0, 0.55, 0.0, 35.0, 5.0, 0.43),
        ("Bangkok_Soft_Clay_Surarak2012", 3000.0, 2000.0, 12000.0, 1.00, 4.0, 23.0, 0.0, 0.61),
        ("Bangkok_Stiff_Clay_Surarak2012", 15000.0, 12000.0, 60000.0, 0.60, 18.0, 28.0, 0.0, 0.53),
    ],
)
def test_t1a_auto_calibration_and_virgin_oedometer(
    name: str,
    E50: float,
    Eoed: float,
    Eur: float,
    m: float,
    c: float,
    phi: float,
    psi: float,
    K0nc: float,
) -> None:
    p = hs.HSParams(E50, Eoed, Eur, m, c, phi, psi, K0nc=K0nc, OCR=1.0)
    assert p.alpha_cap > 0.1
    assert p.beta_cap > 0.0

    # Initialize on the K0nc normal-consolidation line at sig1 = p_ref = 100 kPa:
    sig1_0 = p.p_ref
    sig3_0 = p.K0nc * (sig1_0 + p.a_apex) - p.a_apex
    st = hs.initialize_state_from_stress(
        p, np.array([sig1_0, sig3_0, sig3_0, 0.0, 0.0, 0.0]), ocr=1.0
    )

    # Apply a virgin 1D oedometric compression step (deps_1 > 0, deps_2 = deps_3 = 0):
    deps1 = 1.0e-6
    st1, C_alg, regime = hs.return_map_3d(
        p, st, np.array([deps1, 0.0, 0.0, 0.0, 0.0, 0.0])
    )
    assert regime == "shear_cap", f"{name}: expected simultaneous shear_cap in virgin oedometer"

    # Incremental stress ratio d(sig3)/d(sig1) must match K0nc and d(sig1)/d(eps1) must match Eoed_ref:
    dsig1 = st1.sig[0] - st.sig[0]
    dsig3 = st1.sig[1] - st.sig[1]
    k0_inc = dsig3 / dsig1
    eoed_inc = dsig1 / deps1
    assert abs(k0_inc - p.K0nc) < 2.0e-4, f"{name}: k0_inc={k0_inc} vs K0nc={p.K0nc}"
    assert abs(eoed_inc - p.Eoed_ref) / p.Eoed_ref < 5.0e-4, (
        f"{name}: eoed_inc={eoed_inc} vs Eoed_ref={p.Eoed_ref}"
    )
    assert abs(C_alg[0, 0] - p.Eoed_ref) / p.Eoed_ref < 1.0e-3


def test_t1a_brinkgreve_2010_rd_sweep() -> None:
    """Verify Brinkgreve, Engin & Engin (2010) empirical sand formulas across RD = 30, 50, 80%."""
    alphas = []
    eoeds = []
    for rd in (30.0, 50.0, 80.0):
        r = rd / 100.0
        E50 = 60000.0 * r
        Eoed = 60000.0 * r
        Eur = 180000.0 * r
        G0 = 60000.0 + 68000.0 * r
        m = 0.70 - rd / 320.0
        gam07 = (2.0 - r) * 1.0e-4
        phi = 28.0 + 12.5 * r
        psi = max(0.0, -2.0 + 12.5 * r)
        Rf = 1.0 - rd / 800.0
        p = hs.HSParams(
            E50,
            Eoed,
            Eur,
            m,
            0.0,
            phi,
            psi,
            Rf=Rf,
            use_small_strain=True,
            G0_ref=G0,
            gamma_07=gam07,
        )
        alphas.append(p.alpha_cap)
        eoeds.append(p.Eoed_ref)
    assert eoeds[0] < eoeds[1] < eoeds[2]
    assert all(a > 0.5 for a in alphas)


# --------------------------------------------------------------------------- #
#  Gate T1-B & T1-C: Drained Triaxial Compression (CIDC), Rowe Dilatancy,     #
#  e_max Cutoff, and Exact Coincidence of mc vs mn at TXC and TXE             #
# --------------------------------------------------------------------------- #
def _run_drained_triaxial_step(
    p: hs.HSParams, st: hs.HSState, eps1_target: float, sig3_const: float
) -> Tuple[hs.HSState, str]:
    """Solve for lateral strain eps2 = eps3 such that sig2 = sig3 = sig3_const in CIDC."""
    eps2 = st.eps[1]
    for _ in range(25):
        eps_vec = np.array([eps1_target, eps2, eps2, 0.0, 0.0, 0.0], dtype=float)
        st_try, C_alg, reg = hs.return_map_3d(p, st, eps_vec)
        res = st_try.sig[1] - sig3_const
        if abs(res) < 1.0e-10:
            return st_try, reg
        dres_deps2 = C_alg[1, 1] + C_alg[1, 2]
        eps2 -= res / dres_deps2
    raise RuntimeError("CIDC lateral stress Newton did not converge")


def test_t1b_cidc_hyperbolic_backbone_and_mc_mn_parity() -> None:
    """In CIDC at constant sig3 = p_ref, verify exact hyperbolic curve and mc == mn parity."""
    p_mc = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, K0nc=0.34, OCR=10.0, lode_shape="mc"
    )
    p_mn = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, K0nc=0.34, OCR=10.0, lode_shape="mn"
    )
    sig3_0 = 100.0
    st_mc = hs.initialize_state_from_stress(
        p_mc, np.array([sig3_0, sig3_0, sig3_0, 0.0, 0.0, 0.0]), ocr=10.0
    )
    st_mn = hs.initialize_state_from_stress(
        p_mn, np.array([sig3_0, sig3_0, sig3_0, 0.0, 0.0, 0.0]), ocr=10.0
    )

    qf = sig3_0 * (2.0 * p_mc.sin_phi / (1.0 - p_mc.sin_phi))
    qa = qf / p_mc.Rf

    eps1_grid = np.linspace(0.0, 0.04, 41)[1:]
    for eps1 in eps1_grid:
        st_mc, reg_mc = _run_drained_triaxial_step(p_mc, st_mc, float(eps1), sig3_0)
        st_mn, reg_mn = _run_drained_triaxial_step(p_mn, st_mn, float(eps1), sig3_0)

        # 1. Exact parity between -lodeShape mc and -lodeShape mn in TXC:
        assert np.max(np.abs(st_mc.sig - st_mn.sig)) < 1.0e-9
        assert abs(st_mc.gamma_p - st_mn.gamma_p) < 1.0e-12

        q_num = st_mc.sig[0] - st_mc.sig[1]
        if reg_mc == "shear":
            # Analytical gamma_p(q) at sig3 = p_ref (where s = 1):
            gam_p_ana = (2.0 / p_mc.Ei_ref) * (q_num / (1.0 - q_num / qa)) - (
                2.0 * q_num / p_mc.Eur_ref
            )
            assert abs(st_mc.gamma_p - gam_p_ana) < 1.0e-12
        elif reg_mc == "failure":
            assert abs(q_num - qf) / qf < 1.0e-10


def test_t1c_rowe_stress_dilatancy_and_emax_cutoff() -> None:
    """Verify Rowe stress-dilatancy rate and void-ratio cutoff at e >= e_max."""
    p = hs.HSParams(
        37000.0,
        29600.0,
        90000.0,
        0.50,
        0.0,
        41.0,
        14.0,
        K0nc=0.34,
        OCR=10.0,
        e0=0.55,
        emax=0.551,
    )
    sig3_0 = 100.0
    st = hs.initialize_state_from_stress(
        p, np.array([sig3_0, sig3_0, sig3_0, 0.0, 0.0, 0.0]), ocr=10.0
    )

    hit_cutoff = False
    for eps1 in np.linspace(0.0005, 0.025, 50):
        st_prev = st.copy()
        st, reg = _run_drained_triaxial_step(p, st, float(eps1), sig3_0)
        dgam_p = st.gamma_p - st_prev.gamma_p
        dp = float(np.mean(st.sig[:3]) - np.mean(st_prev.sig[:3]))
        deps_v = float(np.sum(st.eps[:3]) - np.sum(st_prev.eps[:3]))
        deps_v_e = dp / p.Kur_ref
        deps_v_ps = deps_v - deps_v_e

        if st_prev.void_ratio >= p.emax and reg in ("shear", "failure"):
            hit_cutoff = True
            assert abs(deps_v_ps) < 1.0e-12, "Dilatancy must be zero when void_ratio >= emax"
        elif reg == "shear" and dgam_p > 1.0e-7 and st_prev.void_ratio < p.emax:
            q_cur = st.sig[0] - st.sig[1]
            p_cur = float(np.mean(st.sig[:3]))
            s_psim_ana, _, _ = hs._rowe_dilatancy(p, p_cur, q_cur, st_prev.void_ratio)
            assert abs(deps_v_ps / dgam_p + s_psim_ana) < 1.0e-10

    assert hit_cutoff, "Test path must reach void_ratio >= emax dilatancy cutoff"


# --------------------------------------------------------------------------- #
#  Gate T1-D & T1-E: Benz (2007) HSsmall Backbone, Masing Loop & Reversals    #
# --------------------------------------------------------------------------- #
def test_t1d_hssmall_backbone_and_masing_hysteresis() -> None:
    """Verify step-size-independent Hardin-Drnevich backbone and Masing loop."""
    p = hs.HSParams(
        30000.0,
        30000.0,
        90000.0,
        0.55,
        0.0,
        42.0,
        12.0,
        nu_ur=0.20,
        OCR=20.0,
        use_small_strain=True,
        G0_ref=112500.0,
        gamma_07=2.0e-4,
        stiff_stress="p",
    )
    st0 = hs.initialize_state_from_stress(
        p, np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0]), ocr=20.0
    )
    st0.gamma_p = 0.05  # Pre-mobilized shear cone so small-strain cyclic shear stays paraelastic

    # 1. Compare 1-step vs 20-step simple shear to gamma = 1.5 * gamma_07 (below cutoff):
    gam_target = 1.5 * p.gamma_07
    st_1step, _, _ = hs.return_map_3d(
        p, st0, np.array([0.0, 0.0, 0.0, gam_target, 0.0, 0.0])
    )
    st_20step = st0.copy()
    for g in np.linspace(0.0, gam_target, 21)[1:]:
        st_20step, _, _ = hs.return_map_3d(
            p, st_20step, np.array([0.0, 0.0, 0.0, float(g), 0.0, 0.0])
        )
    tau_ana = p.G0_ref * gam_target / (1.0 + 0.385 * gam_target / p.gamma_07)
    assert abs(st_1step.sig[3] - tau_ana) < 1.0e-11
    assert abs(st_20step.sig[3] - tau_ana) < 1.0e-11

    # 2. Unload from +gam_target to -gam_target and reload to +gam_target (Masing 2nd rule):
    st_unl, _, _ = hs.return_map_3d(
        p, st_20step, np.array([0.0, 0.0, 0.0, -gam_target, 0.0, 0.0])
    )
    dtau_ana = -p.G0_ref * (2.0 * gam_target) / (
        1.0 + 0.385 * (2.0 * gam_target) / (2.0 * p.gamma_07)
    )
    assert abs((st_unl.sig[3] - st_20step.sig[3]) - dtau_ana) < 1.0e-11
    assert abs(st_unl.sig[3] + tau_ana) < 1.0e-11


def test_t1e_hssmall_multiaxial_180_and_90_reversals() -> None:
    """Verify 180-deg reversal restores G_0 while partial pi-plane turn partially projects H."""
    p = hs.HSParams(
        30000.0,
        30000.0,
        90000.0,
        0.55,
        0.0,
        42.0,
        12.0,
        nu_ur=0.20,
        OCR=20.0,
        use_small_strain=True,
        G0_ref=112500.0,
        gamma_07=2.0e-4,
        stiff_stress="p",
    )
    st0 = hs.initialize_state_from_stress(
        p, np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0]), ocr=20.0
    )
    st0.gamma_p = 0.05  # Pre-mobilized shear cone so small-strain probes stay paraelastic
    # Pre-strain in pure shear (e11 = +g0/2, e33 = -g0/2) at g0 = 1.0 * gamma_07 (below cutoff):
    g0 = 1.0 * p.gamma_07
    eps_pre = np.array([0.5 * g0, 0.0, -0.5 * g0, 0.0, 0.0, 0.0])
    st1, _, _ = hs.return_map_3d(p, st0, eps_pre)

    dg = 1.0e-8
    # (a) Continue forward (+dg along [0.5, 0, -0.5]):
    _, C_cont, _ = hs.return_map_3d(
        p, st1, eps_pre + np.array([0.5 * dg, 0.0, -0.5 * dg, 0.0, 0.0, 0.0])
    )
    G_cont = 0.25 * (C_cont[0, 0] - C_cont[0, 2] - C_cont[2, 0] + C_cont[2, 2])
    # (b) Reverse 180 deg (-dg along [0.5, 0, -0.5]):
    _, C_rev, _ = hs.return_map_3d(
        p, st1, eps_pre + np.array([-0.5 * dg, 0.0, 0.5 * dg, 0.0, 0.0, 0.0])
    )
    G_rev = 0.25 * (C_rev[0, 0] - C_rev[0, 2] - C_rev[2, 0] + C_rev[2, 2])
    # (c) Turn in pi-plane (de11 < 0 reverses H11 > 0 -> w1=0, while de33 < 0 continues H33 < 0 -> w3=1):
    _, C_turn, _ = hs.return_map_3d(
        p, st1, eps_pre + np.array([-0.5 * dg, dg, -0.5 * dg, 0.0, 0.0, 0.0])
    )
    G_turn = 0.25 * (C_turn[0, 0] - C_turn[0, 2] - C_turn[2, 0] + C_turn[2, 2])

    assert abs(G_rev - p.G0_ref) / p.G0_ref < 5.0e-4
    assert G_cont < G_turn < G_rev


# --------------------------------------------------------------------------- #
#  Gate T2-A & T2-B: Lode Bounding (mc vs mn) on pi-Plane & Plane Strain      #
# --------------------------------------------------------------------------- #
def test_t2a_pi_plane_lode_bounding_mc_vs_mn() -> None:
    """Sweep 25 Lode angles in [-30, +30] deg: verify I1*I2/I3 = m_MN and q_MN >= q_MC."""
    p_mc = hs.HSParams(37000.0, 29600.0, 90000.0, 0.50, 0.0, 35.0, 5.0, lode_shape="mc")
    p_mn = hs.HSParams(37000.0, 29600.0, 90000.0, 0.50, 0.0, 35.0, 5.0, lode_shape="mn")
    s = p_mn.sin_phi
    m_mn_target = (9.0 - s * s) / (1.0 - s * s)
    p_val = 100.0

    for th_deg in np.linspace(-30.0, 30.0, 25):
        th = math.radians(float(th_deg))
        x = math.sin(3.0 * th)
        Th_mc, _ = hs.lode_theta_and_deriv(p_mc, x)
        Th_mn, _ = hs.lode_theta_and_deriv(p_mn, x)
        qf_mc = p_mc.M_c * p_val / Th_mc
        qf_mn = p_mn.M_c * p_val / Th_mn

        # Reconstruct principal stresses on the Matsuoka-Nakai failure locus:
        sig1 = p_val + (2.0 / 3.0) * qf_mn * math.cos(math.pi / 6.0 - th)
        sig2 = p_val - (2.0 / 3.0) * qf_mn * math.sin(th)
        sig3 = p_val - (2.0 / 3.0) * qf_mn * math.cos(math.pi / 6.0 + th)
        I1 = sig1 + sig2 + sig3
        I2 = sig1 * sig2 + sig2 * sig3 + sig3 * sig1
        I3 = sig1 * sig2 * sig3
        assert abs(I1 * I2 / I3 - m_mn_target) < 1.0e-10

        if abs(abs(th_deg) - 30.0) < 1.0e-9:
            assert abs(qf_mn - qf_mc) < 1.0e-11
        else:
            assert qf_mn > qf_mc * 1.001


def test_t2b_plane_strain_compression_bounding_mc_vs_mn() -> None:
    """In Plane-Strain Compression (eps_3 = 0, sig_2 = const), verify mn bounds mc from above."""
    p_mc = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, K0nc=0.34, OCR=10.0, lode_shape="mc"
    )
    p_mn = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, K0nc=0.34, OCR=10.0, lode_shape="mn"
    )
    sig_conf = 100.0
    st_mc = hs.initialize_state_from_stress(
        p_mc, np.array([sig_conf, sig_conf, sig_conf, 0.0, 0.0, 0.0]), ocr=10.0
    )
    st_mn = hs.initialize_state_from_stress(
        p_mn, np.array([sig_conf, sig_conf, sig_conf, 0.0, 0.0, 0.0]), ocr=10.0
    )

    def step_psc(p_mat: hs.HSParams, st_in: hs.HSState, eps1_val: float) -> hs.HSState:
        eps2 = st_in.eps[1]
        for _ in range(50):
            # eps_33 == 0 (Plane Strain in direction 3), sig_22 == sig_conf:
            eps_v = np.array([eps1_val, eps2, 0.0, 0.0, 0.0, 0.0], dtype=float)
            st_try, C_alg, _ = hs.return_map_3d(p_mat, st_in, eps_v)
            res = st_try.sig[1] - sig_conf
            if abs(res) < 1.0e-10:
                return st_try
            deps2 = -res / max(1.0, C_alg[1, 1])
            # Damped line-search to prevent overshooting into tension cutoff near failure:
            alpha = 1.0
            for _ in range(10):
                eps_ls = np.array([eps1_val, eps2 + alpha * deps2, 0.0, 0.0, 0.0, 0.0])
                st_ls, _, _ = hs.return_map_3d(p_mat, st_in, eps_ls)
                if abs(st_ls.sig[1] - sig_conf) < abs(res):
                    break
                alpha *= 0.5
            eps2 += alpha * deps2
        raise RuntimeError("PSC lateral Newton failed")

    for eps1 in np.linspace(0.0005, 0.04, 80):
        st_mc = step_psc(p_mc, st_mc, float(eps1))
        st_mn = step_psc(p_mn, st_mn, float(eps1))
        # Matsuoka-Nakai (mn) must strictly bound Mohr-Coulomb (mc) from above in Plane Strain:
        assert st_mn.sig[0] > st_mc.sig[0]

    phi_ps_mc = math.degrees(
        math.asin((st_mc.sig[0] - st_mc.sig[1]) / (st_mc.sig[0] + st_mc.sig[1]))
    )
    phi_ps_mn = math.degrees(
        math.asin((st_mn.sig[0] - st_mn.sig[1]) / (st_mn.sig[0] + st_mn.sig[1]))
    )
    assert abs(phi_ps_mc - 41.0) < 0.5
    assert phi_ps_mn > phi_ps_mc + 2.5


# --------------------------------------------------------------------------- #
#  Gate T2-C: Central Finite-Difference Algorithmic Tangent Across All Regimes#
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("lode", ["mc", "mn"])
def test_t2c_consistent_tangent_fd_all_regimes(lode: str) -> None:
    p = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 5.0, 38.0, 10.0, K0nc=0.38, lode_shape=lode
    )
    st_ocr = hs.initialize_state_from_stress(
        p, np.array([150.0, 95.0, 65.0, 12.0, -8.0, 5.0]), ocr=15.0
    )
    st_fail = hs.initialize_state_from_stress(
        p, np.array([415.0, 95.0, 90.0, 12.0, -8.0, 5.0]), ocr=15.0
    )
    st_iso = hs.initialize_state_from_stress(
        p, np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0]), ocr=1.0
    )
    st_iso.gamma_p = 0.05
    st_k0 = hs.initialize_state_from_stress(
        p, np.array([100.0, 38.0, 38.0, 0.0, 0.0, 0.0]), ocr=1.0
    )
    st_fail_cap = hs.initialize_state_from_stress(
        p, np.array([415.0, 95.0, 90.0, 12.0, -8.0, 5.0]), ocr=1.0
    )

    cases = [
        ("elastic", st_ocr, np.array([-1.0e-4, 5.0e-5, 5.0e-5, -2.0e-5, 1.0e-5, -1.0e-5])),
        ("shear", st_ocr, np.array([1.2e-3, -4.0e-4, -5.0e-4, 3.0e-4, -1.0e-4, 1.0e-4])),
        ("failure", st_fail, np.array([2.5e-2, -1.2e-2, -1.2e-2, 3.0e-3, -1.0e-3, 1.0e-3])),
        ("cap", st_iso, np.array([2.0e-3, 1.9e-3, 1.8e-3, 1.0e-4, -5.0e-5, 5.0e-5])),
        ("shear_cap", st_k0, np.array([1.5e-3, 0.0, 0.0, 1.0e-4, -5.0e-5, 3.0e-5])),
        ("failure_cap", st_fail_cap, np.array([2.5e-2, -9.0e-3, -9.0e-3, 3.0e-3, -1.0e-3, 1.0e-3])),
        ("tension", st_ocr, np.array([-2.0e-2, -2.0e-2, -2.0e-2, 0.0, 0.0, 0.0])),
    ]

    delta = 1.0e-7
    for expected_reg, st_base, eps_target in cases:
        _, C_alg, reg = hs.return_map_3d(p, st_base, eps_target)
        assert reg == expected_reg
        C_fd = np.zeros((6, 6), dtype=float)
        for j in range(6):
            ep = eps_target.copy()
            ep[j] += delta
            em = eps_target.copy()
            em[j] -= delta
            sp, _, _ = hs.return_map_3d(p, st_base, ep)
            sm, _, _ = hs.return_map_3d(p, st_base, em)
            C_fd[:, j] = (sp.sig - sm.sig) / (2.0 * delta)
        rel_err = np.linalg.norm(C_alg - C_fd) / max(1.0, np.linalg.norm(C_fd))
        assert rel_err < 1.0e-7, f"lode={lode} regime={reg}: FD rel_err={rel_err:.2e}"

    # Also check HSsmall elastoplastic tangent:
    p_hss = hs.HSParams(
        30000.0,
        30000.0,
        90000.0,
        0.55,
        0.0,
        42.0,
        12.0,
        K0nc=0.33,
        use_small_strain=True,
        G0_ref=112500.0,
        gamma_07=2.0e-4,
        lode_shape=lode,
    )
    st_hss0 = hs.initialize_state_from_stress(
        p_hss, np.array([150.0, 90.0, 60.0, 10.0, -5.0, 4.0]), ocr=2.0
    )
    st_hss1, _, _ = hs.return_map_3d(
        p_hss, st_hss0, np.array([1.0e-4, -4.0e-5, -6.0e-5, 8.0e-5, -3.0e-5, 2.0e-5])
    )
    eps_hss2 = np.array([2.2e-4, -8.8e-5, -1.32e-4, 1.76e-4, -6.6e-5, 4.4e-5])
    _, C_alg_hss, reg_hss = hs.return_map_3d(p_hss, st_hss1, eps_hss2)
    assert reg_hss == "shear"
    C_fd_hss = np.zeros((6, 6), dtype=float)
    for j in range(6):
        ep = eps_hss2.copy()
        ep[j] += delta
        em = eps_hss2.copy()
        em[j] -= delta
        sp, _, _ = hs.return_map_3d(p_hss, st_hss1, ep)
        sm, _, _ = hs.return_map_3d(p_hss, st_hss1, em)
        C_fd_hss[:, j] = (sp.sig - sm.sig) / (2.0 * delta)
    rel_err_hss = np.linalg.norm(C_alg_hss - C_fd_hss) / max(1.0, np.linalg.norm(C_fd_hss))
    assert rel_err_hss < 1.0e-7

    # Also verify the nested 1D fallback in _solve_regime4_corner (failure_mode=False)
    # converges to the exact same root and sensitivities as 2D Newton:
    p_n, q_n, x_n, _, _, _ = hs.invariants_and_lode_grad(p, st_k0.sig)
    _, _, k3_n, _, _, _ = hs.lode_factors_and_derivs(p, x_n)
    s_n, _, _ = hs._stiffness_scale_and_derivs(p, p_n, q_n, k3_n)
    Ke_k0 = p.Kur_ref * s_n
    Ge_k0 = p.Gur_ref * s_n
    ok_2d, p_2d, q_2d, dg_2d, pp_2d, sens_2d = hs._solve_regime4_corner(
        p, 110.0, 138.5, 0.35, Ke_k0, Ge_k0, st_k0.gamma_p, st_k0.p_p, st_k0.void_ratio,
        failure_mode=False, force_1d_fallback=False,
    )
    ok_1d, p_1d, q_1d, dg_1d, pp_1d, sens_1d = hs._solve_regime4_corner(
        p, 110.0, 138.5, 0.35, Ke_k0, Ge_k0, st_k0.gamma_p, st_k0.p_p, st_k0.void_ratio,
        failure_mode=False, force_1d_fallback=True,
    )
    assert ok_2d and ok_1d
    assert abs(p_2d - p_1d) < 1.0e-10
    assert abs(q_2d - q_1d) < 1.0e-10
    assert abs(dg_2d - dg_1d) < 1.0e-12
    assert abs(pp_2d - pp_1d) < 1.0e-10
    assert np.max(np.abs(sens_2d - sens_1d)) < 1.0e-6

    # Also verify high-confinement / small-p_ref regime4_corner (ULP scaling regression check):
    p_hi = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 5.0, 38.0, 10.0, p_ref=0.1, K0nc=0.38, lode_shape=lode
    )
    ok_hi2d, p_hi2d, q_hi2d, dg_hi2d, pp_hi2d, sens_hi2d = hs._solve_regime4_corner(
        p_hi, 65000.0, 30000.0, 0.35, p_hi.Kur_ref, p_hi.Gur_ref,
        1.0e-5, 50000.0, 0.65, failure_mode=False, force_1d_fallback=False,
    )
    ok_hi1d, p_hi1d, q_hi1d, dg_hi1d, pp_hi1d, sens_hi1d = hs._solve_regime4_corner(
        p_hi, 65000.0, 30000.0, 0.35, p_hi.Kur_ref, p_hi.Gur_ref,
        1.0e-5, 50000.0, 0.65, failure_mode=False, force_1d_fallback=True,
    )
    assert ok_hi2d and ok_hi1d
    assert abs(p_hi2d - p_hi1d) / p_hi2d < 1.0e-11
    assert abs(q_hi2d - q_hi1d) / q_hi2d < 1.0e-11
    assert abs(dg_hi2d - dg_hi1d) / dg_hi2d < 1.0e-11
    assert abs(pp_hi2d - pp_hi1d) / pp_hi2d < 1.0e-11
    assert np.max(np.abs(sens_hi2d - sens_hi1d)) < 1.0e-5


def test_t1b_cide_triaxial_extension_hyperbolic_backbone_and_dilatancy() -> None:
    """In Drained Triaxial Extension (TXE, theta = -30 deg, x = -1) at constant minor
    principal stress sig3 = p_ref = 100 kPa (with sig1 = sig2 >= sig3), verify:
      1. True minor principal stress sig3_hat = p_hat - k3(-1)*q == sig3_0 = 100 kPa,
         so stiffness scaling s == 1.0 throughout TXE.
      2. Exact hyperbolic curve gamma_p(q) and asymptotic failure q_f = sig3_0 * 2*sin(phi)/(1-sin(phi)).
      3. Exact Rowe stress-dilatancy rate d(eps_v^p,s)/d(gamma^p) = -sin(psi_m) reaching -sin(psi) at failure.
      4. Exact coincidence of -lodeShape mc and -lodeShape mn in TXE.
    """
    p_mc = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, K0nc=0.34, OCR=10.0, lode_shape="mc"
    )
    p_mn = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 0.0, 41.0, 14.0, K0nc=0.34, OCR=10.0, lode_shape="mn"
    )
    sig3_0 = 100.0
    st_mc = hs.initialize_state_from_stress(
        p_mc, np.array([sig3_0, sig3_0, sig3_0, 0.0, 0.0, 0.0]), ocr=10.0
    )
    st_mn = hs.initialize_state_from_stress(
        p_mn, np.array([sig3_0, sig3_0, sig3_0, 0.0, 0.0, 0.0]), ocr=10.0
    )

    # In TXE with constant minor principal stress sig_33 = sig3_0 = 100 kPa and
    # radial compression eps_11 = eps_22 > 0 (so sig_11 = sig_22 = sig3_0 + q > sig_33):
    qf = sig3_0 * (2.0 * p_mc.sin_phi / (1.0 - p_mc.sin_phi))
    qa = qf / p_mc.Rf

    def step_txe_const_sig3(
        p_mat: hs.HSParams, st_in: hs.HSState, eps_rad: float
    ) -> Tuple[hs.HSState, str]:
        eps3 = st_in.eps[2]
        for _ in range(30):
            eps_vec = np.array([eps_rad, eps_rad, eps3, 0.0, 0.0, 0.0], dtype=float)
            st_try, C_alg, reg = hs.return_map_3d(p_mat, st_in, eps_vec)
            res = st_try.sig[2] - sig3_0
            if abs(res) < 1.0e-10:
                return st_try, reg
            eps3 -= res / max(1.0, C_alg[2, 2])
        raise RuntimeError("TXE axial stress Newton did not converge")

    for eps_rad in np.linspace(0.0, 0.025, 41)[1:]:
        st_prev = st_mc.copy()
        st_mc, reg_mc = step_txe_const_sig3(p_mc, st_mc, float(eps_rad))
        st_mn, reg_mn = step_txe_const_sig3(p_mn, st_mn, float(eps_rad))

        assert np.max(np.abs(st_mc.sig - st_mn.sig)) < 1.0e-9
        assert abs(st_mc.gamma_p - st_mn.gamma_p) < 1.0e-12

        p_cur, q_cur, x_cur, _, _, _ = hs.invariants_and_lode_grad(p_mc, st_mc.sig)
        _, _, k3_cur, _, Th_sig_cur, _ = hs.lode_factors_and_derivs(p_mc, x_cur)
        assert abs(x_cur - (-1.0)) < 1.0e-9
        assert abs(k3_cur - (2.0 / 3.0)) < 1.0e-6
        assert abs(Th_sig_cur - 1.0) < 1.0e-6

        s_cur, _, _ = hs._stiffness_scale_and_derivs(p_mc, p_cur, q_cur, k3_cur)
        assert abs(s_cur - 1.0) < 1.0e-6

        if reg_mc == "shear":
            gam_p_ana = (2.0 / p_mc.Ei_ref) * (q_cur / (1.0 - q_cur / qa)) - (
                2.0 * q_cur / p_mc.Eur_ref
            )
            assert abs(st_mc.gamma_p - gam_p_ana) < 1.0e-5
        elif reg_mc == "failure":
            assert abs(q_cur - qf) / qf < 1.0e-5
            dgam_p = st_mc.gamma_p - st_prev.gamma_p
            dp = float(np.mean(st_mc.sig[:3]) - np.mean(st_prev.sig[:3]))
            deps_v = float(np.sum(st_mc.eps[:3]) - np.sum(st_prev.eps[:3]))
            deps_v_ps = deps_v - dp / p_mc.Kur_ref
            assert abs(deps_v_ps / dgam_p + p_mc.sin_psi) < 1.0e-6


def test_t1d_hssmall_non_auxetic_constant_poisson_and_extended_masing() -> None:
    """Verify:
      1. Non-auxetic small-strain elasticity: even when G0_ref / Gur_ref = 5.0 > 3,
         K_e scales with G_e at constant nu_ur = 0.20 (never producing negative nu).
      2. Extended Masing rule: reloading beyond the historical peak gamma_max
         rejoins the virgin Hardin-Drnevich backbone (step-size independent) and
         maintains an exact central-FD consistent tangent across the transition.
    """
    p = hs.HSParams(
        25000.0,
        25000.0,
        72000.0,   # Gur_ref = 72000 / 2.4 = 30000 kPa, Kur_ref = 72000 / 1.8 = 40000 kPa
        0.50,
        0.0,
        40.0,
        10.0,
        nu_ur=0.20,
        OCR=20.0,
        use_small_strain=True,
        G0_ref=150000.0,  # G0_ref / Gur_ref = 5.0 (would be auxetic nu = -0.538 if K were unscaled!)
        gamma_07=2.0e-4,
        stiff_stress="p",
    )
    st0 = hs.initialize_state_from_stress(
        p, np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0]), ocr=20.0
    )
    st0.gamma_p = 0.10  # Keep shear cone inactive during small-strain cyclic checks

    # 1. Check effective Poisson ratio at zero strain (initial tangent) and finite strain (secant):
    _, C0_init, _ = hs.return_map_3d(p, st0, np.zeros(6))
    nu_init = C0_init[0, 1] / (C0_init[0, 0] + C0_init[0, 1])
    assert abs(nu_init - p.nu_ur) < 1.0e-14, f"Expected nu_init == {p.nu_ur}, got {nu_init}"

    st_ax, _, _ = hs.return_map_3d(p, st0, np.array([1.0e-4, 0.0, 0.0, 0.0, 0.0, 0.0]))
    dsig_ax = st_ax.sig - st0.sig
    nu_sec = dsig_ax[1] / (dsig_ax[0] + dsig_ax[1])
    assert abs(nu_sec - p.nu_ur) < 1.0e-14, f"Expected nu_sec == {p.nu_ur}, got {nu_sec}"

    # 2. Extended Masing rule: load to +gam1, unload to -gam1, reload past +gam1 to +gam2:
    gam1 = 0.8 * p.gamma_07
    gam2 = 2.0 * p.gamma_07
    st_load, _, _ = hs.return_map_3d(p, st0, np.array([0.0, 0.0, 0.0, gam1, 0.0, 0.0]))
    st_unl, _, _ = hs.return_map_3d(p, st_load, np.array([0.0, 0.0, 0.0, -gam1, 0.0, 0.0]))

    # Reload from -gam1 to +gam2 in 1 step vs 25 steps:
    st_rel_1step, _, _ = hs.return_map_3d(p, st_unl, np.array([0.0, 0.0, 0.0, gam2, 0.0, 0.0]))
    st_rel_multi = st_unl.copy()
    for g in np.linspace(-gam1, gam2, 26)[1:]:
        st_rel_multi, _, _ = hs.return_map_3d(
            p, st_rel_multi, np.array([0.0, 0.0, 0.0, float(g), 0.0, 0.0])
        )

    tau_virgin_gam2 = p.G0_ref * gam2 / (1.0 + 0.385 * gam2 / p.gamma_07)
    assert abs(st_rel_1step.sig[3] - tau_virgin_gam2) < 1.0e-10
    assert abs(st_rel_multi.sig[3] - tau_virgin_gam2) < 1.0e-10
    assert st_rel_1step.is_virgin and st_rel_multi.is_virgin

    # 3. Central-FD tangent check on a combined volumetric + deviatoric step crossing gamma_hist_max:
    eps_cross = np.array([1.5e-4, 0.5e-4, -0.5e-4, 1.3 * gam1, 0.2 * gam1, -0.1 * gam1])
    _, C_cross, _ = hs.return_map_3d(p, st_unl, eps_cross)
    delta = 1.0e-7
    C_fd = np.zeros((6, 6), dtype=float)
    for j in range(6):
        ep = eps_cross.copy()
        ep[j] += delta
        em = eps_cross.copy()
        em[j] -= delta
        sp, _, _ = hs.return_map_3d(p, st_unl, ep)
        sm, _, _ = hs.return_map_3d(p, st_unl, em)
        C_fd[:, j] = (sp.sig - sm.sig) / (2.0 * delta)
    rel_err = np.linalg.norm(C_cross - C_fd) / max(1.0, np.linalg.norm(C_fd))
    assert rel_err < 1.0e-7


@pytest.mark.parametrize("lode", ["mc", "mn"])
def test_t2d_cap_associated_normality_general_lode_angles(lode: str) -> None:
    """Verify that in REGIME_CAP at general non-triaxial Lode angles (theta != +/- pi/6),
    the plastic strain increment deps^p,c is strictly coaxial with the 3D gradient
    df^c / d(sigma') (zero normality angle error to machine precision).
    """
    p = hs.HSParams(
        37000.0, 29600.0, 90000.0, 0.50, 5.0, 38.0, 10.0, K0nc=0.38, lode_shape=lode
    )
    st_iso = hs.initialize_state_from_stress(
        p, np.array([100.0, 100.0, 100.0, 0.0, 0.0, 0.0]), ocr=1.0
    )
    st_iso.gamma_p = 0.05  # Keep shear cone inactive so pure cap regime is isolated

    w_t = np.array([1.0, 1.0, 1.0, 0.5, 0.5, 0.5])

    # Probe multiple general 3D strain increments spanning plane-strain / pure-shear Lode angles
    eps_targets = [
        np.array([2.0e-3, 1.9e-3, 1.8e-3, 1.0e-4, -5.0e-5, 5.0e-5]),
        np.array([2.2e-3, 2.0e-3, 1.8e-3, 0.0, 0.0, 0.0]),
        np.array([2.1e-3, 1.95e-3, 1.85e-3, 1.5e-4, 8.0e-5, -6.0e-5]),
    ]
    for eps_target in eps_targets:
        st_cap, _, reg = hs.return_map_3d(p, st_iso, eps_target)
        assert reg == "cap"
        p_n, q_n, x_n, _, _, _ = hs.invariants_and_lode_grad(p, st_iso.sig)
        _, _, k3_n, _, _, _ = hs.lode_factors_and_derivs(p, x_n)
        s_n, _, _ = hs._stiffness_scale_and_derivs(p, p_n, q_n, k3_n)
        Ke = p.Kur_ref * s_n
        Ge = p.Gur_ref * s_n
        dsig = st_cap.sig - st_iso.sig
        dp = float(np.mean(dsig[:3]))
        deps_e = np.array(
            [
                dp / (3.0 * Ke) + (dsig[0] - dp) / (2.0 * Ge),
                dp / (3.0 * Ke) + (dsig[1] - dp) / (2.0 * Ge),
                dp / (3.0 * Ke) + (dsig[2] - dp) / (2.0 * Ge),
                dsig[3] / Ge,
                dsig[4] / Ge,
                dsig[5] / Ge,
            ]
        )
        deps_p = eps_target - deps_e
        p_c, q_c, x_c, _, dq_c, _ = hs.invariants_and_lode_grad(p, st_cap.sig)
        assert abs(x_c) < 0.95  # Strictly non-triaxial Lode angle
        p_hat_c = p_c + p.a_apex
        grad_fc = (
            2.0 * p_hat_c * (hs._I6_M / 3.0)
            + (2.0 * q_c / (p.alpha_cap * p.alpha_cap)) * dq_c
        )
        cos_ang = float(np.sum(w_t * deps_p * grad_fc)) / math.sqrt(
            float(np.sum(w_t * deps_p * deps_p)) * float(np.sum(w_t * grad_fc * grad_fc))
        )
        assert abs(1.0 - cos_ang) < 1.0e-14

