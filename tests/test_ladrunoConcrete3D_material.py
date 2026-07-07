"""
Zone-A P0 gate for LadrunoConcrete3D (CDPM2-grade solid concrete) — ADR 31.

P0 is SURFACE-ONLY: there is no nDMaterial yet (return map = P1), so this gate runs the
numpy oracle tests/_testbed/concrete3d_ref.py and asserts the Menetrey-Willam FAILURE surface
is correctly fc-normalized BEFORE any return-map / C++ code (the red-team blocking fix):

  G1  uniaxial compression sigma=(-fc,0,0) lands on f=0 EXACTLY (independent of m0, e)
  G2  uniaxial tension     sigma=(+ft,0,0) lands on f=0 EXACTLY  => fixes m0
  G3  eccentricity Lode-function identity  r(0)*e=1, r(pi/3)=1   (exact)
  G4  Kupfer equibiaxial fcc/fc reproduces target (1.16), recovering e ~ 0.52

The C++ kernel SRC/material/nD/LadrunoConcrete3DKernel.h implements the SAME surface
(yieldF/lodeR/m0Of/invariants); the g++-vs-oracle byte check lands with the return map in P1.
"""
import os
import shutil
import subprocess
import sys

import pytest

HERE = os.path.dirname(__file__)
REPO = os.path.abspath(os.path.join(HERE, os.pardir))
TESTBED = os.path.join(HERE, "_testbed")
sys.path.insert(0, TESTBED)
import concrete3d_ref as ref  # noqa: E402


# CI runs `pytest -m "zone_a"` (.github/workflows/ladruno.yml) — without this marker the WHOLE
# file (surface/return-map/hardening/tangent gates + the C++ kernel byte check) is silently
# deselected and never runs in CI. (Was missing since #240; caught in the PR #249 adversarial review.)
pytestmark = [pytest.mark.zone_a]


CASES = [(30.0, 3.0), (40.0, 3.5), (50.0, 4.0), (25.0, 2.0)]


@pytest.mark.parametrize("fc,ft", CASES)
def test_p0_surface_gate(fc, ft):
    r = ref.run_p0_gate(fc, ft, verbose=False)
    # G1/G2 — on-surface identities, machine-exact
    assert abs(r["G1_uniaxial_compression_f"]) < ref.TOL_ONSURF
    assert abs(r["G2_uniaxial_tension_f"]) < ref.TOL_ONSURF
    # G3 — eccentricity Lode identity, machine-exact
    assert r["G3_ecc_identity_err"] < ref.TOL_RATIO
    # G4 — Kupfer equibiaxial reproduced; e canonical + convex
    assert abs(r["G4_fcc_over_fc"] - 1.16) < 1.0e-6
    assert 0.5 < r["e"] <= 1.0
    assert r["PASS"]


def test_p0_meridian_ratio_rises_with_confinement():
    """Physical sanity: the surface deviatoric section rounds toward von Mises (ratio->1)
    as hydrostatic compression grows."""
    r = ref.run_p0_gate(30.0, 3.0, verbose=False)
    smr = r["surface_meridian_ratio"]
    assert smr[-1.0] < smr[-3.0] < smr[-5.0]
    assert all(0.5 < v < 1.0 for v in smr.values())


# ---------------------------------------------------------------------------
# P1 — semi-implicit Menetrey-Willam return map (perfect plasticity, qh=1):
#   H1 return-to-surface (never outside; lands on it when plastic)
#   H2/H3 uniaxial fc/ft strength
#   H4 confined triaxial fcc(sigma3): driver == analytic surface crossing, monotone up
#   H5 dilatancy knob Df reduces lateral plastic flow
# (hardening qh1/qh2 + ductility measure x(sigma) + damage = later increments)
# ---------------------------------------------------------------------------
def test_p1_return_map_gate():
    r = ref.run_p1_gate(verbose=False)
    # H1 — INDEPENDENT f (recomputed from recorded stresses, not the map's self-report)
    assert r["H1_max_f_signed"] < 1.0e-7        # never outside the surface
    assert r["H1_max_f_on_plastic"] < 1.0e-7    # lands on it when plastic
    assert r["H1_all_converged"]                # every step's Newton converged (fine + coarse)
    assert r["H1_coarse_converged"]
    assert r["H2_err"] < 1.0e-8                  # uniaxial fc (tightened from 5e-3)
    assert r["H3_err"] < 1.0e-8                  # uniaxial ft
    assert r["H4_max_rel_err"] < 1.0e-8          # return-map<->analytic-surface consistency
    assert r["H4_monotone"]                      # fcc rises with sigma3
    assert r["H6_apex_ok"]                       # hydrostatic-tension apex return
    assert r["H5_dilatancy_reduced"]             # Df<1 reduces dilation (ordering-only)
    assert r["PASS"]


def test_p1_confined_strength_near_richart():
    """INDEPENDENT physical check (the only confined gate not sharing the return-map algebra with
    its target): fcc(sigma3) sits just above Richart fc+4.1p, as expected for a model that also
    carries Lode dependence + a softer-than-Richart slope. Band tightened to observed ~1.06-1.08."""
    mp = ref.make_material(30000.0, 0.2, 30.0, 3.0, Df=1.0)
    for p in (1.5, 3.0, 6.0):
        fcc = ref.confined_strength_analytic(mp, p)
        richart = 30.0 + 4.1 * p
        assert 1.0 < fcc / richart < 1.15


def test_p0_dp_limit_e1():
    """DP limit: e=1 => circular deviatoric section => lode_r(theta,1)=1 for all theta
    (the surface degenerates to a smooth Drucker-Prager-like meridian)."""
    for th in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, ref.np.pi / 3.0]:
        assert abs(ref.lode_r(th, 1.0) - 1.0) < 1.0e-14


# ---------------------------------------------------------------------------
# P1-hardening — CDPM2 qh1/qh2/kappa_p/ductility (Grassl et al. 2013 Eq.18,30-36):
#   U1 hardening unit identities; HA reduce-to-P1 (qh0=1,Hp=0); HB failure surface at kappa_p=1;
#   HC pre-peak nonlinear hardening from qh0; HD confinement raises ductility (Eq.32-34).
# (Effective-stress plasticity is monotonic hardening — the peak/softening is DAMAGE = P2.)
# ---------------------------------------------------------------------------
def test_p1_hardening_gate():
    r = ref.run_hardening_gate(verbose=False)
    assert r["U1_ok"]                                # qh1/qh2/ductility identities
    assert r["HA_reduce_to_P1"] < 1.0e-8             # reduces to the verified perfect-plastic map
    assert r["HB_err"] < 0.02                        # failure surface reached at kappa_p=1 (sig=fc)
    assert r["HB_all_conv"]
    assert r["HC_nonlinear"]                         # yields at ~qh0*fc, hardens up monotonically
    assert r["HD_confinement_more_ductile"]          # confinement -> larger strain at kappa_p=1
    assert r["PASS"]


def test_p1_hardening_reduces_surface_to_failure():
    """yield_f with qh1=qh2=1 (Eq.18) is byte-identical to the failure surface (Eq.21)."""
    mp = ref.make_material(30000.0, 0.2, 30.0, 3.0)
    for s in ([-30.0, 0, 0], [-40, -3, -3], [3.0, 0, 0], [-50, -10, -5]):
        sv = ref.principal_to_voigt(*s)
        f_hard = ref.yield_f(sv, 30.0, 3.0, mp["e"], 1.0, 1.0)
        xi, rho, th, *_ = ref.invariants(sv)
        f_fail = ref._yf_inv(xi, rho, ref.lode_r(th, mp["e"]), mp)
        assert abs(f_hard - f_fail) < 1.0e-12


# ---------------------------------------------------------------------------
# P1-tangent — consistent (algorithmic) tangent of the spectral tensor return map:
#   T1 reduce-to-principal; T2 elastic==C_el; T3 NON-symmetric for non-associated flow;
#   T3b associated limit ~20x more symmetric (residual = semi-implicit theta-freeze => Tier-1
#   needs an unsymmetric solver UNCONDITIONALLY); T4 step-stable; T5 frame-objective.
# ---------------------------------------------------------------------------
def test_p1_tangent_gate():
    r = ref.run_tangent_gate(verbose=False)
    assert r["T1_reduce_diag"] < 1.0e-9
    assert r["T2_elastic_err"] < 1.0e-6
    assert r["T3_asym_nonassoc"] > 1.0e-1               # non-associated => non-symmetric tangent
    assert r["T3b_sym_assoc"] < 0.05                    # associated limit ~20x smaller
    assert r["T3_asym_nonassoc"] > 5.0 * r["T3b_sym_assoc"]
    # falsification: associated-limit asymmetry is the SPECTRAL RECOMPOSE (linear in shear), not the
    # theta-freeze (principal-space associated state is machine-symmetric)
    assert r["T3c_assoc_noshear_sym"] < 1.0e-6
    assert 4.0 < r["T3c_shear_linear"] < 6.0           # asymmetry linear in shear
    assert r["T4_taylor_ratio"] > 3.5                  # quadratic-Taylor convergence (~4)
    assert r["T5_objectivity"] < 1.0e-9
    assert r["PASS"]


# ---------------------------------------------------------------------------
# C++ KERNEL <-> oracle byte check (ADR §5 deliverable). Regenerate the oracle numeric-dump
# fixture, compile the standalone g++ self-check of SRC/material/nD/LadrunoConcrete3DKernel.h,
# and diff its return map + analytic consistent tangent against the oracle at the cross-platform
# tolerance floors. Skipped where g++ is unavailable (e.g. the Windows pyd test env); CI (Ubuntu
# Zone-A) runs it for real.
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# P2 — dual scalar DAMAGE (CDPM2 §2.3, Grassl 2013). P2a: tensile wt. Gates (honest, post PR #261
# adversarial review): D0 damage initiates at kappa_p=1; D1 nominal peaks at ft while the P1 effective
# stress stays monotonic (the peak IS damage) and softens to ~0; D2 the crack-band softening-law
# eps_f WIRING (int sig d eps_i *lch == Gf BY CONSTRUCTION — catches eps_f errors, NOT an independent
# objectivity proof). The FE-visible TOTAL dissipation is ~33% lch-dependent (un-regularized effective
# plasticity = a CDPM2 damage-only-regularization characteristic) and is REPORTED (D3), not gated.
# The ω-solve uses a bracketed (bisection-safeguarded) root find — a raw clamped Newton spuriously
# HEALED the cracked material (PR #261 CRITICAL).
# ---------------------------------------------------------------------------
def test_p2_damage_gate():
    r = ref.run_p2_gate(verbose=False)
    assert r["D0_ok"]                       # damage initiates at kappa_p=1 / sig_eff=ft
    assert r["D1_peak_err"] < 0.02          # nominal tensile peak == ft (the damage peak)
    assert r["D1_eff_monotone"]             # P1 effective stress monotonic (no plastic peak)
    assert r["D1_softens"]                  # softens to ~0 with wt -> 1 (no spurious healing)
    assert r["D2_max_rel_err"] < 0.02       # crack-band softening-law wiring: int sig d eps_i *lch == Gf
    assert r["PASS"]


def test_p2b_compression_damage_gate():
    """P2b: compressive damage wc + the alpha_c tension/compression split (CDPM2 Eq.37,46-57).
    C0 alpha_c -> 0 (tension) / 1 (compression) + Eq.37 == eps0 on the failure surface; C0b damage
    initiates at kappa_p=1; C1 nominal compression peaks at fc then softens; C2 the crack-band Gc
    softening-law wiring (Gc/lch BY CONSTRUCTION). FE-visible total non-objectivity (C3) reported."""
    r = ref.run_p2b_gate(verbose=False)
    assert r["C0_ok"]                       # alpha_c split: 0 tension, 1 compression
    assert r["C0_eqstrain_ok"]              # Eq.37 equivalent strain == eps0 on the failure surface
    assert r["C0b_ok"]                      # compression damage initiates at kappa_p=1
    assert r["C1_peak_err"] < 0.03          # nominal compression peak == fc (the damage peak)
    assert r["C1_eff_monotone"]             # P1 effective stress monotonic (no plastic peak)
    assert r["C1_softens"]                  # softens with wc -> 1 (no spurious healing)
    assert r["C2_max_rel_err"] < 0.02       # crack-band softening-law wiring: |sig| d eps_i *lch == Gc
    assert r["PASS"]


def test_p2c_tensor_damage_gate():
    """P2c: the UNIFIED tensor dual-damage update (CDPM2 Eq.1 spectral split) + automatic
    unilateral crack-closure. DT0 split is an exact partition + wt=wc=0 is the identity + the
    UNILATERAL routing is tested DIRECTLY (a compressive principal is carried by (1-wc) ONLY,
    wt-invariant, even under a LIVE wt=0.95) + the physical FLOOR (pure compression does NOT
    spuriously activate wt off the lateral-Newton residual, mirror for tension wc); DT1/DT2 the
    single update reduces to the validated P2a (tension, exact) / P2b (compression, to the
    onset-crossing step) 1D drivers; DT3 the end-to-end tension->compression reversal recovers the
    cracked stiffness via the per-step RE-SPLIT (ADR 4.3 BLOCKING); DT4 the damaged stress is
    frame-objective. (DT5 reports the deferred compression->tension cyclic coupling = P2f.)"""
    r = ref.run_p2c_gate(verbose=False)
    assert r["DT0_split_partition"] < 1.0e-14    # st + sc == sig exactly
    assert r["DT0_identity"] < 1.0e-14           # wt=wc=0 => nominal == effective
    assert r["DT0_unilateral"] < 1.0e-14         # compressive principals via (1-wc), tensile via (1-wt)
    assert r["DT0_compr_wt_invariant"] < 1.0e-14  # compressive entries are wt-INVARIANT (the closed crack)
    assert r["DT0_pure_compression_wt"] < 1.0e-6  # FLOOR: no spurious wt in pure compression (review-fix)
    assert r["DT0_pure_tension_wc"] < 1.0e-6      # FLOOR mirror: no spurious wc in pure tension
    assert r["DT1_sig_maxdiff"] < 1.0e-7         # pure tension reduces EXACTLY to the P2a driver
    assert r["DT1_wt_maxdiff"] < 1.0e-7
    assert r["DT2_sig_maxdiff"] < 1.0e-2         # pure compression matches P2b (to the onset step)
    assert r["DT2_wc_maxdiff"] < 1.0e-2
    assert r["DT3_recovered"]                    # unilateral crack-closure recovers full stiffness (end-to-end)
    assert r["DT4_objectivity"] < 1.0e-9         # damaged stress rotates with the strain
    assert r["PASS"]


def test_p2d_damaged_tangent_gate():
    """P2d: the SINGLE-STEP tensor constitutive update (what the C++ setTrialStrain mirrors) + its
    numerical DAMAGED CONSISTENT TANGENT (the reference the C++ analytic dual-projector tangent is
    FD-checked against). TD0 the single-step update reproduces the P2c path driver; TD1 with no damage
    the tangent reduces to the elastic C / the P1 effective tangent (the (1-w) factor and the
    -sig(x)dw rank-update both vanish pre-onset); TD2 on the softening branch the tangent is DEGRADED
    and INDEFINITE (lambda_min<0, C[0,0]<0) — the Tier-2 IMPL-EX motivation; TD3 non-symmetric for
    non-associated + damaged flow (unsymmetric solver, ADR 4.4); TD4 the update is frame-objective;
    TD5 the tangent stays finite across a load REVERSAL and near an eigenvalue crossing (the hard
    dP_T/dsig cases the C++ analytic tangent must regularize)."""
    r = ref.run_p2d_gate(verbose=False)
    assert r["TD0_tension_maxdiff"] < 1.0e-9        # single-step update == the P2c path driver (tension)
    assert r["TD0_compression_maxdiff"] < 1.0e-6    # ... and compression (eigendecompose-route floor)
    assert r["TD1a_elastic_err"] < 1.0e-6           # no-damage tangent == elastic C
    assert r["TD1b_predamage"]                      # plastic-but-pre-peak step has w_t=w_c=0
    assert r["TD1b_match_effective"] < 1.0e-6       # ... and its tangent == the P1 effective tangent
    assert r["TD2_softening_reached"]
    assert r["TD2_degraded_indefinite"]             # softening tangent: C[0,0]<0 and lambda_min(sym)<0
    assert r["TD3_asym"] > 1.0e-2                    # damaged + non-associated => non-symmetric
    assert r["TD4_objectivity"] < 1.0e-9            # the damaged update is frame-objective
    assert r["TD5_reversal_finite"]                 # tangent finite across a tension->compression reversal
    assert r["TD5_degenerate_finite"]               # ... and near an eigenvalue crossing
    assert r["PASS"]


def test_p2e_analytic_damaged_tangent_gate():
    """P2e: the ANALYTIC dual-projector damaged consistent tangent (the C++ deliverable),
    `C = D_dam : C_eff - sig_t (x) dw_t/deps - sig_c (x) dw_c/deps`, FD-verified against the P2d
    NUMERICAL reference. PE0 the spectral isotropic-function derivative (incl. the l'Hopital
    near-degenerate path); PE1-PE4 analytic == numerical in tension / confined-compression / shear
    non-associated / reversal damaged states (rel ~1e-10); PE5 reduces to the elastic C and the P1
    effective tangent before onset; PE6 at the sigma_lat=0 Macaulay kink (uniaxial-stress compression)
    the analytic tangent is a VALID subgradient — it agrees with the central difference on the loaded
    axial component, and only the ~zero-stress lateral directions differ (the kink, not a bug)."""
    r = ref.run_p2e_gate(verbose=False)
    assert r["PE0_sq"] < 1.0e-5            # spectral dY/dX vs Y=X^2
    assert r["PE0_damage"] < 1.0e-7        # ... and the damage function
    assert r["PE0_degenerate"] < 1.0e-7    # ... near-degenerate eigenvalues (l'Hopital)
    assert r["PE1_tension_rel"] < 1.0e-6 and r["PE1_wt"] > 0.5
    assert r["PE2_compression_rel"] < 1.0e-6 and r["PE2_wc"] > 0.3 and r["PE2_lam_max"] < -0.5
    assert r["PE3_shear_rel"] < 1.0e-6 and r["PE3_wt"] > 0.5
    assert r["PE4_reversal_rel"] < 1.0e-6
    assert r["PE5_elastic_rel"] < 1.0e-9
    assert r["PE5_predamage_w0"] and r["PE5_predamage_rel"] < 1.0e-6
    assert r["PE6_kink_axial_rel"] < 1.0e-3   # the loaded axial component is robust at the kink
    assert r["PASS"]


def test_p2f_gate():
    """P2f: the CDPM2 cyclic beta_c (Eq.50) restored into the compressive-damage plastic driver
    kappa_dc1 (Eq.48) — `beta_c = ft qh2 sqrt(2/3) / (rho_bar sqrt(1+2 Df^2))`. In monotonic compression
    beta_c ~ ft/(fc sqrt(1+2 Df^2)) << 1, so it makes compression markedly MORE DUCTILE than the
    beta_c=1 simplification the P2b/P2c slice used (the chosen 'faithful CDPM2' direction). F1 the
    closed-form value at the uniaxial-compression peak (= ft/(fc sqrt(1+2 Df^2)), in (0,1)); F2 the
    monotonic backbone is STILL valid (peak=fc, softens, effective-stress monotone); F3 NON-tautology —
    the real beta_c changes the post-peak stress by tens of MPa vs beta_c=1 and is a genuine
    state-dependent factor in (0,1). (The analytic damaged tangent now carries the d beta_c/d eps term,
    re-verified by test_p2e_analytic_damaged_tangent_gate; the crack-band Gc wiring re-verified by
    test_p2b_compression_damage_gate. F4 — cyclic omega_c does NOT heal on unload — is now a real gate
    (the monotone-omega_c drive fix, P2g; see test_p2g_monotone_damage_gate).)"""
    r = ref.run_p2f_gate(verbose=False)
    assert r["F1_ok"] and abs(r["F1_bc_peak"] - r["F1_bc_expect"]) < 1.0e-12 and 0.0 < r["F1_bc_peak"] < 1.0
    assert r["F2_ok"] and r["F2_peak_err"] < 0.03 and r["F2_softens"] and r["F2_eff_monotone"]
    assert r["F3_ok"] and r["F3_stress_gap"] > 0.1 * 30.0   # tens of MPa: real beta_c clearly more ductile
    assert 0.0 < r["F3_bc_min"] <= r["F3_bc_max"] < 1.0
    assert r["F4_ok"] and not r["F4_wc_heals_on_unload"]    # P2g: monotone omega_c (no healing on unload)
    assert r["PASS"]


def test_p2g_monotone_damage_gate():
    """P2g: MONOTONE (no-heal) cyclic damage. omega_t/omega_c were re-solved every step against the LIVE
    effective drive stress; the kappa-histories are already monotone, so on an elastic UNLOAD the live
    drive drops and the bracketed solve relaxes omega back (the material spuriously HEALS — #321's F4
    diagnostic). The fix drives each omega with the running MAX of its channel's drive stress, so
    omega = omega(kappa_d) is monotone and the nominal stress unloads along the degraded damage secant
    (1-omega)*sig_bar. On any monotonic path max == live, so this is byte-identical to the pre-P2g
    drivers (DT1/DT2/P2e/P2f are unaffected). G1 tension no-heal + secant unload; G2 reload below the
    previous peak retraces the secant (omega frozen); G3 compression no-heal + secant unload (the F4
    fix); G4 omega never decreases across a tension->compression reversal (persistent crack, no
    cross-heal); G5 the single-step tensor update (the C++ setTrialStrain contract) == the path driver on
    a monotonic path; G6 the UNLOAD damaged tangent is the SPD secant D_dam:C_eff (analytic == numerical
    central diff; the -sig(x)d(omega) rank-update vanished on the frozen channel)."""
    r = ref.run_p2g_gate(verbose=False)
    assert r["G1_ok"] and r["G1_wt_min_diff"] > -1.0e-9 and r["G1_secant_err"] < 1.0e-9 and 0.3 < r["G1_wt_peak"] < 0.999
    assert r["G2_ok"] and r["G2_wt_frozen"] < 1.0e-9
    assert r["G3_ok"] and r["G3_wc_min_diff"] > -1.0e-9 and r["G3_secant_err"] < 1.0e-9 and r["G3_wc_peak"] > 0.1
    assert r["G4_ok"] and r["G4_wt_min_diff"] > -1.0e-9 and r["G4_wc_min_diff"] > -1.0e-9
    assert r["G5_ok"] and r["G5_maxdiff"] < 1.0e-9
    assert r["G6_ok"] and r["G6_tangent_relerr"] < 1.0e-6 and r["G6_lambda_min"] > 0.0 and r["G6_wt_at_unload"] > 0.3
    assert r["PASS"]


def test_p2h_cttemper_gate():
    """P2h: the compression->tension damage-coupling TEMPER (`-ctTemper {none|alphat|proj}`). Per literal
    CDPM2 (Eq.43, no (1-alpha_c)) the tensile damage history accumulates in compression, so a compression
    excursion pre-damages a subsequent tension reload to ~0 (the DT5 diagnostic). The temper scales the
    tensile-history accumulation by a weight w_t: 'none'=1 (faithful, default, byte-identical); 'alphat'=
    1-alpha_c (restores tension after compression AND keeps both monotonic backbones exact); 'proj'=the
    fraction of the plastic-strain increment along positive effective-stress directions (restores tension;
    lightly softens the monotonic tension backbone). H0 'none' kills tension after compression + monotonic
    tension byte-identical; H1 'alphat' restores tension to ~ft with a byte-identical monotonic backbone;
    H2 'proj' restores tension with a <20% monotonic-tension change; H3 both keep omega monotone (P2g
    no-heal preserved); H4 analytic == numerical damaged tangent for both temper modes."""
    r = ref.run_p2h_gate(verbose=False)
    assert r["H0_ok"] and r["H0_none_tac"] < 0.2 * 3.0           # faithful: tension dies after compression
    assert r["H1_ok"] and r["H1_alphat_tac"] > 0.7 * 3.0 and r["H1_mono_maxdiff"] < 1.0e-7
    assert r["H2_ok"] and r["H2_proj_tac"] > 0.5 * 3.0 and r["H2_mono_reldiff"] < 0.2
    assert r["H3_ok"] and r["H3_alphat_wt_mono"] > -1.0e-9 and r["H3_proj_wt_mono"] > -1.0e-9
    assert r["H4_ok"] and r["H4_alphat_tan"] < 1.0e-5 and r["H4_proj_tan"] < 1.0e-4
    assert r["PASS"]


def test_p2i_multiaxial_apportioning_gate():
    """P2i: MULTIAXIAL-DAMAGE APPORTIONING. The tensile omega-solve is driven by E*eps_tilde (the CDPM2
    equivalent strain, Eq.37) instead of the extreme tensile principal; the compressive channel keeps the
    extreme principal (eps_tilde is ft-scaled, so E*eps_tilde could never reach fc to onset omega_c). USER
    decision 2026-06-21. I1 the change is INVISIBLE in uniaxial tension (E*eps_tilde == sig_bar_t), so the
    unified driver still matches the uniaxial reference byte-for-byte (DT1/D1/P2a unaffected). I2 BOTH
    equibiaxial AND triaxial tension drive E*eps_tilde above the uniaxial ft (damage onsets at a lower
    per-principal stress — the CDPM2-consistent envelope the extreme-principal model can't reproduce); the
    ordering is uni<tri<bi (hydrostatic-triaxial is APEX-CAPPED, escalates less than deviatoric biaxial).
    I3 no spurious compression->tension damage (the tension gate keeps omega_t=0 in pure compression).
    I4 the analytic damaged tangent (with its new tensile drive-gradient E*det_deps) == the numerical
    reference at an UNEQUAL biaxial-tension damaged state (off the equal-biaxial eigenvalue degeneracy)."""
    ft = 3.0
    r = ref.run_p2i_gate(verbose=False)
    assert r["I1_reduce_uniaxial"] < 1.0e-9                       # uniaxial byte-identical (invisible change)
    assert abs(r["I2_uni"] - ft) < 1.0e-6 * ft                   # Eq.38 identity E*eps_tilde == ft
    assert r["I2_bi"] > r["I2_uni"] and r["I2_tri"] > r["I2_uni"]  # both escalate above uniaxial
    assert r["I2_bi"] > r["I2_tri"]                              # deviatoric biaxial > apex-capped triaxial
    assert r["I3_pure_compression_wt"] < 1.0e-9                  # no spurious compression->tension damage
    assert r["I4_tangent_rel"] < 1.0e-6                          # analytic == numerical multiaxial tangent
    assert r["PASS"]


def test_p5_confined_fiber_gate():
    """P5: CONFINED-FIBER view (ADR 4.6) — the SAME kernel condensed against a PASSIVE hoop spring, so the
    confinement strength + ductility gain ("Mander by mechanism") EMERGE from the MW cap + non-associated
    dilatancy + ductility measure, with NO pre-baked Mander backbone. The lateral Newton targets the hoop
    residual sig_lat_eff + sig_hoop(eps_lat)=0 instead of the free (sigma3=0) target. F1 reduce-to-plain-1D
    (hoop_K=0 == the free uniaxial-stress driver, BYTE-identical). F2 confined gain: peak strength fcc>fc AND
    the strain-at-peak ductility BOTH grow monotonically with a realistic hoop stiffness (self-mobilized p>0).
    F3 THE HEADLINE: fcc/fc reproduces the Mander circular-hoop formula evaluated at the SELF-MOBILIZED
    p_conf@peak to within 5% for p/fc in [~0.1,0.2] (the mechanism analog of the 3-D active-confinement gate).
    F4 hoop YIELD caps the mobilized pressure (stirrup yields => bounded confinement). Circular/spiral hoops
    only (symmetric eps22=eps33 condensation; rectangular ties are anisotropic — a two-spring extension)."""
    r = ref.run_p5_gate(verbose=False)
    assert r["F1_reduce_1d"] < 1.0e-9                            # hoop_K=0 reduces to the free 1-D driver
    assert r["F2_strength_monotone"]                            # fcc grows with hoop stiffness (confined strength)
    assert r["F2_ductility_monotone"]                          # eps_peak grows with hoop stiffness (ductility)
    assert r["F3_npts"] >= 1 and r["F3_mander_relerr"] < 0.05   # Mander-by-mechanism within 5%
    assert r["F4_caps_below_elastic"]                          # hoop yield bounds the confining pressure
    assert r["PASS"]


def test_p3_implex_gate():
    """P3 Tier-2 IMPL-EX robustness (ADR 4.4), oracle slice (the C++ kernel port is a follow-up build
    PR). IMPL-EX reports an EXPLICIT stress from EXTRAPOLATED internal variables (plastic-strain
    increment + dual damage frozen at the committed-history rate, ratio CLAMPED to [0,R_MAX]) so the
    algorithmic tangent is the degraded-elastic secant `D_dam(w~):C0`, while committing the exact
    IMPLICIT variables. Adversarial-review-hardened (the 6 confirmed findings of #301's review):
    PI1 single-sign SPD across the snap-back where Tier-1 is INDEFINITE; PI2 committed==Tier-1
    byte-exact; PI3 SMOOTH-region O(dt) order across >=3 levels (the global-max overstress is dominated
    by an irreducible onset C0-kink lag and is NOT used as the metric — ALG-1/GAT-2); PI4 error monitor;
    PI5 the secant IS d(sig_rep)/d(deps) (FD) AND the dual-damage SPD limitation is PINNED (mixed-sign
    high-omega is NOT SPD — NUM-1, the contract is conditional not unconditional); PI6 robustness under
    NON-UNIFORM dt via the extrapolation-ratio clamp (ALG-2/NUM-2/NUM-3: no omega over-extrapolation,
    no plastic-strain blow-up, no backward damage on negative dt)."""
    r = ref.run_p3_implex_gate(verbose=False)
    assert r["softening_reached"]
    # PI1 the headline falsification (single-sign): IMPL-EX SPD, Tier-1 indefinite at the same state
    assert r["PI1_implex_SPD"] and r["PI1_implex_min_eig"] > 0.0
    assert r["PI1_tier1_indefinite"] and r["PI1_tier1_min_eig"] < 0.0
    # PI2 committed == Tier-1 to machine precision (no implicit-state corruption)
    assert r["PI2_commit_matches_tier1"]
    assert r["PI2_commit_sig_bar_err"] < 1.0e-12 and r["PI2_commit_wt_err"] < 1.0e-12
    # PI3 smooth-region convergence order >= ~1 across >=3 levels (it is ~2), monotone (not the
    # brittle single-ratio-at-one-N the review flagged); the onset global-max lag is recorded
    assert r["PI3_converges"] and r["PI3_min_order"] > 0.8 and r["PI3_smooth_monotone"]
    assert r["PI4_monitor_ok"]
    # PI5 secant == consistent tangent of the reported stress, AND the SPD scope is honest/pinned
    assert r["PI5_consistent_and_scope_pinned"]
    assert r["PI5_fd_consistency_rel"] < 1.0e-5
    assert r["PI5_singlesign_min_eig"] > 0.0 and r["PI5_mixedsign_min_eig"] < 0.0
    # PI6 the r-clamp makes non-uniform / step-cutting dt robust
    assert r["PI6_robust"]
    assert r["PASS"]


def test_p3_eta_gate():
    """P3 Duvaut-Lions viscoplastic regularization (-eta, ADR 4.4), oracle slice (the C++ kernel port is
    a follow-up build PR). The inviscid (rate-independent) plastic return is relaxed toward the elastic
    trial at the PLASTIC level by the Simo-Hughes closed form sigma=(1-beta)sig_trial+beta sig_inviscid,
    beta=dt/(eta+dt); damage then follows from the relaxed effective stress, and the effective tangent
    blends (1-beta)C_elastic+beta C_inviscid. PV1 eta=0 is BYTE-identical to the inviscid Tier-1 path;
    PV2 the scalar backward-Euler integrator converges (order ~1) to the CONTINUOUS closed-form
    overstress; PV3 the EXACT discrete steady overstress = E*eps_rate*eta (the closed-form 1-D oracle,
    dt-independent); PV4 the tensor kernel's relaxed effective stress IS the (1-beta)trial+beta inviscid
    blend (tension/compression/shear); PV5 the viscous damaged tangent matches its numerical FD (PV5a)
    and, at a genuinely PLASTIC pre-onset step, reduces to the blended effective tangent (PV5b); PV6 the
    overstress norm grows monotonically with eta (the viscous-regularization signature)."""
    r = ref.run_p3_eta_gate(verbose=False)
    # PV1 the headline byte gate: eta=0 recovers the inviscid return EXACTLY (not just to a tolerance)
    assert r["PV1_byte_exact"] and r["PV1_eta0_byte"] == 0.0
    # PV2/PV3 the closed-form 1-D overstress oracle: exact discrete steady state + order-1 transient
    assert r["PV3_exact"] and r["PV3_steady_rel_err"] < 1.0e-10
    assert r["PV2_converges"] and r["PV2_order"] > 0.8
    # PV4 the tensor kernel implements exactly the Duvaut-Lions stress blend
    assert r["PV4_is_blend"] and r["PV4_blend_rel"] < 1.0e-12
    # PV5 the viscous damaged consistent tangent (the C++ deliverable) is correct
    assert r["PV5a_ok"] and r["PV5a_fd_rel"] < 1.0e-5
    assert r["PV5b_ok"] and r["PV5b_plastic"] and r["PV5b_omega"] < 1.0e-9 and r["PV5b_rel"] < 1.0e-9
    # PV6 the regularization signature: more viscosity => more overstress above the rate-independent backbone
    assert r["PV6_monotone"]
    assert r["PASS"]


def test_p2_no_spurious_healing():
    """Regression for the PR #261 adversarial-review CRITICAL: the implicit omega solve must not
    clamp-stall to 0 on a physical softening path (a raw clamped Newton did, so the cracked material
    spontaneously HEALED — stress jumped back to full effective stress). The bracketed solver fixes it.
    Drive the workflow's failing regimes (steep softening: physical params + large lch where eps_f is
    small) and assert wt/wc never collapses 0.3+ -> ~0 while still loading."""
    import numpy as np

    def healed(w):
        return any(w[i] > 0.3 and w[i + 1] < 0.05 for i in range(len(w) - 1))

    # tension: the exact failing case (E=20000,fc=50,ft=4,Gf=0.05,lch=200) + the default-param lch=800
    mp = ref.make_material(20000.0, 0.2, 50.0, 4.0, Df=1.0)
    dt = ref.drive_uniaxial_tension_damaged(mp, np.linspace(0, 0.0006, 400), Gf=0.05, lch=200.0)
    assert not healed(dt["wt"]) and dt["wt"].max() > 0.9
    mp2 = ref.make_material(30000.0, 0.2, 30.0, 3.0, Df=1.0)
    dt2 = ref.drive_uniaxial_tension_damaged(mp2, np.linspace(0, 0.01, 3000), Gf=0.1, lch=800.0)
    assert not healed(dt2["wt"]) and dt2["wt"].max() > 0.9
    # compression: fc=50,Gc=1,lch=800 (steep) — wc must reach ~1 without healing
    dc = ref.drive_uniaxial_compression_damaged(mp, np.linspace(0, -0.05, 3000), Gc=1.0, lch=800.0, As=2.0)
    assert not healed(dc["wc"]) and dc["wc"].max() > 0.9


@pytest.mark.skipif(shutil.which("g++") is None, reason="g++ not available")
def test_cpp_kernel_matches_oracle_dump(tmp_path):
    committed = os.path.join(TESTBED, "concrete3d_oracle_fixture.txt")
    # 1) regenerate the fixture to a TMP path on THIS platform (never overwrite the committed
    #    artifact — else the diff would be self-referential and the test would dirty the tree;
    #    PR #249 review). The g++ check then runs against this SAME-PLATFORM dump.
    fixture = os.path.join(tmp_path, "fixture.txt")
    subprocess.run([sys.executable, os.path.join(TESTBED, "gen_concrete3d_fixture.py"), fixture],
                   check=True, cwd=TESTBED, env=os.environ.copy())
    # 2) THE REAL GATE: compile the self-check (header-only kernel; -I repo root for SRC/) and run it
    #    against the SAME-PLATFORM fresh dump (C++ and oracle compiled/run on one platform => the
    #    precision floors hold exactly).
    # g++ by ABSOLUTE path + env=os.environ.copy(): a bare "g++" is resolved by
    # CreateProcess against the process's LIVE native PATH, which native libs
    # (gmsh.initialize) can replace mid-battery (WinError 2 here); os.environ is
    # the intact startup snapshot the children should inherit.
    gpp = shutil.which("g++")
    exe = os.path.join(tmp_path, "c3dchk.exe")
    src = os.path.join(TESTBED, "concrete3d_kernel_check.cpp")
    subprocess.run([gpp, "-std=c++17", "-O2", "-I", REPO, src, "-o", exe],
                   check=True, cwd=REPO, env=os.environ.copy())
    out = subprocess.run([exe, fixture], cwd=REPO, capture_output=True, text=True,
                         env=os.environ.copy())
    assert out.returncode == 0, f"g++ kernel check failed:\n{out.stdout}\n{out.stderr}"
    assert "KERNEL CHECK: ALL PASS" in out.stdout
    # 3) ROT-GUARD: the committed artifact should still be ~current — but compared NUMERICALLY with a
    #    GENEROUS tolerance, never byte-wise. repr(float) is not bit-stable across platforms, and the
    #    oracle's TANGENT entries are a numerical central-difference (rel_step 1e-6) that AMPLIFIES
    #    last-ULP libm differences ~1e6x (=> ~1e-5 relative cross-platform). 1e-3 still catches real
    #    oracle drift (a genuine equation change shifts values by >>1e-3) while ignoring that noise.
    rt = open(fixture).read().split()
    ct = open(committed).read().split()
    assert len(rt) == len(ct), "committed fixture structure is stale — regenerate + commit"
    for a, b in zip(rt, ct):
        try:
            fa, fb = float(a), float(b)
        except ValueError:
            assert a == b, f"committed fixture token changed ({a!r} vs {b!r}) — regenerate + commit"
            continue
        assert abs(fa - fb) <= 1e-3 * max(1.0, abs(fb)), \
            f"committed fixture numerically stale ({a} vs {b}) — regenerate + commit"
