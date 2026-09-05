"""ADR 90 phase P0 -- Duvaut-Lions viscoplastic regularization ORACLE gate (numpy only, no build).

Two questions, both answered with numbers rather than assertions of intent:

 (1) Is the blend itself right?  PV1..PV6, cloned from the shipped `LadrunoConcrete3D -eta` gate
     (`concrete3d_ref.run_p3_eta_gate`) and run for BOTH candidate architectures over BOTH point
     models: PV1 tau = 0 byte-exact (`== 0.0`, not a tolerance); PV2 first-order convergence to
     the continuous closed form mid-transient; PV3 the EXACT discrete steady overstress
     E*eps_rate*tau, independent of dt; PV4 the blend identity; PV5 the blended tangent against a
     central finite difference; PV6 overstress monotone in tau.

 (2) Is a GENERIC `NDMaterial` wrapper -- which can only see `inner.getStress()`, hence is
     two-track (TT) -- an adequate stand-in for true Duvaut-Lions (TDL), which relaxes the
     internal variables too?  The planning brief calls this OQ3 / decision D2 / risk R1-BLOCKING.
     The oracle found a THEOREM rather than a tolerance: in 1-D, from rest, under MONOTONIC
     loading, TT and TDL are IDENTICAL for any hardening law (see `run_tt_vs_tdl_point`'s
     docstring for the proof).  The tests below therefore assert the identity where it holds
     (linear law, exponential law, J2 proportional, and on the A0 bar) AND assert that it FAILS
     where the proof's hypotheses fail (J2 non-proportional, load/unload/reload) -- the second
     half is what keeps the first from being a tautology.

 (3) A0: the 1-D softening bar.  tau = 0 is the negative control and MUST fail objectivity
     (width == h exactly, dissipated work halving with h); tau > 0 at fixed De must converge.

Zone A, numpy only -- deliberately imports no OpenSees, so it runs on a box with no built fork.
Measured wall time ~12 s (this file); the full sweep behind it is
`python3.12 tests/_testbed/run_a0_sweep.py` (~2 min) and its tables live in
`Ladruno_implementation/_adr90_a0_results.md`.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "_testbed"))
import duvaut_lions_ref as ref                                          # noqa: E402


# --------------------------------------------------------------------------------------------
# 1. PV1-PV6, both variants, both point models
# --------------------------------------------------------------------------------------------
@pytest.mark.zone_a
def test_pv_battery_both_variants():
    """PV1..PV6 for the two-track wrapper AND for true Duvaut-Lions, over the 1-D law (perfect /
    hardening / softening) and the 3-D J2 law.  PV1 is the off-switch gate ADR-59's finding B2
    demands: tau = 0 must reproduce the inviscid material EXACTLY, not to a tolerance."""
    r = ref.run_pv_gate(verbose=False)
    # PV1 -- byte-exact off switch, over a path that yields, softens and unloads
    assert r["PV1_byte_exact"] and r["PV1_max"] == 0.0
    # PV3 -- the closed-form headline: steady overstress = E * eps_rate * tau, dt-independent
    assert r["PV3_exact"]
    assert max(r["PV3_rel_err"].values()) < 1.0e-10
    # PV2 -- the discrete blend is backward Euler of the DL ODE: global order ~1 mid-transient
    assert r["PV2_converges"]
    assert all(v["order"] > 0.8 for v in r["PV2"].values())
    # PV4 -- the update IS (1-beta)*trial + beta*projection
    assert r["PV4_is_blend"] and max(r["PV4"].values()) < 1.0e-12
    # PV5 -- the blended consistent tangent matches a central FD (the C++ deliverable's check)
    assert r["PV5_ok"] and max(r["PV5"].values()) < 1.0e-5
    # PV6 -- more viscosity => more overstress above the rate-independent backbone
    assert r["PV6_monotone"]
    assert r["PASS"]


# --------------------------------------------------------------------------------------------
# 2. TT vs TDL -- the D2 evidence
# --------------------------------------------------------------------------------------------
@pytest.mark.zone_a
def test_tt_equals_tdl_where_the_theorem_holds_and_differs_where_it_does_not():
    """The measured D2 answer.  Identity where the proof applies; a LARGE difference where it
    does not.  The second set of assertions is the non-tautology guard: without them a wrapper
    that simply ignored tau would pass the first set."""
    r = ref.run_tt_vs_tdl_point(verbose=False, nsteps=700)   # reduced for the budget
    # -- identity (round-off only) --------------------------------------------------------
    assert r["perfect_plastic_identical"]                     # perfect plasticity
    assert r["linear_max_d_path_rel"] < 1.0e-10               # linear H/E in [-0.10, +0.10]
    assert r["exp_max_d_path_rel"] < 1.0e-10                  # NONLINEAR exponential law
    assert r["j2_prop_max_d_rel"] < 1.0e-10                   # J2, proportional path
    # -- refutation: the identity is NOT general ------------------------------------------
    assert r["j2_nonprop_max_d_rel"] > 1.0e-3                 # measured 4.4e-2
    assert r["cyclic_max_d_path_rel"] > 1.0e-2                # measured 3.3e-1
    # and the difference grows with De on the non-proportional path
    npr = sorted((x["De"], x["d_rel"]) for x in r["j2_rows"]
                 if x["path"] == "nonprop" and x["H_over_E"] < 0)
    assert npr[0][1] <= npr[1][1] <= npr[2][1]


# --------------------------------------------------------------------------------------------
# 3. A0 bar -- small version of the sweep
# --------------------------------------------------------------------------------------------
@pytest.mark.zone_a
def test_a0_negative_control_tau0_is_not_objective():
    """H1.  With tau = 0 the softening bar localizes into exactly ONE element at every mesh --
    the band width is h, not a material length -- and the dissipated work halves with h.  This is
    the discriminating negative control ADR-31's G5 pairing requires: if this passed objectivity,
    nothing the positive gate says about tau > 0 would mean anything."""
    prev_w = prev_W = None
    for N in (20, 40, 80):
        r = ref.a0_bar(variant="tt", N=N, tau=0.0, T=1.0, nsteps=2000, u_max=2.0,
                       imperfection="one_element")
        assert r["converged"], f"tau=0 N={N} did not converge"
        # width IS the element size, by all three definitions, EXACTLY
        assert r["w1"] == pytest.approx(r["h"], rel=1e-12)
        assert r["w2"] == pytest.approx(r["h"], rel=1e-12)
        assert r["w3"] == pytest.approx(r["h"], rel=1e-12)
        if prev_w is not None:
            assert r["w2"] / prev_w == pytest.approx(0.5, rel=1e-9)
            assert prev_W / r["W_50"] == pytest.approx(2.0, rel=0.02)   # measured 2.000
        prev_w, prev_W = r["w2"], r["W_50"]


@pytest.mark.zone_a
def test_a0_positive_gate_tau_gt_0_converges():
    """H2.  At a FIXED De the width and the load-displacement curve converge under refinement,
    where the tau = 0 control gave an exact factor 0.5 per refinement.  De = 1e-3 with the
    mesh-convergent (graded, fixed-physical-length) imperfection; measured w2 ratios
    N=40/20 = 1.060, N=80/40 = 1.017."""
    prev = None
    for N in (20, 40, 80):
        r = ref.a0_bar(variant="tt", N=N, tau=1.0e-3, T=1.0, nsteps=2000, u_max=2.0,
                       imperfection="fixed_graded")
        assert r["converged"]
        assert r["inelastic_ratio"] > 1.0, "run must be genuinely inelastic (non-tautology)"
        if prev is not None:
            ratio = r["w2"] / prev["w2"]
            assert 0.95 < ratio < 1.10, f"N={N}: w2 ratio {ratio}"
            assert ref.curve_distance(r, prev) < 0.05
        prev = r
    # and the converged width is a strong function of De -- monotone, and NOT the imperfection
    # width (10 mm): 3e-4 -> ~3.6 mm, 1e-3 -> ~42 mm, 3e-3 -> ~91 mm at N = 80.
    ws = [ref.a0_bar(variant="tt", N=80, tau=De, T=1.0, nsteps=2000, u_max=2.0,
                     imperfection="fixed_graded")["w2"] for De in (3.0e-4, 1.0e-3, 3.0e-3)]
    assert ws[0] < ws[1] < ws[2]
    assert ws[0] < 10.0 < ws[1]


@pytest.mark.zone_a
def test_a0_tt_and_tdl_agree_on_the_bar():
    """H4 on the structure.  The bar loads its band monotonically, so the 1-D theorem applies and
    the generic two-track wrapper reproduces true Duvaut-Lions to ~1e-5 (the residue comes from
    elements that yield and then unload, plus the Newton tolerance)."""
    for De, N in ((0.0, 40), (3.0e-4, 40), (1.0e-3, 80)):
        a = ref.a0_bar(variant="tt", N=N, tau=De, T=1.0, nsteps=2000, u_max=2.0,
                       imperfection="fixed_graded" if De else "one_element")
        b = ref.a0_bar(variant="tdl", N=N, tau=De, T=1.0, nsteps=2000, u_max=2.0,
                       imperfection="fixed_graded" if De else "one_element")
        assert a["converged"] and b["converged"]
        assert abs(a["w2"] - b["w2"]) / max(b["w2"], 1e-30) < 1.0e-4
        assert abs(a["W_end"] - b["W_end"]) / max(abs(b["W_end"]), 1e-30) < 1.0e-4
        assert ref.curve_distance(a, b) < 1.0e-4


@pytest.mark.zone_a
def test_a0_two_gates_that_pass_tautologically():
    """The findings the fork records rather than hides -- both are gates that CANNOT fail, so
    neither may be cited as evidence of regularization.

    (a) De collapse at a fixed step count.  beta = dt/(tau+dt) = 1/(1 + nsteps*De) and the strain
        increment u_max/nsteps carry no other dependence on tau or T, so two runs with the same
        (De, nsteps) are BIT-identical whatever (tau, T) produced that De.
    (b) A FLAT fixed-length imperfection (the brief's convention (b) read literally: 2/4/8/16
        elements all at 0.9*sigY).  Every element of the weak zone is exactly as weak as its
        neighbours, so the continuum solution is piecewise constant with the zone boundary on a
        mesh line -- every mesh represents it EXACTLY and the answer is bit-identical in N.
        A graded notch over the same physical length (`fixed_graded`) is the honest version.
    """
    ref_run = ref.a0_bar(variant="tt", N=80, tau=1.0e-3, T=1.0, nsteps=1000, u_max=2.0,
                         imperfection="fixed_graded")
    for tau, T in ((1.0e-2, 10.0), (1.0e-4, 0.1)):
        r = ref.a0_bar(variant="tt", N=80, tau=tau, T=T, nsteps=1000, u_max=2.0,
                       imperfection="fixed_graded")
        assert r["w2"] == ref_run["w2"] and r["P_peak"] == ref_run["P_peak"]   # bit-identical
        assert r["W_end"] == ref_run["W_end"]

    base = ref.a0_bar(variant="tt", N=20, tau=1.0e-3, T=1.0, nsteps=1000, u_max=2.0,
                      imperfection="fixed_flat")
    for N in (40, 80, 160):
        r = ref.a0_bar(variant="tt", N=N, tau=1.0e-3, T=1.0, nsteps=1000, u_max=2.0,
                       imperfection="fixed_flat")
        # identical to round-off (w2 sums a different number of terms per mesh, so the last
        # couple of ulps move; every reported digit agrees)
        assert r["w2"] == pytest.approx(base["w2"], rel=1e-12)
        assert r["P_peak"] == pytest.approx(base["P_peak"], rel=1e-12)
        assert r["W_end"] == pytest.approx(base["W_end"], rel=1e-12)


@pytest.mark.zone_a
def test_a0_large_de_is_frozen_not_regularized():
    """The brief asks for De in {0.01, 0.1, 1}.  Measured on this deck, all three are far past
    the point where the bar still localizes, and De = 1 never softens at all (the load rises
    monotonically to 265 N, 14x the rate-independent peak of 18 N).  "Width converges" there is
    a statement about an essentially elastic bar.  Recorded so the ADR cannot quote it."""
    r = ref.a0_bar(variant="tt", N=80, tau=1.0, T=1.0, nsteps=2000, u_max=2.0,
                   imperfection="fixed_graded")
    assert r["converged"]
    assert r["k_peak"] == len(r["P"]) - 1, "De=1 must have NO post-peak branch at all"
    assert r["P_peak"] > 10.0 * 18.0
    assert r["w2"] == 0.0                       # no post-peak plastic increment to measure
    for De in (0.01, 0.1):
        q = ref.a0_bar(variant="tt", N=80, tau=De, T=1.0, nsteps=2000, u_max=2.0,
                       imperfection="fixed_graded")
        assert q["w3"] == pytest.approx(100.0)  # the "band" is the whole 100 mm bar


@pytest.mark.zone_a
def test_a0_no_snapback_by_construction():
    """The A0 parameter choice, asserted rather than asserted-about.  A softening bar under
    end-displacement control snaps back (and Newton cannot follow it) when the band is shorter
    than L|H|/E.  H = -E/400 puts that threshold at 0.25 mm, a factor 2.5 below the finest
    element (h = 0.625 mm at N = 160), so no leg of the sweep snaps back."""
    L, E, H = 100.0, 20000.0, -50.0
    assert L * abs(H) / E == pytest.approx(0.25)
    assert L / 160 > 2.0 * L * abs(H) / E
    r = ref.a0_bar(variant="tt", N=160, tau=0.0, T=1.0, nsteps=2000, u_max=2.0,
                   imperfection="one_element")
    assert r["converged"]
    assert np.all(np.diff(r["u"]) > 0.0)        # monotone displacement control, no reversal


@pytest.mark.zone_a
def test_band_width_metric_is_calibrated():
    """OQ5 / risk R4.  w2 is threshold-FREE; its within-element variance term h^2/12 is what
    makes a single-element band read EXACTLY h and a k-element top hat read EXACTLY k*h, so the
    number is comparable across meshes.  Without that term a one-element band reads 0."""
    for h in (0.5, 1.0, 2.5):
        p = np.zeros(40)
        p[17] = 3.0
        assert ref.band_widths(p, h)[1] == pytest.approx(h, rel=1e-12)
        for k in (2, 4, 8):
            q = np.zeros(40)
            q[12:12 + k] = 3.0
            assert ref.band_widths(q, h)[1] == pytest.approx(k * h, rel=1e-12)
        assert ref.band_widths(np.zeros(40), h) == (0.0, 0.0, 0.0)


# --------------------------------------------------------------------------------------------
# 4. P0b -- the four measurements the 3-lens adversarial pass showed were missing
#
# Reduced settings (fewer De, fewer steps) so the whole file stays inside its budget; the full
# tables are `python3.12 tests/_testbed/run_a0_sweep.py --p0b` and live in
# `Ladruno_implementation/_adr90_p0b_results.md`.  Thresholds here are order-of-magnitude, so the
# reduction cannot flip a verdict.
# --------------------------------------------------------------------------------------------
@pytest.mark.zone_a
def test_p0b_a_state_dependent_elastic_operator_breaks_the_identity():
    """P0b-(a) / constitutive critic finding 1 (BLOCKING).  The ADR-90 theorem has a hypothesis
    nobody wrote down: a CONSTANT elastic operator.  `psi = sigma + E*alpha` can only advance by
    `E*de` if E is the same number on both tracks -- but a generic wrapper's only elastic operator
    is `inner.getInitialTangent()`, i.e. E at the INNER's stress, while it must predict from its
    OWN stress.  With `E(sigma)` the identity fails on EVERY path, monotonic and proportional
    included, and the error grows monotonically with De.

    SANISAND (G(p), K(p)), Drucker-Prager with G(p), and every hypoelastic inner are in this
    class, so this is not a corner case for the named consumer -- it is the normal case."""
    r = ref.run_pressure_dependent_leg(nsteps=600, De_list=(1.0e-3, 1.0e-2, 0.1, 0.3),
                                       verbose=False)
    assert r["const_max"] < 1.0e-10, "constant E must still satisfy the theorem"
    assert r["var_max"] > 1.0e-2, "a state-dependent E must break it measurably"
    # monotone in De, on the buildable wrapper
    var = sorted((x["De"], x["d_path_rel"]) for x in r["rows"] if x["E_beta"] != 0.0)
    assert all(var[i][1] <= var[i + 1][1] for i in range(len(var) - 1))
    # the obvious repair -- evaluate the predictor modulus at the WRAPPER's own stress -- helps
    # by ~1-2 orders but does NOT close the gap, and is not implementable at the NDMaterial seam
    # (it needs the function E(.), not one sampled number).
    assert r["var_max_ttvp"] < r["var_max"]
    assert r["var_max_ttvp"] > 1.0e-3


@pytest.mark.zone_a
def test_p0b_b_dissipation_goes_negative_on_unloading():
    """P0b-(b).  In true Duvaut-Lions `sigma_bar` is the closest-point projection OF `sigma_vp`,
    so `sigma_vp - sigma_bar` lies in the normal cone and `D_w = sigma_vp . dEps^vp >= 0` follows
    from convexity.  The two-track wrapper projects a DIFFERENT trial, so nothing enforces that
    geometry.  Measured: the increment is genuinely negative on load/unload/reload -- up to
    -0.8 % of the run's cumulative dissipation at De = 0.1 -- while the non-proportional J2 and
    swipe paths stay clean to round-off."""
    r = ref.run_dissipation_gate(nsteps=600, De_list=(1.0e-2, 0.1), verbose=False)
    assert r["any_negative_significant"], "the violation must be reproducible"
    assert all(k.startswith("cyclic1d") for k in r["paths_with_significant_negatives"]), \
        "as measured, only the UNLOADING paths violate it"
    assert r["neg_sum_rel_worst"] < -1.0e-4          # measured -7.9e-3 on the full sweep
    # every run still dissipates net-positive, and the inner's own dissipation is positive
    for v in r["rows"].values():
        assert v["Dw_cum"] > 0.0 and v["D_inner"] > 0.0


@pytest.mark.zone_a
def test_p0b_c_the_converged_width_is_an_imperfection_property():
    """P0b-(c).  At the ONE De where A0 reported both w2 and w3 converging, vary the imperfection.
    The converged width tracks the imperfection amplitude and the zone length, so it is NOT a
    length set by tau: viscosity stops the band from collapsing onto one element, but what it
    collapses TO is chosen by the imperfection field."""
    r = ref.run_imperfection_study(Ns=(80,), amps=(0.05, 0.10, 0.20), lens=(0.05, 0.10, 0.20),
                                   nsteps=1000, verbose=False)
    assert all(x["converged"] for x in r["rows"])
    assert r["w2_span_over_amp"] > 2.0      # measured x4.28 on the full sweep
    assert r["w2_span_over_len"] > 2.0      # measured x3.47
    # w2 is bounded above by L by construction, so "w2 near L" means NO band, not a wide one
    assert all(x["w2"] <= 100.0 + 1e-9 for x in r["rows"])


@pytest.mark.zone_a
def test_p0b_d_blended_acoustic_tensor_restores_ellipticity():
    """P0b-(d).  The only measurement that speaks to WELL-POSEDNESS for the consumer's material
    class; the 1-D A0 bar cannot see it.  On a plane-strain NON-ASSOCIATED softening
    Drucker-Prager point model the inviscid algorithmic tangent loses ellipticity
    (min det Q < 0, Rudnicki & Rice 1975) and the blended tangent regains it below a critical
    BETA -- which means the criterion is on `Delta_t/tau`, not on De: the same De is elliptic or
    not depending on the step count.  Also: TT == TDL on a PROPORTIONAL path even though the flow
    is non-associated (non-associativity alone does not break the theorem), and differs on a
    ROTATING path."""
    r = ref.run_dp_leg(nsteps=400, De_list=(3.0e-4,), n_shear=300, ntheta=121,
                       verbose=False)
    assert r["min_det_inviscid"] < 0.0, "the inviscid non-associated tangent must lose ellipticity"
    assert 0.0 < r["beta_crit"] < 1.0
    # non-associativity alone does NOT break the identity; path rotation does
    assert r["prop_max"] < 1.0e-10
    assert r["rot_max"] > 1.0e-4
    rot = sorted((x["De"], x["d_path_rel"]) for x in r["tt_vs_tdl"] if x["path"] == "rotating")
    assert all(rot[i][1] <= rot[i + 1][1] for i in range(len(rot) - 1))
    # the ellipticity criterion is a BETA criterion: the De it corresponds to scales as 1/nsteps
    des = r["De_crit_at_nsteps"]
    assert des[250] == pytest.approx(4.0 * des[1000], rel=1e-9)
