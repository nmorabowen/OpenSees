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
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "_testbed"))
import concrete3d_ref as ref  # noqa: E402


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
