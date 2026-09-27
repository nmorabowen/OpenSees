"""WP-134 -- the independent SANISAND (DM04) reference integrator, pure-Python
gates.  No OpenSees binary needed (numpy + scipy only); see
Ladruno_implementation/134_sanisand_reference_integrator.md.

Measured wall time: ~17 s on the dev box (11 passed incl. the cross-check, 16.7 s)."""
import math
import os
import sys

import numpy as np
import pytest

pytest.importorskip("scipy")

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                                "Ladruno_scripts"))

from sanisand_reference import (CAMPAIGN, TOYOURA, Control, Options, Params,  # noqa: E402
                                State, e_critical, elastic_moduli, g_lode, integrate,
                                on_yield_state, quantities, yield_f)
from sanisand_reference.crosscheck import benign_states  # noqa: E402
from sanisand_reference.driver import triaxial  # noqa: E402
from sanisand_reference.model import dev, norm  # noqa: E402
from sanisand_reference.ring import REPRODUCER_DEPS, reproducer_state  # noqa: E402


def test_params_follow_the_opensees_argument_order():
    P = Params.from_opensees(CAMPAIGN.as_opensees())
    assert P == CAMPAIGN
    assert (CAMPAIGN.G0, CAMPAIGN.nu, CAMPAIGN.Mc, CAMPAIGN.z_max) == (264.32, 0.312885, 1.3309, 12.5)


def test_lode_interpolation_and_flow_rule_identity():
    c = 0.712
    assert g_lode(1.0, c) == pytest.approx(1.0)
    assert g_lode(-1.0, c) == pytest.approx(c)
    # B - C tr(n^3) = 1 for every unit deviatoric n: n:R' = 1 (DM04 R' design)
    rng = np.random.default_rng(3)
    st = on_yield_state(50.0, 0.5, 0.75, CAMPAIGN)
    for _ in range(20):
        a = rng.normal(size=(3, 3))
        n = dev(a + a.T)
        n /= norm(n)
        c3 = math.sqrt(6.0) * float(np.trace(n @ n @ n))
        k = (1 - c) / c
        g = g_lode(c3, c)
        B = 1 + 1.5 * k * g * c3
        C = 3 * math.sqrt(1.5) * k * g
        assert B - C * float(np.trace(n @ n @ n)) == pytest.approx(1.0, abs=1e-12)


def test_elastic_isotropic_compression_matches_the_closed_form():
    # alpha = 0, r = 0 stays 0 under isotropic compression: f < 0 throughout.
    # With G evaluated at e_init (U4) K = k sqrt(p): p = (sqrt(p0) + k eps_v / 2)^2.
    O = Options(g_void_ratio="initial")
    p0, ev = 20.0, 3.0e-3
    st = State(p0 * np.eye(3), np.zeros((3, 3)), np.zeros((3, 3)), 0.7, np.zeros((3, 3)))
    r = integrate(st, [ev / 3] * 3 + [0, 0, 0], CAMPAIGN, O)
    assert r.status == "ok" and [s["mode"] for s in r.segments] == ["elastic"]
    _, K0 = elastic_moduli(p0, 0.7, CAMPAIGN, O)
    k = K0 / math.sqrt(p0)
    assert r.end["p"] == pytest.approx((math.sqrt(p0) + 0.5 * k * ev) ** 2, rel=1e-8)


def test_plastic_increment_holds_consistency_and_stays_bounded():
    st = on_yield_state(20.0, 0.8, 0.72, CAMPAIGN)
    r = integrate(st, [0, 0, 0, 1e-3, 0, 0], CAMPAIGN, Options())
    assert r.status == "ok"
    assert abs(r.f_end) < 1e-8 * r.end["p"]
    assert r.max_abs_f_plastic < 1e-8 * r.end["p"]
    assert r.max_rho_b < 1.0


@pytest.mark.parametrize("p0", [100.0, 3000.0])
def test_undrained_toyoura_reaches_the_critical_state_line(p0):
    """DM04's e = 0.833 undrained family: dilative at 100 kPa, contractive at
    3000 kPa, both end at the SAME critical state (e fixed => p_c fixed)."""
    r, tab = triaxial(TOYOURA, p0, 0.833, 0.3, drained=False, rtol=1e-8)
    assert r.status == "ok"
    last = tab[-1]
    assert last["eta"] == pytest.approx(TOYOURA.Mc, rel=2e-3)
    assert abs(last["e"] - e_critical(last["p"], TOYOURA)) < 1e-3
    assert last["p"] == pytest.approx(1087.0, rel=5e-3)
    if p0 == 3000.0:   # contractive: q peaks then softens (quasi-steady state)
        assert max(t["q"] for t in tab) > last["q"] * 1.02
        assert min(t["p"] for t in tab) < p0


def test_drained_toyoura_dense_dilates_loose_contracts():
    rd, td = triaxial(TOYOURA, 100.0, 0.735, 0.3, drained=True, rtol=1e-7)
    rl, tl = triaxial(TOYOURA, 100.0, 0.95, 0.3, drained=True, rtol=1e-7)
    assert rd.status == rl.status == "ok"
    evd = sum(td[-1]["eps"][:3])
    evl = sum(tl[-1]["eps"][:3])
    assert evd < -0.05 and evl > 0.01
    assert max(t["eta"] for t in td) > TOYOURA.Mc * 1.1         # peak above M
    assert tl[-1]["eta"] == pytest.approx(TOYOURA.Mc, rel=1e-2)  # loose: at CS
    assert abs(tl[-1]["psi"]) < 5e-3


def test_wp128_smallest_reproducer_stays_inside_the_bounding_surface():
    """sigma = 0.0101 I, alpha = alpha_in = z = 0, one plane-strain d eps_yy = 1e-4.
    C++ ModifiedEuler (campaign) returns rho_b ~ 5.1 in one substep; the
    continuous model gives ~0.25 (WP-134 doc section 6)."""
    r = integrate(reproducer_state(), REPRODUCER_DEPS, CAMPAIGN, Options())
    assert r.status == "ok"
    assert 0.2 < r.max_rho_b < 0.3
    assert r.end["eta"] == pytest.approx(0.535, abs=0.01)


def test_uw_alpha_in_rule_is_a_toggle_that_shows_mechanism_G():
    """Benign admissible start (p 20, TC+shear, eta 0.8 Mc), txLoad 1e-5: the
    stress crosses the m-cone and reloads on the other side.  PAPER rule:
    alpha_in reseated at the onset, integrates.  UW rule (once per increment,
    no in-increment reseat): (alpha - alpha_in):n passes through 0, h -> +1e10
    -> negative, and the rate problem has NO admissible solution (H <= 0 with
    N > 0) -- the integrator stops instead of mapping it to 'elastic'."""
    st = dict(benign_states())["p20_TCshear"]
    de = Control.strain([1e-5, 0, 0, 0, 0, 0])
    rp = integrate(st, de, CAMPAIGN, Options())
    assert rp.status == "ok" and len(rp.reseats) >= 1 and rp.max_rho_b < 1.0
    ru = integrate(st, de, CAMPAIGN, Options.uw_me())
    assert ru.status == "H_nonpositive"


def test_start_outside_the_yield_surface_is_refused_by_default():
    st = on_yield_state(20.0, 0.5, 0.72, CAMPAIGN)
    st.sigma[0, 0] += 1.0                 # push the stress outside the m-cone
    assert yield_f(st.sigma, st.alpha, CAMPAIGN, Options()) > 1e-3
    r = integrate(st, [1e-5, 0, 0, 0, 0, 0], CAMPAIGN, Options())
    assert r.status == "start_outside_yield" and r.t_end == 0.0
