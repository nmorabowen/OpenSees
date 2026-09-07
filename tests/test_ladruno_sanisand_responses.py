"""`psi` and `yieldDistance` responses on `LadrunoSANISAND` (TIMs request
2026-09-07, F4; `_tims_proposed_model_requests_2026-09-07.md`).

Two read-only scalars the model already computes on every update and never
exposed: the state parameter `psi = e - e_c(p')` (the one that feeds M^b and
M^d) and the yield-function value `f = |s - p' alpha| - sqrt(2/3) m p'` (the
signed distance to the cone; negative inside, ~0 on it).  Both are evaluated
from the COMMITTED state at read time, so a recorder sees what fed the last
commit.

Definitions are the model's OWN (`GetStateDependent` / `GetF` verbatim):
`p' = p + p_residual`, floored at 1e-10 for psi.  With the fork's default
`p_r = 0` this is plain `e - e_c(p)`; the `p_r = 1.01` leg proves the response
follows the material's p', not a post-processed p.

Both checks are to machine precision against the SAME recorded quantities
(`stress`, `state[24]`, `alpha`), so they are identities, not calibrations.
The one measured gate is `|f|/(sqrt(2/3) m p')` on the plastic leg -- the
return-map drift -- bounded a decade above what was measured.

Deck: the confine-first cube of `test_ladruno_sanisand.py` (zero free DOFs,
no `Elastic2Plastic` M_c repair), `stdBrick`, 40 deviatoric steps.
"""
import math

import pytest

from _testbed import ops
from test_ladruno_sanisand import (_PARAMS, _P_ATM, _C_N_DEV,
                                   _build_confined, _confine_leg)

pytestmark = [pytest.mark.zone_a]

_ROOT23 = math.sqrt(2.0 / 3.0)
_SMALL = 1.0e-10           # ManzariDafalias::small, the psi p-floor
_M = _PARAMS[9]            # yield-cone size m
_IDENT_TOL = 1.0e-12       # identity checks, relative
# Measured max of |f| / (sqrt(2/3) m p') over the 40 plastic steps, both
# p_r legs: see the PR.  A decade above it; a return map that stops
# correcting drift is what this would catch.
_ON_SURFACE_REL = 1.0e-6


def _read():
    sig = ops.eleResponse(1, 'material', 1, 'stress')     # tension positive
    alpha = ops.eleResponse(1, 'material', 1, 'alpha')
    state = ops.eleResponse(1, 'material', 1, 'state')
    psi = ops.eleResponse(1, 'material', 1, 'psi')
    f = ops.eleResponse(1, 'material', 1, 'yieldDistance')
    assert len(sig) == 6 and len(alpha) == 6 and len(state) == 26
    assert len(psi) == 1 and len(f) == 1, (psi, f)
    return [-v for v in sig], list(alpha), state[24], psi[0], f[0]


def _norm_contr(v):
    return math.sqrt(v[0] ** 2 + v[1] ** 2 + v[2] ** 2
                     + 2.0 * (v[3] ** 2 + v[4] ** 2 + v[5] ** 2))


def _psi_ref(e, p_eff):
    pp = max(p_eff, _SMALL)
    e_c = _PARAMS[6] - _PARAMS[5] * (pp / _P_ATM) ** _PARAMS[7]
    return e - e_c


def _f_ref(s_comp, alpha, p_eff):
    p = (s_comp[0] + s_comp[1] + s_comp[2]) / 3.0
    dev = [s_comp[i] - (p if i < 3 else 0.0) for i in range(6)]
    r = [dev[i] - p_eff * alpha[i] for i in range(6)]
    return _norm_contr(r) - _ROOT23 * _M * p_eff


def _check_identities(pr):
    s, alpha, e, psi, f = _read()
    p = (s[0] + s[1] + s[2]) / 3.0
    assert p > 0.0, p
    p_eff = p + pr
    psi_ref = _psi_ref(e, p_eff)
    assert abs(psi - psi_ref) <= _IDENT_TOL * max(1.0, abs(psi_ref)), \
        (psi, psi_ref)
    f_ref = _f_ref(s, alpha, p_eff)
    scale = _ROOT23 * _M * p_eff
    assert abs(f - f_ref) <= _IDENT_TOL * scale, (f, f_ref, scale)
    return p_eff, f, scale


@pytest.mark.parametrize('pr', [0.0, 1.01])
def test_psi_and_yield_distance_are_the_models_own(pr):
    opts = ('-Presidual', pr)
    _build_confined('LadrunoSANISAND', 1, opts)
    _confine_leg(1)

    # Elastic isotropic state, alpha == 0: f = -sqrt(2/3) m p' exactly, and
    # psi is e - e_c(p') from the recorded void ratio.
    p_eff, f, scale = _check_identities(pr)
    assert abs(f + scale) <= _IDENT_TOL * scale, (f, -scale)
    assert f < 0.0

    # Plastic leg: on the surface to the return map's drift tolerance.
    ops.updateMaterialStage('-material', 1, '-stage', 1)
    worst = 0.0
    for step in range(_C_N_DEV):
        assert ops.analyze(1) == 0, f'deviatoric step {step + 1} failed'
        p_eff, f, scale = _check_identities(pr)
        worst = max(worst, abs(f) / scale)
    assert worst <= _ON_SURFACE_REL, worst
    # and the leg was genuinely plastic -- a state that never left the
    # elastic cone would sit at f = -scale, rel = 1, and fail the line above;
    # this pins the other direction: it did not stay near f = -scale either.
    assert abs(f) < 0.5 * scale


def test_responses_are_absent_on_vanilla_manzari_dafalias():
    """The two names are fork additions; vanilla must not answer them (an
    empty response), so a deck cannot silently read them from the wrong
    class."""
    _build_confined('ManzariDafalias', 1, ())
    _confine_leg(1)
    for name in ('psi', 'yieldDistance'):
        r = ops.eleResponse(1, 'material', 1, name)
        assert len(r) == 0, (name, r)
