"""Material-point drivers on top of `integrate`: strain increments, chains, and
mixed stress/strain-controlled element tests (drained / undrained triaxial,
simple shear).  Axis convention for the element tests: the axial direction is
x (Voigt 0); compression positive throughout."""
from __future__ import annotations

import math

import numpy as np

from .integrator import Control, integrate, path_table
from .model import Options, State, e_critical, norm, t2v, I3


def isotropic_state(p0, e0):
    z = np.zeros((3, 3))
    return State(p0 * I3.copy(), z.copy(), z.copy(), float(e0), z.copy())


def increment(state, deps, P, O=None, **kw):
    """One prescribed strain increment (Voigt, engineering shear)."""
    return integrate(state, Control.strain(deps), P, O, **kw)


def chain(state, deps_list, P, O=None, **kw):
    """Consecutive increments, each from the previous end state.  Stops at the
    first non-ok status.  Returns the list of Results."""
    out = []
    st = state
    for de in deps_list:
        r = integrate(st, de if isinstance(de, Control) else Control.strain(de),
                      P, O, **kw)
        out.append(r)
        if r.status != "ok":
            break
        st = r.state
    return out


def triaxial(P, p0, e0, axial_strain, drained=True, O=None, n_out=200, **kw):
    """Monotonic triaxial compression from an isotropic state (alpha = z = 0).

    drained:   d(sigma_yy) = d(sigma_zz) = 0, shear strains 0, d(eps_xx) given;
    undrained: d(eps_vol) = 0 -> d(eps_yy) = d(eps_zz) = -d(eps_xx)/2.
    One integrate() call over the whole path (the solver chooses its steps);
    returns (Result, table) with the table sampled at the solver's steps."""
    st = isotropic_state(p0, e0)
    if drained:
        ctl = Control((True, False, False, True, True, True),
                      (axial_strain, 0.0, 0.0, 0.0, 0.0, 0.0))
    else:
        ctl = Control.strain([axial_strain, -0.5 * axial_strain, -0.5 * axial_strain,
                              0.0, 0.0, 0.0])
    res = integrate(st, ctl, P, O, **kw)
    tab = path_table(res, P, O)
    if n_out and len(tab) > n_out:
        idx = np.unique(np.linspace(0, len(tab) - 1, n_out).astype(int))
        tab = [tab[i] for i in idx]
    return res, tab


def simple_shear(P, state, gamma, drained=True, O=None, **kw):
    """Simple shear in the x-y plane: gamma_xy given, eps_xx = eps_zz = 0,
    sigma_yy constant (drained) or eps_yy = 0 (undrained / constant volume)."""
    if drained:
        ctl = Control((True, False, True, True, True, True),
                      (0.0, 0.0, 0.0, gamma, 0.0, 0.0))
    else:
        ctl = Control.strain([0.0, 0.0, 0.0, gamma, 0.0, 0.0])
    res = integrate(state, ctl, P, O, **kw)
    return res, path_table(res, P, O)


def csl_distance(P, row):
    """(eta/M_c - 1, psi) for a table row -- both 0 at critical state in
    triaxial compression (g = 1)."""
    return row["eta"] / P.Mc - 1.0, row["e"] - e_critical(row["p"], P)
