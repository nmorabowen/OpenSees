"""Shared adapter for the WP-144 G1 gate tests (Zone B: numpy/scipy, fork-local).

Both oracles expose the same interface (plan §5.1), written independently from the equation
sheet (Ladruno_implementation/144a_norsand_equation_sheet.md). They differ only in two field
names and in how validation is bypassed. This adapter hides that, so a test is written once
and runs against both.

    make_params(oracle, **kw)            common names below -> that oracle's Params (validated)
    make_params(oracle, unchecked=True)  skips validation (for the forced-counterexample tests)
    ORACLES = {"O1": o1_rate, "O2": o2_algo}

Common parameter names (sheet §1.3): p0, kappa_hat, eps_v0, mu0, alpha0, M, N, N_bar, rho,
rho_bar, chi, h, lambda_tilde, v_c0, e0, lambda_c, xi, p_a, c1, c2, csl_mode, zeta, cap.
"""
import os
import sys

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))  # norsand_oracle/

import o1_rate  # noqa: E402
import o2_algo  # noqa: E402

ORACLES = {"O1": o1_rate, "O2": o2_algo}

# common name -> per-oracle field name (only where they differ)
_RENAME = {
    "O1": {},
    "O2": {"lambda_tilde": "lam_tilde", "lambda_c": "lam_c"},
}


class _O2Unchecked(o2_algo.Params):
    """O2 calls params.validate() inside its API; this subclass makes it a no-op."""

    def validate(self):  # noqa: D401
        return self


def make_params(oracle: str, unchecked: bool = False, **kw):
    """Build `oracle`'s Params from common names. Defaults are each oracle's own defaults,
    so tests MUST pass every parameter that matters for the case explicitly."""
    kw = {_RENAME[oracle].get(k, k): v for k, v in kw.items()}
    if oracle == "O1":
        if unchecked:
            kw["check"] = False
        return o1_rate.Params(**kw)
    if oracle == "O2":
        return (_O2Unchecked if unchecked else o2_algo.Params)(**kw)
    raise KeyError(oracle)


# The AB06 §6.1 / K2 parameter set (sheet §14), stated once. rho/rho_bar set per case.
K2_BASE = dict(
    p0=-100.0, kappa_hat=0.01, eps_v0=0.0, mu0=5400.0, alpha0=0.0,
    M=1.2, N=0.4, N_bar=0.2, chi=-3.5, h=280.0,
    lambda_tilde=0.0135, v_c0=1.81, csl_mode="paper", zeta="WW", cap="none",
)


@pytest.fixture(params=["O1", "O2"])
def oracle_name(request):
    return request.param


@pytest.fixture
def oracle(oracle_name):
    return ORACLES[oracle_name]
