"""sanisand_reference -- an INDEPENDENT reference integrator for SANISAND
(Dafalias & Manzari 2004), written from the paper (WP-134).

Oracle for the TIMs ring-state investigation (WP-128) and the WP-129 integrator
fix.  Pure Python (numpy + scipy); never imports OpenSees.  See README.md and
Ladruno_implementation/134_sanisand_reference_integrator.md.

    from sanisand_reference import CAMPAIGN, Options, State, integrate
    res = integrate(state, [0, 1e-4, 0, 0, 0, 0], CAMPAIGN, Options())
"""
from .model import (CAMPAIGN, TOYOURA, Options, Params, State, bounding_report,
                    e_critical, elastic_moduli, g_lode, on_yield_state, quantities,
                    yield_f)
from .integrator import Control, Result, integrate, path_table

__all__ = ["CAMPAIGN", "TOYOURA", "Options", "Params", "State", "Control",
           "Result", "integrate", "path_table", "bounding_report", "e_critical",
           "elastic_moduli", "g_lode", "on_yield_state", "quantities", "yield_f"]
