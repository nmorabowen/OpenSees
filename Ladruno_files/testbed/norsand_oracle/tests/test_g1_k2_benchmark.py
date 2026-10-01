"""G1 gate, K2: AB06 6.1 stress-point localization benchmark (WP-144, Zone B).

Sources: equation sheet 144a section 14 (parameters, loading (S.43), criterion (S.44), the K2
gate form decided at G0) and plan 144 section 5.2 K2.

EVERY expected value and tolerance below is fixed by sheet section 14 (the ordering, the gap
band [2, 6], the nominal bands [19, 25] / [23, 29], the sweep grids) or is the published paper
number (22 / 26, sanity only). None was read from an oracle's output. The paper leaves the
initial pi_i, chi, the exact v_c0 and the crossing criterion unspecified, each of which can move
n by more than one step, so exact n is NOT a gate (sheet section 14, G0 decision 4).

Tests (both oracles, O1 = SciPy Radau rate oracle, O2 = backward-Euler + consistent tangent):
  test_k2_nominal_ordering_and_gap   G1-K2a  gate: n_0.7 < n_1.0 and gap in [2, 6]
  test_k2_nominal_bands              G1-K2b  n_0.7 in [19, 25], n_1.0 in [23, 29]
  test_k2_paper_sanity_22_26         G1-K2c  22 / 26 recorded, NOT asserted
  test_k2_sensitivity_sweep          G1-K2d  27 combinations x both criteria x both oracles,
                                             each: ordering + gap in [2, 6]

Runtime (measured on Esmeralda, serial, WP-144 venv, 2026-09-30): the six nominal tests take
about 8 s together; the 54-case sweep (2 oracles x 27 combinations, two localization runs per
case) takes about 4 min (whole file: 241 s, 60 tests). The sweep is therefore marked
@pytest.mark.slow (registered in norsand_oracle/pytest.ini); deselect with -m "not slow".
"""
from __future__ import annotations

import functools
import itertools

import pytest

from k2_sensitivity_table import (CHI_GRID, NOMINAL, PI_I0_GRID, VC0_GRID, run_case)

# ---- expected values (sheet 144a section 14; written before any oracle output was read) ----
GAP_LO, GAP_HI = 2, 6                  # gap n_1.0 - n_0.7 band [I], section 14 "K2 gate"
NOM_N07_BAND = (19, 25)                # nominal combination, first-step criterion, section 14
NOM_N10_BAND = (23, 29)
PAPER_N07, PAPER_N10 = 22, 26          # AB06 p.1551 [E]; sanity only, never asserted
RHO_07, RHO_10 = (0.7, 0.8), (1.0, 1.0)    # (rho, rho_bar), sheet 14 cases 2 and 1
CRITERIA = ("first", "interp")         # section 14 (d): index 0 = first step, 1 = interpolated

ORACLE_IDS = ("O1", "O2")


def ordering_and_gap_ok(n07, n10) -> bool:
    """The K2 gate of sheet 14: rho 0.7 localizes strictly before rho 1, gap in [2, 6]."""
    if n07 is None or n10 is None:
        return False
    return n07 < n10 and GAP_LO <= (n10 - n07) <= GAP_HI


@functools.lru_cache(maxsize=None)
def _pair(oracle, pi_i0, chi, v_c0):
    """((n_first, n_interp) for rho 0.7/0.8, same for rho = rho_bar = 1)."""
    a = run_case(oracle, pi_i0, chi, v_c0, *RHO_07)
    b = run_case(oracle, pi_i0, chi, v_c0, *RHO_10)
    return a, b


def _nominal(oracle):
    return _pair(oracle, NOMINAL["pi_i0"], NOMINAL["chi"], NOMINAL["v_c0"])


@pytest.mark.parametrize("oracle", ORACLE_IDS)
def test_k2_nominal_ordering_and_gap(oracle, record_property):
    """G1-K2a. Gate (sheet 14): nominal combination, first-step criterion."""
    a, b = _nominal(oracle)
    n07, n10 = a[0], b[0]
    record_property("n_0.7_first", n07)
    record_property("n_1.0_first", n10)
    assert n07 is not None and n10 is not None, f"{oracle}: no localization found (n_0.7={n07}, n_1.0={n10})"
    assert n07 < n10, f"{oracle}: ordering violated, n_0.7={n07} !< n_1.0={n10}"
    gap = n10 - n07
    assert GAP_LO <= gap <= GAP_HI, f"{oracle}: gap {gap} outside [{GAP_LO}, {GAP_HI}]"


@pytest.mark.parametrize("oracle", ORACLE_IDS)
def test_k2_nominal_bands(oracle):
    """G1-K2b. Nominal-combination bands of sheet 14 (first-step criterion)."""
    a, b = _nominal(oracle)
    n07, n10 = a[0], b[0]
    assert n07 is not None and n10 is not None, f"{oracle}: no localization (n_0.7={n07}, n_1.0={n10})"
    assert NOM_N07_BAND[0] <= n07 <= NOM_N07_BAND[1], f"{oracle}: n_0.7={n07} outside {NOM_N07_BAND}"
    assert NOM_N10_BAND[0] <= n10 <= NOM_N10_BAND[1], f"{oracle}: n_1.0={n10} outside {NOM_N10_BAND}"


@pytest.mark.parametrize("oracle", ORACLE_IDS)
def test_k2_paper_sanity_22_26(oracle, record_property):
    """G1-K2c. SANITY ONLY (sheet 14: exact n is not a gate). Records the nominal n against the
    paper's 22 / 26 and prints one line; asserts nothing about the values."""
    a, b = _nominal(oracle)
    for key, val in (("n_0.7_first", a[0]), ("n_1.0_first", b[0]),
                     ("n_0.7_interp", a[1]), ("n_1.0_interp", b[1]),
                     ("paper_n_0.7", PAPER_N07), ("paper_n_1.0", PAPER_N10)):
        record_property(key, val)
    print(f"\nK2 sanity {oracle}: first-step n = ({a[0]}, {b[0]}), interpolated n = "
          f"({a[1]}, {b[1]}); paper ({PAPER_N07}, {PAPER_N10})")


_SWEEP = [(o, pi0, chi, vc) for o in ORACLE_IDS
          for pi0, chi, vc in itertools.product(PI_I0_GRID, CHI_GRID, VC0_GRID)]


@pytest.mark.slow
@pytest.mark.parametrize("oracle,pi_i0,chi,v_c0", _SWEEP,
                         ids=[f"{o}-pi{p}-chi{c}-vc{v}" for o, p, c, v in _SWEEP])
def test_k2_sensitivity_sweep(oracle, pi_i0, chi, v_c0, record_property):
    """G1-K2d. Every combination of sheet 14's swept unknowns, both crossing criteria:
    n_0.7 < n_1.0 and gap in [2, 6]. A miss is a finding against the sheet or the oracle."""
    a, b = _pair(oracle, pi_i0, chi, v_c0)
    failures = []
    for k, crit in enumerate(CRITERIA):
        n07, n10 = a[k], b[k]
        record_property(f"n_0.7_{crit}", n07)
        record_property(f"n_1.0_{crit}", n10)
        if not ordering_and_gap_ok(n07, n10):
            failures.append(f"{crit}: n_0.7={n07}, n_1.0={n10}")
    assert not failures, (f"{oracle} pi_i0={pi_i0} chi={chi} v_c0={v_c0}: ordering/gap [{GAP_LO}, "
                          f"{GAP_HI}] violated for " + "; ".join(failures))
