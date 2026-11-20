"""WP-144 round 3b, G2 (Zone B, SLOW): the floor sensitivity report "F vs F/2" (sheet 9.7, plan 2.7, deck report K5) on the mini strip-footing BVP of
floor_sensitivity.py.  Run with --runslow (about 2-3 minutes: two runs of ~1 minute).

Acceptance (the plan's, restated in the sheet 9.7 "Counted"): the limit load of the footing moves by LESS THAN 2 % when p_min is halved.  Two
non-vacuity statements come with it, both from the sheet: the floor is actually engaged on this model (Gauss points AT the floor, floor events, an
energy source W_f > 0 -- else the report would prove nothing), and the energy source VANISHES with p_min (sum W_f at F/2 is smaller than at F: the sheet's
"the floor's violation is exactly E_f, bounded by W_f, counted and reported, vanishing with p_min").  No refusal may occur (the floor never refuses)
and the footing must reach the full displacement in both runs.

KILLS: a floor that changes the limit load by more than 2 % (M-F9: the floor inside the local Newton, a stiffness regularisation hidden in the tangent), a
floor that is silently inactive on the model, a refusal where the floor should act (M-F1), counters that do not reach the Gauss-point responses (M-F2).
Mutation gate: this test is part of the `-pmin` shell mutants (mutation_gate.md round 3b).
"""
import pytest

import g2_common as G

if G.ops is None:
    pytest.skip(f"opensees.pyd not found/loadable in {G.DIST_BIN}", allow_module_level=True)

import floor_sensitivity as FS  # noqa: E402

pytestmark = [pytest.mark.zone_b, pytest.mark.no_gmsh, pytest.mark.slow]


@pytest.fixture(scope="module")
def runs():
    return {1.0: FS.run_strip(1.0), 0.5: FS.run_strip(0.5)}


def test_limit_load_moves_less_than_two_percent_when_pmin_is_halved(runs):
    a, b = runs[1.0], runs[0.5]
    print("\n" + FS.report(runs))
    assert not a["aborted"] and not b["aborted"], "a run was aborted (a step failed after 6 halvings)"
    assert a["reached"] >= 0.999 * FS.DELTA and b["reached"] >= 0.999 * FS.DELTA, (a["reached"], b["reached"])
    c = FS.compare(a, b)
    assert c["final"] < 0.02 and c["peak"] < 0.02, f"limit load moves {100 * c['final']:.2f} % (final) / {100 * c['peak']:.2f} % (peak) when p_min is halved"


def test_the_floor_is_engaged_on_the_model_and_its_energy_source_vanishes_with_pmin(runs):
    a, b = runs[1.0], runs[0.5]
    assert a["n_at_floor"] >= 4 and a["n_events"] >= 1 and a["sum_Wf"] > 0.0, (a["n_at_floor"], a["n_events"], a["sum_Wf"])
    assert b["sum_Wf"] < a["sum_Wf"], f"W_f does not vanish with p_min: {b['sum_Wf']:.3e} (F/2) vs {a['sum_Wf']:.3e} (F)"
    assert a["R_final"] > 0.0 and b["R_final"] > 0.0
