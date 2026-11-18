"""WP-232k B3 (the AMEND round): S1 one-point legs for the admitted OPTIONS, against 232k-V's
tests/data/mpm/soil/s1_onepoint_options.json (36 legs, d935310b0; each leg carries its `options`,
appended to the deck). Rules (wp-briefs.md "232k", rulings on B3's option legs, 2026-11-18):

- SANISAND IntScheme 129, HS -implex, -smallStrain (HSSmall cyclic shear), the HS deck's tension
  leg: the S1 rule (rel 1e-8 before the FE's own Gauss-point spread passes 1e-8, 1e-3 after).
- SANISAND IntScheme 2: triax and oedo the S1 rule; shear: 1e-8 before the switch step, after it
  every MPM stress inside the FE's Gauss-point envelope widened by 1e-3 relative. The record holds
  GP 1's stress and the GP spread (max - min over the 8 GPs), so the envelope used is
  [sig_GP1 - spread, sig_GP1 + spread] (it contains [min, max]); the state is held to the GP rel
  spread + 1e-3.
- HS -sigmat: the `tension` path is the gated leg (S1 rule); its other paths never go tensile
  and are checked by the S1 rule too.
- SANISAND -implex: with every DOF prescribed it commits the scheme-1 state, so its one-point leg
  is an EQUALITY: at particles, -implex equals the deck (scheme 1) bit for bit, step by step. The
  option's real gate is the boundary-value leg (the K0p column, mode A gimp, against
  k0p_column_capped_implex.json, sigma RMS <= 0.02 rho g l0).
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

sys.path.insert(0, os.path.dirname(__file__))
import ladruno_mpm_soil_s_harness as S  # noqa: E402

pytestmark = pytest.mark.zone_a
DATA = Path(__file__).resolve().parent / "data" / "mpm" / "soil"
REC = DATA / "s1_onepoint_options.json"


def _legs():
    if not REC.exists():
        return []
    return sorted(json.loads(REC.read_text(encoding="utf-8"))["legs"])


def _ops():
    try:
        from _testbed import ops
    except Exception as exc:  # pragma: no cover - environment
        pytest.skip(f"opensees module not importable: {exc}")
    if not hasattr(ops, "mpm"):
        pytest.skip("this opensees build has no `mpm` command")
    return ops


LEGS = _legs()
IMPLEX = [lg for lg in LEGS if lg.startswith("sanisand_implex_")]
PLAIN = [lg for lg in LEGS if lg not in IMPLEX and not (lg.startswith("sanisand_is2_") and
                                                     lg.endswith("_shear"))]
IS2_SHEAR = [lg for lg in LEGS if lg.startswith("sanisand_is2_") and lg.endswith("_shear")]


def _rec():
    return json.loads(REC.read_text(encoding="utf-8"))


@pytest.mark.skipif(not LEGS, reason="s1_onepoint_options.json not committed")
@pytest.mark.parametrize("mode", ["A", "F"])
@pytest.mark.parametrize("leg", PLAIN or ["none"])
def test_option_leg_s1_rule(leg, mode):
    out = S.s1_run(_ops(), leg, mode, rec=_rec())
    bad = [r for r in out["steps"] if not (r["err_sig"] <= r["tol"] and r["err_state"] <= r["tol"])]
    assert not bad, (leg, mode, bad[0], out["max_err"])


def _envelope_run(ops, leg, mode, rec):
    """S1 with the scheme-2 shear rule: per step, after the switch, sigma inside
    [sig_GP1 - spread, sig_GP1 + spread] widened by 1e-3 relative."""
    L = rec["legs"][leg]
    sw = (L.get("fe_self_consistency") or {}).get("first_step_gp_rel_spread_gt_1e8")
    dt = 1.0 if mode == "A" else 1e-4
    cmds, nd = S.s1_commands(rec, leg, mode, dt)
    for name, args in cmds:
        getattr(ops, name)(*args)
    flip = int(L["flip_before_step"])
    bad = []
    for k in range(len(L["dF"])):
        if k + 1 == flip:
            ops.mpm("stage", 1, 1)
        ops.mpm("gridvel", "all", "setAffine", *S.s1_L(rec, leg, k, nd, dt))
        ops.mpm("analyze", 1, dt)
        sig = S._l(ops.mpm("particleData", "sigHat"))
        szz = S._l(ops.mpm("particleData", "sigZZ"))
        st = S._l(ops.mpm("particleResponse", "state"))
        N = len(szz)
        ref = L["steps"][k]
        rs = ref["sig"]
        scale = max(abs(v) for v in rs) or 1.0
        sscale = max(abs(v) for v in ref["state"]) or 1.0
        ns = len(st) // N
        after = sw is not None and k + 1 >= sw
        for p in range(N):
            s6 = S._sig6(sig, szz, nd, p)
            if after:
                w = float(ref.get("gp_spread", 0.0)) + 1e-3 * scale
                ok = all(abs(a - b) <= w for a, b in zip(s6, rs))
                ok = ok and max(abs(st[ns * p + i] - ref["state"][i]) for i in range(ns)) <= \
                    (float(ref.get("gp_rel_spread", 0.0)) + 1e-3) * sscale
            else:
                ok = max(abs(a - b) for a, b in zip(s6, rs)) <= 1e-8 * scale and \
                    max(abs(st[ns * p + i] - ref["state"][i]) for i in range(ns)) <= 1e-8 * sscale
            if not ok:
                bad.append((k + 1, p))
                break
    return bad


@pytest.mark.skipif(not LEGS, reason="s1_onepoint_options.json not committed")
@pytest.mark.parametrize("mode", ["A", "F"])
@pytest.mark.parametrize("leg", IS2_SHEAR or ["none"])
def test_option_is2_shear_envelope(leg, mode):
    bad = _envelope_run(_ops(), leg, mode, _rec())
    assert not bad, (leg, mode, bad[:3])


@pytest.mark.skipif(not LEGS, reason="s1_onepoint_options.json not committed")
@pytest.mark.parametrize("mode", ["A", "F"])
@pytest.mark.parametrize("leg", IMPLEX or ["none"])
def test_option_implex_equals_scheme1(leg, mode):
    """-implex commits the scheme-1 state with every DOF prescribed: bit for bit at particles."""
    ops = _ops()
    rec = _rec()
    base = S.load("s1_onepoint.json")
    deck_leg = leg.replace("sanisand_implex_", "sanisand_")
    L, Lb = rec["legs"][leg], base["legs"][deck_leg]
    dt = 1.0 if mode == "A" else 1e-4
    traces = []
    for r, lg in ((rec, leg), (base, deck_leg)):
        cmds, nd = S.s1_commands(r, lg, mode, dt)
        for name, args in cmds:
            getattr(ops, name)(*args)
        tr = []
        for k in range(len(r["legs"][lg]["dF"])):
            if k + 1 == int(r["legs"][lg]["flip_before_step"]):
                ops.mpm("stage", 1, 1)
            ops.mpm("gridvel", "all", "setAffine", *S.s1_L(r, lg, k, nd, dt))
            ops.mpm("analyze", 1, dt)
            tr.append((S._l(ops.mpm("particleData", "sigHat")),
                       S._l(ops.mpm("particleResponse", "state"))[:26]))
        traces.append(tr)
    assert len(L["dF"]) == len(Lb["dF"])
    first = next((k for k, (a, b) in enumerate(zip(*traces)) if a[0] != b[0]), None)
    assert first is None, f"{leg} {mode}: -implex departs from scheme 1 at step {first + 1}"


@pytest.mark.skipif(not (DATA / "k0p_column_capped_implex.json").exists(),
                    reason="k0p_column_capped_implex.json (232k-V) not committed")
def test_option_implex_boundary_value_leg():
    """-implex's gate: the K0p column (staging only), mode A gimp, against the FE with -implex."""
    import ladruno_mpm_soil_harness as K
    r = K.k0p_implex_bvp(_ops())
    assert r.get("k0p", {}).get("valid", True)
    assert r["gate"]["passed"], (r.get("meter"), r.get("reason"))
