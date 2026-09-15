# ADR-94 addendum (WP F8) — the ASD Drucker-Prager ASSOCIATED-flow wall on the Prandtl deck

**Status:** RESULTS NOTE, 2026-09-14. Branch `wp/100-asd-dp-associated-wall`, PR
[#836](https://github.com/nmorabowen/OpenSees/pull/836). Investigation build
`9c2f964ea` (the `ladruno` tip; `dist/bin` staged at
`.claude/worktrees/release-build-9c2f964`). Predecessor: wp/94f
([#832](https://github.com/nmorabowen/OpenSees/pull/832),
`Ladruno_implementation/_adr94f_results.md`), note 95 §3b.

## 0. The answer

**wp/94f's layer (a) is correct in one direction only, and the other direction is
the associated leg.** It added ADR-97's elastic-metric apex classification to
`Backward_Euler` and **unioned** it with `DruckerPrager_YF::check_apex_region`:

```cpp
bool be_in_apex = yf.check_apex_region(...);            // (p - p_apex) >= eta*q
if (!be_in_apex)
    be_in_apex = cp_apex_region(...);                   // (p - p_apex) >= (K*etabar/G)*q
```

A union takes the **wider** of the two regions. Which one is wider depends on the
flow rule, because the two slopes are `eta` and `K*etabar/G`:

| leg | Euclidean slope `eta` | exact slope `K*etabar/G` | union keeps | verdict |
|---|---|---|---|---|
| psi = 0 (non-associated, the ADR-95 deck) | 0.4457 | **0** (`p >= p_apex`) | the exact one | correct — this is what wp/94f measured |
| psi = phi (**associated**) | 0.4457 | **4.3089** | the **Euclidean** one | **~10x too wide** |

On the associated leg the union therefore apex-projects every trial in the wedge

    eta*q  <=  p - p_apex  <  (K*etabar/G)*q

whose correct return is to the cone **flank**. The committed stress is pinned at
`sigma_apex = p_apex*I` with **no deviator**, and under `tangent_type Continuum`
(what the ADR-95 ASD decks use) that Gauss point additionally reports a **zero
tangent**. **No refusal is issued** — the material reports success — so the
failure is silent, which is why it presents as "the leg walls while still
hardening" rather than as a diagnosable refusal count.

Measured at a single Gauss point on `9c2f964ea`
(`tests/test_f8_asd_dp_associated_apex.py`), associated, the ADR-95 cone
(`eta = 0.445749`, `xi_c = 0.115470`, `p_apex = 0.259047`, `K/G = 9.6667`), trial
at `(p - p_apex)/q = 2.0`:

| | p | sqrt(J2) |
|---|---|---|
| committed by the tip build | **0.2590471** (the apex) | **0.0** |
| closed-form cone return | **-0.5314874** | **0.3523800** |

The closed form is elementary and is the oracle the test pins to 1e-6 relative:
`dgamma = f_tr/(G + K*eta*etabar)`, `q_ret = q_tr - G*dgamma`,
`p_ret = p_tr - K*etabar*dgamma`, valid exactly while `q_ret >= 0` — which *is*
the elastic-metric apex test. ADR-97's own oracle
(`Ladruno_implementation/adr97_oracle/cppm_dp.py`, lines 108-129) already carried
both tests side by side and prints a MISMATCH when they disagree; nothing new had
to be derived.

## 1. The fix

For yield functions declaring `yf_apex_elastic_metric` (Drucker-Prager only) the
elastic-metric test **replaces** the Euclidean one in `Backward_Euler` instead of
widening it — which is what the trait means, and what `Closest_Point` has always
done. At `etabar = 0` the two regions coincide (the elastic-metric answer is
`p > p_apex`, which strictly contains the Euclidean cone), so **wp/94f's result is
untouched by construction, not by luck**. Every other yield function still calls
`check_apex_region` and is unchanged. The flank-first fallback (wp/94f layer (b))
is left exactly as it was and remains the safety net behind the classification.

One file, one branch: `SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h`
(plus the stale union claim corrected in `YieldFunctions/DruckerPrager_YF.h`).
Marked `// Ladruno (ADR-94 addendum, F8)`.

## 2. The deck measurement — the ADR-95 R3 gate's own associated control

The ADR-95 collapse gate (`tests/test_r3_prandtl_collapse_gate.py`) is now
material-pluggable: `_define_material()` carries the UW-to-ASD cone/flow mapping,
`_run_leg(material="UW", gp_probe=None)` keeps every existing gate leg building
exactly what it built before, and
`Ladruno_files/testbed/hypo_bearing/r3_assoc_probe.py` drives one leg with a
per-step Gauss-point census in the `(p, sqrt(J2))` half-plane. The associated leg
is run at `h0 = 1.0` (`CONTROL_H0 = 0.5` is the gate's own rung; `h0 = 1.0` is the
cheapest and is where the TIMs act reported the wall).

<!-- F8-RESULTS-TABLE -->

## 3. What is NOT claimed

<!-- F8-CAVEATS -->
