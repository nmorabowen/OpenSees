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

A union takes the **wider** of the two regions, and the exact slope `K*etabar/G`
overtakes the Euclidean `eta` as soon as

    etabar > eta * G / K

On this deck `G/K = 0.10345` and `eta = 0.4457`, so the crossover is at
`etabar = 0.0461` — **psi ~ 2.3 deg**. The union is therefore unsafe not "under
associated flow" but under **essentially any dilatancy**; `etabar = 0` is the one
case where it is safe, and it is the one case wp/94f measured.

| leg | Euclidean slope `eta` | exact slope `K*etabar/G` | union keeps | verdict |
|---|---|---|---|---|
| psi = 0 (the ADR-95 deck wp/94f measured) | 0.4457 | **0** (`p >= p_apex`) | the exact one | correct |
| psi ~ 2.3 deg (`etabar = eta*G/K`) | 0.4457 | 0.4457 | either | the crossover |
| **psi ~ phi/2** (`etabar = eta/2`) | 0.4457 | **2.1545** | the **Euclidean** one | **4.8x too wide** |
| psi = phi (**associated**) | 0.4457 | **4.3089** | the **Euclidean** one | **~10x too wide** |

Above the crossover the union therefore apex-projects every trial in the wedge

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

All four legs below are the **same mesh, same deck, same session**
(200 `LadrunoBrick -formulation bbar`, 1386 DOF, `system Pardiso`,
`SUBDIV_BUDGET = 80`, `WALL_BUDGET = 3600 s`, push to s/B = 0.15). Exact
`q_u = q0*N_q = 138.907 kPa`.

| leg | build | q_max (kPa) | ratio | mode | tail % | ds/floor | CAPACITY | failed / subdiv | **relaxed steps** | wall s |
|---|---|---|---|---|---|---|---|---|---|---|
| UW associated (reference) | `9c2f964ea` | 268.75 | **1.9348** | TARGET | 0.120 | 4400 | **yes** | 0 / 0 | **0 / 329** | 63 |
| **ASD associated, PRE-fix** | `9c2f964ea` | 226.51 | 1.6307 | **BUDGET** | 8.270 | 50 | **NO** | 898 / 81 | **256 / 560** | 894 |
| **ASD associated, POST-fix** | `3324485f7` | 268.38 | **1.9321** | BUDGET | 0.163 | 2500 | **yes** | 1528 / 81 | **638 / 687 (92.9 %)** | 1122 |
| ASD psi = 0, PRE-fix | `9c2f964ea` | 150.71 | 1.0850 | TARGET | 0.001 | 5000 | yes | 0 / 0 | **0 / 329** | 72 |
| ASD psi = 0, POST-fix | `3324485f7` | 150.71 | 1.0850 | TARGET | 0.001 | 5000 | yes | 0 / 0 | **0 / 329** | 107 |
| **UW psi = 0** (the gate's own leg, refactor + fix regression) | `3324485f7` | 150.71 | **1.0850** | TARGET | 0.001 | 4400 | yes | 0 / 0 | 41 |

The last row is the check that the material factoring and the `gp_probe` hook
change nothing the gate builds, and that the C++ change does not reach the
vanilla `DruckerPrager`: the gate's own `h1.0_nonassoc` record is **1.0849**
(module docstring) / 1.0850 on the repaired material (ADR-95 §4).

`BUDGET` is the gate's own "capacity WITH A NAMED ALLOWANCE" — the mode all three
of the gate's own legs are RECORDED in (its module table: BUDGET at ds/floor
2500 / 1250 / 800), though in this session, on this box, both reference legs ran
on to `TARGET` at ds/floor 4400. The load had been flat for the last tenth of the
run with the step still 2500x the floor. **The post-fix ASD associated leg is a
CAPACITY at 1.9321 against the UW reference's 1.9348 — 0.14 % apart.**

**The `relaxed steps` column is not decoration.** 638 of the post-fix ASD
associated leg's 687 converged steps (92.9 %) needed the ladder's THIRD rung —
`KrylovNewton` at 10x the `NormUnbalance` tolerance, 60 iterations — where the UW
associated leg and both psi = 0 legs needed **zero**. The gate records `nrelax`
precisely so this cannot pass unnoticed, and it is the honest qualifier on
"agree to 0.075 %": the ASD leg reaches the same answer, on a looser
tolerance, most of the way. Part of the asymmetry is deck, not material: the ASD
decks run `strict_convergence 1` and `n_max_iterations 100` (inherited from
`asd_path_diag.py`, the settings ADR-95 measured on) and the vanilla
`DruckerPrager` has no equivalent switch, so a step the ASD material refuses is a
step the UW material would have silently accepted. That is also where the 1528
failed ladder attempts come from.

### The two implementations now follow the same PATH, not merely the same peak

`q` (kPa) at matched settlement, both legs of this session:

| s/B | UW assoc | ASD assoc POST-fix | rel | ASD assoc PRE-fix |
|---|---|---|---|---|
| 0.0050 | 126.68 | 126.68 | +0.000 % | 130.08 |
| 0.0169 | 226.19 | 226.26 | +0.027 % | *(walled at 0.01685)* |
| 0.0200 | 233.04 | 233.14 | +0.043 % | — |
| 0.0400 | 255.97 | 256.08 | +0.045 % | — |
| 0.0600 | 261.63 | 261.83 | +0.075 % | — |
| 0.0800 | 264.18 | 264.34 | +0.062 % | — |
| 0.1000 | 265.93 | 266.08 | +0.056 % | — |
| 0.1366 | 268.24 | 268.38 | +0.052 % | — |

Worst over the whole common range: **0.075 %**.

**How large the stress error was, and why that is not the point.** Over its own
(short) range the PRE-fix leg's load-settlement curve was never far from UW's:
worst 5.17 % at s/B 0.00202, and 0.141 % at its terminal point. The
misclassification's stress error at a handful of Gauss points is small in
absolute terms (the apex sits at 0.259 kPa in a 200 kPa field). What it destroys
is the **iteration**: an apex projection is a different, non-smooth map with a
zero `Continuum` tangent, so the outer Newton loses its quadratic convergence at
exactly the Gauss points the mechanism is forming around. The symptom is a
controller death, not a wrong number — which is why "the load path looks fine"
is not evidence that the return map is.

**The pre-fix failure is SILENT.** `grep -c "rejecting step"` over the whole
894-second pre-fix associated log is **0**: no refusal, no NaN, no apex message.
Every Gauss point reported success; 12 → 32 of them were simply pinned at
`p = p_apex = 0.259047` with `sqrt(J2) = 0` while the surrounding field ran to
−187 kPa. What failed was the outer Newton — 898 failed attempts, the step
ground down to 50x the floor over the final tenth of the run, and the load was
still climbing at 8.3 % of the initial tangent when the subdivision budget went.
That is why this reads as "the element walls" rather than as a material defect.

**The psi = 0 leg is BYTE-IDENTICAL across the fix**: 329 rows, all of
`(s_m, s_over_B, q_kPa, ds_mm, relaxed)` equal, `q_max` 150.709859 both sides.
That is the construction argument made experimentally — at `etabar = 0` the
elastic-metric region contains the Euclidean one, so replacing the union by the
exact test cannot change anything.

### Gate battery, post-fix build `3324485f7`

| file | result |
|---|---|
| `tests/test_f8_asd_dp_associated_apex.py` (new) | **4/4** — and **1 failed / 3 passed on `9c2f964ea`**, so it gates the fix |
| `tests/test_adr94f_asd_apex_fallback.py` | **4/4** — wp/94f's own cases unchanged |
| `tests/test_adr97_p4_inertness.py` | **10/10** — `Backward_Euler` still BYTE-IDENTICAL on all 23 baseline decks |
| `tests/test_adr94c_numerics.py` + `test_adr94_redblue_numerics.py` | **11/11** |
| total | **29 passed** (pre-fix control on `9c2f964ea`: 25 passed, the new file excluded) |
| `tests/test_r3_prandtl_asd_associated.py` (new, slow tier) | see Sec. 4 |

## 3b. Review round 1 — what adversarial review found and what changed

Verdict MERGE-OK conditional. The classification fix itself held: 43 sweep rows
across psi and nu (0.2 / 0.49) with zero misclassifications, byte-identical for
every yield function without the trait, 232 other ASDPlastic tests green, and
1.9321 confirmed a genuine plateau. Four SHOULD-FIXes and five nits, all applied.

**(1) The condition is not "associated".** Restated everywhere as
`etabar > eta*G/K` — psi ~ 2.3 deg on this deck — with the reviewer's
`etabar = eta/2` row (exact slope 2.1545 vs Euclidean 0.4457, a 4.8x wedge)
added as two `zone_a` gate cases: the flank row against the closed form, and the
apex row at ratio 2.1760. Both passed on the round-0 build already; they are
there because the *framing* was wrong, not the fix.

**(3) A real outcome regression, found and fixed.** Narrowing the apex region
routes near-boundary trials into the flank scalar Newton, whose `dPhi/dlambda`
carries the pinned vanilla `df/dk = -1` cohesion term that `f` does not contain.
With cohesion SOFTENING (`ScalarLinearHardeningParameter = -20000`) that Newton
**converges** — `rc = 0`, `|f| ~ 1e-7`, no exhaustion — onto a state whose
deviator points OPPOSITE the trial deviator. Measured, `etabar = eta`, one Gauss
point, identical at `strict_convergence` 0 and 1:

| (p - p_apex)/q | HS = 0 (correct) | HS = -20000, pre-guard |
|---|---|---|
| 3.0 | p -0.189103, q 0.199762, `s_zz-s_xx` **+0.346** | p -0.478023, q 0.328548, `s_zz-s_xx` **-0.569** |
| 4.3089 | the apex, q 1.0e-15 | p -0.838550, q 0.489254, `s_zz-s_xx` **-0.847** |

Pre-fix those trials were apex-projected (q = 0), so the round-0 PR made them
worse. Layer (b) could not catch it: its other trigger `be_exhausted` is
computed **only** under `strict_convergence` (off by default), and this state
does not exhaust — it converges to the wrong root.

Fixed with a geometric guard after the flank loop, `yf_apex_elastic_metric`
scope: a Drucker-Prager return is a non-negative radial scaling of the trial
deviator plus a pressure change, so `dot(dev_ret, dev_tr) < 0` is inadmissible
for any parameters. It is reported through the existing `be_flank_failed`
channel, so it inherits layer (b)'s apex fallback and the existing fail-loud
refusal without adding an exit.

**The guard needed a size test, and finding that cost the round-1 build.** The
raw sign test regressed wp/94f's own zero-dilatancy acceptance path: that path
walks up the cone in equal steps (q = 0.0924, 0.0693, 0.0462, 0.0231, ...) and
lands EXACTLY on the vertex on its 5th, where `dev_ret` is zero to round-off and
the sign of the dot product is a coin flip — the leg refused at step 5 with
`codes = [0,0,0,0,-3]`. The guard now also requires `||dev_ret|| > tol_yf`, the
stress-unit scale the rest of the integrator measures in (ADR-94 M5). In the
softening reproducer `||dev_ret|| ~ 0.46`; in the round-off case ~1e-17.

**(6) Stale comments.** `DruckerPrager_YF.h`'s wp/97b note ("used only by
Backward_Euler ... two integrators, two answers") is corrected: for Drucker-Prager
`Backward_Euler` no longer calls `check_apex_region` at all, so the member is
dead code for that YF. The residual asymmetry is stated in its place —
`Closest_Point` applies the elastic-metric test to EVERY `yf_has_apex` YF while
`Backward_Euler` applies it only to the opted-in ones, so MohrCoulomb /
HoekBrown / TensionCutoff still get two answers — and the vanilla-ledger row's
"still used by every non-opted-in path" now says "except DP".

**(8) Measurement honesty — the `relaxed` column.** Added to the table in Sec. 2,
to the slow gate's module docstring and to the PR: the ASD associated leg needed
the ladder's third rung on **638 of 687** steps against **0** for UW and both
psi = 0 legs, and the ASD decks carry `strict_convergence 1` /
`n_max_iterations 100` with no UW equivalent.

**Nits.** (4) `apex_stress()` ignores the back stress and `cp_apex_region` tests
`dev(sigma)` rather than `s - alpha` — recorded in `LEDGER_quirks` as
pinned-not-fixed (identical before and after F8; measured `rc = -3` under strict,
`q = 0.393` inadmissible under non-strict with `alpha0` nonzero), and the
`be_apex_project` comment claiming every `yf_has_apex` YF is perfectly plastic is
corrected — Drucker-Prager opts in WITH hardening template parameters. (5)
`adr97_oracle/baselines/dump_hist.py` now records `ladrunoBuild()` in each dump
(recorded, not asserted: a baseline is a different build by construction; what
was missing was provenance). (9) the "all three end BUDGET" sentence is softened
to say that is the gate's RECORDED mode, while both reference legs here ran on to
TARGET. (10) `n_eucl_wedge` is documented as diagnostic-only and unable to see
the defect: it is evaluated on COMMITTED stresses, so the only thing that can land
in it is apex round-off — UW 125-137 vs ASD 0 is round-off, not mechanics. (11)
the cheap gate's docstring said ~15 s; measured 0.71 s for 12 tests.

**Not changed, and why.** `tests/test_adr94_matrix.py` regenerates the tracked
`_adr94_matrix.md` on every run (a wart already recorded in the ADR-94
implementation log) and its regeneration also strips the file's provenance
preamble, so the regenerated file was REVERTED rather than committed. Its cell
diff is worth reading though: the only cells that moved are the two
`Numerical_Algorithmic_*` columns of VonMises / MohrCoulomb / HoekBrown, i.e.
ADR-97 P4's re-point (#829) drifting against a table last regenerated at
`3622d6214` — **no Drucker-Prager cell moved**, which is independent evidence
that this WP's change is confined to Drucker-Prager.

### Gate battery, round-1 build

| file | result |
|---|---|
| `tests/test_f8_asd_dp_associated_apex.py` | **12/12** in 0.71 s (was 4; +2 for the `eta/2` reframing, +6 for the softening guard and its non-softening controls) |
| `tests/test_adr94f_asd_apex_fallback.py` | **4/4** |
| `tests/test_adr97_p4_inertness.py` | **10/10** — still byte-identical on all 23 baseline decks |
| `tests/test_adr94c_numerics.py` + `test_adr94_redblue_numerics.py` | **11/11** |
| the five files together | **37 passed in 10.16 s** |
| `pytest -k "adr84 or adr94 or adr95 or adr97 or asdplastic or f8"` | **273 passed, 9 skipped** (+ the pre-existing `test_adr94_matrix.py` cwd-relative-path error, which passes when run from `tests/`) |

## 3c. Review round 2 — MERGE-OK, two follow-ups closed

Re-verification returned **MERGE-OK**: the gate cases reproduce the reviewer's
table exactly, the flip guard fixes the reported softening defect, the vertex
floor was verified at 14 resolutions and lands at the predicted `||r_ret||`
crossover, inertness 10/10, 205 other ASDPlastic tests green. Two follow-ups:

**(1) Two C++ comment blocks still said "associated"** — `ASDPlasticMaterial3D.h`
enumerated only `etabar = 0` and `etabar = eta`, and `DruckerPrager_YF.h` said
"Under ASSOCIATED flow it is the narrower one". Both now state the crossover
`etabar > eta*G/K` (psi ~ 2.3 deg on this deck) and carry the `etabar = eta/2`
row. Comment-only.

**(2) The apex region is in the RELATIVE deviator — fixed, and the attribution
is one layer up from where the review put it.** Drucker-Prager's surface is
written in `r = dev(sigma) - alpha`, so both the flip guard AND the apex
classification have to be. The review named the guard; the measured refusal is
actually the apex **classification** reaching the state first — the message is
`predictor-classified apex: ... |f(sigma_apex)| = 0.0346`, which is exactly
`sqrt(J2(alpha))`. `cp_apex_region` tested the flip on `dev(sigma)`, so it called
a cone state APEX; `be_apex_project` was then asked for a vertex
`DruckerPrager_YF::apex_stress()` cannot supply (it ignores `alpha`), its own
`|f(sigma_apex)| <= tol_yf` guard fired, and the step was refused. Both were
fixed in the same variable; the guard's premise was wrong too and would have been
the next thing to bite. Reviewer's reproducer,
`alpha = (0.02, 0.02, -0.04, 0, 0, 0)`, `Ht = 0`, associated:

| trial `(q_tr, p_tr)` | exact return (closed form) | round 1 | round 2 |
|---|---|---|---|
| (0.05, 0.5103) | `sqrt(J2(r))` 0.0173156018, `p` 0.2202010508 | **rc -3** (strict 0 and 1) | **rc 0**, matches to 1e-6 |
| (0.02, 0.3810) | `sqrt(J2(r))` 0.0173206061, `p` 0.2201898240 | **rc -3** (strict 0 and 1) | **rc 0**, matches to 1e-6 |

The two trials do **not** share a return point — they differ by 5.0e-6 in
`sqrt(J2(r))`. The review quoted one pair (0.017321, 0.220190); that is trial 2's,
and the test now computes each trial's own closed form in-test rather than
transcribing a constant, with a separate provenance case asserting that the
in-test oracle reproduces the reviewer's quoted figure to 5e-7.

**This refusal is older than F8 and is not the union's doing.** The Euclidean
`check_apex_region` over-classifies the same trial (`p - p_apex = 0.1220` against
`eta*q_rel = 0.0244`), so every arrangement since wp/94c made the apex projection
live — Euclidean alone, Euclidean OR elastic-metric, or elastic-metric alone —
reaches it. What round 2 changed is the test the integrator now relies on.

**What is still pinned:** `apex_stress()` ignores `alpha`. The vertex of this
surface is `alpha + p_apex*I`, not `p_apex*I`, so a state that genuinely IS in
the apex region with a nonzero back stress is still refused rather than
projected. Refusing is the safe half of wrong, and the fix is a one-line change
to a vanilla yield function that also moves `Closest_Point`'s apex return — its
own WP, its own gate.

In the guard, both `r` are formed against the back stress current at their own
state (trial internal variables for the returned stress, committed ones for the
elastic predictor) — the pairing the yield function itself uses; with `Ht = 0`
they coincide. The back stress is reached through
`YieldFunctionType::internal_variables_t`'s first element, which for every yield
function declaring `yf_apex_elastic_metric` IS the back stress, so no new
accessor was needed, and `cp_apex_region`'s correction is behind the same
`if constexpr` so every other family is byte-identical there. The
`||r_ret|| > tol_yf` floor is unchanged, and the guard comment now says plainly
that a flip smaller than `f_absolute_tol` is invisible **by construction** — that
is the integrator's own `f` scale, not a tuning knob.

Regression checks re-measured on the round-2 build: wp/94f's psi = 0 walk
`codes = [0]*10` ending at the apex; the `HS = -20000` rows now commit the apex
(`q = 0`) instead of a flipped deviator, at both strict settings; the `HS = 0`
and `HS = +2000` rows reproduce their pinned values to the last printed digit.

| gate | round 1 | round 2 |
|---|---|---|
| `tests/test_f8_asd_dp_associated_apex.py` | 12/12 | **17/17** |
| the five ASDPlastic files | 37 | **42 passed in 6.43 s** |
| `-k "adr84 or adr94 or adr95 or adr97 or asdplastic or f8"` | 273 passed, 9 skipped | **278 passed, 9 skipped** |

The `Closest_Point` gates (ADR-97 P1/P2/P3/P5/P6) are inside that sweep and stay
green: `cp_apex_region` is shared with them, and with a zero back stress the
correction is arithmetically a no-op.

## 4. Artifacts

All under `Ladruno_files/testbed/hypo_bearing/`, prefix `f8_`:
`f8_r3_h1.0_assoc_{uw,asd_PRE,asd_POST}.csv` (the load-settlement curves),
`f8_r3_h1.0_nonassoc_asd_{PRE,POST}.csv` (the byte-identical pair),
`f8_census_assoc_{uw,asd_PRE,asd_POST}.csv` (the Gauss-point census),
`f8_asd_assoc_h1.0_{PRE,POST}.log` (trimmed, with the refusal/NaN counts of the
full logs in their headers). Driver: `r3_assoc_probe.py`.

## 3. What is NOT claimed

* **The `h0 = 0.5` rung (`CONTROL_H0`) was not re-run on ASD.** The gate's own
  associated control lives at `h0 = 0.5`; this WP measured `h0 = 1.0`, the
  cheapest rung and the one where the wall was reported. The ASD-vs-UW agreement
  is therefore established at one resolution, not across the sequence.
* **The associated collapse load itself is not a validated capacity in the
  physical sense.** `psi = phi` on a bounded mesh is the strong upper solution
  and the gate has always treated the associated leg as a falsification control,
  not as an answer to compare against Prandtl. What is claimed here is that the
  two IMPLEMENTATIONS of one cone now agree on it.
* **The zero apex tangent is untouched.** wp/94c chose `Stiffness = 0` at the
  apex under `tangent_type Continuum` deliberately ("the honest continuum
  operator at a perfectly plastic apex is ZERO"), and `Secant` — the default —
  blends it with the elastic operator. This WP only changed WHICH states are
  classified as apex, not what happens to one that is. On the post-fix leg **412
  of 1600 Gauss points are apex-pinned at s/B 0.1225** and the leg still advances
  to a plateau (the UW reference carries 436 at s/B 0.148), so the zero tangent
  is not by itself a wall — that was the alternative hypothesis and it is not
  supported.
* **The associated leg is still the expensive one.** Post-fix it takes 1122 s and
  1528 failed ladder attempts against the UW material's 63 s and zero, and it
  ends on `BUDGET` rather than `TARGET`. What changed is that it advances: tail
  0.163 % vs 8.270 %, terminal step 2500x the floor vs 50x, s/B 0.1366 vs 0.0169.
  Why the ASD path costs more Newton work than UW's on the same cone is not
  answered here.
* **Nothing is claimed about `Closest_Point`.** It has always classified in the
  elastic metric; this WP makes `Backward_Euler` agree with it, which is what
  wp/94f said it was doing.
* The number the task brief quoted for the campaign's associated UW leg (1.60)
  is the **`h0 = 0.5`** measurement recorded in note 95 §4 (1.6026). At
  `h0 = 1.0`, measured here for the first time, UW associated reads **1.9348**.
