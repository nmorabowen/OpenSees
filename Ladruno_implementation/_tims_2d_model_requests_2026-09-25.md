# TIMs 2D model — consolidated requests to the fork, 2026-09-25 (F18–F23)

Intake file, written by the TIMs Workbench `2d-model` act (`C:\Users\nmb\Dropbox\UANDES EC\TIMs Workbench`, branch `work/ape/2d-model`, `Tries/2d-model/`) and committed here the way `_adr92_tims_request_2026-09-05.md` was. It supersedes every open item of the act's earlier prompts (the 2026-09-15, -15b, -15c, -18 and -18 Bézier requests, which were pasted, not committed). Attachments sit beside it in `_tims_2d_model_requests_2026-09-25/`.

**How to work it.** Read `CLAUDE.md` first and hold it. One `wp/<n>-<slug>` branch and a draft PR per item, on day one. Build with `Ladruno_scripts\build.bat`. Run the affected pytest suites before and after and paste the output. `// Ladruno` on every touched upstream line. Ledgers, guide and banner in the same PR. **Never run the Workbench's models and never edit anything in the Workbench**: reproduce on your own decks and on the attached material states. Verify every line number and every number cited here against this checkout before relying on it — the act read them on `79e062367`. End each PR with what was verified, what was not, and which ledger rows were added. Answer in a `*_tims_report.md` as before.

---

## 0. Where the earlier requests stand (the act's reading)

| item | subject | fork answer | act's status |
|---|---|---|---|
| F10 | IMPL-EX recipe on the strip | ADR-92 F10 verdict (26d5c607f) | read; recipe does not transfer to this deck (bare `-implex` aborts at s/B 0.0004 on the edge Gauss point) |
| F10b | a point `ModifiedEuler` cannot integrate under small increments | — | **open → folded into F18 and F21 below** |
| F11 | threaded state-determination loop | WP-107, #843 | read; SANISAND refused from the threaded loop → **continued as F19** |
| F12 | `IntScheme 2` qualification | WP-105, #844 | read and accepted (qualified per increment, refuted as the BVP integrator) → **F18(c) builds on it** |
| F13 | `-pRe` elastic-only floor | WP-106, #842 | read; on the strip it does not move the wall (0.0294 against 0.0288) |
| F14 | `-flipAlphaIn init` default | WP-112, #849 | closed |
| F15 | `GetElastoPlasticTangent` defects | WP-110, #847 | closed; the tangent families still part (§1.5) |
| F16 | `"tangent"` response | WP-110 wired `"tangent"` | first half closed; the `"tangentEP"` half **open → F20(c)** |
| F17 | BezierTri6 -bbar, non-associated | WP-114, #848 | closed; the fixed element tracks the quad within 1 % |

## 1. What the act measured since (context for every ask)

**The deck.** Plane-strain strip footing, `B = 1.5` m, `15B × 12B`, self-weight, 7.65 kPa surcharge outside the footprint, 18.4 kN/m footing weight, Jaky K0 through `nu*`, rigid footing by `LadrunoKinematicCoupling`, displacement push. `LadrunoQuad -bbar`, B/8 uniform: 2 430 elements, 9 720 Gauss points, 4 860 DOF. `system Pardiso`, `NormUnbalance` at 1e-5 of the reference load, Newton → NewtonLineSearch → KrylovNewton (the last at tol ×10). `LadrunoSANISAND`, the campaign set (G0 264.32, e_init 0.6944, Mc 1.3309, c 0.71, λc 0.027, e0 0.83, ξ 0.45, h0 1.3, ch 0.968, nb 3.5, A0 0.05, nd 5.75, zmax 12.5, cz 1100), `IntScheme 1`, `TanType 0`, `-flipAlphaIn init`, `-Pmin 0.0101`, `-maxSubsteps 2000`, `-Presidual 0`, `-honorTolR 0`. Build `79e062367`, Esmeralda (8 threads) unless noted.

**1.1 The deck is sound; the material is not.** On the identical deck, only the material changing:
- UW `DruckerPrager`, ψ = 0, at 33 / 38 / 43 / 48°: every leg reaches s/B 0.15 in 35–139 s with 0–3 subdivisions, carrying up to 1 262 kPa.
- `PressureDependMultiYield` at a 33° plane-strain cone: a clean limit point, 417.6 kPa at s/B 0.116, 38 min.
- `LadrunoSANISAND`, every variant — the full set; `zmax 0`; `nb 0, zmax 0`; `nb = nd = zmax = 0` (a critical-state sand): **every one stops on the step floor** ("ds below floor") between s/B 0.026 and 0.041 (0.048 with a 25 kPa surcharge, §1.4), carrying 279–814 kPa, after 3–5 h. The full set stops at 650 kPa, half what DP at 48° carries without subdividing.

**1.2 The profile** (the fork's own `profiler start -deep`, every case to s/B 0.002, same deck; ms per **accepted** Newton iteration):

| material | update | formTangent | formUnbalance | linearSolve | total | rungs N/LS/K |
|---|---:|---:|---:|---:|---:|---|
| linear elastic | 2.8 | 10.2 | 5.0 | 14.3 | 50.2 | 15/0/0 |
| DruckerPrager | 16.0 | 10.4 | 5.6 | 8.2 | 51.4 | 15/0/0 |
| PDMY, 38° plane-strain cone | 9.7 | 24.1 | 40.3 | 7.4 | 86.8 | 11/4/0 |
| LadrunoSANISAND, nb = nd = zmax = 0 | **1 721.8** | 19.1 | 6.9 | 6.8 | **2 010.8** | 8/9/4 |

85.6 % of SANISAND's wall is inside `Domain::update()`. The solve is never the constraint (7 % with a real soil). Substep snapshot at one committed step over 9 720 points: total 229 649, mean 23.6, max 666 — a snapshot, because no response accumulates.

**1.3 Where it stops.** Per-Gauss-point field dumps every 10 steps to the wall:
- Full set, last dump before the wall (s/B 0.036): the ten highest η (1.90–2.00, on the bounding surface) all sit in the **top row of elements outside the footing**, x/B 0.69–2.06, at p' ≈ 3.5 kPa. Not the corner.
- Critical-state set (s/B 0.041): the extremes lie along the shear band, 8 of 10 on one side of a symmetric deck; the lowest p' (3.39 kPa) at the surface just outside the edge.
- The explicit lane's failure state (attached, §4) puts the worst point at x = 0.844 m, y = −0.094 m (just outside the edge at x = 0.75 m) with **p' = 0.352 kPa and η = 12.87**.

**1.4 What does not move it.** `-Presidual 1.01` (the vanilla default) and `5.05` kPa: the wall lands at s/B 0.0263 and 0.0343, inside the baseline's own run-to-run range; at 5.05 kPa the extreme points simply move further out along the free surface (x/B 2.1–2.7). `-pRe` (F13): no effect. A **real** surcharge of 25 kPa relieves the ring, moves the extremes to the edge column going down, and postpones the wall to 0.048 — on a different, stronger deck. The act's reading: the shallow ring is where the wall *shows*; the cause is the material's integration at a high stress ratio and a few kPa of confinement.

**1.5 Tangents.** On `79e062367`, held at the declared 1e-5 with no relaxed rung, `TanType 0` reproduces itself to 0.5 %; `TanType 1` and `2` still give different curves (at s/B 0.009: 248.9 / 257.5 / 325.2 kPa) and stop earlier (0.0102, 0.0091). Only `TanType 0` is dependable here.

**1.6 Reproducibility.** The same leg, same spec, same build, run twice on different nodes: the wall at s/B 0.0278 (650.8 kPa) and 0.0365 (790.6 kPa); the critical-state set at 0.0320 and 0.0413. Only the node and the MKL thread schedule changed.

## 2. What the act read in the source (verify each)

- `SRC/material/nD/UWmaterials/ManzariDafalias.cpp:1490` — `TolE = mHonorTolRInME ? mTolR : 1e-4`: under `-honorTolR 0` the substepper runs on the hardcoded 1e-4, not the declared `TolR`.
- The substep error (the block around `:1866-1876` on this checkout): `err = ‖dσ₂ − dσ₁‖ / (2‖σ‖)` when `‖σ‖ ≥ 0.5`, and the absolute `‖dσ₂ − dσ₁‖` below it — **relative to the stress itself**, with an abrupt switch at 0.5 kPa. At p' ≈ 3.5 kPa that asks for an absolute stress error of about 1e-3 kPa.
- `:1766` and `:1839` — each substep forms two full 6×6 elastoplastic tangents (`aCep1`, `aCep2`). They are the two stages of Heun's method, so they are not redundant; but each stage only needs the stress increment `Cep:dε`, a 6-vector.
- SANISAND's plastic modulus is proportional to p and its elastic moduli to √p: at a few kPa with η on the bounding surface the rate equations are stiff, and an explicit scheme is stability-limited there, not accuracy-limited.
- F12 (`Ladruno_files/testbed/hypo_bearing/adr92_f12/F12_intscheme2_verdict.md`): `BackwardEuler_CPPM` is 4–30× more accurate and 4–13× cheaper **per given increment**, its `TanType 2` is the algorithmic tangent (§5.1), but under a global Newton an off-path trial iterate sends it through up to 2⁹ recursive halvings (§4, §5.2; 12–134 s per failing step), and nothing reaches the element as a refusal (§5.3).
- `Ladruno_implementation/107_ladruno_openmp_element_loop.md:308` — SANISAND refused from the threaded loop: `IntScheme 1` segfaults despite a clean static audit; `IntScheme 2` has shared static work arrays in `NewtonIter()`. The implex guide documents `implexRefusals` and `avgImplexError` as process-wide (`LadrunoSANISAND_implex_guide.md:315-330`).
- No `OPS_PROFILE_SCOPE` exists inside `ManzariDafalias.cpp` or the PDMY family.
- `SRC/material/nD/soil/PressureDependMultiYield03.cpp:275-282` hard-codes `ei, cs1, cs2, cs3`.

---

## 3. Asks

### F18 — integrate SANISAND at low confinement (the core ask)

(a) **An error norm with an absolute floor** in `ModifiedEuler`: `err = ‖dσ₂ − dσ₁‖ / max(2‖σ‖, σ_ref)`, with `σ_ref` a new flag (e.g. `-errFloor <kPa>`, default reproducing today's behaviour byte-identically). Report, on the attached ring states (§4) and on a single point driven to η/M^b → 1 at p₀ = 2, 5, 20 kPa: substeps per increment and the stress error against a reference integration (`IntScheme 45` at a tight `TolR`, or `ModifiedEuler` at 1e-8) for σ_ref ∈ {0, 0.1, 1, 5} kPa. The floor is only acceptable if the error it admits is smaller than the global Newton's own tolerance at that point; say where that holds.

(b) **Rate-form stages.** Inside the substep loop compute `dσ = C:dε − Λ·C:m` directly instead of forming `aCep1`/`aCep2`, and form the 6×6 tangent once at the end for `TanType` 1/2. Results must match today's to round-off; report the per-substep cost before and after (F20's scopes).

(c) **Make `IntScheme 2` usable under a global Newton**, starting from F12's diagnosis. When `BackwardEuler_CPPM` cannot return a trial iterate, **refuse at once with the refusal code the element forwards** (F7's roster) instead of recursing 2⁹ times, so the global step is cut in milliseconds rather than minutes; add a line search or a better start to the local Newton; remove the static work arrays (also F19's blocker). Then rerun F12's bearing deck: is the global Newton quadratic with `TanType 2` on the solution path, and how deep does the leg go against `IntScheme 1` for the same wall clock?

(d) **The per-point fallback** (F10b(b), still open): when `ModifiedEuler` hits `-maxSubsteps`, hand that point's increment to `BackwardEuler_CPPM` and refuse only if that also fails — behind a flag, with the one-element test.

(e) **Say plainly if the ring state is one no integrator can take** (F10b(c)): the attached b8 point at p' 0.352 kPa, η 12.87. If η can exceed M^b by a factor of six at a committed state, that is a question about the committed state, not only the integrator — trace how it got there.

**Test:** (a) byte-identical at the default; (b) round-off identical; (c) F12's bearing-deck comparison, with iterations per step; the ring states integrated under each variant with substep counts. **Ledger:** quirks rows (the hardcoded TolE, the relative error with the 0.5 kPa switch), implementations rows, the SANISAND guide's integration section.

### F19 — SANISAND in the threaded state-determination loop

WP-107 measured `Domain::update()` at 51.2 % of the step on this deck shape and refuses SANISAND because it segfaults. **First deliverable, before any code: the inventory** of shared mutable state reachable from `setTrialStrain` / `commitState` in `ManzariDafalias` and `LadrunoSANISAND` — file-scope and static buffers, the process-wide ledgers (`implexRefusals`, `avgImplexError`), the warning budgets, anything in the `Pmin`/flip bookkeeping — with every write site. Then make them per-instance, thread-local with a reduction at the end of the phase, or locked where a count must stay global. **Test:** a deck of your own with a few thousand SANISAND points at 1/2/4/8 threads: identical committed curves to round-off, and the speed-up table. **Ledger:** implementations row; the 107 note's refusal table updated.

### F20 — make the cost measurable

(a) A per-instance **cumulative** substep counter (since `revertToStart`), plus the count of the last `setTrialStrain` and the number of cap hits, exposed like `implexRefusals`, stating which columns are per point. The explicit lane's failure dump (§4) read `substeps` = 0 at every point right after a failed `analyze` — the counter a post-mortem needs is exactly the one that is lost.
(b) `OPS_PROFILE_SCOPE`s inside the integration: elastic predictor, `GetStateDependent`, each tangent stage, `Stress_Correction`, the CPPM Newton. Report the split at a ring state and at a deep state.
(c) F16's second half: a `"tangentEP"` response returning the continuum elastoplastic tangent at the committed state regardless of `TanType`, checked against a numerical tangent at a plastic point.

### F21 — replay a dumped material state

The act cannot hand the fork its model, but it can hand it material states (§4). Today there is no way to put a `LadrunoSANISAND` point into a given state (σ, α, α_in, fabric z, e) and drive it with a strain increment. **Ask:** a material-point facility — a `setParameter`/response, a `-initState <file>` or a small Python helper over a one-element deck with zero free DOF — that loads such a state and applies a prescribed `dε`, and reports substeps, the error per substep (a trace of `T`, `dT`, `err`: F10b(a)) and the returned stress. Use it on the attachments for F18.

### F22 — a deterministic mode

§1.6: the wall moved 30 % between two identical runs. Document and expose a deterministic mode for the implicit static path: at least MKL's conditional numerical reproducibility for Pardiso (`MKL_CBWR`, `iparm[33]`), and a note of what else is order-dependent (the threaded loop once F19 lands, any OpenMP reduction). **Test:** the same deck twice on 8 threads, byte-identical curves. **Ledger:** guide paragraph; a banner line when the mode is on.

### F23 — PDMY, two small items

(a) Expose `PressureDependMultiYield03`'s `ei, cs1, cs2, cs3` as optional arguments defaulting to today's values (byte-identical), with a quirks row.
(b) A guide note with the act's evidence: in drained plane-strain compression, ten dense-sand candidates across PDMY01/02/03 (with and without a retuned critical-state curve, `dilationParam3` 0.6–2.0) all failed a saturation gate — the volumetric strain rate at the end never fell below 27 % of its peak, usually above 95 % — and on the strip the dense set drove p' under the footing from 19.7 to 1 652 kPa without a plateau. The act's reading, to be checked: PDMY's dilation brake is a switch on void ratio (`PressureDependMultiYield.cpp:2189-2214`) that a genuinely dense sand does not reach. No code change asked beyond (a).

---

## 4. Attachments (`_tims_2d_model_requests_2026-09-25/`)

- `ring_points_b8.csv`, `ring_points_b16.csv` — 40 points each with p' < 10 kPa, sorted by η/M^b (compression-side M^b = Mc·exp(−nb·ψ), g(θ) = 1), with the full state: σ (6, Voigt, kPa, compression negative as OpenSees stores it), α (6), α_in (6), fabric z (6), e, ψ, element, Gauss point, centroid x, y (m). Source: the explicit central-difference lane's failure dumps, Esmeralda jobs 147068 (B/8, failed at s/B 0.0207) and 147069 (B/16, s/B 0.0127), build `48c0e99bc`, read with the model in memory right after `analyze` returned −2; parameters as §1 with `nu = 0.312885` (the K0 substitution).
- `README.md` in the same folder states the column order.

## 5. Order of work

F20(a) and F21 first (a day, and they make everything after measurable), then F18(a) and (b), then F18(c)–(e), then F19. F22 and F23 are independent and small.
