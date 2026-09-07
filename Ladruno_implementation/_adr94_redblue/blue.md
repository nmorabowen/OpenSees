# ADR-94 R2 — blue verdict

Build `52314165a`. New evidence: `tests/test_adr94_redblue_blue.py` (2 tests,
0.2s, coexists cleanly with the 4 baseline + R1 files, 13/13 combined).

Per-H: (a) results vs cost, (b) existing mitigation, (c) where red overstates.

**H1 (static `Stiffness`).** (a) RESULTS UNAFFECTED, COST ONLY. Measured
directly: a 2-cube heterogeneous VonMises model (one driven plastic, one
elastic — H1 fully active) converges (`NormDispIncr` 1e-9) to stresses
matching each cube solved alone to 1e-6 relative — the residual is built from
each element's own committed stress, never from the shared static, so a
converged step cannot be wrong. Cost: 201 total Newton iterations together vs
124 (106+18) solved separately — a **62% iteration penalty**, not a wrong
answer. (b) Single-material meshes do NOT make H1 ordering-invariant in
general — only a mesh with genuinely IDENTICAL GP states does (test 2:
0-iteration-delta, bit-identical, order-swap invariant). Cerro Lindo's
horseshoe-cavity model has real stress gradients, so H1 is live there too —
but per (a) it cost iterations, not the reported defect (missing tensile cap /
inadmissible corner states, ADR-84 §1), which were real physics/algorithm
bugs, already fixed. (c) **This is where red is most wrong**: H1 is billed as
a blocker that corrupts results. It corrupts the ASSEMBLED TANGENT, which is
a Newton-direction artifact — every existing Cerro Lindo/M-series result
obtained with `analyze()==0` is exact regardless of H1, because convergence
already certifies the residual. The real cost is iteration count and, in the
worst case, non-convergence (`rc=-3`) — annoying, loud, and already visible in
existing logs, not a silent-wrong-answer class like H7/H10b.

**H4 (revert no-op).** (a) Only bites a workflow that reads material state
between a Newton failure/cutback and the next trial, or calls `ops.reset()`
expecting a real restart — a fixed-step-schedule deck never does either. (b)
`TenNodeTetrahedron`'s eleResponse self-heals, masking it further there. (c)
Red is right this is real but overstates universality — production decks
that don't adaptive-step or inspect mid-cutback state never exercise it.

**H5 / H7 (unguarded elastic-exit, BE inconsistency fallback).** (a) Both are
silent-accept, so RESULTS can be wrong by construction — not cost-only. But
H7's trigger needs `H_iso` steeper than `2G`, a softening regime MC/MCTC with
`H=0` (Cerro Lindo's actual usage) cannot reach. H5's non-default-integrator
sites only bite decks that explicitly select FE/FE_sub/ME/RK45 — Cerro Lindo
and the ADR-84 battery use `Backward_Euler` throughout. (b)
`strict_convergence 1` + `Backward_Euler` (the shipped default pair) already
gates the one path real decks use. (c) Red correctly calls these
major/blocker; blue's disagreement is narrower: severity for TODAY's actual
users is bounded by "don't use softening or non-BE integrators without
`strict_convergence`" — documentable, not a rewrite.

**H8 (`Backward_Euler_LineSearch`).** (a) Cost/robustness only for anyone who
selects it, and it is strictly worse than the default it's meant to harden.
(b) Mitigation = don't select it (document) or refuse it — cheaper than
fixing three independent bugs in an option nobody needs.

**H9 (ME/RK45 dead drift check).** (a) Cost/silent-drift only for ME/RK45
users; `return_to_yield_surface` is the working correction. BE (default
integrator) users are unaffected.

**H10 (HB drift, DP apex NaN).** (a) H10a changes results only in NET TENSION
on rock — compression-dominated tunnel/slope work (Cerro Lindo's actual
geometry) is bit-identical per R1-C's control test. H10b needs a pure
hydrostatic-tension path through the apex — a soil-unloading corner case, not
the general case. (c) Red is right these are confirmed; blue's point is
scope — both are upstream YF issues gated on coordination (D4), not fork
regressions to rush.

**H12 (cout diagnostics).** (a) Perf/hygiene only, zero effect on results.

**H13 (silent parser typos).** (a) Does not corrupt a CORRECT deck's results,
but makes a WRONG deck (typo) silently run with defaults — a **guardrail
gap**, not evidence any existing, already-validated deck is wrong today.

**H6 tangent finding.** (a) Convergence cost only — stresses are exact
(REFUTED accuracy claim on VM; BE lands on the closest-point answer when the
return direction doesn't rotate). Shipped default `Secant` costs 5.3x
`Continuum`'s iterations. (b) ADR-84 §9.4 already gives this exact guidance
(`tangent_type Continuum`) for MCTC from independent prior work — confirmation
of standing advice, not a new gap. (c) Changing the *default* is what ADR-84 §3's
bit-identical strategy explicitly avoided doing for the YF; the same logic
applies to the tangent default — a global default flip is upstream-facing
(changes every deck's iteration signature) where a one-line doc addition does
not.

## Blue-ranked action list

**Do now, opt-in, fork-only (cheap, no coordination):**
1. Document `tangent_type Continuum` as the recommended default for
   ASDPlasticMaterial3D generally (H6), not just MCTC (already ADR-84 §9.4).
2. Document `strict_convergence 1` + `Backward_Euler` as the required pairing
   for any softening (`H_iso<0`) or non-ideal-plastic MC/MCTC/VM deck (H5/H7).
3. Document/refuse `Backward_Euler_LineSearch` until fixed — it has no
   legitimate use case today (H8): a one-line refusal in the parser is
   cheaper and safer than the 3-bug fix, and is fork-local.
4. Add a loud parser error for unknown `Begin_Integration_Options` tokens and
   unknown model-parameter names (H13) — small, local, no behavior change for
   any correctly-spelled existing deck, and closes the fork's own recorded
   worst failure mode (ADR-84 §9.1).

**Worth fixing now, fork-local, but needs a build + regression pass (medium
risk, still opt-in via existing flags where possible):**
5. Gate H7's fallback under `strict_convergence` (currently the "eighth
   silent-accept site" that ADR-84 P2a intended to close).
6. H4 `revertToLastCommit`/`revertToStart` — real fix, bounded blast radius,
   needed before any adaptive-stepping/augmented-Lagrangian work leans on
   ASDP state fidelity.

**Send upstream / coordinate with jaabell first — do not fork-diverge:**
7. H1 (per-instance `Stiffness`) — touches every specialization's ABI-ish
   static layout; also the ADR-75b threading blocker, so it should land as
   one deliberate change with a full-battery bit-identical check, the same
   discipline ADR-84 §3 used for the YF, not a quick patch under this ADR.
8. H10a (HB composite) and H10b (DP apex, both the stub AND the dead call
   site) — both are joint-authorship YF/framework code; guessing a fix here
   risks silently diverging from jaabell's tree the way H10a already has.
9. H9 dead drift checks, H15 (BE vs ME/RK45 elasticity evaluation) — low
   urgency since `return_to_yield_surface`/BE-default already mitigate; batch
   with H1's framework pass.

**Documentation-only, no code change warranted:**
10. H2 (dangling `c_str()` — works by SSO luck, one-line fix but zero
    observed impact; fine to note and defer), H14 (`getCopy` `first_step`,
    latent under all current usage), H3 (SP/MP scope, already D3-decided).
