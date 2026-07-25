---
title: ADR-76 session handoff — what shipped, what is wrong with the ADR, what is loose
project: Ladruno
status: handoff
owner: nmora
tags: [handoff, adr, adr-76]
---

# ADR-76 handoff (2026-07-25)

## Shipped and verified

| Commit | Change | Gate |
|---|---|---|
| `1a1502bb7` | `OPS_ModifiedNewton` parses **all** options (`-initial -factoronce` used to drop one), + camelCase/`-Initial`/`-Secant` spellings, factor resets, unknown-token warning | `adr76_smoke/` 11/11 |
| `804c4bda3` | **LAPACK singular matrix reported SUCCESS.** `return -info+1;` == `(-info)+1`, so `info==1` returned 0. X had been copied from B before the call, so callers consumed the LOAD VECTOR as a displacement. Fixed in BandGen/FullGen/BandSPD | `lapack_singular_regression/` 8/8 |
| `0ce0b08fb` | Doc self-consistency; generated artifacts stripped from history | — |
| `074aff56f` | `OPS_Algorithm` returns −1 on a null (closes the silent-swallow for EVERY algorithm factory); two false claims corrected | rebuild + both decks |

## THE ADR IS NOT IN A GOOD STATE — this is the main job

`76_ladruno_tangent_reuse_adr.md` was written, then amended repeatedly in place as
five adversarial reviews landed. Amendments were applied **locally** and nothing
downstream was swept. A reader meets, in order: a §3 table, 150 lines of §4
design, a §5 that §8 calls wrong, a §6 battery of seven tests where six gate
withdrawn work, a §7 of questions §8 answers — and only then §8, at line ~600.

Patched already (do not redo): title, §1 lead, §3 table rows, §4 SUPERSEDED
banner, the user-facing rule, the §6 over-claim + `8/8`→`11/11`, the 9.4e-10
mislabel, the `:155`→`:166` citation, both ledger contradictions.

**Still wrong — needs a REWRITE, not more patching:**

- **§4** — collapse to one paragraph + move the body to an appendix. §4.5's
  P0/P1/P2/P3 rollout is the most actionable-looking stale text in the file.
  **Keep §4.4's default-`false` box** — it is the best-argued passage and is
  correctly reconciled with §2.
- **§5** — §8.1 states `-krylov` is *antagonistic*, not complementary. §5 says the
  opposite and carries no marker. Fold what survives into §8.3, delete the rest.
- **§6** — split. Test 4 + its table is a real shipped result: keep as the
  acceptance record. Tests 1/2/3/5/6/7 move under §8 as "what would have had to
  be proven".
- **§7** — §8 answers every question in it. Fold into §8.3.
- **§9 provenance** — cites `ISSUE — Newton-initial refactorization.md`, which is
  **not in the repo**. Either commit it or mark it external/unarchived.
- ADR-76 never mentions the LAPACK fix, though both ledgers call it an "ADR-76
  audit spin-off". It needs its own section or its own ADR.
- One-way links: `40_ladruno_performance_adr` (which ADR-76 `amends:`) and
  `75_ladruno_sparse_direct_strategy_adr` contain no reference back.

## Loose threads (nowhere else durable)

1. **`Linear.cpp:110` opens the `formTangent` profiler scope OUTSIDE the
   `factorOnce` guard**, so the call count is `nSteps` regardless of whether an
   assembly happened. `ModifiedNewton.cpp:198` gets it right — a 1-line fix with
   a known-good sibling. Every `Linear -factorOnce` profile ever taken overstates
   assemblies. **Not in any ledger.**
2. **PR2 — `factorOnce` has no `domainChanged` reset** in `Linear`,
   `ModifiedNewton`, `ExpressNewton`. Designed but unstarted; see §8.3 + the
   quirks entry. The correct shape is an `invalidateTangent()` helper called from
   `domainChanged()` **and** every negative return (the early-return holes are
   not accompanied by a domain change). Also harden `sendSelf` — all three
   serialize the mutable latch. **Covers only the topology subset**: Δt change,
   staged construction, `updateMaterialStage`, `region -rayleigh` and modal
   damping never reach `domainChanged()`, so the ledger warning must survive.
3. **`CorotCrdTransf3d::T` is `static`** (`.h:135`) — cross-element aliasing in
   vanilla. Flagged in quirks, **not queued in `upstream_pr_campaign.md`**.
4. **ADR-75 PARDISO profiler scopes** — ADR-76 §5 claims ADR-75 tracks this. It
   does not. Track it or drop the claim.
5. **Neither deck is in CI.** `lapack_singular_model.tcl` is dependency-free and
   self-verifying — wiring it into Zone-A is ~5 lines and the highest-value item
   left. Gate on the terminal marker, never the exit code (`OpenSees.exe` exits
   **0** after a Tcl parse error — recorded quirk).
6. The smoke checker asserts mutual equality but **no absolute anchor**; 7.0 is
   analytically available and free. A regression driving every case to the same
   wrong answer passes.
7. Two dormant threaded solvers (`BandSPDLinThreadSolver.cpp:189`,
   `ThreadedSuperLU`) keep the same singular-as-success defect via `return info;`.
   Unbuilt, so consequence-only — worth one ledger line.
8. Five `LEDGER_vanilla_files` rows have `PR = —`; fill after merge.

## Corrections to earlier findings — do not re-propagate

- **`LadrunoIMKBeam`/`LadrunoIMKBeam2d` are CLEAN.** An earlier session note called
  them dirty; they use `getInitialTangent()` via `ladrunoIMKInitBlock` and no
  current geometry. Never committed, but do not repeat it.
- **`LadrunoDispBeamColumn2d/3d` are clean EXCEPT** under `-damping` (the
  multiplier reads domain time + the `OPS_GetStaticAnalysis()` singleton) or a
  `CorotCrdTransf3d`. Not intrinsically dirty.
- **Do NOT freeze `Ki` at the reference configuration.** This was implemented and
  reverted. `betaK0 * getInitialStiff()` enters the **residual**
  (`Element::getRayleighDampingForces`), so it changes converged answers; a frozen
  `K0` does not annihilate rotated rigid-body modes and injects a spurious axial
  damping force ≈ `β·EA·θ̇·sinθ`; and `-initial` then diverges past ~2-8° of chord
  rotation. If the aliasing is fixed at all, fix it in `CorotCrdTransf3d`.
