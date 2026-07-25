# ADR-75 P1e — `system Pardiso -krylov <L>` (factorization-preconditioned CGS): MEASURED

**Verdict: this is the reuse lever that works under FULL NEWTON — the case P1a's `factored` gate
cannot touch — and its payoff GROWS with model size on the production thread count.**
**1.51× vs direct PARDISO at 50.7k DOF / 4 threads; 1.57× when stacked with P1d's half-storage.**
Tip displacement bit-identical throughout.

Date: 2026-07-25 · worktree `opensees-pr-strategy-9c17b8` · same binary for every row.

## Why this exists

P1a's factorization reuse gates the numeric phase on the SOE `factored` flag. That flag answers
"has A been re-assembled since the last solve", so it can only *skip* work when the tangent is
unchanged — it pays under `ModifiedNewton` / `Initial` / `Krylov` / IMPL-EX and pays **nothing**
under full Newton, which is what most decks actually run.

`iparm[3]` is the complementary axis. It keeps the previous L/U and uses it as a **preconditioner**
for a tangent that *has* changed: phase 23 runs CGS (Intel's `K=1`, `mtype 11`) or CG (`K=2`,
`mtype 2`, SPD only) and falls back to a real refactorization by itself if the iteration stalls.
P1a and P1e therefore cover disjoint algorithm families; neither subsumes the other.

## What was measured

Model: Lane B — n³ `LadrunoBrick` + `LadrunoJ2`, implicit static, 15 `LoadControl` steps,
**`algorithm Newton` (full)**, `numberer RCM`, `test NormDispIncr 1e-7`.
Harness: `p1e_krylov_bench.py`, median of 3 **interleaved** rounds (A/B/A/B…).
Correctness gate first: `p1e_smoke.py`; Tcl path in `p1d_tcl_smoke.tcl`.

Metric is **end-to-end `analyze()` wall time**, which includes assembly and element
state-determination — work `-krylov` does not touch. The solve-phase speedup is therefore
*larger* than every ratio below. Back-of-envelope: ADR-40b's ~66% `linearSolve` at 11.5k DOF would
put the 11.5k/4-thread row at ~1.4× on the solve alone — **but treat that as illustrative, not
measured**: ADR-40b's fraction is a single-threaded one, and the solve fraction is precisely what
threading shrinks, so applying it to a 4-thread run overstates the implied solve speedup. A real
number needs a factor-vs-triangular-solve profiler scope, which does not exist yet (ADR-40 Phase-0
gap).

Repeats are 3 (2 for the 60-step run), below P1d's 5–7. The headline 25³ rung was re-run in a
separate process as a check (§1).

## 1. Time — the two axes pull in opposite directions

| MKL threads | grid | DOF | direct (s) | `-krylov 6` (s) | speedup |
|---|---|---|---|---|---|
| 1 | 15³ | 11,520 | 14.259 | 9.426 | **1.51×** |
| 1 | 20³ | 26,460 | 46.825 | 24.143 | **1.94×** |
| 4 | 15³ | 11,520 | 9.099 | 7.430 | **1.22×** |
| 4 | 20³ | 26,460 | 24.606 | 18.954 | **1.30×** |
| 4 | 25³ | 50,700 | 60.258 | 39.964 | **1.51×** |
| 4 | 25³ | 50,700 | 60.416 | 39.872 | **1.52×** (independent re-run) |

The headline rung was run twice, in separate processes: **1.508× and 1.515×**, 0.5% apart. P1d's
caution about a volatile unsymmetric-PARDISO anchor does not bite here.
Raw data: `p1e_krylov_t1_n15-20.json`, `p1e_krylov_t4_n20_compose.json`, `p1e_krylov_t4_n25.json`.

**Threading erodes the win; model size restores it.** Numeric factorization parallelizes well; the
triangular solves inside each CGS iteration do not. So 1→4 threads makes the *direct* path
disproportionately cheaper and the 11.5k win collapses from 1.51× to 1.22×. That is not the end
state: at 4 threads the trend in N is **1.22× → 1.30× → 1.51×**, still rising at the largest rung
measured.

**This is the ADR-40b trap in miniature** — a verdict taken at 11.5k ("marginal, 1.22×") inverts by
50.7k. Do not quote the small-model number.

## 2. It composes with P1d's half-storage (sublinearly)

20³ / 26,460 DOF / 4 threads, all four rows from one interleaved sweep:

| config | median s | vs direct |
|---|---|---|
| direct | 24.459 | 1.00× |
| `-krylov 6` | 18.851 | 1.30× |
| `-matrixType 1` (P1d) | 18.224 | 1.34× |
| **`-matrixType 1 -krylov 6`** | **15.615** | **1.57×** |

1.30 × 1.34 = 1.74 if the levers were independent; the measured 1.57× is short of that, as
expected — both attack the same factorization cost, one by making each factorization cheaper and
one by doing fewer of them. **Legal only for `-matrixType 1`** (see risk 1).

Incidental reproducibility check: the `direct` anchor came out 24.606 s and 24.459 s in two
independent sweeps (0.6%). P1d had warned that the unsymmetric PARDISO row was its volatile term;
it was stable here.

## 3. The stopping tolerance `L` barely matters

At 1 thread / 20³: `-krylov 3` 23.074 s, `-krylov 6` 24.143 s, `-krylov 8` 24.517 s — a 6% spread
across five orders of magnitude of requested accuracy.

The `-stats` trace explains it: CGS converges in **1 iteration during the elastic steps and 3–5
once the model yields**, so a tighter criterion costs at most a couple of extra matvecs. `L` is not
a delicate choice, and there is no accuracy argument for the loose end — `-krylov 3` and
`-krylov 8` produce the same tip displacement. **Recommend `L = 6`** (Intel's own example value)
over the marginally faster `3`: the extra three digits are nearly free and keep the inexact solve
well inside the `1e-7` outer Newton tolerance.

## 4. Accuracy

Tip displacement relative error vs the direct PARDISO row: **0.0** at 15³ and 25³, **1.5e-16** (one
ULP) at 20³, at 1 and 4 threads alike. Under `ModifiedNewton` the `-krylov` run is **bit-identical**
to the direct `ModifiedNewton` run.

**Zero CGS fallbacks occurred on this model** — every phase-23 call was answered by the iteration
and `cgs_error` was never raised. That is a property of Lane B's well-conditioned J2 tangent and
must NOT be generalized; see risk 2.

## 5. The stale-preconditioner guard, tested adversarially

`factorsCurrent` exists because a CGS *win* returns a correct answer while leaving the stored L/U
one tangent behind, so the phase-33 shortcut must be forbidden afterwards. Failure mode is
`error = 0` plus a wrong answer, so the guard was verified by **deliberately removing it,
rebuilding, and re-measuring**.

It is load-bearing — the answer changes. But two obvious checks were measured **not** to catch it,
and that is worth recording so nobody re-derives them:

| check | with guard | guard removed | discriminates? |
|---|---|---|---|
| tip ux vs the UmfPack/Newton reference | 5.21e-10 | 5.33e-10 | **no** — both "ok" at any sane tolerance |
| ModifiedNewton iteration count | 93 | 94 (1.01×) | **no** |
| ...same, swept to perfect plasticity (`Hiso` 2000→0) | — | 1.03× throughout | **no** |
| **-krylov vs DIRECT under the same ModifiedNewton** | **0.0** | **1.2e-11** | **yes** |

Why the first two fail: a solve against stale factors is still a descent direction, so Newton's
residual test launders it into a slower-but-convergent quasi-Newton. And Lane B's tangent barely
moves — which is also why CGS converges in 1–5 iterations — so even at `Hiso = 0` the stale factors
stay nearly as good as the current ones.

**Consequence for the risk assessment:** on a Newton-family path the hazard is self-correcting and
mild. The exposure is on paths with **no outer iteration to launder it** — a single linear solve, a
sensitivity/DDM right-hand side, an eigen shift-invert — where a stale-factor solve is simply the
wrong answer, returned cleanly. The guard is justified by Intel's documented contract (a CGS win
does not recompute the factors), not by a measured Lane-B failure. `p1e_smoke.py` now asserts the
row that does discriminate, at a 1e-13 threshold.

## 5-bis. Preconditioner aging: hypothesised, tested, REFUTED

CGS never refactors while it keeps winning, so the preconditioner is still the *step-1*
factorization while the tangent walks away from it. The `-stats` trace appeared to show exactly
that drift — iterations creeping 1 → 3 → 4 → 5 across 15 steps — which would mean a 15-step
benchmark **overstates** the steady-state win. It does not:

| steps | 20³ / 4 threads, direct → `-krylov 6` | speedup |
|---|---|---|
| 15 | 24.459 → 18.851 s | 1.30× |
| 60 | 90.234 → 67.514 s | **1.34×** |

Run-length probe on the 10³ model, counting every phase-23 outcome:

| steps | solves | CGS wins | fallbacks | iterations (first 8 → last 8) |
|---|---|---|---|---|
| 15 | 39 | 39 | 0 | `1,1,1,1,1,1,1,1` → `3,4,5,5,4,5,5,5` |
| 60 | 147 | 147 | 0 | `1,…` → `4,5,5,5,4,5,5,5` |
| 150 | 340 | 340 | 0 | `1,…` → `5,5,4,5,5,5,5,5` |

**The iteration count plateaus at 4–5 and stops growing.** The drift is the elastic→plastic
transition, not unbounded aging: once Lane B reaches a steady plastic state the tangent stops
moving away from the step-1 factorization. No refresh policy is needed, and the 15-step figures
are not short-run artefacts.

## 5-ter. The fallback branch is UNTESTED — the top residual risk

Across every run above — 150 steps, 340 solves, `-krylov` 6/8/9, hardening swept to perfect
plasticity, load doubled to a tip displacement of 55.5 (deeply plastic, 18 CGS iterations) —
**PARDISO never once gave up.** Zero fallbacks, `cgs_error` never raised.

So the `iparm[19] < 0` branch has **never executed**: its decode (`-it_cgs*10 - cgs_error`), its
`cgs_error 5` advice, and its perturbed-pivot check are unverified at runtime. This is stated
plainly rather than papered over.

**Mitigation, by inspection:** the only *correctness*-relevant action in that branch is
`factorsCurrent = true`, and it is set unconditionally — it does not depend on the decode. That
assignment is right for both reachable cases: `iparm[19] < 0` (Intel: "the factors L and U are
recomputed for the matrix A") and `iparm[19] == 0` (no CGS attempted ⇒ phase 23 was a plain
factorize-and-solve). A wrong decode would therefore corrupt a *diagnostic string*, not an answer.
The residual exposure is a misleading message and an unproven `cgs_error 5` warning.

## 6. Risks and limits (read before enabling this by default)

1. **No CGS mode for `mtype -2`.** Intel documents `K=2` for symmetric **positive definite** only.
   `-matrixType 2` — the correct symmetric choice for a softening or buckling tangent, because such
   a tangent is genuinely indefinite — therefore **cannot use this lever at all**; it warns and
   falls back to direct. So the model class that most needs solver speed (nonlinear, softening) is
   exactly the one where `-krylov` is either unavailable (`-matrixType 2`) or facing a much harder
   preconditioning problem (`-matrixType 0`). The 1.57× composition result is an **SPD-only**
   number.
2. **Untested near a limit point, and the fallback path has never fired** (§5-ter). Lane B hardens
   monotonically. A tangent approaching singularity is where CGS should start failing, and the
   fallback cost (per Intel, up to about half a factorization wasted, then the factorization
   anyway) has **not been measured** — attempts to provoke it by tightening `L` to 9 and going to
   perfect plasticity under double load all failed to produce a single fallback.
3. **Inexact solves inside a nonlinear loop can cost outer iterations.** Not observed here — the
   answers are bit-identical and iteration counts match — but it is the mechanism that made MUMPS
   BLR *slower* in P2b. Judge any new model on total wall time, never on solve time alone.
4. **Not measured above 50.7k DOF**, and the trend was still rising there.

## 7. Recommendation

Ship **opt-in, off by default**, as `-matrixType` and `-BLR` are; Intel's own position is that
"other values are only recommended for an advanced user". The default costs nothing — `-krylov`
absent is byte-identical to P1d.

Worth enabling when: full Newton (or any tangent-refreshing algorithm), unsymmetric or SPD,
≳25k DOF, smooth material path. Worth *disabling* the moment `cgs_error 5` appears.
