# `algorithm Newton -initial` re-assembles and re-factorizes an invariant matrix every iteration

**Reported by:** TIMs project (SFIM continuum reference model)
**Build:** Ladruno binary of 2026-07-25 02:29, built from `8a07e7642` (MKL PARDISO enabled) — the hash the
executable reports at startup.
**Verified against:** `origin/ladruno` at `bdf8adf9f` (which contains `8a07e7642`; the three commits
between them are two merges and a docs change, touching neither `SRC/analysis` nor `SRC/system_of_eqn`).
Working tree confirmed level with origin, not a stale checkout.
**Severity:** performance, not correctness. Converged results are unaffected.

---

## 1. Summary

`algorithm Newton -initial` performs a full tangent assembly **and** a full numeric
factorization on **every iteration**, of a matrix that in most configurations cannot have
changed since the previous iteration. In a static analysis the matrix in question is
constant for the entire analysis, yet it is rebuilt and re-factorized at every iteration of
every step.

The engine already has the machinery to avoid the factorization — every sparse solver in the
fork gates its numeric phase on the SOE `factored` flag (`UmfpackGenLinSolver.cpp`,
`PARDISOGenLinSolver.cpp:228`, `MumpsSolver`). The flag is cleared because the matrix is
re-assembled. The missing piece is one level up, at the assembly layer: nothing asks whether
the assembly was necessary.

---

## 2. Where it happens

`SRC/analysis/algorithm/equiSolnAlgo/NewtonRaphson.cpp:163-190` — `formTangent` sits inside
the iteration `do`-loop. The only tangent flag given special treatment is
`INITIAL_THEN_CURRENT_TANGENT` (`:166`), which correctly forms the initial tangent on
iteration 0 and the current tangent thereafter. `INITIAL_TANGENT` falls through to the
generic branch and is re-formed unconditionally on every pass.

Each re-formation costs, in order:

1. `theSOE->zeroA()` and a full loop over every `FE_Element` calling `getTangent()` and
   `addA()` — `IncrementalIntegrator::formTangent`, `:100-127`;
2. the clearing of the SOE `factored` flag, which forces the solver to redo its numeric
   factorization on the next `solve()`.

Contrast `ModifiedNewton::solveCurrentStep` (`ModifiedNewton.cpp:139-149`), which forms the
tangent **once before** the loop.

## 3. Why the work is redundant — and where it is not

This is the part that must not be over-claimed; the answer differs by integrator family.

**Static analyses — redundant, always.**
`StaticIntegrator::formEleTangent` (`:69-83`) under `INITIAL_TANGENT` is exactly
`zeroTangent(); addKiToTang();`. The assembled matrix is a pure function of the elements'
initial stiffness. It does not depend on the current state, on the load factor, or on the
iteration index. It is invariant not merely within a step but for the whole analysis, and it
is nonetheless assembled and factorized on every iteration. A 1 000-step pushover at 5
iterations/step performs 5 000 assemblies and 5 000 factorizations of one matrix. Since
`algorithm Newton -initial` is a standard robustness fallback in pushover work, this is a
well-travelled path, not an exotic one.

**Transient analyses — conditional.**
`Newmark::formEleTangent` (`:295-298`) under `INITIAL_TANGENT` builds

```
K_eff = c1*Ki + c2*C + c3*M          c1=1, c2=gamma/(beta*dt), c3=1/(beta*dt^2)
```

`c1..c3` depend only on `dt`; `M` is constant; but

```
C = alphaM*M + betaK*Kt + betaK0*Ki + betaKc*Kc        (Element::getDamp, SRC/element/Element.cpp:211-231)
```

so when `betaK != 0` the damping term carries the **current** tangent `Kt`, the matrix
genuinely changes between iterations, and the re-formation is legitimate work.

Therefore: **redundant for all static analyses; redundant for transient analyses with
`betaK == 0`; genuine work when `betaK != 0`.** Any fix must respect that split — an
unconditional "form once per step" would silently alter the iteration for stiffness-
proportional damping.

A related observation, offered because it surprised us and may be worth documenting: with
`rayleigh $a $b 0 0` the `-initial` flag is very nearly inert. In our model
(`Newmark(0.5,0.25)`, `dt=1e-3`, `zeta=0.5`, periods 0.5 s / 0.2 s) we get `c1 = 1`,
`c2 = 2000`, `b = 0.02274`, so the effective tangent is

```
1.0 * Ki  +  45.5 * Kt  +  4.02e6 * M      (INITIAL_TANGENT)
46.5 * Kt              +  4.02e6 * M      (CURRENT_TANGENT)
```

`-initial` exchanges about **2 %** of the stiffness content; the rest re-enters through the
damping regardless. Users selecting `-initial` for robustness under stiffness-proportional
Rayleigh damping are very likely not getting the initial-stiffness iteration they think they
are, while paying full price for it.

## 4. Measured cost

SFIM shallow-footing model: 38 984 DOF, 9 248 `SSPbrickUP` (u-p) + `ShellMITC4`,
`Newmark(0.5,0.25)`, MKL PARDISO, 8 threads, i7-13700KF. Unit costs measured by differencing
profiled `Newton` against `ModifiedNewton` runs (both perform 2 triangular solves per step;
`Newton` performs one extra assembly and one extra factorization):

| Operation | Cost |
|---|---:|
| One tangent assembly | 325 ms |
| One numeric factorization (PARDISO, 8 threads) | 255 ms |
| One triangular substitution | 32 ms |
| Whole step, `ModifiedNewton` | 934 ms |

The two spellings were then run directly against each other on this build, same model, same
solver, same thread count, same tolerance, 25 profiled push steps:

| `algorithm` | `formTangent` calls/step | iters/step | s/step | settlement (m) |
|---|---:|---:|---:|---|
| `ModifiedNewton -initial` | 1.00 | 2.00 | **0.936** | −0.000974399 |
| `Newton -initial` | 2.00 | 2.00 | **1.504** | −0.000974399 |

Identical converged settlement, identical iteration count, **1.61× the cost** — 568 ms per
step spent re-assembling and re-factorizing to arrive at the same place. On our previous
UMFPACK build, note 21 of our own performance series measured `Newton -initial` at 6.56 s/step
against `ModifiedNewton` at 3.88, the same ratio the operation counts predict.

This is a transient run with `betaK != 0`, i.e. the case of §3 in which the re-formed matrix
genuinely does differ between iterations. Even here — where the extra work is *defensible* —
it buys nothing: the same answer in the same number of iterations. In the static case, where
the matrix is provably identical, there is not even a theoretical argument for it.

## 5. What we would like

Ranked by value over risk.

**R1 — Document the trap. (zero risk, immediate)**
State in the `algorithm` documentation that `Newton -initial` re-forms and re-factorizes at
every iteration, and that `ModifiedNewton -initial` provides an initial-stiffness iteration at
one assembly and one factorization per step. Users reach for `Newton -initial` as *the* robust
fallback without knowing it is the expensive spelling of it.

**R2 — The real fix: let the engine know the matrix has not changed.**
Introduce a monotone *tangent version* counter on the SOE or `AnalysisModel`, incremented by
`domainChanged()`, by a change in `dt`, by anything altering initial stiffness
(`updateMaterialStage`, material/section replacement), and — when any element carries
`betaK != 0` — by state update/commit. Any algorithm can then compare the version it last
assembled against the current one and skip `formTangent` when they match. This subsumes the
special case, is not specific to `-initial`, and would additionally allow a static
`Newton -initial` pushover to factorize exactly **once** for the whole analysis. It also
composes with the existing `factored` flag rather than duplicating it.

**R3 — Cheap intermediate, if R2 is too large.**
In `NewtonRaphson`, when `tangent == INITIAL_TANGENT`, skip `formTangent` for iterations > 0
*only when the integrator declares the initial tangent state-independent*. Add a virtual
`IncrementalIntegrator::initialTangentIsInvariant()` returning `true` for `StaticIntegrator`,
and for transient integrators `true` iff no element in the model carries `betaK != 0`. Please
do **not** make this unconditional: with `betaK != 0` it changes the iteration silently.

**R4 — Parser: allow `-initial` and `-factoronce` together.**
`OPS_ModifiedNewton` (`ModifiedNewton.cpp:53-86`) reads exactly one option string through an
`else if` chain, so `algorithm ModifiedNewton -initial -factoronce` silently honours only
`-initial` and drops the rest. The combination is the natural way to express "initial
stiffness, assembled and factorized once", which for a static analysis is arguably the
cheapest robust algorithm the framework can offer. Loop over the remaining arguments instead
of reading one.

## 6. Acceptance tests we would find convincing

1. **Static pushover, `Newton -initial`** — `formTangent` call count for the analysis equals
   1 (or the number of domain changes), and the displacement history is bit-identical to the
   current build.
2. **Transient, `betaK == 0`, `Newton -initial`** — `formTangent` calls == 1 per step;
   history identical to the current build.
3. **Transient, `betaK != 0`, `Newton -initial`** — behaviour **unchanged**, i.e.
   `formTangent` still called per iteration. This is the guard against an over-eager
   optimisation, and is the test we would most want to see.
4. **`algorithm ModifiedNewton -initial -factoronce`** — both options take effect.

Note for whoever picks this up: because `K_f` determines only the *path* to equilibrium and
not the equilibrium itself, none of this can change a converged result. It changes iteration
counts and cost. Tests 1 and 2 assert identity because in those configurations the matrix is
provably unchanged, so even the path must be preserved.

## 7. Unrelated, noticed in passing

- ~~`PARDISOGenLinSolver.cpp` carries no profiler scopes, where `UmfpackGenLinSolver.cpp` has
  `soe.symbolic` / `soe.factor` / `soe.trisolve`. On a PARDISO run `linearSolve` is a single
  opaque block, so factorization cannot be separated from substitution without differencing
  two algorithms. Instrumentation parity would be welcome.~~
  ✅ **FIXED 2026-07-27, [#667](https://github.com/nmorabowen/OpenSees/pull/667)** — this report
  is what put it on the list. PARDISO now has the three UmfPack names plus `soe.cgs` (phase 23),
  `dc.s.fill`/`dc.s.verify`, and a deep `soe.addA`.
- ~~The profiler HDF5 run attributes record `threads=1` and `nElem=0` regardless of
  the actual run configuration (observed with `MKL_NUM_THREADS=8` and 9 325 elements).~~
  ✅ **FIXED 2026-07-27 (ADR-75 P1i).** Root cause: `Profiler::buildMeta()` set
  `m.threads = threads_.size()` — the count of threads registered with the *profiler*, which is 1
  for any single-threaded command layer no matter what MKL is doing; `nElem`/`nNode` were promised
  by a comment in the same function and populated by nobody. Both reports in this section were
  correct and both are now closed.
