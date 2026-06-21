---
title: "ADR-30 P6 — near-singular LᵀML condition gate + frozen-Ccr staleness guard"
project: Ladruno
type: design
status: design (pre-code)
related:
  - "[[30_ladruno_explicit_constraint_projection_adr]]"
  - "[[projection_handler_handoff]]"
tags: [adr, constraints, explicit-dynamics, projection, robustness, condition-number]
---

# P6 — two robustness guards for LadrunoProjection

Backlog items folded into one phase (both are warn/refuse-only — NO change to the projection
math, so low-risk). Closes ADR-30 open item **O1**.

## 3a — near-singular `LᵀML` condition-number gate (projector)

**Gap.** `LadrunoConstraintProjector::buildMass()` builds `LtML = LᵀML` (SPD nRet×nRet) per group
and guards singularity with a SINGLE test solve `tmp.Solve(e0,x0)` — which only fails on an
**exact** zero pivot. A *near*-singular `LtML` (a barely-dependent retained direction, a tiny
relative mass, a near-redundant constraint) passes the solve but makes `project()` amplify
round-off into a garbage acceleration — silently wrong results. ADR O1.

**Fix.** Replace the test-solve with a real condition-number estimate. `LtML` is symmetric (SPD
for positive masses + full-rank L), small (nRet typically 1–6, ≤ a few tens for a big diaphragm),
so a self-contained **cyclic Jacobi** symmetric eigensolve (no LAPACK-symbol dependency — matches
the fork's OpenSees-free-kernel style; sidesteps the `dsygv_`-missing gotcha) gives all eigenvalues
cheaply. Then `cond = λmax / λmin`:
- `λmin ≤ 0` (numerically `≤ ε·λmax`) **or** `cond > COND_REFUSE (1e12)` → **refuse** with a named
  error (supersedes the exact-pivot test-solve; the old message about "a retained direction has no
  mass / redundant constraint" is retained and extended with the measured cond).
- `cond > COND_WARN (1e8)` → **warn once per buildMass** (the projection still runs; the user is
  told the constraint set is ill-conditioned).
- Keep the existing per-DOF massless guard (it fires earlier with a more specific recipe).

Self-contained: a `static int symEigJacobi(const Matrix& A, Vector& eig)` in the projector .cpp.
No API change; `buildMass` already returns int.

**Test.** Build a group whose `LtML` is near-singular (two retained DOFs coupled so the effective
mass matrix is near rank-1, e.g. an equalDOF chain through a near-massless intermediary, or a
direct unit-test of the eigensolver + a constraint model with a mass ratio ~1e14) → refused with
the cond in the message; a normal equalDOF/diaphragm group → passes; a moderately ill-conditioned
one → warns but runs.

## 3b — frozen-`Ccr` runtime staleness guard (handler)

**Gap.** The P2 transport bakes the rigidLink-beam / rigidDiaphragm lever arms into `L` from
`Ccr = mpPtr->getConstraint()` ONCE (small-rotation, frozen — the same limitation `Transformation`
carries). For accumulated master rotation ≳ 0.1 rad the lever-arm coupling is stale and the tie
drifts, silently. We can't refresh `Ccr` cheaply (that's the RBE2 element / finite-rotation route),
but we CAN warn the user.

**Fix.** Detect the at-risk couplings at `handle()` and watch the master rotation at runtime.
- **Flag at transport:** in the MP edge loop (where `c = Ccr(i,j)` for retained DOF `rDofs(j)` of
  master `rNode` is read), a coupling is staleness-prone iff `c != 0` AND `rDofs(j)` is a
  **rotational** DOF of the master. Rotational test by OpenSees convention from `ndm =
  node->getCrds().Size()` and `ndf = node->getNumberDOF()`: `(ndm==2 && dof>=2) || (ndm==3 &&
  ndf>=6 && dof>=3)`. Record the master node tag (deduped), its rotational local dofs, and capture
  `θ0` = the master's committed rotational displacement at build time.
- **Warn at runtime:** in `applyLoad()` (already the per-step handler hook, P4b), for each monitored
  master read the current rotational disp `θ`, compute `drift = max_k |θ_k − θ0_k|`, and if
  `drift > STALE_WARN_RAD (0.1)` and not-yet-warned → emit a ONE-TIME warning naming the node and
  drift, pointing to RBE2 (`LadrunoKinematicCoupling`) / `Transformation` for finite rotation.
  Early-out when there are no monitors or all have warned.

New handler state: `struct RotMonitor { int nodeTag; std::vector<int> rotDof; std::vector<double>
theta0; bool warned; }; std::vector<RotMonitor> rotMonitors;` cleared on every rebuild (mirrors the
other per-build vectors). equalDOF (Ccr = I) and rigidLink-bar (translation-only) are NOT flagged.

**Test.** A rigidLink-beam (or rigidDiaphragm) whose master node is driven to a large rotation
(>0.1 rad) → the staleness warning fires once; the same model kept at small rotation → no warning.
(Captured via a stderr/opserr capture, like the other refusal/warn tests; or assert the run
completes and the tie error grows — the *behavior* is unchanged, only the diagnostic is added.)

## Scope / bookkeeping
- Files: `SRC/analysis/handler/LadrunoConstraintProjector.cpp` (3a), `LadrunoProjectionHandler.{cpp,h}`
  (3b). Both fork-authored → ZERO vanilla. No new classTag. No banner change.
- Tests: `tests/test_adr30_projection_p6.py`.
- Ledgers: impl-row note + a LEDGER_quirks entry (cond gate thresholds; Ccr staleness is a warn not
  a fix). Handoff §P6 + remove items from the backlog. ADR O1 marked closed.
- Verification gate: per the new gate-policy memory, this touches the projector core math path
  (the cond gate decides whether project() runs), so a focused review/verify is warranted — but it
  is guard-only (no change to the a_proj computation), so a lighter pass + the test battery suffice.
