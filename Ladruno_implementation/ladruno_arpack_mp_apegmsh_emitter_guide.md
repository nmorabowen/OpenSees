---
title: Distributed ARPACK modal emitter — implementation guide for apeGmsh
project: Ladruno
status: ready to implement — the fork dependency is MERGED (PR #668)
priority: medium
adr: ADR-1000 §34 · apeGmsh ADR 0077
tags:
  - apegmsh
  - modal
  - mpi
  - emitter
  - arpack
  - mumps
---

# Distributed ARPACK modal emitter — implementation guide for apeGmsh

**What we are asking for:** a second Tier-1 backend in apeGmsh's parallel modal
surface — a **partitioned** deck that runs plain `eigen` over
`MumpsParallelSOE` under `OpenSeesMP`, alongside the existing replicated FEAST
deck.

**Who this is for:** whoever owns emission and the modal result surface in
apeGmsh. No fork C++ knowledge is needed. The normative sources are
[[1000_ladruno_cms_adr]] §34 (the wiring and its scope holes) and apeGmsh
ADR 0077 (the Tier 0 / Tier 1 structure this extends).

**Unlike [[LadrunoCMS_apegmsh_emitter_guide]], this one is not on hold.** The
fork dependency is merged and live-verified; the only thing standing between
apeGmsh and this capability is the emitter.

---

## 0. Why this exists, in three sentences

Until 2026-07-27 the only *correct* distributed modal path in this ecosystem
was FEAST, because plain `eigen` under `OpenSeesMP` built its `ArpackSOE`
unwired — `processID` stayed −1, the `M*v` merge never fired, and a
partitioned deck returned an empty spectrum on rank 0 while rank 1 deadlocked.
apeGmsh ADR 0077 correctly refuted that route as "never emitted" (finding F1).
**Fork PR #668 wired it**, so the refutation no longer holds against a fork
build carrying that commit, and apeGmsh's ADR needs amending along with the
emitter.

**Why you would want it even though FEAST works:** FEAST is *replicated* —
every rank assembles the full `(K, M)` CSR, and only the factorization is
distributed. This route is *partitioned*: each rank holds only its slice of
the model **and** the factorization is distributed. For the memory-bound case
that is the entire reason to go parallel, this strictly dominates. FEAST keeps
two advantages it does not lose: band-targeting (ask for a frequency window
rather than the lowest N) and `-certify` completeness.

---

## 1. The deck we need emitted

A **partitioned** deck — the ordinary `_emit_partitioned` output, `if
{[getPID]==K}` blocks and all — plus this preamble and solve. Contrast with
`modal_deck`'s FEAST path, which deliberately emits **flat** (ADR 0077 INV-3);
here partitions are the point, not an obstacle.

```tcl
# ... partitioned model blocks, as emitted today ...

constraints Transformation      ;# forced, see §2
numberer ParallelPlain          ;# the existing partitioned default
system Mumps                    ;# LOAD-BEARING here, see §2

set _lam [eigen $numModes]      ;# ONE capture, never a second [eigen ...]
if {[getPID] == 0} {
    set _fp [open eigenvalues.out w]; puts $_fp $_lam; close $_fp
}
```

Then the per-rank mode-shape harvest of §3.

### Runtime requirement to state in the docstring

This deck needs an `OpenSeesMP.exe` built from `ladruno` at or after
`5a522b03b` (PR #668). Against an older binary it does not degrade — it
**hangs or returns an empty spectrum**. There is no runtime feature probe;
treat the build as a deployment precondition and say so where a user will read
it.

---

## 2. Invariants — each of these is load-bearing

- **INV-A1 — `system Mumps` is genuinely load-bearing.** This is the exact
  opposite of the FEAST deck, where ADR 0077 INV-4 records that the `system`
  line plays no part in the solve (the RCI kernel owns its own dmumps). Here
  the fork wires the `ArpackSOE` collectives **only when the analysis
  `LinearSOE` is `MumpsParallelSOE`**. Emit any other system and the wiring
  silently does not engage — back to the broken path. Do not let a user
  override this on a partitioned modal deck; fail loudly instead.
- **INV-A2 — partitioned, not flat.** The whole point. Do **not** reuse the
  FEAST path's `supports_partitions = False` seam.
- **INV-A3 — `constraints Transformation`, forced.** Same reasoning as ADR 0077
  INV-4: `Penalty` pollutes M with penalty mass and `Lagrange` injects
  zero-mass DOFs, either of which fabricates spurious modes. The existing
  auto-emit only fires when MP constraints exist, so force it here.
- **INV-A4 — a single captured `eigen`.** Identical to ADR 0077 INV-5. A second
  `[eigen ...]` inside the rank-0 write-out is a redundant distributed solve
  *and* a rank-0-only collective ⇒ deadlock.
- **INV-A5 — no `modalProperties`.** Still MPI-blind upstream; #668 does not
  change that. The existing `ParallelModalResult` guard (raises on
  `participation_factors` / `mass_ratios`) applies unchanged.
- **INV-A6 — nodal-mass owner discipline.** See §5. This one has no code guard
  anywhere in the stack and apeGmsh is the right place to enforce it.

---

## 3. Mode-shape harvest — this is where you cannot copy FEAST

**The FEAST harvest does not transfer.** ADR 0077 P3 gets to use an ordinary
rank-0 recorder precisely *because* the FEAST deck is replicated — rank 0 holds
every node, so one recorder captures the whole field. Here each rank holds only
its partition, so a rank-0 recorder would silently capture **only rank 0's
slice** and `mode_shape_field` would return a partial mode with no error.

What we need instead — per-rank harvest plus a client-side merge:

1. Each rank writes its **own** sidecar, `mode_shapes_rank<P>.json`, carrying
   the node tags **it owns**, in its own recorder column order, plus `ndf` and
   `ndm` (keep the existing key names so the reader stays close to the FEAST
   one).
2. Each rank creates one recorder per found mode over its own nodes:
   `mode_shape_<k>_rank<P>.out`.
3. Same mechanics as FEAST P3, and for the same source-verified reasons:
   create the recorders **after** the captured solve (the eigen dataFlag reads
   `Node::getEigenvectors()` only at record time, and `Domain::addRecorder`
   does not auto-fire), fire them with a single `record`, close with `remove
   recorders`.
4. `ParallelModalResult.from_job` gains a partitioned branch: glob the per-rank
   sidecars, concatenate, and sort by node tag into the same
   `(n_nodes, ndf)` layout the replicated path produces.

### The boundary-node check you get for free

Nodes on a partition boundary are defined on **both** touching ranks, so they
appear in two sidecars. Do not silently take the first. **Compare them** — the
eigenvector components must agree to solver tolerance, because both ranks
solved the same global problem. This is a genuine correctness assertion:
disagreement means the wiring is not live (each rank ran a private Lanczos),
which is exactly the failure mode this whole slice exists to prevent. Raise on
mismatch, with the offending node tag and the two values.

That check is the cheapest available proof that the run was really distributed,
and it costs one comparison per shared node.

---

## 4. Proposed API

Extend the existing surface rather than adding a parallel one. The `target`
seam in ADR 0077 INV-7 is about *runtime* (tcl vs pymp); this is a different
axis — the *solver* — so it wants its own parameter:

```python
apeSees.modal_deck(
    path,
    *,
    band=(f_min, f_max),      # FEAST only
    num_modes=None,            # ARPACK only
    solver="feast",            # "feast" | "arpack"
    certify=False,             # FEAST only
    target="tcl",
    out="eigenvalues.out",
)
```

Rules, all of which should fail loudly rather than guess:

- `solver="feast"` requires `band` and rejects `num_modes`; `solver="arpack"`
  requires `num_modes` and rejects `band` and `certify`.
- `solver="arpack"` on a model with **no partitions** is an error — this route
  has no reason to exist unpartitioned, and silently emitting a one-rank deck
  would hide the mistake. Point the message at Tier 0.
- `solver="arpack"` with `target="pymp"` raises, same as FEAST's (the modern
  interpreter still builds its `ArpackSOE` bare — ADR-1000 §34.4; #668 is
  Tcl-only).

Emitter method: `TclEmitter.eigen_parallel(...)`, sibling to
`eigen_feast_parallel`. Keep the docstring discipline of the existing one — it
is the reason the FEAST path is auditable a month later.

---

## 5. Traps — do not rediscover these

Every item here was measured or source-verified during #668 and its adversarial
review. Full detail in ADR-1000 §34.

- **Nodal mass on shared boundary nodes double-counts, silently.** The `M*v`
  merge sums per-rank contributions, so a boundary node given mass on both
  ranks contributes twice. Because the Tcl `eigen` path always has `shift = 0`,
  the `addM` side-channel quick-returns and **K stays exactly right while M is
  wrong** — so nothing ever surfaces as an error. You get a plausible spectrum,
  biased low. **apeGmsh should own this**: it already knows partition
  membership, so assign each node's mass to exactly one owning rank at emission
  (the same owner rule already used for nodal loads) rather than trusting the
  user.
- **`modalDamping` is silently wrong on a partitioned deck** — the modal
  projection is a rank-local partial sum, so damped response changes with rank
  count, with no warning. apeGmsh already refuses `ops.damping.modal` under MPI
  (`apesees.py`); **keep that guard** — this route makes the preceding `eigen`
  succeed, which removes the accidental protection that a failing eigen used to
  provide.
- **A partitioned deck cannot switch to a serial system.** With `numberer
  ParallelPlain`, emitting `system UmfPack` fails `LinearSOE::setSize()` inside
  `domainChanged` *before* `eigen` runs. Relevant if you ever emit a
  multi-analysis deck that mixes modal and static stages.
- **Exit codes lie under MPI.** A Tcl error exits **0** (`g3TclMain` discards
  its exit code; `mpiParameterMain` returns 0 unconditionally), and a rank that
  `exit`s while a peer is deadlocked blocks forever in `MPI_Finalize`. If job
  orchestration keys success off the process exit code, it will report PASS on
  an aborted run. Key off harvested artifacts instead, under a timeout.
- **No feature probe exists.** Wrong binary ⇒ hang or empty spectrum, not a
  clean error.

---

## 6. Acceptance criteria

What we will check this against:

1. **Deck text** — partitioned blocks present (contrast the FEAST deck's flat
   emit); forced `Transformation`; `ParallelPlain`; `Mumps`; exactly one
   `[eigen ...]` capture; no `modalProperties`; per-rank shape recorders and
   sidecars.
2. **Numerical** — on a model small enough to also run serially: distributed
   spectrum == single-process `eigen` oracle to solver tolerance, at np = 1, 2
   and 4. The fork's own smoke, `Ladruno_scripts/verify_arpack_mp_mumps.tcl`,
   is the reference pattern (analytic fixed-free chain + mode-shape ratios).
3. **Mode-shape merge** — `mode_shape_field(k)` returns the **full** field, and
   the boundary-node agreement check of §3 is exercised by a test with a real
   shared boundary.
4. **Mutation** — the discipline that caught the weak test in #668, and worth
   repeating here: break the emitter deliberately (emit `system UmfPack`
   instead of `Mumps`, or a flat deck instead of partitioned) and confirm the
   test suite fails. A test that passes both ways is not a test. If a
   mutation cannot be made to fail, say so rather than papering over it.
5. **Guardrails** — `solver="arpack"` without partitions raises; with `band`
   raises; with `target="pymp"` raises.

---

## 7. Related apeGmsh work this implies

- **ADR 0077 needs amending, and it is the prerequisite.** It currently lists
  this route under *Rejected alternatives* — "REFUTED (F1) … never emitted"
  (lines 118 and 244). That was correct against vanilla and is now false
  against a #668 build. The F1 appendix also describes the reduction as gated
  by `#ifdef _PARALLEL_PROCESSING`; the operative gate is the **runtime**
  `processID`, which is why the class's own collectives were present but
  dormant. Amend both, and add the backend to the Tier-1 decision.
- **Tier 0 stays the default.** Nothing here displaces it: serial ARPACK is
  still the fastest path at every size measured, and it remains the **only**
  route to correct participation factors and effective modal mass. This
  backend is for the model that does not fit, and the docstring should say so
  in those words.

---

## 8. Status of the fork side

| item | state |
|---|---|
| `ArpackSOE::setProcessID` / `setChannels` | merged, PR #668 |
| MP `eigen` wiring, creation + reuse paths | merged, gated on `MumpsParallelSOE` |
| `verify_arpack_mp_mumps.tcl` (4 checks) | merged; serial/np2/np4 ≤1.96e-15 |
| mutation-tested against a gate-deleted build | yes — old smoke false-passed, current one fails |
| openseespy / PyMP | **not fixed** — same F1 defect, out of scope |
| `modalProperties` under MPI | **still blind** — upstream, unchanged |
