---
title: LadrunoCMS partitioned emitter — implementation request for the apeGmsh team
project: Ladruno
status: request-on-hold
priority: low
adr: ADR-1000
tags:
  - apegmsh
  - cms
  - modal
  - partitioning
  - mpi
  - emitter
---

# LadrunoCMS partitioned emitter — implementation request for the apeGmsh team

**What we are asking for:** an apeGmsh emission mode that writes a **physically
partitioned** OpenSees deck (one fragment per MPI rank) plus an auditable
`cms_partition_manifest.json`, for use with the fork's `eigen -ladrunoCMS`
hierarchical modal solver.

**Who this is for:** whoever owns emission in
[apeGmsh](https://github.com/nmorabowen/apeGmsh). Nothing here needs fork C++
knowledge. The normative source is ADR-1000 §17.5 ([[1000_ladruno_cms_adr]]);
this document is the apeGmsh-facing restatement, with the acceptance criteria
we will check against.

---

## 0. STAND DOWN — read this before you plan anything

**Correction, 2026-07-26, same day this document was written. Do not schedule
this work yet.**

An earlier version of this section told you to treat the emitter as a
**production dependency**. That was based on a decision to ship CMS to
production in exchange for its known slowness. **New measurements taken hours
later have undermined the reason CMS exists**, and the fork is now likely to
close the lane rather than ship it. You should not spend a sprint on a
dependency of something that may be shelved.

**What changed.** Measured on an identical mesh, standard solver vs CMS on four
ranks (both computing the same six eigenvalues to ~1e-12):

| mesh | n | standard | CMS, 4 ranks | ratio |
|---|---:|---:|---:|---:|
| 20x20 | 3 200 | 0.0294 s | 0.2335 s | 7.9x slower |
| 40x40 | 12 800 | 0.2061 s | 107.90 s | **523x slower** |
| 60x60 | 28 800 | 0.6480 s | did not finish | — |

CMS does not merely lose on speed, it **scales worse**: 4x the DOFs costs the
standard solver 7x and costs CMS 462x.

**Why that kills the rationale.** The justification for CMS was never speed — it
was *memory capacity*: it is the route to models too large for one node, because
no rank holds the whole problem. Two independent findings undercut that:

1. **The capacity route already exists in this codebase without CMS.** The
   standard `eigen` command drives shift-invert through whatever `LinearSOE` the
   analysis owns (`ArpackSolver` calls `theSOE->solve()`), `ArpackSOE` is
   parallel-aware, and `MumpsParallelSOE` is present and reachable from the MP
   interpreter. So `eigen` over a distributed MUMPS system already spreads the
   factorization — the dominant memory object — across ranks, using code that
   has been hardened for two decades.
2. **CMS's memory growth is asymptotically backwards.** Craig-Bampton
   materialises a DENSE constraint-mode block of size
   `n_interior x n_boundary`, which grows faster than the sparse factor fill it
   is supposed to avoid. Consistent with the real building model, where CMS used
   **more** total memory than the sequential solver, not less.

**What we are asking of you right now: nothing.** No emitter work, no manifest,
no `coarse_of_fine`. Keep the document for reference.

**One item is still worth having, and is cheap:** a **second valid 4-way
partition of Building 1A** (item 1 of §6). It is a data deliverable, not code,
and it unblocks a verification gate regardless of what happens to CMS.

**If CMS is revived** it will most likely be for a different purpose —
*superelements*: reducing a substructure once and reusing it across many
analyses (e.g. a foundation/soil block across a suite of ground motions), which
is what Craig-Bampton is actually good at. That use case would need much of this
same contract, so the document stays as written below.

---

## 1. Background in one paragraph

CMS (component mode synthesis) reduces a modal problem by splitting the model
into subdomains, reducing each one, and assembling a much smaller pencil. The
fork's implementation is two-level and MPI-parallel: `T2 → S2 → T1 → S1`, where
level 2 is the fine subdomain (one per MPI rank) and level 1 is a coarse group of
fine subdomains. In `physical` mode each rank builds **only its own part of the
model** — its own nodes and elements — and the global matrices are never
assembled anywhere. That is what a partitioned deck has to express.

---

## 1.5 What already exists on your side (verified against source)

Checked against the apeGmsh working tree at `334cc70e` before writing this, so
the ask is scoped to what is actually missing rather than to what we assumed.

**Already implemented — please reuse, do not rebuild:**

| Capability | Where |
|---|---|
| `g.mesh.partitioning.partition(n_parts, *, weights=None, backend=None) -> PartitionInfo` | `mesh/_mesh_partitioning.py`. Must run after `g.mesh.generation.generate()`. A `partition_explicit` variant also exists. |
| `build_element_partition_owner(fem)` → `{fem_eid: rank}` | `opensees/_internal/build.py` |
| `build_node_partition_owners(fem)` → incident ranks per node, plus `NodePartitionOwners.primary_owner` / `primary_owner_map` giving `{node_id: min(rank)}` | `opensees/_internal/build.py` |
| `ops.tcl(path, ..., per_rank=True)` — driver deck + one fragment per rank (your ADR 0061) | `opensees/apesees.py`. Note `split=True` and `per_rank=True` are mutually exclusive. |
| **Additive quantities already emitted once** — both `mass` lines and pattern `load` lines are bucketed by PRIMARY owner and emit on that rank only | `apesees.py` (`_bucket_mass_targets_by_rank`, `load_nodes=primary_nodes`) |
| Per-rank `fix` lines and owned-set caching | `apesees.py` |

That covers most of §2.1 already. **We are not asking you to redo any of it** —
we are asking you to make it *auditable*, *refusing*, and *CMS-aware*.

**Genuinely missing — this is the work:**

1. **`cms_partition_manifest.json`** (§2.2). Nothing like it exists today;
   `PartitionInfo` carries only `n_parts`, `elements_per_partition` and
   `weights_per_partition`, so it is not a starting point for the manifest.
2. **`coarse_of_fine`** (§2.2). CMS-specific, absent — nothing in apeGmsh has a
   notion of grouping partitions into coarse groups.
3. **The CMS preflight refusals** (§3). Validators exist for other purposes
   (e.g. `_run_staged_bc_validators`), but not this contract, and not the v1
   unsupported-construction list.
4. **A second valid 4-way Building-1A partition** (§6) — a data deliverable, not
   code.

If any of the "already implemented" rows is less complete than it looks from the
outside — particularly the additive-once guarantee under **staged** models — tell
us, because the fork side is currently assuming it holds.

---

## 2. What to emit

### 2.1 Per-rank deck fragments

One fragment per rank, plus a global driver. Your existing
`ops.tcl(..., per_rank=True)` path is the intended vehicle — per §1.5 it already
satisfies most of this table. The table is the **contract we will verify**, not a
list of things to write from scratch; treat it as the checklist to confirm the
existing path against.

Requirements per fragment:

| Requirement | Meaning |
|---|---|
| **Owned elements only** | Each element is emitted by **exactly one** rank. Never two, never zero. |
| **Nodes where needed** | A node is emitted on every rank incident to one of that rank's elements. Interface nodes therefore appear on ≥2 ranks — that is expected and required. |
| **Additive quantities once** | Nodal **mass** and **loads** emitted **only** on their single primary owner rank; summing across ranks must reproduce the monolithic model exactly. **Already implemented** (§1.5) — confirm it also holds for staged models. This is the easiest thing to get wrong and the hardest to see. |
| **Fixities / constraints where needed** | Emitted on the ranks that need them; restraint semantics must match the monolithic model. Building 1A has 1 333 base nodes and they must survive the split with the same restraints. |
| **Same global numbering** | Copies of one shared node must resolve to the **same global equation id** on every incident rank. |

### 2.2 The manifest

A single `cms_partition_manifest.json` next to the fragments. ADR-1000 §17.5
requires at minimum:

```jsonc
{
  "schema": "ladruno-cms-partition-manifest",
  "schema_version": 1,
  "femdata": {
    "hash": "<stable hash of the FEMData snapshot>",
    "identity": "<human-readable model identity/label>"
  },
  "global": { "num_nodes": 11841, "num_elements": 27360 },
  "num_partitions": 4,
  "coarse_of_fine": [0, 0, 1, 1],
  "ranks": [
    {
      "rank": 0,
      "num_elements": 6840,
      "elements_hash": "<hash of the SORTED owned element tags>",
      "num_nodes": 3102,
      "nodes_hash": "<hash of the SORTED emitted node tags>",
      "fragment_hash": "<hash of the emitted fragment file>"
    }
    // ... one entry per rank
  ],
  "additive_owners": {
    "mass": { "<nodeTag>": 0 },
    "loads": { "<nodeTag>": 0 }
  },
  "shared_nodes": { "<nodeTag>": [0, 1] }
}
```

Notes on the fields that are less obvious:

- **`coarse_of_fine`** — the CMS coarse grouping: for each fine subdomain (=rank)
  index, which coarse group it belongs to. Length = `num_partitions`. For 4 ranks
  with `p=2, m=2` this is `[0,0,1,1]`. **Group members should be topologically
  adjacent** — a coarse group is meant to be a contiguous chunk of structure, not
  an arbitrary set of ranks. Ranks in a group need not be contiguous *rank
  numbers* (the fork tests a non-contiguous grouping explicitly), but they should
  be neighbours in the model.
- **Hashes** — any stable algorithm, as long as it is documented in the manifest
  and reproducible. They exist so a rerun can be proven identical, and so we can
  tell "the deck changed" from "the solver changed".
- **`additive_owners`** — one owner rank per nodal mass and per nodal load. This
  is the audit trail for the "emitted once" rule.

The manifest is **generation evidence, not proof**. We verify the algebra
independently on the fork side (see §5); the manifest tells us *what was
intended* so that when the algebra disagrees we know which side is wrong.

---

## 3. Preflight: fail loud, before any MPI call

The emitter must **refuse to emit** — not warn — when it cannot honour the
contract. A wrong partitioned deck produces a plausible-looking wrong answer,
which is far worse than a refusal.

Refuse when:

- element coverage is not exact and unique (any element in two ranks, or in none);
- a node's incident-rank set is inconsistent with the elements actually emitted;
- an additive load or mass has no owner, or more than one;
- the model uses a construction whose partitioned semantics are not defined (below).

### 3.1 Explicitly out of scope for v1

ADR-1000 §17.5 excludes these from the first physical version. If a model uses
one, the emitter must reject it rather than emit something approximate:

- embedded reinforcement ties;
- auto structural rebar;
- `g.embed`;
- fork contact;
- zero-length / node-pair explicit elements.

Building 1A uses none of them. Any new model must pass the same preflight.

---

## 4. What the deck is driven with

For context only — this is the fork side, you do not implement it:

```text
eigen -ladrunoCMS -domainMode physical
      -hierarchy logical -level1 <p> -level2 <m>
      -modesL2 <k2> -modesL1 <k1>
      <numModes>
```

with `p * m == number of ranks`. `-verifyAssembly` is **not** valid in `physical`
mode (it is a parser error there, not a silently ignored option). The
non-degenerate reference case is 4 ranks with `p=2, m=2`; for 6 ranks we prefer
`p=2, m=3`.

---

## 5. How we will verify what you emit

So you know what "done" means. These are the P3a rows of the CMS P3 execution
plan ([[ladruno_cms_p3_execution_plan]]):

| Check | Assertion |
|---|---|
| Exact element coverage | union of per-rank element sets == monolithic set; count matches; zero missing or extra |
| Unique interior nodes | every interior node tag appears in exactly one rank |
| Interface nodes | a shared node appears on ≥2 ranks, all geometrically incident; never on a non-incident rank |
| Additive once | each nodal mass/load has exactly one primary owner; the sum over ranks equals the monolithic value |
| Deck equivalence | the monolithic-with-guards deck and the `per_rank=True` deck give the same eigenvalues and eigenvectors |
| **Neg: duplicate element** | a manifest with an element in two ranks → fail-loud preflight, no solve |
| **Neg: ownerless element** | an element in zero ranks → fail-loud |
| **Neg: inconsistent owner set** | a node whose incident ranks disagree → fail-loud |

The negative rows matter as much as the positive ones: we need to see the
emitter **refuse**, so please make the refusals testable (a distinct exception
type or error code is ideal).

Independently of the manifest, the fork already proves the algebra is a genuine
split — `tests/ladruno_cms_assembly_check.cpp` compares distributed `K·x` and
`M·x` against a serial reference to 1e-12 and has negative controls for a
replicated pencil and a double-owned element. Your manifest does not need to
prove correctness; it needs to make the *inputs* auditable.

---

## 6. The two things that would unblock us fastest

If you want the smallest useful first delivery, in priority order:

1. **A second valid 4-way partition of Building 1A** (any different, valid
   partition). The P3e invariance gate needs two partitions and only one exists.
   This alone unblocks a gate.
2. **The manifest + preflight** for the existing partitioning path — i.e. make
   the current capability auditable and refusing, before adding anything new.

Everything else (arbitrary `n_parts`, generalised grouping) can follow.

---

## 7. Open questions for your side

We do not want to over-specify your internals. Please push back on any of these:

1. Is `coarse_of_fine` something the emitter should **compute** (from mesh
   adjacency) or something the caller should **pass in**? We have no strong
   preference; we need it recorded in the manifest either way.
2. `ops.tcl(..., per_rank=True)` looks like the right vehicle to us (it already
   does the owner bucketing) — but would you rather the manifest ride that call,
   or come from a separate CMS-specific entry point that wraps it? We have no
   preference; you own the API shape.
3. Which hash do you want to standardise on, and over exactly what byte stream?
   We will read whatever you document.
4. Does the openseespy path need the same treatment, or is Tcl enough for now?
   The fork harness drives `OpenSeesMP` with a Tcl deck today.

---

## 8. Status and cross-references

- Fork-side CMS status, gate by gate: [[1000_ladruno_cms_adr]] (§17 defines the
  physical-distribution contract; §20–§25 record the 2026-07-26 session that
  closed Part 0, P3b and P3d).
- The remaining P3 matrix, including the P3a rows above:
  [[ladruno_cms_p3_execution_plan]].
- The performance picture and why the win — if it exists — is memory rather than
  speed: [[ladruno_cms_p4_execution_plan]].
- The apeGmsh-facing feature table: [[ladruno_apegmsh_contract]].

**Fork-side state at the time of writing (2026-07-26):** `LADRUNO_CMS=ON` builds
and all five standalone numerical checks pass at 2 and 4 ranks, plus two
end-to-end `eigen -ladrunoCMS` smokes. P3a and P3e are blocked on exactly the
work in this document.
