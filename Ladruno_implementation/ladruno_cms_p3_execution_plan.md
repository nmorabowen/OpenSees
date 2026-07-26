# LadrunoCMS P3 — physical-distribution verification: execution plan

**Status:** P3 NOT complete. One physical 4-rank Building 1A gate passed
(`rho_max=9.84323e-9`, eigenvalue error `4.47408e-12`, MAC `0.9999999999997988`;
see [[ladruno_cms_building_1A_acceptance]]). This document is the executable
plan for the remaining P3 matrix from ADR §17.8 ([[1000_ladruno_cms_adr]]) plus
the "remaining gates" in the acceptance report. It is a spec, not a claim of
completion — every item here requires **running** CMS across MPI ranks.

**Authored 2026-07-23 without a runtime build** (CMS is MP-only + `LADRUNO_CMS`
OFF by default; the authoring box could not produce a stable `LADRUNO_CMS=ON`
build). Expected values / assertions below are derived from the ADR contract and
the existing serial oracles, and MUST be confirmed on first real run — treat any
expected number as "assert to be pinned", not "known good".

## Prerequisites (the runtime environment P3 needs)

1. A working `LADRUNO_CMS=ON` MP build with the oneAPI MPI + MUMPS toolchain.
   Local: `set LADRUNO_CMS_BUILD=1 & Ladruno_scripts\build.bat` builds+runs the
   standalone checks (see [[../Ladruno_internal/BUILD_GOTCHAS]] #7). CI: the
   nightly self-hosted Zone-B lane (PR #612).
   **Status 2026-07-26: the local half is now MET** — two defects blocked it (a
   static/dynamic MKL link collision and a `build.bat` batch-parse bug, see
   [[../Ladruno_internal/BUILD_GOTCHAS]] #10); with both fixed all four checks
   pass under `mpiexec -n 4`. **The CI half is NOT met and never was:** the repo
   has zero self-hosted runners registered, so `zone-b-nightly` is cancelled on
   every scheduled run. Treat the local box as the only harness until a runner
   is registered ([[1000_ladruno_cms_adr]] §20).
2. The frozen Building-1A Gmsh harness (11 841 nodes, 27 360 elements) +
   `building_1A_cms_physical_run.ipynb` and its partition emitter.
3. A **second** valid 4-way Gmsh partition of Building 1A (for the P3e invariance
   gate) — currently only one 4-way partition exists.

Extension points already in-tree: `tests/ladruno_cms_{mumps,lanczos,hierarchy,`
`subspace,topology}_check.cpp` (standalone, mpiexec), `tests/`
`ladruno_cms_physical_smoke.tcl`, `tests/ladruno_cms_openseesmp_smoke.{tcl,py}`,
`tests/_testbed/ladruno_cms_reference.py` (numpy serial oracle).

---

## Part 0 — the 2-rank MUMPS ordering failure (blocks 2-rank support)

> **RESOLVED 2026-07-26 — but not by the fix proposed below, and it was never a
> 2-rank problem.** The 2-rank fixture this part asked for was built
> (`testDistributedMumpsAtScale` / `testSerialMumpsAtScale` in
> `tests/ladruno_cms_mumps_check.cpp`) and reproduces the exact failure with a
> DENSE order-12000 pencil. Measured: it fails at **np=4 as well as np=2**; a
> *sparse* order-64000 pencil passes at both; and **`ICNTL(28)=1` does not help**
> — so the "MUMPS chose parallel analysis" diagnosis below is wrong. The trigger
> is the dense pattern: `ICNTL(7)=7` selects PORD, whose multisector ordering
> cannot form stages on a dense matrix. Fix shipped: **candidate 2,
> `ICNTL(7)=0` (AMD)**, applied at BOTH factorization sites — including the
> serial `MumpsSPD` (`MPI_COMM_SELF`) path that the rank-local Craig-Bampton
> pencils actually use, which was outside this part's original scope and fails
> identically on one rank. Full evidence: [[1000_ladruno_cms_adr]] §21. The
> analysis below is kept as the historical hypothesis.

**Symptom** (scalability/RESULTS.md decision 3): 2-rank physical CMS fails in the
MUMPS analysis phase after ~124 s having consumed **15.6 GiB**, with
`orderMinPriority: no valid number of stages in multisector`. np=4 and np=6 work.

**Diagnosis.** `LadrunoCMSMumps.cpp:302-306` (distributed factorize) sets
`ICNTL(7)=7` (auto ordering) and `ICNTL(18)=3` (distributed assembled input) but
leaves **`ICNTL(28)` at its default 0 (auto), permitting MUMPS to pick *parallel*
analysis**. "multisector" is the parallel/PORD multisector ordering; it fails to
form valid stages for the np=2 process/matrix geometry, and the 15.6 GiB spent in
analysis (not factorization) is the signature of the parallel analysis path.

**Fix (ranked — validate at BOTH np=2 AND np=4; the gate is "np=2 fixed AND np=4
not regressed"):**

1. **`ICNTL(28)=1` (force sequential analysis)** — set `impl_->control.icntl[27]
   = 1;` in the distributed block (~line 306). MUMPS gathers structure to the
   host and orders sequentially, skipping the parallel multisector path. Fully
   compatible with `ICNTL(18)=3`; negligible cost for CMS's small reduced pencils
   (`r_D <= denseMax`). **Primary candidate.**
2. If sequential PORD still trips multisector: pin `ICNTL(7)=0` (AMD) — no
   multisector concept, always built into MUMPS.
3. `ICNTL(7)=5` (METIS) ONLY if the MUMPS build compiled the METIS ordering
   interface (the fork links its own METIS to the *library*, not necessarily into
   MUMPS — verify before relying on this).

**Validation protocol:** run `ladruno_cms_hierarchy_check` and a small physical
deck at np=2 (must now pass) and np=4 (must match current eigenvalues/residuals
byte-for-byte, or to 1e-12). Only then lift the "2-rank unsupported" note in
RESULTS.md. Do NOT ship the ICNTL change without the np=4 non-regression run.

> **Correction 2026-07-26 — the np=2 half of that protocol is currently vacuous.**
> `testDistributedMumps` ([ladruno_cms_mumps_check.cpp:89](../tests/ladruno_cms_mumps_check.cpp))
> and `checkDistributedFourRankFlow` ([ladruno_cms_hierarchy_check.cpp:290](../tests/ladruno_cms_hierarchy_check.cpp))
> both `return` immediately when `size != 4`. At np=2 the distributed
> factorization — the exact path that fails — is never entered, so the check
> passes without exercising anything. **Prerequisite for Part 0: a 2-rank
> distributed fixture (or a physical deck run at np=2).** Until one exists there
> is no way to satisfy this gate, so the `ICNTL(28)` change stays unapplied.
> See [[1000_ladruno_cms_adr]] §20.5.

---

## Part P3a — emitter + partition manifest

Harness: extend `topology_check` (has `OpenSeesLIB`) + a Python emitter test on
`cms_partition_manifest.json`.

| Check | Assertion |
|---|---|
| Exact element coverage | union of per-rank element sets == monolithic set; count == 27 360 (Building 1A); zero missing/extra |
| Unique interior nodes | every interior node tag appears in exactly one rank |
| Interface nodes on incident ranks only | a shared node's tag appears in ≥2 ranks, all of which are geometrically incident; never in a non-incident rank |
| Additive loads/masses once | each nodal mass / load has exactly one primary owner rank; sum over ranks == monolithic |
| **Neg: duplicate element** | manifest with an element in two ranks → fail-loud preflight, no solve |
| **Neg: ownerless element** | element in zero ranks → fail-loud |
| **Neg: inconsistent owner set** | a node whose incident ranks disagree on ownership → fail-loud |
| Deck equivalence | monolithic-with-guards deck == `per_rank=True` deck (same eigenvalues/vectors) |
| Preflight for unsupported | every not-yet-supported construction (non-plain load, unsupported constraint) rejected loudly before any MPI collective |

## Part P3b — distributed numbering + assembly (the Kx/Mx oracle)

Harness: new `ladruno_cms_assembly_check.cpp` (mpiexec -n 4), + reuse
`reference.py` serial matrices.

> **LARGELY DONE 2026-07-26** — `tests/ladruno_cms_assembly_check.cpp` exists and
> passes at np=2, np=3 and np=4 (size-generic by design; it announces the skip
> loudly at np=1). Covered: sparse local IDs, contribution locality, shared
> equations exist and agree across incident ranks, the **Kx and Mx oracles** (3
> deterministic + 3 fixed-seed random probes, 1e-12), an entrywise sum check, and
> both negative controls. **Not covered, and still open:** (a) the collective
> global-dimension row — it needs a real `Domain`, not a fixture; (b) the negative
> controls prove *the oracle* detects a double count, **not** that the production
> `-verifyAssembly signature|full` guard fails loud in the `EigenSolver` path,
> which also needs a `Domain`. See [[1000_ladruno_cms_adr]] §22.

| Check | Assertion |
|---|---|
| Shared-node global IDs | copies of one shared node carry the SAME global equation id on all incident ranks |
| Sparse local IDs | local ids are non-contiguous (a strict subset of `[0,n)`) |
| Collective global dimension | `n` obtained by collective, equals the monolithic DOF count (63 048 for Building 1A) |
| Contribution locality | each rank's captured `addA`/`addM` records == ONLY its owned elements' contributions |
| **Kx oracle** | distributed `K·x` == explicit serial `K·x` for ≥3 deterministic vectors (e1, ones, ramp) AND ≥3 fixed-seed random vectors, to 1e-12 |
| **Mx oracle** | same for `M·x` |
| **Neg: replicated input** | a rank fed the full (replicated) matrix instead of its split → detected + fail-loud (double-count guard) |
| **Neg: corrupt owner set** | a perturbed owner set → assembly signature mismatch, fail-loud, no silent wrong answer |

This is the single most important P3 gate: it proves "physically distributed"
means the algebra is a genuine split, not a replicated matrix.

## Part P3c — small physical CMS (correctness)

Harness: extend `hierarchy_check` + a small physical deck (e.g. the ZCB80-class
fixture), mpiexec -n 4.

- `physical` vs `replicatedReference` vs serial LAPACK — same spectra.
- Full basis → round-off residual; truncated basis → original-pencil residual + MAC.
- Lumped AND consistent mass, including `M_IB != 0` (coupling block preserved).
- Massless coordinates + static reconstruction (the `Z` condensation path).
- Repeated solves + `domainChanged` between solves (no state leak; identical result).
- Collective failure injected on one rank → all ranks return together, **no deadlock**.

## Part P3d — non-degenerate hierarchy (p=2, m=2)

Harness: `hierarchy_check` fixtures, mpiexec -n 4.

- 4 ranks, `p=2, m=2` (the mandatory non-degenerate case).
- Three interface topologies: level-2-only, level-1-only, combined.
- **Rank/partition-permutation invariance**: permuting which rank owns which
  partition leaves eigenvalues, residuals, AND the subspace (MAC≈1) unchanged.
- Level ablation (omit T1) allowed ONLY as a labelled diagnostic, never an
  accepted config.

## Part P3e — Building 1A physically distributed (the real gate)

Harness: `building_1A_cms_physical_run.ipynb` + a SECOND 4-way partition.

- One rank fragment per partition emitted from the same `FEMData`.
- Zero duplicated/missing elements; 1 333 base nodes preserved with compatible
  restraint semantics.
- The same 8 modes as P2c within eigenvalue / residual / MAC tolerance.
- `rho_max <= 1e-8` on the original pencil.
- **Invariance across ≥2 valid 4-way partitions** (the new requirement).
- 2 repeated runs (reproducibility).
- Per-phase time + peak-RSS-per-rank report — **without** yet claiming speedup
  (that is P4).

---

## P3 exit gate

Zero double-count, zero cross-rank divergence, P2c reproduced under genuine
physical distribution, 2-rank supported, and every negative control fails loud.
Only then does P3 close and P4 (performance regime + the dense→sparse workspace
work, see the P4 task) becomes the remaining gate to `shipped`. The nightly
Zone-B lane is the standing runtime harness for all of the above.
