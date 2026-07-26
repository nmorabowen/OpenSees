---
title: ADR-75 session handoff — sparse-direct solver lane (PARDISO desktop / MUMPS cluster)
project: Ladruno
status: CLOSED for this session — P0/P1/P1c/P1d/P1e/P1f/P1g/P2/P2b + ADR-75b ALL MERGED to ladruno (tip `103d35e`). Desktop solver lane is DONE. Remaining work is cluster-gated or Lane-3; see "What to do next". Pick up from here in a new session.
owner: nmora
relates: 75_ladruno_sparse_direct_strategy_adr
---

# ADR-75 handoff (session of 2026-07-24/25)

Everything below is **merged to `ladruno`**. Pick up from "What to do next".

## What shipped

| PR | Merge | What |
|---|---|---|
| [#617](https://github.com/nmorabowen/OpenSees/pull/617) | `e8fda0e` | ADR-75 scoping (3 lanes + Abaqus/Kratos/LS-DYNA precedent) |
| [#618](https://github.com/nmorabowen/OpenSees/pull/618) | `57c5b52` | Adversarial-review corrections to the ADR |
| [#619](https://github.com/nmorabowen/OpenSees/pull/619) | `3d36da5` | **P0** trade study (portfolio confirmed) + Lane-B bench harness |
| [#620](https://github.com/nmorabowen/OpenSees/pull/620) | `434b841` | Lane-B baseline measured (UmfPack 22.711 s = P1 gate) |
| [#621](https://github.com/nmorabowen/OpenSees/pull/621) | `c9494db` | Clarify the baseline's MUMPS row is serial-only |
| [#622](https://github.com/nmorabowen/OpenSees/pull/622) | `f5c43db` | **P1a** PARDISO compiles + factorization reuse |
| [#623](https://github.com/nmorabowen/OpenSees/pull/623) | `4144a65` | **P1b** wire `system Pardiso`, threaded MKL — gate PASSED |
| [#624](https://github.com/nmorabowen/OpenSees/pull/624) | `6921f93` | **P1c** scaling study + capability wall |
| [#625](https://github.com/nmorabowen/OpenSees/pull/625) | `2bb348c` | **P2** MUMPS `-BLR` + oneAPI-2026 toolchain fixes |
| [#626](https://github.com/nmorabowen/OpenSees/pull/626) | `434b1a6` | **P2b** `-stats` (INFOG/RINFOG) + BLR memory verdict |
| [#630](https://github.com/nmorabowen/OpenSees/pull/630) | `222dca8` | **P1d** symmetric PARDISO `-matrixType` + `-stats` + Tcl verb (+ adversarial review) |
| [#633](https://github.com/nmorabowen/OpenSees/pull/633) | `b5b7dc1` | **P1e** `-krylov` factorization-preconditioned CGS (`iparm[3]`) — the full-Newton reuse lever (+ adversarial review) |
| [#635](https://github.com/nmorabowen/OpenSees/pull/635) | `1f5cde8` | **P1e follow-up** softening/limit-point sweep — CGS never falls back; post-peak branch caveat |
| [#636](https://github.com/nmorabowen/OpenSees/pull/636) | `580cfb5` | **P1f** `addA` binary search (1.09-1.10x) + **threaded PARDISO is NOT byte-reproducible** (corrects P1) + apeGmsh recipe |
| [#637](https://github.com/nmorabowen/OpenSees/pull/637) | `24a0535` | docs — record P1e-followup + P1f merge commits |
| [#632](https://github.com/nmorabowen/OpenSees/pull/632) | `8f3d1a4` | **ADR-75b — Lane 3 OPENED**: determinism CI policy + scatter remedy SETTLED, stage L3-0 measured, **3 adversarial passes** (the 3rd a 64-agent workflow that reversed the staging recommendation) |
| [#638](https://github.com/nmorabowen/OpenSees/pull/638) | `103d35e` | **P1g** UmfPack `addA` binary search — the other half of P1f's finding |

## Headline results (all measured, not estimated)

**Desktop PARDISO (`system Pardiso`) — the big win.**
- Lane B, same binary: **1.71× UmfPack @4 threads**, 1.76× @8, **1.19× even single-threaded**.
- **The win COMPOUNDS with size: 1.61× (11.5k DOF) → 2.15× (26k) → 3.40× (51k)** — UmfPack scales
  ~O(N²), PARDISO ~O(N^1.45).
- **Capability, not just speed: UmfPack OOM'd at 86,490 DOF; PARDISO solved it in 30.4 s and
  136,080 DOF in 68.6 s.** Largest single-machine model raised ≥1.6× in DOF; true ceiling untested.
- **Bit-identical to UmfPack at every thread count** (accuracy). ⚠ **But NOT byte-reproducible
  run-to-run when threaded** — corrected by P1f ([#636](https://github.com/nmorabowen/OpenSees/pull/636)):
  10 runs of one binary at 4 threads give **2 distinct answers (5/5 split, ~1 ULP)**; 1 thread is
  10/10 identical. The original "no FP-determinism problem in the solver lane" came from one run per
  thread count. **Pin `MKL_NUM_THREADS=1` for any byte-identical gate** (or wire `iparm[33]` CNR).
- **4 threads is the practical default** (1→8 only buys 1.58×; memory-bandwidth-bound).

**Symmetric PARDISO (`-matrixType 2`) — the memory lever BLR was supposed to be.** (P1d, new)
- **1.94-1.96x UmfPack** @4 threads (1.62-1.63x @1 thread), reproducible to +-1% over two sweeps.
  Versus *unsymmetric PARDISO* it is **~1.25x** (1.24-1.35x): that ratio is noisier because the
  unsym PARDISO row is the least reproducible measurement in the set, not the symmetric one.
- **Peak memory -41.8%** (105.60 -> 61.48 MB), stored nnz -49.3%, factor nnz -47.3%. `-matrixType 1`
  (SPD) is -46.4% peak but **fails on any indefinite tangent** - prefer 2.
- **Exact**, not approximate: bit-identical answers, so unlike BLR it is legal on oracle paths.
- **Refutes the `SparseSYM`-based worry** that symmetric ≠ better — that was SparseSYM's
  implementation, not symmetric storage.
- **Default stays unsymmetric on purpose** (half-storage silently solves the wrong system on a
  genuinely unsymmetric tangent — contact / non-associated / follower loads / `LadrunoUP`).

**MUMPS BLR (`-BLR <eps>`) — shipped, but did NOT pay at any size reachable here.**
- At ~32k DOF/np2: `eps=1e-4` → factor entries −21.8%, BLR flops −45%, **but peak MB/proc only −4.6%**,
  and **3.17× slower**. `eps=1e-9` is strictly worse (+8.4% MB/proc, 1.70× slower).
- **Peak memory is dominated by the active frontal/working space, not stored factors** — so "BLR saves
  memory" is largely false of the number that decides whether a model fits (`INFOG(21)`).
- Nonlinear trap: an inexact solve degrades Newton ⇒ more iterations ⇒ looser eps can be *slower*.
- This bounds the **small end only**. BLR's payoff is governed by **front size** (≈ nested-dissection
  top separator, ~N^(2/3) DOF; front memory ~N^(4/3)) — i.e. by **model size, not rank count**.
  ~20 MB fronts here vs ~1.7 GB at 1M DOF, ~36 GB at 10M DOF.

## Build state (IMPORTANT)

- This worktree's `dist/` has **both** rebuilt modules: serial `dist/bin/opensees.pyd` **with PARDISO**
  (477 pardiso symbols, was 0) and `dist/openseesmp/openseesmp.pyd` **with `-BLR`/`-stats`**.
- **The SHARED checkout `C:\Users\nmb\Documents\Github\OpenSees\dist\` still has the OLD binaries** —
  rebuild there before anyone expects `system Pardiso` outside this worktree.
- ⚠ **oneAPI 2026.1 ships NO Fortran** (no `ifx.exe`, no `ifconsol.lib`). An update re-pointed
  `compiler\latest` 2025.1 → 2026.1 mid-session and broke MP builds three ways. **All three fixed** in
  `Ladruno_scripts/setup_env.bat` + `build.bat`; written up as **`BUILD_GOTCHAS.md §8`**. Serial builds
  are unaffected, which is why it hides until an MP build.

## ⇢ START HERE (new session, 2026-07-25 close-out)

**Everything below the line is merged.** `ladruno` tip = `103d35e`. The **desktop solver lane is
DONE** — `system Pardiso` with factorization reuse, symmetric half-storage (`-matrixType 0|1|2`),
`-krylov` CGS reuse, `-stats`, and now a binary-search `addA` in both desktop SOEs. Nothing in
ADR-75's P-series is outstanding.

**The four things left, and what gates each:**

> **⚠ SUPERSEDED IN PART — see "§ 2026-07-25 (b) session" below.** **The cluster is UP** (it was not
> re-checked last session, only assumed down). **Item 4 is ANSWERED** (park Lane 3, ADR-75b §12).
> **Item 1 was never "~10 min, no code"** — the Tcl `system Mumps` ladder could not parse `-BLR` or
> `-stats` at all, and dropped them *silently*; fixed as **P2h**, but the cluster binary still needs
> a rebuild before the measurement is possible. The table below is kept as written for provenance.

| # | Item | Gated on | Effort |
|---|---|---|---|
| **1** | **Cluster BLR validation** — production deck, `-stats`, BLR on/off, compare `INFOG(21)` | **the cluster being up** (it was down all session) | ~10 min, no code |
| **2** | **MUMPS `-matrixType 2`** on the cluster — same symmetric lever P1d measured at −41.8% peak memory on the desktop; untouched on MUMPS, strong prior | the cluster | small |
| **3** | **Rebuild the SHARED checkout** so `system Pardiso` works outside a worktree | **owner decision** — an agent must not switch the shared checkout's branch under the other ~7 worktrees ([[ladruno-build-in-worktree-not-shared-checkout]]) | one build |
| **4** | **Lane 3 (ADR-75b)** — threaded element assembly | see below; **genuinely open, not merely unstarted** | large |

## § 2026-07-25 (b) session — THE CLUSTER IS UP; item 4 answered; item 1's premise was wrong

### 1. The cluster is UP. Check it, do not inherit "it was down".
`ssh esmeralda` answers on the first try, **BatchMode, no prompt**. Verified 2026-07-25:

| fact | value |
|---|---|
| host | `esmeralda` — `~/.ssh/config` → `200.55.208.141:6023`, user `nmorabowen` |
| uptime | **33 days** — so it was *also* up during the session that recorded it as down |
| login node | 24 cores / 31 GB / 749 GB free on `/mnt/deadmanschest` |
| scheduler | **SLURM at `/opt/slurm/bin`, NOT on `PATH`** — `which sbatch` fails; `/opt/slurm/bin/sinfo` works |
| nodes | 18 in `computes`: **3 idle** (`node[1-2,11]`), 12 `mix`, 3 `down*` |
| queue | busy with another user's jobs (`pxpalaci`, many `PD`) — expect to wait, or target the idle nodes |
| github | the login node reaches `github.com`; `git fetch origin ladruno` pulled `d4f68cf4a` fine |

**The likely reason it read as "down" before: `which sbatch` returns nothing.** That is not a dead
cluster, it is a `PATH` that omits `/opt/slurm/bin`. Build/run guide: `Ladruno_internal/02_esmeralda_linux_build_guide.md`.

### 2. Item 1 was NOT "~10 min, no code" — it was unmeasurable, and would have returned a FALSE NEGATIVE
`OpenSeesMP` is **Tcl**-driven, and the Tcl `system Mumps` ladder (`SRC/tcl/commands.cpp`) carried
**only** `-ICNTL14`/`-ICNTL7`/`-matrixType`. `-BLR` and `-stats` were wired into
`SRC/interpreter/OpenSeesCommands.cpp` (Python) **only**, and the Tcl ladder's `else currentArg++`
arm discarded unknown tokens **silently**. So the exact command the handoff prescribed —
`system Mumps -ICNTL14 200 -stats` vs `... -BLR 1e-8 -stats` — would have run **full-rank both
times, with no stats**, and "`INFOG(21)` did not move" would have been written down as the answer.

- ✅ **Fixed as P2h (this PR):** `-BLR`/`-stats`/`-ICNTL35`/`-CNTL7` in the Tcl ladder, 5-arg ctor on
  both parallel arms, a missing-value bounds check (it also read `argv[currentArg+1]` past the end),
  and the silent skip replaced by a warning. **Compile-verified on all three `#ifdef` arms with a
  negative control; NOT executed** — no Tcl+MUMPS binary was reachable (see 3).
- **`-matrixType` (item 2) needs NO code** — it is vanilla and already present in *both* ladders.
  Item 2 really is just a run, once a current binary exists.

### 3. The cluster binary is 460 commits stale — this is now item 1's real gate
`~/ladruno_build_test/OpenSees` was at **`e503ce4c0` (#580, built 18 Jul)** — predates `-BLR`
(#625) entirely. Its tree is **clean** (no local commits, no modified tracked files; only untracked
`feast_l2_profile/*.json` artifacts), so updating it is safe. It has since been fetched to
`origin/ladruno` = `d4f68cf4a`; **the working tree was left untouched at `e503ce4c0`.**

⚠ **That checkout is on the local branch `ladruno-p5-build`, NOT `ladruno`** — it is ADR-74's cluster
tree ("the Esmeralda tree (`ladruno-p5-build` @ #580)", [[74_ladruno_parallel_numberer_adr]] §Cluster
note), even though `02_esmeralda_linux_build_guide.md` §2 still describes the clone as being on
`ladruno`. The staged script therefore uses `git checkout --detach <ref>` on purpose: it moves the
build tree without rewriting or clobbering `ladruno-p5-build`, so ADR-74's reference point survives.
**Do not `git reset --hard origin/ladruno` on that branch** (which is what the build guide §6 says)
unless you have first checked that ADR-74 is done with it.

**A rebuild script is staged and ready at `esmeralda:~/adr75_rebuild.sh <ref>`** — fetch + detached
checkout + `conan install` + the proven MP cmake line (`/opt/openmpi` wrappers, `/mnt/nfshare` MUMPS
+ ScaLAPACK) + `--target OpenSeesMP -j20`, logging to `~/ladruno_build_test/adr75_build.log` and
ending with the `ldd` link check. **It was NOT run:** this environment's tool policy blocked every
`ssh` that builds or executes on esmeralda (read-only inspection was fine), and the owner elected to
**park the cluster lane** rather than change the permission this session.

**So item 1/2 now read:** rebuild esmeralda at a ref containing P2h → run the two-way (and
`-matrixType 0` vs `2`) comparison → compare `INFOG(21)`, total wall, and the answer vs full-rank.
**Capture ADR-75b's G-L3 fractions in the same run** (see 4) — it is the same deck and the same cost.

### 4. Item 4 (Lane 3) is ANSWERED: **park it.** — [[75b_ladruno_threaded_assembly_adr]] §12
`§11 q8` is closed. ADR-75b status → **`PARKED — measurement-gated (G-L3)`**, and **prerequisite
work is de-authorized too**, which is the part that changed. The argument in one line:

> **The reason Lane 3 was attractive no longer describes any target it ranks.** §2.1's contribution
> was that loop A has *no FP reduction*, so threading it is bit-identical and free of the determinism
> policy. After L3-0b shipped, loop A is **rank 3**; the two above it (loop C 46.33%, loop B 61.56%)
> are **both reducing loops**. What is left is the expensive version of the lane — P-a…P-e in full,
> the §3 ordered/fast policy, an OpenMP build, and (q9) a **1–4 GB serial** memory regression, to buy
> a **≤2.78×** 4T ceiling — against a gate measured at **11 520 DOF**, ~12× below the 136k DOF
> ADR-75 §1 already calls non-representative, in the direction ADR-75 §1 says will shrink it.

**G-L3, the one measurement that reopens it:** on a **≥500k DOF** production deck, profile with
**coarse** `OPS_PROF_FLAGS="-perStep"` (not `-deep`) and report **each individual assembly loop's %
of step** (loop A `Domain::update`, loop B `formTangent`, loop C `formUnbalance`) plus solve %.
**Largest SINGLE loop ≥40% ⇒ re-authorize Lane 3 at that loop. No loop reaches 40% ⇒ close Lane 3**
as Amdahl-irrelevant and send the residual back to the solver lanes. One run, no code.

> **⚠ The gate is on a SINGLE LOOP's fraction, not the aggregate kernel fraction** — the first draft
> of ADR-75b §12 got this wrong and it reversed the answer. Lane D reads **kernel 85.30% but
> gate FAIL**, because loop A was 38.95%. Also note **coarse mode cannot give the kernel-vs-scatter
> split inside a loop** (`elemByType_` is deep-mode only, `Profiler.h:149`); that is a separate,
> later measurement and belongs on a reduced-size deck, since `-deep` at production scale is both
> prohibitive and distortive. Both points are written up in ADR-75b §12.4.

### 5. Item 3 is unchanged and is the owner's
Confirmed with the owner this session: **nmora rebuilds the shared checkout, not an agent.** It sits
on `ladruno` at `2b06dba0f`; `d4f68cf4a` is a clean fast-forward (the ~14 attached worktrees are all
on their own `claude/*` or `up/*` branches, so no branch switch is involved) — but it stays the
owner's call. Left untouched.

---

**On item 4, read `75b_ladruno_threaded_assembly_adr.md` §7 before planning anything.** Its own
measurement now argues *against* starting where the first draft said to:
- **Work removal outranks threading.** L3-0b(i) `-commitSolveState` on explicit decks is **−29.6%
  wall, bit-identical, shipped, zero code** — and enabling it drops Lane D's loop-A fraction from
  55.7% to **38.95%, i.e. FAILING the >40% gate**, while loop C rises to 46.3%. L3-0b(ii), the `addA`
  search, is now **shipped on both solvers** (P1f + P1g) and bought 1.03–1.10×.
- **So no lane currently passes the gate for loop A on a correctly-configured, fork-owned deck.**
  L3-1 is **de-authorized as the first threading stage**.
- **The prerequisites are heavier than first assumed** (§5.4, from the 64-agent review): the
  re-entrancy inventory was directory-scoped and missed `FE_Element::theMatrices` (a class-wide pool
  handed straight to `addA` — a **100% collision** no element allowlist can fix), `LadrunoBrick`'s 16
  `update()` statics, the **shared mutable loop iterator** (the loop is not even index-addressable),
  `Node`'s lazily-allocating trial *getters*, `SRC/coordTransformation/`, and the **profiler
  instrument itself**.
- **Honest question to answer first (§11 q8):** is Lane 3 still worth it after the cheap work-removal
  landed? Re-measure the gate *after* L3-0b, not before.

**Method rules this session paid for — apply them:**
1. **"Deterministic" is a claim about a DISTRIBUTION.** Threaded PARDISO gives 2 distinct answers over
   10 runs at 4 threads (1 thread: 10/10 identical). The original "bit-identical at every thread
   count" came from n=1 per thread count. **Pin `MKL_NUM_THREADS=1` for any byte-identical gate.**
2. **Check an invariant where you DEPEND on it**, not where it happens to be established today.
3. **Prefer running the omitted configuration to re-reading your argument.** Every real defect found
   across three adversarial passes was a *confounded comparison* or a *selectively scoped statistic* —
   never a logic error.
4. **A phase baseline in a dated report may have been optimized away by a later ADR.** This bit twice
   (ADR-67 P-NEW-1, then P-NEW-2 two bullets away in the same file).

---

## What to do next (ranked)

1. **Cluster validation of BLR — the open question, ~10 min, NO code needed.** Run a production deck
   twice and compare: `system Mumps -ICNTL14 200 -stats` vs `... -BLR 1e-8 -stats`. Compare
   **`INFOG(21)` (MB/proc)**, total wall (must include the whole Newton loop), and the answer vs
   full-rank. Start near `1e-8`. If `INFOG(21)` doesn't move, BLR is not the memory lever at your scale
   either and the next candidate is **MUMPS out-of-core `ICNTL(22)`** (trades memory for I/O, not
   accuracy — not yet exposed by the fork; small addition next to `-BLR`).
2. ~~**P1b — symmetric PARDISO (`mtype ±2`)** (briefed in [[_adr75_p1b_brief]]).~~
   ✅ **DONE 2026-07-25, shipped as P1d ([#630](https://github.com/nmorabowen/OpenSees/pull/630))** —
   `phase1/RESULTS_p1d_symmetric.md`. Won on **both** axes: 1.94-1.96x vs UmfPack @4T (~1.25x vs
   unsymmetric PARDISO) and **-41.8% peak memory**, bit-identical answers. The brief's "must be measured, not assumed" caution
   was right to insist on the measurement and wrong about the direction — `SparseSYM`'s 2.10×
   slowdown was *its* implementation quality, not a property of symmetric storage.
   **Follow-on, now carrying a strong prior: exercise MUMPS `-matrixType 2` on the cluster** — same
   lever, still untouched there, and it composes with item 1.
3. ~~**Lane 3 — threaded element assembly.**~~ ✅ **OPENED 2026-07-25 as
   [[75b_ladruno_threaded_assembly_adr]]; policy SETTLED and stage L3-0 MEASURED.**
   - **Determinism policy settled (§3):** *per loop class.* `Domain::update` has **no FP reduction at
     all**, so threading it is bit-identical and the existing oracles gate it unchanged. For the
     reducing loops (`addA`/`addB`), **ordered mode is the CI gate AND the default when threading is
     on; fast/atomic mode is opt-in, 1e-12-gated, and forbidden on oracle paths** — LS-DYNA's
     `ncpu=-N` *shape* with the default **inverted**, because this fork's QA is exact, not
     tolerance-based (and BLR already set that precedent).
   - **Scatter remedy settled (§4), not relitigated:** freeze-sparsity + atomic scatter. New checked
     fact — **OpenSees already freezes the sparsity** (`zeroA` never touches `colA`/`rowStartA`), so
     the whole race is one `A[k] += m(i,j)`. Ordered mode = **gather** over that same frozen CSR
     (race-free without atomics, exact), which is why de-statication must produce **per-element**
     buffers, not `thread_local` ones.
   - **L3-0 measured** (`Ladruno_files/testbed/perf/lane3/RESULTS_l3a_update_scope.md`): **the Lane-1
     solver win is what makes Lane 3 worth doing** — same Lane-B model, element-kernel fraction
     **35.8% (UmfPack) → 74.9% (`Pardiso -matrixType 2` @4T)**, solve 55.9% → 16.9%; the >40% gate
     fails under UmfPack and passes hugely under PARDISO. Lanes A/D element-bound regardless
     (81.9%/86.8%). `ForceBeamColumn2d::getTangent` = **0.12 µs/ele** (confirms the regression hazard).
   - **⚠ Next step REVERSED by a 64-agent max-effort review (3rd pass).** The draft said "next =
     L3-1, thread `Domain::update` on Lane D (54.4%)". That 54.4% is **half redundant work**: it sums
     two `elem.update` sites and ADR-67 P-NEW-2's **shipped `-commitSolveState`** deletes the second.
     Re-measured (3 interleaved A/B rounds, `disp_z` bit-identical in all 6): step 24 555 → 17 278 ms
     (**−29.6%**), **loop A 55.74% → 38.95% = FAILS the >40% gate**, loop C 32.77% → **46.33%**.
     ⇒ **L3-1 is de-authorized as the first threading stage; no lane passes the gate for loop A on a
     correctly-configured, fork-owned deck.**
   - **NEXT = L3-0b, work removal, before any threading** (ADR-40's standing order, which the ADR had
     inverted): (i) **`-commitSolveState` on explicit decks** — shipped, bit-identical, −29.6%;
     (ii) ~~replace `addA`'s O(idSize² × rowlen) search~~ ✅ **FULLY SHIPPED, both solvers:** PARDISO
     as P1f ([#636](https://github.com/nmorabowen/OpenSees/pull/636), 1.09-1.10×) and UmfPack as P1g
     ([#638](https://github.com/nmorabowen/OpenSees/pull/638), 1.03-1.06×). L3-0 named it and it
     shipped within a day, needing no re-entrancy audit, no determinism policy and no OpenMP build.
     **Note the asymmetry, it is the useful part:** the UmfPack win is *smaller and shrinks with
     size* where PARDISO's grows — `addA` pays in **inverse proportion to how expensive the solve
     already is**, and UmfPack scales ~O(N²). Nothing left on this item.
   - **⚠ P1f also corrected a claim ADR-75b's determinism policy rested on:** threaded PARDISO is
     **not byte-reproducible run-to-run** (4T: 2 distinct answers over 10 runs of one binary, ~1 ULP;
     1T: 10/10 identical). So the "existing byte-identical oracles gate it unchanged" story holds
     **only with `MKL_NUM_THREADS=1`**, and ADR-75b §7's acceptance protocol — which said
     "bit-identical at 1/2/4/8 threads", i.e. **n=1 per thread count** — had the same flaw that
     produced the original wrong claim. Both fixed (§3 P-6, §7 item 4): N≥10 repeats, each config must
     reproduce **itself**, and pin MKL to 1 thread.
   - **Then, before ANY threading code:** the review found the re-entrancy inventory was
     directory-scoped and missed the hazards that actually race — **`FE_Element::theMatrices` is a
     class-wide pool** (one Matrix shared by every same-numDOF FE_Element, handed straight to `addA`
     ⇒ a **100% collision** no element allowlist can fix), `LadrunoBrick::update()`'s 16 statics, the
     **shared mutable loop iterator** (so the loop is not even index-addressable), `Node`'s lazily
     allocating trial **getters**, `SRC/coordTransformation/`, and **the profiler instrument itself**.
     ADR-75b §5.4 lists all seven; §10's vanilla-file plan roughly tripled as a result.
   - **Also corrected:** the taxonomy had 3 loops; there are **five** — the transient path assembles
     in `TransientIntegrator.cpp`, and both paths carry **DOF_Group `addA`/`addB` loops** carrying the
     nodal mass and nodal loads, so a gather keyed on element matrices would silently drop them.
4. **Rebuild the shared checkout** so `system Pardiso` (+ `-krylov`, `-matrixType`, P1f/P1g `addA`)
   is available outside a worktree. **Now fully unblocked — every ADR-75 PR is merged.**
   **Owner: nmora, once P1d merges** — decided 2026-07-25 rather than have an agent switch the shared
   checkout's branch under the other 6 worktrees. One build from a clean `ladruno` gets P1b + P1d.
5. ~~**Version-control the perf skill.**~~ ✅ **DONE 2026-07-25 — but NOT in the fork.** The fork's
   `.gitignore` has a blanket `.claude/` rule ("session scratch — never commit"), so committing there
   needed a negation exception on an upstream file. Owner chose the pattern every other skill uses:
   its own private repo — `C:\Users\nmb\Documents\Github\Performance Skill\` →
   **`github.com/nmorabowen/opensees-performance-skills`** (manuals stay in Seafile). The install at
   `~/.claude/skills/opensees-performance/` remains; the repo is now the source of truth.
   ⚠ The global `~/.claude/CLAUDE.md` still says it lives in the fork — that line is stale.

## Gotchas banked this session (don't rediscover)

- `theSOE` is a **member of `OpenSeesCommands`** (`OpenSeesCommands.h:178`), not a file-scope global —
  a free `OPS_*` factory cannot assign it. The serial `#else` branch of `OPS_MumpsSolver()` still does,
  which is a **third proof serial MUMPS has never compiled** (with missing `libseq` + CMake gating).
- `~LinearSOE()` **deletes the solver**, so by the time a solver destructor runs the SOE has already
  freed `A/B/X/rowStartA/colA`. Never dereference `theSOE` in a solver destructor (cache what you need).
- The `system Mumps` option loop was `> 1` on the premise **every option takes a value**; a bare flag
  (`-stats`) at the end was silently dropped. Now `> 0`.
- `setup_env.bat` must **not** `setlocal` (it exports env to its caller) ⇒ **no delayed expansion
  `!VAR!`**; write those blocks flat with labels.
- CMake caches the **resolved absolute** compiler path; `-DCMAKE_Fortran_COMPILER=ifx` does **not**
  override a stale cached path. `build.bat` now deletes such a cache.
- A BLR run that gives a **bit-identical answer at a loose eps means BLR never engaged** — don't read it
  as a passing test.
- `sweep_p1.sh` originally `mv`'d the locked baseline JSON and destroyed it (twice). Now `cp`.

### Banked in the P1d session

- **A mistyped `system` option costs you the SOLVER, silently.** A factory that `return 0`s leaves
  `theSOE` null and OpenSees falls back to the default **`ProfileSPDLinSOE`** with only a warning —
  on a 3D solid mesh that is catastrophically slow, so the run looks *slow*, not *wrong*. Trigger
  here: openseespy `system('Pardiso','-matrixType','2')` — a **STRING** value fails
  `OPS_GetIntInput`. Cost 25 minutes of a dark bench. Two lessons: pass ints, and never let a
  bench script `grep` away the log that would have shown the warning.
- **MKL 2026.1 renamed its compute DLLs `.2.dll` → `.3.dll`** (scalapack/blacs stayed `.2`), so
  `build.bat`'s hardcoded `if exist`-guarded staging list copied *nothing* and the build still
  reported success — failing later as `ImportError: DLL load failed while importing opensees`.
  Fixed concurrently and more thoroughly by #627 (base-name globbing + stale-DLL purge); P1d only adds `mkl_avx10`. **`BUILD_GOTCHAS.md §9`**; same family as §8.
- **PARDISO symmetric needs its own `iparm`, not just `mtype`.** Scaling/matching
  (`iparm[10]`/`iparm[12]`) must be **0** for `mtype ±2` — MKL applies them as an *unsymmetric*
  permutation, which is the classic "symmetric PARDISO returns garbage" report. Use `iparm[9]=8`
  and `iparm[20]=1` (Bunch-Kaufman).
- **Derive, don't set.** `mtype` is read from the SOE's `matType` at the symbolic phase with no
  setter, so storage format and factorization mode cannot drift apart.
- `system Pardiso` was Python-only until P1d: `SRC/tcl/commands.cpp` has a **completely separate
  `system` if-ladder** from `SRC/interpreter/OpenSeesCommands.cpp`. Wiring one does not wire the
  other — check both when adding a `system` verb.

## Reproduce the benches
```bash
cd Ladruno_files/testbed/perf/phase1
./sweep_p1.sh                                    # UmfPack vs PARDISO, 1/2/4/8 threads
./sweep_p1d.sh                                   # + symmetric (-matrixType 1/2), 1/4 threads
OPS_PYD=<worktree>/dist/bin python3.12 -S p1d_memory.py    # -stats peak memory, sym vs unsym
<worktree>/dist/bin/OpenSees.exe p1d_tcl_smoke.tcl         # Tcl `system Pardiso` parity
OPS_SCALE_SIZES=15,20,25,30 MKL_NUM_THREADS=4 python3.12 -S laneB_scaling.py
BLR_NX=20 BLR_NY=20 BLR_NZ=24 ../../../../dist/openseesmp/mpiexec.exe -n 2 \
  python3.12 -S mp_blr_smoke.py ../../../../dist/openseesmp
```
Results: `RESULTS_laneB_baseline.md`, `RESULTS_p1_pardiso.md`, `RESULTS_p1c_scaling.md`,
`RESULTS_p2_blr.md`, `RESULTS_p1d_symmetric.md`.

⚠ Pass numeric option values as **ints**, not strings (`'-matrixType', 2`) — see the P1d gotchas.

## Standing context
The production regime is **huge solid nonlinear models** (not fiber frames). ADR-75 §1 carries the
correction: Lane B is the *primary* lane, the 66% solve fraction is a **floor** that grows with N, and
no gate measured at ≤136k DOF should be treated as production-representative.
