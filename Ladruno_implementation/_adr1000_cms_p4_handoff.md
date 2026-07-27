---
title: ADR-1000 CMS — P4 session handoff (why CMS scales worse, and the two experiments that settle it)
project: Ladruno
status: handoff
priority: high
adr: ADR-1000
tags:
  - cms
  - modal
  - performance
  - p4
  - handoff
---

# ADR-1000 CMS — P4 handoff

**Read this first, then [[1000_ladruno_cms_adr]] §27–§30.** Everything below is
measured, not assumed; where something is an assumption it says so.

---

## 1. The one thing to know

CMS is **7.9x slower than the standard solver at n=3200 and 523x slower at
n=12800** on an identical mesh (§30). It does not merely lose on speed — **it
scales worse**. Quadrupling the DOFs costs Arpack+UmfPack 7x and costs CMS 462x.

That matters because the *accepted production rationale for CMS is memory
capacity* — the model too large for one node. If cost grows near `O(n^4)`, the
large model CMS exists to make possible would take an unusable amount of time.
A capacity win is worth nothing if the answer never arrives.

**Your job is to find out whether that is a tuning artefact or an algorithmic
property.** Two experiments, below. Do experiment A first: it is one run and it
discriminates.

---

## 1b. Scope note added after this handoff was written

**The final verdict is deferred to an Esmeralda campaign (ADR §31).** Your two
experiments still matter — they decide whether the 523x is tuning or algorithmic,
and that answer shapes how the cluster campaign is configured — but do NOT close
P4 on desktop evidence alone. All of it was taken at 4 ranks, which is where the
source paper's own data predicts CMS cannot win. §31 has the campaign design and
the exit criterion.

## 2. The hypothesis

CMS parameters were held **fixed** (`k2 = -modesL2 12`, `k1 = -modesL1 24`) at
every problem size in the §30 ladder. As the per-rank subdomain `m` grows, a
12-mode fixed-interface basis is relatively poorer, so:

- the **local fixed-interface Lanczos** (T2) needs more iterations and restarts —
  it was already seen exhausting `maximumRestarts` as `m` grew (§27.4, §29), and
- the **subspace refinement** has to correct a worse starting space.

So the suspicion is: *the T2 Lanczos cost grows super-linearly in `m` at constant
`k2`*, and that — not memory — is P4's real obstacle.

**If true**, the fix is a `k2` heuristic (scale the retained fixed-interface
modes with the subdomain, e.g. `k2 ~ sqrt(m)` or a fraction of `m`), plus
possibly an automatic default. **If false** — the cost is super-linear even with
`k2` scaled — then CMS's formulation is the problem at these sizes and P4 should
be closed honestly as "no winning regime exists", which the ADR explicitly
permits (§10, and the P4 plan's exit gate).

---

## 3. Experiment A — profile the 40x40 case (ONE run, do this first)

The 40x40 point is where the ratio explodes (107.90 s vs 0.2061 s). The T2 and
Lanczos instrumentation is already in place, so a single `-verbose` run tells you
where those 107 seconds go.

```bash
cd Ladruno_implementation/cms_profile
# env: see section 5 for the TCL_LIBRARY trap
set LADRUNO_CMS_PROFILE_WIDTH=40
set LADRUNO_CMS_PROFILE_HEIGHT=40
set LADRUNO_CMS_PROFILE_RESTARTS=4000
set LADRUNO_CMS_PROFILE_MAXITER=200000
mpiexec -n 4 <build>/OpenSeesMP.exe cms_profile_sheet.tcl
```

Read these lines:

- `LadrunoCMS timing:` — assembly / hierarchy / refinement / total
- `LadrunoCMS hierarchy phases [s]:` — partition / T2 / S2 / T1 / S1 / backSub / publication
- `LadrunoCMS T2 breakdown [s]:` — factorize / constraintModes / condensation / **lanczos** / reconstruct / scatter / congruence
- `LadrunoCMS T2 Lanczos [s]:` — rayleighRitz / orthonormalize / operator / residual **and the counts `(ritzCalls= opApplications= restarts=)`**
- `LadrunoCMS dimensions:` — `n`, `r2`, `r*`, `n/r2`, `r2/r*`

**What each outcome means:**

| Observation | Reading |
|---|---|
| `restarts` large, `opApplications` huge | Lanczos is thrashing — supports the hypothesis; go to experiment B |
| `refinement` >> `hierarchy` | the poor CMS starting space is the cost, not T2 itself — still supports the k2 story, but the fix is `k2`/`q`, not the Lanczos |
| `rayleighRitz` dominant with modest counts | the `j^2` projection is the cost — see §6 "known remaining inefficiency" |
| everything modest and total still ~108 s | something is wrong with the measurement, not the solver — re-check before believing it |

For reference, the **20x20** numbers (from §29, same deck) are:
`total 0.2461 | hierarchy 0.1194 | T2 0.0918 of which lanczos ~0.0753 |
ritzCalls 15 opApplications 120 restarts 1 | n=3200 r2=208 r*=88`.
Compare like-for-like: the interesting quantity is how `opApplications` and
`restarts` grow from 20x20 to 40x40 relative to `m`.

---

## 4. Experiment B — the ladder with `k2` scaled to the subdomain

Re-run the §30 comparison, but **scale `k2` with the subdomain** instead of
holding it at 12. Suggested rule to start from: `k2 = max(12, round(sqrt(m)))`,
`k1 = 2*k2`, where `m` is the per-rank interior count (the run prints it as
`interior m=` in the `T2 shape` line).

For the sheet at 4 strips, `m ≈ 2 * width * height - b`:

| mesh | approx m | `sqrt(m)` | suggested k2 |
|---|---:|---:|---:|
| 20x20 | 760 | 28 | 28 |
| 40x40 | 3 120 | 56 | 56 |
| 60x60 | 7 080 | 84 | 84 |

Use `cms_compare_arpack.tcl` (it builds the SAME mesh at np=1 and np=4, so the
comparison is honest) and add `-modesL2 / -modesL1` per size. **≥3 repeats per
point** — the P4 plan's statistical discipline, and this box has real variance.

**Acceptance / decision:**

- **Ratio flattens** (roughly constant CMS/standard as `n` grows) → it is tuning.
  Ship an automatic `k2` heuristic, re-run, and update §30. P4 may then be
  reachable on the memory axis.
- **Ratio still degrades** → algorithmic. Write it up honestly, and recommend
  closing P4 as "no winning regime at reachable scales" per the P4 exit gate.
  Do **not** keep optimising in the hope it turns around; that is what §26 and
  §29 already taught (see §7).

Note `denseMax` may need raising as `k2` grows: `r*` grows with `k2`, and the
final dense solve is capped by `-denseMax` (the deck passes 4000).

---

## 5. Environment — the traps that will cost you an hour each

1. **`TCL_LIBRARY`.** A freshly built `OpenSeesMP.exe` run straight out of the
   build tree starts with no commands (`invalid command name "wipe"`). Point it
   at the conan Tcl:
   `set TCL_LIBRARY=C:\Users\nmb\.conan2\p\b\tcl1fa6686758830\p\lib\tcl8.6`
   (find it with `dir /s /b %USERPROFILE%\.conan2\init.tcl`). The curated `dist\`
   resolves this; the raw build tree does not.
2. **`cmd` eats a trailing space into `set`.** `set VAR=400 && ...` assigns
   `"400 "`, which the option parser rejects as a non-integer. Write
   `set VAR=400&& ...` with no space.
3. **Build:** `set LADRUNO_CMS_BUILD=1 & Ladruno_scripts\build.bat` builds and
   runs the five standalone checks and exits. For the MP binary:
   `cmake --build build\build\Release --target OpenSeesMP -j 8` after
   `call Ladruno_scripts\setup_env.bat`.
4. **Stale processes lock the exe.** A killed `mpiexec` leaves
   `OpenSeesMP.exe` / `hydra_pmi_proxy` running and the next link fails with
   `LNK1104: cannot open file`. Kill them before rebuilding.
5. **Build in a worktree**, never the shared checkout (see the global notes).
6. Do not `sed -i` a `.bat` — it strips CRLF and `cmd` silently miscompiles
   multi-line blocks (`BUILD_GOTCHAS` #11).

---

## 6. What already exists (do not rebuild it)

**Decks** — `Ladruno_implementation/cms_profile/`:
- `cms_profile_chain.tcl` — 1-D chain. **Friendliest possible topology** (1-DOF
  interfaces): do not draw conclusions about S2/T1 from it. It has an analytic
  correctness anchor.
- `cms_profile_sheet.tcl` — 2-D quad sheet, real interface `b = 2*height`. This
  is the one to profile with.
- `cms_compare_arpack.tcl` — same mesh at np=1 (standard) and np=4 (CMS). Use
  this for any ratio claim.

**Instrumentation** (all under `-verbose`): hierarchy phases, T2 breakdown,
Lanczos breakdown + counts, `n/r2`, `r2/r*`, per-rank peak RSS.

**Options:** `-maxRestarts` is now exposed (§28; default 20). `-maxIter` is the
operator-application budget. The two interact — raising restarts without raising
`maxIter` just moves the failure.

**Known remaining inefficiency (measured, not fixed):** `rayleighRitz` is the
largest Lanczos term after §29's work (~0.029 s at 20x20). Its cost is the `j^2`
projection dot products, which could be made incremental — only new rows/columns
change as the basis grows — but that is surgery on the restart path and was
deliberately not attempted without evidence it matters at scale. **Experiment A
gives that evidence.**

---

## 7. Two lessons this lane already paid for — please do not re-buy them

1. **§26:** the obvious-looking dense buffer was refactored (O(dim²) → O(nnz),
   real memory win) and it moved wall clock **0%**. Measured A/B, ±3%.
2. **§29:** the next "obvious" target was `phi`/`psi`. Profiling first showed
   they are **0.317 MiB and 0.1% of T2** — the Lanczos was 80%. Optimising them
   would have been optimising noise. The actual fix (batching the MUMPS inverse
   action, caching `M*operatorBasis`) gave **1.65x**.

**Measure before you change anything.** Both wins and both non-wins in this lane
came from that rule, and every violation of it cost a session.

Also: a **1-D chain profile does not transfer**. §27 read `r2/r* = 1.06` from the
chain and concluded level 1 buys nothing; on the sheet it is **2.36**. Always
profile on the sheet or a real model.

---

## 8. State of the lane

- **Shipped this session:** Part 0 (MUMPS ordering, §21), P3b assembly oracle
  (§22), P3d complete (§23–§25), P4 §4 partial memory refactor (§26), P4 §1
  instrumentation (§27), `-maxRestarts` (§28), T2 Lanczos 1.65x (§29), the P4
  re-measurement (§30).
- **P3 open:** P3a and P3e are **blocked on apeGmsh** — see
  [[LadrunoCMS_apegmsh_emitter_guide]]. P3c extras open.
- **P4 open and now looking worse**, per §30.
- **CI:** the nightly self-hosted lane still does not run — the repo has **zero
  self-hosted runners**. The local box is the only harness.
- **Building 1A is still not in the repo**, so no measurement here is the real
  P4 gate.
