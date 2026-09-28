# TIMs 2D-model requests F18–F23 — implementation plan and orchestration

Plan for `_tims_2d_model_requests_2026-09-25.md` (the intake, carried on WP-127's
branch so it lands on `ladruno` with the first WP). Written 2026-09-26 after a review
pass that verified every §2 source citation against `fb1afe58b` (no commit touched
`SRC/material/nD/UWmaterials/` since the act's `79e062367`, so all line numbers hold).

## 0. Review findings the plan is built on

- **A — attachment sign convention.** The CSV README says "compression negative as
  OpenSees stores it". On all 80 rows `p_kPa == +tr(σ)/3` and every normal stress is
  ≥ 0: the columns are the **internal, compression-positive** `mSigma`. The replay
  facility (F21) takes an explicit convention switch and the report tells the act.
- **B — η = √(3/2)‖α‖ at the ring points.** With `m = 0.005` the stress rides on the
  back-stress; at element 1950 gp 2/3 (the p' 0.352 kPa, η 12.87 point) **α itself is
  ~6× outside the bounding surface** (M^b ≈ 2.10) while gp 1 of the same element sits
  at η 1.54. Hypothesis: explicit-substep overshoot of the stiff α law
  dα ∝ h(α^b − α) when h → ∞ as (α − α_in):n → 0 (just after a flip). If it holds,
  F18(a) (error floor) treats a symptom; the cure is α-stability. **F18(e) is
  therefore moved ahead of F18(a)/(b).**
- **C — silent accept at `dT_min`.** `ManzariDafalias.cpp:1884-1896`: a substep that
  FAILS the error test at `dT == dT_min` is accepted with elastic tangent, a radial
  clamp to `m_Mc` (not M^b, compression side only) and a re-derived α. Uncounted.
  F20(a) counts it.
- **D — F18(b) is not round-off neutral for `TanType 2`.**
  `aCep_Consistent = ½(aCep1+aCep2)·(aD·aCep_Consistent + T·mIImix)` (`:1916-1918`)
  consumes the per-stage 6×6s. Rate-form stages must keep stage tangents when
  `TanType == 2` (or the claim is limited to stresses + TanType 0/1).
- **Minor** — b8 row 1859/2 has tr(α) = 2.3e-3 (others ~1e-10). Replay projects α, z
  to deviatoric and warns.

## 1. Work packages

All: branch `wp/<n>-<slug>` from fresh `origin/ladruno`, **draft PR day one with
`--base ladruno`**, commit continuously, `// Ladruno` on every touched upstream line,
ledgers + guide + banner in the same PR, pytest before/after pasted in the PR, PR ends
with verified / not verified / ledger rows added. Agents never merge and never flip
to ready; the owner does. A WP that needs an unmerged sibling's code is cut from that
sibling's tip but STILL opens with `--base ladruno` (global CLAUDE.md stacked-PR rule)
and says so in its body.

| WP | Items | Depends on | Build | Kind |
|---|---|---|---|---|
| **127** `sanisand-replay-counters` | F20(a) cumulative substeps / last-call / `-maxSubsteps` cap hits / **accept-at-dT_min** count (finding C), per instance, since `revertToStart`; F21 state replay (load σ, α, α_in, z, e with explicit sign convention; drive dε; report substeps, per-substep `T, dT, err` trace, returned stress) + Python helper; intake doc + this plan | — | yes | C++ (vanilla MD + LadrunoSANISAND) |
| **128** `sanisand-ring-trace` | F18(e): replay the 80 ring states; test finding B (α overshoot mechanism); say plainly whether any integrator can take the b8 worst point and how a committed η/M^b ≈ 6 arises; baseline substeps/errors vs a reference (`IntScheme 45` tight / ME 1e-8) — the numbers F18(a) needs | 127 binary | no (uses 127's `dist`) | investigation + report |
| **129** `sanisand-errfloor-rateform` | F18(a) `-errFloor` (byte-identical default) + expose the hard-coded `TolE 1e-4` as a quirks row; F18(b) rate-form stages honouring finding D; F20(b) `OPS_PROFILE_SCOPE`s; F20(c) `"tangentEP"` vs numerical tangent; plus the α-stability fix if 128 confirms B | 127, 128 verdict | yes | C++ |
| **130** `sanisand-cppm-under-newton` | F18(c) CPPM refuses at once with the element-forwarded refusal code, local line search / better start, **de-static `NewtonIter`**; F18(d) per-point ME→CPPM fallback behind a flag; rerun F12 bearing deck | 127 | yes | C++ |
| **131** `sanisand-threaded` | F19. **Step 1 (now, no build): the inventory** of shared mutable state reachable from `setTrialStrain`/`commitState` in MD + LadrunoSANISAND with every write site. Step 2 (after 130): per-instance / thread-local / locked; 1/2/4/8-thread identity + speed-up | 130 for code | step 2 | doc → C++ |
| **132** `deterministic-mode` | F22: MKL CNR (`MKL_CBWR`) + Pardiso `iparm[33]`, a documented list of remaining order-dependence; runtime banner/print when on (the splash list is compile-time) | — | yes | C++ + guide |
| **133** `pdmy03-cs-params` | F23(a) optional `ei, cs1, cs2, cs3` (byte-identical default) + quirks row; F23(b) guide note (dilation brake = void-ratio switch, `PressureDependMultiYield.cpp:2189-2214`) | — | yes | C++ (vanilla, additive) |
| **reply** | `_tims_2d_model_report.md` answering F18–F23 + findings A–D, assembled as WPs land | all | — | doc |

## 2. Orchestration

- **Wave 1 (parallel):** 127, 131-step-1, 132, 133.
- **Wave 2:** 128 as soon as 127's build is green (runs on 127's `dist`).
- **Wave 3:** 129 (after 128's verdict) and 130 (after 127) — 130 can overlap 129.
- **Wave 4:** 131 step 2 after 130; the reply doc.

**Build mutex.** One full `build.bat` at a time on this box (16 threads, 28 GB). An
agent acquires the lock by atomically creating the directory
`%LOCALAPPDATA%\Temp\ladruno_build.lock` (`mkdir` fails if held), writes its WP id
into `owner.txt` inside it, builds, and removes the directory when the build ends
(success or failure). Waiting agents poll every few minutes. Tests may run
concurrently with another agent's build.

**Test runner.** CPython 3.12 at
`C:\Users\nmora\AppData\Local\Python\pythoncore-3.12-64\python.exe`, launched with
`-S` + manual paths to the WP worktree's own `dist\bin`, and assert
`opensees.__file__` is that worktree's pyd (a boot `.pth` otherwise preloads a stale
pyd from another worktree).

**Adversarial gate.** Required for 129 (new math, vanilla), 130 (vanilla, numerics),
131 step 2 (shared state). Not for 127/132/133 beyond the test battery + `/code-review`.
