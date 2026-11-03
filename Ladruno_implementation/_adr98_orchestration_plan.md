---
title: "ADR 98 orchestration plan — agents, models, effort, work packages, and the orchestrator's protocol"
project: Ladruno
status: draft — awaiting the owner's approval of §1 (go/no-go item 1, §6)
orchestrator: "Fable (house precedent: ADR-90, ADR-94) — the owner's session; coordinates, synthesizes, does not merge"
companion: 98_ladruno_mpm_scoping_adr
updated: 2026-09-26
---

# ADR 98 orchestration plan

ADR-87 requires that "agent assignment and effort level per WP are decided at
kickoff and recorded". This plan does that for ADR-98, in the format of
[[_adr90_orchestration_plan]] and [[94_asdplastic_review_plan]]. Nothing here is
approved until the owner signs §6.

## 0. Two things that are deliberately different

- **Only M0 is fully assigned.** ADR-98 has no architecture winner yet: A, F and C
  are finalists and M0b decides. M1–M7 change shape with the winner (F adds an
  engine-owned explicit loop; A adds a particle element; C adds materialization),
  so their rows in §2 are **provisional** and get re-decided at the M0c kickoff.
- **"Effort" is not yet operational for working agents.** The `Agent` tool takes a
  `model` only; reasoning effort comes from an agent definition's frontmatter, and
  this repo has none (`.gitignore:151-153` ignores `.claude/*` except `skills/`; no
  user-level `~/.claude/agents/` exists). Until go/no-go item 2 is settled, effort in
  this plan means two things, as in ADR-90: **the `/code-review` level per PR** (fully
  operational), and **recorded intent for working agents** (what the definition
  would say). **The orchestrator is the exception:** it is the owner's main session,
  not a subagent, so its model and effort are ordinary session settings and are
  fully operational.

## 1. Model and effort selection

| Role | Model · effort | Rationale |
|---|---|---|
| **Orchestrator** — sequencing, decisions, synthesis, checkpoints, PR state; launches builds and long runs. **The owner's main session**, so both settings are real | **Fable · high**, raised to **max** at the three decision points — adjudicating the 98a re-review, the 98e verdict, the M0c decisions — and back to high after | house precedent (ADR-90, ADR-94 R6); holds the only context spanning ADR-86, ADR-88, ADR-98 and the measurements. High by default because most orchestration is sequencing, launching and reading results, and polling a build at max burns tokens for nothing; max where one wrong call commits months of work. Effort is not the whole defence: both ADR-86's and ADR-98's blockers were missed by an author session running at max, and both were *cross-document* misses — the procedural fix is §3 step 8 |
| **Adversarial critics** — the owed ADR-98 re-review, and each family PR | **Fable · high** (lead) + **Opus · high** ×2 lenses, background, independent prompts | Fable found the blocker in both ADR-86 and ADR-98; critics must not share the author's context — ADR-98's author is an Opus session, so the Fable lead is a different model, and the Opus lenses run in fresh contexts |
| **Numpy oracles** (M0a; V3/V4 at M2; V8 wiring at M4) | **Opus · high** | correctness-critical numerics; precedent: ADR-90's constitutive oracle, ADR-71 P0 |
| **C++ numerics** — §6-0 re-seed contract (M0-pre), F's explicit engine (M0b-iii), every M1–M7 core | **Opus · high**, **one agent per WP**, continued via `SendMessage` across review rounds | vanilla integrators, byte-identity, P2G/G2P, tangents; precedent: ADR-90 WP-C |
| **Measurement harnesses** — the synthetic churn deck (M0b-i), re-seed timing (M0b-ii), sweep matrices | **Sonnet · medium** | a harness written once and run many times (ADR-94 rationale); the orchestrator launches the runs |
| **Interpretation critic** on the M0b numbers and on validation decks (V6, V7) | **Opus · high** | physics interpretation decides the winner; precedent: ADR-90 WP-D's physics-lens critic |
| **Wiring** — ledgers, verification manifest, guide, banner, stale docs, apeGmsh seeding v0 | **Sonnet · medium** | well-specified, pattern-following; precedent: ADR-90 WP-B / WP-E |
| **Read-only sweeps** during any WP | **Explore (Sonnet)**, "very thorough" | keeps the orchestrator's context lean |
| **Out-of-family verdict** (ADR-87 D4) at family level | **SuperGrok**, run manually by the owner from `reviews/handoff_adr98.md` | the citable review; same-family panels are bug-finding tools |
| Haiku | **not used** | nothing here is mechanical enough to be safe (ADR-90, ADR-94) |

**Rationale in one line:** Opus where mechanics judgment decides the answer; Sonnet
where a harness or wiring is written once and run many times; Fable to orchestrate
and to attack.

## 2. Work packages

Branches follow ADR-90's lettered scheme under ADR-87 D9: one branch per WP, a draft
PR on day one, one deliberate merge by the owner.

### Fully assigned (M0)

| WP | Branch | Content | Working agent | `/code-review` | Gate |
|---|---|---|---|---|---|
| **98a** | `wp/98a-mpm-scoping-adr` | ADR-98 + this plan + the three `LEDGER_quirks` rows (W1, W4, W6); the **owed re-review** of the revised ADR | orchestrator; critics per §1 | `low` (docs) | re-review folded; owner sign-off on §6 |
| **98b** | `wp/98b-integrator-reseed` | **M0-pre** — `TransientIntegrator::reseedFromCommitted()` in Newmark / HHT / GeneralizedAlpha / CentralDifference / CentralDifferenceLadruno / DynamicRelaxation; the `Domain` flag; the drivers; `LEDGER_vanilla_files` rows | **Opus · high** (one agent) + **Sonnet · medium** for ledger wiring | `high` | byte-identical when the flag is never raised; the W1 falsifier deck turns green |
| **98c** | `wp/98c-mpm-oracle` | **M0a** — numpy 1D/2D MPM: explicit MUSL (the F engine's reference) + implicit Newton; linear / GIMP / quadratic B-spline; V1, V2, V5 | **Opus · high**, under the ADR-71 independence protocol (reads nothing under `SRC/`) | `medium` | closed forms reproduced; convergence orders measured |
| **98d** | `wp/98d-mpm-measurements` | **M0b** — (i) churn deck on the existing fork: add/remove N = 10³…10⁵ elements between `analyze` calls, `system Diagonal` and PARDISO; (ii) re-seed cost from 98b; (iii) F's explicit engine as the 98c oracle's C++ twin, elastic and SANISAND particles | (i)(ii) **Sonnet · medium**; (iii) **Opus · high** — a *different* agent from 98c's author, meeting it only at a frozen interface contract | `medium` (i, ii) · `high` (iii) | three numbers written into ADR-98 |
| **98e** | `wp/98e-mpm-verdict` | **M0b synthesis + M0c** — declare the winner among A / F / C; decision memo for contact, parallel, D2 re-confirmation | **orchestrator** (reads the result files, not transcripts) + **Opus · high** interpretation critic | `low` (docs) | owner sign-off; M1+ rows below re-decided |

### Provisional (re-decided at the 98e kickoff)

| Phase | Working agent | `/code-review` | Extra review |
|---|---|---|---|
| M1 — winner, 2D + 3D, explicit + implicit, particle recorder, both interpreters, manifest, mutation gate, guide | **Opus · high** core; **Sonnet · medium** wiring; apeGmsh seeding v0 as a separate apeGmsh WP (**Sonnet · medium**) | `high` | **D4 out-of-family verdict** (first family PR) |
| M2 — GIMP + quadratic B-spline; V3, V4 | **Opus · high** core; **Opus · high** oracle (independent agent) | `high` | — |
| M3 — locking, undrained; V6, V7 | **Opus · high** | `high` | **Opus** physics-lens critic |
| M4 — two-phase particles via `LadrunoUPKernel.h`; V8 | **Opus · high** | `high` | — |
| M5 — FE → MPM handoff (depends on ADR-88 S1–S2) | fork side **Opus · high**; apeGmsh transfer side **Sonnet · medium**, with the identity-round-trip gate interpreted by **Opus** | `high` | — |
| M6 — contact (multi-velocity-field MPM–MPM; FE–MPM penalty) | **Opus · high** | `high` | **D4 out-of-family verdict** (contact is where this fork's blockers have hidden) |
| M7 — parallel, if IN | **Opus · high** | `high` | **Opus** critic on determinism and W6 |
| Skill `material-point-method` (ADR-86 D3) | content **Opus · medium** from ADR-98 §2–§8; repo + Seafile scaffolding and `resync-skills.bat` **Sonnet · medium** | — | — |

## 3. Orchestrator protocol

1. Open every WP as a draft PR on day one (ADR-87 D9); never merge; never
   `gh pr merge --auto`.
2. **Launch builds and long measurements yourself**, via WMI, watched by a bounded
   foreground loop on the artifact's own mtime — never by the working agent, never a
   detached background watcher (ADR-97 verdict §7.1; ADR-94 closeout).
3. Read result files, not transcripts (ADR-94 R6).
4. Review loop: `/code-review <level>` → findings → `SendMessage` to the **same**
   agent, context intact → re-review. Critics are separate background agents with
   independent prompts (ADR-90).
5. If a session stalls, resume it with the **finisher pattern** — a new session
   inspects the worktree's `git log`, `git status` and last build log — as the
   default path, not the fallback (ADR-97 verdict §7.2).
6. Keep the oracle author and the C++ author different agents that meet only at a
   frozen interface contract (ADR-71 P0 independence protocol).
7. After each WP, update ADR-98's implementation log (§14) with who ran what, on
   which build, with which result.
8. **Before accepting or writing any "works with X" / "X is missing" claim**, open the
   code that X actually runs through and re-read the decision and implementation-log
   sections of the companion ADRs. Both ADR-86's and ADR-98's blockers were facts
   already written down elsewhere in the fork; the author had read one of them the
   same session.

## 4. Standing rules handed to every agent (verbatim in the prompt)

1. Work only in your own worktree; never write in the shared checkout.
2. **Licence rule:** never open the source of GPL / LGPL / unlicensed codes (Anura3D,
   Karamelo, AMPLE, GeoTaichi, Matter, NairnMPM) — papers and documentation only.
   Kratos (BSD-4) code may enter the fork with acknowledgement but **never apeGmsh**.
3. **W1:** transient integrators ignore nodal state written between steps unless
   `domainChanged()` runs or the §6-0 re-seed is raised. Never rely on a node write.
4. **W4:** only `recorder ladruno` may watch entities that can be removed.
5. Ledgers (`LEDGER_vanilla_files`, `LEDGER_implementations`, `LEDGER_quirks`) are
   updated in the same PR as the change.
6. Stage edits as scripts, so work survives a session restart (ADR-94 closeout).
7. Run a test file once before pushing any edit to it, however trivial.
8. MSVC green is not GCC green; cross-platform float pins use ≥ 1e-6 relative, never
   1e-9 (ADR-94 closeout).
9. Every number you report names the build hash (`ladrunoBuild`), the deck, the
   output path and the date.
10. Commit trailer per the harness; no `--no-verify`; never merge.

## 5. Sequencing

98a first (the re-review must clear before anything is built). Then 98b and 98c in
parallel — they share nothing. 98d(i) needs neither and can start with 98b;
98d(ii) needs 98b; 98d(iii) needs 98c's frozen interface. 98e closes M0. M1+ waits
for the 98e decisions. No calendar is given: under ADR-87 an unmeasured duration is
a claim.

## 6. Go/no-go items for the owner before kickoff

1. **Approve the §1 model/effort table** — or reassign.
2. **Make "effort" operational, or not.** Either add agent definitions (e.g. a
   numerics role, a harness role, a critic role) under `.claude/agents/` with model
   and effort in their frontmatter — which needs a `!.claude/agents/` exception in
   `.gitignore`, as WP-115 did for `skills/` — or keep effort as recorded intent plus
   `/code-review` levels, as ADR-90 did.
3. Approve the re-review panel for 98a (Fable lead + two Opus lenses).
4. Approve 98b's vanilla footprint (six integrators and the drivers).
5. Confirm the D4 out-of-family verdicts at M1 and M6.
6. **Choose the orchestrator session.** Either switch the session that wrote ADR-98
   (an Opus 5.5 session, currently at max effort) to Fable — the transcript and its
   context are kept — or start a fresh Fable session that reads ADR-98, this plan and
   §15's review log (the finisher pattern), which is independent of the author's
   context. Either way, set effort to high and raise it at the decision points.
