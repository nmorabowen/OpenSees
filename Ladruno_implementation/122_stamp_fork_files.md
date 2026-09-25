# WP-122 — Stamp the 31 unstamped fork files so the quirk lint scans them

Revision 1. Not adversarially reviewed: a comment-only change, proven by fingerprint (below); per
`feedback_adversarial_gate_when`, the full gate is for novel math, core/vanilla code or weak coverage.

Status: **complete; PR #857 marked ready (2026-09-25). The owner merges.**

Scoped 2026-09-25. Branch `wp/122-stamp-fork-files`, cut from `ladruno` @ `bc5c33453` (not stacked on
WP-120 #855). Discharges WP-120 R1 and its open question 2.

## Problem

`ci/check_quirk_patterns.py` scans only files carrying the `LADRUNO-HEADER-START` stamp. WP-120 compared
the stamp with git (files on `ladruno` that are neither on `upstream/master` nor at the merge-base
`e1237189a`). It found **31 fork-added files, 8,042 lines, without the stamp**: invisible to L1/L2/L4, the
WP-116 failure mode ×31. The same comparison found `Ladruno_scripts/stamp_headers.py` GLOBS out of step
with the tree: 10 stamped files missing from it, 5 entries matching no file.

## Shape

1. **GLOBS.** Add the 31 files (exact paths; `ASDPlasticMaterial3D/` is a vanilla directory, so the kit
   headers are listed individually or by `HoekBrown_*` / `StiffSoil_*`). Add the 10 stamped-but-missing
   files. Remove the 5 dead `ExplicitBathe{SMS,SMSConsistent,LNVD,LNVDSMS,LNVDSMSConsistent}.*` entries
   (#419 deleted those files).
   *Accept:* GLOBS resolves to exactly stamped ∪ the 31 (268 files), and no entry matches nothing.
   **Met** (268; 0 dead; 0 stamped files outside GLOBS).
2. **Stamp.** `python Ladruno_scripts/stamp_headers.py` stamped 36 files. 31 are new stamps. 5 had a
   hand-written 5-line stub between the markers, which is now the canonical block: `LadrunoAutoPenaltyReduce.{h,cpp}`,
   `LadrunoContactAbort.{h,cpp}`, `LadrunoParallelBuild.cpp`. `--check`: "All 268 authored files carry a
   current header."
3. **Prove comment-only.** `wp122_stamp/fingerprint.py` hashes every C/C++ file under `SRC`. Each hash covers
   the code with comments stripped by the lint's `clean()`, whitespace-normalised, plus every preprocessor
   line and every string literal outside comments.
   *Accept:* all hashes identical before and after. **Met:** 3,570 / 3,570 unchanged. The `SRC` diff is
   726 insertions, 0 deletions, and every added line is a `//` comment or blank. No rebuild was needed for
   evidence; Zone-A builds it on ready.
4. **Lint the widened scope, then triage.** `check_quirk_patterns.py`: 0 findings. That zero was checked,
   not assumed (step 5).
5. **Scope acceptance (mutation).** `wp122_stamp/scope_acceptance.py` plants the same L2 incident (a
   singleton `instance()` with no `clearAll` reset and no waiver) into `SRC/utility/LadrunoThreads.cpp` in
   `git archive` trees:

   | Tree | Expected | Result |
   |---|---|---|
   | this branch, unplanted | clean | PASS (0 findings) |
   | `origin/ladruno` + plant (file unstamped) | lint silent on the file | PASS (0 findings, the incident missed) |
   | this branch + plant (file stamped) | lint flags it | PASS (1 finding, rc = 1) |

6. **Ledgers.** A `LEDGER_implementations` row. `WORKFLOW_GOTCHAS.md` §5 gains why the stamp matters (it
   is the lint's scope) and the two traps behind the drift.
7. **CI gate (added at the owner's request, 2026-09-25).** A new `static-gates` step, "header stamp covers
   GLOBS (WP-122)", runs `python Ladruno_scripts/stamp_headers.py --check`. It sits just before the quirk
   lint, because a red stamp means the lint's scope is wrong. The required job name is unchanged.
   *Accept:* `wp122_stamp/stamp_gate_acceptance.py` runs `--check` in `git archive` copies of the branch
   (each tree runs its own copy of the script). GLOBS globs are case-sensitive on the Linux runner, so they
   were also re-matched with `fnmatchcase` against `git ls-files`: 268 files, 0 globs differing from the
   Windows match. **Met:**

   | Tree | Expected | Result |
   |---|---|---|
   | branch, unchanged | exit 0 | PASS ("All 268 authored files carry a current header.") |
   | stamp block removed from `LadrunoMassCache.h` | exit 1, file named | PASS |
   | stamp block stale (credit line edited) | exit 1, file named | PASS |

## Results

**Why the lint found nothing:** no rule has anything to read in the 31 files yet. They contain 0 calls to
`getRayleighDampingForces()` (L1), 0 `static X& instance()` singletons (L2), and 0 `Element` subclasses (so
none of L4/L5/L6). The value is forward-looking: the next edit that adds a Rayleigh term, a registry or an
element hook to one of these files is now checked.

**Process-wide mutable statics in the 31 files.** These are outside every rule; WP-115 rejected "lint every
mutable static". Each was checked by hand:

| Static | Kind | Reachable from a parallel region? |
|---|---|---|
| `CriticalTimeStep.cpp:360-361` `augSizeMismatchWarned`, `augSelfReportWarned` | one-shot warning latches | no — Δt estimation runs serially in the integrator |
| `LadrunoLoadControl.cpp:350` `noSPTermNotes` | note counter | no — integrator, serial |
| `LadrunoThreads.cpp:45` `ladrunoNumThreads` | thread-count setting | no — set from the command layer |
| `MohrCoulombTensionCutoff_YF.h:331` `warned_T_above_apex`, `:654` `apex_event_count` | latch and a non-atomic counter inside the yield function | **not today.** `ASDPlasticMaterial3D` does not override `ladrunoThreadSafeUpdate()` (`Material.h:72` default `false`), so the WP-107 allowlist runs any model containing it serially (`Domain.cpp:2638`). Allowlisting ASDPlastic later would make these two a data race: audit them first |

All six survive `wipe`, so a latch warns once per process rather than once per model. That matches the
one-shot-latch convention WP-115 accepted.

## Rejected approaches

- **Wildcard GLOBS for the ASDPlastic kit** (`ASDPlasticMaterial3D/**/*.h`). That directory is vanilla
  (Abell's upstream framework); a wildcard would stamp upstream files. Exact paths, plus two prefix globs
  (`HoekBrown_*`, `StiffSoil_*`) that match only fork-added files today.
- **Leave the 5 stub blocks as they were.** They were outside GLOBS, so `--check` never saw them. Bringing
  them under GLOBS makes the script own them; the canonical block replaces the stub.
- **Skip the ASDPlastic kit headers as "José's code".** They are fork-added (not on `upstream/master`),
  and José Abell is one of the four credited authors. The upstream campaign
  (`upstream_pr_campaign.md` §3, package 1.5) re-stamps headers in its own scrub step anyway.

## Open questions

1. **The second trap is still open.** CI now runs `--check` (step 7), which catches a GLOBS file that
   loses its stamp. It cannot see a fork file that was never added to GLOBS; that is how all 31 escaped.
   Further options: (b) fail if a file whose path contains `ladruno`/`Ladruno` lacks the stamp (would have
   caught 14 of the 31, nearly free); (c) WP-120's `inventory.py` against `upstream/master`. That is exact,
   but CI would have to fetch the upstream remote. Owner's decision, together with WP-120 (ii).
2. **Before allowlisting ASDPlasticMaterial3D for the WP-107 threaded loop**, make
   `MohrCoulombTensionCutoff_YF.h`'s two statics atomic or per-instance.
