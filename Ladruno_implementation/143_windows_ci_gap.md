# WP-143 — the Windows-only test gap in CI

Status: DONE (draft PR #884). Found by WP-136 (#870).

## The gap (measured 2026-09-27)

- 9 `zone_a` test files skip on non-Windows (`sys.platform != "win32"`), and
  Zone-A runs on `ubuntu-latest`, so PR CI never runs them:
  `test_ladruno_sanisand_flip_determinism.py`, `test_ladruno_sanisand_replay_counters.py`,
  `test_wp132_deterministic_pardiso.py`, `test_pardiso_stats.py`, `test_pardiso_solver.py`,
  `test_pardiso_asym_rearm.py`, `test_feastEigen.py`, `test_adr97_p4_inertness.py`,
  `test_wire_venv_pth_override.py`.
- The self-hosted Windows nightly jobs (`zone-b-nightly`, `cross-tier-nightly`)
  were `cancelled` on all 100 scheduled runs from 2026-06-20 to 2026-09-27. The repo has
  zero registered runners. These 9 files have had NO CI for at least three months.
- `cross-tier-nightly` has no build step: even with a runner online it tests
  whatever `opensees.pyd` is installed on the box, not the checked-out commit.

## Scope

1. Make portable what is not MKL-specific. Each of the 9 files is triaged: legs
   that need MKL/Pardiso/FEAST/Windows stay gated; the rest run on any solver
   so Zone-A covers them. Done for flip-determinism: its hold-count and warning legs
   now run everywhere; the thread leg stays Windows-only.
2. Quirk lint **L8** (L7 is WP-140, #880): a `zone_a` test with a win32 skip must
   carry `# ci-coverage: <where it runs>`, or the lint fails.

Out of scope (owner): register the `ladruno-perf` runner; add a build step to the
nightly jobs; an optional PR-triggered Windows job.

## Result

Inventory after this WP (`python ci/check_quirk_patterns.py --list-waivers | grep ci-coverage`):

| kind | files |
|---|---|
| local-only (MKL-only, no CI until a Windows job builds and runs them) | `test_pardiso_solver`, `test_pardiso_stats`, `test_pardiso_asym_rearm`, `test_wp132_deterministic_pardiso`, `test_feastEigen` |
| partial (run on Ubuntu at a floor; the strict leg is Windows-only) | `test_adr97_p4_inertness`, `test_ladruno_sanisand_replay_counters`, `test_ladruno_sanisand_flip_determinism` |
| portable | `test_adr74_numberer_1` (taskkill vs kill cleanup) |

The five local-only files cannot be ported: Pardiso and FEAST are MKL. Only the
owner items close them.

## Verification

- L8 flagged exactly the 9 files before any annotation, and reports 0 findings after.
  Self-tests: 67/67, including 16 L8 cases covering each gate form, the ternary
  exemption, text-only mentions, unknown kind, short reason, stale annotation and an
  unparseable file.
- Flip file on Windows (build `3cc41f12a`), Pardiso path: 3/3 pass.
- Flip file on Windows, portable path (`LADRUNO_FLIP_SYSTEM=FullGeneral`): 2 pass, and
  the thread leg is skipped with its reason.
- Break-test on the portable path (default mapped to `vanilla`): the hold leg fails
  (spread 0.576) and so does the warning leg. So the legs that now run on Ubuntu are a
  real F14 gate.
- Ubuntu, via `workflow_dispatch` run 36371695499, Zone-A green: 2749 passed,
  149 skipped, 14:43.
  - Flip hold + warning legs PASSED; the thread leg SKIPPED.
  - `adr97_p4` and `replay_counters` ran and passed.
  - FEAST SKIPPED, as annotated.
- Static gates on the draft (Ubuntu): lint 0 findings `[L1..L6,L8]`, self-tests 67 passed.
