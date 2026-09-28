# WP-143 — the Windows-only test gap in CI

Status: IN PROGRESS (draft PR). Found by WP-136 (#870).

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
   so Zone-A covers them. First target: flip-determinism's hold-count and warning
   legs (the thread leg stays Windows-only). That waits for #870 to merge, to avoid
   editing the same file on two branches.
2. Quirk lint **L8** (L7 is WP-140, #880): a `zone_a` test with a win32 skip must
   carry `# ci-coverage: <where it runs>`, or the lint fails.

Out of scope (owner): register the `ladruno-perf` runner; add a build step to the
nightly jobs; an optional PR-triggered Windows job.
