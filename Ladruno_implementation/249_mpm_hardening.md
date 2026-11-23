---
title: WP-249 — MPM hardening before phase 7 (#1095, #1110, #1103, #1106)
project: Ladruno
tags:
  - mpm
  - wp-249
---

# WP-249 — MPM hardening before phase 7

Program #1075. A fix/hardening WP queued by the owner on 2026-11-22, to land before
phase 7 (232m, multi-body contact). It is not a §9 phase of ADR 232.

Scope, in priority order:

1. **#1095** — the mode-A crash on a step-cut retry (`Matrix::operator*(Vector):
   incompatible sizes`, then 0xC0000005). Root cause, fix, and a deterministic test that
   forces cuts across a connectivity change.
2. **#1110** — check whether 232l's `mpm traction` is the single-phase surcharge carrier
   #1110 asks for, in both modes; demonstrate it on 232k's S2 column, or record the gap.
3. **#1103, #1106** — the small, local items of the 232i and 232j code-review follow-ups.

Out of scope: #1107, #1124–#1127, #1141, contact.

## Status

Ready for review (2026-11-23). #1095 stays open with a guard test and no fix (owner's call,
2026-11-23); #1110 narrowed; the small #1103/#1106 items done.

## #1095 — not reproduced; guard test only

- **Code reading** (cut/rollback path, mode A). A step that fails before commit is cut and
  retried with nothing to undo: `prepareStep` rebuilds supports, active set, node diff, P2G,
  repoint and SPs from the committed particles, and `StaticAnalysis::analyze` reverts the
  Domain. DruckerPrager3D's `revertToLastCommit` is a no-op and its Matrix×Vector products have
  fixed sizes. LadrunoMPMPoint does not opt into threaded element loops, so at 8 threads only
  MKL/PARDISO is threaded: threads change the numbers, not the code path. No size mismatch was
  found on the path.
- **Reproduction attempts** on the #1095 configuration (N = 1e5, DP, quasi-static, PARDISO,
  renormalize, 8 threads): `c45af1407` 2/2 clean; `9db5803f6` (the build that crashed) 2/2
  clean; `2236d1ca1` built with debug symbols (`/Z7`, `/DEBUG:FULL`) under a Win32 catcher
  with breakpoints on both `incompatible sizes` prints in `Matrix.cpp` and a minidump on an
  access violation, along r2's cut path: no breakpoint, access violation or dump in either
  repeat; one completed, the other hit the 2 h timeout during step 8's last cut. Also clean at
  N = 1e4, 1 thread. On `9db5803f6` and the symbols build, step 7 exhausts `mpm cut 2 4` and
  raises a normal `OpenSeesError`.
- **Reading.** The 1-in-3 crash looks like a corruption that depends on the run, not a
  deterministic defect on the cut path. #1095 stays open (owner, 2026-11-23).
- **Guard.** `tests/test_ladruno_mpm_wp249_hardening.py::test_cut_retry_across_connectivity_changes_3d`
  forces failed attempts (`mpm debug failAttempt … afterTrial`) on two steps whose own
  preparation re-points the particles (a connectivity change at the failed attempt, asserted on
  the reference run's first sub-step; the per-call retry count pins which step is armed) on a 3D
  DP block, and checks the retried run equals the direct sub-steps (rel 1e-12), with none of
  the #1095 messages. It passes on the base too: **no mutation evidence**; it guards the path,
  it fixes nothing.

## #1110 — checked, narrowed

`mpm traction` is a grid-face load (nodal forces λ t ∫hat_a on a grid boundary row), GIMP/linear
only. Single-phase particles reach it in mode A (general path, ndf = ndm) and, by the code, in
mode F (route S). On S2 with the grid top at the soil top and gravity as the unscaled particle
body force: σ_yy max / ρgl0 = 0.108 (HS) and 0.033 (SANISAND) against the gate's 0.03, from
0.27 with the gated top-row body force; the residual sits in the top cell. No gated number
changes. #1110 stays open for B-spline, a non-grid surface and the top-cell residual. Probe:
`Ladruno_scripts/mpm_249/`. Quirk: `ledger/quirks/WP-249-mpm-traction-is-a-grid-face-load.md`.

## #1103 and #1106

| Item | Disposition |
|---|---|
| #1103-2 restart pack v1 | Done: the reader accepts version 1 (no body force); `mpm debug packVersion 1` writes it and falls back to the current layout for a body force or a B-spline basis. Two tests; both fail on the base, which refuses `packVersion 1` (neither reads back which version was written; the B-spline fallback has no test; no pack from the v20261116 binary is used -- the layout was checked against its source). The fallback covers only the body force and B-spline: an anti-locking flag or a v4 gridvel field is dropped by a v1 pack, as by 2 and 3. |
| #1103-1, -4 | Left: they touch gated B-spline paths. |
| #1103-3 | Left: unreachable today. |
| #1106-2 `failDump` switch | Done: cleared at mode-A teardown, re-asserted from the live host; `mpm info failDump`. Test; mutation (reset removed) fails it. |
| #1106-4 | Done: `addFbarTangent` scratch is `thread_local`. |
| #1106-5 | Done: `detN`'s 3×3 branch is `ladruno_hypo::det3` (bit-identical). |
| #1106-1 | Left: gated anti-locking numerics. |
| #1106-3 | Left: belongs with #1097. |

## Verification

Five targets built at `2236d1ca1` (stamps match on exe and pyd). WP-249 tests 4/4; MPM suite
(`tests/test_ladruno_mpm_*.py -m "not zone_b"`) 514 passed, 15 xfailed, 61 skipped.

## Review of record

- **Fable critic (other family): AMEND, applied.** (1) The cut-retry test armed steps whose
  own preparation did not re-point: the "across a connectivity change" assertion passed on later
  sub-steps. Re-aimed at steps 2 and 11, whose preparation re-points (asserted on the
  reference's first sub-step), with a per-call retry count `[3, 3]` that pins the armed step (an
  off-by-one in the re-aim, step 12, made the runs take different step sequences and differ by
  8.6e-5; the count now catches it). The re-aimed guard passes on `c45af1407` too. (2) Wording:
  "regression" → guard, "bit for bit" → rel 1e-12. (3) The v1 fallback comment narrowed to what
  it covers. (4) Pack-test limits recorded above. (5) `mpm info failDump` reports the process
  switch, not the host's path (kept; documented). (6) The traction quirk's mode-F half corrected.
- **`/code-review high`:** 4 findings. The timed-out symbols run counted as clean: fixed. The
  `failDump` report, the v1 readback and the v20261116 fixture: kept as notes (C++ + rebuild, or
  the old binary).
