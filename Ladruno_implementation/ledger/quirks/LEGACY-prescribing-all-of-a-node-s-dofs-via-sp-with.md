---
wp: LEGACY
title: "Prescribing ALL of a node's DOFs via sp with the Transformation handler ⇒ 0 free equations ⇒ process terminates"
legacy_seq: 56
---
### Prescribing ALL of a node's DOFs via `sp` with the Transformation handler ⇒ 0 free equations ⇒ process terminates
- **FIXED 2026-08-13** -- see the `FullGenLinSOE::getX - vectX == 0` entry below for the
  root cause (six SOEs, null size-0 `Vector` wrappers) and the fix. This entry is kept
  because its *symptom* description (pytest aborting with no summary, looking like a
  hang) is the one you are most likely to search for.
- **Bites:** a static test that imposes both DOFs of the only free node via `ops.sp`
  (with the other node fully `fix`ed) leaves the system with **zero unknowns**.
  Under `constraints('Transformation')` the solve does not return an error code —
  it **terminates the process** mid-`analyze()` (no Python traceback, exit 0),
  which under pytest aborts the whole run with no summary (looks like a hang).
- **Why:** the Transformation handler condenses out the constrained DOFs; with none
  left the assembled system is degenerate and the path hits a hard exit rather than
  a graceful failure (same family as the MPCO `exit(-1)` kernel-kill pattern).
- **Workaround:** use `constraints('Penalty', 1e14, 1e14)` for fully-prescribed
  configurations (DOFs are retained and penalised, so ≥1 equation remains), and
  read element force via `eleForce` rather than `nodeReaction`. Learned 2026-06-02.
