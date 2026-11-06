---
wp: WP-130
title: "CPPM at low p and at strain reversals branches on ROUND-OFF -- a one-ulp input change flips the path; never pin such decks' numbers cross-platform (WP-130)"
legacy_seq: 516
---
### CPPM at low p and at strain reversals branches on ROUND-OFF -- a one-ulp input change flips the path; never pin such decks' numbers cross-platform (WP-130)
- **Bites:** the WP-130 byte-id decks passed on Windows and failed on the first Linux Zone-A run
  (every earlier branch run was a draft, where Zone-A is skipped). `ls_ps_s2_cyc` (p ~ 1-3 kPa,
  `-Pmin 0.0101`) takes the full 2^9 halving ladder plus explicit at almost every step. WHICH
  sub-increments converge is decided by round-off: the Windows/Linux census matches through row 23
  (fixed tangent: row 14) and diverges at the next reversal. On ONE Linux binary a 1e-15 relative
  change of the confinement strain moves row 15's tangent from 41874.9 to 61786.6, the Windows value;
  1e-11 moves state entries by up to 1.4x. The opt-out free decks' global Newton wanders too (iteration
  counts 9 vs 8, 31 vs 18). And a homogeneous quad's 4 GPs differ only by round-off, so which of them
  refuses is platform-dependent (MSVC all 4; GCC GP 2-4, 3-4 or only 3). No code defect: the
  pre-WP-130 source built on Linux is bit-identical to the WP-130 opt-out there; `-Wall -Wextra
  -Wuninitialized -Wmaybe-uninitialized -Wsequence-point` is clean.
- **Rule:** in these regimes pin INVARIANTS (opt-out == pre-change on the same platform) or STRUCTURE
  (rows, rc, "the ladder ran", GP-summed or any-GP census), not trajectories. Compare stable decks
  per entry (1e-6 relative + 1e-10 of the row's quantity group), not against a deck-wide floor: a
  `1e-6 * max|entry of the deck|` floor is dominated by the tangent (~0.3 absolute on ~1 kPa
  stresses) and makes the stress comparison vacuous.
- **Status (WP-130, #868):** `test_ladruno_sanisand_cppm_newton.py` compares md3d_s2, ls3d_s2,
  ls3d_s2_big, ls_ps_s2 (and the fixed-tangent free decks) per entry on every platform (smallest
  measured margin 78x); the round-off decks keep their numeric pins on win32 only, with a structural
  check elsewhere; the halvings test reads all 4 GPs. Related: the c < 7/9 Lode non-convexity makes
  extension paths round-off-selected too (WP-151 memo §6.3).
