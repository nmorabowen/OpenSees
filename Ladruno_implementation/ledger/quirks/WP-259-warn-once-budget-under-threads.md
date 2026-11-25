---
wp: WP-259
title: "A process-wide warn budget under threads: claim with one fetch_add, print by the claimed index"
date: 2026-11-25
---
### A warn-once / "first N lines" budget under threads overshoots or races — claim it with one `fetch_add` (WP-259)

- **Symptom.** A material that prints "the first N failures, then a suppression line" keeps the
  count in a process-wide static (`LadrunoConcrete3D`'s `s_returnMapWarnPrinted`, budget 20).
  Written as `if (count < N) print; if (++count == N) suppress;` it is a plain read-modify-write:
  under the WP-107 / WP-192 loops two workers lose increments (the suppression line can be
  skipped or printed twice). Making the counter `std::atomic` but keeping the separate check and
  increment (WP-258's DruckerPrager NormCep budget) is race-free but can print up to N + T − 1
  lines with T threads.
- **Fix.** One claim: `const long k = budget.fetch_add(1);` then `if (k < N) print;` and
  `if (k + 1 == N) print the suppression line;`. Serially this is the old sequence line for line
  (k is the old count); threaded, exactly N lines and one suppression line are printed, in
  whatever order the claims land. Each print still goes through `LadrunoOutputGuard`.
- **Gate.** `tests/test_wp259_concrete_threads.py` `3dcr` pushes a stdBrick prism past the crush
  so return maps fail in volume on every thread count, and asserts the 1-thread run printed
  exactly 20 lines plus the suppression line and every threaded run the same counts.
- **Reset.** A `wipe`-time reset is `store(0)` on the main thread, between loops.
