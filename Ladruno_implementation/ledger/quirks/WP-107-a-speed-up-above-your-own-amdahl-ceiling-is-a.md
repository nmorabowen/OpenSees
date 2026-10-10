---
wp: WP-107
title: "A speed-up above your own Amdahl ceiling is a measurement defect, not a result (WP-107, red-team S4)"
legacy_seq: 457
---
### A speed-up above your own Amdahl ceiling is a measurement defect, not a result (WP-107, red-team S4)

WP-107 reported 1.11x at 8 threads on a deck whose loop-A fraction it had itself
measured at 7.80 %, i.e. a computed ceiling of **1.08x**. The number was above the
ceiling in the same document and nobody flagged it. Re-measured on an idle box:
1.03x / 1.03x / **0.94x** at 2/4/8 threads — 8 threads is a *regression*.

Two traps in one:

1. **Check every measured speed-up against the ceiling you derived.** Exceeding it
   is proof the measurement is noise (or that the fraction is wrong); it is never
   good news.
2. **"The box was busy, so this is a lower bound" is a sign error as often as not.**
   Contention inflates the serial baseline *and* suppresses the oversubscription
   penalty of the wide thread counts. Whether the bias helps or hurts depends on
   the deck, so the caveat is never a substitute for re-running idle.
