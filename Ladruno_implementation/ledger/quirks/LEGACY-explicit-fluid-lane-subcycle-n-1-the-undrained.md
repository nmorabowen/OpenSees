---
wp: LEGACY
title: "Explicit-fluid lane + -subcycle N>1: the UNDRAINED CFL binds the SYNC interval N*dt — the implicit lane's large-N freedom does NOT transfer"
legacy_seq: 199
---
### Explicit-fluid lane + `-subcycle N>1`: the UNDRAINED CFL binds the SYNC interval N*dt — the implicit lane's large-N freedom does NOT transfer
- **Bites:** carrying the E7.3a intuition ("all N <= 50 stable, error ~ N^1.2") from the implicit `-fsL zero` lane to `-fluidUpdate explicit`: N=4 at 0.4x the undrained pencil (sync interval 1.6x) diverges in ~400 steps on the e72 column — where the SAME dt at N=1 is bounded for 60k steps. Measured and toy-twinned (C++ step 395 / toy step 408).
- **Why:** with the fluid explicit, the coupled staggered stability is set by the frozen-force interval = the sync interval N*dt (the undrained stiffening acts once per sync); the implicit-at-commit fluid absorbed exactly that stiffening, which is where E7.3a's freedom came from. Subcycling also has no purpose on the explicit lane — the fluid step is an axpy, there is no solve to amortize.
- **Workaround/status (2026-07-19, ADR-73 P3b):** `-subcycle auto` under `-fluidUpdate explicit` resolves N=1 with a notice; manual N>1 prints a loud one-time sync-CFL warning (keep N*dt within the N=1 margin). Battery s5 pins the expected-diverge demo, the bounded sync-0.4x leg, and the auto->1 notice. ADR §12 P3b item 9.
