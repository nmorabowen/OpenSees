---
wp: LEGACY
title: "A central-difference FD tangent of ManzariDafalias straddles the alpha_in reversal reset — use a one-sided, direction-safe difference"
legacy_seq: 468
---
### A central-difference FD tangent of `ManzariDafalias` straddles the `alpha_in` reversal reset — use a one-sided, direction-safe difference
- **Bites:** `(sigma(eps + h e_j) - sigma(eps - h e_j)) / 2h` at a plastic state gives off-diagonal entries 40-100 % off any analytic tangent, at every `h`, and looks like a tangent bug (WP-110 F15a first pass).
- **Why:** `integrate()` resets `alpha_in := alpha_n` when `(alpha_n - alpha_in):Ce:d_eps < 0`. Whether that fires depends on the SIGN of the probe increment, so `+h` and `-h` land on two different internal states (read back: `alpha_in(+h) != alpha_in(-h)` in all six directions) and the central difference differences across a kink, not along one branch.
- **Workaround/status:** per direction, pick the sign whose run leaves `alpha_in` equal to the committed base state, and difference one-sidedly against the base stress (h = 1e-8 matched the corrected tangent to 0.02-0.17 %). Also make the last load increment tiny, because both integrators take G and K at the START of an increment while the FD samples the committed state. Implemented in `tests/test_manzari_ep_tangent_gate.py`.
