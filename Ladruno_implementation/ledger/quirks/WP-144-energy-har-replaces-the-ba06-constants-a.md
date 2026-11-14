---
wp: WP-144
title: "-energy HAR REPLACES the BA06 constants: a constant of the other energy is REFUSED, never ignored; -p_a is ONE flag; the reference pressure moves from |p0| to…"
date: 2026-11-06
---
### `-energy HAR` REPLACES the BA06 constants: a constant of the other energy is REFUSED, never ignored; `-p_a` is ONE flag; the reference pressure moves from `|p0|` to `p_a` (WP-144 round 3, sheet 2.4)
- **Bites:** (1) HAR replaces `alpha0`; `-p0 -kappa_hat -eps_v0 -mu0 -alpha0` given with `-energy HAR` would be read by nothing, so a deck that carries both silently runs a different model from the one it documents. The shell refuses them (code 201), and `-k -g -n -G0 -nu -e_ref` with BA06 (202), both stiffness forms (203), one of a pair (204), an out-of-range DM04 input (206); the kernel's `validate()` refuses `k <= 0`, `g <= 0`, `n` outside `[0, 1)`, `p_a <= 0`, `p_min < 0` (21 - 25). (2) `p_a` is one flag shared by the HAR energy and the fork CSL: with `-energy HAR -csl fork` the energy's and the CSL's `p_a` are the same number (TIMs: 101 kPa, the campaign `Patm 101`; 101.325 is withdrawn). HAR requires `-p_a` (no unit-blind default). Under BA06 with the paper CSL it is unused (a note is printed). (3) Every scaling that used `|p0|` (`F_tol`, `r4`, the F0 band of the initial state, the `p_min` default) uses `p_ref`: `|p0|` (BA06) or `p_a` (HAR, with `p0 := -p_a`).
  - A FD tangent test cannot see a HAR deck silently running BA06 (a consistent swap passes its own FD check): the kill is an expected value (K1 13.1h, 13.11, 13.13, 13.14) or parity with O2 running the HAR law.
  - The DM04 mapping (`-G0 -nu [-e_ref]`, e_ref defaults to `v0 - 1`) is exact on the isotropic axis only; off the axis HAR's moduli differ from DM04's (G ratio 1.61 at `eta = M`). It is a convenience of the parser: the material stores `k` and `g`.
- **Rule:** one flag, one meaning; any parameter a model does not read must be refused when given, not accepted and dropped. The echo prints the energy, `p_ref` and the floor so a deck shows what it runs.
- **Workaround/status:** IMPLEMENTED (WP-144 round 3, `parseLadrunoNorSand` + kernel `validate`). Final gate pending.
