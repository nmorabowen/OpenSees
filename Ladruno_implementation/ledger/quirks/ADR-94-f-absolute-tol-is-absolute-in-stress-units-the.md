---
wp: ADR-94
title: "f_absolute_tol is absolute in stress units — the unit system decides whether strict_convergence refuses (ADR-94 M5)"
legacy_seq: 391
---
### `f_absolute_tol` is absolute in stress units — the unit system decides whether `strict_convergence` refuses (ADR-94 M5)
- **Bites:** the same Mohr-Coulomb problem completes 20/20 in kPa at the default `1e-6` and is refused on step 1 in Pa. `|Phi|` scales with σy (VM), c·cosφ (MC), σci·s^a (HB, ~5 MPa at 50 MPa rock), four decades across the catalogue before units. Tightening to `1e-10` refuses both.
- **Rule:** quote units next to every tolerance; size `f_absolute_tol` to the YF's own strength scale.
- **Status (wp/94c):** FIXED (opt-in). `f_relative_tol` makes the yield tolerance `max(f_absolute_tol, f_relative_tol * yf.strength_scale())`, with `strength_scale()` supplied by the yield function itself (`sqrt(2/3)*sigma_y` for VM, `xi_c` for DP, `c*cos(phi)` for MC/MCTC, `sigma_ci*s^a` for HB, the cohesion-like term for StiffSoil; a YF that declares none returns 0 and is unaffected). Default is **0, i.e. OFF and byte-identical**, so this is a switch you have to reach for; the shipped default is still absolute and still unit-dependent. `tests/test_adr94c_numerics::test_C4_f_relative_tol_makes_the_verdict_unit_independent` runs the same MC problem in two unit systems 1e9 apart and gets 20/20 both times with the option on.
- **The reproducer moved, and that is worth knowing:** wp/94c's shear-slot fix (below) also corrected `Backward_Euler`'s consistency scalar for every Voigt-convention YF, and the MC deck that used to be refused at x1000 now completes even with `f_absolute_tol 0`. The unit gap needed to reproduce M5 grew from x1e3 to x1e9. The defect is unchanged in kind; a faster-converging return map just hid it further out. Do not read "my deck passes now" as "the tolerance is scale-free".
