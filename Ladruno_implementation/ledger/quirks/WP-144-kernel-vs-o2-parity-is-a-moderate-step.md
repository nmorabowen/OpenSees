---
wp: WP-144
title: "Kernel-vs-O2 PARITY is a MODERATE-STEP statement: at ~1e-2 strain per step both codes decide on round-off, and the iteration-count equality gate holds on the F…"
date: 2026-11-06
---
### Kernel-vs-O2 PARITY is a MODERATE-STEP statement: at ~1e-2 strain per step both codes decide on round-off, and the iteration-count equality gate holds on the FIXED paths only (WP-144)
- **Bites:** the kernel reproduces O2's algorithm expression by expression, so state, stress and tangent agree to ~1e-12..1e-14 on the P1a path set. That statement has a domain.
  - **Step size.** Parity holds for steps up to about 3e-3 strain. At 1e-2 a step is decided by quantities that are tied to round-off in BOTH codes: the line-search accept/backtrack, the pi_i scan bracket, the substep ladder. An ulp-level difference in the eigen-decomposition (LAPACK `eigh` in O2, Jacobi in the kernel) can then take the two codes down different, equally valid branches (substep depth, a backtrack), and the states stop agreeing to 1e-14.
  - **Iteration counts.** The gate "local / pi_i iteration counts equal" is meaningful only where the path is fixed by construction (the P1a fixed paths, moderate steps). On randomised or large-step paths a count differs for a non-bug reason (a branch tied at round-off), so counting mismatches there would be noise. It is reported, not gated.
- **Rule:**
  - State the domain next to every parity number (moderate step, <= ~3e-3). Do not stretch a parity gate to 1e-2 steps by loosening the tolerance: a loose gate stops seeing the real defects (see the mutant gate).
  - For large steps compare against the TRUTH (O1, the rate oracle) with a convergence argument, not kernel against O2 to 1e-10.
  - Gate iteration-count equality on the fixed paths only; on the others print it.
- **Workaround/status:** ✅ Documented and in force in `Ladruno_files/testbed/norsand_oracle/kernel_parity/test_kernel_parity.py` (WP-144 P1a, 2026-10-01); a refusal PATH (kernel refuses where O2 refuses) is compared by reason, not by iteration count.
