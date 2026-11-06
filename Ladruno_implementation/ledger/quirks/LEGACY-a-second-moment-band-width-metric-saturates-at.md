---
wp: LEGACY
title: "A second-moment band-width metric SATURATES at the specimen length — w2 near L means \"no band\", not \"a wide band\""
legacy_seq: 351
---
### A second-moment band-width metric SATURATES at the specimen length — `w2` near `L` means "no band", not "a wide band"
- **Bites:** you adopt a threshold-free band width `w2 = sqrt(12*Var)` over the plastic-strain profile (calibrated so a one-element band reads exactly `h` and a k-element top hat reads exactly `k*h`), sweep the regularization parameter, and read off a beautifully monotone `w2` = 3.5 -> 42 -> 91 mm on a 100 mm bar. The last two are not band widths. A **uniform** profile returns `w2 = L` exactly, so `w2` is bounded in `[h, L]` **by construction**, and the metric reports "the whole specimen is yielding" and "there is a 91 mm band" with the same number.
- **Why:** the second moment of a distribution supported on `[0, L]` is maximised by the uniform distribution, and the `sqrt(12*Var)` normalisation was chosen precisely so the uniform case maps to `L`. Convergence of `w2` toward `L` under refinement therefore *looks* like convergence of a width when it is convergence to a homogeneous, non-localizing solution.
- **Workaround:** always report `w2/h` (the band-resolution floor — below ~3 the "band" is two elements and the number is meaningless) **and** `w2/L` (above ~0.5 there is no band); gate the declared operating point on `w2/h >= 3`. Report a threshold metric (FWHM) alongside and expect disagreement — on the same runs the threshold count and the FWHM differed by up to 40x, and the FWHM failed to converge at parameter points where `w2` converged. Learned 2026-09-05, ADR-90 P0b leg (c), [[_adr90_p0b_results]] §3.
