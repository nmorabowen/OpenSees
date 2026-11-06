---
wp: WP-127
title: "The TIMs ring-point CSVs carry the INTERNAL, compression-POSITIVE mSigma, although their README says \"compression negative\" (finding A, WP-127)"
legacy_seq: 482
---
### The TIMs ring-point CSVs carry the INTERNAL, compression-POSITIVE `mSigma`, although their README says "compression negative" (finding A, WP-127)
- **Bites:** `_tims_2d_model_requests_2026-09-25/ring_points_b{8,16}.csv`: on all 80 rows `p_kPa == +tr(sigma)/3` and every normal stress is >= 0. `ManzariDafalias` stores `mSigma` compression-positive; only the wrappers' `getStress()` flips it (`eleResponse ... stress` is tension-positive). A dump taken from internal members (or a post-processor that negated twice) looks exactly like the README's claim until you check `p` against the trace. `alpha`, `alpha_in` and `z` are ratios and are NOT flipped by the wrappers in either direction. b8 row 1859/2 also carries `tr(alpha) = 2.3e-3` (the rest ~1e-10).
- **Workaround/status:** `ladrunoSANISANDReplay` has **no default convention** — `-convention compressionPositive|tensionPositive` is required — and projects alpha/alpha_in/z to deviatoric with a warning above round-off (the pre-projection traces are returned). `sanisand_replay.check_sign_convention()` verifies a dump row. Tell the act (WP-127 plan, reply doc).
