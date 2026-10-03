---
wp: LEGACY
title: "ManzariDafalias's elastic G uses the INITIAL void ratio m_e_init, not the current one — and the correct line is sitting right there, commented out"
legacy_seq: 338
---
### `ManzariDafalias`'s elastic `G` uses the INITIAL void ratio `m_e_init`, not the current one — and the correct line is sitting right there, commented out

- **What the paper says.** Dafalias & Manzari (2004) p.623: `G = G0*p_at*(2.97-e)^2/(1+e)*(p/p_at)^0.5`, with `e` the **current** void ratio.
- **What the code does.** All three `ManzariDafalias::GetElasticModuli` overloads use `m_e_init` — six lines, two per overload (`mElastFlag == 0` and `else`). The current void ratio is passed in as the parameter `en` and was **never read**. The corrected expression is present in the file, inside the disabled `/* ... */` block three lines above the first site. `b0` elsewhere in the same file uses the current `e`, and `SAniSandMS` uses `en`. So the file disagrees with itself, with its own commented-out reference, and with its sibling.
- **Why it is NOT simply fixed.** It moves a **calibrated** quantity. Any `G0` fitted against the frozen form absorbs the error (the reporting project's `G0 = 264.32` did), so "correcting" `G` without refitting `G0` just preserves one error by means of another — and does it silently, to every existing deck. ADR-86 D9: a calibration-moving error takes a **flag seam**, not a bugfix.
- **How it ships (ADR-86 PR-2 commit 7).** Protected `bool mUseCurrentVoidRatioInG` on the base, false in all four constructors, read at each site as `const double& eG = mUseCurrentVoidRatioInG ? en : m_e_init;`. Vanilla is bit-identical (proven on a low-`p` triaxial fingerprint AND a Ramberg-Osgood Tcl deck), and **nothing is wired to the flag yet** — turning it on is a separate, deliberate decision that must come with a `G0` refit.
- **The trap next to it.** `ManzariDafaliasRO` shadows two of those three overloads with identical signatures. Adding a `bool` member is safe; adding `virtual` would convert the shadows into overrides and start running Ramberg-Osgood elasticity inside every base integrator. See the dedicated `virtual` entry above — the seam pattern exists precisely because `virtual` is off the table here.
- **Useful A/B fact if you ever test this.** RO does **not** shadow the `(sigma, en, K, G, const double& D)` overload, and `commitState:464` calls it on every commit — so a plain RO deck exercises RO's own elasticity *and* a base `GetElasticModuli` body in the same run. That is what makes an RO deck a real test of a base-side edit rather than a null one.
- **Learned:** 2026-08-27.
