# WP-144 — LadrunoNORSAND (near-future plan, not started)

**Status: PLANNED.** Recorded 2026-09-27 on the owner's decision after the TIMs F18–F23 program. No code yet. The ND tags are reserved in `LEDGER_implementations.md` (**33023** base, **33024** 3D, **33025** PlaneStrain); they are not in `classTags.h` until the code lands, following the ADR-90 precedent.

## Why

The TIMs strip-footing "wall" was traced to **integration defects**, not to the model: the α escape triggered by the stale α_in, enabled by the stress-only error, the Λ<0-as-elastic branch and the frozen moduli; see [[128_sanisand_ring_trace]] and [[134_sanisand_reference_integrator]]. SAS-ME (WP-129) and CPPM (WP-130) fix that for the calibrated DM04 SANISAND. What remains is a structural weakness: DM04's back-stress α has **no energy function**, so non-negative dissipation is not guaranteed, and its α_in memory and h ∝ 1/((α−α_in):n) singularity are the exact sources of the bugs we fixed.

The model survey ([[_sand_model_survey_2026-09-27]]) ranks **NorSand in the Borja & Andrade (2006, CMAME) form** first for this problem:
- a conservative **hyperelastic** energy;
- **proven non-negative plastic dissipation** when N̄ ≤ N;
- an implicit return map with 3 local unknowns and a closed-form tangent;
- state = (σ, ψ, an image-stress hardening variable), with **no back-stress, no reversal memory and no h singularity**.

**Honest limits:**
- **No true variational (minimum-principle) structure** for a dilatant sand: the admissible non-associated flow gives a non-symmetric tangent. "Admissible, implicit, consistent tangent" is the target, not "variational".
- **No fabric, weak cyclic behaviour.** It is the monotonic/footing model; DM04 remains the cyclic/SSI model, and HySand (WP-145) is the long-term consistent cyclic candidate.

## Design (from the survey; verify before building)

1. **Critical state line: power law, in DM04's form** e_c = e0 − λc (p/Patm)^ξ, NOT NorSand's original logarithmic CSL. The log form sends ψ → −∞ as p' → 0 (unbounded dilatancy at the free surface: the TIMs ring). With DM04's form, e0, λc and ξ carry over from TIMs' calibration.
2. **3-invariant Lode dependence** for plane-strain strength. **Open risk:** Borja & Andrade's dissipation proof is 2-invariant; the 3-invariant extension must be re-derived (survey "could not verify"; adversarial gate).
3. **Hyperelastic** stiffness with G ∝ √p-type pressure dependence from a stored-energy function (check convexity at the ring's stress ratio η ≈ 2.1; Houlsby/Amorosi/Rojas-type energies have limits there).
4. **Low-confinement rule** stated up front: a documented p' floor (Pmin ≤ 0.5 kPa), projected and **counted, not refused**, with the floor-sensitivity report (limit load at F and F/2, accepted if Δ < ~2 %).
5. **Standalone LadrunoJ2-style class** (not the ASDPlasticMaterial3D kit: that has no hyperelastic energy or ψ-driven hardening), with 3D and PlaneStrain wrappers like the SANISAND family; refusal through the WP-99 commit latch from day one (the material checklist's bounded-work and refusal-plumbing rules).

## Plan (the survey's estimate: ~5–6 engineer-weeks)

| Phase | Content | Effort |
|---|---|---|
| P0 | Python **kernel oracle** from the paper (+ the curved CSL + 3-invariant extension), with the dissipation check D ≥ 0 at every step; re-derive the 3-invariant proof | ~1 week (+0.5–1 for the proof) |
| P1 | C++ class + wrappers, implicit return map, closed-form tangent, parser, wire, responses | 1.5–2 weeks |
| P2 | Tests: kernel parity vs the oracle, FD tangent, D ≥ 0 census, byte-identity of everything else, refusal on discarding elements, bounded work | ~1 week |
| P3 | Calibration from TIMs' existing data (the CSL, Mc, G(√p), ν and e transfer; refit χ from peak dilatancy vs ψ, H0/Hy and N/N̄ from the same drained triaxials); single-point gates vs DM04 and PM4Sand; the strip deck (B/8 and B/16) | ~1 week |
| Gate | Adversarial review (new maths); banner, ledgers, guide | — |

**Start condition:** after WP-138 shows what the calibrated DM04 + SAS-ME/CPPM reaches on the strip deck. Build NorSand if DM04 still leaves TIMs short of a defensible limit load, or when TIMs wants a model that is consistent by construction for design use.

## Sources

The survey ([[_sand_model_survey_2026-09-27]]) holds the citations and the E / E-sec / R / I tags. Key: Borja & Andrade 2006 (CMAME); Jefferies 1993 and later (NorSand); the TIMs intake ([[_tims_2d_model_requests_2026-09-25]]).
