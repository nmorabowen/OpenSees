---
wp: LEGACY
title: "ManzariDafalias::ModifiedEuler's low-p clamp wrote a p that contradicted the stress it had just built — and got away with it because the store is DEAD"
legacy_seq: 336
---
### `ManzariDafalias::ModifiedEuler`'s low-`p` clamp wrote a `p` that contradicted the stress it had just built — and got away with it because the store is DEAD

- **Symptom:** none, ever. That is the whole point of the entry.
- **What it was:** the clamp is `if (p < m_Pmin + m_Presidual) { NextStress = GetDevPart(NextStress) + m_Pmin*mI1; p = m_Pmin; }`. Everywhere else in that function `p` denotes `one3*GetTrace(sigma) + m_Presidual`, and the rebuilt tensor has `one3*GetTrace(NextStress) == m_Pmin`, so the consistent scalar is `m_Pmin + m_Presidual` — **1.0201 against the 0.0101 that was written, a factor of 101** at the vanilla constants. `Stress_Correction()` writes `m_Pmin + m_Presidual` at its analogous clamp; `RungeKutta45()` recomputes from the rebuilt stress. `ModifiedEuler` was the only one of the three that disagreed with itself.
- **Why nobody ever saw it:** the store is dead. `p` is a local, the `while (T < 1.0)` loop immediately below always runs, and its first contact with `p` is an unconditional recompute — nothing reads the clamped value. A 7-leg fingerprint (including three legs with `-Pmin` raised until the clamp branch is genuinely entered, one of which completes 1200 steps) is **byte-identical** before and after the repair.
- **The lesson worth keeping:** *dead stores hide wrong values, they do not make them right.* This block is exactly where PR-2 added the clamp's user-facing diagnostic, and it is exactly where the next reader of `p` will be added. A latent wrong value sitting one inserted line above a recompute is a trap primed for whoever touches the block next — and it would have surfaced as a factor-of-101 error in a *diagnostic*, i.e. as a lie about the model rather than a crash.
- **Corollary for auditing this family:** when you find two sites that "disagree", check reachability before you rank them. Here the disagreement was real and the blast radius was zero; in [[LEDGER_vanilla_files]] commit 4 of the same PR the reverse held (three identical-looking sites, load-bearing at exactly one).
- **Not made uniform, on purpose:** `Stress_Correction()` also collapses the deviator and lands one `m_Presidual` further inside the floor. Those are recovery-strategy choices on a LIVE return value, not internal contradictions — ADR-86 D9 opinion, not error.
- **Learned:** 2026-08-27.
