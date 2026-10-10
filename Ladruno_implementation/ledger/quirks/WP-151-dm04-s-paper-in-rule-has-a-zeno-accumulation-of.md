---
wp: WP-151
title: "DM04's paper α_in rule has a Zeno accumulation of re-seats near the peak: b:n → 0⁺ and |dα/dt| → ∞ in FINITE pseudo-time — that, not H ≤ 0 alone, is the SAS-ME…"
legacy_seq: 530
---
### DM04's paper α_in rule has a Zeno accumulation of re-seats near the peak: b:n → 0⁺ and |dα/dt| → ∞ in FINITE pseudo-time — that, not H ≤ 0 alone, is the SAS-ME `loadingNonPosH` wall (WP-151)
- **Bites:** The WP-138 footing (E_B, E_D, E_B16) walls on `loadingNonPosH`, with 10 M α_in re-seats and 88 M rejected reversals in E_B. The exact oracle (WP-134, Radau) shows why:
  - After a re-seat, h = ∞ makes α slide along b.
  - With b nearly normal to n, that slide rotates n on the thin cone (√(2/3)m = 0.004 for m = 0.005). (α−α_in):n turns negative again, and the next re-seat fires.
  - The intervals shrink geometrically (∝ (b:n)², ×0.04 per re-seat) while b:n → 0 from ABOVE and |dα/dt| ∝ 1/b:n → ∞. H ≤ 0 is only the b:n < 0 exit.
  - On the five real wall states × 64 trials, today's SAS-ME refuses exactly the 102/320 trials the exact oracle cannot integrate, one-to-one.
  - At c = 0.71 every refuser has |α−α_in| of about one cone radius (a re-seat a moment ago) and ρ_b > 1 > ρ_α (n at an extension-side Lode angle).
  - A c = 0.80 footing walls too, with compression-side refusers (cos3θ ≈ +0.65) and the same sequence. The singular set belongs to DM04, not to the calibration.
  - Chen, Ghorbani, Zhang & Kodikara (2022, §3.9.1; verified in Chen's published-works thesis, doi:10.26180/23639730.v1, Ch. 3) report a SANISAND04 plane-strain footing on loose sand that aborts when (α − α_in):n drops suddenly to 0, and worse with finer steps or a tighter tolerance. It is the same singular factor; b:n is not analysed there.
- **Rule:** Neither piece alone cures it:
  - a floor on h alone: 97/320 still chatter, and the C++ discretizes that into `-maxSubsteps`;
  - a floor gated on b:n ≤ 0: it misses the b:n → 0⁺ side (102/320);
  - a re-seat threshold alone: h stays 1e10 in its band (102/320).

  Use the everywhere floor AND the hysteresis together: two re-seats then need a finite α travel, so they
  cannot accumulate.
  - At footing scale add the softening CAP too. Past the old wall, b:n < 0 post-peak points near a reversal
    make even the floored h drive H ≤ 0.
  - Floor + hysteresis without the cap walls at s/B 0.0525. The floor alone and the hysteresis alone wall
    EARLIER than DM04.
- **Workaround/status:** WP-151's opt-in flags `-sasHFloor 1 -sasReseatHyst 1 -sasSoftCap 0.5` (SAS-ME only; default OFF and byte-identical) give 0/320.
  - Small-strain cost: about 1 % softer at γ ~ 1e-5; take G0 from the elastic range (memo §6.4).
  - Full set on the footing: 0 NonPosH at s/B 0.054, still hardening, where E_B walls at 0.0508.
  - Owner (relayed 2026-09-28): merge with the full set recommended. κ stays an owner/TIMs choice, and the
    κ 0.25 / 0.75 legs are pending.
  - On the footing (Esmeralda, 2026-09-28) they pass E_B's onset with 0 `loadingNonPosH` on B/8, B/16 and B/4. The B/4 leg reaches s/B 0.127, 2.5× E_B's wall.
  - q–s stays within ±0.22 % of E_B below s/B 0.03.
  - Floor alone gives `maxSubsteps` instead; hysteresis alone keeps NonPosH, as predicted.
  - **Do not "fix" it by changing the Lode parameter c:** a c = 0.80 footing walls too (s/B 0.048).
  - [[151_sanisand_reseat_singularity]] §2.5, §9.
