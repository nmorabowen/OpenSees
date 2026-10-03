---
wp: WP-130
title: "Two WPs indexed their sendSelf/recvSelf blocks at the SAME offset (35 + LMS_COUNT) -- a textual merge conflicted only in sendSelf, recvSelf AUTO-MERGED (WP-130…"
legacy_seq: 526
---
### Two WPs indexed their sendSelf/recvSelf blocks at the SAME offset (35 + LMS_COUNT) -- a textual merge conflicted only in sendSelf, recvSelf AUTO-MERGED (WP-130 x WP-129, review #868 item 1)
- **Bites:** WP-129 (SAS-ME options + sasStats) and WP-130 (CPPM options) each appended their block
  "after the census" at `35 + LMS_COUNT`. git conflicted in sendSelf and on the vector size, but
  the two recvSelf READ blocks merged cleanly and read the same slots twice: a restored or
  MP-received IntScheme-2 point got mCPPMOnFail = -1 (refuse) and mCPPMHalvings = 0 from WP-129's
  default SAS values, unclamped. Nothing failed loudly.
- **Fixed (WP-130 merge, #868):** the layout is derived from named constants in
  `LadrunoSANISAND.h` (`LWIRE_CENSUS`, `LWIRE_CPPM`, `LWIRE_CPPM_N`, `LWIRE_SAS`,
  `LWIRE_SAS_OPT_N`, `LWIRE_SIZE`), documented above `LadrunoSANISAND::sendSelf`; received CPPM
  options are clamped (setLadrunoCPPMOptions' rule). Pinned by
  `test_wire_round_trip_both_blocks_after_the_129_merge` (non-default options of BOTH blocks saved,
  restored into a DEFAULT-built skeleton; options, both censuses and the next two steps exact).
  **Rule: a new wire block gets its own named offset, never "35 + LMS_COUNT + k".**
