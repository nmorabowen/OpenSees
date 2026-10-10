---
wp: WP-130
title: "A WP-130 flag given with its DEFAULT value on a deck where it cannot act was accepted; refuse still ground the full ladder (WP-130 review r1)"
legacy_seq: 518
---
### A WP-130 flag given with its DEFAULT value on a deck where it cannot act was accepted; `refuse` still ground the full ladder (WP-130 review r1)
- **Fixed (WP-130, #868):** every flag is validated on having been GIVEN (`-cppmOnFail explicit`
  on IntScheme 1, `-cppmTangent` under TanType 0/1, `-meFallback off` without a cap: refused).
  `-cppmOnFail refuse` without `-cppmHalvings` bounds the ladder at 3 (echoed): vanilla's 9 cost up
  to 1.1 s per refused update. After a RESCUED ME cap hit `lastCapHit` is 2 (capHits still counts
  it). The WP-99 latch warning names its cause (companion vs CPPM refusal) and says it is sticky
  until `reset()`.
