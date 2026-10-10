---
wp: LEGACY
title: "HB_sigma_ci is not a parameter; it is HB_sigci"
legacy_seq: 413
---
### `HB_sigma_ci` is not a parameter; it is `HB_sigci`

ADR-97 P1 and P2 both wrote `HB_sigma_ci` in their Hoek-Brown refusal tests. Under
the ADR-94 contract a missing model parameter is an ERROR, so those decks were
rejected for the wrong reason and the refusal assertions never exercised the family
gate at all. A refusal test that is not ALSO checked in the positive direction (the
same deck must CONSTRUCT under an integrator that does support it) cannot tell the
two apart. Every gate-6 row in `tests/test_adr97_p3_hoekbrown.py` is checked both
ways for exactly this reason.
