---
wp: WP-133
title: "pAtm is a STATIC member of PDMY01/02/03 — the last material created sets the atmospheric pressure for every material of that class (WP-133, found by reading)"
legacy_seq: 508
---
### `pAtm` is a STATIC member of PDMY01/02/03 — the last material created sets the atmospheric pressure for every material of that class (WP-133, found by reading)
- **Bites:** each constructor ends with `pAtm = atm;` on `static double pAtm`. Two PDMY03 materials with different `$pa` (e.g. one in kPa, one in Pa, or a sensitivity study) silently share the last one's value in every pressure normalisation (`isCriticalState`, the contraction/dilation `(p/pa)` factors). Found by reading during WP-133; not exercised by a test.
- **Workaround/status:** keep one `$pa` per process for each PDMY class (the usual case). Not fixed (vanilla; a per-material array would follow the three-edit rule of the entry above).
