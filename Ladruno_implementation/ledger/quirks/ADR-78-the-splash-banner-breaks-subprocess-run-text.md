---
wp: ADR-78
title: "The splash banner breaks subprocess.run(text=True) test harnesses on cp1252 consoles (ADR 78)"
legacy_seq: 242
---
## The splash banner breaks `subprocess.run(text=True)` test harnesses on cp1252 consoles (ADR 78)

- **Bites:** tests that spawn a python child with `capture_output=True, text=True`
  and import the engine (`test_ladruno_up_element_th.py` winding gate,
  `test_ladruno_up_mp_smoke.py` serial roundtrip) die with
  `TypeError: argument of type 'NoneType' is not iterable` on a cp1252-locale
  Windows box: the child prints the splash banner, whose UTF-8 art/feature text
  contains bytes cp1252 cannot map (e.g. the superscript minus in `F0⁻¹`,
  byte 0x81, present since the StagedDefGrad banner line), the reader thread
  raises `UnicodeDecodeError`, and `proc.stdout` comes back `None`. Looks like a
  physics regression; is an encoding trap. Predates ADR 78 (verified: the byte
  is in the committed `banner_features.txt`).
- **Workaround/status:** set `LADRUNO_OPENSEES_QUIET=1` in the child env (both
  tests pass then), or pass `encoding="utf-8", errors="replace"` to
  `subprocess.run`. Proper fix (harden the harnesses) spun off as its own task.
  Keep NEW banner-feature lines ASCII-only. *2026-07-28 (ADR 78).*
