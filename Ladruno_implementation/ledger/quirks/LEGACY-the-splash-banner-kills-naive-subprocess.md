---
wp: LEGACY
title: "The splash banner kills naive subprocess capture on Windows — UnicodeDecodeError (cp1252) at ~offset 3209, and a broad except turns it into a silent \"build unk…"
legacy_seq: 277
---
### The splash banner kills naive `subprocess` capture on Windows — `UnicodeDecodeError` (cp1252) at ~offset 3209, and a broad `except` turns it into a silent "build unknown"
- **Bites:** `subprocess.run([opensees...], capture_output=True, text=True)` on a cp1252-locale Windows box raises `UnicodeDecodeError` decoding the banner's UTF-8 box-drawing glyphs; any wrapper with a broad `except` around it silently degrades to "hash unknown". Measured cost: the TIMs provenance harness lost two runs to exactly this while pinning engine identity after the T1 probe incident (2026-08-10) — the one context where you MUST capture the banner, because the `Ladruno OpenSees build: <hash>` line lives inside it.
- **Why:** the banner is emitted as UTF-8; Python's text-mode subprocess decodes the child's stdout with the locale default (`cp1252` on most Windows), which has no mapping for the box-drawing bytes.
- **Rule:** capture engine output with `encoding="utf-8", errors="replace"`, never bare `text=True`. `LADRUNO_OPENSEES_QUIET=1` sidesteps the decode but ALSO suppresses the build-hash line — so it is the wrong tool when the capture is FOR provenance.
- **Workaround/status:** documented; the real fix is a machine-readable build-stamp query (Tcl + Python) so provenance never scrapes the banner — proposed as follow-up. *Learned 2026-08-10 (TIMs T1 provenance incident).*
