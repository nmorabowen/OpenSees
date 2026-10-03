---
wp: LEGACY
title: "contact_dump byte-identity is a WITHIN-TOOLCHAIN-SESSION observable -- cross-session hashes do not compare"
legacy_seq: 315
---
### `contact_dump` byte-identity is a WITHIN-TOOLCHAIN-SESSION observable -- cross-session hashes do not compare
- **Bites:** anyone diffing a `contact_dump` hash recorded by an earlier session/build machine against a fresh build and freezing a phase over the mismatch.
- **Why:** T2 recorded `da73b6f8...7782` at source tree bce163ccf; a fresh full build of the IDENTICAL tree (harness unchanged since T0) hashes `B0F8F770...81E4`, twice, byte-identical within the session. The T2 artifacts no longer exist so the byte-level cause cannot be narrowed below {independent MUMPS/conan rebuild, toolchain drift between builds}; either way the gate's semantics were always same-session pre/post-change (the ADR's "both ladrunoBuild stamps recorded" clause).
- **Workaround/status (2026-08-18, ADR-85 T3):** every phase RE-CAPTURES its baseline at its own tip with its own binary (twice, hashes must match) before touching C++; the recorded hash is scoped to the session that measured it.
