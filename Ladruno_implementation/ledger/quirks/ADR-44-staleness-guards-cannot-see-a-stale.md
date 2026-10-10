---
wp: ADR-44
title: "Staleness guards cannot see a stale DomainModalProperties that reproduces the SAME spectrum — write guard tests with unique stiffnesses (ADR 44 P3)"
legacy_seq: 176
---
## Staleness guards cannot see a stale `DomainModalProperties` that reproduces the SAME spectrum — write guard tests with unique stiffnesses (ADR 44 P3)

`DomainModalProperties` survives `wipe()` (the [[#`wipe()` does NOT recreate the
Domain — new domain-level state MUST be reset in `Domain::clearAll()` (ADR 46 P1)|
clearAll leak]] family). The P1a/P2/P3 staleness guards compare eigenvalue count +
element-wise values between the Domain and the snapshot — so a
`wipe(); rebuild-IDENTICAL-model; eigen` sequence leaves a stale-but-equal snapshot
the guard legitimately CANNOT distinguish (same spectrum ⇒ same Γ/Vscale up to sign
⇒ numerically the same answer, so it is also harmless). The trap is in TESTS: a
`guard_no_modalproperties` pytest that rebuilds the same `m,k` as any earlier test
in the file will NOT raise. Give guard-test models a stiffness unique within the
file (`test_ladrunoRandomResponse.py` uses k=512 for exactly this reason).
