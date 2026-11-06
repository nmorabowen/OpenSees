---
wp: WP-107
title: "A process-wide static bool \"warn once\" latch makes every model after the first MUTE (WP-107, red-team S3)"
legacy_seq: 456
---
### A process-wide `static bool` "warn once" latch makes every model after the first MUTE (WP-107, red-team S3)

The pattern is everywhere in this codebase and it is wrong whenever the message
describes a **per-model** decision rather than a per-process one:

```cpp
static bool warned = false;
if (!warned) { warned = true; opserr << "WARNING ..."; }
```

WP-107's threaded element loop re-audits the domain on *every* `Domain::update()`,
so `wipe` + rebuild, a runtime `element`, and `remove element` were all handled
correctly — but with three such latches, only the FIRST model in the process ever
said what it had decided. Measured: clean deck announces THREADED; `wipe` + a deck
with a bad element 7 refuses correctly and says so; `wipe` + a deck with a bad
element 33 refuses **silently**; `wipe` + a clean deck threads **silently**. A run
that quietly went serial and a run that stayed threaded then look identical, which
is the exact confusion the message exists to prevent. A one-process-per-run bench
driver hides it; pytest, apeGmsh and any in-process parameter study do not.

**Rule:** latch per *(owning object, outcome, identifying tag, parameter)*, not per
process, and give the owning object a generation counter that its mutators bump
(`Domain::addElement` / `removeElement` / `clearAll` here). Keep the steady-state
quiet so the message cannot become per-iteration spam — both failure modes are
real, and the test file asserts both directions.
