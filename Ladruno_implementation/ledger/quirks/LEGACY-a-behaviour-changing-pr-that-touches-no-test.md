---
wp: LEGACY
title: "A behaviour-changing PR that touches no test file leaves the suite asserting the OLD contract, and the failures look like unrelated breakage"
legacy_seq: 292
---
### A behaviour-changing PR that touches no test file leaves the suite asserting the OLD contract, and the failures look like unrelated breakage
- **Bites:** ADR-78 P1 (#730/#731) converted fifteen silent contact degradations into hard aborts. `git show --stat` on both commits: only `CMakeLists.txt`, `LadrunoContactAbort.{cpp,h}`, `LadrunoContactHandler.cpp`. Four tests still asserting the pre-P1 "skip loudly and keep going" contract went red on `ladruno` and stayed there, discovered later by an unrelated regression sweep.
- **Why it evades the usual guards:** the tests fail on `assert ops.analyze(...) == 0` with a `-1`, which reads as "the model stopped converging" — generic breakage — rather than "the contract this test encodes was deliberately replaced". Nothing links the failure back to the PR that caused it. The test NAMES were the only surviving record of the old contract (`_skips`, `_skipped_loudly`, `_inert_`).
- **What catches it:** when a PR changes a contract, grep the suite for tests whose NAME encodes the old one before merging. And when re-greening such tests later, decide per test whether the broken precondition WAS the subject: if it was, invert the assertion (repairing the deck deletes the gate); if it was incidental, repair the deck (inverting turns a physics test into a refusal test). The two are not interchangeable, and choosing by whichever is less typing loses coverage either way.
