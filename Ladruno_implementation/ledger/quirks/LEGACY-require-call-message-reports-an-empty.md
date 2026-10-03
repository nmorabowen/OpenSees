---
wp: LEGACY
title: "require(call(...), \"...\" + message) reports an EMPTY diagnostic — argument evaluation order is unspecified"
legacy_seq: 236
---
### `require(call(...), "..." + message)` reports an EMPTY diagnostic — argument evaluation order is unspecified
- **Bites:** every CMS standalone check uses `require(bool, std::string)` and the natural idiom `require(doThing(..., message) == 0, "doThing failed: " + message)`. C++ does not specify the order in which function arguments are evaluated, so MSVC builds the message string **before** running the call — capturing `message` while it is still empty. The failure prints `FAIL: doThing failed:` with nothing after the colon, which reads like the callee returned no diagnostic and sends you looking in the wrong place. Cost the P3d work a full debug cycle: the real message was `invalid distributed hierarchy input on at least one rank`, which points straight at the cause.
- **Why it is not a correctness bug:** the *condition* is still evaluated correctly, so nothing passes that should fail. Only the diagnostic is lost — and only on the failure path, which is exactly when you need it.
- **Workaround/status:** call first, store the bool, then `require(ok, "... " + message)`. Swept 2026-07-26: `REQUIRE_CALL(status, text)` evaluates the call first and only then builds the diagnostic; 11 sites converted across `assembly` (4), `lanczos` (4), `mumps` (2) and `topology` (1). Sites that already computed the bool into a variable were left alone -- they never had the problem. **Use `REQUIRE_CALL` for any new `require(someCall(...), "..." + message)`.** *2026-07-26 (ADR-1000 P3d).*
