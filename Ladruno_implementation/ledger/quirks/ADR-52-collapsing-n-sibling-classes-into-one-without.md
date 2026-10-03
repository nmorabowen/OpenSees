---
wp: ADR-52
title: "Collapsing N sibling classes into one without breaking serialization (ADR-52 W1-E2)"
legacy_seq: 135
---
## Collapsing N sibling classes into one without breaking serialization (ADR-52 W1-E2)

**Pattern (reusable):** when you fold a family of integrator/element/material *classTags* into ONE
class selected by flags, the trap is `sendSelf`/`recvSelf` and the object brokers — a saved DB or an
MPI rank reconstructs an object by its **classTag**, via `FEM_ObjectBroker::getNewX(classTag)`, then
calls `recvSelf` on it. If you delete the retired tags or point them at a default-constructed unified
object, the flag state is lost and the reconstructed object behaves wrong (silently — no error).

**What works (ExplicitBathe 6→1, #419):**
1. **Keep all N retired `#define`s** in `classTags.h` (annotate DEPRECATED; never free/reuse — a
   future feature grabbing 33009 would alias a saved model). The collapse frees ZERO tags.
2. **Derive the classTag from the flag combo** in the ctor (`tagForFlags(lnvd,sms,consistent)`),
   passing it to the base `TransientIntegrator(classTag)`. So `getClassTag()` of a `-lnvd -sms` object
   IS the legacy `ExplicitBatheLNVDSMS` tag ⇒ `sendSelf` writes the right tag with no extra field.
3. **Route every retired tag through one static factory** `X::makeForBroker(classTag)` that decodes
   the flags from the tag (`flagsForTag`) and constructs the unified object pre-set. Both brokers
   (`FEM_ObjectBrokerAllClasses.cpp` + `runtime/.../TclPackageClassBroker.cpp`) use a 6-way
   case-fallthrough → one `return X::makeForBroker(classTag);`.
4. **Send a FIXED-SIZE param superset** in `sendSelf`/`recvSelf` (union of all variants' payloads).
   The flags do NOT travel in the payload — they are implied by the classTag the broker already used.
   recvSelf only fills the numerics; it must reset transient bookkeeping (the SMS injected-map,
   energy-registry, warn-once flags).
- **Validity caveat:** the flag↔tag map must be a true bijection over the *valid* combos. ExplicitBathe
  has 6 valid combos (`consistent` implies `sms`), matching the 6 legacy tags exactly — so `tagForFlags`
  defends `if (!sms) consistent=false`. If your flag space has combos with NO legacy tag, you need a NEW
  tag for them (and a manifest row), not a silent collapse.
- **Architecture gotcha that fed this:** the "6 siblings" were really a **2-base × 3-mode lattice**
  (`ExplicitBathe`/`ExplicitBatheLNVD` bases, each ×{none,sms-lumped,sms-consistent}). The two bases'
  `newStep`/`update`/`commit` were **byte-identical**; only LNVD's `formUnbalance` override + the SMS
  `domainChanged` injection differed. Diff the candidate classes' hot paths BEFORE assuming a merge is
  behavior-preserving — here it was, so the byte-identity tests (assert on disp/vel/accel, not stderr)
  hold. Keep each retired OPS_ parser verbatim (exact historical positional grammar) one release; each
  just constructs the unified class with fixed flags.
