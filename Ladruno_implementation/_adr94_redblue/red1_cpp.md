# ADR-94 R2 — RED lane 1 (C++ state machine / API contract)

Read-only, build `52314165a`. No C++ edited. Line numbers are
`SRC/material/nD/ASDPlasticMaterial3D/ASDPlasticMaterial3D.h` unless noted.
Reproducer (structural, Zone-A, <1s): `tests/test_adr94_redblue_cpp.py`.

## Verdict

The state machine is **unsafe to rely on as shipped**: the class-static
sharing R1-A found in `Stiffness` (H1) is not an isolated bug, it is a
*pattern* repeated at three more layers — getters, model-parameter I/O, and
the YF/PF functor layer — plus two independent silent-failure sites R1
did not name (`getInitialTangent()`'s side effect, `sendSelf`/`recvSelf`
claiming success while doing nothing).

## Findings

**F1 (new, blocker-adjacent) — `getInitialTangent()` is a getter with a
global side effect that clobbers `Stiffness`.** Lines 708–718: it computes
`Eelastic` from `CommitStress` and does `Stiffness = Eelastic;` before
copying to the return buffer. `Stiffness` is the same class-static H1
already showed is shared by every GP/element/tag of one YF×PF×EL combo. Any
code path that calls `getInitialTangent()` (initial-stiffness iteration,
`-initial` algorithms, diagnostics) on *any* instance overwrites the tangent
every other instance's later plain `getTangent()` (698–705, a pure copy)
will read — independent of, and in addition to, H1's "last GP integrated"
mechanism. A "read-only" API call mutates shared state.

**F2 (new, broader-than-H1 threading blocker) — the YF/PF static return
buffers are scoped to the *functor type*, not the material combo.** Every
`YieldFunctionType`/`PlasticFlowType` returns its vector results through a
private `static VoigtVector vv_out` (confirmed in all 7 PF headers and
`DruckerPrager_YF.h:71,73,117,119,128`). `df_dsigma_ij`/`apex_stress` share
one `vv_out` per YF **type parameters only** (e.g.
`DruckerPrager_YF<Alpha,Cohesion>`), so two *different* full
`ASDPlasticMaterial3D<E,Y,P,tag>` combos — different `EL`, different
material tag, even unrelated elements — that happen to reuse the same YF
type share the same buffer. This is H1's defect one layer deeper and with a
**wider** sharing key than `Stiffness` (which R1-A scoped to "every GP,
element and tag of one combo" — this is every GP/element/tag/combo that
picks that YF or PF). BE's loop at 2251/2254 binds `n`/`m` as
`const VoigtVector&` into these statics; within one serial BE call the
values are consumed before the next `df_dsigma_ij` call (checked against
`DruckerPrager_YF.h`'s own `hardening()`, which does not re-enter
`df_dsigma_ij`, so no live corruption today), but this is a threading
blocker strictly worse than the one R1-A recorded for ADR-75b, since it
crosses combo/tag boundaries that `Stiffness` does not.

**F3 (new) — `sendSelf`/`recvSelf` claim SUCCESS (return 0) while doing
nothing** (1223–1237): both print `"...not implemented!!!"` to `cerr` and
`return 0`. Contrast with `revertToStart()` (774–778), which at least
returns `-1`. Any parallel (`OpenSeesMP`/`SP` domain migration) or database
(`save`/`restore`) path involving this material silently loses state and
reports success — no caller has any signal to check, unlike H4's `-1` which
is merely *ignored*.

**F4 (H4 root cause, one level earlier than R1-B cited)** —
`Domain::revertToStart()` itself (`SRC/domain/domain/Domain.cpp:2357-2361`)
calls `elePtr->revertToStart();` in a `while` loop and never inspects the
return value, before `OPS_resetModel` (R1-B's cite) even gets a chance. The
swallow is at the framework's own domain-sweep, not just the Tcl/Python
command layer — makes H4 unfixable from the material side alone.

**F5 — `setParameter`/`updateParameter` id map mixes Commit/Trial
inconsistently.** id 1 (`stress`) writes both `CommitStress` **and**
`TrialStress`; ids 7/8 (`K02D`/`K03D`) and 16–21
(`commitStressIncrementXX..XZ`) write only `CommitStress`, leaving
`TrialStress` stale until the next `setTrialStrain`. All of these set
`stress_set_externally = true`, which (line 232) suppresses the
`first_step` `InitialP0` seed in `setTrialStrain` — so a `K02D`/commit-only
update *can* bypass geostatic `InitialP0` initialization while leaving
`TrialStress` unset, a state the integrators (which start from
`CommitStress`) may not tolerate consistently across id choices. Separately,
`setParameter` (832–844) has its material-tag match commented out to
`if (true)`, so the "does this argv target my tag" guard is dead — currently
masked because the framework only calls `setParameter` on already-tag-
matched materials, but it means the object's own contract does not enforce
what its code implies it should.

**F6 (confirms/extends H14)** — `getCopy()`/`getCopy(type)` (780-829) never
copy `first_step`; `ASDP_TAG` is `#define ASDP_TAG this->getTag()`, so
clones constructed with `(ASDP_TAG)` share the *user's own material tag* —
by design the per-tag `INT_OPT_*`/`GLOBAL_*` maps (4100-4111, static per
full `<E,Y,P,tag>` combo) are correctly keyed and don't collide across
different tags. The factory (`OPS_AllASDPlasticMaterial3Ds.cpp:227-262`)
constructs one full throwaway instance **per registered specialization**
(46×) on every `nDMaterial ASDPlasticMaterial3D <tag> ...` call just to
probe `getYFName()`/`getPFName()`/`getELName()` strings and discards all but
the match — no prototype/`getCopy()` step at parse time; the object handed
back to the parser is the real, tag-bearing instance, and every element's
own `getCopy()` (H14) clones **that**. No `list` mode leak confirmed — `list`
only prints, `instance` stays `nullptr` when `yf_type=="list"`.

## What R1 missed

R1-A/R1-B treated `Stiffness` (H1) and `revertToStart`/`revertToLastCommit`
(H4) as the state-machine defects. Neither lane looked at (a) getters with
side effects (F1), (b) the YF/PF functor static buffers, which are the same
defect pattern with a *wider* sharing scope than `Stiffness` (F2), (c)
`sendSelf`/`recvSelf`'s false-success contract (F3), or (d) that the H4
swallow originates in `Domain::revertToStart()` itself, not just the Tcl/Py
command layer (F4).

## Severity ranking (serial, multi-element model)

1. **H1 + F1 (tied, blocker)** — wrong assembled tangent / clobbered shared
   tangent. Changes iteration count and, if it triggers non-convergence
   cutback, load history — potentially changes **results**, not just
   diagnostics.
2. **H13 (blocker, confirmed by R1-B)** — silent misconfiguration; changes
   **results** with no signal at all.
3. **F3 (major)** — false success on `sendSelf`/`recvSelf`; changes results
   silently, but only exercised under MPI/database, not the default serial
   single-process path this ADR otherwise covers.
4. **H4/F4, H7, H8, H9, H5 (major, per R1)** — unchanged from R1's calls.
5. **F5 (major/doc-gap)** — parameter-ordering trap; changes results only
   for staged/geostatic decks using `K02D`/`K03D`/commit-increment ids.
6. **F2 (major for future threading, doc-only today)** — no live corruption
   found on serial execution; strictly enlarges the ADR-75b blast radius
   R1-A already flagged.
7. **F6/H14, H2 (minor/latent)** — unchanged from R1.

Only H1/F1, H13, and F3 change **results**; F2 and F6 are latent/threading;
the rest are diagnostics- or configuration-shaped as R1 already classified.
