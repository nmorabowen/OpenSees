# ADR 39 — P1 design: skeleton + handler injection + zero-force hybrid

> Pre-code design for P1. **Revised after the adversarial design gate** (Workflow
> wv35sge9t, 4 reviewers + synth → SALVAGEABLE-WITH-CHANGES; all fixes folded in,
> flagged `[GATE]`). Grounded in real OpenSees source. Parent:
> `39_ladruno_contact_domain_adr.md`. Loop: `_adr39_loop_state.md`.

## P1 goal

Prove the **hybrid plumbing + handler injection + commit/revert lifecycle** with
a **zero-force** contact adapter — NO narrow phase. Split into **P1a** (prove
hybrid+injection, ≈1 vanilla edit) then **P1b** (Domain surface + lifecycle
hooks). `[GATE]`

## Grounded interface facts (verified in source)

- `ConstraintHandler::handle(nodesNumberedLast)` is the **only** FE_Element/
  DOF_Group factory; called inside `domainChanged()` after `clearAll()`
  (`DirectIntegrationAnalysis.cpp:412,423`). `update()` is NOT in the step loop.
- **`ConstraintHandler::applyLoad()` IS called every step** via
  `AnalysisModel::applyLoadDomain`/`updateDomain` (`AnalysisModel.cpp:567,604`) —
  but PRE-solve on TRIAL (rejectable) state. ⇒ right place for the **broad phase**
  (P2.5), WRONG place for state **commit**. `[GATE MINOR-8]`
- Adapter reading `Node::getTrialDisp()` in `getResidual` sees the CURRENT iterate
  in BOTH families: explicit sets trial via `setResponse`+`updateDomain` in
  `newStep` before `formUnbalance` (`CentralDifferenceLadruno.cpp:509-518`);
  implicit updates trial each Newton iter (`Newmark.cpp:482-505`,
  `NewtonRaphson.cpp:200-215`). `[GATE Q-P1-1 RESOLVED]`
- `FE_Element` bare ctor `(tag, numDOF_Group, ndof)` → `myEle==0`. **With
  `myEle==0` the base `getResidual`/`getTangent` `exit(-1)` and `zeroResidual`/
  `addMtoTang`/etc. early-return** (`FE_Element.cpp:167-194,296,323`). ⇒ the
  adapter MUST override `getResidual` AND `getTangent` AND own its own
  `Vector`/`Matrix` buffers. `[GATE MAJOR-4]`
- **Explicit DOES call `getTangent(this)`** (`IncrementalIntegrator::formTangent`
  `:110-124`, unconditional) — it just delegates to `formEleTangent`, which under
  CDL is mass-only (`CentralDifferenceLadruno.cpp:222-227`, no `addKtToTang`) ⇒
  zero contact LHS. The earlier "explicit never calls getTangent" was imprecise.
  `[GATE MAJOR-4]`
- Implicit assembles the tangent as `addKtToTang(c1)` (`Newmark.cpp:289-294`) ⇒ the
  adapter's `getTangent` must return **`c1·K_c`** (read the integrator combo
  factors), not raw `K_c`. `[GATE MAJOR-3]`

## Design decision: self-contained, STATELESS-VIEW adapter

**The contact FE_Element computes its own narrow phase in `getResidual`/
`getTangent` (like a real Element) — no per-step trigger, no vanilla integrator
edit.** The adapter is a **stateless VIEW**: all path-dependent pair state
(gap0, friction slip/stick) lives on the **Domain-owned `LadrunoContactDomain`**,
re-bound to fresh adapters by pair-key at every `handle()`. `[GATE Q-P1-4]`

```cpp
class LadrunoContactFE : public FE_Element {     // ELE_TAG 33015
  Vector resid;  Matrix tang;                    // OWN buffers (myEle==0)
  // P1a: empty connectivity -> FE_Element(tag,0,0); resid size 0, tang 0x0
  const Vector& getResidual(Integrator*) override {            // both families
     resid.Zero();
     for (key in myPairKeys) { st = contactDomain->pair(key);  // VIEW: read state
        project (Node::getTrialDisp); if (gap<0) scatter F into resid over myID; }
     return resid;                               // P1: zero
  }
  const Matrix& getTangent(Integrator* I) override {           // implicit only
     tang.Zero(); ... K_c ...; tang *= I->getCfactors...c1;    // [GATE MAJOR-3]
     return tang;
  }
  void addMtoTang(double) override {}            // contact pairs carry no mass
};
```

## Lifecycle (the gate's BLOCKER fixes)

- **Commit:** `Domain::commit()` is the single integrator-agnostic choke point
  (`Domain.cpp:2157-2186`). Add `// Ladruno: if (theContactDomain)
  theContactDomain->commit();`. (Domain iterates only Nodes/Elements, never
  AnalysisModel FE_Elements ⇒ the adapter can't self-commit.) `[GATE BLOCKER-1]`
- **Revert:** `Domain::revertToLastCommit()` (`Domain.cpp:2188-2214`) — failed
  implicit steps call it (`DirectIntegrationAnalysis.cpp:201,219,230,260`). Add
  `// Ladruno: theContactDomain->revertToLastCommit();` else rejected-trial
  friction state leaks into the retry. **Was missing.** `[GATE BLOCKER-2]`
- **clearAll/rebuild:** `domainChanged` does `AnalysisModel::clearAll()` (destroys
  adapters) NOT `Domain::clearAll()`, so the Domain-owned `LadrunoContactDomain`
  (and its pair state) SURVIVES. New adapters re-bind to it by pair-key in
  `handle()`. `[GATE Q-P1-4]`
- **Serialization** rides `LadrunoContactDomain::sendSelf/recvSelf` (or `Node`
  slots à la `projTieForce`), never the FE. P1 single-proc stubs `→0`.

## Handler: REPLICATE PlainHandler, do not compose `[GATE Q-P1-2 / MAJOR-5]`

`LadrunoContactHandler::handle()` mirrors `LadrunoProjectionHandler::handle()`
(`:148-218`): replicate the PlainHandler DOF_Group + element-FE loop, **track
`numFe` locally** (PlainHandler returns `count3` ≠ `numFe`, `:299`, so a composer
can't know the next FE tag — and `addFE_Element` silently drops duplicate tags),
then append the contact adapter(s) with tags continuing above `numFe`. Override
`doneNumberingDOF` to chain `this->ConstraintHandler::doneNumberingDOF()`.
Numbering stays valid: contact FEs reuse existing node DOF_Groups (zero new
equations); `setID()` auto-maps their equation IDs.

## P1a — minimal (build first)

Files: `SRC/domain/contact/LadrunoContactFE.{h,cpp}` (empty-connectivity zero
adapter, owns buffers), `SRC/analysis/handler/LadrunoContactHandler.{h,cpp}`
(replicate-PlainHandler + inject one empty adapter, hardcoded — no ContactDomain
yet). Vanilla: the `constraints LadrunoContact` branch in `OPS_ConstraintHandler`
(`OpenSeesCommands.cpp`). classTags: `HANDLER_TAG_LadrunoContactHandler 33002`
(after 33001), `ELE_TAG_LadrunoContactFE 33015` (ELE band stops at 33014).
`sendSelf/recvSelf → 0` stub. CMake: add `SRC/domain/contact/CMakeLists.txt` +
the handler file to the analysis target.

**P1a acceptance (`tests/test_adr39_contact_p1.py`, zone_a):**
- 2-block model runs under `CentralDifferenceLadruno` AND `Newmark` with
  `constraints LadrunoContact`; displacements **BITWISE** identical to no-contact
  (empty connectivity ⇒ graph-neutral). `[GATE Q-P1-3: empty conn + bitwise]`
- An `equalDOF` is still enforced through the handler (delegation works).
- `clearAll`/rebuild roundtrip clean, incl. a forced unrelated `domainChanged`
  (`remove element` mid-run) → contact behaviour unchanged. `[GATE MINOR-11]`

## P1b — Domain surface + lifecycle hooks (on proven plumbing)

`LadrunoContactDomain*` on `Domain` (nullable → byte-identical when null) +
`LadrunoContactSurface`; the `Domain::commit()`/`revertToLastCommit()` `// Ladruno`
hooks (wired, no-op at zero force); `contactSurface`/`contact` parser + Tcl/Py
wrapper registrations. P2 switches the adapter to **declared per-segment
connectivity** + the 1e-12 (not bitwise) gate. `[GATE Q-P1-5: per-segment at P2]`

## Dead wiring removed `[GATE MINOR-10]`

`getExplicitCriticalTimeStep()` on the adapter is UNREACHABLE — `CriticalTimeStep`
scans `Domain->getElements()` only (`CriticalTimeStep.cpp:299,310`), never
AnalysisModel FE adapters. Dropped from the P1 sketch; P3 routes contact `dt_cr`
via `LadrunoContactDomain` (new CFL seam or a participating Domain element).

## Resolved open questions (all closed by the gate)

- Q-P1-1 trigger → self-contained `getResidual` correct; no per-step trigger
  (broad-phase trigger = `applyLoad`, only needed at P2.5).
- Q-P1-2 delegation → REPLICATE, track `numFe` locally.
- Q-P1-3 byte-identical → EMPTY connectivity + **bitwise** gate (P2 → declared + 1e-12).
- Q-P1-4 commit → Domain-owned state + `Domain::commit()`/`revertToLastCommit()`
  hooks; adapters are stateless views.
- Q-P1-5 connectivity → moot in P1 (empty); per-segment at P2.
