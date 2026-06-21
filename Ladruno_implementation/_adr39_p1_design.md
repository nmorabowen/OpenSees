# ADR 39 — P1 design: skeleton + handler injection + zero-force hybrid

> Pre-code design for P1. Goes through the adversarial DESIGN gate BEFORE C++
> (mirrors ADR-30 Gate-A). Grounded in the real OpenSees interfaces read 2026-06-21.
> Parent: `39_ladruno_contact_domain_adr.md`. Loop: `_adr39_loop_state.md`.

## P1 goal

Prove the **hybrid plumbing + handler injection** with a **zero-force** contact
adapter — NO narrow phase yet. Acceptance: the same model runs under
`CentralDifferenceLadruno` (explicit) AND `Newmark` (implicit) with contact
defined-but-inactive and gives the **same answer as no-contact** (tolerance, see
Q-P1-3), and the handler delegates the non-contact constraints correctly.

## Grounded interface facts (verified in source)

- `ConstraintHandler` (abstract): `handle(nodesNumberedLast)` is the **only**
  FE_Element/DOF_Group factory; called inside `domainChanged()` after
  `clearAll()` (`DirectIntegrationAnalysis.cpp:413,423`). `setLinks(Domain&,
  AnalysisModel&, Integrator&)`, `update()`, `applyLoad()`, `doneNumberingDOF()`,
  `clearAll()` are lifecycle hooks.
- **`ConstraintHandler::update()` is NOT called in the transient step loop** —
  the loop is `newStep → solveCurrentStep → commit`
  (`DirectIntegrationAnalysis.cpp:215-258`). ⇒ the per-step contact trigger
  CANNOT be the handler's `update()`.
- `FE_Element` has a **bare constructor** `FE_Element(int tag, int numDOF_Group,
  int ndof)` for adapters not backed by a Domain `Element`. Overridable:
  `getResidual(Integrator*)`, `getTangent(Integrator*)`, `addMtoTang`, `getID`,
  `getDOFtags`. A normal `Element` FE reads node trial disp inside its
  `getResistingForce` — an adapter can do the same.
- `LadrunoProjectionHandler` is the in-fork `ConstraintHandler` precedent
  (HANDLER_TAG 33001).

## Design decision: self-contained adapter (narrow phase in getResidual)

**The contact FE_Element computes its own narrow phase when the assembly asks
for its residual/tangent — exactly like a real Element.** No per-step global
trigger, no vanilla integrator edit.

```
class LadrunoContactFE : public FE_Element {
  // connectivity = conservative static superset of nodes its pairs can couple
  // (P1: one adapter per master segment; nodes = segment nodes ∪ candidate slaves)
  const Vector& getResidual(Integrator*) override {
     resid.Zero();
     for (pair in myPairs) {                 // P1: myPairs empty / returns zero
        project slave onto master seg (read Node::getTrialDisp);
        if (gap < 0) { F = -kn*gap*n (+friction); scatter F into resid over myID; }
     }
     return resid;                           // explicit: this is all that's used
  }
  const Matrix& getTangent(Integrator*) override { ... K_c ... }  // implicit only
  void addMtoTang(double) override {}        // contact pairs carry no mass
  double getExplicitCriticalTimeStep() override { return active? 2*sqrt(m/kn):-1; }
};
```

- **Explicit:** `formEleTangent` is mass-only (verified) ⇒ `getTangent` never
  called; `getResidual` injects `F_c`. Self-contained.
- **Implicit:** `getResidual` + `getTangent` both used; Newton reads current trial
  disp each iteration. Self-contained.
- **Friction state** (path-dependent) commits in the adapter's `commitState`-equiv
  (FE_Element has no commitState; the OWNING `LadrunoContactPair` state is committed
  via a hook the handler calls at commit — Q-P1-4).

## Component sketch (P1 scope = skeleton, zero force)

1. **`LadrunoContactDomain`** (Domain owns `LadrunoContactDomain*`, nullable →
   byte-identical when null). Holds surfaces + (P2.5) bucket grid + the adapter
   set. P1: builds adapters from a **brute-force** candidate pairing.
2. **`LadrunoContactSurface`** — from a node-set (slave) or element-face-set
   (master). P1: stores node tags + segment connectivity; no projection.
3. **`LadrunoContactFE : FE_Element`** — the adapter (above). P1: `getResidual` /
   `getTangent` return **zero**; connectivity declared but force = 0.
4. **`LadrunoContactHandler : ConstraintHandler`** — `handle()`:
   (a) do the PlainHandler work (DOF_Groups for all nodes, FE_Elements for all
   Domain elements) — by **composing an owned base handler** (`thePlain->handle()`)
   or replicating it; then
   (b) ask `theContactDomain` for its adapter set and `theModel->addFE_Element(fe)`
   for each. `clearAll()` deletes the adapters. `commit` hook commits pair state.

## Open questions for the gate

> [!question] Q-P1-1 (per-step trigger)
> With the self-contained-adapter design, is a per-step ContactDomain trigger
> needed at all for P1/P2, or only at epoch re-emit (P2.5)? Confirm the adapter
> reading `Node::getTrialDisp()` inside `getResidual` sees the CURRENT iterate in
> both explicit (`newStep`-set trial) and implicit (Newton-updated trial). Any
> case where the residual is requested with stale trial disp?

> [!question] Q-P1-2 (handler delegation)
> Compose an owned `PlainHandler` and call its `handle()`, then add contact FEs?
> Or replicate PlainHandler's DOF_Group/FE loop? Composition risk: the owned
> handler's `setLinks`, numbering (`doneNumberingDOF`), and `clearAll` must chain
> correctly. Does adding contact FEs AFTER the base `handle()` but before
> `doneNumberingDOF` keep the numbering valid (contact FEs reuse existing node
> DOFs — no NEW equations)?

> [!question] Q-P1-3 (zero-force is NOT byte-identical — connectivity changes the graph)
> A contact adapter declares connectivity (graph edges between surface nodes) even
> at zero force. That changes the `DOF_Graph` → the numberer may produce a
> different equation order / bandwidth → results equal but **NOT byte-identical**.
> So the P1 gate must assert **same answer to ~1e-12**, NOT bitwise. OR: at
> zero-force/inactive, the adapter declares EMPTY connectivity (getID size 0) so
> the graph is untouched — is an empty-connectivity FE_Element legal (does
> addFE_Element / DOF_Graph tolerate it)? Decide which.

> [!question] Q-P1-4 (pair state commit)
> `FE_Element` has no `commitState`. Friction/gap state lives in
> `LadrunoContactPair`. Who commits it each step? Options: the handler's
> `applyLoad`/a commit hook; or the adapter registers a `Recorder`-like callback;
> or ContactDomain commits all pairs when the Domain commits. Which is clean +
> serialization-safe?

> [!question] Q-P1-5 (conservative-static connectivity, P1 form)
> P1 uses brute-force pairing (every slave × every master segment). Connectivity
> per adapter = ? Simplest: ONE adapter over ALL surface nodes (dense but trivially
> correct for P1's small tests). Per-segment adapters come with bucket sort (P2.5).
> Is the one-big-adapter P1 simplification acceptable, or does it mask a real
> numbering/assembly issue that per-segment would expose?

## P1 acceptance gates (tests)

- Zero-force regression: 2-block model, contact defined, runs under CDL **and**
  Newmark, displacements match no-contact to 1e-12 (Q-P1-3 tolerance).
- Handler delegation: a model with an `equalDOF` + the contact handler →
  the equalDOF is still enforced (delegation works).
- `addFE_Element` of the adapter succeeds; `clearAll`/rebuild roundtrip clean
  (no leak, no dangling).
- Rigid-body translation of both bodies → ΣF_c = 0 (trivially true at zero force;
  real test lands in P2).
