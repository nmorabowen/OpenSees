# ADR-30 P1 design (equalDOF core) — pseudocode for Gate-A

Working doc; the input the Gate-A math panel attacks BEFORE any C++. Scope = P1:
`equalDOF` only, the handler + projector + CDL hook. Spec: ADR 30.

## Components
1. `LadrunoProjectionConsumer` (abstract interface, 1 method)
   ```
   class LadrunoProjectionConsumer {
     virtual void setConstraintProjector(LadrunoConstraintProjector*) = 0;
   };
   ```
   `CentralDifferenceLadruno` implements it (stores the pointer; non-owning).

2. `LadrunoConstraintGraph` — builds the undirected DOF-DOF graph from the MPs and
   returns connected components (groups). For equalDOF each MP is one edge
   (retainedDOF — constrainedDOF). Group = set of node/dof pairs + which are retained.

3. `LadrunoConstraintProjector` — per group g:
   - `ID allEqn_g`   equation numbers of every DOF in g, stacked = rows of L
   - `ID retIdx_g`   which rows of allEqn_g are RETAINED (the free master DOFs)
   - `Matrix C_g`    constrained = C_g · retained  (equalDOF ⇒ C_g = I, one retained)
   - `Matrix L_g`    = [ I ; C_g ] (allEqn_g.Size() × retIdx_g.Size())
   - `Matrix LtML_g` factored (dense, tiny) — built lazily from diag(M)
   - `Vector m_g`    the group's per-row lumped mass (diag of M_g)

4. `LadrunoProjectionHandler` (HANDLER_TAG 33001).

## Handler.handle()  (Plain-style assembly, MPs NOT eliminated)
```
for each Node n:
    dof = new DOF_Group(n)
    set all its ID slots to -2                       # free, to be numbered
    for each homogeneous SP on n: set that slot -1   # excluded from equations
    # *** KEY DEVIATION FROM PlainHandler: do NOT set MP-constrained slots to -4. ***
    #     A slave DOF keeps its own equation and its own diagonal mass; the
    #     projection enforces the tie at the acceleration level, not by elimination.
    addDOF_Group(dof)
setNumEqn(count of -2 slots)
create FE_Elements for all elements (as PlainHandler)

# --- build group structure (node/dof, not yet equation numbers) ---
graph = LadrunoConstraintGraph(all MP_Constraints)
groups = graph.connectedComponents()
runDiagnostics(groups)        # §5.2, see below — errors at handle() time
store groups (node/dof lists + C) for doneNumberingDOF()
return count3
```

### runDiagnostics (the T5 battery, equalDOF-relevant subset)
- slave DOF under TWO MPs                         → named error (both tags)
- chain: a slave DOF is retained by another MP    → named error (both tags)
- SP-fixed DOF that is also an MP slave            → named error
- redundant cycle (u1=u2,u2=u3,u3=u1)              → per-group rank(L) < cols ⇒ error
- zero free DOFs after exclusion                   → clean error (not segfault)
- massless retained OR massless slave DOF          → error w/ recipe
                                                     (P0 proved the SOE won't catch it)
- duplicate SP on one DOF                          → warning (Auto parity)

## Handler.doneNumberingDOF()   (equation numbers now known)
```
for each group g:
    allEqn_g  = [ node.dof_group.getID()[localdof] for each (node,dof) in g ]
    drop rows whose eqn < 0 (SP-fixed retained ⇒ delete that column of L, §2.4)
    retIdx_g  = rows that are retained masters
    C_g, L_g  = build from g (equalDOF ⇒ C_g = I)
theProjector = new LadrunoConstraintProjector(groups' eqn data)
# push to the integrator via the interface (NOT a concrete downcast):
consumer = dynamic_cast<LadrunoProjectionConsumer*>( getIntegratorPtr() )
if consumer == 0:
    opserr << "LadrunoProjectionHandler requires an explicit projection-aware "
              "integrator (e.g. CentralDifferenceLadruno); use constraints "
              "Transformation/Penalty for implicit analyses\n";
    return -1                                   # hard, early, named
consumer->setConstraintProjector(theProjector)
return 0
```

## Projector.buildMass(LinearSOE* soe)   (lazy, once per domainChanged)
The integrator calls this right after its first M-only `formTangent()` (mass is
constant between domain changes). Reads the per-equation lumped mass.
```
for each group g, each row r with eqn e = allEqn_g[r]:
    m_g[r] = diagonalMass(e)        # see MASS SOURCE below
    if m_g[r] <= 0 for a row that L couples into a retained dir: error (massless)
    LtML_g = factor( L_gᵀ diag(m_g) L_g )
```
**MASS SOURCE (open, for Gate-A):** candidate routes to diag(M) at equation e —
(a) read DiagonalSOE's stored diagonal A[e] (needs a small accessor/friend; only
valid for a Diagonal SOE — acceptable since the handler is explicit-only); (b)
assemble from Node masses via DOF_Group getID()↔getMass() (SOE-agnostic). Pick in
impl; correctness is identical (both are the lumped diagonal).

## Projector.project(Vector& a)    (in-place, per group)
```
for each group g:
    araw = a[allEqn_g]                       # gather (size n_all)
    rhs  = L_gᵀ ( diag(m_g) · araw )         # size n_ret
    ar   = LtML_g.solve(rhs)                 # (LᵀML)⁻¹ Lᵀ M araw
    afull= L_g · ar                          # size n_all  (retained AND constrained)
    a[allEqn_g] = afull                       # scatter back (overwrites BOTH)
```
equalDOF degenerate check: n_ret=1, L=[1;1], m=(m_r,m_c) ⇒ ar=(m_r araw_r+m_c araw_c)
/(m_r+m_c), af=(ar,ar). = Theory Eq. 28.1. ✔

## CDL hook (CentralDifferenceLadruno.cpp)
New members: `LadrunoConstraintProjector* theProjector=0; bool massBuilt=false;
Vector Aproj;`  + implement the consumer interface.

1. domainChanged(): `massBuilt=false;` (mass cache invalid). IC compliance check:
   for each group, verify `C·u_r − u_c ≈ 0` and `C·v_r − v_c ≈ 0` from committed
   state; warn+abort by default, project if `-projectICs` (O2).
2. newStep() starter block, AFTER `*Aprev = theLinSOE->getX();` (line 446):
   ```
   if (theProjector && !massBuilt) { theProjector->buildMass(theLinSOE); massBuilt=true; }
   if (theProjector) theProjector->project(*Aprev);     # project a0
   ```
   THEN the existing `Vhalf->addVector(1.0,*Aprev,-0.5*deltaT)` seeds v_{-1/2} from
   the PROJECTED a0. (Ordering matters: project before the seed.)
3. update(U) main path, before the Vfull computation (~line 507):
   ```
   Aproj = U;
   if (theProjector) theProjector->project(Aproj);
   # use Aproj (not U) for Vfull, setResponse, and Aprev:
   *Vfull = *Vhalf;  Vfull->addVector(1.0, Aproj, 0.5*deltaT);
   theModel->setResponse(*Ut, *Vfull, Aproj);
   ...
   *Aprev = Aproj;
   ```
   (Circuit-breaker max|a| check runs on U as before, BEFORE projection.)

## Correctness claims for Gate-A to attack
- C1 project() = M-orthogonal projection onto range(L): `a_proj = L(LᵀML)⁻¹LᵀM a_raw`,
  idempotent (P²=P), self-adjoint in the M-inner-product.
- C2 momentum: `Lᵀ M (a_proj − a_raw) = LᵀM a_proj − LᵀM a_raw = LᵀML(LᵀML)⁻¹LᵀM a_raw
  − LᵀM a_raw = 0`. Retained-conjugate momentum conserved exactly, ANY mass.
- C3 manifold invariance: leap-frog is linear; if u0,v0 comply and every a (incl a0)
  is projected, then `C a_r=a_c ∀steps ⇒ C v_r=v_c, C u_r=u_c ∀t` to machine eps
  (constant C). No drift term.
- C4 Δt-neutrality: range(L) is a subspace ⇒ Rayleigh-quotient interlacing ⇒
  ω_max(constrained) ≤ ω_max(unconstrained). The per-element dt_cr stays conservative.
- C5 equalDOF reduces to the mass-weighted average (Theory Eq. 28.1).
- C6 the slave DOF keeping its own equation + diagonal mass is consistent: its
  a_raw is overwritten by project(); its mass enters LᵀML (the m_c term). No double
  count (the element/inertia is assembled once per DOF, Plain-style).

---
# Gate-A resolutions (BINDING — supersede the above where they conflict)

Panel `wf_09671e08` (15 agents, 4 lenses). **Algebra core C1/C2/C5 verified SOUND**
(M-orthogonal projector, momentum, equalDOF→mass-weighted average; operator is
correctly general-C so P2 is not blocked). All blockers were in the seam/hook. Fixes:

**R1 [BLOCKER] IC check must use the MP offset (incremental) convention.** OpenSees
MP enforces `u_c − Uc0 = C(u_r − Ur0)`, NOT `u_c = C u_r`. Precompute per group the
constant offset `δ = Uc0 − C·Ur0` (from `MP_Constraint::getConstrainedDOFsInitial
Displacement()` / `getRetainedDOFsInitialDisplacement()`); check `(u_c − C·u_r) ≈ δ`.
Velocity check carries NO offset: `C·v_r − v_c ≈ 0`. The acceleration projector is
offset-agnostic (offsets differentiate away; `a_c=C a_r` exact) → C1–C6 unchanged.
*Consequence for tests:* this handler matches the Penalty/Lagrange INCREMENTAL
convention; vs a Transformation reference there is a constant `C·Ur0−Uc0` displacement
offset unless constraints are built at a coincident/zero-offset config — build T2/T3
ties at zero offset (or subtract the constant) so a correct run is not misread as a bug.

**R2 [BLOCKER] Mass source = route (a) ONLY: read the assembled Diagonal-SOE diagonal.**
Route (b) (Node::getMass) is WRONG — it omits element-contributed mass (element `-rho`/
`-cMass` mass enters via `FE_Element::addMtoTang`, never `Node::getMass`). Only the SOE
diagonal equals the `M` the integrator inverts (required for C1/C2/C6). `buildMass`:
`dynamic_cast<DiagonalSOE*>(soe)`; null → hard named error ("requires system Diagonal").
**Timing (CRITICAL):** `DiagonalDirectSolver::solve()` overwrites `A[e]` with `1/aii`
in place (verified, DiagonalDirectSolver.cpp:124) and sets `isAfactored`. So `buildMass`
MUST read the diagonal **between `formTangent` (CDL line 434) and `solve` (line 442)** in
the starter, where `A[e] == m`. Add a one-line accessor to DiagonalSOE
(`const double *getDiagonalA() const { return A; } // Ladruno`, LEDGER_vanilla row).

**R3 [MAJOR] update() hook — `*Aprev` (line 533) is the load-bearing write.** Project
`U` once at the TOP of `update()`, after the NaN/Inf breaker (496–501, runs on raw U) and
before line 507, into the `Aproj` member. Consume `Aproj` at ALL THREE sites: `Vfull`
(507–508), `setResponse` (511), AND **`*Aprev = Aproj;` at line 533** (was `*Aprev = U;`).
The main leap-frog advance `Vhalf += dt*Aprev; Ut += dt*Vhalf` lives in newStep() (452–453)
and consumes `*Aprev` — so manifold invariance (C3) rides entirely on line 533 being the
projected accel. KE/Vfull breakers unchanged (Vfull already projected).

**R4 [MAJOR] SP-fixed retained DOF → zero-overwrite, not column-deletion-to-empty.** If
every retained column feeding a slave row is SP-fixed (deleted), the slave is NOT free —
it must take `a_c = C·0 = 0`. Carry a per-group "fixed-row" set; `project()` first writes
`a[eqn_c] = 0` for orphaned slave rows, then runs the M-orthogonal solve over the remaining
genuinely-free sub-group (may be empty for pure equalDOF — fine). Prune `allEqn_g`/`m_g` in
lockstep with dropped L rows so gather never indexes a negative eqn.

**R5 [MAJOR] Consistent-mass guard.** The projector weight `diag(M)[e]` is only the true
inertia under LUMPED mass. At buildMass, for each tied (retained/constrained) DOF scan its
incident FE_Elements' `getMass()`; if any has a non-negligible off-diagonal touching a tied
equation (`max|M(p,q)|, p≠q > tol·max|M(p,p)|`), **refuse** with a named error ("requires
lumped mass; element <tag> assembles consistent mass on a tied DOF; use -cMass 0 or
constraints Transformation"). Mirrors the RBE2 massless-scan precedent + the P0 SOE-trust
finding. (Scope-limited to tied DOFs → bounded cost.)

**R6 [MAJOR] IC check: cover velocity, every domainChanged, null-projector safe.** Run the
R1 check on BOTH u and v from committed state, every `domainChanged()`, default warn+ABORT
(return −1, named group/node/dof). If `theProjector==0` at check time (the
`DirectIntegrationAnalysis::setIntegrator`-after-stamp path bypasses `doneNumberingDOF`),
emit a WARNING (never silent skip) + LEDGER_quirks row. `-projectICs` projects BOTH u0 and
v0 on the INCREMENT before the line-286–290 seed so `v_{-1/2}` is born on-manifold.

**R7 [MINOR→fold] Graph = one edge per (retained,constrained) DOF PAIR (row of the MP),
not one per MP** (`equalDOF R C 1 2 3` = one MP, three pairs). Retain the retained/
constrained ROLE per DOF so the chain + SP∩slave diagnostics are decidable; build L from
the FULL edge set (no spanning-tree) so the redundant-cycle `rank(L)<cols` check can fire.

**R8 [MINOR] project() entry invariant.** Enforce `massBuilt ⇒ buildMass ran` at the top of
`project()` (error if not), don't infer it from firstStep co-occurrence.

**R9 [MINOR] Diagnostic phase-split.** handle()-time (node/dof only): slave-under-2-MPs,
chain, SP∩slave, redundant-cycle rank, zero-free-DOF. buildMass-time (needs assembled
diagonal): massless (rank-deficient `LᵀML`), consistent-mass guard. Demote the old
handle()-time massless Node-mass scan to a non-fatal hint.

**R10 [doc] Sharpen C4.** `a_proj = L(LᵀML)⁻¹Lᵀr` (the M⁻¹ from the solve cancels), so the
projected leap-frog IS the central difference of the Galerkin-reduced pencil `(LᵀKL, LᵀML)`,
whose spectrum interlaces the full pencil's ⇒ `ω_max` can only drop. The M-weight is
*required* for this cancellation.
