# WP-232f — mode A interface (232f-C, round 1: the freeze candidate)

Status: **proposed** by 232f-C, 2026-11-14; binding once the orchestrator says "frozen". Authority:
ADR `Ladruno_implementation/232_ladruno_mpm_general_driver_adr.md` (§4, §5.1, §5.3, §5.4, §6, §8,
§9). Anchors are quoted from `origin/ladruno` @ `f0b9a1df6`; line numbers drift. 232f-V's
reference scripts were not read. Decisions are numbered (D*n.m*) so a review can cite them.

## 1. The particle element `LadrunoMPMPoint` (`ELE_TAG_LadrunoMPMPoint` = 33026, provisional)

- **D1.1 What it is.** One element per particle: one integration point, one owned `NDMaterial`
  copy (ownership moves from the host to the element at the first mode-A step). Not
  user-creatable: no parser, no `OpenSeesElementCommands` / `TclElementCommands` row; the host
  creates it. Source in `SRC/element/ladrunoMPMPoint/`, so the element lints, which walk
  `SRC/element` (`ci/check_quirk_patterns.py:648`), see it.
- **D1.2 State.** The kinematic truth stays in the host's `ladruno_mpm::Particles`
  (`LadrunoMPMCore.h:287`: x, v, V, F, R, eps, sigHat, sig, sig0, m, V0, lp), plus a new `a`
  (implicit acceleration). The element holds its index p and writes there at G2P, so
  `mpm particleData`, `mpm output` and the reset snapshot work unchanged. The element itself
  holds: p; a host link through a shared liveness token (an element orphaned by
  `remove loadPattern` fails loudly instead of touching freed memory); the step frame (support
  node tags, S_a and ∇S_a at the step-start x_p, inactive flags, projected v̂_a, â_a); the trial
  cache (F_rel, F₁, R₁, J₁, Δε̂, σ, f, K); a state: ARMED, COMMITTED, RESTORED or ORPHAN.
- **D1.3 Nodes = the support at the step-start position.** uGIMP puts node i in the support iff
  |x_p − x_i| < h + l_p (`LadrunoMPMCore.h:91`, `} else if (a < h + lp) {`; zero entries are
  dropped, `:204`, `if (S == 0.0 && dS == 0.0) continue;`). With 0 < l_p ≤ h/2 that open
  interval is 2(h + l_p) ≤ 3h long: **2 or 3 nodes per axis, so 4–9 nodes in 2D and 8–27 in 3D**
  (the scan is `for (int i = c - 1; i <= c + 2; i++)`, `:200`); the linear hat gives 4 / 8.
  Every grid node has ndf = ndm, so numDOF = ndm × nodes and changes from step to step.
- **D1.4 Re-pointing in place.** `repoint(tags, S, ∇S, inactive, v̂, â)` from the host; the
  element refetches each `Node*` with `theDomain->getNode(tag)` (null = host bug, loud), resizes
  its scratch, drops its initial stiffness (`LadrunoShell::dropKi`) and resets the base Rayleigh
  scratch: `index` is protected (`Element.h:257`, `int index, nodeIndex;`) and keyed to the
  numDOF of its first allocation (`Element.cpp:191-201`), so the element sets `index = -1` and
  calls `Element::setRayleighDampingFactors` again. The handler then rebuilds every `FE_Element`
  from `getExternalNodes()` at the forced `domainChanged()`.
- **D1.5 Kinematics** (`update()`, ARMED only). The step reference is the reset grid (u = 0
  everywhere after P2G). F_rel = I + Σ u_a ⊗ ∇S_a; F₁ = F_rel·F_n; then
  `hypoIncrementFromF(Fn, F1, deps, R1, &J1)` (`LadrunoHypoKernel.h:276`) or
  `hypoIncrementFromFPlaneStrain` (`:340`); ε̂₁ = ε̂_n + Δε̂; `setTrialStrain(ε̂₁)`;
  σ̂ = σ̂₀ + σ_mat; σ = R₁σ̂R₁ᵀ (`pushStress6`, `pushStressPlaneStrain` `:401`) — the slice-1
  algebra. A kernel refusal, or a trial centre x_p + Σ S_a u_a outside the real grid (interface
  A7), makes `update()` return −1. `Domain::update` keeps that code
  (`Domain.cpp:2537`, `ladrunoKeepWorstCode(ok, theEle->update());`), so the step fails
  **before** commit and is cut and retried (§3).
- **D1.6 Force and tangent** — LadrunoBrick's hypo assembly with the particle as its single
  point (`LadrunoBrick.cpp:2062`, `resid(3*a + i) += fi * dv;`; `:2078`,
  `pushTangent6(R1, Dmat, Dcur)`; `:2116`, the geometric term). g_a = ∇S_a·F_rel⁻¹,
  V = J₁V₀ (as mode F, `LadrunoMPMCore.h:672`); f_a = V σ g_a − S_a m_p λ b (gravity as a body
  force, W5); K_ab = V B_aᵀ c B_b + V (g_a·σ·g_b) I, c the pushed tangent. Plane strain needs a
  new `pushTangentPlaneStrain` in `LadrunoHypoKernel.h` (the {xx, yy, xy} block of `bond6` for an
  in-plane R), g++-checked against the 6×6 embedding. As in the shipped `-geom hypo` lane, K is
  not the exact linearisation of the Hughes–Winget increment. `getInitialStiff` is the same
  assembly at u = 0 with the material's initial tangent, rebuilt at every re-point. A
  non-associated material gives an unsymmetric K: the user picks a general solver.
- **D1.7 Mass, one model:** lumped, M_aa = S_a m_p I, the row sum of the consistent N m_p Nᵀ
  (Σ S_a = 1). `getMass` returns it; `getResistingForceIncInertia` adds it with
  `LadrunoShell::addNodalInertia(…, diagOnly = true, …)`, snapshotting the static residual before
  the Rayleigh forces (lint).
- **D1.8 Rayleigh:** `alphaM` and `betaK` admitted; `betaK0` and `betaKc` refused at
  `mpm analyze` (they name a stiffness of an earlier support). `commitState` calls
  `Element::commitState()` first (lint; `Kc` stays null). The element overrides no threading
  opt-in (`Element.h:143`, `:155`, default `false`): off every WP-107 / WP-192 path. Ground-motion
  patterns are refused at `mpm analyze` in phase 1.
- **D1.9 `commitState()` — G2P, never failing (W11).** `Element::commitState()`, material
  commit, then a host counter. **G2P is batched:** the commit of the N-th particle element in
  this `Domain::commit` runs G2P for all particles, unless `ladrunoPendingCommitRefusals() > 0`
  (`LadrunoMaterialStatus.h:162`, a read that does not clear) — then nothing moves and
  `Domain::commit` aborts (`Domain.cpp:2343`). Particles therefore stay unmoved on a refusal, as
  in mode F, and G2P still precedes the recorders (`Domain.cpp:2420`). G2P: x += Σ S_a u_a; F, R,
  V, ε̂, σ̂, σ take their trial values; v_p ← α_F(v_p + Σ S_a(v_a − v̂_a)) + (1 − α_F)Σ S_a v_a,
  a_p likewise with a (α_F = `mpm flip`); quasi-static leaves v_p, a_p alone. Then COMMITTED:
  `update()` returns −1 loudly ("a step not prepared by `mpm analyze`"), so a plain `analyze`
  cannot run the particles; `getResistingForce*` return the converged cache (reactions).

## 2. Grid nodes as Domain nodes

- **D2.1 Lifetime.** The Domain holds exactly the active support nodes of the current step, so an
  unsupported lattice node never gets a `DOF_Group` (W3, `PlainHandler.cpp:120-132`). A node
  that stays active keeps its `Node` and tag; one that leaves is retired (its host SPs removed
  first — removal cleans nothing, W5 — then `removeNode` and `delete`); one that enters is new.
  Coordinates x₀ + i h, ghost layer included (interface A5).
- **D2.2 Tags.** Nodes: `nodeBase` + a counter that only increases — a retired tag is never
  reused in the Domain's life (W5); int overflow is a loud refusal. Elements: `eleBase` + p
  (phase 1 never removes a particle). Both bases default to 1 + the largest tag in use at the
  first mode-A step; `mpm tagBase <node> <ele>` overrides. A later collision fails loudly at add.
- **D2.3 Empty support — phase-1 pick: exclude.** Active iff m_a > c_m·max m_p (mode F's rule,
  `LadrunoMPMCore.h:544`; `mpm massCutoff`, default 1e-10). An inactive support node is neither
  a Domain node nor in the particle's connectivity: it is held at zero increment, which is mode
  F's inactive-node semantics, with no renormalisation of S. A particle left with no active node
  is refused (W15: a zero-DOF element exits the process, `FE_Element.cpp:126-129`).
  **Switch for M0b:** `mpm emptyNodes exclude | ghost <γ>`. Ghost keeps every support node and
  adds a spring γ k_p I to K_aa (and γ k_p u_a to f_a) at each inactive node, with
  k_p = M_p⁰ V_p / h² (M_p⁰: the initial normal modulus, slice 1's `normalModulus`). **What it
  changes:** exclude is the γ → ∞ limit (the dropped nodes do not move); ghost lets them move
  against a spring, changing the answer by O(γ) and only through nodes below the cutoff; γ → 0
  is the unregularised, possibly singular, system.

## 3. One mode-A step (the driver lives in the host)

1. Refusals (§4, once per `mpm analyze`); an inadmissible centre is refused (A2).
2. Supports at the committed x_p; active set; node diff (retire / create); re-point every element
   (ARMED).
3. P2G onto Domain nodes: m_a; v̂_a = Σ S m_p v_p / m_a; â_a likewise (quasi-static: 0); then
   `setTrialDisp(0)`, `setTrialVel(v̂)`, `setTrialAccel(â)`, `commitState()` — the committed disp
   is zero (`Node.cpp:1108`).
4. Grid conditions rebuilt as host-pattern SPs (§5).
5. `theDomain->domainChange()` from C++, **unconditionally**, never through the user command
   (W4): `analyzeStep` re-seeds only when the stamp moved (`DirectIntegrationAnalysis.cpp:206`,
   `int stamp = the_Domain->hasDomainChanged();`), and that is where Newmark re-reads the
   committed disp, vel and accel (`Newmark.cpp:398`, `getCommittedDisp()`; W1).
6. t₀ = Domain time. Transient: `dia->analyzeStep(dt)` — never `analyze`, which sub-steps
   (W10). Quasi-static: `sa->analyze(1, false)` (W16).
7. Return ≥ 0: accepted. G2P must have run once for all N (else a loud abort); reactions (§5);
   the particle writer, after `Domain::commit`, as in mode F.
8. Return < 0 and no G2P in this attempt: `analyzeStep` has already reverted the Domain and the
   integrator (`DirectIntegrationAnalysis.cpp:229-236`). Cut dt by `<factor>` (quasi-static:
   `LoadControl::setDeltaLambda`, `LoadControl.h:61`; no other static integrator is admitted),
   re-seed from the unchanged particles and retry; the sub-steps fill the original interval.
   `mpm cut <factor> <maxCuts>` (default 2, 4); exhausted → `mpm analyze` fails with the particles
   at the last accepted step. Nothing to undo.
9. Return < 0 after a commit began (`Domain::commit` returned −1 on a refusal → −4): the run is
   aborted (W17); `mpm reset` clears it.
10. After the last step: `Domain::flushRecorders()`.

## 4. Admission; the two bridges

| Mode (`mpm mode …`) | Needs | Admitted | Refused, loudly |
|---|---|---|---|
| `implicit` | a transient analysis | `getClassTag()` ∈ {Newmark 1, HHT 2, GeneralizedAlpha 44}, as `LadrunoStaggeredDriver.cpp:160` | LadrunoHHT 33013, LadrunoGeneralizedAlpha 33014, every other tag; `VariableTransient` by `dynamic_cast<VariableTimeStepDirectIntegrationAnalysis*>` (both engines store it as the transient analysis: `commands.cpp:3386`, `OpenSeesCommands.cpp:884`); `PFEMAnalysis` likewise |
| `quasiStatic` | a static analysis | `LoadControl` (6) | DisplacementControl (7) — on a grid node by ADR §5.3, and phase 1 has no free FE node to drive (free FE DOFs stay refused); ArcLength (8), ArcLength1 (13), LadrunoArcLength (33004), every other static integrator |
| both | — | Plain handler for homogeneous conditions | a prescribed-velocity condition under `PlainHandler`, which warns and assumes homogeneous (`PlainHandler.cpp:92`); `betaK0` / `betaKc`; ground-motion patterns; `mpm wall`, `mpm damping local` (mode F only) |

`mpm mode explicit | implicit | quasiStatic` (default `explicit`: slice 1 unchanged) is fixed at
the first step of an epoch; changing it needs `mpm reset`. Mode F's Δt refusal does not apply to
mode A. The free-FE-DOF check skips the host's own grid nodes.

**Bridges.** `struct LadrunoMPMAnalysisBridge { StaticAnalysis*; DirectIntegrationAnalysis*;
StaticIntegrator*; TransientIntegrator*; ConstraintHandler*; }` and a shared
`int OPS_LadrunoMPMWith(const LadrunoMPMAnalysisBridge *)`. Classic Tcl: the existing `ladrunoMPM`
(`commands.cpp:765`) fills it from `theHandler`, `theStaticAnalysis`, `theStaticIntegrator`,
`theTransientIntegrator`, `theTransientAnalysis` (`:651`, `:655-656`, `:687-688`). Python: a new
`OPS_LadrunoMPMDL()` in `OpenSeesCommands.cpp` fills it from `cmds->getHandler()`, `getStaticAnalysis()`,
… (`OpenSeesCommands.h:95-131`), the `OPS_LadrunoStaggeredAnalyze` pattern (`:4085`). The
`LadrunoCommandTable.h:113` row becomes `LADRUNO_COMMAND("mpm", OPS_LadrunoMPMDL, ladrunoMPM)`.

## 5. Commands added to `mpm`

| Command | Meaning |
|---|---|
| `mpm mode explicit\|implicit\|quasiStatic` | §4 |
| `mpm analyze <n> <dt>` (implicit) / `mpm analyze <n>` (quasiStatic) | §3; quasi-static steps are the user's `LoadControl` increments |
| `mpm velocityBC <id> box <lo..> <hi..> -v <v..> [-comps <xyz>] [-moving]` | every active node in the closed box gets Δu_d = v_d (t − t₀) on the listed components: a host-pattern `SP_Constraint(node, d, v_d, false)` whose value the host's `applyLoad` sets with `applyConstraint(t − t₀)` (`valueC = loadFactor*valueR;`, `SP_Constraint.cpp:335`). `-moving`: the box translates with v. The rough rigid footing is `-comps xy` with v = (0, −v_f) |
| `mpm fix <face> <comps>` (mode A) | homogeneous host SPs on the face nodes and their ghost layer — mode F's face set (`LadrunoMPMCore.h:794`) |
| `mpm reaction <id>` | Σ `getReaction()` over the condition's nodes at the last accepted step, after `Domain::calculateNodalReactions(flag)` (`Domain.cpp:4589`; flag 0 quasi-static, 1 transient). ndm doubles (classic Tcl: full precision, the slice-1 double-query path) |
| `mpm emptyNodes exclude\|ghost <γ>`, `mpm cut <f> <k>`, `mpm tagBase <n> <e>` | §2, §3 |
| `mpm info acceptedSteps\|g2p\|retries\|domainChanges\|connectivityChanges\|gridNodes\|mode` | g2p = particle G2P updates (= N × accepted steps); connectivityChanges = steps in which any support changed |
| `mpm debug failAttempt <k>` | the first attempt of the k-th next step fails before commit (every element's `update()` returns −1 once) |
| `mpm debug skipReseed 0\|1` | when 1, `domainChange()` is skipped on a step whose node set and every connectivity are unchanged (with GIMP that is stricter than "no cell change": the support edge is at h + l_p). P2G still writes the nodes, so the integrator keeps last step's U, U̇, Ü — the W1 failure, made visible |

## 6. Restart (database)

- The slice-1 refusal (`LadrunoMPMPattern.cpp:185-197`) stays in mode F and is lifted in mode A.
- **Host `sendSelf`:** grid, parameters, faces, velocity conditions, mode and options,
  `Particles` (with `a`) and the reset snapshot, epoch, the abort flag, `nodeBase`, the node
  counter, `eleBase`, the lattice → tag map, its `TimeSeries` (class tag + db tag, as
  `LoadPattern::sendSelf` does). Not sent: the writer (re-issue `mpm output`), the per-step SPs
  (rebuilt each step).
- **Element `sendSelf`:** tag, p, host tag, ndm, connectivity, `alphaM`, `betaK`, material class
  tag, db tag and its own `sendSelf`. Save only between `mpm analyze` calls (COMMITTED).
- **Restore:** the stamp moves every step, so `Domain::recvSelf` takes its clear-and-rebuild
  branch (`Domain.cpp:4094`) through the broker. A RESTORED element answers the `update()` that
  `Domain::addElement` makes (`:520`) with 0, then behaves as COMMITTED; the host re-links its
  elements (tag eleBase + p, class-tag checked) at the next `mpm analyze`. `recvSelf` into a live
  element drops every cache.
- **Broker:** `getNewElement` (`FEM_ObjectBrokerAllClasses.cpp:935`) ← `ELE_TAG_LadrunoMPMPoint`;
  `getNewLoadPattern` (`:2771`) ← `PATTERN_TAG_LadrunoMPM` (the overlay case `:2793` is the model).

## 7. Mutation

The existing family `MPM` (`Ladruno_mutation.h:140-147`; `CMakeLists.txt:244-251`), no new
family. **ZERO** drops the element's V σ g_a term (body force and inertia stay) — the mode-A
counterpart of `LadrunoMPMCore.h:532-535`. **SCALE** multiplies the batched G2P's position
increment by 1.5 — the counterpart of `:581`. The mode-A tests join the MPM suite in
`Ladruno_scripts/mutation_gate.py`. The ordering hazards stay guarded by tests (re-seed,
synchrony).

## 8. Files

**New (fork):** `SRC/element/ladrunoMPMPoint/{LadrunoMPMPoint.h, LadrunoMPMPoint.cpp,
CMakeLists.txt}`; `SRC/domain/pattern/ladrunoMPM/{LadrunoMPMImplicit.cpp, LadrunoMPMBridge.h}`;
`tests/test_ladruno_mpm_implicit_{engine,refusals,restart}.py` (zone_a, classic Tcl and Python);
a g++ check for `pushTangentPlaneStrain`; ledger fragments `WP-232-mpm-implicit-mode-a.md`
(implementations, vanilla, quirks as found). **Edited (fork):** `LadrunoMPMPattern.{h,cpp}`,
`OPS_LadrunoMPM.cpp`, `LadrunoMPMCore.h` (`Particles::a`), `LadrunoMPMIO.cpp` (`MODE` implicit /
quasi-static with `KIND` transient / static, already pinned by 232a §2.3), `LadrunoHypoKernel.h`,
`LadrunoCommandTable.h`, `ladruno_mpm_harness.py`, `mutation_gate.py`, the guide, manifest row
(`feature_class: element`), banner line, `stamp_headers.py` GLOBS. **Upstream — each a
`// Ladruno` marker and a vanilla row:** `SRC/classTags.h`; `FEM_ObjectBrokerAllClasses.cpp`;
`SRC/element/CMakeLists.txt` (`add_subdirectory`); `SRC/interpreter/OpenSeesCommands.{cpp,h}`;
`SRC/tcl/commands.cpp`.

## 9. Gate rows → observables

| Row (wp-briefs "232f") | Observable |
|---|---|
| cantilever QS / transient | `mpm particleData x` of the cases file's tip particle(s), minus the initial x; transient: after each `mpm analyze 1 dt` → first peak, period |
| footing | `mpm reaction <id>` per step → q = −R_y / B (thickness 1) |
| self-weight column, QS | `mpm particleData sig` against ρ₀g(l₀ − Y), Y from the initial `x` |
| plate with a hole | `mpm particleData sig` at the cases file's hole-edge particles; σ_net from `mpm reaction` or the applied traction |
| 1D bar, Newmark | `mpm particleData x` / `v` history → frequency |
| failed-step retry | `mpm debug failAttempt k`; `mpm info g2p` = N × `acceptedSteps`; `retries`; `particleData x v F` against a run given the cut schedule |
| no-cell-crossing re-seed | `mpm info connectivityChanges` proves the step; normal vs `mpm debug skipReseed 1`; `particleData` |
| recorder synchrony | `mpm output -every 1` file against `particleData x sig` read after each step |
| restart | `database File`, `save n`, `wipe` + `restore n`, rebuild the analysis, continue; `particleData` rel 1e-10, both engines |
| refusals | each case of §4: `mpm analyze` raises (Tcl error / Python exception) with the named rule |

## 10. Proposed ADR amendments (none applied)

1. §5.3 / §5.4: G2P runs in `commitState()` **batched at the last particle commit**, so a refusal
   (W17) leaves particles unmoved in mode A too — the per-element form would move siblings first.
2. §5.3: mode-A grid conditions are `SP_Constraint`s owned by the host pattern; a non-homogeneous
   one needs a non-Plain handler.
3. §5.3: quasi-static phase 1 admits `LoadControl` only (no FE seam, so no admissible
   displacement-control node); retry through `setDeltaLambda`.
4. §6 / element guide: `betaK0` / `betaKc` refused for particle elements; `alphaM`, `betaK` kept.
5. W4: node churn changes `LadrunoRecorder`'s topology fingerprint, so `recorder ladruno` would open
   a stage per step — mode-A output is `mpm output` only (guide).

## 11. Open questions for the orchestrator

1. Mode-A G2P blend: FLIP increments with `mpm flip` (default 1), as written — or a PIC default?
2. Mode-A default `massCutoff`: shared 1e-10 (tail nodes ill-conditioned in quasi-static) or a
   mode-A default such as 1e-6? M0b measures either.
3. `ELE_TAG` 33026 and `SRC/element/ladrunoMPMPoint/` (one upstream CMake line).
4. Restart test by `wipe` + `restore` in one process, or a second process as well?
