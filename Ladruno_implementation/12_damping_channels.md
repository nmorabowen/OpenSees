---
title: Damping channels in OpenSees — complete reference
project: Ladruno
status: reference
priority: medium
owner: nmora
tags:
  - reference
  - solver
  - damping
  - explicit
---

# Damping channels in OpenSees — complete reference

> [!summary]
> There are **six** distinct mechanisms by which damping enters the OpenSees
> equations of motion. They live on different objects, combine by **different
> rules** (some overwrite, some add), and behave differently under implicit vs
> explicit integration. This doc is grounded in the source (line refs below) and
> in a runnable verification battery, `tests/test_damping_channels.py`
> (8 passed, 1 xfail). Two behaviours here were **only** settled by testing —
> they are not visible from reading the code.

Everything feeds the same equation:

$$M\ddot u + \underbrace{C}_{\text{damping}}\dot u + R(u) = P(t)$$

`C` is **never** stored as one global matrix. It is assembled on the fly inside
the integrator's `formTangent` / `formUnbalance` from per-object contributions.
That assembly is the key to understanding double-assignment.

---

## The six channels at a glance

| # | Channel | Set by | Stored on | C contribution | Combine rule |
|---|---------|--------|-----------|----------------|--------------|
| 1 | **Element Rayleigh** | `rayleigh`, `region -rayleigh`, `setElementRayleighFactors` | each `Element` (`alphaM,betaK,betaK0,betaKc`) | αM·Mₑ + βK·Kₜ + βK0·K₀ + βKc·K_c | **OVERWRITE** per element |
| 2 | **Nodal Rayleigh** | `rayleigh ... -node`, domain/region rayleigh | each `Node` (`alphaM` only) | αM·M_node | **OVERWRITE** per node |
| 3 | **Modal damping** | `modalDamping`, `modalDampingQ` | `Domain` → integrator | 2ζᵢωᵢ via M-orthonormal φᵢ | **ADDITIVE** (on top of 1,2) |
| 4 | **Element `Damping` objects** | `damping` cmd + `-damp` on element/region | each supporting `Element` | frequency-band viscous, in residual | **ADDITIVE** (overwrites only other `-damp`) |
| 5 | **Material / dashpot** | `Viscous`, `ViscousDamper`, Maxwell… | inside constitutive law | enters `R(u̇)`, **not** C | additive (it's real force) |
| 6 | **Numerical / algorithmic** | HHT-α, generalized-α, bulk viscosity | integrator | not C — α-dissipation | n/a |

Mental model: **channels 1 & 2 are "Rayleigh" — per-object scalars that get
overwritten.** Channels 3 & 4 are *separate objects* that *add* to whatever
Rayleigh you already set. That asymmetry is the #1 source of "I accidentally
double-damped / under-damped" confusion.

---

## 1. Element Rayleigh damping

[`Element::setRayleighDampingFactors`](../SRC/element/Element.cpp) (Element.cpp:116)
stores four scalars — **plain assignment**, so it **overwrites**, never accumulates.
[`Element::getDamp()`](../SRC/element/Element.cpp) (Element.cpp:211) builds:

```cpp
if (alphaM != 0.0) M->addMatrix(0.0, this->getMass(),         alphaM);
if (betaK  != 0.0) M->addMatrix(1.0, this->getTangentStiff(), betaK);   // current tangent
if (betaK0 != 0.0) M->addMatrix(1.0, this->getInitialStiff(), betaK0);  // initial
if (betaKc != 0.0) M->addMatrix(1.0, *Kc,                     betaKc);  // last-committed
```

So `Cₑ = αM·Mₑ + βK·Kₜ(t) + βK0·K₀ + βKc·K_c`. The four β's pick *which* stiffness.
**For nonlinear runs prefer `betaK0` (initial)** — tangent `betaK` can go zero/negative
on softening and destabilize.

Three entry points, all routing to the same per-element setter:
- `rayleigh αM βK βK0 βKc` → [`OPS_rayleighDamping`](../SRC/interpreter/OpenSeesMiscCommands.cpp) (line 119). With no flag it calls [`Domain::setRayleighDampingFactors`](../SRC/domain/domain/Domain.cpp) (Domain.cpp:2119) which **loops every element and node**.
- `rayleigh αM βK βK0 βKc -ele $tag…` / `-node $tag…` — targeted.
- `region $r … -rayleigh αM βK βK0` → [`MeshRegion::setRayleighDampingFactors`](../SRC/domain/region/MeshRegion.cpp) (MeshRegion.cpp:289).

> [!warning] Overwrite, not additive — VERIFIED
> ```tcl
> rayleigh 0.0 0.0 0.002 0.0              ;# domain: every element gets βK0=0.002
> region 1 -ele 5 -rayleigh 0.0 0.0 0.006 0.0  ;# element 5 OVERWRITTEN to 0.006
> ```
> Element 5 ends at βK0=0.006 — the domain value is **gone**, not summed.
> `test_element_rayleigh_overwrites_domain` measures ζ≈0.06 (not 0.08).

---

## 2. Nodal Rayleigh damping (mass-proportional only)

[`Node::setRayleighDampingFactor`](../SRC/domain/node/Node.cpp) (Node.cpp:1221) stores a
**single** `alphaM`; a node has no stiffness, so there is no β term.
[`Node::getDamp()`](../SRC/domain/node/Node.cpp) (Node.cpp:1227) returns `αM·M_node`.

A domain/region `rayleigh` pushes `alphaM` onto **both** nodes and elements. Whichever
object carries the mass (nodal lumped vs element consistent) is the one whose
mass-proportional damping bites — generally no double-count, *unless you hand-set both
nodal and element mass for the same DOF*.

---

## 3. Modal damping — scope & semantics

Set by [`modalDamping`](../SRC/interpreter/OpenSeesCommands.cpp) (line 3215) /
[`modalDampingQ`](../SRC/interpreter/OpenSeesCommands.cpp) (line 3276), stored on the
`Domain`, consumed by [`IncrementalIntegrator`](../SRC/analysis/integrator/IncrementalIntegrator.cpp).

> [!important] Hard prerequisite: `eigen N` must run first.
> The factor vector is sized to the number of computed eigenvalues; with none, you get
> a warning and **zero** damping.

**What it computes** (IncrementalIntegrator.cpp:529–554, using M-orthonormalized φᵢ from
`setupModal` :483–494):

$$C = \sum_i 2\,\zeta_i\,\omega_i\,(M\phi_i)(M\phi_i)^T$$

This is a **true modal damping matrix**, not Rayleigh — it assigns exactly ζᵢ to mode i,
no spillover, no frequency-dependence artifacts.

### Scope
- **Only the modes you computed are damped.** `eigen 3` + `modalDamping 0.05` → modes 4+
  have **zero** modal damping (FE high-frequency content is undamped by this channel; often
  paired with a little stiffness Rayleigh to mop up).
- **Global only** — no per-region/per-element modal damping.
- **ADDITIVE with Rayleigh** — `rayleigh` + `modalDamping` both apply, summed
  (`test_rayleigh_plus_modal_is_additive`: 0.03+0.04→0.07). Easy to over-damp.

### `modalDamping` vs `modalDampingQ` — the difference is one flag
`setModalDampingFactors(&v, inclMatrix)` — `true` vs `false`:
- **`modalDamping`** → `inclModalMatrix=true` → damping **matrix added to the tangent**
  ([`addModalDampingMatrix`](../SRC/analysis/integrator/IncrementalIntegrator.cpp) :563,
  called in `formTangent` TransientIntegrator.cpp:88) **and** force added to RHS. Consistent implicit.
- **`modalDampingQ`** → `inclModalMatrix=false` → **only the force** is added to the
  unbalance (TransientIntegrator.cpp:135). The tangent is *not* modified.

> [!bug] modalDampingQ applies damping with the WRONG SIGN — confirmed structural
> In free vibration, `modalDampingQ ζ` *amplifies* the response: measured ζ ≈ **−ζ_target**
> under **both** Newton and Linear, while `modalDamping` gives the correct +ζ.
>
> **The error is Δt-independent**, which rules out a lag/explicit-stability artifact:
>
> | steps/period | `modalDamping` (matrix) | `modalDampingQ` (force-only) |
> |---|---|---|
> | 100 | +0.04996 | −0.04972 |
> | 200 | +0.04999 | −0.04989 |
> | 400 | +0.05000 | −0.04995 |
> | 800 | +0.05000 | −0.04998 |
>
> Q converges to −0.05, not toward +0.05 → a **sign inconsistency**, not numerics.
>
> **Confirmed pure sign inversion:** `modalDampingQ(−0.05)` damps correctly at **+0.05003**.
> So the force-only path injects `+Cv` energy instead of dissipating `−Cv`.
>
> **Where it is NOT:** the residual sign convention is fine —
> `FE_Element::addRIncInertiaToResidual` does `theResidual += −1.0·getResistingForceIncInertia`
> (FE_Element.cpp:517), and the modal `−Cv` from `addModalDampingForce`/`setB`
> (IncrementalIntegrator.cpp:502/556) matches it. The two earlier force variants (:303, :349)
> are commented out. Both `modalDamping` and `modalDampingQ` add that SAME `−Cv` force in
> `formUnbalance` (TransientIntegrator.cpp:135); the ONLY difference is `modalDamping` *also*
> adds `+c·C` to the tangent (`addModalDampingMatrix` :563, `formTangent` :88).
>
> **Conclusion:** a velocity-proportional force applied through the *residual only* (no tangent
> `C`) is integrated by Newmark with the opposite effective sign — i.e. `modalDampingQ` as a
> standalone force-only mode appears to have never worked; only `modalDamping` (matrix + force)
> is correct. A naive "flip the force sign" would break `modalDamping` (shared force), so a fix
> must special-case `inclMatrix==false`. Encoded as a `strict` xfail
> (`test_modalDampingQ_force_only_matches_matrix`). **Use `modalDamping`, never `modalDampingQ`.**
> Candidate genuine upstream bug — surfaced, NOT patched (change-budget rules).
> *(This corrects an earlier assumption that the two forms give identical results.)*

---

## 4. Element `Damping` objects (frequency-band, "modern")

A separate hierarchy in [`SRC/damping/`](../SRC/damping/Damping.h) (Yuli Huang et al.,
2020–2022; *J. Earthquake Eng.* 2022). **Not Rayleigh** — targets a near-constant ζ over a
**frequency band** via a rational-filter approximation, evolved as committed state.

Created with the `damping` command ([`TclCommand_addDamping`](../SRC/damping/TclDampingCommand.cpp)):
- **`Uniform`** — `damping Uniform $tag $zeta $freq1 $freq2` — constant ζ across the band.
  `<-activateTime>/<-deactivateTime>/<-fact tsTag>` switch it on/off in time (e.g. no damping
  during gravity staging).
- **`SecStif`** — secant/committed-stiffness-proportional.
- **`URD` / `URDbeta`** — user-defined multi-point ζ(freq) least-squares fit.

**Attaches** via `-damp $tag` on the element command, or `region $r -damp $tag`
([`MeshRegion::setDamping`](../SRC/domain/region/MeshRegion.cpp) :324). Only elements that
override `setDamping` support it — base [`Element::setDamping`](../SRC/element/Element.cpp)
(:204) just warns. Supporting set: ElasticBeam2d/3d, DispBeamColumn, ForceBeamColumn, Brick,
FourNodeQuad, the Shell family, ZeroLength.

**Math** (e.g. [ElasticBeam2d](../SRC/element/elasticBeamColumn/ElasticBeam2d.cpp) :709, :1007):
the element calls `theDamping->update(q)`, adds `getDampingForce()` to its resisting force, and
**multiplies its stiffness by `getStiffnessMultiplier()`** — so it modifies both residual and
tangent. **ADDITIVE** with Rayleigh; a region `-damp` overwrites an element's prior `-damp`.

---

## 5. Material / dashpot damping

- **Rate-dependent uniaxials** — `Viscous` (σ=C·ε̇^α), `ViscousDamper` (Maxwell), Kelvin-Voigt,
  viscoelastic/viscoplastic. Damping is intrinsic to the stress-strain law → enters `R(u̇)`,
  **not** C. **This is the idiomatic way to damp a `zeroLength`**: put a `Viscous` material on a
  direction. Localized, physical, and **explicit-safe** (no critical-Δt penalty beyond the
  dashpot's own stiffness).
- **Material nonlinearity + element βK** — the element's `βK·K` uses the material tangent, so
  `betaK` (current) damping drifts as the material yields → another reason to prefer `betaK0`.

---

## 6. Numerical / algorithmic damping

HHT-α, generalized-α, KR-α, Bathe, bulk viscosity — controlled dissipation of high-frequency
modes via the integrator, **not** part of C. Out of scope here; see
[[04_explicit_dynamics_and_energy_balance]] and [[Ladruno_explicit_roadmap]].

---

## Implicit vs explicit assembly

**Implicit** ([`TransientIntegrator::formTangent`](../SRC/analysis/integrator/TransientIntegrator.cpp) :67):
1. zero A;
2. if `modalDamping` (incl-matrix) → `addModalDampingMatrix` (ch. 3);
3. DOF_Group tangents → `c1·K + c2·C_node + c3·M_node` (nodal Rayleigh, ch. 2);
4. FE_Element tangents → `c1·K + c2·C_ele + c3·M`, with the `Damping`-object stiffness multiplier already folded into K (ch. 1 & 4).

`c2,c3` are the integrator velocity/accel coefficients (Newmark γ/(βΔt), 1/(βΔt²)). C is in the
Jacobian → all six channels work and Newton converges. **All channels supported.**

**Explicit** ([`CentralDifference`](../SRC/analysis/integrator/CentralDifference.cpp) :155):
the "tangent" is the effective mass, `M_eff = c3·M + c2·C` with `c2=0.5/Δt, c3=1/Δt²`. Caveats:

> [!danger] βK is a trap in explicit
> Stiffness-proportional Rayleigh injects full stiffness into the effective system *and shrinks
> the critical time step* (Δt_cr ∝ 1/ω_max). **Prefer mass-proportional αM** in explicit.
> See [[../Ladruno_implementation/04_explicit_dynamics_and_energy_balance]] and the
> CentralDifferenceLadruno work. The `CentralDifference` integrator ctor can itself set
> domain-wide Rayleigh at `domainChanged` (CentralDifference.cpp:181) — another place Rayleigh
> originates, overwriting element values.

- Modal-matrix form fights the diagonal-mass advantage; the force-only `modalDampingQ` shape is
  the natural explicit fit — **but see the anti-damping bug above**.
- Material/dashpot damping (ch. 5) is cleanest for explicit.

---

## ⚠️ ZeroLength gotcha: `-doRayleigh` defaults OFF — VERIFIED

> [!warning]
> A `zeroLength` (and `zeroLengthSection`) element does **NOT** contribute
> **stiffness-proportional** Rayleigh damping (`betaK`/`betaK0`/`betaKc`) unless built with
> `-doRayleigh 1`. The default `-doRayleigh 0` yields **exactly zero** stiffness damping.
> Mass-proportional `alphaM` (which lives on the node) works regardless.

```tcl
element zeroLength 1 1 2 -mat 1 -dir 1                ;# doRayleigh=0 → βK0 does NOTHING
element zeroLength 1 1 2 -mat 1 -dir 1 -doRayleigh 1  ;# βK0 now damps
```

Pinned by `test_zeroLength_doRayleigh_default_off` (ζ≈0 with default) and
`test_betaK0_realises_target_zeta` (ζ=0.05 with the flag). Most other elements default-on;
zeroLength is the trap. Logged in `LEDGER_quirks.md`.

---

## Interaction / double-assignment cheat-sheet

| Combination | Result |
|---|---|
| domain `rayleigh` + region/element `rayleigh` (same ele) | **OVERWRITE** — last call wins, NOT summed |
| nodal αM + element αM (same DOF) | both apply; double-counts only if the mass itself is double-counted |
| `rayleigh` + `modalDamping` | **ADDITIVE** — both in C |
| `rayleigh` + element `-damp` object | **ADDITIVE** — independent channels |
| `modalDamping` then `modalDampingQ` | last call wins (same Domain slot, `inclMatrix` flips) |
| `Viscous` zeroLength + `rayleigh` on it | **ADDITIVE** (one is R(u̇), other is C·u̇) |
| `betaK` on zeroLength, default `-doRayleigh 0` | **ZERO** — needs `-doRayleigh 1` |
| `betaK*` on zeroLength + `-doRayleigh 1` | works |
| `eigen` not run before `modalDamping` | **silently no damping** (warning only) |
| `betaK` (tangent) on softening material | tangent can vanish/negate → instability; use `betaK0` |
| `modalDampingQ` (force-only), free-vib SDOF | **anti-damps in test** — flagged bug |

---

## Verification battery

`tests/test_damping_channels.py` — free-vibration log-decrement, implicit Newmark
average-acceleration (no algorithmic damping to contaminate ζ). **8 passed, 1 xfail.**

| Test | Asserts |
|------|---------|
| `test_alphaM_realises_target_zeta` | nodal αM → ζ=αM/(2ω) |
| `test_betaK0_realises_target_zeta` | element βK0 (`-doRayleigh 1`) → ζ=βK0·ω/2 |
| `test_zeroLength_doRayleigh_default_off` | default `-doRayleigh 0` → ζ≈0 (regression guard) |
| `test_element_rayleigh_overwrites_domain` | domain-then-`-ele` overwrites (0.06, not 0.08) |
| `test_rayleigh_plus_modal_is_additive` | Rayleigh + modal → sum (0.07) |
| `test_modal_per_mode_value` | 2-DOF: each mode damps at its own factor |
| `test_modal_scope_zero_mode_is_undamped` | mode with factor 0 stays undamped (no spillover) |
| `test_modalDamping_matrix_realises_target` | modalDamping (matrix) → ζ=0.05 |
| `test_modalDampingQ_force_only_matches_matrix` | **xfail** — modalDampingQ anti-damps |

> [!note] Test environment
> The built `dist/bin/opensees.pyd` is linked against **python312.dll** — run with
> `…\pythoncore-3.12-64\python.exe` (has numpy + pytest). `opensees_venv` is Python 3.11 and
> **cannot** load the pyd. Bootstrap = `os.add_dll_directory(DIST); sys.path.insert(0, DIST)`;
> no oneAPI setvars needed to *run*. See the memory note `project_opensees_test_env`.
>
> ```
> pythoncore-3.12-64\python.exe -m pytest tests\test_damping_channels.py -v
> ```

---

## Open items

> [!success] modalDampingQ root cause — characterized
> **Resolved at behavioral level:** confirmed *pure sign inversion* (`modalDampingQ(−ζ)` → +ζ),
> Δt-independent, both Newton & Linear. The `−Cv` force is correctly signed for the residual
> convention; the failure is that applying it through the residual *without* the `+c·C` tangent
> term (which only `modalDamping` adds) integrates with the opposite effective sign. So the
> standalone force-only mode appears to have never worked.
> **Decision (2026-06-01): DOCUMENT-ONLY.** No fork patch, no upstream report — the
> `modalDamping` (matrix) form fully covers the use case, so the workaround is sufficient.
> Pinned as a `strict` xfail so a future upstream fix is detected automatically (the xfail
> would start XPASSing). If revisited, a fix must special-case `inclMatrix==false`; a blanket
> force-sign flip is wrong (it would break `modalDamping`, which shares the force).

> [!question] Coverage gaps
> - Element `Damping` objects (ch. 4) are not yet in the battery — add a `Uniform` band-damping case.
> - Explicit-side: a CentralDifference βK critical-Δt collapse demo (theory known from CDL; not yet a pinned test).
