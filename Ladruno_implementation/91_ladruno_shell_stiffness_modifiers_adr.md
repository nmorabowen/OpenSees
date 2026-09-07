# ADR 91 — ETABS-style shell section stiffness modifiers

**Status:** Accepted (WP-91, branch `wp/91-shell-modifiers`)
**Author:** N. Mora-Bowen
**Class tag:** `SEC_TAG_LadrunoShellModifier 33000`
**Related:** ADR-87 (PR/branch protocol)

---

## 1. Context

Elastic shell models of RC shear walls need *cracked-section* stiffness. Code practice
(ACI 318-25 §6.6.3.1.1, ASCE 41 Table 10-5) reduces in-plane and out-of-plane wall
stiffness by different factors — typically 0.35–0.70 — and the reduction is **not
isotropic**: a wall may keep its vertical (gravity-carrying) membrane stiffness while its
in-plane shear and its bending stiffness are cracked away.

ETABS expresses this with ten *area section property modifiers* applied directly to the
shell's 8-component section stiffness matrix:

| ETABS | meaning | OpenSees plate resultant |
|---|---|---|
| f11, f22 | membrane direct | `SECTION_RESPONSE_FXX` (0), `FYY` (1) |
| f12 | membrane **shear** | `FXY` (2) |
| m11, m22 | bending | `MXX` (3), `MYY` (4) |
| m12 | twisting | `MXY` (5) |
| v13, v23 | transverse shear | `VXZ` (6), `VYZ` (7) |
| mass | inertia | section rho |
| weight | self-weight | *(see §5 — deliberately not implemented)* |

The mapping onto the OpenSees plate resultant ordering is 1:1.

**What the fork has today:** `ElasticMembranePlateSection` carries exactly ONE modifier,
`Ep_mod` (upstream contribution, "Out-of-Plane stiffness modifier added by Pearl Ranchal /
Degenkolb"). It is a single scalar substituted for `E` in the bending modulus and in the
transverse-shear correction. It cannot express f11 != f22, it cannot touch the membrane
block at all, and it is hard-wired into that one elastic class — `LayeredShell`,
`PlateFiber`, and the Nunez membrane sections have nothing.

## 2. Decision

Add **one** fork class: `LadrunoShellModifierSection`, a `SectionForceDeformation`
**decorator** that wraps any order-8 plate section and applies eight stiffness modifiers
plus a mass modifier.

Decorator, not a new elastic section, because:

- it works with `LayeredShell`, `PlateFiber`, `RCLayeredMembraneSection` and anything
  future, not just the one elastic class;
- the vanilla footprint is registration lines only — no upstream class is restructured;
- every shell element in the fork is already compatible. `ASDShellQ4.cpp:149` fetches a
  `SectionForceDeformation*` by tag and copies it; `ShellMITC4` / `DKGQ` / `NLDKGQ` /
  `ASDShellT3` only ever call `getStressResultant`, `getSectionTangent`, `getRho`. There is
  no `dynamic_cast` to a concrete section anywhere in `SRC/element/shell/`.

`ParallelSection` is the in-repo precedent for a wrapper section (notably its
`sendSelf`/`recvSelf` handling of a nested section).

## 3. The mathematical contract

Let `m[0..7] = (f11, f22, f12, m11, m22, m12, v13, v23)` and `s[i] = sqrt(m[i])`,
`S = diag(s)`. The wrapper is a **congruence transform** of the inner section:

```
eps_inner = S * eps                       (what the inner section is shown)
sigma     = S * sigma_inner(eps_inner)    (what the element is handed back)
D'        = S * D * S            i.e.  D'(i,j) = sqrt(m_i * m_j) * D(i,j)
rho'      = massMod * rho
```

`getSectionDeformation()` returns the **outer** strain.

Properties that motivated this form:

1. **SPD-preserving.** `D' = (s s^T) . D` is a Hadamard product with a rank-1 PSD matrix,
   so by the Schur product theorem it preserves definiteness for any non-negative
   modifiers.
2. **Exact consistent tangent for arbitrary nonlinear inner sections.** It is a true change
   of variables, so `dsigma/deps` is exactly `S D S` — no convergence penalty, and the
   decorator is transparent to `LayeredShell` and friends.
3. **Energy-consistent.** `W = 1/2 eps^T S D S eps`.
4. **Reduces to the obvious answer.** f11=f22=f12=c gives exactly `c * D_membrane`,
   matching the "scale E" intuition; and m11=m22=m12=v13=v23=r reproduces the upstream
   `Ep_mod=r` behavior.

**Sign note.** `ElasticMembranePlateSection` stores its plate block negative
(`tangent(3,3) = -D`, curvature sign convention). `S` is positive on both sides, so the
sign structure is preserved. This is not a bug; do not "fix" it.

## 4. Rejected alternative: diagonal-only scaling

Scaling only `D(i,i)` and leaving the Poisson coupling `D(0,1)` at full value **destroys
positive definiteness**. With f11 = f22 = 0.1 and nu = 0.2, the membrane block becomes

```
M * [[0.1, 0.2], [0.2, 0.1]]      -> eigenvalues of opposite sign
```

which is exactly the cracked-wall regime this feature exists to serve. Rejected.

The accepted trade-off: under the congruence, cross terms move as `sqrt(m_i m_j)`, so
f11=1 / f22=0.01 drops the nu-coupling by 10x, not 100x. If literal ETABS-bit-matching on
the coupling term is ever demanded, add a `-diagonal` opt-in **with an SPD guard that
refuses** rather than silently emitting an indefinite section. Demand-driven; not built.

## 5. Scope and non-goals

**No weight modifier.** Owner decision, and independently forced by the code: in
`ShellMITC4.cpp:1725-1755` the self-weight body force is `momentum = appliedB; momentum *=
rhoH` with `rhoH = materialPointers[i]->getRho()` — the *same* rho that builds the mass
matrix. A `weightMod` could therefore not be independent of `massMod`, and shipping an
argument that silently aliases another one is worse than not shipping it. Scale self-weight
at the load level instead (the `eleLoad -type -selfWeight` factors, or the gravity pattern).

Also out of scope: DDM/sensitivity methods (base defaults retained); exposing the modifiers
themselves as `setParameter` targets; aggregate shorthand flags (`-membrane`, `-bending`) —
they add a parse-order trap for no ETABS-parity gain.

## 6. Named refusals

- **R1** inner section `getOrder() != 8`
- **R2** inner section's `getType()` is not exactly `{FXX, FYY, FXY, MXX, MYY, MXY, VXZ, VYZ}`
  (an order-8 section with other response codes is not a plate section)
- **R3** any modifier `< 0`
- **R4** inner section tag not found

A modifier of exactly `0.0` is **allowed** (ETABS permits it) but warns once, naming the
mode made singular.

## 7. Interface

```
section LadrunoShellModifier $tag $innerSecTag \
    [-f11 v] [-f22 v] [-f12 v] [-m11 v] [-m22 v] [-m12 v] [-v13 v] [-v23 v] [-mass v]
```

All nine optional, default `1.0`, order-independent. All-defaults is byte-identical to the
inner section — that is gate G1.

Modifiers are per-section-tag. ETABS's per-area-object overwrite is expressed by defining
one wrapped section per distinct modifier set.

## 8. Implementation map

New: `SRC/material/section/LadrunoShellModifierSection.{h,cpp}`

Registration (vanilla files, each edit marked `// Ladruno ...`):
`SRC/classTags.h`; `SRC/material/section/{CMakeLists.txt,Makefile}`;
`SRC/interpreter/OpenSeesSectionCommands.cpp`;
`SRC/material/section/TclModelBuilderSectionCommand.cpp` **and**
`SRC/runtime/commands/modeling/section.cpp` (this fork has two Tcl paths);
`SRC/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp`;
`SRC/runtime/runtime/TclPackageClassBroker.cpp`.

`33000` is free in the SEC_TAG registry (highest vanilla `SEC_TAG` is 7703). Per-registry
banding means 33000 is independently used by `NUMBERER_TAG`, `RECORDER_TAGS` and
`INTEGRATOR_TAGS` — not a collision.

Recorders need no work: `SectionForceDeformation::setResponse` handles `force`/`deformation`
generically, so `section force` reports the **modified** resultants — which is the ETABS
reporting convention.

## 9. Verification

| Gate | Assertion |
|---|---|
| G1 identity | all modifiers 1.0 -> element K and R identical to the bare inner section (`ShellMITC4` + `ASDShellQ4`) |
| G2 congruence | single element, imposed uniform section deformation, 8 resultants vs numpy `S D S eps`, all nine modifiers distinct |
| G3 analytic | uniform membrane 0.5 exactly doubles in-plane tip deflection; uniform bending 0.5 exactly doubles out-of-plane |
| G4 SPD | f11=f22=0.01, nu=0.25 -> all global eigenvalues positive; numpy companion shows diagonal-only would be indefinite |
| G5 `Ep_mod` parity | m11=m22=m12=v13=v23=r reproduces upstream `Ep_mod=r` to round-off |
| G6 mass | massMod=0.5 -> frequencies scale by sqrt(2), stiffness untouched |
| G7 refusals | R1-R4 each error cleanly; modifier 0.0 accepted |
| G8 passthrough | `LayeredShell` wrapped at 1.0 -> identical nonlinear run |
| G9 frame equivalence | flexure-controlled cantilever (L/d=10) built as a frame (`A,I x0.35`) and as a shell mesh (`f11,f22,f12 x0.35`): both soften by exactly 1/0.35 to 8 sig figs, and the shell/frame ratio is identical gross and cracked -- the modifier scales THROUGH the discretisation gap |
| G10 in-plane no-op | `m11,m22,m12` and `v13,v23` change in-plane flexure by 0.0000%: out-of-plane modifiers are a silent no-op on a wall loaded in its own plane |

## 10. ETABS cross-validation

Reference model driven through **apeETABS** (`github.com/nmorabowen/apeETABS`, ETABS 22
OAPI via comtypes). Its read stack is complete; its `Assign.modifiers` write path was a
documented stub and is being implemented as a prerequisite so the reference model is
scripted and rerunnable rather than hand-built.

Note there is **no ETABS MCP connector** — the OAPI reference is the `etabs-oapi` skill
bundled in `apeETABS/.claude/skills/`.

The ETABS OAPI area modifier array carries 10 entries including `weight`; the campaign pins
`weight = 1.0` throughout, consistent with §5.

## 11. Open questions

- **OQ-1** Does ETABS apply the modifier to the Poisson coupling term as `sqrt(f11 f22)`
  (congruence) or leave it untouched? The CSI manual does not say. §4 argues the congruence
  is the only numerically admissible reading; the cross-validation in §10 is what will
  actually settle it. If ETABS turns out to differ, the discrepancy is confined to the
  off-diagonal membrane term and shows up only when f11 != f22.
- **OQ-2** Whether the same decorator should be offered for frame sections (ETABS also
  modifies Area/As2/As3/J/I22/I33). Not in WP-91.
- **OQ-3** The scaled work buffers (`static Vector sigma`, `static Matrix D`, and the
  function-local `static Vector eIn` in `setTrialSectionDeformation`) follow the upstream
  section idiom — `ElasticMembranePlateSection` does the same. Deliberate: a per-instance
  8x8 `D` would cost ~512 B per integration point, i.e. hundreds of MB on a large shell
  model. It is nevertheless the static-idiom hazard catalogued by ADR-75b Lane 3
  (threaded assembly). If the element loop is ever threaded, this class must be converted
  along with the rest of the section library, not before it.
