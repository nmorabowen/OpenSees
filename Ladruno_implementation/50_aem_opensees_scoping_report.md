# Applied Element Method (AEM) — scoping report against OpenSees / Ladruno

**Author:** research synthesis for N. Mora Bowen (Ladruno fork)
**Date:** 2026-06-23
**Status:** scoping / literature review (not an ADR — informs a future decision)
**Method:** deep-research web harness (19 sources, 72 claims, 25 adversarially
verified) + direct reading of primary sources in the personal Seafile library
(JSCE/JNDS originals + ELS Theoretical Manual V9 + Munjiza FDEM book).

> **Verification legend.** ✅ = adversarially verified (web) or read directly from a
> primary source in Seafile. 🟡 = credible lead, primary-source-plausible but the
> automated verifier abstained (rate-limited) — *not* refuted, just not machine-
> confirmed this run. ❌ = genuinely refuted by a verifier. See §5 for the audit.

---

## 0. Executive summary

The **Applied Element Method (AEM)** discretizes a body into **rigid virtual
elements** connected at discrete edge points by pairs of **normal + shear springs**.
All deformation lives in the springs; elements only translate and rotate (3 DOF in
2D). It is a descendant of Kawai's **Rigid-Body-Spring Model (RBSM)**, distinct from
**DEM** (Cundall–Strack, particles), **DDA** (Shi, blocks + contact constraints),
and **FDEM** (Munjiza, FE mesh + discrete fracture/contact). Its selling point is a
*single* formulation that tracks a structure continuously through the **continuum
stage** (elastic → cracking → yield → buckling) into the **discrete stage**
(spring rupture → element separation → debris fall → collision/recontact) **without
remeshing and without a geometric stiffness matrix**. This is what makes it the
engine of choice for progressive-collapse / blast / demolition simulation
(commercial product: **Extreme Loading® for Structures, ELS**, by Applied Science
International — the AEM authors' company).

**For the Ladruno fork**, the honest finding is: AEM is *not* a contact feature you
bolt onto FEM — it is a **different discretization**. But the *capability* AEM is
prized for (separation + recontact + collapse) is reachable from an OpenSees base
through a **hybrid FEM + element-removal + contact** route, and there is a real
OpenSees precedent (Talaat & Mosalam's element-removal work, §4). Your existing
NTS + mortar contact is a genuine head-start on the hardest missing piece —
**but** the truly missing capability is *dynamic broad-phase collision discovery*
(finding contact pairs at runtime among hundreds of separated bodies), which
neither standard OpenSees nor the current Ladruno contact lane provides.

---

## 1. AEM theory & formulation

### 1.1 Discretization and spring stiffness ✅

The structure is divided into small **rigid elements**. Two adjacent elements are
connected at discrete points along their common edge by a **pair of springs** — one
**normal** (Kn) and one **shear** (Ks). Each spring represents the stiffness of a
material *volume* of dimensions `d × T × a`:

```
Kn = E · d · T / a        Ks = G · d · T / a
```

where `d` = spacing between springs, `T` = element thickness, `a` = length of the
representative area (the element dimension normal to the edge to the spring line),
`E` = Young's modulus, `G` = shear modulus.

- **2D element → 3 DOF** (two translations + one rotation): pure rigid-body motion.
  The element *shape never changes*; every internal deformation is a spring
  deformation. The *assembly* is deformable. *(Verified ✅ against the primary
  source, Meguro & Tagel-Din 2002, JNDS 24(1); corroborated by Wikipedia and the
  Kharel formulation notes. 2-0 vote.)*
- **3D element → 6 DOF**, spring sets of **one normal + two shear** springs, giving a
  **12×12** pair stiffness matrix. *(🟡 — the web verifier abstained on the 3D
  details; **confirmed by direct reading of the ELS Theoretical Manual V9 §2**,
  which uses 8-node hexahedral solid elements with exactly this spring layout.)*

When reinforcement crosses a spring line, the **rebar stiffness is added** to the
matrix-material spring stiffness at that location (smeared at the spring, not a
separate truss). *(ELS Manual §2.5.2, §3.5.2.)*

### 1.2 Stiffness-matrix assembly ✅

The pair stiffness matrix is built by applying a **unit displacement to each DOF**
(holding the rest fixed) and summing the surrounding-spring force contributions —
the classic direct-stiffness construction, but with the spring geometry (`θ`, `α`,
lever arm `L`) baked into a trigonometric 6×6 (2D) given as Eq. (2) of Meguro &
Tagel-Din 2000. The **global matrix is sparse and banded** like FEM, so a standard
sparse solver applies. Poisson's-ratio coupling was added by Tagel-Din & Meguro
(1998). *(Read directly from Meguro & Tagel-Din 2000, JSCE 17(1), Eq. 1–2.)*

### 1.3 Geometric nonlinearity — no geometric stiffness matrix ✅

AEM handles large displacement **without** a geometric (initial-stress) stiffness
matrix. Per increment it (1) updates element geometry from the computed
displacement increment, (2) **rotates the spring force vectors** into the new
configuration, and (3) carries a **geometrical residual** back into the next
iteration:

```
K · ΔU = Δf + R_m + R_G,     R_G = f − F_m
```

with `R_m` the material/state residual and `R_G = applied − assembled-spring force`.
A **small-deformation assumption holds within each increment**; geometric
nonlinearity accumulates across increments. *(Verified ✅ verbatim against Meguro &
Tagel-Din 2002, JNDS 24(1), §3.1, Eq. 3–4. 2-0 vote.)*

> **Caveats the authors themselves state:** constant load direction (no follower
> loads), small load increments required, and load-control diverges past a limit
> point (needs displacement/arc-length control). Relevant when judging robustness
> vs OpenSees' arc-length + line-search machinery.

### 1.4 Separation criterion ✅ (Seafile)

A spring (hence the local material connection) **fails / separates** when its strain
reaches a **separation strain** threshold (ELS Manual §3.1.8). On the concrete side
this couples to tensile cracking + a **failure-softening factor** (§3.1.18) and a
combined normal–shear **failure envelope** (§3.1.14–3.1.16, Fig 3-13); on the rebar
side to **bar rupture** (§3.1.16). Once *all* springs on an edge fail, the two
elements are **physically separated** and thereafter interact only through contact.
This is the qualitative difference from FEM element-erosion: separation is a
**spring-level, directional, partially-reversible** event, not a binary
kill-the-element flag. *(Read from ELS Theoretical Manual V9 §3; Fig 3-9 "elements
separate and re-contact again".)*

### 1.5 Automatic contact / recontact / collision 🟡✅(Seafile)

After separation, AEM **automatically re-detects contact** between bodies and
reconnects them with **contact springs** (normal + shear) generated on the fly.
ELS Manual §6 enumerates three discrete contact types:

- **Corner-to-Face** (§6.2.1) — a corner of one element pressing a face of another;
- **Edge-to-Edge** (§6.2.2);
- **Corner-to-Ground** (§6.2.3) — debris hitting the rigid ground.

with **contact stiffness** factors (§6.3), **energy dissipation during contact**
(§6.4, an unloading-stiffness factor governing rebound — Fig 6-7/6-8), and a
**per-contact-type stable time step** (§6.5). *(The web verifier could only reach a
NAFEMS webinar for this and abstained 🟡; **the ELS Theoretical Manual V9 §6 in
Seafile confirms the mechanism directly** ✅.)* This is the capability that
separates AEM from element-deletion FEM — the **recontact loop is intrinsic**.

### 1.6 Material models at springs ✅ (Seafile)

Springs carry full uniaxial constitutive laws: elastic, tension-only, bilinear,
multi-linear, reinforced-concrete (Maekawa-type concrete + a reinforcing-bar model
with Bauschinger), composite steel-concrete, masonry, bearing, and a user-defined
material hook (UDMM). *(ELS Manual §3.2–3.13.)* Conceptually close to OpenSees
`uniaxialMaterial` objects living on the springs.

### 1.7 Time integration (dynamic collapse) ✅ (Seafile)

Dynamic collapse is run by **implicit Newmark-β step-by-step integration** with a
per-step Newton solve (load/displacement control), lumped mass, and Rayleigh-type
damping; the contact time-step limiter (§6.5) keeps collision stable. Earthquake,
blast (pressure-time), and element demolish/construct loading are all staged in
**loading stages** (§4.1–4.5). *(ELS Manual §4, §2.6.2.)* Note this is **implicit**
— unlike LS-DYNA-style explicit erosion FEM.

### 1.8 Annotated bibliography

**Foundational AEM (primary — citations verified against the PDFs in Seafile):**

1. **Meguro, K. & Tagel-Din, H. (2000).** *Applied Element Method for Structural
   Analysis: Theory and Application for Linear Materials.* Structural Eng./
   Earthquake Eng., JSCE, **17(1)**, 21s–35s (= J. Struct. Mech. Earthquake Eng.,
   JSCE, No. 647/I-51). — The origin paper: springs, Kn/Ks, 3-DOF rigid elements,
   Poisson coupling, RC + confinement. ≈174 citations.
2. **Tagel-Din, H. & Meguro, K. (2000).** *Applied Element Method for Simulation of
   Nonlinear Materials: Theory and Application for RC Structures.* Structural Eng./
   Earthquake Eng., JSCE, **17(2)**, 137s–148s. — Nonlinear/cracking companion.
3. **Tagel-Din, H. & Meguro, K. (2000).** *Applied Element Method for Dynamic Large
   Deformation Analysis of Structures.* Structural Eng./Earthquake Eng., JSCE,
   **17(2)**, 215s–224s. — Dynamic large-deformation lane.
4. **Meguro, K. & Tagel-Din, H. S. (2002).** *Applied Element Method Used for Large
   Displacement Structural Analysis.* J. Natural Disaster Science, **24(1)**, 25–34.
   — Buckling/post-buckling; the "no geometric stiffness matrix" result; the I–VI
   application matrix. [PDF: jsnds.org/jnds/24_1_3.pdf]
5. **Meguro, K. & Tagel-Din, H. (2001).** *Applied Element Simulation of RC
   Structures under Cyclic Loading.* J. Structural Engineering, ASCE, **127(11)**,
   1295–1305. — The widely-cited ASCE cyclic-RC paper.
6. **Tagel-Din, H. (1998).** *A New Efficient Method for Nonlinear, Large
   Deformation and Collapse Analysis of Structures.* PhD thesis, Univ. of Tokyo. —
   The originating thesis.
7. **Applied Science International, LLC (2004–2022).** *Extreme Loading® for
   Structures — Theoretical Manual, Version 9.* Durham, NC. — The authoritative
   engineering reference for the production AEM: separation strain (§3.1.8),
   contact types (§6), material models (§3), time integration (§4). *(Seafile,
   proprietary — not on the open web.)*

**Lineage (so the families stay distinct):**

8. **Kawai, T. (1978).** *New discrete models and their application to seismic
   response analysis of structures.* Nuclear Engineering and Design, **48(1)**,
   207–229. — **RBSM**, AEM's direct ancestor.
9. **Cundall, P. A. & Strack, O. D. L. (1979).** *A discrete numerical model for
   granular assemblies.* Géotechnique, **29(1)**, 47–65. — **DEM** origin.
10. **Shi, G.-H. (1988/1992).** *Discontinuous Deformation Analysis* (PhD, UC
    Berkeley; Engineering Computations **9(2)**, 157–168). — **DDA**.
11. **Munjiza, A. (2004).** *The Combined Finite-Discrete Element Method.* Wiley. —
    **FDEM** reference book. *(Seafile: FEM/Teoria/.)*
12. **Meguro, K. & Hakuno, M. (1989).** *Fracture analyses of concrete structures by
    the modified distinct element method.* Structural Eng./Earthquake Eng., JSCE,
    **6(2)**, 283s–294s. — The Extended/Modified DEM (EDEM) AEM was reacting to.

**Review / follow-on (from the ELS Manual reference chain + web):**

13. Lincy Christy, D., Pillai, T. M. M., Nagarajan, P. (2018). *Analysis of masonry
    using AEM* — IOP Conf. Ser. Mater. Sci. Eng. **330**, 012117. (review-grade)
14. Grunwald, C., et al. (Fraunhofer EMI, 2018) — see §2. Salem, H., et al.
    (2011–2021 series); Sasani & Sagiroglu (2008); Park, Suk, Hong (2009);
    Helmy, Salem, Mourad (2012/2013); Attia & Salem (2021). (RC/steel collapse
    validation chain cited wholesale in the ELS Manual references.)

A fuller raw catalog of the Seafile AEM holdings is staged at
`_aem_seafile_bibliography.md`.

---

## 2. AEM vs FEM / OpenSees — capability comparison

### 2.1 The collapse spectrum (ELS Manual Fig 1-1) ✅ (Seafile)

A collapsing structure passes through a **continuum stage** then a **discrete
stage**:

```
Linear → Crack/Yield/Crush → Buckling/Post-buckling │ Separation → Debris (rigid) → Collision
└──────────────── continuum ───────────────────────┘ └──────────────── discrete ──────────────┘
```

The ELS manual's own framing: FEM is **accurate** in the continuum range but
**"not automated"** at element separation and **"time consuming"** at collision;
AEM spans the whole spectrum in one formulation. *(Read from ELS Manual Fig 1-1.)*
This is the marketing claim — see §2.4 for the verified nuance.

### 2.2 What AEM genuinely does that standard FEM does not

| Capability | AEM | Standard nonlinear FEM (incl. OpenSees) |
|---|---|---|
| Element separation | Intrinsic (spring rupture) ✅ | Element **erosion/deletion** (mass/energy lost; binary) |
| Recontact after separation | **Automatic** contact-spring regeneration ✅(Seafile) | Requires a general contact algorithm + broad-phase search; rare in implicit structural FEM |
| Debris–debris / debris–ground collision | Built-in (corner-face/edge-edge/corner-ground) ✅(Seafile) | Needs self-contact + collision detection; expensive |
| Geometric nonlinearity | No geometric stiffness matrix ✅ | Corotational / total-Lagrangian + geometric stiffness |
| Crack path | Emergent along spring lines (mesh-biased) | Smeared, XFEM, cohesive, or phase-field |
| Mesh dependence | Crack/separation **biased to element edges** | Smeared softening needs regularization; XFEM mesh-independent |

### 2.3 Verified AEM-vs-FEM benchmarks ✅

- **Grunwald, C., et al. (Fraunhofer EMI, 2018).** *Reliability of collapse
  simulation — comparing the finite and the applied element method at different
  levels.* Engineering Structures (Elsevier), PII S0141029618301986. A **controlled
  FEM-vs-AEM benchmark** against experiments at **three levels** — small-scale
  connections, mid-size building elements, full complex buildings — under blast and
  earthquake. *(Study design verified 3-0 ✅.)*
- **(2025) ELS-vs-SAP2000 RC-frame study**, J. Building Engineering, PII
  S2352012425001729. A **five-span five-bay RC frame**, **54 scenarios** (2 methods
  × 6 slab models × 3 thicknesses × 3 damping ratios), three-stage corner-column
  removal at 3-s intervals — a current, rigorous AEM(ELS)-vs-FEM(SAP2000)
  head-to-head. *(Verified 3-0 ✅.)* Useful as an external reference to validate any
  OpenSees-based AEM-like capability against.

### 2.4 Honest nuance — what was *refuted* ❌

The claim that **"FEM fails to produce realistic collapsed shapes at full scale
while AEM succeeds across all scales"** was **refuted (0-2)**. The Fraunhofer study's
*existence and design* are confirmed; its *ranking conclusion is not* — do not
repeat "AEM beats FEM at collapse" as established fact. The defensible statement is:
AEM provides a **single-formulation, lower-setup-cost** path through separation +
collision that conventional implicit FEM reaches only with significant additional
machinery — a **workflow / automation** advantage, not a proven accuracy supremacy.

### 2.5 Family boundaries (taxonomy)

- **RBSM** (Kawai): rigid elements + interface springs, *small deformation* — AEM's
  parent. AEM ≈ RBSM + large-displacement + separation + recontact.
- **DEM** (Cundall–Strack): rigid *particles*, contact-governed, explicit; great for
  granular/debris but not a continuum stress field.
- **DDA** (Shi): deformable *blocks* with rigorous contact constraints (no
  interpenetration via open-close iteration); geomechanics roots.
- **FDEM** (Munjiza): a *real FE mesh* per body + cohesive fracture + DEM contact —
  the most "FEM-faithful" cousin; closest in spirit to a hybrid OpenSees route.
  AEM trades FDEM's mesh fidelity for spring simplicity and speed.

---

## 3. Feasibility for the Ladruno / OpenSees fork

> **Verification status:** the web harness could **not** machine-verify any
> Angle-3 claim this run (every OpenSees-specific lead hit rate-limiting and
> abstained 🟡 — *not* refuted). The analysis below is grounded in (a) the ELS
> Manual's spec of what AEM needs, (b) the Talaat & Mosalam OpenSees precedent
> (a real, well-known body of work — primary-source check recommended), and
> (c) OpenSees `Domain` architecture. Treat as an informed assessment, not a
> verified one.

### 3.1 Core finding: AEM is a different discretization, not a contact feature

You cannot turn OpenSees into AEM by adding a contact element. AEM replaces the FE
shape-function machinery with rigid-element + edge-spring kinematics. **But** the
*capability* you'd want from AEM — separation + recontact + progressive collapse —
is reachable as a **hybrid: nonlinear FEM (continuum stage) + runtime element
removal (separation) + general contact with broad-phase collision (discrete
stage)**. That hybrid is exactly the FDEM philosophy and is more natural on an
OpenSees base than a from-scratch AEM kernel.

### 3.2 The three architectural gaps

1. **Runtime topology change.** The OpenSees `Domain` is not designed for elements/
   nodes appearing and disappearing mid-analysis. Removing an element means:
   deleting it + its loads/constraints from the `LoadPattern`s, updating lumped
   nodal mass, and recursively cleaning up **dangling nodes** and **floating
   (rigid-body-unstable) sub-assemblies**. **Precedent 🟡:** *Talaat, M. & Mosalam,
   K. M. (2009), "Modeling progressive collapse in reinforced concrete buildings
   using direct element removal," Earthquake Engng Struct. Dyn., DOI
   10.1002/eqe.898* (+ **PEER 2007/10**) implemented exactly this in OpenSees via a
   `Recorder`-derived **element-removal class** invoked after each converged step
   (it has `Domain` access, checks per-element removal criteria, and mutates the
   `Domain`). The sudden release of the removed element's internal forces is
   imposed as **dynamic equilibrium at the end nodes** (correct transient
   acceleration) — superior to the "multiply stiffness by ε" trick (ill-conditioning
   + wrong dynamics). This is **the** OpenSees-native template to study first.
2. **Dynamic broad-phase collision discovery — the real missing piece.** Your NTS +
   mortar contact assumes **predefined contact pairs** (a known master/slave
   surface set). AEM/collapse needs to **discover** new contact pairs at runtime
   among hundreds of separated, tumbling bodies — i.e. a **broad-phase** stage
   (spatial hashing / sweep-and-prune / BVH) feeding a **narrow-phase** (corner-
   face / edge-edge / corner-ground, à la ELS §6). **Neither standard OpenSees nor
   the current Ladruno contact lane has this.** Notably, Talaat & Mosalam
   **explicitly scoped general collision *out*** (only a simplified "soft-impact"
   model, their Appendix B) — pinpointing this as the precise gap an AEM-like
   extension must fill.
3. **Solver robustness through separation events.** Spring rupture / element removal
   injects sharp transients; an implicit Newmark solve needs energy-consistent
   release + likely event-to-event stepping. Your viscous contact stabilization
   (`-visc`) and the fork's robust-solve track are relevant here.

### 3.3 Where Ladruno is already ahead

Your **mortar + NTS contact with friction and viscous stabilization** is a real
head-start on **narrow-phase** contact and on the constitutive side of collision
(normal/tangential enforcement, energy dissipation ≈ ELS §6.3–6.4). The missing
layer above it is **broad-phase discovery + contact-pair lifecycle management**
(create/destroy contact pairs as bodies separate and collide).

### 3.4 Prior art to verify before any build (open leads 🟡)

- **Talaat & Mosalam (2009 EQE / PEER 2007/10)** — element removal in OpenSees. The
  canonical precedent; the OpenSees wiki "Element Removal During Simulation" page
  documents the user-facing feature. **Verify the primary sources directly.**
- **Lu, X. (Tsinghua) FEM-DEM coupling** — failed FE quads removed at runtime and
  *replaced by granular DEM particles* (mass/volume-equivalent; nodal velocities
  interpolated to particles), implemented on MSC.MARC via `UACTIVE`/`FORCDT` user
  subroutines. A concrete model of the FEM→DEM runtime transition — though *not*
  OpenSees, and *not* AEM.
- Any **OpenSees + DEM (e.g. YADE/LIGGGHTS) co-simulation** — unverified; worth a
  literature pass if a DEM-debris route is of interest.

### 3.5 Recommendation (scoping-level)

1. **Do not** attempt a from-scratch AEM kernel inside OpenSees — wrong altitude;
   it fights the `Domain`/element abstraction.
2. **Do** treat "AEM-like collapse" as a **3-capability program** on the existing
   FEM + contact base: (i) **runtime element removal** (port/adapt the Talaat–
   Mosalam pattern), (ii) **broad-phase collision discovery** + contact-pair
   lifecycle (the genuinely new subsystem), (iii) **separation-event-robust
   integration** (build on the `-visc` / robust-solve work).
3. **Sequence:** element removal first (self-contained, real precedent, immediately
   useful for progressive-collapse-by-column-removal studies — and it dovetails
   with the existing `49_ladruno_integrator_study` work), broad-phase collision
   second (the hard, novel part), full debris recontact last.

---

## 4. Open questions / next steps

- **Primary-source-verify** Talaat & Mosalam (EQE 2009 + PEER 2007/10): exact class
  design, dangling-node/floating-element cleanup, mass/constraint update. (The
  feasibility analysis rests on this; the automated run could not confirm it.)
- Read **ELS Manual §6 in full** (Read pp. 105–112) to extract the actual
  corner-face / edge-edge contact geometry + the stable-time-step formula — the
  spec for a narrow-phase implementation.
- Confirm the **3D AEM** element-matrix details from a primary Meguro/Tagel-Din 3D
  paper (only 2D is web-verified; ELS Manual §2 corroborates 3D qualitatively).
- Decide scope: **column-removal progressive collapse** (needs only capability (i),
  achievable near-term) vs **full debris-field collapse** (needs (i)+(ii)+(iii), a
  multi-quarter program).

## 5. Methodology & verification audit

- Web harness: 5 search angles, 19 sources fetched, 72 claims extracted, 25
  adversarially verified (3-vote, 2/3 to kill). **4 confirmed, 21 "killed."**
- **Critical caveat:** the bulk of the 21 "killed" were `0-0 (3 abstain)` — the
  verifier agents hit **API rate-limiting**, so the votes never completed. These are
  **unverified, not refuted** (this includes *all* of Angle 3's OpenSees leads).
  Only **3 claims were genuinely content-refuted** (§2.4 ranking; one HW thesis
  validation detail; one blast-validation overreach).
- **Seafile primary sources materially upgrade the report**: the JNDS/JSCE
  originals confirm §1.1–1.3 first-hand, and the **ELS Theoretical Manual V9
  resolves §1.4–1.7 and §2.1** (separation strain, contact types, materials, time
  integration) that the web pass left at 🟡 — these are the strongest-tier evidence
  in this document.
- **Source quality:** lean on the JSCE/JNDS/ASCE originals + ELS Manual (primary);
  Engineering Structures / J. Building Engineering benchmarks (peer-reviewed); treat
  Wikipedia/blog/NAFEMS-webinar as corroboration only.

> **Bibliography gap to close (no rush):** the report would be stronger on the
> contact-algorithm depth (§3.2 broad/narrow phase) with **Wriggers,
> *Computational Contact Mechanics* (Springer, 2006)** and **Laursen,
> *Computational Contact and Impact Mechanics* (Springer, 2002)** — neither is
> currently in the Seafile library (only Kikuchi & Oden 1988).
