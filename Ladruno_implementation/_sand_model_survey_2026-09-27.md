# Sand constitutive models for the TIMs strip footing: physics, thermodynamic consistency, numerics

Survey for the Ladruno fork's TIMs soil work. Written 2026-09-27. Read-only research: web search and fetch
(abstracts, open PDFs, publisher metadata) plus the fork's own notes on `origin/ladruno`. No external code
was downloaded or run, and no repository or git state was changed. It builds on the earlier
`sanisand_external_survey.md` (implementations and integration practice for SANISAND) and does not repeat it.

**Tags on every claim**
- **E**: I read the source (full text, the relevant pages, or an authoritative abstract or metadata record).
- **E-sec**: known through a source that cites the original.
- **R**: from memory; I did not open the source in this session. Treat as unverified.
- **I**: my own inference, including arithmetic on published equations.

Many publisher pages (ScienceDirect, Wiley, Springer, ICE, cdnsciencepub) returned HTTP 403. Where the only
evidence was an abstract, the tag is E and the claim is limited to what the abstract says.

---

## 0. Executive summary

1. **The TIMs wall is an integration failure, not a model failure.** Changing the model does not fix it by
   itself. WP-128 and WP-134 traced the escapes of the integrated state to discrete defects in
   `ModifiedEuler`:
   - mechanism **F**: Λ < 0 is treated as elastic, which gives err = 0;
   - **E**: the error test looks at stress only;
   - **U9**: K and G are frozen over the increment;
   - **G**: α_in is reseated once per increment.

   The paper equations, integrated exactly from admissible ring starts, never escape (0 of 960 runs)
   **[E, local 134 §6.5]**. So the near-term question is whether DM04 is *physically* acceptable at the
   ring once it is integrated properly. The model question is a longer-term one.
2. **DM04 SANISAND's weaknesses are structural, but they are not the cause of the wall.**
   - It has no free energy:
     - The elasticity is hypoelastic, G ∝ √p with constant ν. That is non-conservative in closed cycles
       **[E-sec, Zytynski et al. via Houlsby et al. 2005]**.
     - The back-stress law dα = ⟨L⟩(2/3)h(α^b − α) with h = b0/((α−α_in):n) comes from no potential
       **[I]**.
   - Non-negative dissipation is therefore not guaranteed by construction. The fabric z and the α_in
     memory make this worse **[I]**.
   - On the physics side, its **power-law CSL keeps ψ, M^b and M^d bounded as p' → 0**. For the campaign set
     at e = 0.6944: ψ(p→0) = −0.136, M^b = 2.14 (a 52° triaxial friction angle) and M^d = 0.61 **[I, arithmetic
     on DM04 E4/E9]**. So DM04 represents the low-confinement ring with a finite, if high, strength. That is
     physically defensible.
3. **No realistic sand model is variational.**
   - Variational (incremental-potential) constitutive updates exist for associative critical-state
     plasticity, for example Ortiz & Pandolfi's variational Cam-clay **[E]**.
   - Sand needs non-associated dilatancy.
   - The thermodynamically admissible ways to get it are hyperplasticity with stress-dependent (frictional)
     dissipation or a kinematic constraint (Collins & Houlsby 1997; Houlsby 2019), or explicit restrictions
     such as Borja & Andrade's N̄ ≤ N.
   - These guarantee D ≥ 0 but in general give a **non-symmetric** tangent and no incremental minimum
     principle **[I; framework E-sec/E]**.
   - "Thermodynamically admissible + robust implicit return map + consistent tangent" is achievable.
     "Variational" is not, for a dilatant sand.
4. **Best physically consistent candidate for a monotonic dense-sand footing: NorSand in the Borja & Andrade
   (2006) form.**
   - Hyperelastic Houlsby-type energy (conservative).
   - State parameter ψ with image state.
   - Dissipation σ:ε̇ᵖ ≥ 0 proved under N̄ ≤ N.
   - Reduced dissipation inequality satisfied while p_i < 0.
   - Fully implicit return map in strain invariants with 3 local unknowns, and a closed-form consistent
     tangent **[E, read §2.1–2.7]**.
   - Widely used in practice: FLAC/3DEC/PFC built-in **[E]**, PLAXIS **[E]**, open GPL VUMAT **[E]**.
   - It needs two changes before TIMs use:
     - replace the logarithmic CSL, under which ψ → −∞ as p → 0 and dilatancy is unbounded at the free
       surface, with the curved or power-law CSL that the Itasca version already offers **[I; option E]**;
     - add a Lode-angle dependence for plane-strain strength.
5. **Ranked shortlist (§6):**
   1. NorSand-BA hyperelastic, implicit.
   2. DM04 SANISAND with the new SAS-ME/CPPM integrators and a hyperelastic upgrade.
   3. HySand, the 2026 multisurface hyperplastic sand model: the most rigorous thermodynamics, the best
      cyclic behaviour, and the least evidence and code.
6. **Near term.** Keep the calibrated DM04 with SAS-ME/CPPM and a documented p'-floor rule. This is
   defensible (§7.1). Even HySand's authors needed a **10 kPa surface surcharge "to ensure convergence at low
   stress levels"** in their 3D FE **[E]**. So a floor or surcharge is established practice even for a
   thermodynamically rigorous model.

   Known weaknesses of the near-term route:
   - hypoelastic energy non-conservation (small under monotonic loading);
   - the α law has no potential;
   - the UW `D_factor` sigmoid is a non-paper change that acts exactly in the ring (p < 5.05 kPa);
   - softening makes the post-peak response mesh-dependent without regularisation **[E, Gao et al. 2021]**.
7. **Longer term.** Implement a NorSand-BA-type class (curved CSL, 3-invariant, hyperelastic, implicit, with
   a Python kernel oracle first). About **5–6 engineer-weeks** including calibration from the same triaxial
   data **[I]**. Calibrate against DM04 on the same data, then on the strip deck. Keep HySand as a research
   option for the later cyclic SSI work.

---

## 1. The problem the model must serve (what matters for TIMs)

Source: the local intake and notes, **E**:
- Plane-strain rigid strip footing on dense sand, B = 1.5 m, pushed under displacement control; the target is
  the limit load or plateau.
- A 7.65 kPa surcharge outside the footprint.
- The DM04 campaign set: e_init 0.6944, e0 0.83, λc 0.027, ξ 0.45, Mc 1.3309, c 0.71, nb 3.5, nd 5.75, A0 0.05,
  h0 1.3, z_max 12.5, and the rest.
- The failure zone is a low-confinement ring at the free surface just outside the footing edge, with
  p' ≈ 0.3–10 kPa.
- PDMY reaches no plateau because its "dilation brake" is a critical-line crossing detector, not a
  state-dependent dilatancy **[E, local 133]**.

Physics requirements, ranked for a *monotonic limit load*:

| # | requirement | why | tag |
|---|---|---|---|
| P1 | Critical state and a state parameter ψ, so that dilatancy → 0 at critical state | a plateau and residual strength need it; PDMY's failure shows it | E (local 133) |
| P2 | Stress–dilatancy with state-dependent peak and softening of dense sand | N_γ depends strongly on density and stress level | E (Loukidis & Salgado 2011 abstract); I |
| P3 | Pressure-dependent stiffness | settlement before the peak | R |
| P4 | Consistent behaviour as p' → 0 | the ring. Real sand has no tensile strength, zero stiffness at p' = 0, and very high peak friction at low p' | R (Bolton 1986); I |
| P5 | Plane-strain strength (Lode dependence) | a strip footing is plane strain; φ_ps exceeds φ_tc | R |
| P6 | Fabric / inherent anisotropy | second-order for monotonic capacity, but neglecting it "may lead to significant overestimation" of N_γ for bedded sand | E-sec (Gao et al. snippets) |
| P7 | Cyclic capability | secondary, for later SSI | — |

**A property of DM04 worth stating explicitly [I].**
- Every DM04 modulus scales with √p: G and K through E2, and K_p ∝ p·h with h ∝ p^(−1/2). The rate equations
  are therefore pressure-homogeneous, and the dimensionless stiffness ratios stay bounded as p → 0.
- What shrinks is the elastic strain scale p/G ∝ √p. A given strain increment is a much larger *relative*
  load step at 3 kPa than at 100 kPa: a factor of about 5.8 in p/G.
- WP-134 found that the escape criterion is a *pressure ratio per increment*, Δp/p ≈ 6–8 **[E, local 134
  §6.4]**. That is consistent with this scaling. It is why the ring is where a stability-limited explicit
  scheme breaks.

---

## 2. Comparison table

Ratings: ++ strong, + adequate, 0 weak or partial, − poor or absent. Evidence tags as defined above.

| Model | (1) Physics for TIMs | (2) Consistency | (3) Numerics | (4) Calibration from DM04 data | (5) Availability / fork effort | (6) Footing evidence |
|---|---|---|---|---|---|---|
| **DM04 SANISAND** (Dafalias & Manzari 2004) | **++** ψ, CSL, M^b/M^d, Lode g(θ), fabric z. At p→0 the power-law CSL bounds ψ (M^b 2.14 on the campaign set) [I]. Cyclic ++ [E-sec] | **−** Hypoelastic, non-conservative [E-sec]. α law and α_in memory have no potential; D ≥ 0 not guaranteed [I]. The UW α_in rule gives h < 0, a *repelling* α law, in 146 of 480 exact runs [E, local 134] | **0 → +** Stiff, zero elastic range (m = 0.005), h singular at (α−α_in):n = 0 [E, local]. Robust with SAS-ME / CPPM + refusals [E, local plan]. Implicit line-search + complex-step Jacobian shown on a strip footing [E, Zhang et al. 2026 abstract] | **++** already calibrated | **++** in the fork (`LadrunoSANISAND`, `ManzariDafalias`); WP-129/130 in flight [E, local] | **+** Zhang 2026 strip footing (implicit) [E]; Roy et al. 2020 modified DM04 surface footing vs centrifuge, with "challenges" across strain levels [E]; Chaloulos, Papadimitriou & Dafalias 2019 strip footing [E metadata; content R] |
| DM04 + hyperelasticity (Lan, Liu & Zhao 2023 on SANISAND-MS; Irani et al. 2024 on DM) | as DM04 | **0** Elastic part conservative; no elastic ratcheting over 100 cycles [E abstracts]. Plastic part unchanged [I] | as DM04; hyperelastic tangent symmetric [I] | ++ plus the energy exponent [I] | + a modest change to the elastic block [I] | none found |
| SANISAND with cap (Taiebat & Dafalias 2008) | adds a cap for plastic compression under constant η; not needed for a monotonic footing [R/I] | − as DM04 [I] | 0 | + | 0 not in OpenSees [R] | none found |
| SANISAND-MS (Liu, Abell, Pisanò 2019) | cyclic ratcheting via a memory surface; monotonic ≈ DM04 [R] | − (hypo) / 0 (hyper, Lan 2023) [E] | 0; OpenSees uses RK4 with error control [E-sec, prior survey] | + | ++ `SAniSandMS` is in the fork's source [E, local tree] | suction caisson (cyclic) [E-sec] |
| SANISAND-F / -Z / -H | fabric evolution / zero elastic range / high pressure [E-sec] | − [I] | Z: implicit DIRK with substepping (2023) [E-sec] | 0 extra parameters | − | none found |
| Anisotropic critical state theory, ACST (Li & Dafalias 2012; Gao et al.) | ++ adds fabric anisotropy to ψ-based dilatancy [E-sec] | − (bounding surface, hypo) [I] | 0; **nonlocal void-ratio regularisation gives mesh-independent strip-footing curves** [E, Gao, Li & Lu 2021] | 0 needs anisotropy tests [I] | − | **++** strip footing, bearing capacity vs fabric [E metadata / E-sec] |
| **NorSand, Borja & Andrade 2006 form** | **+** ψ and image state, D* = αψ_i, softening on the dry side [E]. p→0: the **log CSL makes ψ unbounded** [I]; needs a curved CSL. No fabric; cyclic limited [R] | **+** Hyperelastic Houlsby-type energy, conservative [E]. D^p = σ:ε̇ᵖ ≥ 0 when N̄ ≤ N [E]. Reduced dissipation holds while p_i < 0 [E]. Additive free-energy split holds only after dropping an O(Δε_s^p²) term [E] | **++** Implicit return map: 3 unknowns (ε_v^e, ε_s^e, Δλ) plus a scalar sub-Newton on p_i, and a closed-form algorithmic tangent [E]. Non-symmetric when non-associated [I]. Needs a cap on Q near η → 0 [E] | **+** CSL, M, G transfer; H, χ, N refit from the same drained triaxials [I] | 0 not in OpenSees [I]; GPL VUMAT (explicit) [E]; FLAC/PLAXIS built-in (closed) [E] | 0/+ plane-strain localization of dense sand [E]; PLAXIS landslide with Prandtl-type bands [E, Woudstra 2021]; no bearing-capacity validation found |
| NorSand (Jefferies 1993; Jefferies & Been; Itasca) | as above plus principal-stress-rotation and cyclic terms [E, Itasca doc; E-sec Cheng & Jefferies 2020] | 0/+ built on the Cam-clay work equation [R]; elasticity power law "G = G_ref(p/p_ref)^m", not stated as hyperelastic [E] | not documented [E, doc silent] | + | − closed (Itasca); VUMAT GPL [E] | as above |
| **HySand** (Simonin, Houlsby & Byrne 2026) | **+** multisurface Matsuoka–Nakai, density lines (B, Γ, Δ), density- and anisotropy-dependent dilation, consolidation mechanism, peak and softening [E] | **++** Gibbs energy + yield surfaces (hyperplastic); convex yield surfaces imply D ≥ 0 [E] | + ABAQUS and PLAXIS UMATs exist [E]; **used a 10 kPa surcharge for convergence at low stress** [E]; integration in Simonin 2023 and Houlsby 2025a (unread) | 0 14 parameters, calibrated on monotonic + cyclic triaxials [E]; the CSL maps from DM04, the rest is refit [I] | − no public code found [E, none found]; multisurface state is N × (6+) variables [I] | 0 monopiles only (confidential validation) [E] |
| Other hyperplastic sand models (Golchin & Lashkari 2014; Orakci et al. 2026 preprint; 2023 C&G state-dependent dilatancy angle) | + ψ-linked dilatancy [E abstracts] | ++ by construction [E abstracts] | Orakci: "robust return mapping" near the apex [E] | 0 | − | none found |
| Continuous / critical-state hyperplasticity (Einav & Puzrin 2003–04; Coombs & Crouch 2011; Coombs 2017) | clay-oriented critical state; no sand ψ [E-sec/R] | ++ [E] | ++ implicit + consistent tangent for 3-invariant CS hyperplasticity [E, Coombs & Crouch 2011 abstract] | − (clay) | − | none |
| Manzari–Dafalias 1997 / Li–Dafalias 2000 / Li 2002 | ++ first ψ-dependent dilatancy (Li & Dafalias 2000) [R] | − [I] | 0; Li 2002 was integrated explicit vs implicit by Carow & Rackwitz 2021 [E, prior survey] | + | − | Loukidis & Salgado's DM-type two-surface model: N_γ [E] |
| Severn-Trent sand (Gajo & Wood 1999) | + MC failure, CSL, ψ-dependent strength and stiffness, kinematic hardening [E-sec] | − [I] | + "normalized stress space" eases implementation [E-sec] | + | − | footings and walls listed [E-sec] |
| CASM (Yu 1998) | + unified clay/sand, ψ [R] | 0 [R] | + [R] | + | − | none found |
| **Hypoplasticity** (von Wolffersdorff 1996 + Niemunis–Herle IGS; Fuentes ISA) | ++ density, pressure, CSL asymptotes; cyclic with IGS/ISA [E-sec]. p→0: stiffness → 0, needs a small-stress fix [R] | **− by design:** incrementally non-linear rate law with no free energy. A dissipation inequality can be imposed afterwards (Toll 2011) [E metadata; E-sec abstract] | + explicit adaptive RK is standard [E-sec, prior survey]; no yield surface to return to | 0 different parameter set (h_s, n, e_d0, e_c0, e_i0, α, β) [R] | 0 numgeo, soilmodels UMATs (registration) [E-sec] | **+** Tejchman & Herle 1999 "class A" plane-strain footings with polar hypoplasticity vs model tests [E-sec] |
| **PM4Sand** (Boulanger & Ziotopoulou) | + DR-based, plane strain, liquefaction-oriented; floors p at 0.5 kPa or 0.005·p_in [E, prior survey] | − [I] | 0 explicit (FLAC) / ME in OpenSees [E, local header] | 0 DR-centric, not a CSL fit [R] | **++** in the fork (`PM4Sand.cpp`) [E, local] | none found for monotonic footings |
| PDMY01/02/03 | − no saturation of dilation [E, local 133] | − | + | 0 | ++ | the TIMs deck: 417.6 kPa limit point at a 33° cone [E, local intake] |
| Breakage mechanics (Einav 2007) | crushing at high p; irrelevant in the ring [I] | ++ [R] | + [R] | − gradation tests | − | none |
| H-directional micromechanics (Nicot & Darve) | micro-directional; cost [R] | + [R] | 0 [R] | − | − | none |
| GSH / "thermodynamic CS model for sands" (Granular Matter 2025) | CS behaviour "without yield surface" [E abstract] | ++ [E abstract] | unknown | − | − | none |

---

## 3. SANISAND family

### 3.1 Physics

- **Critical state.** DM04 is critical-state-compatible through ψ-dependent bounding and dilatancy images
  (E9 in WP-134's table) **[E, local 134 §2]**.
- **Dense-sand behaviour.** Peak, softening towards the CSL, and phase transformation are all reproduced by
  WP-134's independent reference **[E, local 134 §6.1]**.
- **Low p'.**
  - With e_c = e0 − λc (p/pa)^ξ, ψ tends to e − e0 as p → 0. On the campaign set that gives
    M^b ≤ 1.3309·e^(3.5·0.136) = 2.14 and M^d = 0.61 **[I]**. The fork's intake quotes M^b ≈ 2.10 at the ring
    **[E, local 127 §0]**.
  - Strength stays finite and dilatancy stays bounded, which is the physically right behaviour: peak φ rises
    at low stress, but not without limit **[R, Bolton 1986]**.
  - Two terms are genuinely singular: G ∝ √p → 0, and b0 ∝ p^(−½) → ∞ **[E, E2/E10]**. In the continuum this
    is harmless because of homogeneity, but it makes explicit integration stiff (§1).
- **UW `D_factor`.** The UW implementation multiplies D by a sigmoid for p < 0.05·P_atm = 5.05 kPa. This is
  **not in DM04**, and it acts across most of the ring **[E, prior survey §4.3; local 134 U1]**. It is a
  physics choice: less dilatancy at the surface. It should be a documented decision, not an inherited
  default **[I]**.
- **Plane strain.** Handled by the Lode interpolation g(θ, c) **[E]**. c = 0.71 sets the extension/compression
  ratio. The plane-strain M then comes out between Mc and Me **[I]**.

### 3.2 Consistency

- **Elasticity.** Hypoelastic with ν constant and G ∝ √p. Such laws can generate energy in closed loops
  **[E-sec, Zytynski et al., via Houlsby et al. 2005 and search summaries]**. Hyperelastic replacements for
  DM-type models exist:
  - Lan, Liu & Zhao 2023 (SANISAND-MS), with volumetric–deviatoric coupling: no elastic ratcheting, and
    stress-induced anisotropy **[E, abstract]**;
  - Irani et al. 2024, who assessed energy potentials under Dafalias–Manzari plasticity: reversible under
    hyperelasticity, stress accumulation over 100 cycles under hypoelasticity **[E, abstract]**.
- **Plastic part.** No free energy and dissipation pair is published for DM04's α law, α_in memory, or fabric
  z **[I; none found]**.
  - Houlsby & Richards (2023) show that bounding-surface-type models can be built *within* hyperplasticity
    **[E, metadata + E-sec summary]**. That is a different model, not a proof for DM04.
  - The concrete danger is the UW once-per-increment α_in rule. Integrated exactly, it produced h < 0 in
    146 of 480 ring runs, turning dα ∝ h(α^b − α) into a repelling law, plus 27 rate problems with no
    admissible solution **[E, local 134 §6.5, §7]**. The paper's event-driven reseat avoids this entirely
    (0 of 960) **[E, local 134]**. That is a model-definition issue as much as an integration one.
- **No thermodynamic reformulation of SANISAND was found.**

### 3.3 Numerics

- Covered in depth by WP-128/134 and the earlier survey.
- Established remedies:
  - SAS-type error control over (σ, α, z), refusal instead of force-accept, a stiffness step limit, per-stage
    K and G **[E, local]**;
  - implicit alternatives: line-search Newton with a complex-step Jacobian for SANISAND-04 in ABAQUS,
    validated on a **strip footing** bearing-capacity problem (Zhang, Gao, Lu, Zhou, Lai & Du 2026)
    **[E, abstract]**.

  Nothing here is variational.

### 3.4 Footing evidence

- **Zhang et al. 2026.** Strip footing with implicit SANISAND-04 **[E, abstract]**.
- **Roy, Chow, O'Loughlin, Randolph & Whyte 2020 (CGJ).** A modified DM04 applied to a surface circular footing
  and a plate anchor against centrifuge tests: "potential" shown, with challenges across strain levels
  **[E, abstract]**.
- **Chaloulos, Papadimitriou & Dafalias 2019 (JGGE).** Fabric effects on strip footings with an anisotropic
  SANISAND-type model **[E metadata; content R]**.
- **Loukidis & Salgado 2011 (Géotechnique 61(2)).** FE with a two-surface (DM-type) critical-state model,
  which includes non-associated flow, softening and anisotropy. They evaluate N_γ and s_γ and propose a
  friction-angle selection rule **[E, abstract]**.
- **Gao, Li & Lu 2021.** With ψ-softening, strip-footing load–displacement and shear-band width are
  **mesh-dependent unless regularised**. A nonlocal void ratio gives mesh-independent curves when h < ℓ
  **[E, abstract]**.

---

## 4. Hyperplasticity family

### 4.1 Framework

- Two scalar potentials, an energy (Helmholtz or Gibbs) and a dissipation or yield function, define the
  whole model **[E-sec, Houlsby & Puzrin 2006 via search]**.
- Non-associated frictional flow is admissible when the dissipation depends on stress, or through a
  kinematic dilation constraint. Collins & Houlsby 1997 is the origin, and Houlsby 2019 recasts it in convex
  analysis with Fenchel duals **[E-sec / E, abstract]**.
- The yield surface in *true* stress then differs from that in dissipative (generalised) stress, which is the
  "shift stress" idea **[R; the term itself was not confirmed in anything I read]**.
- Consequence **[I]**: D ≥ 0 is guaranteed, but the incremental problem has no minimum principle, so the
  algorithmic tangent is generally non-symmetric.
- **Stored plastic work** (Collins 2005; Collins & Muhunthan 2003) explains dilatancy and the effect of
  density within this framework **[E-sec]**. It is the natural home for a state parameter.

### 4.2 Models

- **HySand** (Simonin, Houlsby & Byrne 2026, Géotechnique 76(13), doi 10.1680/jgeot.25.00091) **[E]**:
  - 14 parameters;
  - N Matsuoka–Nakai-type yield surfaces combined with consolidation surfaces;
  - dilation depending on density and anisotropy through an internal anisotropy variable;
  - density measured against loosest (B), critical (Γ) and densest (Δ) lines;
  - pressure-dependent Gibbs energy g ∝ p^(2−m)-type with m = 0.73 in the Karlsruhe calibration
    **[E, ISFOG 2025-301 Table 1–2]**;
  - reproduces the peak and softening of dense sand and undrained phase transformation on Karlsruhe fine
    sand **[E, SEG23 abstract; ISFOG figures]**;
  - implemented as ABAQUS and PLAXIS 3D user subroutines **[E, ISFOG 2025-506]**;
  - in the monopile FE, **K0 = 1 and a 10 kPa surcharge "to ensure convergence of the analysis at low stress
    levels"** **[E, ISFOG 2025-301 §4]**;
  - integration details are in Simonin (2023, DPhil) and Houlsby (2025a) **[E-sec]**, which I did not read;
  - no public code found.
- **Golchin & Lashkari 2014 (IJSS).** A critical-state sand model with elastic–plastic coupling
  (G ∝ p^χ) in a hyperelastic framework with bounding-surface plasticity **[E-sec, abstract via search]**.
  Its plastic part is bounding-surface, so its full thermodynamic status is **[unverified]**.
- **Orakci, Anoyatis, Chow & François 2026 (engrXiv preprint).** Adds a *third potential* to Ziegler's two
  for non-associativity, uses a ψ-linked dilatancy angle and a phase-transformation condition, "robust
  return mapping" near the apex, and triaxial and DSS validation **[E, abstract]**. Not peer-reviewed.
- **2023 C&G hyperplastic model with a state-dependent induced dilatancy angle and crushing** **[E-sec,
  snippet]**.
- **Continuous hyperplastic critical state (CHCS)** (Einav & Puzrin 2004) and the Einav–Puzrin–Houlsby 2003
  numerical studies (single, multiple and continuous yield-surface fields): smooth response, continuous
  memory, "guaranteed to obey the laws of thermodynamics" **[E-sec]**. These are clay-oriented.
- **Coombs & Crouch 2011 (CMAME 200).** Implicit stress integration and consistent tangents for 3-invariant
  critical-state hyperplasticity **[E, abstract]**. This is the numerics template for any hyperplastic CS
  model.

### 4.3 Verdict

- Hyperplasticity is the only family whose consistency is *structural* rather than argued case by case.
- For a monotonic footing, its sand models are young: HySand was published in 2026, with no footing studies
  and no public code.
- Multisurface formulations are expensive: N surfaces, each with tensor internal variables **[I]**.
- It is the right long-term home for cyclic SSI. It is not the fastest route to a credible TIMs limit load
  **[I]**.

---

## 5. NorSand

### 5.1 Physics

Jefferies' NorSand **[E, Itasca doc]**:
- Cam-clay-like "bullet" yield surface, sized by the image stress p_i.
- State parameter ψ = e − e_c. Hardening H = H0 − H_y ψ, towards a limit p_i,max set by D_min = χψ_i.
- Flow D^p = M_i − η.
- An internal cap so that unloading to very low p yields.
- Options for principal-stress-rotation softening.
- CSL either as Γ, λ (log) or C1–C3 (curved).

It captures the density-dependent peak and softening that P1 and P2 need. It has no fabric. Cyclic ability
comes from later add-ons (NorSand-PSR, NorSand-VT UMAT) **[E-sec / E, GitHub listing]**.

**p' → 0 [I].**
- With the log CSL, e_c → ∞ as p → 0, so ψ → −∞. The limiting dilatancy D* = αψ_i is then unbounded, which
  gives unbounded peak strength at the free surface. That is physically wrong, and a numerical hazard in
  exactly the TIMs ring.
- The curved (power-law) CSL, as in DM04 or Itasca's C1–C3 option, bounds it.
- In Borja & Andrade's hyperelastic energy, p = p0·exp(ω)·[…]. p reaches 0 only asymptotically, so the elastic
  law has no tension branch. μ0 > 0 adds a constant shear-stiffness floor that regularises the p → 0 limit
  **[I from their eqs. 2.3, 2.51]**.

### 5.2 Consistency

Borja & Andrade 2006, CMAME 195:5115, §2.1–2.5 **[E, read]**:
- Stored energy Ψᵉ(εᵉ_v, εᵉ_s) of the Houlsby 1985 type, "conservative" in a closed elastic loop.
- Non-associated volumetric flow through Q with N̄.
- **Plastic dissipation σ:ε̇ᵖ ≥ 0 provided N̄ ≤ N.**
- With p_i taken as the stress-like variable conjugate to ε_s^p, the reduced dissipation inequality
  σ:ε̇ᵖ − p_i ε̇_s^p ≥ 0 holds, because ε̇_s^p = λ̇ ≥ 0 and p_i < 0 (compression negative).
- They show that the free energy is **not** exactly additive: Ψᵖ depends on εᵉ through p_i*. They drop the
  O(Δε_s^p²) coupling term in the stress.

So the model is thermodynamically admissible up to that stated second-order approximation. This is the
strongest consistency statement I found for any *practice-grade* sand model **[I]**.

Jefferies' own (Itasca) NorSand uses a power-law G and does not state an energy **[E]**. Its admissibility
therefore rests on the work equation **[R, Jefferies 1997]**.

### 5.3 Numerics

Borja & Andrade 2006, §2.6–2.7 **[E]**:
- Classical return map in strain invariants.
- Local unknowns (εᵉ_v, εᵉ_s, Δλ), with p_i solved by a nested scalar Newton, or as a fourth unknown in one
  loop with "about the same" efficiency.
- A closed-form algorithmic tangent from the converged local Jacobian (their eq. 2.75). It is non-symmetric
  when N̄ ≠ N **[I]**.
- They also used B-bar near the critical state.
- A plastic-potential cap for η < χM (χ ≈ 0.1) avoids a corner or negative η at the hydrostatic axis. That is
  a known failure mode, relevant to the ring where η swings **[E, Remark 2]**.

Other implementations:
- Itasca documents no integration scheme **[E]**.
- The INL `GranularFlowModels` VUMAT (GPL-2.0, Abaqus/Explicit) states "explicit (sub-stepping)/implicit
  (return mapping)" updates **[E]**.

### 5.4 Calibration and evidence

- Calibration uses the same drained and undrained triaxials plus the CSL **[E-sec, Loukidis/Jefferies
  lineage; Woudstra 2021]**.
- In PLAXIS, dense soils showed indefinite undrained hardening, which needed a cavitation cut-off
  **[E, Woudstra 2021 abstract]**. That is irrelevant for the drained strip.
- No published strip-footing bearing-capacity validation was found **[I, none found]**.

---

## 6. Ranked shortlist (top 3) and the trade-off

| rank | model | why | the trade-off, plainly |
|---|---|---|---|
| **1** | **NorSand, Borja & Andrade hyperelastic implicit form**, with a curved/power-law CSL and a Lode-dependent M(θ) | Proven dissipation ≥ 0 (N̄ ≤ N) and a conservative hyperelastic energy **[E]**. A small, smooth implicit local problem with a closed-form tangent **[E]**. ψ-based peak and softening **[E]**. None of DM04's α back-stress, α_in memory, h ∝ 1/((α−α_in):n) singularity or zero-elastic-range cone, which are exactly the TIMs pathologies **[I]** | You give up fabric and good cyclic behaviour. Monotonic only in its clean form **[R/I]**. The published consistency proof is 2-invariant; 3-invariant M(θ) must be re-checked **[I]**. Not in OpenSees; about 5–6 weeks to build **[I]**. Footing evidence is thin **[I]** |
| **2** | **DM04 SANISAND with SAS-ME/CPPM, the paper α_in rule, and a hyperelastic upgrade** | Already calibrated and in the fork. The best cyclic path for SSI later. The oracle exists (WP-134). Implicit SANISAND-04 has been shown on a strip footing **[E, Zhang 2026]** | No free energy for the plastic part, so D ≥ 0 is not guaranteed **[I]**. The ring stays numerically stiff **[E, local]**. `D_factor` is a non-paper low-p physics change **[E, local]** |
| **3** | **HySand** (hyperplastic multisurface) | The strongest consistency by construction, with density/CSL, dilation and cyclic behaviour in one parameter set **[E]**. ABAQUS/PLAXIS UMATs exist **[E]** | Newest (2026), no public code, no footing evidence. Multisurface cost. Its authors still needed a 10 kPa surcharge at low stress **[E]**. Research risk; about 10–14 weeks **[I]** |

Not shortlisted:
- **PM4Sand** is the practice benchmark and is already in the fork, so run it as a *check* on the TIMs deck.
  It is plane-strain, liquefaction-oriented and not thermodynamic **[E/I]**.
- **Hypoplasticity** has good footing evidence (Tejchman & Herle) but is non-thermodynamic by design **[E-sec]**.
- **ACST** is the best *physics* for anisotropy and has the regularised strip-footing evidence **[E]**, but its
  thermodynamics are no better than DM04's.

---

## 7. Recommendation for TIMs

### 7.1 Near term: keep DM04 with SAS-ME / CPPM and a documented p'-floor rule

**Verdict: defensible**, on four conditions.

1. **Model-level α_in rule.** Adopt the paper's event-driven reseat, not UW's once-per-increment rule.
   - WP-134 shows the UW rule is itself a source of h < 0 and repelling α, even under exact integration
     **[E, local 134]**.
   - Without this, "better integrators" integrate a flawed rate problem accurately.
2. **Floor rule, stated as a regularisation and not as physics.**
   - Precedent: PM4Sand floors p at 0.5 kPa or 0.005·p_in, numgeo at 0.5 kPa, and the Tamagnini/Mašín UMAT
     adds 0–5 kPa of artificial cohesion **[E, prior survey]**. HySand's FE used a 10 kPa surcharge
     **[E]**.
   - Proposed rule **[I]**:
     - (a) The *elastic* floor (`-Pmin`) is at most 0.5 kPa.
     - (b) The *strength* shift (`-Presidual`) is 0 by default, and at most 1 kPa if used.
     - (c) Report the limit load at floor F and F/2. Accept the floor if the load changes by less than
       about 2 %.
     - (d) Report the number of Gauss points at the floor at the limit state.
   - Scale check **[I]**: an apparent cohesion c ≈ p_r·tanφ ≈ 0.7 kPa per kPa of p_r, times N_c ≈ 30–75 at
     φ 35–45°, gives about 20–50 kPa per kPa of p_r. Against about 650 kPa, **p_r = 1 kPa can move the load
     by up to about 3–8 %**. That is why (b) and (c) matter.
   - The act's own data: `-Presidual` 1.01 and 5.05 kPa did not move the wall outside run-to-run scatter
     **[E, local intake §1.4]**. That was under the defective integrator, so it must be re-measured.
3. **Decide `D_factor` explicitly.** The UW low-p dilatancy sigmoid, active below 5.05 kPa, is not DM04
   **[E]**. Run with it on and off, and report its effect on the ring and the load **[I]**.
4. **Expect mesh dependence after the peak.** ψ-softening localises. Report B/8 and B/16. Treat the post-peak
   plateau as mesh-dependent unless a regulariser is added, for example a Gao-type nonlocal void ratio
   **[E, Gao, Li & Lu 2021]**. The peak is less sensitive than the softening branch **[R]**.

**Known weaknesses of the near-term route, stated plainly:**
- No energy function for α, so dissipation positivity is not guaranteed. In practice it is watched, not
  proved **[I]**.
- The elasticity is hypoelastic: a small error under monotonic loading, a real one for cyclic SSI later
  **[E-sec / I]**. The optional fix is a Houlsby-type energy (n = 0.5), as Lan et al. 2023 and Irani et al.
  2024 did for DM-type models **[E]**.
- m = 0.005 makes the elastic range essentially zero, so almost every increment is plastic and stiff at low
  p **[E, local; I]**.
- p → 0 remains singular in G and b0 and is only regularised by the floor **[E/I]**.

### 7.2 Longer term: implement a NorSand-BA-type model ("LadrunoNorSand"; name and class tag to be reserved)

**Specification [I, from the sources above]:**
- Hyperelastic energy of the Borja & Andrade (Houlsby 1985) form, or Houlsby, Amorosi & Rojas (2005) with
  n ≈ 0.5 to match DM04's √p. Check the convexity limit at high stress ratio **[R]**.
- Curved CSL e_c = e0 − λc (p/pa)^ξ, **the same form as DM04**. This bounds ψ at p → 0, and its three
  parameters transfer unchanged.
- M(θ) from Jefferies or Matsuoka–Nakai-type Lode dependence (plane strain). Re-derive the dissipation proof
  of B&A §2.2 with θ-dependence.
- Image-state hardening and D* = χψ_i, with the internal flat cap and B&A's Q-cap near η → 0.
- Return map in invariants with an algorithmic tangent (B&A §2.6–2.7), with refusal codes forwarded to the
  element (F7 roster). Optional: substepping on local-Newton failure (Carow & Rackwitz; arXiv 2504.17476
  **[E, prior survey]**).
- Kernel-oracle doctrine: a Python reference (SciPy Radau, like WP-134) before any C++. FD-tangent gate,
  one-element drivers, a dissipation monitor (σ:ε̇ᵖ − p_i ε̇_s^p ≥ 0 asserted in tests), and a mutation gate.

**Fit to the fork [I].** A standalone class in the `LadrunoJ2` style is preferred over the
ASDPlasticMaterial3D kit.
- The kit has linear, Duncan-Chang and StiffSoil elasticity, and MC/DP/HB/StiffSoil yield functions and flow
  directions **[E, local tree]**, but no hyperelastic energy and no ψ-driven image hardening.
- Adding those as kit components is possible but touches its templates.

**Effort [I]:**

| stage | time |
|---|---|
| Python oracle | ~1 week |
| C++ class + parser + recorder responses | ~1.5–2 weeks |
| 3-invariant / plane-strain extension and proof check | ~0.5–1 week |
| Test battery (FD tangent, patch, one-element, dissipation) | ~1 week |
| Calibration and strip-deck campaign | ~1 week |

**Total about 5–6 engineer-weeks.** Adversarial gate required (new maths).

**Calibration plan from TIMs' existing DM04 data [I]:**
1. **Transfer directly:** the CSL (e0 0.83, λc 0.027, ξ 0.45); Mc = M_tc 1.3309; G from G0 with the √p
   exponent; ν; the initial e.
2. **Refit from the same drained triaxials (dense and loose):**
   - χ_tc from the peak-dilatancy-versus-ψ plot (D_min = χψ). This is the analogue of DM04's nd/A0.
   - H0 and H_y from the pre-peak stiffness and the strain to peak.
   - N and N̄ from the volumetric curves.
3. **Undrained tests** as a check on ψ-dependence only. They are not needed for the drained strip.
4. **Cross-model gate:** single-point drained and undrained triaxial and plane-strain compression, NorSand
   against DM04 (the oracle-integrated `uw_model`). The peak q and the ε_v curves must agree within the
   triaxial scatter.
5. **BVP gate:** the strip deck at B/8 and B/16 with NorSand and DM04 (SAS-ME), plus PM4Sand as a practice
   check. Compare the peak, the post-peak and the mechanism. Also compare against Loukidis–Salgado-type N_γ
   at the deck's density and stress level **[E, abstract only; the numbers need the paper]**.

Keep **HySand** (or the 2026 third-potential hyperplastic model) on the list for the cyclic SSI phase, when
its code and integration papers become readable.

---

## 8. Established facts vs inference (summary)

**Established [E]:**
- Borja & Andrade's NorSand variant: hyperelastic energy, D ≥ 0 under N̄ ≤ N, implicit return map and
  algorithmic tangent, the free-energy coupling approximation, the Q-cap.
- HySand: formulation, parameters, UMATs, and the 10 kPa surcharge.
- Variational Cam-clay is associative.
- Zhang 2026: SANISAND-04 implicit, strip footing.
- Roy 2020: DM04 footing challenges.
- Gao 2021: nonlocal fix for mesh dependence.
- Irani 2024 and Lan 2023: hyperelastic DM-type.
- The fork's available models.
- WP-134's integration findings.

**Inference [I]:**
- DM04's plastic part has no potential and positivity is not guaranteed.
- Non-symmetric tangents in non-associated hyperplasticity.
- The ln-CSL unboundedness of ψ at p → 0.
- The floor-sensitivity estimate.
- The effort estimates.
- The ranking itself.

---

## Sources

**Local (fork, read-only, `origin/ladruno`):**
- `Ladruno_implementation/_tims_2d_model_requests_2026-09-25.md`, `127_tims_2d_requests_plan.md`,
  `128_sanisand_ring_trace.md` (incl. its correction note), `134_sanisand_reference_integrator.md`,
  `133_pdmy_notes.md`
- scratchpad `sanisand_external_survey.md`
- the SRC tree listing (`UWmaterials/`, `UANDESmaterials/SAniSandMS*`, `ASDPlasticMaterial3D/`, `PM4Sand.h`)

**Read (E):**
- Borja & Andrade (2006). Critical state plasticity Part VI. CMAME 195:5115–5140. Author PDF:
  https://geomechanics.civil.northwestern.edu/Papers_files/ccvi.pdf
- Simonin, Houlsby & Byrne (2026). The HySand hyperplasticity constitutive model for sand: theory.
  Géotechnique 76(13). https://doi.org/10.1680/jgeot.25.00091
- Saberi, Simonin, Houlsby & Byrne (2025). ISFOG 2025.
  https://www.issmge.org/uploads/publications/132/133/ISFOG2025-301.pdf
- Byrne et al. (2025). PICASO. ISFOG 2025. https://www.issmge.org/uploads/publications/132/133/ISFOG2025-506.pdf
- Simonin, Houlsby & Byrne (2023). SEG23. https://proceedings.open.tudelft.nl/seg23/article/view/631
- Ortiz & Pandolfi (2004). A variational Cam-clay theory of plasticity. CMAME 193:2645.
  https://doi.org/10.1016/j.cma.2003.08.008 (via https://authors.library.caltech.edu/records/zfc70-0zg50)
- Houlsby (2019). Frictional plasticity in a convex analytical setting. Open Geomechanics.
  https://opengeomechanics.centre-mersenne.org/articles/OGEO_2019__1__A3_0/ (abstract)
- Zhang, Gao, Lu, Zhou, Lai & Du (2026). IJNAG. https://doi.org/10.1002/nag.70354 (abstract)
- Roy, Chow, O'Loughlin, Randolph & Whyte (2020). CGJ. https://doi.org/10.1139/cgj-2019-0841 (abstract)
- Loukidis & Salgado (2011). Géotechnique 61(2):107. https://doi.org/10.1680/geot.8.P.150.3771 (abstract)
- Gao, Li & Lu (2021/22). Acta Geotech. 17:427. https://doi.org/10.1007/s11440-021-01236-3 (abstract)
- Gao, Lu & Du (2020). JEM 146(8). https://doi.org/10.1061/(ASCE)EM.1943-7889.0001814 (metadata)
- Chaloulos, Papadimitriou & Dafalias (2019). JGGE 145(10). https://doi.org/10.1061/(ASCE)GT.1943-5606.0002082
  (metadata)
- Irani et al. (2024). IJNAG 49(1). https://doi.org/10.1002/nag.3852 (abstract)
- Lan, Liu & Zhao (2023). C&G 159. https://www.sciencedirect.com/science/article/pii/S0266352X23001854
  (abstract via ADS/search)
- Toll (2011). The dissipation inequality in hypoplasticity. Acta Mech. 221:39.
  https://doi.org/10.1007/s00707-011-0487-x (metadata)
- Houlsby & Richards (2023). C&G 156:105143. https://doi.org/10.1016/j.compgeo.2022.105143 (metadata)
- Orakci, Anoyatis, Chow & François (2026). engrXiv preprint. https://engrxiv.org/preprint/view/6998/version/9080
- Itasca NorSand documentation. https://docs.itascacg.com/itasca940/common/models/norsand/doc/modelnorsand.html
- Woudstra (2021). TU Delft MSc. https://repository.tudelft.nl/record/uuid:dc29fd0a-6e8a-4f94-92d1-cd7b9a66c4fa
- INL GranularFlowModels. https://github.com/idaholab/GranularFlowModels
- Coombs & Crouch (2011). CMAME 200:2297. https://www.sciencedirect.com/science/article/abs/pii/S0045782511001320
  (abstract)
- Granular Matter (2025). A thermodynamic critical state model for sands. https://doi.org/10.1007/s10035-024-01492-6
  (abstract via search)
- Cheng & Jefferies (2020). Geo-Congress. https://ascelibrary.org/doi/10.1061/9780784482810.002 (abstract via
  search)

**Secondary / recollection (E-sec / R):**
- Collins & Houlsby (1997). Proc. R. Soc. A 453:1975.
- Collins (2005). Géotechnique 55(5):373. https://doi.org/10.1680/geot.2005.55.5.373
- Collins & Muhunthan (2003). Géotechnique 53:611.
- Houlsby, Amorosi & Rojas (2005). Géotechnique 55(5):383.
- Einav & Puzrin (2004), CHCS, IJSS. https://www.sciencedirect.com/science/article/abs/pii/S0020768303005134
- Einav, Puzrin & Houlsby (2003). IJNAG. https://onlinelibrary.wiley.com/doi/abs/10.1002/nag.303
- Golchin & Lashkari (2014). IJSS. https://www.sciencedirect.com/science/article/pii/S0020768314001383
- Gajo & Wood (1999). Géotechnique 49(5):595 and IJNAG 23:925.
- Tejchman & Herle (1999). Soils Found. https://www.sciencedirect.com/science/article/abs/pii/S0038080620311318
- Dafalias & Manzari (2004). JEM 130(6):622.
- Taiebat & Dafalias (2008). IJNAG 32.
- Li & Dafalias (2000). Géotechnique 50(4).
- Li & Dafalias (2012). JEM 138(3).
- Jefferies (1993). Géotechnique 43(1).
- Jefferies (1997). Géotechnique 47(5):1037.
- Bolton (1986). Géotechnique 36(1).
- Yu (1998). CASM, IJNAG.
- Einav (2007). JMPS, breakage mechanics.
- Houlsby & Puzrin (2006). *Principles of Hyperplasticity*, Springer.
- A thermodynamically consistent variational neural update framework for geomaterials with non-associated
  flow (C&G 2026). https://www.sciencedirect.com/science/article/pii/S0266352X26006725 (403)

---

## Could not verify

- The **full text of any paywalled paper** (HTTP 403). This includes Loukidis & Salgado's N_γ numbers, whether
  they floored low stress, and their mesh study; Zhang et al. 2026's strip-footing details (p-floor,
  plateau); and Roy et al. 2020's footing results.
- **HySand's integration scheme and low-stress handling** beyond the ISFOG surcharge remark (Simonin 2023
  DPhil; Houlsby 2025a), and whether any HySand code is public.
- The **"shift stress" terminology** and its exact definition in Collins & Houlsby 1997. The mechanism
  (stress-dependent dissipation leading to non-associativity) is E-sec. The term is R.
- Whether **Jefferies' NorSand** (Itasca/PLAXIS) uses a hyperelastic energy. The Itasca documentation gives a
  power-law G only.
- Whether **3-invariant NorSand** keeps B&A's dissipation proof. It is 2-invariant in the source; re-derivation
  is needed.
- **Houlsby, Amorosi & Rojas' convexity limit** for n = 0.5 energies at high stress ratio. This matters because
  the ring reaches η ≈ M^b ≈ 2.1.
- **Tejchman & Herle's agreement** with the model tests, and hypoplasticity's small-stress fix (only
  summaries and recollection).
- The **2026 variational neural update** for non-associated geomaterials. It could change the "no sand model
  is variational" statement; I read only its title.
- **PM4Sand's monotonic drained dense-sand performance** at footing scale. No study found.
- No **thermodynamic reformulation of DM04/SANISAND** was found. Absence of evidence, not evidence of absence.
