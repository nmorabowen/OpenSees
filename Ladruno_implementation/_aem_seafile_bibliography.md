# AEM — Seafile local bibliography (primary sources)

Staging catalog of Applied Element Method (AEM) references found in the personal
Seafile library, to be merged into the deep-research report on *AEM vs OpenSees*.
These are **primary sources not on the open web** (proprietary manuals, JSCE/JNDS
papers, the FDEM book) and should anchor the theory + feasibility sections.

Library root: `C:\nmb\My Libraries\Libros\Estructuras\`

## Foundational AEM papers (verified citations — read directly)

1. **Meguro, K. & Tagel-Din, H. (2000).** "Applied Element Method for Structural
   Analysis: Theory and Application for Linear Materials." *Structural Eng./
   Earthquake Eng., JSCE*, Vol. 17, No. 1, 21s–35s (J. Struct. Mech. Earthquake
   Eng., JSCE, No. 647/I-51), April 2000. — **The foundational AEM paper** (~174
   citations). Springs in normal + tangential directions; Poisson's-ratio effect;
   RC + confinement; small-deformation linear-elastic formulation.
   File: `AEM/Applied_element_method_for_structural_analysis_The.pdf`
   - Spring stiffness: `Kn = E·d·T/a`, `Ks = G·d·T/a` (d = spring spacing,
     T = thickness, a = representative-area length). Each spring = stiffness of a
     volume (d × T × a). Rebar stiffness added to matrix stiffness at the spring.
   - 2D element = rigid body, **3 DOF**; element shape fixed, all deformation lives
     in the springs; the *assembly* is deformable. Stiffness matrix (Eq. 2) built
     by unit-displacement per DOF, summing surrounding-spring contributions.

2. **Meguro, K. & Tagel-Din, H. S. (2002).** "Applied Element Method Used for
   Large Displacement Structural Analysis." *Journal of Natural Disaster Science*,
   Vol. 24, No. 1, pp. 25–34. — Large-displacement / buckling / post-buckling;
   **no geometric stiffness matrix needed**.
   File: `AEM/Applied Element Method Used for Large Displacement Structural Analysis.pdf`
   - "Areas of application" matrix (I–VI) maps the AEM development series:
     I small-def elastic monotonic (Tagel-Din & Meguro 1998); II small-def
     nonlinear monotonic (Meguro & Tagel-Din 1998); III cyclic (Meguro &
     Tagel-Din 1997); IV large-def elastic [this paper]; V dynamic monotonic
     (Tagel-Din & Meguro 1999); VI dynamic cyclic (Tagel-Din & Meguro 2000).

## ELS Theoretical Manual V9 (authoritative, proprietary — read TOC + intro)

**Applied Science International, LLC (2004–2022).** *Extreme Loading® for
Structures — Theoretical Manual, Version 9.* Durham, NC.
File: `AEM/ELS/Theoretical Manual V9.pdf` (+ Blast Manual V9, User-defined Material
Manual V9, ELS-OpenFOAM User Guide V9, Viewer Manual V9, Quick Start V9).

The single most useful source for the **feasibility** angle — it lays out exactly
the machinery OpenSees lacks. Structure of the formulation (TOC, p.iv–v):
- **§2 Element overview**: 8-node hexahedron solids, link members, dash-pot
  dampers; connectivity (FEM vs AEM); DOF + rigid-body option; **connectivity
  springs** (matrix springs, reinforcement springs, stiffness §2.5.3); thin-
  structure shear-locking fix; equilibrium equations; effective modal mass.
- **§3 Material models & failure criteria**: *separation strain* (§3.1.8),
  *friction coefficient* (§3.1.9), residual shear strength, rupture of reinforcing
  bars, failure-softening factor; concrete + reinforcing-bar models; masonry;
  **Delaunay/Voronoi meshing** (§Fig 3-25). Fig 3-9 "elements separate and
  re-contact again" is the separation↔recontact cycle.
- **§4 Loading**: load/displacement control; dynamic time step; earthquake; blast
  scenario; pressure control; **element demolish & construction** (§4.4.7 —
  runtime topology change); temperature.
- **§6 Element Contact** — the broad-phase/collision core OpenSees has no analogue
  for: contact types **Corner-Face, Edge-Edge, Corner-Ground** (§6.2); contact
  stiffness (§6.3); **energy dissipation during contact** (§6.4); per-contact-type
  time step (§6.5). Unloading-stiffness factor governs rebound (Fig 6-7/6-8).
- Narrative (p.2–4): collapse passes through a **continuum stage → discrete
  stage**; Fig 1-1 "Collapse History": Linear → Crack/Yield/Crush → Buckling/
  Post-buckling (continuum) → **Element Separation → Debris as Rigid Bodies →
  Collision** (discrete). FEM flagged "not automated" at separation and "time
  consuming" at collision; AEM spans the whole range.
- Embedded validation reference chain (>50 experimental/theoretical comparisons):
  Tagel-Din & Meguro 2000, Meguro & Tagel-Din 2001, Tagel-Din 2002, Meguro &
  Tagel-Din 2003, Tagel-Din & Rahman 2004, Sasani & Sagiroglu 2008, Sasani 2008,
  Park et al. 2009, Wibowo et al. 2009, Helmy et al. 2009, Galal & El-Sawy 2010,
  Khalil 2011, Salem et al. 2011, Helmy et al. 2012/2013, Salem & Helmy 2014,
  Salem et al. 2014/2016, Attia et al. 2017, Zerin et al. 2017, Alanani et al.
  2020, El-desoqi et al. 2020, Attia & Salem 2021.

## Other AEM-folder papers (not yet read in detail — titles/filenames)

- `AEM-ULTIMATE ANALYSIS OF PROGRESSIVE COLLAPSE.pdf`
- `AEM simulation of collapse.pdf`
- `Seismic Pregressive Collapse of Reinforced Concrete Frames AEM.pdf`
- `Simulation_of_Concrete-Frame_Collapse_due_to_Dynam.pdf`
- `Analysis_of_small_scale_RC_building_subjected_to_s.pdf`
- `VORONOI APPLIED ELEMENT METHOD FOR STRUCTURAL ANALYSIS.pdf` (+ `_ANAL` variant)
- `AEM Masonry Retrofit.pdf`, `aem masonry.pdf`,
  `SIMULATION_OF_BRICK_MASONRY_WALL_BEHAVIOR_UNDER_IN.pdf`
- `AEM ground surface deformation.pdf`, `ASSESSMENT_OF_SEISMIC_DAMAGE_TO_RAILWAY_STRUCTURES.pdf`
- `Lincy_Christy_2018_IOP_Conf._Ser._Mater._Sci._Eng._330_012117.pdf` (review)
- `B2006-ERS-Meguro.pdf`, `SF-Applied-Element-Method-Tagel-Rahman.pdf`,
  `PP-band-State-of-Art-Report.pdf` (masonry retrofit, ICUS)
- ELS sample set (`AEM/ELS/Samples/*.ELS + .pdf`) — ~70 worked verification models.
- `AEM/AEM APPLICATION (YOU TUBE)/*.flv` — demo videos (blast, demolition,
  progressive collapse, impact, seismic).

## Related-family + comparison references (outside the AEM folder)

- **FDEM (Munjiza):** `FEM/Teoria/The Combined Finite-Discrete Element Method -
  Ante Munjiza.pdf` — the FDEM reference book; the rigorous combined finite-
  discrete cousin of AEM (contact detection, contact interaction, fracture).
- **Contact mechanics:** `Civilax/Analysis & Design Books/Contact Problems in
  Elasticity.pdf` (Kikuchi & Oden); `…/Finite Element Analysis in Beam to Beam
  Contact.pdf`; `…/Structural Impact.pdf` (Jones).
- **OpenSees + progressive collapse:** `Matrix/Erratum-Krylov-subspace-
  accelerated-newton-algorithm-Application-to-dynamic-progressive-collapse-
  simulation-of-frames-JSE.pdf` — KSA-Newton applied to frame progressive-collapse
  dynamics (directly relevant to the OpenSees feasibility angle).
- **Progressive collapse (general):** `Simica/Best Practices for Reducing the
  Potential for Progressive Collapse in Buildings.pdf` (NIST); `Tesis/Sezen -
  Progressive Collapse.pdf`; `Articulos/Progressive Failure.pdf`.
- **Blast loading:** `Blast/[Needham] Blast Waves.pdf`; `AEM/ELS/Blast Manual V9.pdf`.

> Note: the PBD-folder "collapse" papers (FEMA P695 / ATC-63, IDA, fragility) are
> *statistical sideways-collapse risk*, NOT physical element-separation collapse —
> different sense of "collapse"; keep distinct in the report.
