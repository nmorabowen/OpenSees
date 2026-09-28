# R1 - Literature survey: the h singularity at load reversal in SANISAND-type bounding-surface models

> **Post-hoc verification (main session, 2026-09-28).**
> - Chen et al. (2022), §3.9.1 (the SANISAND04 footing that aborts as (α − α_in):n drops to 0, worse with
>   finer steps or a tighter STOL) is re-read in Chen (2023), Monash PhD thesis, doi:10.26180/23639730.v1,
>   pp. 3-37–3-39. That thesis declares its Chapter 3 to be the published C&G 2022 paper, with sections not
>   renumbered. The journal PDF itself was not accessed.
> - The quotations the WP-151 memo §3 takes from this survey were re-checked against the saved primary texts.
>   PM4Sand v3.3 p. 25: "to avoid division by zero" and the numerical-stability sentence are verbatim.
>   Ghorbani et al. (2023) says "spurious oscillations" (without "numerical").
>   Jeremić et al. (2008) says "the denominator of Equation (32a) becomes negative".
>   The PM4Sand init-guard sentence below is the agent's paraphrase, not a quotation.


Survey date: 2026-09-28. Read-only; no repository file was modified.

Tags:
- [E] - I opened the source in this session.
- [E, SA] - a parallel sub-agent opened it in this session; I did not re-open it.
- [E-abstract] - only an abstract or metadata record was read.
- [E-sec] - confirmed only through a named secondary source.
- [S] - search snippet only.
- [E-absence] - "not found in any source opened in this session" (the three sub-agents checked independently).
- [R] - my recollection, not verified.
- [I] - my own inference.

Notation: x = (alpha - alpha_in):n. b:n = (alpha^b_theta - alpha):n. H = Kp + 2G(...) - K D n:r (the loading-index denominator).

## Executive summary (what the literature offers, ranked by how established it is)
1. **Additive regularization of the reversal distance.** This is the most established approach in production codes.
   - PM4Sand (since v2, 2012) and PM4Silt, in FLAC and OpenSees: Kp = G h0 sqrt(b:n) / (exp(x_app) - 1 + C_gamma1), with C_gamma1 = h0/200. The manuals say the constant is there to avoid a division by zero.
   - Pisano's SANISAND-MS PLAXIS UDSM: h = b0 / (|x| + 0.001), citing PM4Sand.
   - Result: a bounded Kp at every re-seat. [E]
2. **Remove alpha_in from h (memoryless, bounded h).**
   - Taiebat & Dafalias (2008): h = b0 / [(3/2)((b_ref - b):n)^2], because alpha_in is hard to implement implicitly. [E]
   - MD97-type b_ref forms. [E-sec]
   - Chen, Ghorbani, Zhang & Kodikara (2022, C&G 152:105008): a softplus ratio, h = (b0/theta) ln(1 + e^(b:n)) / ln(1 + e^((alpha - alpha^b_theta+pi):n)), plus tanh smoothing of Kp across reversals. This paper also shows a SANISAND04 **footing** BVP aborting because x drops suddenly to 0, getting worse with smaller steps or tighter tolerance, and fixes it. It is the single most relevant source found. [E]
3. **Hysteretic / threshold re-seat and memory.** Established for overshooting.
   - Dafalias (1986) and SANISAND-Z (2016): plastic-strain-threshold weighting of alpha_in. [E-sec]
   - PM4Sand: alpha_in^app, alpha_in^p and C_rev. [E]
   - SANISAND-F: unloading alone does not re-seat. [E]
   - LiPa: formal/informal reversals. [E-sec]
   - Ghorbani et al. (2023, Comput. Mech.): a positive scalar floor x = J^r m_q, designed for spurious reversals from numerical oscillations; it reduced iterations and CPU in FE. [E]
   - P2PSand: `ratio-reverse` = 0.02. [E]
   - Caveat: Chen et al. show threshold schemes can still fall back to h = infinity and degrade as steps shrink. [E]
4. **Kp sign handling.**
   - PM4Sand v3.3 sets Kp = 0 outside the bounding surface and reports that this improved numerical stability. v2/v3 had allowed a negative Kp. [E]
   - The SANISAND family intends negative Kp (softening), and SANISAND-MSf calls alpha outside the bounding surface standard DM04 behaviour. [E]
   - The OpenSees DM04 ports leave Kp unclamped and turn dGamma < 0 into an elastic step. [E]
5. **Implementation-only caps, undocumented in any paper.** [E]
   - OpenSees UW ManzariDafalias: h = 1e10 only if |x| < 1e-10, otherwise a signed h.
   - SAniSandMS: x >= 1e-10 and h <= 1e7.
   - Legacy UCD DM04: x >= 1e-10, re-seating inside iterations.
6. **Conceptual alternatives.** [E]
   - Hashiguchi's subloading surface: a continuous normal-yield ratio, no memory reset, and a strain-space loading criterion.
   - SANISAND-MSf's cancellation of a singular hM denominator by |.| and Macaulay factors.

**Gap.** No source treats the DM04 re-seat (h = infinity) coinciding with b:n <= 0, i.e. the infinity * (negative) or 0/0 limit. [E-absence]
- The nearest analyses are Chen et al. (2022), on reversal during softening for MD97 (h < 0, Kp << 0, spurious stiffening), and PM4Sand, which has C_gamma1 and the Kp = 0 floor as separate provisions.
- [I] The literature-backed combination for our model is a bounded distance (additive epsilon, positive floor, or memoryless softplus), plus a Kp floor or cap (Kp >= 0, or at least H > 0), plus reversal detection on the committed state with hysteresis.

## Q1. DM04 itself

Citation: Dafalias, Y.F. & Manzari, M.T. (2004). "Simple plasticity sand model accounting for fabric change effects." *J. Eng. Mech.* 130(6):622-634. doi:10.1061/(ASCE)0733-9399(2004)130:6(622).
Access: the ASCE full text is paywalled. ResearchGate returned 403, the Semantic Scholar page came back empty, and IA Scholar showed an anti-scrape challenge (not attempted). **I did not read DM04 itself.** Everything below is [E-sec] from papers that quote or reproduce it (several co-authored by Dafalias) or [R]. **DM04's own equation numbers are unverified.**

- Multiaxial hardening, as reproduced by several independent secondary sources:
  - h = b0 / [(alpha - alpha_in):n], with b0 = G0 h0 (1 - ch e)(p/p_at)^(-1/2)
  - Kp = (2/3) p h (alpha^b_theta - alpha):n
  - d(alpha) = <L> (2/3) h (alpha^b_theta - alpha)
  - [E-sec] Sources: Jeremic, Cheng, Taiebat & Dafalias (2008), IJNAG 32(13):1635-1660 [vol./pages R], eqs. (31)-(32b) (sokocalo.engr.ucdavis.edu/~jeremic/wwwpublications/CV-J20.pdf); Duque et al. (2022), Acta Geotech. 17:2235-2257, Table 5; the OpenSees ManzariDafalias doc page (opensees.github.io/OpenSeesDocumentation/.../ManzariDafalias.html); PM4Sand v3.3 eqs. (33)-(36).
- alpha_in rule, as reproduced by Jeremic et al. (2008): alpha_in is the value of alpha at the initiation of a new loading process, and it is updated to the current alpha when the denominator (alpha - alpha_in):n becomes negative. [E-sec] PM4Sand v3.3 Sec. 2.5 describes the same rule as the Dafalias (1986) bounding-surface practice that DM04 adopted. [E-sec] SANISAND-MSf (Yang, Taiebat & Dafalias 2022, Table 1 and text) repeats it: alpha_in is updated when the denominator of h becomes negative, following the rules of Dafalias (1986). [E-sec, co-authored by Dafalias]
- h = infinity at the initiation of loading (by design):
  - Taiebat & Dafalias (2008, IJNAG 32:915-948, p. 926) describe DM04 as using |alpha - alpha_in| in the denominator of h (triaxial form), which gives an infinite h at initiation, hence an infinite Kp and a zero loading index L. [E-sec, written by Dafalias]
  - Taiebat & Dafalias (2015, 6ICEGE, SANISAND-Z preview) say the same for the r_in form. It gives Kp -> infinity at initiation, which yields a smooth elastic-plastic transition, but it also creates "the so-called overshooting response upon reverse loading/immediate reloading", known since Dafalias (1975). [E - for that paper's statement]
  - [R] DM04's triaxial form is h = b0/|alpha - alpha_in|, with alpha_in updated at each load increment reversal.
- Overshooting and small reversals in DM04 itself: not verified. Dafalias' later papers (ICEGE 2015; SANISAND-Z 2016, see Q2) introduce the Dafalias (1986) overshooting remedy as something new for the SANISAND family. [E] [I] That suggests DM04 used the plain re-seat without a remedy.
- b:n < 0 (alpha outside the bounding image): no DM04 statement verified. Within the family it is explicitly allowed as softening:
  - Taiebat & Dafalias (2008, p. 926) note that alpha^b varies with the state parameter psi, so alpha can momentarily "cross" a decreasing alpha^b during monotonic loading; alpha^b - alpha < 0 then gives d(alpha) < 0 (softening). [E]
  - Yang, Taiebat & Dafalias (2022, SANISAND-MSf, p. 232) state that during dilation and softening alpha lies outside the bounding surface, and call this a standard feature of DM04. [E] [I] So the collision in our footing - a re-seat (h = infinity) at a moment when (alpha^b_theta - alpha):n <= 0 (alpha outside the Lode-dependent image in the NEW direction n) - joins two behaviours the model intends (infinite h at initiation; negative Kp outside the image). Their product (+infinity) x (negative), or 0/0 as both terms -> 0, is not addressed in any source I opened.
- Loading-index denominator H = Kp + 2G(B - C tr(n^3)) - K D n:r (OpenSees UW form; PM4Sand eq. 31 gives the B = 1, C = 0 version). A strain-driven step has a plastic solution only if H > 0. [E-sec for the equation; I for the well-posedness remark]

## Q2. SANISAND successors

Notation: x = (alpha - alpha_in):n and b:n = (alpha^b_theta - alpha):n.

### 2.1 Dafalias, Papadimitriou & Li (2004) - DPL04
"Sand plasticity model accounting for inherent fabric anisotropy." *J. Eng. Mech.* 130(11):1319-1333. doi:10.1061/(ASCE)0733-9399(2004)130:11(1319).
- Full text: paywalled (ASCE 403). Abstract read via OpenAlex. [E-abstract, SA]
- The model adds A = F:n dependence of the plastic modulus and of the CSL location to a pre-existing stress-ratio driven bounding-surface model. [E-abstract, SA]
- h form, alpha_in rule, regularization: not verified. My recollection is a DM04-type h = b0/x with an A-dependent factor. [R]

### 2.2 Taiebat & Dafalias (2008) - SANISAND (closed yield surface)
"SANISAND: Simple anisotropic sand plasticity model." IJNAG 32(8):915-948. doi:10.1002/nag.651. Read in full: https://escholarship.org/content/qt06r6c3cw/qt06r6c3cw.pdf [E]
- h form:
  - Triaxial, eq. (17a): h = b0 / (b_ref - s(alpha^b - alpha))^2, with b_ref = alpha^b_c + alpha^b_e, the diameter of the bounding surface.
  - Eq. (17b): b0 = G0 h0 (1 - ch e)(p_at/p)^(1/2).
  - Multiaxial, eq. (35): h = b0 / [ (3/2) ((b_ref - (a^b - a)) : n)^2 ], with b_ref = sqrt(2/3) b_ref n.
  - [E]
- Regularization: alpha_in is ELIMINATED. The authors say the DM04 |alpha - alpha_in| gives an infinite h at initiation, hence an infinite Kp and a zero L, which would block p0 changes under constant-eta loading. The new form also simplifies implementation, because it removes an updatable discrete memory parameter that is difficult to implement implicitly (p. 926). [E]
  - The squared fixed-reference denominator is never negative. It vanishes only when alpha sits at the opposite bounding-surface image, e.g. a full reversal from the bounding surface. [I]
  - For b:n < 0 the denominator exceeds b_ref^2, so Kp stays finite and negative (bounded softening). [I]
- Softening is explicitly allowed: alpha^b - alpha < 0, from alpha crossing a shrinking alpha^b, gives a negative rate. [E]
- |eta - alpha| is needed to avoid a wrong sign (pp. 926-927). [E]
- Reversal rule: none needed. [E]

### 2.3 Li & Dafalias (2012) ACST, and Gao, Zhao, Li & Dafalias (2014)
- Li, X.S. & Dafalias, Y.F. (2012). "Anisotropic critical state theory: role of fabric." *J. Eng. Mech.* 138(3):263-275. doi:10.1061/(ASCE)EM.1943-7889.0000324. Paywalled. Abstract via OpenAlex. [E-abstract, SA]
  - ACST is a theory (F, A, dilatancy state line). Its illustrations use a triaxial model. SANISAND-F describes that model as having no kinematic hardening. [E-sec, SA]
  - No h/alpha_in to report; not verified beyond that.
- Gao, Z., Zhao, J., Li, X.S. & Dafalias, Y.F. (2014). "A critical state sand plasticity model accounting for fabric evolution." IJNAG 38(4):370-390. doi:10.1002/nag.2211. The journal version is paywalled. The companion conference paper was read (https://eprints.gla.ac.uk/112972/1/112972.pdf). [E, SA]
  - Isotropic-hardening, monotonic model.
  - No alpha_in and no reversal rule.

### 2.4 Dafalias & Taiebat (2016) - SANISAND-Z (zero elastic range)
"SANISAND-Z: zero elastic range sand plasticity model." *Geotechnique* 66(12):999-1013. doi:10.1680/jgeot.15.P.271.
- The journal is paywalled. The abstract was read via Crossref. [E-abstract, SA]
- I read the conference precursor: Taiebat, M. & Dafalias, Y.F. (2015). "A true zero elastic range sand plasticity model." 6ICEGE, Christchurch. https://www.issmge.org/uploads/publications/59/60/196.00_Taiebat.pdf [E]
- Form and loading definition:
  - Eq. (6): Kp = (2/3) p h [(r^b - r):n] / [(r - r_in):n], with h = G0 h0 (1 - ch e)(p/p_at)^(-0.5), i.e. DM04's b0.
  - The image point is defined by a stress-rate-dependent mapping (eq. 4), so L never becomes negative. [E]
  - A new plastic loading event is declared when (r - r_in):n <= 0, and r_in is re-seated so loading restarts with infinite Kp and no purely elastic phase. [E]
- Overshooting remedy: the paper applies the one from Dafalias (1986). A threshold on previous cumulative plastic strain weights the r_in update between its past and current values. [E] The journal abstract confirms an r_in updating scheme that avoids overshooting. [E-abstract, SA]
- The exact formula is reproduced by Ghorbani et al. (2023, Comput. Mech.), eqs. (11)-(12) [E-sec]:
  - alpha_in^(i+1) = alpha^(i) + m_q (alpha_in^r - alpha^(i))
  - m_q = <1 - (eps_q^p / eps_bar_q^p)^j>, with j = 1 by default and eps_bar = 0.01%, exponent placement per Chen et al. 2022 eq. (3-21) (see Q6.1 and Q6.1b)
- Denominator regularization: none. [E - conference version]

### 2.5 Petalas, Dafalias & Papadimitriou - SANISAND-F (2020) and SANISAND-FN (2019)
- SANISAND-F: "SANISAND-F: Sand constitutive model with evolving fabric anisotropy." *Int. J. Solids Struct.* 188-189:12-31 (2020). doi:10.1016/j.ijsolstr.2019.09.005. Read in full: https://research.chalmers.se/publication/512957/file/512957_Fulltext.pdf [E]
  - d(alpha) = <L> H (alpha^b_theta - alpha) (eq. 16).
  - **H = (2/3) h(e,p,A) / <(alpha - alpha_in):n>** (eq. 17) - a Macaulay bracket IN THE DENOMINATOR.
  - h = G0 h1 exp(h2 A)(e^-1 - ch)^2 (p/p_at)^(-1/2) (eq. 18).
  - Kp = p H (alpha^b_theta - alpha):n (eq. 19).
  - [E]
- The re-seat rule (pp. 16-17) [E]:
  - When <.> becomes zero or negative, n points against the recent loading direction alpha - alpha_in, and a new plastic loading process begins: alpha_in := alpha, and H restarts infinite.
  - This can happen through UNLOADING-reloading, but also as n rotates progressively during continued loading without L ever being < 0.
  - With continuous rotation of n, the update happens at the moment x = 0, "without any discontinuity" of H, which approaches infinity.
  - An elastic unloading (L < 0) by itself does NOT re-seat alpha_in. The update happens only if the new n gives x <= 0. This avoids an unjustified infinite Kp when reloading starts at, or very near, the unloading point.
  - [I] This is a reversal rule tied to the direction test only, not to the loading/unloading event, and it is intended to suppress spurious re-seats.
- Softening: Kp can become zero and then negative when alpha moves outside the bounding surface, because alpha^b depends on zeta (eq. 12a). [E]
- A re-seat at the moment b:n < 0 is not discussed. [E - absence, SA confirmed]
- SANISAND-FN: Petalas, Dafalias & Papadimitriou (2019). "SANISAND-FN: An evolving fabric-based sand model accounting for stress principal axes rotation." IJNAG 43(1):97-123. doi:10.1002/nag.2855. Abstract only [E-abstract, SA]: it modifies dilatancy and Kp through a non-coaxiality measure. That it inherits eq. (17) is [I].
- Petalas & Dafalias (2019), implicit integration of SANISAND-Z, *Comput. Geotech.* 112:386-402: not read.

### 2.6 SANISAND-MS (Liu, Abell, Diambra & Pisano) and variants
- Liu, H.Y., Abell, J.A., Diambra, A. & Pisano, F. (2019). "Modelling the cyclic ratcheting of sands through memory-enhanced bounding surface plasticity." *Geotechnique* 69(9):783-800. doi:10.1680/jgeot.17.P.307. Read the TU Delft accepted manuscript (copy saved by a sub-agent). [E]
  - Kp = (2/3) p h (r^b - r):n (eq. 8).
  - **h = b0/[(r - r_in):n] * exp[mu0 (p/p_atm)^0.5 (b^M/b_ref)^2]** (eq. 11).
  - b^M = (r^M - r):n; b_ref = (r^b - r^b_+):n, the opposite projection, which is always > 0 (eq. 12).
  - r_in is r at the onset of load reversal, updated to the current r whenever (r - r_in):n < 0 (footnote). [E]
  - No denominator regularization in the paper. [E]
- Liu (2020) TU Delft PhD thesis, Ch. 3 [E, SA]:
  - Overshooting in small unload-reload cycles is acknowledged and deliberately left untreated (p. 40, citing Kan & Taiebat 2014 and Dafalias & Taiebat 2016).
  - The r (stress-ratio) form can push part of the yield locus outside the bounding surface and cause artificial softening (p. 47). That is why the back-stress (alpha) form was adopted in OpenSees and later chapters.
- Liu & Pisano (2019). "Prediction of oedometer terminal densities through a memory-enhanced cyclic model for sand." *Geotechnique Letters* 9(2):81-88. Same h (GL eq. 2). [E, SA]
- SANISAND-MSu: Liu, Diambra, Abell & Pisano (2020). "Memory-enhanced plasticity modeling of sand behavior under undrained cyclic loading." JGGE 146(11):04020122. doi:10.1061/(ASCE)GT.1943-5606.0002362. Accepted manuscript read. [E, SA]
  - h = b0/x * exp[mu0 (p/p_atm)^0.5 (b^M/b_ref)^(w1) (1/eta)^(w2)] (eq. 16), with eta lower-bounded by m.
  - A Macaulay bracket on (alpha^b - alpha):n was added to the memory-surface shrinkage because that quantity becomes negative during softening (eq. 4).
  - The authors flag that the h^M denominator can vanish, "rare but possible", and leave it to the implementation (eq. 9 discussion).
- OpenSees SAniSandMS code guards: see Q5.2. [E]

### 2.7 Yang, Taiebat & Dafalias (2022) - SANISAND-MSf
"SANISAND-MSf: a sand plasticity model with memory surface and semifluidised state." *Geotechnique* 72(3):227-246. doi:10.1680/jgeot.19.P.363. Read in full: https://escholarship.org/uc/item/0227z2t1 [E]
- Authors: exactly three - Ming Yang, Mahdi Taiebat, Yannis F. Dafalias - per the cover page and running head. There is no "Mohammadi-Haji". [E]
- h:
  - Table 1: h = b0/x, with b0 = G0 h0 (1 - ch e)(p/p_at)^(-1/2). An alternative b0 multiplies this by a Lode-angle factor g(theta,c) raised to a constant n_g. The sign of that exponent is lost in the text extraction and was not verified.
  - **Eq. (11): h = [b0/x] * exp[ mu0/(||alpha_in||^u + eps) * (b^M/b_ref)^w ], defaults w = 2, eps = 0.01.** eps guards only the 1/||alpha_in||^u factor, NOT x. [E]
- alpha_in rule: updated "when the denominator of h becomes negative, as per the rules discussed by Dafalias (1986)" (p. 229). No overshooting rule is given in the paper. [E]
- The singular-denominator traps and their cancellation fix, eqs. (9)-(10): see Q6.3. [E]
- Semifluidised-state factor: h0 = h0'[(1 - <1 - p/p_th>)^(x_l) + f_l], f_l = 0.01 (Table 2). [E, SA]

### 2.8 Barrero, Taiebat & Dafalias (2020) - SANISAND-SF
"Modeling cyclic shearing of sands in the semifluidized state." IJNAG 44(3):371-388. doi:10.1002/nag.3007.
- Abstract only (Wiley PDF 403) [E-abstract, SA]: a strain-liquefaction factor at low p reduces plastic shear stiffness and dilatancy.
- As adopted in MSf: l_dot = <L>[c_l <1 - p_r> (1 - l)^(n_l)] - c_r l |eps_v_dot|. [E-sec]
- h/alpha_in: plain DM04 [R]. Ghorbani et al. (2023) list Barrero et al. among works preventing drastic Kp changes at spurious reversals [E-sec]; the mechanism is not verified.

### 2.9 Other variants and reversal-surface relatives
- SANISAND-FM: Zeng, Taiebat & Dafalias (2026), *Geotechnique* 76:787-804, doi:10.1680/jgeot.25.00540. ACST plus a memory surface. [E-abstract, SA]
- Reyes, Taiebat & Dafalias (2025), JGGE, doi:10.1061/JGGEFK.GTENG-12658. SANISAND-MSf for non-zero mean shear. [E-abstract, SA]
- Papadimitriou, Chaloulos & Dafalias (2019), *Acta Geotech.* 14:253-277. Reversal surfaces: the last stress reversal point is the projection centre. [E-abstract]
- Limnaiou & Papadimitriou (2023), *Acta Geotech.* 18:235-263. Reversal surfaces with an overshooting-safe reversal-point update. [E-abstract]; formal/informal rule per Q6.1 [E-sec]
- Chen, C. et al. (2024; online 2023), SANISAND-Z implicit integration with sub-stepping, *Comput. Geotech.* 165:105899: not read. [SA metadata only]

## Q3. PM4Sand / PM4Silt (Boulanger & Ziotopoulou)

Sources opened (full PDFs, text extracted and equation pages rendered and read):
- PM4Sand v2: Boulanger, R.W. & Ziotopoulou, K. (2012). *PM4Sand (Version 2): A sand plasticity model for earthquake engineering applications.* Report UCD/CGM-12/01, UC Davis (rev. 1). https://faculty.engineering.ucdavis.edu/boulanger/wp-content/uploads/sites/71/2014/09/Boulanger_Ziotopoulou_Sand_Model_CGM-12-01_2012_rev1.pdf
- PM4Sand v3: Boulanger & Ziotopoulou (2015). *PM4Sand (Version 3)...* Report UCD/CGM-15/01, March 2015. https://faculty.engineering.ucdavis.edu/boulanger/wp-content/uploads/sites/71/2014/09/Boulanger_Ziotopoulou_PM4Sand_Model_CGM-15-01_2015.pdf
- PM4Sand v3.3: Boulanger & Ziotopoulou (2023). *PM4Sand (Version 3.3)...* Report UCD/CGM-23/01. https://itasca-software.s3.amazonaws.com/udm-library/Boulanger_Ziotopoulou_PM4Sand_v3.3_CGM-23-01.pdf
- PM4Silt v2.1: Boulanger & Ziotopoulou (2023). *PM4Silt (Version 2.1)...* Report UCD/CGM-23/02 (rev. June 2023). https://itasca-software.s3.amazonaws.com/udm-library/Boulanger_Ziotopoulou_PM4Silt_Model_CGM-23-02-rev-2.pdf
- NOT read: v3.1 manual (UCD/CGM-17/01, 2017) - pm4sand.engr.ucdavis.edu is behind a Cloudflare challenge, the Boulanger page links a box.com share, Scribd shows only a preview. Ziotopoulou & Boulanger (2016) SDEE 84:269-283, doi:10.1016/j.soildyn.2016.02.013 - paywalled (Elsevier), not read.

### 3.1 Reversal detection and the alpha_in family (v3.3, Sec. 2.5, pp. 20-21)
- A reversal is identified, following traditional bounding-surface practice, whenever (alpha - alpha_in):n < 0. This is v3.3 eq. (19), identical to PM4Silt v2.1 eq. (22). [E]
- The manual states the problem plainly. Small load-reversal cycles reset alpha_in, raise Kp, and make the response overly stiff after a small reversal. It calls this a well-known bounding-surface problem with several partial remedies. [E] PM4Silt v2.1 cites Dafalias & Taiebat (2016) and Duque et al. (2021) for it. [E]
- On a reversal: alpha_in^p := alpha_in (previous initial), alpha_in := alpha (current). [E]
- alpha_in^app (apparent initial back-stress ratio), tracked per component over the whole history: for positive loading directions, the minimum value the component ever had, but not less than 0; for negative directions, the maximum ever had, but not more than 0. Stated rationale: avoid over-stiffening after small unload-reload cycles on an otherwise monotonic branch without tracking the history through many cycles. [E]
- Kp dependency, v3.3 eq. (20): if (alpha - alpha_in^p):n < 0 then Kp = f(alpha_in^true, alpha_in^app), else Kp = f(alpha_in^app). Kp is controlled first (inversely) by (alpha - alpha_in^true):n (large initial stiffness), then increasingly by alpha_in^app, and only by alpha_in^app once loading passes the previous reversal point. [E]
- History of the rule:
  - v2 (2012, Sec. 2.5) already constrained the re-seat. For positive loading, the new alpha_in = max(minimum value alpha_in has ever had, 0), mirrored for negative loading, so small reversals after much stronger loading cannot produce an overly stiff response. [E]
  - v3 (2015) lists changes to the initial back-stress-ratio tracking logic among its revisions. It introduced alpha_in^app with a hard switch, v3 eq. (20): if (alpha - alpha_in^p):n < 0, alpha_in = alpha_in^true, else alpha_in = alpha_in^app. [E]
  - v3.3 replaces the switch by the multiplicative factor C_rev (below). [E]
  - Ziotopoulou & Boulanger (2016) is the journal paper for the v3 changes. [R - not opened]
- Initialization guard (added in v3.2, documented in v3.3 Sec. 1.1 and Sec. 2.5): alpha_in at initialization is scaled so that its stress ratio is <= 0.9 Mb. Otherwise, for an initial state above the bounding surface, Kp = 0 (Mcur > Mb) and D = 0 (alpha - alpha_in = 0), so the stresses cannot change. [E] [I] This is the same kind of degenerate zero-distance pairing as our failure, fixed by keeping alpha_in off the current state.

### 3.2 Plastic modulus (v3.3, Sec. 2.7, pp. 24-25)
- DM04 restated (v3.3 eqs. 33-35): d(alpha) = <L>(2/3) h (alpha^b - alpha) (33); Kp = (2/3) p h (alpha^b - alpha):n (34); h = (3/2) Kp / [p (alpha^b - alpha):n] (35). Loading index, eq. (31): L = [2G n:de - n:r K d(eps_v)] / [Kp + 2G - K D n:r]. [E] (PM4Sand computes Kp first and then back-computes h from (35), so its h is singular at b:n = 0, not at reversal. [I])
- DM04 in PM4Sand notation, eq. (36): Kp = (2/3) G h0 [(1+e)/(2.97-e)^2 (1 - Ch e)] * [(alpha^b - alpha):n] / [(alpha - alpha_in):n]. [E]
- PM4Sand form, eq. (37): Kp = G h0 * [(alpha^b - alpha):n]^0.5 / { exp[(alpha - alpha_in^app):n] - 1 + C_gamma1 } * C_rev. [E]
- Eq. (38): C_rev = [(alpha - alpha_in^app):n] / [(alpha - alpha_in^true):n] for (alpha - alpha_in^p):n <= 0; C_rev = 1 otherwise. [E] PM4Silt v2.1 eqs. (40)-(41) are identical. [E]
- The C_gamma1 statement (v3.3 p. 25; the same text appears in v2 p. 23 and v3 p. 24): C_gamma1 prevents division by zero and slightly affects small-strain nonlinearity and damping. With C_gamma1 = 0, Kp is infinite at the start of each loading cycle because (alpha - alpha_in):n = 0, and nonlinearity appears only once that distance grows enough to bring Kp/G down to about 100-200. C_gamma1 = h0/200 was found to give a reasonable response. [E] [I] Near the reversal the denominator is exp(x) - 1 + C_gamma1 ~ x + C_gamma1, i.e. an additive regularization of the DM04 1/x term with epsilon = h0/200. So Kp_max ~ G h0 sqrt(b:n) / C_gamma1 = 200 G sqrt(b:n) when C_rev = 1.
- Sign of Kp outside the bounding surface - the rule changed between versions:
  - v2 (2012) eq. (37) and v3 (2015) eq. (38): for loose-of-critical states with (alpha^b - alpha):n < 0, the signs are modified to allow negative Kp: Kp = G h0 * ( -[-(alpha^b - alpha):n]^0.5 ) / { exp[(alpha - alpha_in):n] - 1 + C_gamma1 }. [E]
  - v3.3 (2023) p. 25: outside the bounding surface, Kp is set to zero rather than allowed to go negative. The manual's reason, quoted: "This restriction on the plastic modulus improved numerical stability" (the same sentence adds that it had little effect on computed responses). [E] Which version introduced this change could not be pinned down (the v3.1 manual was not accessible). The v3.3 revision notes say the v3.1 report added formulation and implementation clarifications prompted by user questions. [E]
  - PM4Silt v2.1 p. 22: the stress ratio is precluded from being outside the greater of the bounding and dilatancy surfaces in the implementation. [E]
- [I] As printed, eq. (38) does NOT bound Kp in the reload branch.
  - Right after a genuine reversal, x_true = 0 while x_app > 0, so C_rev -> infinity and Kp -> infinity (the intended stiff restart).
  - If alpha_in^app = alpha_in^true = alpha, C_rev = 0/0. In the OpenSees port this happens when the xy back-stress changes sign (see 3.4).
  - The OpenSees port evaluates C_rev as (<x_app> + C_gamma1)/(<x_true> + C_gamma1), which stays finite.
  - The FLAC DLL source was not seen.
  - A Kp = 0 floor for b:n < 0 must be applied before this multiplication to avoid 0 * infinity.

### 3.3 FLAC implementation details relevant to reversal robustness (v3.3 Sec. 3.1-3.2, Table 3.1, pp. 59-66)
- Explicit integration without sub-stepping. [E]
- In each zone, stresses and internal variables are averaged over the four sub-zones at the end of each step. The averaged state gets a drift correction (alpha projected along the zone-averaged stress-ratio direction). If the averaged stress ratio lies outside the bounding surface, it is projected back along the normal to the bounding surface. [E]
- D and Kp are recomputed from the averaged end-of-step state and used by all sub-zones in the NEXT step, so they lag one step. The manual notes other elasto-plastic models in FLAC use the same approach. [E]
- Reversal detection runs once per step on the averaged zone quantities: Table 3.1 step 11 computes (alpha - alpha_in):n for the averaged zone and decides whether a reversal happened, and step 17 then updates the initial back-stress ratios. [E]
- Version 1 gave each sub-zone its own D and Kp. It sometimes produced unusual deformation modes, attributed to sub-zones loaded in opposite directions (e.g. one strongly contractive after a reversal while others dilated). [E] [I] Detecting reversals on a committed, averaged state is a practical defence against iterate-level chatter.

### 3.4 OpenSees PM4Sand port (Chen & Arduino, UW) - source read locally [E]
File: SRC/material/nD/UWmaterials/PM4Sand.cpp. These are vanilla lines with no fork edits. I re-checked them against upstream master (Cg1 at line 2402, C_rev factor at 2426, Kp clip at 2465-2466; https://raw.githubusercontent.com/OpenSees/OpenSees/master/SRC/material/nD/UWmaterials/PM4Sand.cpp).
- GetStateDependent(): Cg1 = h0/200; x_app = Macaulay((alpha - alpha_in):n); x_true = Macaulay((alpha - alpha_in_true):n).
  - If |b:n| < 1e-10, h = 1e10 (the comment says this avoids a division by zero).
  - Else if (alpha - alpha_in_p):n <= 0: h = 1.5 G h0 / p / (exp(x_app) - 1 + Cg1) / sqrt(|b:n|) * Cka / (...) * (x_app + Cg1)/(x_true + Cg1).
  - Otherwise the same without the last factor.
  - Then Kp = (2/3) h p (b:n). [E]
  - [I] This means Kp = G h0 sign(b:n) sqrt(|b:n|) / (...), i.e. v3 eq. (38) negative Kp in the dilative branch.
- In the contraction branch ((alpha^dr - alpha):n > 0), Kp = max(0, Kp), with the comment "bound K_p to non-negative, following flac practice". [E]
- integrate(): reversal if (alpha - alpha_in_true):n_tr < 0, where n_tr is the yield normal at the elastic trial stress sigma_n + Ce:(eps - eps_n), evaluated against the committed alpha. Then alpha_in_p := alpha_in, alpha_in_true := alpha. alpha_in is set from the min/max history component-wise if alpha_xy * alpha_in_p,xy > 0 (same shear sign), else alpha_in := alpha. [E]

## Q4. Classic bounding-surface overshooting and memory rules

Most originals in this group are paywalled (ASME, Springer, ASCE, ICE, Elsevier). Much of this section is therefore [E-sec] or [R]; the secondary source is named each time.

### 4.1 Dafalias & Popov (1975, 1976)
- Dafalias, Y.F. & Popov, E.P. (1975). "A model of nonlinearly hardening materials for complex loading." *Acta Mech.* 21(3):173-192. doi:10.1007/BF01181053.
- Dafalias, Y.F. & Popov, E.P. (1976). "Plastic internal variables formalism of cyclic plasticity." *J. Appl. Mech.* 43(4):645-651. doi:10.1115/1.3423948.
- Metadata only [E-abstract, SA]; full texts not read.
- Plastic modulus: Kp = Kp_bar + h(delta) * delta / (delta_in - delta). delta is the distance from the stress point on the yield surface to its bounding image; delta_in is delta at the start of the current plastic loading process. Kp is infinite at re-yield by design, and delta_in is reset at each new plastic loading. [R] [I] The denominator delta_in - delta is the exact analogue of DM04's (alpha - alpha_in):n.
- Taiebat & Dafalias (2015) state that overshooting caused by the alpha_in update has been known since the method's inception in Dafalias (1975). [E]
- Metals remedy (Petersson & Popov 1977, *J. Eng. Mech. Div.* 103(4):611-627): interpolate the size of intermediate surfaces with a weight W that depends on the cumulative plastic strain up to the last reversal and the increment since. [E-sec, via Minagawa, Nishiwaki & Masuda 1987, JSCE Struct. Eng./Earthq. Eng. 4(2):361s-370s; SA]

### 4.2 Dafalias (1986)
"Bounding surface plasticity. I: Mathematical foundation and hypoplasticity." *J. Eng. Mech.* 112(9):966-987. doi:10.1061/(ASCE)0733-9399(1986)112:9(966). Paywalled, not read.
- The alpha_in memory variable and the reversal test (alpha - alpha_in):n < 0 are the traditional bounding-surface practice that DM04 adopted (PM4Sand v3.3, Sec. 2.5). [E-sec]
- Kp = infinity at initiation is intended. [E-sec, ICEGE 2015]
- Overshooting remedy [E-sec - Taiebat & Dafalias 2015; Ghorbani et al. 2023 eqs. 11-12; Chen et al. 2022 eqs. 3-20/3-21]:
  - A threshold on the cumulative plastic strain accumulated during the reversed segment weights the new alpha_in between its previous value and the current alpha.
  - Reconstructed: alpha_in_new = m alpha_in_old + (1 - m) alpha_current, with m = <1 - (eps_q^p / eps_bar_q^p)^j>.
  - SANISAND-Z used j = 1 and eps_bar = 0.01% (per Chen et al.).
- Dafalias (1986) also pointed out that the Fardis et al. (1983) d_min remedy makes Kp discontinuous. [E-sec, Chen et al. 2022]

### 4.3 Mroz, Norris & Zienkiewicz (1978, 1979)
- (1978) IJNAG 2(3):203-221, doi:10.1002/nag.1610020303; (1979) *Geotechnique* 29(1):1-34, doi:10.1680/geot.1979.29.1.1.
- Abstracts only [E-abstract, SA]: a field of hardening moduli on nested surfaces, with cyclic application.
- Piecewise-linear response with no smooth transition at a reversal. [E-sec, Minagawa et al. 1987 and Chen 2023 Sec. 2.4, SA]
- Taiebat & Dafalias (2015) group the stress-reversal models of Mroz et al. (1979) and Mroz & Zienkiewicz (1984) with models that have elastic neutral loading. [E]
- Projection-centre (homology-centre) jump and a finite modulus at the jump. [R]

### 4.4 Kan & Taiebat (2014)
E-Kan, M. & Taiebat, H.A. (2014). "On implementation of bounding surface plasticity models with no overshooting effect in solving boundary value problems." *Computers and Geotechnics* 55:103-116. doi:10.1016/j.compgeo.2013.08.006.
- Citation verified via Crossref [E-abstract, SA]. **The authors are M. E-Kan and H.A. Taiebat (UNSW); Khalili is NOT an author.** Companion paper: Kan, Taiebat & Khalili (2014), *Int. J. Geomech.* 14(2):239-253, doi:10.1061/(ASCE)GM.1943-5622.0000307, a mapping rule based on the last stress-reversal point [E-abstract, SA].
- Full text blocked (ScienceDirect/ResearchGate 403).
- Remedy, from a search snippet only [S]: "clouds" of loading surfaces with a strain margin inside which small unload-reload cycles do not move the homology centre, controlled by a threshold on accumulated plastic shear strain. Implemented in the UNSW bounding-surface model in an explicit finite-difference code and tested on monotonic and dynamic BVPs.
- Findings Chen et al. (2022) attribute to it [E-sec]:
  - Overshooting occurs under both implicit and explicit integration.
  - A smaller oscillation amplitude makes it worse.
  - Montans & Borja's virtual bounding surface suits Masing-type models.
- Newton-iterate noise: not verifiable. [I] Their host code was explicit, so there were no global Newton iterates.

### 4.5 Chen, Ghorbani, Zhang & Kodikara (2022) - see Q6.1b
This is the paper the brief cited as "C&G 152:105008": "Stress overshooting solution for soil plasticity models". Read via the author's thesis chapter. [E]

### 4.6 Papadimitriou & Bouckovalas (2002) - NTUA-SAND
"Plasticity model for sand under small and large cyclic strains: a multiaxial formulation." *Soil Dyn. Earthq. Eng.* 22(3):191-204. doi:10.1016/S0267-7261(02)00009-X.
- Paywalled; the NTUA DSpace record has metadata only. [E-abstract, SA]
- Triaxial predecessor: Papadimitriou, Bouckovalas & Dafalias (2001), JGGE 127(11):973-983. It combines a bounding surface with a Ramberg-Osgood small-strain formulation, and its Kp depends on accumulated plastic volumetric strain (fabric). [E-abstract, SA]
- The Itasca NTUA-SAND UDM info sheet lists a state variable "rijLR", the stress-ratio tensor at the last load reversal. [E, SA; https://itasca-software.s3.amazonaws.com/udm-library/NTUA-SAND_ShortInfo_0.pdf]
- TD08 notes that Papadimitriou et al. used an h different from the alpha_in form. [E]
- Exact Kp, trigger and bounding at reversal: not verified.

### 4.7 Andrianopoulos, Papadimitriou & Bouckovalas (2010a, b; 2005)
- (2010a) "Bounding surface plasticity model for the seismic liquefaction analysis of geostructures." *SDEE* 30(10):895-911. doi:10.1016/j.soildyn.2010.04.001.
- (2010b) "Explicit integration of bounding surface model for the analysis of earthquake soil liquefaction." IJNAG 34(15):1586-1614. doi:10.1002/nag.875.
- Abstracts read (University of Thessaly repository ir.lib.uth.gr/xmlui/handle/11615/25611 and /25612; Crossref) [E-abstract, SA]:
  - (a) A discontinuously relocatable projection centre tied to the LAST load-reversal point serves both the mapping and the reference of the Ramberg-Osgood "elastic" nonlinearity (no elastic region).
  - (b) Explicit integration with automatic error control and sub-stepping, checked with iso-error maps including Lode-angle changes and on a VELACS BVP.
- The Itasca UDM uses modified Euler (Sloan et al. 2001) with STOL and T_min. [E, SA]
- Reversal detection inside sub-steps and the relocation rule: not verified.
- Andrianopoulos et al. (2005), "Bounding surface models of sands: pitfalls of mapping rules for cyclic loading", 11th IACMAG pp. 241-248. [E-sec, Pisano & Jeremic 2014]

### 4.8 Loukidis & Salgado (2009)
"Modeling sand response using two-surface plasticity." *Computers and Geotechnics* 36(1-2):166-186. doi:10.1016/j.compgeo.2008.02.009.
- Not accessible; nothing verified. Treatment at reversal: unknown. [No claim made.]

### 4.9 Li (2002)
"A sand model with state-dependent dilatancy." *Geotechnique* 52(3):173-186. doi:10.1680/geot.52.3.173.41008.
- Abstract only. [E-abstract, SA]
- The projection centre jumps to the stress point at a reversal, and Kp depends on the ratio rho_bar/rho. Just after the jump rho = 0, so Kp is infinite only for that instant. [R] (This is the brief's own description; I could not verify it.)

### 4.10 Manzari & Dafalias (1997)
"A critical state two-surface plasticity model for sands." *Geotechnique* 47(2):255-272. doi:10.1680/geot.1997.47.2.255.
- Abstract only. [E-abstract, SA]
- **h does NOT depend on alpha_in:** h_MD = h0 |(alpha^b_theta - alpha):n| / [ b_ref - |(alpha^b_theta - alpha):n| ], with b_ref = (alpha^b_theta - alpha^b_theta+pi):n. [E-sec, Chen et al. 2022 eq. 3-18, read by me]
- TD08 also refers to MD97's alternative h choice. [E]
- Consequences, per Chen et al. [E-sec]:
  - Kp jumps when n flips.
  - After a reversal during softening, b_ref < |b:n|, so h < 0 and Kp << 0 - a spuriously stiff response and overshooting.
  - h is 0 on the bounding surface. [I]

### 4.11 Wang, Dafalias & Shen (1990)
"Bounding surface hypoplasticity model for sand." *J. Eng. Mech.* 116(5):983-1001. doi:10.1061/(ASCE)0733-9399(1990)116:5(983).
- Abstract only [E-abstract, SA]: a hypoplastic model in which the loading and plastic-strain-rate directions depend on the stress-rate direction.
- It is one of the works that applied the stress-rate-dependent mapping rule to sands. [E, Taiebat & Dafalias 2015]
- Kp as a distance ratio measured from a projection centre at the last reversal, later inherited by Li (2002). [R]

### 4.12 Other items found
- Petalas, A.L. (2025/26). IJNAG 50(2):597-612. Abstract only [E-abstract, SA]: the memory loss in MD97/DM04 is a constitutive choice, not something intrinsic to bounding-surface plasticity.
- Zhang, W., Lim, K., Ghahari, S.F., Arduino, P. & Taciroglu, E. (2021). IJNAG 45(8):1091-1119, doi:10.1002/nag.3194. A bounding-surface model in Abaqus with an "overshooting correction scheme" [E-abstract, SA]. Chen et al. class it as a plastic-shear-strain threshold scheme. [E-sec]
- Carow & Rackwitz (2021), *Comput. Geotech.* 140:104206: SANISAND-Z with backward Euler and a damped Newton adaptive trial step. [S]
- Chen, C. et al. (2024; online 2023), *Comput. Geotech.* 165:105899: implicit SANISAND-Z with sub-stepping. [metadata only]
- Montans, F.J. & Borja, R.I. (2002). IJNME 55(10):1129-1166: implicit J2 bounding-surface plasticity with a virtual bounding surface. [E-sec, Ghorbani 2023 and Chen 2022]
- Tseng & Lee (1983), *J. Eng. Mech.* 109(3):795-810: an elastic-unloading chord criterion. [E-sec, Chen 2022]

## Q5. Implementations

### 5.1 OpenSees ManzariDafalias (UW; Ghofrani & Arduino) - source read [E]
File SRC/material/nD/UWmaterials/ManzariDafalias.cpp (local fork copy). The lines quoted carry no `// Ladruno` marker, i.e. they are the upstream code. I re-checked each of them against upstream master, https://raw.githubusercontent.com/OpenSees/OpenSees/master/SRC/material/nD/UWmaterials/ManzariDafalias.cpp (5175 lines, fetched 2026-09-28) [E]:
- reversal test: lines ~864-876
- ModifiedEuler neutral-loading and negative-dGamma branches: ~1330-1345
- GetStateDependent h cap: ~4719-4722
- parser default IntScheme = 1: line 93

Findings:
- h, GetStateDependent(): x = (alpha - alpha_in):n. If |x| < small (small = 1e-10), h = 1.0e10; else h = b0 / x. [E] Note the fabs: a NEGATIVE x of magnitude above 1e-10 gives a NEGATIVE h, with no floor. [E] [I] x < 0 can occur inside a step or sub-step because alpha_in is frozen for the whole step while n rotates. Kp = (2/3) p h (b:n) then takes the sign of x*(b:n).
- Re-seat, integrate(): trialDirection = Ce : (eps - eps_n), the elastic trial stress increment from the committed state. If (alpha_n - alpha_in,n) : trialDirection < 0, alpha_in := alpha_n; otherwise alpha_in := alpha_in,n. [E] The code comment says it assumes a fully elastic step and checks whether the new stress direction differs strongly from the path, measured from the yield-surface centre. [E] [I] Because it contracts with the deviatoric (alpha_n - alpha_in,n), the test is effectively (alpha_n - alpha_in,n) : de < 0, a strain-space test.
  - It is re-evaluated at every Newton iterate from the committed (alpha_n, alpha_in,n). It can therefore flip between iterates, but it never accumulates.
  - Within a step, alpha_in is held fixed for all sub-steps. No re-seat occurs inside the integrator.
- H <= 0 handling:
  - In ModifiedEuler (default IntScheme 1; the parser default is oData[0] = 1), if |H| < 1e-10 the sub-step is treated as neutral loading: d(sigma) = 0, d(alpha) = 0, and the whole strain increment is booked as plastic. [E]
  - If dGamma < -1e-10 (e.g. H < 0 with a loading numerator), the code prints (debug only) "dGamma cannot be negative!", sets dGamma = 0, takes an elastic stress increment, translates alpha with the stress ratio (d(alpha) = d(s/p)), and switches to the elastic tangent. [E]
  - ForwardEuler sets H = +1e-10 if |H| < 1e-10, and passes a negative H through Macaulay(dGamma), which gives an elastic step. [E] A TODO comment in the code admits the H = 0 case is not handled properly.
  - [I] So vanilla UW silently converts the H <= 0 / Kp -> -infinity case into an elastic or rigid-translation step. That is an implementation patch, not a model rule.
- Implicit scheme (IntScheme 2) Jacobians use their own thresholds [E]:
  - GetJacobian (~4195): if |x| <= 1e-10, then h = 1e10, x = 1e-10, and the dh terms are dropped.
  - NewtonSol (~2925): x -> 1e-10 if |x| <= 1e-10, else |x|.
  - NewtonSol2 (~3332): x -> 1e-4 if x <= 1e-3, with dh dropped.
  - [I] So the linearization near a re-seat is inconsistent across code paths.
- Incidental bug (IntScheme 5, ForwardEuler only): `Vector r(6); if (p > small) Vector r = ...;` shadows r, which stays 0, so the K D n:r term drops out of that scheme. [E; also flagged by a sub-agent]

### 5.2 OpenSees SAniSandMS (Liu, Abell, Diambra & Pisano) - source read [E]
File SRC/material/nD/UANDESmaterials/SAniSandMS.cpp (vanilla, no fork markers). The upstream master copy is identical, 3023 lines; the lines below are at 1119-1123 and 2652-2663 (https://raw.githubusercontent.com/OpenSees/OpenSees/master/SRC/material/nD/UANDESmaterials/SAniSandMS.cpp) [E]. The OpenSees doc page credits H. Liu, J.A. Abell, A. Diambra and F. Pisano and describes RK4(5) with error control (opensees.github.io/.../SAniSandMS.html). [E]
- h: x = (alpha - alpha_in):n. If x < small (1e-10), x = small. This is a ONE-SIDED floor, so negative x is lifted to +1e-10. Then h = min(1.0e7, b0/x * exp(mu0 sqrt(p/p_at) (bM_distance/bref)^2)). [E]
- Memory-surface hardening: b_rM_rin = (r_alphaM - alpha_in):n, floored at 1e-7; hM = min(1.0e10, 0.5 b0/b_rM_rin + 0.5/sqrt(2/3) Z <-D>/b_bM). [E]
- Re-seat: the same trial-elastic-increment test as UW MD, (alpha_n - alpha_in,n) : [Ce:(eps - eps_n)] < 0, giving alpha_in := alpha_n. [E]
- Loading-index denominator: `if (fabs(temp4) < small) temp4 = small;` at every RK stage. [E]
- [I] Net effect: h is always positive and capped at 1e7. A re-seat with b:n < 0 therefore gives a large but finite negative Kp, about (2/3) p 1e7 (b:n). H can still go negative.

### 5.3 OpenSees PM4Sand / PM4Silt
- PM4Sand is covered in Q3.4. [E]
- PM4Silt.cpp follows the same pattern (Cg1 = h0/200; Kp >= 0 clip in the contraction branch). Its RungeKutta4 also resets Kp whenever the dGamma denominator is < 0. [E, SA]
- Chen, L. & Arduino, P. (2021). *Implementation, verification, and validation of the PM4Sand model in OpenSees.* PEER Report 2021/02. https://peer.berkeley.edu/sites/default/files/2021_chen_final.pdf. It gives the reversal as (alpha - alpha_in):n_trial < 0 (eq. 2.30), says C_gamma1 is added to keep the denominator from becoming zero, and sets Kp = 0 when b:n < 0. [E, SA]

### 5.4 Pisano SANISAND-MS PLAXIS UDSM (TU Delft; open source)
https://github.com/FedericoPisano/SANISAND-MS-UDSM (commit 205c13b)
- README eq. (read by me) [E]:
  - h = b0 / ( |(alpha - alpha_in):n| + C_eps ) * exp[mu0 (p/p_atm)^0.5 (b^M/b_ref)^2], with C_eps proportional to h0/100 and C_eps = 0.001, "see Boulanger and Ziotopoulou (2017)"
  - hM uses the same |.| + C_eps denominator, plus a sgn[...] variant
- Source (read by a sub-agent) [E, SA]:
  - StatedepSub.f90 L1243-1249: abs value plus 0.001
  - exp exponent capped at 9.2103 (about 1e4)
  - h = min(1e7, h), or min(user h_Max, h)
  - Integration.f90: re-seat only if (alpha - alpha_in):n < -tol with tol = 1e-6 (L83). n is taken at the yield-entry point, or at the current stress. The check runs only once plastic loading is established.
  - No Kp clamp observed.
- [I] This is the second production code, after PM4Sand, with an additive denominator constant. It explicitly borrows the idea from PM4Sand. The |x| + C_eps form also makes h positive and bounded for x < 0 inside a step.

### 5.5 Itasca P2PSand (Cheng & Detournay) - FLAC/FLAC3D online doc
https://docs.itascacg.com/itasca900/common/models/p2psand/doc/modelp2psand.html (the FLAC3D 7.0 page is identical). Read by me in the sub-agent's saved copy. [E]
- Kp, eq. (15): Kp = (2/3) h0 f1 f2 Dr G [(alpha^b_theta - alpha):n] / [(alpha - alpha_in):n]. The text notes Kp = 0 on the bounding surface. No regularizing constant or negative-Kp clamp is documented. [E]
- Property `ratio-reverse`: the minimum change of the back-stress ratio for the path to count as a reverse path, default 0.02. [E] Which metric is thresholded is not stated.
- [I] This is the only explicit reversal-recognition hysteresis found in a commercial code.
- The paper - Cheng, Z. & Detournay, C. (2021). "Formulation, validation and application of a practice-oriented two-surface plasticity sand model." *Computers and Geotechnics* 132:103984 - is open access (per OpenAlex, SA) but ScienceDirect blocked retrieval. Not read.

### 5.6 Itasca DM04 UDM (Z. Cheng)
- Manuals v3.2 (2019) and v4.0 (2023) come in the UDM zips at itasca-software.s3.amazonaws.com/udm-library/ (PDFs only; the DLLs were not downloaded). They list no h/Kp regularization parameter and no reversal parameter. The auxiliary properties are kcut (a low-pressure cut-off, default 0.01), flag-ini, and flag-origin (resets fabric to the origin "to avoid possible overshooting"). [E, SA]
- Cheng, Dafalias & Manzari (2013), Itasca symposium paper: not reachable (ResearchGate 403). The paper number is also unverified.

### 5.7 Legacy UC Davis DM04 (NewTemplate3Dep, 2005; Z. Cheng, B. Jeremic, M. Taiebat) - still carried in xara
https://github.com/peer-open-source/xara/blob/13cc747ffaab23ebea4debaa60268699e358105d/SRC/material/ucsd/NewTemplate3Dep/DM04_alpha_Eij.cpp#L199-L215 [E, SA]
- First call: alpha_in := alpha and h = 1e10 * b0.
- Afterwards, on the CURRENT (iterate) state: a_in = (alpha - alpha_in):n; if a_in < 0, re-seat alpha_in := alpha; if a_in < 1e-10, set a_in = 1e-10; then h = b0/a_in.
- [I] The re-seat mutates alpha_in during a hardening evaluation, i.e. inside Newton iterations and sub-steps - the most fragile variant found.

### 5.8 numgeo (Machacek et al.)
- The reference pages for Sanisand, Sanisand-2, Sanisand-F and Sanisand-MSf (https://j-machacek.github.io/numgeo/2026-09/reference/material/mechanical/sanisand.html and siblings) document no h cap, no regularization and no alpha_in rule. [E, SA]
- A 1D site-response tutorial uses an undocumented option `update_alpha, 1` (alpha updating strategy, detect load reversal). [E, SA]
- Sanisand-2 was implemented with M. Taiebat and S. Zeng (UBC). The benchmark used OpenSees results, and one drift option mentions "the OpenSees drift correction". [E, SA] [I] Likely OpenSees heritage.
- Machacek, Staubach, Tafili, Zachert & Wichtmann (2021), *Comput. Geotech.* 138:104276 (CC-BY) was not read: ScienceDirect returned 403 and the tuprints copy sits behind a bot challenge that was not attempted.

### 5.9 Others
- Real-ESSI sanisand2004 (Jeremic, UC Davis): the DSL manual exposes no h parameter and recommends explicit integration with strain sub-increments < 1e-4. The lecture notes (Sec. 104.6.11, eq. 104.402) give only the textbook rule (update alpha_in when the denominator becomes negative). The source is not public. [E, SA]
- soilmodels.com SANISAND UMAT/PLAXIS (Martinelli, Miriano & Tamagnini; updated by Masin): code not opened. [E-sec, soilmodels.com page]
- MATLAB port cgl-sd/SANISand04 (GitHub): its 1e10 branch is overridden by h = min(1e7, b0/x), so x < 0 still gives a negative h. The reversal rule is the OpenSees one. [E, SA]
- jghorbani2/SANISAND_rep (C++): HPARA = BREF / max(1e-15, x), where BREF is that code's b0-type factor, not TD08's b_ref. No cap. Re-seat if (alpha - alpha_in):n < 0 at the start-of-increment stress. [E, SA]
- Yu, Wang & Zhang (2020), "three integration schemes for SANISAND-04" (J-STAGE JGS Special Publication 8(3)): nothing on reversal or singularity. [E, SA]

## Q6. Explicit discussions: spurious reversals, thresholds, singular modulus with softening

### 6.1 The key paper: spurious reversals from numerical oscillations (open access, read in full)
Ghorbani, J., Chen, L., Kodikara, J., Carter, J.P. & McCartney, J.S. (2023). "Memory repositioning in soil plasticity models used in contact problems." *Computational Mechanics* 71:385-408. doi:10.1007/s00466-022-02245-z. https://link.springer.com/content/pdf/10.1007/s00466-022-02245-z.pdf [E]
- Model: MUD. On saturation its hardening law degenerates to the SANISAND one (their conclusions):
  - h = h0 G0 (1 - ch e)(p'/p_atm)^(-1/2) / [(alpha_k - alpha_in)^T n] (eq. 8)
  - Re-seat: if (alpha_k - alpha_in)^T n < 0 at the start of step i+1, then alpha_in^(i+1) = alpha_k^(i) (eq. 9)
  - Kp = (2/3) p' h (alpha^b - eta)^T n (eq. 10)
  - [E]
- Problem statement:
  - "Spurious numerical oscillations" trigger tiny reversal events followed by reloading, so the model becomes unrealistically stiff and overshoots. [E]
  - In dynamic and contact problems the oscillations cannot be avoided, and analyses terminate early. [E]
  - Their Fig. 1 example: a single -0.00017 axial-strain increment at eps_a = 0.08 in an undrained triaxial test produces a drastic jump in q. [E]
  - When the induced reversal occurs, (alpha - alpha_in):n < 0 is observed first and then reset to zero by eq. (9), giving a very large Kp (their Fig. 2 discussion). [E]
- Their literature map: some remedies prevent drastic Kp changes when spurious oscillations trigger a reversal - Dafalias 1986; Kan & Taiebat 2014; Dafalias & Taiebat 2016; Barrero et al. 2020; Limnaiou & Papadimitriou 2022; Duque et al. 2022; Ziotopoulou & Boulanger 2016. Montans & Borja (2002) use a virtual bounding surface (J2 bounding-surface only). Typically a reversal counts as true only if the deviatoric plastic strain accumulated during it exceeds a threshold (their ref. [15] = SANISAND-Z). [E for what they state]
- Scheme of Dafalias & Taiebat (2016), as reproduced (eqs. 11-12):
  - alpha_in^(i+1) = alpha_k^(i) + m_q (alpha_in^r - alpha_k^(i)) (eq. 11)
  - m_q = <1 - (eps_q^p / eps_bar_q^p)^j> (eq. 12; the exponent placement is ambiguous in their typesetting, but with the default j = 1 it makes no difference)
  - eps_q^p = sqrt(2/3 e_q^p : e_q^p), the deviatoric plastic strain accumulated during the reversal step (i); eps_bar_q^p is a threshold
  - A trivial unloading (m_q -> 1) leaves the memory essentially unchanged. A large one (m_q -> 0) reverts to the plain re-seat, with an infinite Kp.
  - It needs an extra stored tensor, alpha_in^r.
  - [E-sec for SANISAND-Z; E for the reproduction]
- Stated limitation: repositioning cannot guarantee (alpha - alpha_in):n > 0 afterwards (they say this was highlighted in [15]). When it fails, the model must degenerate to eq. (9). [E]
- Limnaiou & Papadimitriou (2022), as reproduced: every reversal is treated as informal (alpha_in NOT updated) until the difference between the accumulated stress ratio and alpha_in exceeds a prescribed tolerance; only then is it formal and repositioned. [E-sec]
- Their new scheme: J = (alpha_k - alpha_in)^T n; after a reversal set J_1^(i) = J^r m_q (eq. 16). J^r is a scalar state: the last positive value of J before the current loading process was reversed. Equivalently alpha_in^(i) = alpha_k^(i) - J^r m_q n^(i) (eq. 17). [E]
  - For trivial reversals (oscillations) the distance, and so h, stays finite and continuous.
  - For true reversals (eps_q^p >= eps_bar_q^p, so m_q = 0) the original infinite-Kp re-seat is recovered.
  - Only a scalar is stored.
- Integration: explicit modified-Euler/Euler with error control (STOL). The reversal test (alpha_k - alpha_in)^T n < 0 runs inside the sub-stepping. An odd/even counter "Lstepcount" alternates between the reversal step (store eps_q^p) and the reloading step (reposition). [E]
- Documented failure mode with no known remedy: a tiny first reversal is NOT detected ((alpha - alpha_in):n stays > 0), but the following reloading IS flagged as a reversal. Both schemes then misfire. The authors call it rare and manageable through STOL and step size. [E]
- FE contact-impact results: both repositioning schemes reduced iterations and CPU time. Larger contact-penalty coefficients produced stronger oscillations. [E]
- [I] Directly relevant: (i) it is the only source found that frames the alpha_in re-seat as a numerical-robustness problem in BVPs; (ii) the scalar-floor idea (keep (alpha - alpha_in):n >= J^r m_q > 0 after a non-genuine reversal) removes the h = infinity spike without changing the model for genuine reversals.

### 6.1b The most directly relevant analysis: reversal during softening, and a SANISAND04 footing BVP that aborts
Chen, L., Ghorbani, J., Zhang, C. & Kodikara, J. (2022). "Stress overshooting solution for soil plasticity models." *Computers and Geotechnics* 152:105008. doi:10.1016/j.compgeo.2022.105008.
- Read as Chapter 3 of Chen, L. (2023), *Modelling of hydro-mechanical shakedown and ratcheting of unsaturated granular materials*, PhD thesis, Monash University. The chapter reproduces the published paper. PDF saved by a sub-agent; I rendered and read pp. 3-5 to 3-46. [E] Equation numbers below are the THESIS numbers (3-xx); the paper's own numbering was not checked.
- Review of remedies (p. 3-5) [E]:
  - Fardis et al. (1983): store the minimum distance ratio. This creates a Kp discontinuity.
  - Montans & Borja (2002): a virtual bounding surface, suited to Masing-type models.
  - Tseng & Lee (1983): an elastic-unloading chord-length test, suited to metals.
  - Accumulated-plastic-shear-strain thresholds (Dafalias 1986; Dafalias & Taiebat 2016; Kan & Taiebat 2014; Zhang et al. 2021). A reversal counts as genuine only above the threshold. These remain prone to Kp discontinuities.
- **Reversal during softening (Sec. 3.6.1)** [E]:
  - MD97 hardening, eq. (3-18): h_MD = h0 |(alpha^b_theta - alpha):n| / [ b_ref - |(alpha^b_theta - alpha):n| ], with b_ref = (alpha^b_theta - alpha^b_theta+pi):n.
  - After softening, where (alpha^b - alpha):n < 0, a load reversal can give b_ref < |(alpha^b - alpha):n|. Then h_MD < 0 and Kp = (2/3) p' h_MD (alpha^b - alpha):n < 0, far below zero, with <d(lambda)> = 0. The model becomes unrealistically stiff and the stress overshoots.
  - Their drained-triaxial demonstration: loose Karlsruhe sand, reversal at 12% axial strain.
  - [I] This is the only published analysis found of a reversal coinciding with alpha outside the bounding image. Theirs is the MD97 variant of our failure mode.
- **SANISAND04, eq. (3-19): h_SAN = G0 h0 (1 - ch e)(p'/p_a)^(-1/2) / [(alpha - alpha_in):n].** A reversal-reloading drives (alpha - alpha_in):n to 0 and h to infinity. In their triaxial "Path A" with small oscillation amplitudes (0.2%, 0.05%, 0.001%), overshooting gets WORSE as the oscillation amplitude shrinks. The same holds for MD97 (Figs. 3-5, 3-7). [E] They note Kan & Taiebat (2014) reached the same conclusion. [E-sec]
- **SANISAND-ZO (Dafalias & Taiebat 2016 weighting)** - their reconstruction, eqs. (3-20)-(3-21) [E]:
  - alpha_in^(i+1) = m alpha_in^(i-1) + (1 - m) alpha^(i)
  - m = <1 - (eps_q^p(i) / eps_bar_q^p)^j> (3-21), with j = 1 by default; Chen et al. used the Dafalias & Taiebat (2016) values j = 1 and eps_bar = 0.01%. Their failure demo used strain increments of 3e-6 versus 1e-7.
  - Their Cases 2 and 3 (Fig. 3-8): after the weighted update, (alpha - alpha_in):n can still be < 0. Per SANISAND-Z the model then re-seats fully, which again gives h = infinity. Chen et al. say Dafalias & Taiebat acknowledge Case 2.
  - With smaller strain increments or tighter stress tolerance the overshooting gets WORSE (Figs. 3-9, 3-10). Their explanation: h can take two very different values (finite and infinite) at the same state, and small sub-steps resolve that discontinuity instead of stepping over it.
- **Footing BVP (Sec. 3.9.1)** [E]:
  - Setup: 2D plane-strain flexible footing on loose Karlsruhe sand (e0 = 0.98), static, their coupled FE code (288 quadratic elements).
  - SANISAND04 with time increments equivalent to 40, 400 and 4000 steps: the coarse run completes. With finer steps, overshooting becomes pronounced and the ANALYSIS ABORTS, "caused primarily by the sudden reduction of (alpha - alpha_in):n to zero". Tighter STOL also makes the aborts worse.
  - SANISAND-ZO fixes the coarser cases but fails at the finest increment.
  - Their proposed model completes all cases.
- **Their proposed remedy - memoryless, strictly positive, bounded h** (Sec. 3.7) [E]:
  - h_p = (b0/theta) * ln(1 + exp((alpha^b_theta - alpha):n)) / ln(1 + exp((alpha - alpha^b_theta+pi):n))   (3-22)
  - b0 = G0 h0 (1 - ch e)(p'/p_a)^(-1/2)   (3-23)
  - theta = tanh(omega1 * integral|d eps_q| + omega2)   (3-24), with omega1 = 1e4 and omega2 either 0 (elastic-like start of virgin shearing) or 1e6 (feature off)
  - Inside h, alpha_in is replaced by the opposite bounding image alpha^b_theta+pi, so h itself needs no reversal detection.
  - The softplus ln(1 + e^x) keeps h > 0, and it never gives h = infinity or Kp = infinity at a reversal.
  - A smoothing step handles the remaining Kp discontinuities (eqs. 3-25 to 3-27):
    - Kp = (2/3) p' b_h1:n and d(alpha) = (2/3) <d(lambda)> b_h1
    - b_h1 = (||b_h|| F + ||b_h||^r (1 - F)) n_h, with b_h = h_p (alpha^b_theta - alpha)
    - ||b_h||^r is the value at the last reversal, still detected by the DM04 test (alpha - alpha_in):n < 0
    - F = tanh(omega1 (alpha - alpha_in):n), so Kp is continuous across a reversal
    - So a DM04-style alpha_in is still tracked, but only for this continuity blend, not for the magnitude of h.
  - No new material parameters. They re-checked it against Karlsruhe and Toyoura data (supplementary material, not read).
  - [I] The numerator softplus keeps h > 0 even when (alpha^b - alpha):n < 0. The sign of Kp still comes from the (alpha^b - alpha):n factor in Kp = (2/3) p h b:n, so softening survives, but finite: |Kp| <= (2/3) p (b0/theta) ln(1 + e^(b:n)) |b:n| / ln 2.
- Companion conference paper: Chen, L., Ghorbani, J., Zhang, C. & Kodikara, J. (2022/2023). "A robust solution to address overshooting in bounding surface plasticity models." IACMAG 2022, *Lecture Notes in Civil Engineering* 288:71-78. Paywalled; abstract via search snippet only. [E-sec]

### 6.2 Other explicit statements on reversal-detection robustness
- Pisano, F. & Jeremic, B. (2014). "Simulating stiffness degradation and damping in soils via a simple visco-elastic-plastic model." *Soil Dyn. Earthq. Eng.* 63:98-109. Preprint read: https://sokocalo.engr.ucdavis.edu/~jeremic/wwwpublications/CV-J32.pdf [E]
  - The Borja & Amies (1994) unloading criterion - a reversal whenever the hardening distance kappa starts to increase - proved not robust under irregular loading. Near the bounding surface kappa is small and "can be easily corrupted even by numerical inaccuracies", which triggers spurious unloading. [E]
  - They adopt (sigma - sigma_0) : dr < 0, equivalently (sigma - sigma_0) : n_dev < 0 (eq. 20). They note it produces overshooting (Dafalias 1986) and point to E-Kan & Taiebat (2014) for remediation, and to Andrianopoulos et al. (2005, 2010) for robustness alternatives. [E]
- Jeremic, B., Cheng, Z., Taiebat, M. & Dafalias, Y.F. (2008). "Numerical simulation of fully saturated porous materials." IJNAG 32(13):1635-1660 [volume/pages R]. Preprint: https://sokocalo.engr.ucdavis.edu/~jeremic/wwwpublications/CV-J20.pdf [E]
  - With explicit integration, a load-reversal step large enough to miss the elastic region lands on the opposite side of the yield surface. The step is then evaluated with derivatives from the wrong side and can be completely erroneous. [E]
  - [I] This is acute for a thin cone (m = 0.005), where almost any reversal step jumps across the cone.
- Andrianopoulos, K.I., Papadimitriou, A.G. & Bouckovalas, G.D. (2005). "Bounding surface models of sands: pitfalls of mapping rules for cyclic loading." Proc. 11th IACMAG, pp. 241-248. Not opened [E-sec via Pisano & Jeremic 2014].
- Duque, J., Yang, M., Fuentes, W., Masin, D. & Taiebat, M. (2022). "Characteristic limitations of advanced plasticity and hypoplasticity models for cyclic loading of sands." *Acta Geotech.* 17:2235-2257. doi:10.1007/s11440-021-01418-z (eScholarship copy read: https://escholarship.org/content/qt1vh8z8hp/qt1vh8z8hp.pdf) [E]
  - Their Limitation 1 is overshooting after reverse loading followed by immediate reloading. For DM04 and SANISAND-MSf it is attributed to the discrete memory variable being updated at any stress reversal. [E]
  - The bounding surface limits how far the overshoot can go. [E]
  - Suggested remedy: memorize recent loading history - Dafalias (1986), with implementation details in Dafalias & Taiebat (2016) - added to SANISAND-MSf in their Fig. 3. [E] (Their ref. [6] mis-cites "Dafalias 1986" as Mech. Res. Commun. 13(6); DM04 pages are also mis-cited. [E])
- Tafili, M., Duque, J., Masin, D. & Wichtmann, T. (2024). "Repercussion of overshooting effects on elemental and finite-element simulations." *Int. J. Geomech.* 24(3). doi:10.1061/IJGNAI.GMENG-8842. Abstract only [E-abstract]: overshooting is among the most serious limitations of the models studied (one bounding-surface model) and has a major impact on FE simulations. Full text paywalled.
- Limnaiou, T.G. & Papadimitriou, A.G. (2023). "Bounding surface plasticity model with reversal surfaces for the monotonic and cyclic shearing of sands." *Acta Geotech.* 18:235-263. doi:10.1007/s11440-022-01529-1. Abstract read (Springer landing page; full text paywalled) [E-abstract]
  - SANISAND-type model with no (small) yield surface. The last stress reversal point defines both the elastic and the plastic strain rates, and the paper emphasises updating the reversal point to avoid overshooting.
  - Companion paper: the SDEE (2022) verification in BVPs with an explicit finite-difference code (ScienceDirect PII S0267726122002433; abstract via search only [E-sec]).
- Papadimitriou, A.G., Chaloulos, Y.K. & Dafalias, Y.F. (2019). "A fabric-based sand plasticity model with reversal surfaces within anisotropic critical state theory." *Acta Geotech.* 14:253-277. doi:10.1007/s11440-018-0751-5. Abstract read [E-abstract]: the last stress reversal point is the projection centre for the image stress on the bounding surface (monotonic focus). Full text paywalled.

### 6.3 Model-intrinsic precedents for removing a singular denominator (from the Dafalias school)
- SANISAND-MSf: Yang, M., Taiebat, M. & Dafalias, Y.F. (2022). "SANISAND-MSf: a sand plasticity model with memory surface and semifluidised state." *Geotechnique* 72(3):227-246. doi:10.1680/jgeot.19.P.363. eScholarship copy: https://escholarship.org/uc/item/0227z2t1 [E]. Their pp. 232-233 describe two traps in the memory-surface (MS) hardening:
  - (i) During dilation and softening, alpha - and hence the MS image - lies outside the bounding surface (BS). The MS size could then go negative.
  - (ii) Solving the consistency condition for hM puts (alpha^b_theta - alpha^M_theta):n in the denominator. That quantity can vanish when part of the MS lies outside the BS, giving an infinite hM, which they warn can cause serious numerical problems in implementations. [E] They note that Liu et al. (2019) contains such a singular case. [E]
  - The fix, eq. (9), is to reformulate the evolution law rather than cap it. Macaulay brackets go on the kinematic term. The dilation term is multiplied by |(alpha^b_theta - alpha^M_theta):n|, so the would-be zero denominator cancels in the closed-form hM, eq. (10), which uses a sgn(...) term. [E]
  - [I] This is the closest published analogue of a model-intrinsic cure for a 0/0 in a SANISAND hardening law: multiply the evolving quantity by the factor that would otherwise divide.
- SANISAND 2008 (Taiebat & Dafalias, IJNAG 32:915-948) removed alpha_in altogether, eqs. (17a)/(35), see Q2. The stated motivation includes the ease of implicit implementation. [E]
- PM4Sand's C_gamma1 (Q3) is the additive regularization used in production codes. [E]

### 6.4 Reversal-recognition rules found (hysteresis and thresholds), weakest to strongest
1. **Plain sign test**: (alpha - alpha_in):n < 0. Used by DM04, PM4Sand eq. (19), SANISAND-MS/MSf and Liu 2019. [E]
2. **Where and when the test is evaluated**:
   - Legacy UCD DM04 (xara): on the current iterate, inside the hardening function. [E, SA]
   - OpenSees UW MD and SAniSandMS: once per call, from the committed state, using the elastic trial increment direction Ce:d(eps). [E]
   - PM4Sand/FLAC: once per step, on the zone-averaged committed state, with Kp lagged one step. [E]
   - Chen et al. 2022: inside every sub-step rate evaluation (used only to reset the smoothing memory). [E]
3. **Tolerance on the sign test**:
   - Pisano UDSM: re-seat only if x < -1e-6, with n taken at yield entry. [E, SA]
   - P2PSand: `ratio-reverse` = 0.02, a minimum change of the back-stress ratio before a path counts as reversed. [E]
4. **Unloading alone is not a reversal**: SANISAND-F re-seats only if the NEW n gives x <= 0. An elastic unloading/reloading near the same point therefore keeps alpha_in. [E]
5. **Plastic-strain threshold with weighted repositioning**:
   - Dafalias (1986) / SANISAND-Z: m_q. [E-sec]
   - Kan & Taiebat (2014): a strain margin (loading-surface cloud). [S]
   - Zhang et al. (2021). [E-abstract]
   - Known failure: after the weighted update x can still be <= 0, which forces a full re-seat and h = infinity (Chen et al. 2022 Cases 2-3; Ghorbani et al. 2023). The problem gets worse with finer increments and tighter tolerances (Chen et al.). [E]
6. **Formal/informal reversals with a tolerance on the accumulated stress-ratio change** (Limnaiou & Papadimitriou 2022). [E-sec]
7. **Guaranteed-positive repositioning**: x := J^r m_q > 0 after trivial reversals (Ghorbani et al. 2023). [E]
8. **Memoryless h (no reversal needed for h)**: TD08 eq. (35); Chen et al. 2022 eq. (3-22); MD97. [E / E-sec]

### 6.5 What was searched for and NOT found (absence claims, over the sources opened)
- No source analyses the DM04-specific coincidence of a re-seat (x -> 0, h -> infinity) with (alpha^b_theta - alpha):n <= 0, whether as the product infinity * (negative) or as the 0/0 limit. [E-absence; confirmed independently by the three sub-agents]
  - Closest: Chen et al. (2022), reversal during softening for MD97 (h < 0, Kp << 0).
  - Closest: PM4Sand, where C_gamma1 and the Kp = 0 floor coexist as separate provisions.
  - Closest: SANISAND-MSf, which cancels a singular hM denominator, but for the memory surface, not h.
- No source treats alpha_in re-seats caused by GLOBAL Newton iterates.
  - Ghorbani et al. (2023) treat oscillations in the time/load history.
  - Pisano & Jeremic (2014) treat numerical inaccuracy near the bounding surface.
  - Chen et al. (2022) and Kan & Taiebat (2014) treat oscillation amplitude, step size and STOL.
- A Macaulay-bracketed denominator exists only in SANISAND-F's H = (2/3) h / <x>. It is a formal device that makes H = +infinity at the re-seat, NOT a regularization. [E]
- No published minimum value of x was found.
  - Code floors exist: 1e-10 (OpenSees UW, SAniSandMS, legacy UCD) and 1e-15 (jghorbani2).
  - Additive constants exist: C_gamma1 = h0/200 (PM4Sand) and C_eps = 0.001 (Pisano UDSM).

## Q7. Hashiguchi subloading surface (conceptual alternative)

Source opened: Hashiguchi, K. (2015). "Complete formulation of the subloading surface model." *Proc. VI Int. Conf. on Computational Methods for Coupled Problems in Science and Engineering (COUPLED PROBLEMS 2015)*, B. Schrefler, E. Onate & M. Papadrakakis (eds.), CIMNE, pp. 837-848. https://upcommons.upc.edu/bitstreams/cca5560b-2368-4c94-8aff-5a01b6671316/download [E]

- The subloading surface f(sigma_bar) = R F(H) (eq. 11) always passes through the current stress and is similar to the normal-yield surface about a similarity centre (the "elastic core"). R, with 0 <= R <= 1, is the normal-yield ratio, so no purely elastic domain exists. [E]
- R evolves as dR/dt = U(R) ||d^p|| (eq. 16), with U(R) -> +infinity as R -> 0, U > 0 for R < 1, U = 0 at R = 1 and U < 0 for R > 1 (eq. 17). Hashiguchi uses U(R) = u cot(pi R / 2) (eq. 18). [E] A value of R above 1, produced by a finite step, is therefore pulled back automatically. [E]
- The consistency condition (eqs. 40-41) gives a plastic modulus M^p (eq. 43) that contains a U(R)/R term. That term blows up only as R -> 0, i.e. when the stress sits exactly at the elastic core. It does not blow up at a reversal. [E for the presence of the term; I for the reading]
- At a reversal the stress does not jump back to a projection centre. It keeps its current R (R decreases elastically during unloading and grows again on reloading), so the modulus varies continuously. The advantages Hashiguchi lists include a smooth elastic-plastic transition that always meets the smoothness condition, and no need for a yield judgment. [E]
- Loading criterion (eq. 48): d^p != 0 when n:E:d > 0, i.e. in strain space with the elastic trial rate. Hashiguchi notes that a criterion written with the stress-rate plastic multiplier cannot be used for softening. [E] [I] This is the same structural point as the DM04 denominator H = Kp + 2G - K D n:r: the strain-driven problem is well posed only while H > 0.
- Over-stiff reloading after reversal (the Masing-type behaviour) is handled by translating the elastic core c and by making the effective u depend on n : n_hat_c (reloading vs reverse loading; Sec. 3.11). There is no discrete memory reset. [E]
- [I] What to take for DM04: a reversal-free "distance" variable that evolves continuously (a subloading-type R, or a continuously relaxed alpha_in) removes both the reset discontinuity and the singular modulus. The cost is changing the model's small-strain and reloading response, so it would need recalibration.

## Access notes: what was paywalled or unreachable (no paywall, captcha or bot challenge was bypassed)
- **Paywalled, not read:**
  - Dafalias & Manzari (2004) JEM (ASCE)
  - Ziotopoulou & Boulanger (2016) SDEE (Elsevier)
  - Limnaiou & Papadimitriou (2023) Acta Geotech. (abstract only) and their SDEE 2022 verification paper
  - Papadimitriou, Chaloulos & Dafalias (2019) Acta Geotech. (abstract only)
  - Tafili et al. (2024) IJG (abstract only)
  - Chen, Ghorbani, Zhang & Kodikara (2022) IACMAG chapter (snippet only)
  - Marinelli et al. (2026) SANISAND-MSf implementation chapter (abstract stub only)
- **Blocked, not read:**
  - PM4Sand v3.1 manual: Cloudflare challenge on pm4sand.engr.ucdavis.edu; box.com share; Scribd preview only.
  - Zenodo "IC MAGE Model 11 - SANISAND": files restricted.
  - ResearchGate: HTTP 403.
  - IA Scholar: anti-scrape puzzle, not attempted.
  - eScholarship blocks curl (403), but WebFetch retrieved the PDFs.
  - Cheng & Detournay (2021) and Machacek et al. (2021): open access, but ScienceDirect or the repository returned 403 or a bot challenge (SA).
  - Cheng, Dafalias & Manzari (2013) Itasca symposium paper: ResearchGate 403 (SA).
  - Ghofrani's UW thesis: not located (SA).
- **Also paywalled or blocked, reported by the sub-agents (SA):**
  - DPL04 (ASCE 403)
  - Li & Dafalias (2012) and Gao et al. (2014), journal versions
  - SANISAND-Z (2016), Geotechnique version
  - SANISAND-FN full text (Wiley/Durham, Cloudflare)
  - Barrero et al. (2020) full text (Wiley 403)
  - UBC theses (Yang; Reyes): repository flagged "unusual activity"
  - Dafalias & Popov (1975, 1976)
  - Dafalias (1986)
  - Mroz et al. (1978, 1979)
  - Papadimitriou & Bouckovalas (2002)
  - Andrianopoulos et al. (2010a, b) full texts
  - Li (2002); Wang et al. (1990); MD97
  - Kan & Taiebat (2014) (ScienceDirect/ResearchGate 403)
  - Loukidis & Salgado (2009) (no abstract retrieved)
  - Dafalias, Petalas & Feigenbaum (2023) IJSS (CC-BY, but behind a bot check)
  - An unidentified 2018 IJNAG paper at par.nsf.gov (connection refused)
- **Sub-agent tag convention:** "[E, SA]" marks a source that a sub-agent opened in full in this session and that I did not re-open myself. "[E-abstract]" means only an authoritative abstract or metadata record was read. "[S]" means only a search-engine snippet, because the page returned 403. Every other [E] was opened by me.

## Final table

x = (alpha - alpha_in):n. b:n = (alpha^b_theta - alpha):n. "BS" = bounding surface.

| Model / code | h / Kp form (source eq.) | Denominator regularization (param, value) | alpha_in / reversal rule (threshold? memory?) | Kp sign handling outside BS | Tag |
|---|---|---|---|---|---|
| Dafalias & Popov 1975/76 | Kp = Kp_bar + h delta/(delta_in - delta) | none; Kp = infinity at re-yield (intended) | delta_in reset at each new plastic loading; metals remedy Petersson & Popov 1977 (weight W on cumulative plastic strain) | n/a (metals, hardening) | R; E-sec |
| Dafalias 1986 | alpha_in-type (as adopted by DM04) | none; infinity intended | x < 0; overshooting remedy: alpha_in weighted by m = <1 - (eps_q^p/eps_bar)^j> | not verified | E-sec |
| Mroz, Norris & Zienkiewicz 1978/79 | field of moduli on nested surfaces | finite, piecewise | centre jump at reversal (R) | not verified | E-abstract; R |
| Manzari & Dafalias 1997 | h = h0 \|b:n\| / (b_ref - \|b:n\|), b_ref = (alpha^b_theta - alpha^b_theta+pi):n (Chen eq. 3-18) | no alpha_in; singular only if \|b:n\| -> b_ref | memoryless; Kp jumps when n flips | reversal during softening gives h < 0 and Kp << 0 (spurious stiffening) | E-sec (Chen 2022) |
| Wang, Dafalias & Shen 1990 / Li 2002 | rho_bar/rho from a projection centre at the last reversal | none; infinite only at the instant of the jump | centre jumps to the stress point | not verified | R; E-abstract |
| **DM04** | h = b0/x; Kp = (2/3) p h b:n | none; h = infinity at re-seat (intended) | alpha_in := alpha when x < 0; no threshold, no memory | negative Kp intended (softening); product -> -infinity at a re-seat with b:n < 0 (I) | E-sec |
| DPL04 | DM04-type with A = F:n dependence | unknown | unknown (DM04, R) | unknown | E-abstract; R |
| SANISAND 2008 (TD08) | h = b0 / [(3/2)((b_ref - b):n)^2] (eq. 35; triaxial 17a) | alpha_in ELIMINATED; squared fixed reference | none needed | finite negative Kp (softening allowed) | E |
| SANISAND-Z 2015/16 | Kp = (2/3) p h [(r^b - r):n] / [(r - r_in):n] (eq. 6) | none | r_in := r when (r - r_in):n <= 0; Dafalias-1986 threshold weighting (eps_bar = 0.01%, j = 1 per Chen) | image defined inside/on/outside BS (abstract) | E (conf.); E-sec |
| SANISAND-F 2020 | H = (2/3) h / <x> (17); Kp = p H b:n (19) | Macaulay makes H = +infinity at x <= 0 (not a regularization) | re-seat when x <= 0, including n-rotation during loading; elastic unloading alone does NOT re-seat | Kp -> 0, then negative outside BS | E |
| SANISAND-MS 2019 | h = b0/((r - r_in):n) exp[mu0 (p/pa)^0.5 (b^M/b_ref)^2] (11) | none in paper | r_in := r when (r - r_in):n < 0; overshooting acknowledged, untreated | r-form can give artificial softening (thesis) | E; E, SA |
| SANISAND-MSu 2020 | eq. 16 = MS form x (1/eta)^w2 | eta >= m only | alpha_in at the stress-increment reversal | Macaulay on b:n in MS shrinkage; singular hM flagged, left to the code | E, SA |
| SANISAND-MSf 2022 | h = b0/x exp[mu0/(\|\|alpha_in\|\|^u + eps) (b^M/b_ref)^w] (11) | eps = 0.01 guards only \|\|alpha_in\|\|; \|.\|, <.>, sgn cancel the singular hM (9-10) | x < 0 (per Dafalias 1986); no overshooting rule in paper | alpha outside BS in softening is standard DM04 behaviour | E |
| PM4Sand v3.3 / PM4Silt v2.1 | Kp = G h0 sqrt(b:n) / (exp(x_app) - 1 + C_gamma1) * C_rev (37-38) | C_gamma1 = h0/200 | x < 0; alpha_in^p, alpha_in^true, alpha_in^app (component min/max, all history), C_rev; alpha_in at init limited to 0.9 Mb | Kp = 0 outside BS (v3.3); negative in v2/v3; PM4Silt keeps the stress ratio inside max(BS, DS) | E |
| Chen, Ghorbani, Zhang & Kodikara 2022 | h_p = (b0/theta) softplus(b:n) / softplus((alpha - alpha^b_theta+pi):n) (3-22); Kp via tanh-smoothed b_h1 (3-25..27) | alpha_in removed from h; softplus > 0 | x < 0 only resets the smoothing memory \|\|b_h\|\|^r | finite negative Kp (sign from b:n) | E (thesis ch. 3) |
| Ghorbani et al. 2023 (MUD) | SANISAND-type h = b0'/x (8) | after a trivial reversal x := J^r m_q > 0 (17) | x < 0 plus plastic-strain threshold m_q; scalar memory J^r | not addressed | E |
| Kan & Taiebat 2014 (UNSW) | Kp = H_f + H_b, H_f -> infinity at reversal | not verified | plastic-shear-strain margin (loading-surface cloud) | not verified | S; E-sec |
| Limnaiou & Papadimitriou 2022/23 | reversal surfaces, no yield surface | not verified | formal/informal reversals with a tolerance | not verified | E-abstract; E-sec |
| NTUA-SAND / Andrianopoulos 2010 | projection from the last-reversal point (rijLR) | not verified | relocatable centre at the last reversal | not verified | E-abstract; E, SA |
| Itasca P2PSand | Kp = (2/3) h0 f1 f2 Dr G b:n / x (15) | none documented | `ratio-reverse` = 0.02 (minimum back-stress-ratio change) | Kp = 0 on BS; no clamp documented | E |
| OpenSees UW ManzariDafalias | h = b0/x | h = 1e10 iff \|x\| < 1e-10; x < 0 gives h < 0 | (alpha_n - alpha_in,n):(Ce:d(eps)) < 0, per call, committed state | no clamp; dGamma < 0 -> elastic step with alpha dragged; \|H\| < 1e-10 -> neutral-loading branch | E |
| OpenSees SAniSandMS | h = min(1e7, b0/max(x, 1e-10) exp(...)) | floor 1e-10; cap 1e7 | as UW MD | no clamp; \|H\| < 1e-10 -> 1e-10 | E |
| OpenSees PM4Sand/PM4Silt | v3.1-type, C_rev = (<x_app> + Cg1)/(<x_true> + Cg1) | Cg1 = h0/200; Macaulay on x | (alpha - alpha_in^true):n_tr < 0 (trial stress); app/p memory | Kp >= 0 only in the contraction branch | E |
| Pisano SANISAND-MS UDSM (PLAXIS) | h = b0/(\|x\| + C_eps) exp(<= 9.21), h <= 1e7 or h_Max | C_eps = 0.001 (~ h0/100) | x < -1e-6; n at yield entry | none observed | E; E, SA |
| Legacy UCD DM04 (xara) | h = b0/max(x, 1e-10); first call h = 1e10 b0 | floor 1e-10 | re-seat on the current iterate inside the hardening function | none | E, SA |
| Itasca DM04 UDM / numgeo | not documented | not documented | not documented (numgeo: undocumented `update_alpha` option) | not documented | E, SA |
| Hashiguchi subloading surface | M^p with a U(R)/R term (43) | finite except R -> 0 (elastic core) | no discrete memory; strain-space loading n:E:d > 0 (48) | R > 1 pulled back (U < 0) | E |

