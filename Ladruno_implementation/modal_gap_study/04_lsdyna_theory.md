# LS-DYNA modal / eigenvalue / frequency-domain theory dossier

Source manuals (PDFs on this machine, `…/Lit Review/Explicit/LS-DYNA/`):
- **Theory** = `LS-DYNA Theory.pdf` (680 pp).
- **Vol I** = `LS-DYNA_Manual_Vol_I_R16.pdf` (keyword reference, R16; cited by PDF page).

Citations give manual + section/keyword + PDF page. Quotes are short paraphrases; the
manuals are copyrighted. This dossier feeds the OpenSees "Ladruno" modal-feature gap study.

---

## 1. Implicit eigensolver (the core)

**Method.** The default eigensolver is **Block Shift-and-Invert Lanczos**, taken from
**BCSLIB-EXT** (Boeing's Extreme Mathematical Library) — Theory §36, PDF p.629
("Block Shift and Invert Lanczos eigensolver from BCSLIB-EXT"). Confirmed in Vol I
`*CONTROL_IMPLICIT_EIGENVALUE`, Remark 2 (PDF p.299 of that keyword block = manual
p.12-299): the code "is from BCSLIB-EXT, Boeing's Extreme Mathematical Library."

**Generalized problem.** Solves `K Φ = M Φ Λ` (Theory §36, p.629), with `K`,`M` the
assembled stiffness/mass, `Φ` mode shapes, `Λ` eigenvalues. Lanczos in its native form
works the *standard* problem `A Φ = Φ Λ` with only matrix–vector products and converges
fastest to *extreme* eigenvalues. To make the lowest structural modes "extreme," the
problem is recast in **shift-invert** form:

> `(K − σM)⁻¹ M Φ = Φ Θ`, with `θ_i = 1/(λ_i − σ)` (Theory §36, p.629).

Eigenvalues near the shift `σ` map to the largest `θ_i`, so Lanczos finds them quickly.
Each shift requires a **sparse direct factorization of `(K − σM)`** (the expensive step).

**Shift strategy + Sturm/mode count.** BCSLIB-EXT chooses a *sequence* of shifts `σ_i`
automatically. Each factorization yields the **matrix inertia** (count of negative
diagonals = number of eigenvalues left of `σ_i`) — i.e. a **Sturm-sequence count**. From
inertia at two shifts the algorithm knows exactly how many eigenvalues lie in an interval
and can certify that none were missed (Theory §36, p.629: inertia "tells the algorithm how
many eigenvalues are to the left of any given σ_i … determine if all of the eigenvalues in
that interval have been computed"). This is why it is called robust.

**I/O-bound, not CPU-bound.** Theory §36 (p.629) warns the SMP parallel speed-up is limited
because "a vast amount of data … has to be stored on I/O files"; wall-clock is governed by
the I/O subsystem as much as CPU.

**Underlying sparse direct solver** (Theory §35, pp.625-628): all 5 LS-DYNA direct solvers
are **multifrontal** (Duff & Reid 1983), with sparsity orderings **MMD** (small problems)
or **METIS** (large/3-D). Solver 4/5 = multifrontal with **distributed-memory parallelism**;
Solver 6 = BCSLIB-EXT (shared-memory). Inertia/Sturm counting rides on this factorization.

> **↔ OpenSees contrast.** OpenSees ships only a *standard symmetric* Lanczos/`-genBandArpack`
> path (ARPACK shift-invert via a banded/`UmfPack` factorization); it has **no built-in
> inertia/Sturm interval certification** that all modes in `[f_min,f_max]` were captured —
> the user must request `numEigen` and trust the count.

---

## 2. `*CONTROL_IMPLICIT_EIGENVALUE` (Vol I, PDF pp.1603-1615)

Activated with `IMFLAG=1` on `*CONTROL_IMPLICIT_GENERAL` + nonzero `NEIG` (Remark 1, p.1611).

**Card 1 (p.1604):** `NEIG, CENTER, LFLAG, LFTEND, RFLAG, RHTEND, EIGMTH, SHFSCL`.
- `NEIG` — number of modes (required). `NEIG>0` ⇒ compute at t=0 then terminate.
- `CENTER` — center frequency; finds the `NEIG` modes nearest `CENTER` (Hz; Remark 10).
- `LFLAG/LFTEND`, `RFLAG/RHTEND` — left/right **frequency-band** endpoints `[f_min,f_max]`.
  If the band has two finite endpoints and `NEIG` is set large, LS-DYNA auto-counts the
  eigenvalues in the band (via inertia) and reduces `NEIG` to that number (Remark 2, p.1612).
- `SHFSCL` — optional manual initial Lanczos shift (Remark 9, p.1614).

**`EIGMTH` — eigenvalue extraction method (p.1604-1605):**
| EIGMTH | Method | Notes |
|---|---|---|
| 2 (default) | Block Shift-Invert Lanczos (BCSLIB-EXT) | the workhorse |
| 3 | Lanczos with `[M]=[I]` | debugging only |
| 5 | as 3 + **Dynamic Terms** | |
| 6 | as 2 + **Dynamic Terms** | |
| 101 | **MCMS** (AMLS-type, SMP only, R11+) | thousands of modes for NVH |
| 102 | **LOBPCG** (preconditioned, BLR precond.) | few modes, low memory; SMP R11, MPP R14 |
| 103 | **Fast Lanczos** | thousands of *approximate* modes in **MPP** for NVH |
| 111 | **Sectoral Symmetry** | cyclic-symmetry reduction by harmonic index |

**THE COMPLEX / DAMPED EIGEN OPTION (key).** There is **no `NEIG<0` complex flag** — that's
intermittent analysis. Complex modes are triggered **on a different card**: set
**`LCPACK=3` on `*CONTROL_IMPLICIT_SOLVER`** (Remark 3, p.1612). Then:
- LS-DYNA treats the eigenproblem as **nonsymmetric** and **automatically switches the
  eigensolver to ARPACK** ("By setting LCPACK to 3, ARPACK will automatically be chosen").
- **All damping terms** in the model are added to the **first-order terms** of the
  eigenproblem ⇒ the quadratic/damped eigenproblem `(λ²M + λC + K)φ = 0`, yielding
  **complex eigenpairs** (complex frequencies + complex mode shapes, i.e. damped natural
  frequencies and modal damping ratios).
- **Implicit Rotational Dynamics** also forces ARPACK because it injects first-order
  (gyroscopic/Coriolis) terms — see Theory §36.1 (rotating systems), pp.630-631: the
  rotational body force contributes to both `C` and `K`, the inertial damping is **not
  symmetric**, so eigenvalues/vectors are **complex**; LS-DYNA forms
  `C^R = C + MΩ`, `K^R = K + K^σ + MΩ²` and solves the resulting nonsymmetric problem.
- When `LCPACK=3`, the band/center fields `CENTER, LFLAG, LFTEND, RFLAG, RHTEND, EIGMTH,
  SHFSCL` are **ignored** (Remark 3, p.1612).
- **Critical limitation:** for **eigenvalue analysis `LCPACK=3` is implemented in SMP
  only** (`*CONTROL_IMPLICIT_SOLVER` Remark 8, p.1685: "For eigenvalue analysis, LCPACK=3
  is only implemented in SMP."). MPP `LCPACK=3` works for static/dynamic implicit but
  **not** for the complex eigensolve.

**Intermittent (in-transient) eigen analysis (p.1604, p.1611).** `NEIG<0` ⇒
**intermittent eigenvalue analysis**: the curve ID `(−NEIG)` schedules eigen extractions
*during* a transient or dynamic-relaxation run, so modes track evolving geometry, stress,
contact, material state. Works with implicit (`IMFLAG=1`) or explicit (`IMFLAG=6`) transient.
Each extraction emits its own numbered `d3eigv`/`eigout` database. `EIGMSCL` (Card 4) chooses
original vs mass-scaled mass for the extraction (R15+ fix).

**Other useful flags:** `MSTRES` (compute modal stresses), `EVDUMP` (write eigenvectors —
ASCII `>0` / binary `<0`; in MPP one file per rank `Eigen_Vectors.xxxx`), `ROTSCL` (rotational
DOF inertia scaling, default 1/1000), and outputs to `d3eigv` plot DB + `eigout` summary
incl. **modal participation factors and modal effective mass tables** (Remark 11, p.1614).

> **↔ OpenSees contrast.** OpenSees `eigen` solves only the **undamped symmetric** real
> problem; it has **no complex/quadratic eigensolver** — no damped natural frequencies, no
> complex mode shapes, no gyroscopic/rotating-frame modes. This (LCPACK=3 → ARPACK
> nonsymmetric quadratic eigen) is the single most distinctive LS-DYNA capability for a
> Ladruno modal feature to target. Modal effective-mass tables are also not native.

---

## 3. Frequency domain / NVH — `*FREQUENCY_DOMAIN_*` (Vol I)

All of the modal frequency-domain procedures **require a preceding eigen extraction**
(`*CONTROL_IMPLICIT_EIGENVALUE` → `d3eigv`); they are **modal-superposition** by default.

### 3a. FRF — `*FREQUENCY_DOMAIN_FRF` (PDF pp.2763-2777)
Frequency response functions from nodal/base excitation. Explicit NOTE on the card: natural
frequencies & mode shapes are needed, so `*CONTROL_IMPLICIT_EIGENVALUE` **must** be present
(p.2763). Excitation `VAD1`: base vel/accel/disp, nodal force, pressure, enforced motion by
large-mass, torque, base angular motions (p.2765-2766). Modal damping via `DAMPF`/`LCDAM`
(ζ, mode- or frequency-dependent) or Rayleigh `DMPMAS/DMPSTF` (Card 2, p.2768). Modal
truncation via `FNMAX, MDMIN, MDMAX`. Output band `FMIN..FMAX`, `NFREQ`, linear/log/biased
spacing; amplitude-phase or real-imaginary pairs (Card 4, p.2770). Restart from `d3eigv`.

### 3b. SSD — `*FREQUENCY_DOMAIN_SSD` (PDF pp.2826-2840)
Steady-state dynamics: response to a spectrum of harmonic loads. Options BLANK/`DIRECT`/
`DIRECT_FREQUENCY_DEPENDENT`/`FRF`/`ERP`/`SUBCASE` (p.2826).
- **Default = modal superposition** (`MDMIN..MDMAX`, `FNMIN..FNMAX`, `MEMORY` in/out-of-core).
- **`DIRECT` / `DIRECT_FREQUENCY_DEPENDENT` = direct method** (solve the complex dynamic
  stiffness `(−ω²M + iωC + K)` at each frequency, no modes) — **"only available for SMP"**
  (p.2826). The standalone **`*CONTROL_IMPLICIT_SSD_DIRECT`** (PDF p.1688/12-376) requests a
  **direct complex SSD** (`ISSFLG`, `FMIN/FMAX/NFREQ`, structural `LOSS` factor) — double
  precision SMP only.
- `ERP` option = Equivalent Radiated Power (radiated sound power proxy, Rayleigh-integral or
  classic), with fluid `RO`, `C` (p.2837-2839).

### 3c. Random vibration — `*FREQUENCY_DOMAIN_RANDOM_VIBRATION` (PDF pp.2785-2801)
PSD-input → response **PSD and RMS** by **modal superposition**. `VAFLAG` loading: base
acceleration, random pressure, plane/shock/progressive/reverberant/TBL waves, nodal force,
base vel/disp, enforced motion by large mass (p.2789). `METHOD`: auto / modal superposition /
**modal acceleration** / **modal truncation augmentation** (p.2789). Damping: modal or
broadband (`DMPTYP`). Inputs are auto-PSD (`NAPSD`) and **cross-PSD** (`NCPSD`, real/imag or
mag/phase) load definitions (p.2796-2800). `FATIGUE` option computes fatigue life under random
vibration (Steinberg/Dirlik-type via S-N curves). Outputs absolute/relative PSD & RMS.

### 3d. Response spectrum — `*FREQUENCY_DOMAIN_RESPONSE_SPECTRUM` (PDF pp.2802-2806)
Peak response to an input response spectrum (e.g. earthquake) by **modal superposition** +
**mode-combination rule** `MCOMB`: SRSS(0), NRC grouping(1), **CQC(2)**, double-sum
Rosenblueth-Elorduy(3), **NRL-SUM(4)** incl. closely-spaced-modes variants(−4/−14), Gupta-
Cordero(5/6), Rosenblueth(7), ABS(8), or weighted blend of two rules(99) (p.2804).
Multi-point/multi-directional combination `MPRS`: SRSS or **100-40-40** (Newman) rule
(p.2805). `DDAM` option = U.S. Navy Dynamic Design Analysis Method (shipboard shock / UNDEX).
`MISSING_MASS_CORRECTION` option (ZPA + static unit-g file) recovers truncated high-freq mass.

> **↔ OpenSees contrast.** OpenSees has **no native frequency-domain / NVH stack**: no FRF,
> no steady-state harmonic (SSD), no random-vibration PSD→RMS, no built-in response-spectrum
> modal-combination (CQC/SRSS) — users hand-roll these in Python on top of `eigen` + modal
> participation. A modal-superposition FRF/SSD/CQC layer is the natural Ladruno target;
> direct-method SSD (complex dynamic stiffness per frequency) is a heavier second tier.

---

## 4. Parallel (MPP) eigensolve — `*CONTROL_IMPLICIT_SOLVER` (Vol I, PDF pp.1676-1687)

The shift-invert Lanczos eigensolve's cost is the repeated **sparse factorization of
`(K − σM)`**, so MPP eigen scalability = MPP linear-solver scalability. `LSOLVR` (Card 1,
p.1677):
- **`LSOLVR=7` (default) = parallel multifrontal sparse solver** — **distributed-memory
  (MPP)** factorization; this is what distributes the shift-invert solves across MPI ranks.
- `LSOLVR=2` = older parallel multifrontal (deprecated).
- `LSOLVR=30` = **MUMPS** external parallel direct/hybrid solver (CERFACS/ENS-Lyon;
  Amestoy-Duff-Koster-L'Excellent distributed multifrontal), coupled since R11.1
  (Remark 2, p.1683). Can be pure direct or, with `RPARM5>0`, a Block-Low-Rank approximate
  factorization used as a CG preconditioner.
- `LSOLVR=6` = BCSLIB-EXT direct (SMP only). `LSOLVR=22-26` = iterative CG/GMRES.

**Distributed ordering** (`ORDER`, p.1677 + Remark 6, p.1684): MMD / **METIS (serial)** /
**ParMETIS (MPP, R14+)** / **LS-GPart (MPP, Ansys, R11+)**. For large MPI-rank counts a
*parallel* ordering (ParMETIS / LS-GPart) is recommended because serial METIS becomes a
bottleneck — directly relevant to distributed eigensolves.

Note (p.1685): `LPRINT` and `ORDER` from this card **also apply to the eigensolution
software** — the eigensolve reuses the implicit linear-solver factorization/ordering stack.
BLR low-rank factorization (`IBLROPT`, `ISINGLE`, `RPARM5`) reduces factor storage/time
(refs: Amestoy et al. on BLR multifrontal, p.1687) — important because eigen factorizations
are memory- and I/O-bound (cf. Theory §36 I/O caveat).

**MPP availability by eigen method:**
- `EIGMTH=2` Block Shift-Invert Lanczos: runs in MPP via the `LSOLVR=7`/MUMPS distributed
  factorization (BCSLIB-EXT's own threading is SMP/I/O-bound).
- `EIGMTH=103` **Fast Lanczos**: the **MPP** path for thousands of approximate NVH modes
  (assumes the `LSOLVR=2` multifrontal factorization stays in memory between shifts).
- `EIGMTH=101` **MCMS**: **SMP only**. `EIGMTH=102` **LOBPCG**: SMP R11, **MPP R14+**.
- **`LCPACK=3` complex/nonsymmetric eigen: SMP only** (p.1685) — the damped/complex modal
  capability does **not** run distributed.

> **↔ OpenSees contrast.** Distributed-memory modal analysis in OpenSees (OpenSeesMP) is not
> a first-class facility: `eigen` runs on the (assembled) system and there is no distributed
> shift-invert-Lanczos riding a MUMPS/ParMETIS distributed multifrontal factorization the way
> LS-DYNA's `LSOLVR=7`/`30` + ParMETIS stack does. A Ladruno MPP modal feature would need to
> couple the eigensolve to a distributed sparse factorization (MUMPS is already a MUMPS-class
> dependency in the OpenSees ecosystem), mirroring LS-DYNA's design.

---

## Cross-cutting takeaways for the Ladruno modal feature

1. **Robustness via inertia/Sturm interval certification** (Theory §36) is the feature that
   makes "give me all modes in `[f_min,f_max]`" trustworthy — OpenSees lacks it.
2. **Complex/damped modal = `LCPACK=3` → ARPACK nonsymmetric quadratic eigen** (Vol I
   `*CONTROL_IMPLICIT_EIGENVALUE` Remark 3 + `*CONTROL_IMPLICIT_SOLVER` Remark 8). SMP-only
   even in LS-DYNA. This is the highest-value, highest-novelty target.
3. **Frequency-domain stack (FRF/SSD/random/response-spectrum)** is *all modal superposition*
   layered on the eigensolve, plus a direct (per-frequency complex-stiffness) fallback for SSD.
4. **MPP eigen = distributed multifrontal factorization** (`LSOLVR=7`/MUMPS=30 + ParMETIS),
   not a special distributed Lanczos; the parallelism lives in the linear solve and ordering.
