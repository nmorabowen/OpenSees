---
title: ADR-97 P2 — principal-space closest point for MohrCoulomb / MCTC (report)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P2 (`wp/97c-cp-principal`) — implementation report

**PR** [#824](https://github.com/nmorabowen/OpenSees/pull/824) (draft, based on
`wp/97b-cp-smooth` = P1 [#819](https://github.com/nmorabowen/OpenSees/pull/819),
itself based on `wp/97a-plan-oracles` [#817](https://github.com/nmorabowen/OpenSees/pull/817))
· **plan** [[97_ladruno_asdp_closest_point_adr]] · **P1** [[reviews/adr97_p1_report]]
· **oracle** `Ladruno_implementation/adr97_oracle/cppm_mc.py`

## 1. What shipped

`integration_method Closest_Point` + `tangent_type Algorithmic` for the
**Mohr-Coulomb family**, through a **principal-stress-space multi-surface
closest-point return** (Clausen, Damkilde & Andersen 2006/2007) rather than P1's
smooth 6D Newton.

The shipped 6D Mohr-Coulomb yield function is the exact Lode-angle form
`f = A(theta) sqrt(J2) + I1 sin(phi)/3 - c cos(phi)`, whose gradient carries a
`1/cos(3 theta)` that the header dodges with a `|theta| >= 29 deg`
Drucker-Prager substitution and, by default, a central difference over the six
raw Voigt slots. A coupled Newton on that gradient cannot converge quadratically
at a corner — which is where Mohr-Coulomb models live. In principal space the
same surface is six **planes**: on the sorted sextant `s1 >= s2 >= s3` (tension
positive) one plane, two corner **lines**, one vertex. Every return is then a
closed-form linear projection in the elastic metric — **no Newton anywhere** —
and the Koiter tangent is constant on each region.

| piece | how |
|---|---|
| region selection | Clausen's **boundary planes** through the apex, spanned by each edge direction and the face return direction `D m_13`. Never a `dLambda >= 0` active-set search: at the apex with `psi < phi` the three Koiter multipliers are **not** all positive (ADR-97 P0 header finding 6). |
| boundary-plane signs | taken **analytically** — any point strictly inside the open face is `apex + b1 ell_1 + b2 ell_2` with `b1, b2 > 0` and `n_1 . ell_1 == 0`, so `sgn_1 = sign(n_1 . ell_2)` and symmetrically. This replaces the oracle's dimensional `apex - 20.0` reference point (an ADR-94 M5 unit trap). Verified identical to the oracle's calibration over **60 (phi, psi, c) combinations and 13094 trial states, 0 mismatches**. |
| edge directions | `null(A) = a_i x a_j` (cross product) instead of the oracle's SVD, oriented by the same `ell(0) >= ell(2)` rule. Agrees with the SVD to **2.2e-16**. |
| back-transform | `C = Rs T Rs^-1 E` with `Rs` the Voigt image of `E -> Q E Q^T` (`Q` the TRIAL eigenvectors), `T`'s normal block `dy/dx` and `T`'s shear slots the eigenprojection **rotation term** `(y_i - y_j)/(x_i - x_j)`, with the l'Hôpital limit `dy_i/dx_i - dy_i/dx_j` on a degenerate trial eigenvalue and a threshold **relative** to `strength_scale()` (ADR-94 M5). `Rs^-1` is built as the Voigt image of `Q^T . Q`, not inverted numerically (agrees with `inv(Rs)` to 2.4e-15). |
| plastic strain | `d eps^p = E^-1 (sigma_tr - sigma_ret)` — convention safe (it never touches the engineering-shear factor of the principal flow vectors) and exact on every region including the vertex. |
| MohrCoulombTensionCutoff | the trial is offered **first** to ADR-84's `special_return` hook (cutoff face / Rankine edge / MC∩TC corner / compound corner / apex), whose raw Koiter `stiffness_return` maps to `Algorithmic` and whose `SR_QUALITY_FALLBACK` is refused under `strict_convergence` — the same contract `Backward_Euler` has. **No ADR-84 geometry is re-derived.** When the hook declines (cutoff inactive at the trial, or an MC-dominant Stage-3c trial), the plain-MC principal return takes over and the **composite** `f` is re-checked. |

`Backward_Euler` and every YF/PF path it executes are untouched (ADR-97 D1).

**Build.** `1072c27ae` — the last commit that changes anything under `SRC/`
before the mutation gate. `ops.ladrunoBuild()` on `dist/bin/opensees.pyd`
reports `1072c27ae7048d776da8154ef66c8bbdee84b36d`.

## 2. Support count — 22 of 46, not the plan's 31

The plan wrote "support rises to 31 (YF {VM,DP,MC} x PF {VM,DP,MC} = 30, + MCTC)".
That assumed every cross combination is a registered, verifiable specialization.
Reading `ASD_material_definitions.cpp`: `MohrCoulomb_YF` appears in **7**
specializations and `MohrCoulomb_PF` in **6**, and in only **one** of them are
both of the Mohr-Coulomb family.

The principal-space map is valid only when the yield surface **and** the flow
potential are both the piecewise-linear Mohr-Coulomb ones. A mixed pairing has
no verified map at all: the principal return assumes linearity of both, and P1's
smooth 6D map cannot use Mohr-Coulomb's Lode-angle gradient. So P2 ships:

| specialization | status |
|---|---|
| `MohrCoulomb_YF<Null>` x `MohrCoulomb_PF<Null>` | **supported** |
| `MohrCoulombTensionCutoff_YF<Null>` x `..._PF<Null>` | **supported** |
| `MohrCoulomb_YF` x `VonMises_PF` (Null / TensorLinear / AF) | refused |
| `MohrCoulomb_YF` x `DruckerPrager_PF` (2) | refused |
| `MohrCoulomb_YF` x `HoekBrown_PF` | refused |
| `VonMises_YF` x `MohrCoulomb_PF` (2), `DruckerPrager_YF` x `MohrCoulomb_PF` (2), `HoekBrown_YF` x `MohrCoulomb_PF` | refused |

**20 (P1) + 2 = 22 of 46.** The parse-time refusal names the mixed pairing
explicitly and says why. Enforced at compile time by
`ladruno_cp_principal_family`, which is non-zero only when the YF's and PF's
family markers are **equal and non-zero** *and* every internal variable is inert
(the map has no `q`-row, so a live hardening law would be silently ignored).

## 3. Gate-by-gate results

### Gate 1a — the regions, against the P0 oracle (`tests/test_adr97_p2_principal.py`)

One step of prescribed total strain `eps = E^-1 sigma_tr`, so the elastic
predictor IS the oracle's trial state. `E = 30000, nu = 0.25, phi = 30,
psi = 10, c = 10` — `cppm_mc.py`'s own constants.

| trial | region | rel. error vs the oracle | committed \|f\| | `cp_iterations` |
|---|---|---|---|---|
| face, no shear | face | **1.5e-16** | 5.3e-15 | 1 |
| face, sheared (rotation term live) | face | **9.1e-16** | 5.3e-15 | 1 |
| edge `s1 == s2` | line1 | **0.0** | 1.8e-15 | 1 |
| edge `s2 == s3` | line2 | **8.0e-16** | 5.9e-14 | 1 |
| apex, hydrostatic tension | apex | **0.0** | 0.0 | 1 |
| apex, slightly deviatoric | apex | **0.0** | 0.0 | 1 |

Pinned at 1e-10 relative, measured at the 1e-16 floor. That pin is tighter than
the fork's usual cross-platform 1e-6 and it is justified: the return is
closed-form linear algebra whose only platform-variable step is a 3x3 symmetric
eigen-decomposition, and every one of these trials is either non-degenerate or
returns a state that is **invariant under the eigenvector ambiguity** (the
degenerate pair comes back equal). Everything that passes through a Newton
(`Backward_Euler` comparisons) is pinned at 1e-6.

`cp_iterations` reads **1** on every plastic step and **0** on an elastic one:
on this path it is a region count, not an iteration count.

### Gate 1b — admissible at every commit

| path | steps | plastic | worst committed `f` | `cp_iterations` |
|---|---|---|---|---|
| triaxial compression | 10 | 5 | **+1.2e-14** | {0, 1} |
| simple shear | 10 | 8 | **+1.8e-15** | {0, 1} |
| rotating principal directions (2 legs) | 20 | 16 | **+2.0e-14** | {0, 1} |

The rotating path is checked to actually rotate: the first eigenvector of the
committed stress moves **20.70 deg** between the legs.

### Gate 1c — MohrCoulombTensionCutoff, against ADR-84's own oracles

| deck | result |
|---|---|
| hydrostatic tension, 3x the cap | `sigma = 24.7 * I` **exactly** (cap error 0.0); CP vs BE gap **0.0**, i.e. bit identical |
| uniaxial tension, 4x the yield strain | `sigma_zz = 24.7` exactly (the Rankine face return); CP vs BE gap **0.0** |
| confined compression (cutoff inactive ⇒ fall-through to the plain-MC return) | worst `f_MC` **+2.8e-14** on a strength scale of 93.97; worst `f_TC` **-217** (the cutoff is satisfied with room to spare) |

The two zero gaps are the load-bearing result: `Closest_Point` calls ADR-84's
hook rather than re-deriving its geometry, so wherever the hook governs the two
integrators commit the same doubles.

### Gate 2 — the consistent tangent

`Algorithmic` against a central difference of the **binary's own assembled
internal force**. Two rigs, because no single one reaches both regions:

| rig | region | rel_err | notes |
|---|---|---|---|
| `fd_check` oedometric (`rig="uniaxial"`, `nu = 0.15`) | edge, **degenerate eigenvalue** | **2.88e-11** (`rel_fro` 2.15e-11) | `s1 == s2` exactly ⇒ the l'Hôpital branch |
| free-node, face axis aligned | face, separation 0.190 | **8.40e-09** | rotation term live |
| free-node, face sheared | face, separation 0.271 | **1.08e-08** | rotation term live |

The **free-node rig** is new in this WP and local to the test file: the
homogeneous strain field of the P1 driver is prescribed on seven of the eight
nodes and node 7 is left free. Seven prescribed nodes mean there is no limit
point at any stress level, so any state can be reached — in particular a genuine
face state with separated principal stresses, which neither of
`fd_tangent_driver`'s rigs can produce (both of them land on `s1 == s2`).

**Negative control.** Not a number, an outcome: on the same free-node face rig
`Backward_Euler`/`Continuum` and `Backward_Euler`/`Secant` both fail to converge
(`analyze -> -3`), where `Closest_Point`/`Algorithmic` converges and reproduces
the assembled tangent to 8.4e-9. ADR-94 M3 measured those two operators at 57 %
and 80 % against a central difference of the material's own response; at that
error a Newton on a nearly singular perfectly plastic element simply stops.

**Iteration contrast** (confined cube, load driven, 6 steps, `ops.testIter()`):

| deck | `Closest_Point`/`Algorithmic` | `Backward_Euler`/`Continuum` | `Backward_Euler`/`Secant` |
|---|---|---|---|
| MC oedometric (`nu = 0.15`) | 2, 2, 3, 3, 3, 3 = **16** | 2, 2, 13, 15, 15, 15 = 62 (**3.9x**) | 2, 2, 10, 12, 12, 12 = 50 (3.1x) |
| ADR-84 MCTC (kPa) | 3, 3, 4, 4, 4, 4 = **22** | 10, 13, 13, 13, 12, 13 = 74 (**3.4x**) | — |

### Gate 4 — `Backward_Euler` inertness, and a finding

`tests/test_adr97_p4_inertness.py` re-run on this binary: **23 decks / 282
committed-stress rows BYTE-IDENTICAL** against the pre-ADR-97 baseline
`3622d6214`, each in a fresh subprocess. The MohrCoulomb, MCTC and
MC-from-the-MCTC-file decks are among them, so the four headers P2 edits are
covered directly.

`Closest_Point` is **step-size independent**, as a closed-form projection onto a
fixed surface from a proportional path must be:

| N | rel. error vs the oracle |
|---|---|
| 1 | 1.5e-16 |
| 4 | 8.9e-16 |
| 10 | 4.4e-16 |
| 40 | 5.2e-16 |

And at the apex `Closest_Point` **equals** `Backward_Euler` to **2.1e-16** — the
one region where the shipped gradient cannot be wrong, because the integrator
does not use it.

**The finding.** Mohr-Coulomb's flow direction is CONSTANT inside a sextant (the
surface is piecewise linear), so on a deck whose iterates stay in one sextant the
Ortiz-Simo cutting plane and the closest point are the SAME point and the two
maps must agree. Measured, they do not — unless `MC_ds > 0`:

| `MC_ds` | branch taken by `MohrCoulomb_YF::df_dsigma_ij` / `MohrCoulomb_PF` | `Backward_Euler` vs the closest point |
|---|---|---|
| `0` (what every deck in this repo passes) | **analytic** Lode-angle `c1/c2/c3` | **2.87e-01** |
| `1e-8` | central difference of its own `f` / `g` | 1.03e-10 |
| `1e-6` | ″ | 1.29e-12 |
| `1e-4` | ″ | **2.11e-14** |

`f` is exactly linear in principal stress, so a central difference of it over the
Voigt slots is the EXACT 6D gradient — and its accuracy *improves* with a larger
step, which is the signature of differencing a linear function. Two independent
references (this closed-form return, and the header's own difference of its own
`f`) therefore agree with each other and against the shipped analytic
coefficients. `Backward_Euler` through the analytic branch is 4.4e-2 (N = 1) to
2.9e-1 (N = 10) away from the exact return and is step-size dependent.

**Recorded, not fixed.** Fixing it changes `Backward_Euler`, which D1 keeps
byte-identical. It is pinned in BOTH directions by
`test_gate4_backward_euler_agrees_only_through_its_own_finite_difference`, so it
cannot rot silently: the test turns red the day either half changes. Every
Mohr-Coulomb deck in `tests/` passes `MC_ds 0.0`, i.e. runs on the analytic
branch; the ADR-84 MCTC battery is largely insulated because `special_return`
does not use that gradient at all.

### Gate 6 — fail loud (12 tests in the P2 file, 13 in the P1 file)

* the **six mixed pairings** the generator registers are still refused at parse
  time (`MohrCoulomb_YF` x {`VonMises_PF` Null, `VonMises_PF` TensorLinear,
  `DruckerPrager_PF`, `HoekBrown_PF`}, {`VonMises_YF`, `DruckerPrager_YF`} x
  `MohrCoulomb_PF`);
* `HoekBrown` (P3) still refused; `MohrCoulomb` and `MCTC` now **accepted**
  (the positive half of the matrix);
* `tangent_type Algorithmic` still refused with `Backward_Euler` (D2 survives P2);
* `MC_phi == 0` (Tresca: the apex is at infinity and Clausen's boundary planes,
  which pass through it, are undefined) is **refused loudly** — `analyze -> -3`
  on a `LadrunoBrick`, which checks the refusal sentinel;
* the **degenerate-eigenvalue** path: exactly hydrostatic, hydrostatic + 1e-12
  deviator, hydrostatic + 1e-7 deviator all commit the finite admissible vertex
  `17.32050808 * I` — no NaN. This is the class of state ADR-94 B4 used to
  commit NaN through;
* `strict_convergence 1` is **byte-inert** on a converging MC deck (rel gap
  exactly 0.0).

`tests/test_adr97_p6_failloud.py::test_closest_point_is_refused_for_unconverted_families`
asserted that MohrCoulomb was refused. That half is now **inverted**, said so in
its own docstring, and replaced by the mixed-pairing and HoekBrown assertions.

### Gate 5 — mutation

See [[reviews/adr97_p2_mutation]].

## 4. Found while implementing

Four things worth the next agent's time; all four are in [[LEDGER_quirks]].

1. **The oedometric deck at `nu = 0.25` with `phi = 30` never yields.** The
   elastic `K0 = nu/(1-nu) = 1/3` coincides EXACTLY with the Mohr-Coulomb
   compression meridian `(1-sin phi)/(1+sin phi) = 1/3`, so `f` is identically
   `-c cos(phi)` at every stress level. A tangent gate built on that rig is
   silently vacuous — it measures the ELASTIC tangent and passes. Every rig in
   the P2 file now asserts the state is plastic before measuring anything.
2. **`yf_tolerance()` is not a usable admissibility tolerance for a stress
   reassembled from a spectral decomposition.** It defaults to the ABSOLUTE
   `f_absolute_tol = 1e-6` (`f_relative_tol` defaults to 0). On the ADR-84 MCTC
   deck (kPa, `|sigma| ~ 5.4e3`) the header's own `f` recomputed from
   `Q diag(y) Q^T` came back at 3.4e-6 — 6e-10 RELATIVE, i.e. round-off — and
   the first version of this map refused the step. The round-off is amplified
   because an EDGE return lands exactly on a corner where the Lode angle is ill
   conditioned (`dtheta/dJ3 ~ 1/cos(3 theta)`, and `dA/dtheta` is not stationary
   there). The fix is two checks: an EXACT one in principal space against all
   three surfaces of the sextant (perfectly conditioned; this is the one that
   catches a coding error) and the header's composite `f` at a tolerance
   relative to `max(strength_scale, |sigma_ret|)`.
3. **The shipped analytic Lode-angle gradient is wrong** (gate 4 above).
4. **A load-driven rig cannot reach a Mohr-Coulomb FACE state.** Both of
   `fd_tangent_driver`'s rigs put `s1 == s2` (the oedometric one by construction,
   the 12-DOF one because a lateral load takes it past its limit point first),
   so a face-region tangent measurement needs a kinematically over-determined
   rig — hence the free-node rig above.

## 5. Open

* The analytic `c1/c2/c3` Lode-angle branch of `MohrCoulomb_YF` /
  `MohrCoulomb_PF` (and their MCTC twins) — a separate PR, because it changes
  `Backward_Euler`. It should be opened with the gate-4 measurement above as its
  warrant and the P0 oracle as its reference.
* `MC_phi == 0` (Tresca) is refused rather than supported. The boundary planes
  can be re-anchored on a point of each edge LINE instead of the apex, which
  generalises; not done here because no registered deck uses it and a loud
  refusal is honest.
* The **apex tangent is rank 0** (the vertex does not move), so an element every
  one of whose Gauss points sits at the apex is singular. That is the true state
  of affairs — the same statement P1 makes for the Drucker-Prager apex — but it
  means the free-node FD rig cannot measure the apex tangent, and the P2 gates
  do not try.
* Mixed YF/PF pairings (six specializations) and `HoekBrown` (P3), `StiffSoil`
  (P5).
* `Closest_Point` is still **not** the default and `Algorithmic` is **not** the
  default tangent (D1; the flip is gated on P6).
