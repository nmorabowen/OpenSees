---
title: "WP-159 — Smoothed (C1) mortar contact law: -smoothN g0, -smoothT r"
project: Ladruno
type: ADR (amends the ADR-41 C2/C3 mortar normal law and return map; status row in the ADR-48 capstone)
status: "PR #903 open (not merged) — opt-in, default byte-identical; the R3 pile blocker R3-N1 is NOT closed by it (section 5)"
owner: nmora
related:
  - "[[41_ladruno_mortar_alm_contact_adr]] (the C2 normal law and the C3 return map this smooths)"
  - "[[155_pile_contact_r05]] (-augment never, -adjust, -gapOffset, -maxGap; the definitions stream v4)"
  - "[[157_mortar_friction_pair_state]] (per-pair friction state; the onset weight is per pair)"
  - "[[158_mortar_tangent_diagnosis_consistanttan]] (open item: 'residual kink at lift-off'; -consistanttan)"
  - "[[48_ladruno_contact_capstone_adr]] (status-of-record row)"
  - "apeGmsh piles-validation ladder/R3_interface/VERDICT.md (R3-N1), ladder/R0_8_fork_slice/BUILDER.md"
tags: [adr, contact, mortar, friction, smoothing, regularization, pile, wp-159]
updated: 2026-10-02
---

# WP-159 — Smoothed (C1) mortar contact law

> [!summary] The short version
> Two opt-in options on a 3D `-mortar` contact make the contact residual C1:
>
> | Flag | Smooths | Law |
> |---|---|---|
> | `-smoothN <g0>` | the normal kink at gap 0, and the jump of the friction traction at lift-off | `P = 0` (x ≤ −S), `(x+S)²/(4S)` (abs(x) < S), `x` (x ≥ S), with `x = −εN·ḡ`, `S = εN·g0`; the friction traction is scaled by the C1 weight `χ = t²(3−2t)`, `t = (x+S)/(2S)` clamped to [0, 1] |
> | `-smoothT <r>` | the stick/slip corner of the return map | `min(ρ, cap)` becomes `ρ − (ρ−cap+δ)²/(4δ)` over abs(ρ−cap) < δ = r·cap |
>
> The consistent tangent is FD-checked by a numpy oracle (4000 states, 7.6e-9) and on the binary
> (`printA` against a central FD of `printB`: every state to the FD noise, ≤ 4.5e-6). Off, the build is
> byte-identical: the 61-file contact battery gives 382/382 on both binaries and all 145,862 post-`analyze` hashes match.
> **On the R3 pile deck the law helps but does not unblock R3**: the bonded lateral push gets from
> 0.83 mm to 6.3 mm, but α-cohesion lateral and axial fail in S1 or the first step. The measured
> reason is §5: once friction is engaged from the reference, the *shipped* law fails α S1 the same
> way. R3's α S1 passes on the shipped law only because the engagement origin latches at the first
> Newton iterate. The remaining blocker is the mortar Tresca slip on the faceted interface, not the
> normal kink.

## 1. Question

R3-N1 (piles-validation `ladder/R3_interface/VERDICT.md`): implicit Newton cannot carry the 3D mortar
interface through active-set changes. Lift-off behind a laterally loaded pile and a moving slip front
both stall or diverge. All the unbalance sits at contact nodes. It persists with `-consistanttan`,
every εT, matched polygons, line search and Krylov, and small steps. ADR-158 left it open as "a
residual kink at lift-off". The shipped normal law `p = min(0, λ + εN·ḡ)` has a kink at ḡ = 0, and
a cohesive cone `μN + c` makes the friction traction **jump** from up to c to 0 when a node opens.
Does a C1 law let Newton converge through these changes?

## 2. Decisions

### D1 — the normal law: a quadratic onset over a gap band ±g0 (`-smoothN g0`)

Per slave node of a pair, with `x = −(λ + εN·ḡ)` and `S = εN·g0`:
`P = 0` for x ≤ −S, `(x+S)²/(4S)` for abs(x) < S, `x` for x ≥ S; `p = −P`. The pressure is C1 with a
Lipschitz slope `P′ = (x+S)/(2S)`, never tensile, exactly the shipped ramp outside the band, and
within `S/4` of it inside (the maximum, at x = 0). The tangent replaces the 0/1 active mask by `P′`:
`K = εN·P′·b bᵀ/a ⊗ n⊗n`.

**Rejected: softplus** `p = εN·s·ln(1+e^(−g/s))`. Its pressure is positive at every gap, so every
paired node would be "in contact" for friction (engaged, with a tangential bond) at any distance.
The quadratic onset has compact support. **Rejected: a one-sided band** (onset on the closed side,
zero pressure at ḡ = 0). Under `-adjust` every node starts at ḡ = 0, where the one-sided law has
zero stiffness: the pile floats in the first iterate (the ADR-155 G-9 problem, which its closed-branch
trick cannot fix for a C1 law, whose slope at the onset is exactly 0). The centred band starts with
half the penalty stiffness. **The price**: under `-adjust` the reference state carries `P = S/4`.
Pick `g0` so that `εN·g0/4` is small against the working pressures.

### D2 — the friction onset: a C1 weight over the same band (part of `-smoothN`)

The friction traction of the return map (unchanged inside) is multiplied by
`χ = t²(3−2t)`, `t = (x+S)/(2S)` clamped to [0, 1]. Outside the band (x ≥ S) χ = 1: the shipped
friction. The tangent is `χ·K_ss(dN/dz = εN·P′) + εN·χ′·tF ⊗ n` (the second term is non-symmetric
and is assembled only under `-consistanttan`, the ADR-158 D2 rule).

**The recipe rule this forces: `g0 ≥ c/εN` for a cohesive interface.** Releasing a bond c over the
band is a local softening of about `0.75·c/g0` (shear-normal). When that exceeds εN the
force-controlled path has a limit point and no Newton variant follows it (oracle G3, pure cohesion at
g0 = 1e-5 with c/εN = 2e-5: a snap). For the R3 α bands (c = 7–32 kPa, εN = 1e7) that is g0 ≥ 3e-6;
for the artificial "bond" (c = 1e4) it is g0 ≥ 1e-3, where the S/4 prestress is 2.5 MPa.

**Rejected (measured, then reverted): spread the onset over the pressure range [0, max(S, c)].** It
bounds the softening for any c, but it is not the shipped law in the limit g0 → 0 for a cohesive
interface (closed nodes with P < c lose part of their cohesion), and on R3 it made α S1 worse
(a 2-cycle at 4.6e3 kN).

### D3 — the rounded stick/slip corner (`-smoothT r`, 0 < r < 1)

`LadrunoFrictionKernel::frictionReturnMapSmooth` / `frictionTangentBlockSmooth`, added beside the
shipped functions (untouched). The traction magnitude `min(ρ, cap)` becomes the C1 blend
`ρ − (ρ−cap+δ)²/(4δ)` over abs(ρ−cap) < δ = r·cap, with the plastic slip that keeps
`T = kt(gT − gpT_trial)` (so the trial state stays a pure function of committed state).
`dφ/dcap` is the *total* derivative (δ moves with cap); the first oracle missed that and the FD gate
caught it (3.6e-2 → 7.6e-9). It needs friction (refused without `-mu`/`-cohesion`).

### D4 — interactions, refusals

- **ALM: refused.** Both options need `-augment never`. With the Uzawa update λ ← −P(x), a node whose
  pressure is below εN·g0 converges to an *open* gap held by a positive pressure, so the smoothing
  would not vanish as the multiplier converges. Defining that is not worth it; the pile recipe is
  already `-augment never`. The setter refuses `commit` and `request` by name.
- `-adjust`, `-gapOffset`: compose (the gap shift is applied in `mortarActive`, before the law).
  `-gapOffset +g0` moves the reference to the band edge (zero pressure, but also zero stiffness; D1).
- `-maxGap`: unaffected (handle-time pairing).
- ADR-157 per-pair state: χ and the rounded corner act per (slave node, facet pair). The
  engagement origin `gT0` is still captured at the first evaluation with P > 0, now at the outer
  band edge.
- Refused by name: NTS (no `-mortar`), `-tie`, `-soft` (SOFT=2 explicit) and `-visc` (both keep the
  unsmoothed active set), `g0 ≤ 0`, `r ∉ (0, 1)`; a 2D pair draws a handle-time FATAL.
- **Definitions stream v4 → v5** (+2 mortar slots, tail-appended). A v4 stream (every ADR-155..158
  database) still reads, as smoothN = smoothT = 0.

## 3. Verification

**Oracle** (`contact_prototypes/proto_adr159_smooth_normal.py`, exit 0):

| Gate | Result |
|---|---|
| G1 law: C1 at the edges (normal and corner), P ≥ 0, ramp outside the band, sup abs(P − ramp) = S/4 | PASS |
| G2 one-node tangent vs central FD, 4000 random states (open / band / closed, stick / blend / slip, cohesion / Coulomb / Tresca), with and without `-smoothT` | worst 7.6e-9 |
| G3 crease node (two facet pairs, ±20°) pushed to full slip with lift-off, force control | Coulomb and mixed cones: ≤ 5 iterations per step at g0 = 1e-6…1e-4; pure cohesion snaps at g0 = 1e-5 < c/εN (D2); shipped 2–3 |
| G4 g0 → 0 on a closed case | error 3.7, 0.54, 4.9e-5, 1e-15 for g0 = 1e-4, 3e-5, 1e-5, 3e-6 (penetration 1e-5) |

A one-node toy does **not** reproduce the shipped stall: the stall needs the mesh coupling of the
pile deck.

**Binary FD** (scratch probe, one slave quad on springs over a master quad, `printA` against the
central FD of `printB` at `setNodeDisp -commit` states, friction engaged at zero slip first):
stick, uniform and differential slip, one node open, in band, slip plus band: every case ≤ 4.5e-6
with and without `-smoothN`/`-smoothT`. Near the cap the shipped tangent is off by 32 % (the kink);
with `-smoothT 0.2` it is 4.1e-6.

**Tests** (`tests/test_adr159_mortar_smooth_contact.py`, 18 cases; every behavioural case fails on
the base binary, which refuses the flags):

- (law) a stiff block on springs: in the band and beyond it the equilibrium matches the smoothed law
  to 1e-6, beyond the band it equals the shipped law, pulled off past the band the master reaction
  is exactly 0;
- (a) a cohesive block with a lateral pull, pressed then pulled off past separation (force control,
  12 steps): ≤ 8 iterations per step (3, 3, 3, 3, 4, 4, 6, 2, 1, 1, 1, 1), zero reaction once open.
  A rigid block flips its four nodes at once, so the shipped law converges here too (1–2);
- (b) the ADR-158 solid roof (μ = 0.2, `-consistanttan`) with a band far inside the penetration:
  the analytic descent to 2e-4 in ≤ 6 iterations per step;
- (c) g0 → 0 on a closed case converges monotonically to the shipped state, and equals it once the
  band no longer reaches the penetration;
- (db) save → wipe → restore reproduces a smoothed contact exactly;
- (ref) the nine refusals of D4, and the recipe is accepted.

**Battery.** The 61-file contact battery (the ADR-157 list plus the ADR-158 test), MKL_NUM_THREADS=1,
with the byte-dump plugin (`contact_prototypes/bytedump_plugin.py`): **382/382 on `fc75db7f3` and
382/382 on this build; all 145,862 post-`analyze` snapshot hashes (282 tests) identical.** No flag
given => `setMortarSmoothing` is never called and the 3D mortar code path is the shipped one.

**Linux.** Esmeralda, gcc 11.4, the CI-like OpenSeesPy build (conan, bundled BLAS/LAPACK, MKL
disabled): `test_adr159` + `test_adr157` + `test_adr158`, 33 passed.

## 4. R3 pile deck (P20, M1, εN = 1e7, `-adjust -augment never -maxGap 0.1 -consistanttan`, Pardiso, MKL_NUM_THREADS=1)

| Run | Flags | S1 | Lateral / axial |
|---|---|---|---|
| bond, shipped (`M1_b_bond`, R3) | — | strict, 11 its | 0.83 mm after a cut, then fails |
| bond | `-smoothN 1e-5 -smoothT 0.1` | loose (strict stalls at 0.13 kN) | 0.83, 1.67, 2.5 mm in 5, 8, 6 its; 3.3, 4.2, 5.0 mm with cuts; fails at 6.3 mm (diverges) |
| bond | `-smoothN 1e-6 -smoothT 0.1` | strict, 6 its | 0.83 mm (7 its), fails at 1.67 mm |
| bond, skin rotated 3.75° + neighbour masters | `-smoothN 1e-5 -smoothT 0.1` | strict 6 its, then a quadratic zero pass (2.6e-3 → 1.1e-5) | 0.83, 1.67 mm (10, 7 its); fails at 2.1 mm in a fixed 5-cycle. Shipped: fails at 0.42 mm |
| α lateral | `-smoothN 1e-5 -smoothT 0.1` | loose | fails in the first step, norms growing. Shipped: fails at 0.42 mm |
| α axial | `-smoothN 1e-5 -smoothT 0.1` | loose | diverges in the first 1 mm step. Shipped: 1–3 mm in 7–11 its, fails at 4 mm |
| α S1 only, **shipped law** | `-gapOffset -1e-6` (engaged from the reference) | **fails** (norms growing, 340 kN) | — |

Lateral did **not** reach 0.1 D in any configuration.

## 5. What R3 says about R3-N1

1. **The normal kink is part of it, not all of it.** The bonded lateral push goes 7.6× further, and
   the rotated-skin S1 becomes quadratic. The failures that remain are fixed cycles (a 5-cycle at
   9.3e1 … 9.8e0 kN) and slow growth, not the lift-off chatter.
2. **The α lanes are blocked by mortar Tresca slip, under either law.** The shipped α S1 passes only
   because `gT0` latches at the first Newton iterate, which forgives the gravity settlement slip.
   Engage the nodes from the reference (a 1 µm interference) and the shipped law fails α S1 with
   the same signature as the smoothed one. The smoothed law engages every node at the reference
   (P = S/4 > 0), so it inherits that harder problem.
3. **The bond release is a softening, whatever the smoothing.** A bond c released over a band of
   width w costs about c/w of shear-normal stiffness. D2's rule `g0 ≥ c/εN` keeps it below εN.

## 6. Follow-ups

- The mortar Tresca slip on the faceted interface (the α S1 / axial blocker): FD-probe the R3 S1
  tangent at the reference-engaged state, the ADR-158 way, before adding any term.
- The `gT0` latch inside a Newton step is history-dependent: a node's stick origin depends on which
  iterate first touched. That is a separate, deeper issue than the law (all three lanes have it).
- Geometric clip non-smoothness (ADR-158 D3) remains on a mismatched polygon; the rotated-skin twin
  separates it from the law.
- NTS is not smoothed (its kernel is separate; refused by name).
