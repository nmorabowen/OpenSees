---
title: relativeMassDamping — mass-proportional Rayleigh damping relative to a reference velocity
project: Ladruno
tags:
  - damping
  - ssi
  - energy-balance
---

# `relativeMassDamping`

Mass-proportional Rayleigh damping that acts on the velocity **relative to a
prescribed reference velocity** instead of the absolute velocity:

    f_D,i = alphaM * M_i * (v_i - r * v_ref(t))        (nodal masses only)

Motivation: in an absolute-coordinate soil-structure model (input through the
soil boundary: absorbing boundaries, DRM) the structure's `alphaM * M * v_abs`
contains a drag `alphaM * M * v_ff` towards the fixed frame that a relative-
coordinate (fixed-base or spring) model does not have. With `v_ref` = the free-
field velocity the mass term damps the motion relative to the free field, the
damping of the relative formulation. The stiffness terms (betaK, betaK0, betaKc)
are untouched; for a structure `K * (rigid translation) = 0`, so they are
already relative.

## Command

```tcl
relativeMassDamping <tag> -vel <dof> <tsTag> [-vel <dof> <tsTag> ...] \
                    ( -region <regTag> ... | -node <nodeTag> ... | -allMassNodes )
```

| Item | Meaning |
|---|---|
| `-vel dof tsTag` | the reference velocity of dof `dof` (1-based) is TimeSeries `tsTag` evaluated at the domain time. A copy of the series is taken. Unlisted dofs (rotations) have `v_ref = 0` |
| `-region` / `-node` / `-allMassNodes` | the target set, exactly one. Only nodes that carry a **nodal** mass on this process are attached; the others are skipped |
| return | `nAttached nAlphaZero  Σ αM·m_dof1  Σ αM·m_dof2 ...` (per `-vel`, this process). The rank sums add to `alpha * structural mass`, the audit |

`alphaM` keeps coming from `rayleigh` / `region -rayleigh`; the command only
sets the frame. Issue it **after the last `loadConst -time 0.0`** (where a
compensating load pattern would sit), so gravity stages are untouched. Delays
(`v_ff` lagging the input by H/Vs) are a property of the series:
`timeSeries Path ... -startTime tau`.

Python: `ops.relativeMassDamping(tag, '-vel', 1, ts, '-region', r)`. The name is
one string, `OPS_RelativeMassDampingCmdName` (`RelativeMassDamping.cpp`).

## What changes, what does not

| Path | Change |
|---|---|
| Residual — `Node::getUnbalancedLoadIncInertia` | `-alphaM M (v - v_ref)` when a reference is attached; original statement otherwise |
| Tangent — `Node::getDamp`, `c2 * C` | **unchanged**: `v_ref` depends on time only, so `dR/dU` is the same |
| Reactions flag 2 (`rayleighForces`, Ladruno recorder `rayleighForce`, MPCO `RAYLEIGH_FORCE`), flag 0, `getResponse(RayleighForces)` (`nodalRayleighForces`) | relative force |
| EnergyBalance `DW_nod` | power `v^T C (v - v_ref)` (closure exact) |
| Everything else, and every node without a reference | byte-identical (all edits guarded by `massDampingRef == 0` with the original statement in the unguarded branch) |

## Energy bookkeeping

`v^T alphaM M (v - v_ref) = (v - v_ref)^T alphaM M (v - v_ref) + v_ref^T alphaM M (v - v_ref)`

`DW_nod` books the total (the closure `ULW = KE + IE + DW` holds as without the
feature). The total can be negative where the reference drives the structure,
so it is not a dissipation figure. `recorder EnergyBalance -v2 -relativeMassSplit`
adds two advisory columns before `RES`:

| Column | Definition |
|---|---|
| `DW_rel` | `∫ (v - v_ref)^T alphaM M (v - v_ref) dt`, non-negative, non-decreasing |
| `W_ref` | `∫ v_ref^T alphaM M (v - v_ref) dt`, the work delivered to the reference frame |

`DW_rel + W_ref = DW_nod` to round-off; neither enters `RES`. The columns are
also present, without the flag, when a reference exists at the recorder's first
record (a recorder defined before gravity records before the command is
issued, so pass the flag in that case). Trapezoidal, same first-record rule as
the accumulators.

## Integrator coverage

Supported: integrators that form the nodal unbalance through
`TransientIntegrator::formNodUnbalance` → `Node::getUnbalancedLoadIncInertia`
(Newmark, HHT, LadrunoHHT without sensitivity, the central-difference family).

Refused with a FATAL line: integrators that form the nodal damping force as
`getDamp() * Udot` through `DOF_Group::addD_Force` — the `*_TP` family and the
sensitivity branches of Newmark / LadrunoHHT / LadrunoGeneralizedAlpha /
PFEMIntegrator — because they would silently return the absolute-frame drag.
The guard exits with `_Exit(-1)` (an `exit()` from inside the analysis loop runs
static destructors and ends in a heap fault on Windows).

Explicit integrators evaluate the velocity at a half step; the reference is
evaluated at the domain time of the residual — a first-order inconsistency to
assess before using the feature with them.

## Parallel (flat-per-rank OpenSeesMP)

The reference term is `alphaM * M * (...)`: it exists only where the nodal mass
exists. A shared node whose mass is written on one rank contributes once; the
copies without mass contribute zero whatever the reference. No rank placement
is needed. Verified on 2 and 4 ranks (the shared nodes on partition borders
carry their mass on one rank) against the serial run: same deviation as the
compensating-load oracle under the same partitions (MUMPS vs BandGeneral).

## Limitations

- Nodal masses only; elements with their own mass keep `alphaM M v` through
  the element Rayleigh forces.
- Not serialized (`sendSelf`), like `alphaM` itself: re-issue on restart.
- `Node(const Node&)` does not copy the reference (OpenSeesSP partitioning).
- The reference cache is shared by the nodes of one object; the nodal unbalance
  is formed serially. A threaded nodal loop needs the cache filled before it.

## Tests

- Zone A: `tests/test_relative_mass_damping.py` (classic Tcl): SDOF against the
  closed Newmark recursion and against the compensating-load oracle; Rayleigh
  force responses; energy split; zero reference bitwise equal to no command;
  `HHT_TP` refused; audit return and refusals; openseespy registration.
- Zone B: `Ladruno_files/testbed/relative_mass_damping/` — the SSI drag toy
  (needs `drag_toy.py` from the Paper 1 repository, `DRAG_TOY_DIR`), the
  2/4-rank MPI comparison (`run_mp.bat`), the energy-split and guard script, and
  `energy_viscous_zerolength.py` (the pre-existing EnergyBalance double count of
  viscous zeroLength dashpots, see `LEDGER_quirks.md`).
