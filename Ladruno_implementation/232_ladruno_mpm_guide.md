---
title: "WP-232 — the explicit MPM engine (mode F): user guide"
project: Ladruno
tags: [mpm, wp-232, guide]
---

# The explicit MPM engine (`mpm`) — user guide

Slice 1 of WP-232 ([[232_ladruno_mpm_general_driver_adr]]): a native explicit Material Point
Method engine, hosted as a `LoadPattern` (`LadrunoMPMPattern`, class tag 33027) and driven by
the `mpm` command in classic Tcl (`OpenSees.exe`) and Python (`opensees.pyd`). The numerics are
frozen by `tests/mpm_interface.md` (WP-232b); this guide says how to drive them.

## A minimal run

```tcl
model basic -ndm 2 -ndf 2                      ;# ndm 2 = plane strain (thickness 1), or 3
nDMaterial ElasticIsotropic 1 1.0e7 0.3
mpm grid 1 -origin 0 0 -h 0.1 -cells 50 10     ;# the host is load pattern 1
mpm body 1 box 0.5 0.0 1.5 0.4 -ppc 2 -mat 1 -rho 1000
mpm wall ymin 0.3                              ;# stationary Coulomb wall on the grid face
mpm gravity 5.0 -8.66
timeSeries Path 2 -time {0 0.2 1e9} -values {0 1 1}
mpm timeSeries 2                               ;# body force = factor(t) * gravity
set dtc [mpm criticalTimeStep]
mpm analyze 5000 2.0e-4
puts [mpm particleData x]
```

Python is the same with `ops.mpm('grid', 1, '-origin', 0.0, 0.0, ...)`.

## Commands

| Command | What it does |
|---|---|
| `mpm grid <tag> -origin <x0..> -h <h> -cells <n..> [-basis gimp\|linear]` | creates the host (load pattern `<tag>`); one per Domain in v1. `gimp` (uGIMP) is the production basis; `linear` is the cell-crossing contrast |
| `mpm body <id> box <lo..> <hi..> -ppc <n> -mat <tag> -rho <ρ>` | lattice seeding (interface §5): `n` particles per cell per axis at `x0 + (a+½)h/n`, `V0 = (h/n)^ndm`, `lp = h/(2n)`; also `disc cx cy R` (2D) and `ball cx cy cz R` (3D). Ids run 0, 1, 2, … in seeding order. Returns the count |
| `mpm shift <body> <dx..>` | translate a seeded body, before the first step |
| `mpm particle <id> <body> <mat> <ρ> <V0> -x .. -lp .. [-v ..] [-s ..]` | one explicit particle (`-s`: initial stress, Voigt, tensor shear) |
| `mpm velocity <body> uniform <v..>` / `sine <comp> <axis> <v0> <β> <x_ref>` | initial velocity of a body |
| `mpm stress <body> uniform <s..>` / `geostatic <y_top> <ρ> <g> <K0>` | initial stress of a body (2D geostatic: σyy = −ρg(y_top − y), σxx = K0 σyy) |
| `mpm fix <face> <comps>` | zero the listed grid-velocity components on a face (`xmin … zmax`; `comps` ⊂ `xyz`). A roller is the face's normal component |
| `mpm wall <face> <μ>` | stationary rigid wall with Coulomb friction μ (no penetration, friction once per step) |
| `mpm gravity <g..>` | body force per unit mass, signed |
| `mpm timeSeries <tsTag>` | the factor applied to gravity (default 1) |
| `mpm damping local <α>` | local non-viscous (FLAC) damping, α ∈ [0, 1) — for quasi-static runs |
| `mpm flip <α_F>` | FLIP fraction: 1 pure FLIP (default), 0 pure PIC |
| `mpm massCutoff <c_m>` / `mpm safety <s>` | node active iff `m_i > c_m · max m_p` (default 1e-10); Δt refused iff `Δt > s · Δt_c` (default 0.5) |
| `mpm criticalTimeStep` | `Δt_c = h / max_p (c_p + |v_p|)`, `c_p = sqrt(M_p/ρ_p)` from the current tangent |
| `mpm analyze <n> <Δt>` | n explicit steps; recorders run at each `Domain::commit`, then flush |
| `mpm reset` | particles and materials back to the state at the first `mpm analyze`; counts an epoch (`mpm info epoch`) |
| `mpm info numParticles\|epoch\|time\|summary` | `summary` = [mass, momentum.., max \|v\|, kinetic energy] |
| `mpm particleData <field>` | flat list in storage order: `id body mat m V0 V x v lp F R sig sigHat epsHat sig0` |
| `mpm debug shape <x..> -lp <lp..>` / `shape1d <h> <lp> <r>` / `activeNodes` | unit-level hooks (interface Q4) |

## One step

1. refuse if a particle centre left the real grid (interface A2; the error names the lowest such id) or `Δt > s·Δt_c`;
2. `Domain::applyLoad(t+Δt)`: the domain time advances and the host reads its series factor;
3. P2G of mass, momentum and forces over the particles' support nodes; nodes below the cutoff drop;
4. damping, grid momentum update, grid conditions (walls, then fixes);
5. MUSL: particle velocity (PIC/FLIP blend) and position from the grid, a second projection,
   ΔF and the incremental hypo-kernel entry (Hughes–Winget, unrotated frame), the material
   update — with no active element (`ops_TheActiveElement = 0`);
6. commit the particle materials (a declared refusal aborts the run);
7. particles take their trial state; 8. `Domain::commit()` (recorders); 9. profiler rows.

The engine grid carries one ghost node layer beyond every face (interface A5), so a particle's
GIMP support is never cut: the partition of unity holds right up to the face, and a body slides
along a face without drag. A face's `fix` or `wall` acts on its face nodes and its ghost layer
alike; a face with no condition just has free ghost nodes. Node ids count the ghost layer
(interface §2): `id = (i+1) + (nx+3)·((j+1) + (ny+3)·(k+1))`, with the ghost layer at index −1.

## Materials

| Material | v1 |
|---|---|
| `ElasticIsotropic` | admitted; reproduces the interface's hypoelastic law, with initial stress as an offset |
| `LadrunoJ2` | admitted; `-autoRegularization` refused (it reads an element length; there is none at particles) |
| `DruckerPrager` (UW) | admitted; stage switch `materialState` through `parameter` (below) |
| anything else | refused at `mpm body` / `mpm particle`, naming the class tag |

**Initial stress** (`mpm stress`, `-s`) is accepted on `ElasticIsotropic` only. On LadrunoJ2 or
DruckerPrager it is refused: an additive offset would bypass the yield check, and v1 has no
yield-consistent seeding (interface Q2, 232c's design).

**Stages.** `updateMaterialStage` walks Domain elements only and never reaches particles
(ADR W12). Use the host as the parameter target — the spelling is `loadPattern` in both engines:

```tcl
parameter 1 loadPattern 1 materialState
updateParameter 1 0           ;# DruckerPrager elastic
```

## What `reset` does

`reset` (the Domain's `revertToStart`) reaches a load pattern only as a re-`applyLoad`, so the
particles **stay where they are**. `mpm reset` puts them back. For a clean restart use both.
`wipe` and `remove loadPattern <tag>` delete the host and its particles.

## Refused in v1

- a Domain with a free FE DOF (no FE–MPM seam yet): fix every node or remove the FE model;
- database `save` / `restore` while the host exists (no particle restart path yet);
- `mpm` in OpenSeesSP / OpenSeesMP;
- adding particles after the first step (until `mpm reset`).

## Output

Stock recorders cannot see particles. Use `mpm particleData`. The `.ladruno` particle group
(Contract 2, `FORMAT_MINOR_VERSION` 2) arrives with phase c2.

## Verification

- `tests/test_mpm_engine_core_cpp.py` — the pure core (`LadrunoMPMCore.h`) and the kernel entry,
  standalone with g++, against every case of `tests/mpm_cases.txt` and
  `tests/hypo_incremental_cases.txt` (232b's independent numpy oracle).
- `tests/test_ladruno_mpm_engine.py` — every case through `OpenSees.exe` **and** `opensees.pyd`,
  and the two engines against each other. The MPM mutation family's suite.
- `tests/test_ladruno_mpm_lifecycle.py` — refusals, lifecycle, admitted materials, the stage path,
  thread-count determinism.
- Known: the 20%-compression self-weight column misses its frozen 0.03 (0.072, as in the oracle;
  uGIMP cell crossing) — a strict xfail per interface A4.
