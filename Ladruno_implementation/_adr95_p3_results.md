# ADR-95 P3 -- VERDICT (2026-09-07, build cf239c9d)

**H1 event TRANSFERS to Lagrange (tet10) and Bernstein-std, NOT (in range) to
Bernstein-bbar.** `TenNodeTetrahedron` FLOORs at s/B=0.00161, `BezierTet10 std`
at 0.00168 -- ~7x EARLIER than H20uri's 0.01114 -- with a corner GP (branch=3)
1-2 stations before the wall at the IDENTICAL location (element 534, GP 3,
x=-2.32, z=-0.25), I1 just above T=0.777 kPa, detAmin collapsing 1e6-1e8x off a
-0.02..-0.08 ambient floor (matches P1's "not the event" baseline exactly).
`BezierTet10 -bbar`: ZERO branch-2/3 GPs s/B 0.001-0.02, reaches TARGET (still
hardening, tail 6.1%) -- volumetric relief removes the corner onset entirely
over the tested range, same direction bbar/uri worked on the hex family.

Matched s/B=0.008 (N/A = walled first; DOF: hex from P1, tet 2583 nodes/7749
DOF/1200 tet10 all 3 legs, `build_mesh_tet10_prandtl.py`, 1.66x H20/5.6x H8):
h8bbar 1386 DOF q=110.09 ratio 0.7925 | h20uri 4659 DOF q=94.32 ratio 0.6790 |
beziertet10bbar 7749 DOF q=99.01 ratio **0.7128** | beziertet10 7749 DOF N/A
(walled 0.00168, ratio 0.410) | tet10 7749 DOF N/A (walled 0.00161, ratio 0.391).
`bearing_mesh_tet10.npz` (ADR-79 square-footing box) is a DIFFERENT problem,
verified unusable and NOT reused; a new plane-strain mesh was built instead.
Reading: Bezier-bbar (0.713) sits BETWEEN h20uri (0.679) and h8bbar (0.793) at
matched settlement -- deficit B is small here, not note 82's ~1.4x (different
allowance regime); the two walled legs cannot even enter this table -- their
deficit IS the early wall, not a lagging capacity.

---

## 1. What was built

- `bearing_mesh_tet10.npz` inspected first: x,y in [-10,10], z in [-8,0], 49-node
  SQUARE footing -- ADR-79's PDMY box, not this plane-strain strip. Confirmed
  unusable for the Prandtl-Reissner oracle (no q0*Nq shape-factor-free answer for
  a square footing) before writing anything else.
- `build_mesh_tet10_prandtl.py` (new): plane-strain strip, one element thick in y
  (THICK=0.5), x/z graded on `h20_prandtl.strip_mesh`'s own block boundaries
  (XLIM=30, ZBOT=-20, B_FOOT=2, R_GRADE=1.35) at h0=1.0 -- 9 graded + 2 uniform +
  9 graded in x, 7 graded + 3 uniform in z, 200 hex-shaped cells (matches the
  H20/H8 h0=1.0 legs), each split into 6 structured tets by gmsh (transfinite
  volume, no recombine). Verified: volume 600.000000 m^3 exact (rel 0), gmsh
  TET10 edge slots 4..9 match `(0,1)(1,2)(0,2)(0,3)(2,3)(1,3)` to 1e-9 (straight
  sides), which is the SAME order `BezierTet10`/`TenNodeTetrahedron` use (both
  source files say BezierTet10 was built to match TenNodeTetrahedron's order --
  confirmed by grep, not assumed) -- so gmsh connectivity needs no permutation.
  2583 nodes (462 vertex + 2121 mid-edge), 1200 tet10, 7749 DOF. Written to
  `bearing_mesh_tet10_prandtl.npz` (a NEW file; the ADR-79 npz was left untouched).
  Built with `C:\Users\nmb\venv\opensees_env\Scripts\python.exe build_mesh_tet10_prandtl.py`.
- `tet_path_diag.py` (new): the same deck as `h20_prandtl.py`/`quad_path_diag.py`
  (phi_txc=20, nu=0.45, rho_bar=0, SY=0.2, non-associated) built on the tet10
  mesh for `--elem tet10 | beziertet10 | beziertet10bbar`. Reuses
  `quad_path_diag.py`'s `sample_branch` / `sample_tangent` / `tangent_health` /
  `mobilisation` BY IMPORT (that file is untouched); the `_branch.csv` column
  layout and `_branch.npz` `stNNN_<field>` station convention are byte-identical
  to `quad_path_diag.py`'s, just under a `tpd_` prefix.

### The load-basis trap (verified against source, not assumed)

`SRC/element/bezierTetrahedron/OPS_BezierTet10.cpp`'s own docstring: Bernstein
face functions each integrate to A/6, so a uniform traction is **q*A/6 on ALL SIX
face nodes** (3 vertices + 3 mid-edges) for BezierTet10 (either formulation).
`TenNodeTetrahedron` gets the standard T6-consistent rule instead: **0 at the 3
vertices, q*A/3 at the 3 mid-edges**. `tet_path_diag.consistent_surcharge_tet`
computes each basis's own vector from the real face geometry and asserts the sum
against q0*A_top=30.0 m^2 to **1e-9** in an ELASTIC pre-step (control 3, reaction
resultant): all three legs PASS at machine precision (errors 0 to 3e-13 %). Node
ordering, edge convention, and this load split were all verified against the C++
source before any run, per the plan's explicit warning that this exact mismatch
already caused a Bezier deck to diverge at first yield (note in `OPS_BezierTet10.cpp`).

## 2. Runs (sequential, build `cf239c9d`, `ADR95_DIST` staged copy, `ladrunoBuild()`
checked at the top of every leg)

```
py -3.12 tet_path_diag.py --elem beziertet10bbar --branch --cond-at 5e-4 --cond-every 25 --sfrac 0.02 --budget 200 --tmax 2700 --suffix _p3
py -3.12 tet_path_diag.py --elem tet10           --branch --cond-at 5e-4 --cond-every 25 --sfrac 0.02 --budget 200 --tmax 2700 --suffix _p3
py -3.12 tet_path_diag.py --elem beziertet10     --branch --cond-at 5e-4 --cond-every 25 --sfrac 0.02 --budget 200 --tmax 2700 --suffix _p3
```
`--cond-every 25` (vs `quad_path_diag`'s default 4) was chosen after a smoke test
showed each dense `FullGeneral` SVD sample at 4785 free DOF costs ~20-25 s;
`--sfrac 0.02` bounds the push past H20uri's known wall (0.0111) while keeping
each leg inside the wall-clock cap. Logs: `tpd_{elem}_p3_run.log`. Outputs:
`tpd_{elem}_h1.0_p3.{csv,json}`, `tpd_{elem}_h1.0_p3_branch.{csv,npz}`.

### 2a. `beziertet10bbar` -- MODE=TARGET (did not wall)

215 steps, 0 subdivisions, 384 s. q_max = 135.10 kPa = **0.9726** of exact; end
s/B = 0.02000 of 0.02 (reached the imposed cap); tail dq/ds = 6.13 % of initial
(**STILL HARDENING**, not a plateau -- same "wall was never tested" caveat as
every quadratic leg in this campaign; 0.02 was a wall-clock choice, not the
element's true reach). Branch census at 9 stations (s/B 0.001 to 0.02, plus the
final at 0.02): **branch(0,1,2,3) never left (~4300-4700, ~350-700, 0, 0)** --
zero f2/corner GPs at EVERY sampled station. detAmin_min stayed in [-0.079,
-0.064] throughout (the ubiquitous non-associated-plasticity negative baseline
P1 already showed is not the event). At s/B=0.02: 4 of 4785 negative eigenvalues
of the symmetric tangent part -- present, small, and co-existing with a
perfectly healthy TARGET termination (H2 strong-form stays dead here too).
**H1 event: NOT OBSERVED in the tested range.**

### 2b. `tet10` (TenNodeTetrahedron, Lagrange) -- MODE=FLOOR

77 steps, 19 subdivisions, 635 s. q_max = 54.28 kPa = 0.3908 of exact; end s/B =
**0.00161** (8% of the 0.02 cap). Branch census (4 stations):

| s/B | branch(0,1,2,3) | detAmin_min | corner GP (if any) |
|---|---|---|---|
| 0.00070 | (4790,10,0,0) | -0.0219 | -- |
| 0.00122 | (4529,271,0,0) | -0.0686 | -- (last healthy station) |
| 0.00159 | (4699,100,0,**1**) | -1.53e7 | ele 534 gp 3, x=-2.319 z=-0.250, I1=0.7935 |
| 0.00161 (wall) | (4567,232,0,**1**) | -1.76e8 | same GP, I1=0.7935 |

I1 at the corner GP = 0.7935 kPa, just above T = sqrt(2/3)*SY/rho = **0.7771**
kPa (rho from the same `h20_prandtl.alpha_from_phi_txc`/`rho=sqrt(2)*alpha`
formula P1 used). detAmin at that GP is -1.5e7 to -1.8e8; every OTHER GP at the
same two stations stays in [-0.063, -0.056] -- an **8-9 order-of-magnitude**
outlier, the same signature as H20's "elsewhere -0.06" baseline. The last
healthy station (0.00122) had zero corner GPs; the corner GP appears one
sampled station later and the leg FLOORs 2e-5 of s/B after that. 5 of 4785
negative tangent eigenvalues at the wall.

### 2c. `beziertet10` (BezierTet10 std, Bernstein, no bbar) -- MODE=FLOOR

34 steps, 10 subdivisions, 530 s. q_max = 56.94 kPa = 0.4099 of exact; end s/B =
**0.00168**. Only 2 branch stations landed before the wall (coarse cadence, fast
FLOOR): s/B=0.00100 -- (4642,158,0,0), detAmin_min -0.0392, zero corner GPs;
s/B=0.00168 (wall) -- (4548,251,0,**1**), detAmin_min **-1.43e6**, corner GP at
the SAME location as the tet10 leg (**ele 534, gp 3**, x=-2.319, z=-0.250),
I1=1.4525 kPa (further above T=0.7771 than the Lagrange leg's 0.7935 -- std has
no volumetric relief, so the mean-stress rise into the cutoff is faster and
overshoots further before the tangent is sampled). Non-corner detAmin at the
wall station stays at -0.057. 4 of 4785 negative tangent eigenvalues.

## 3. Reading

The corner-branch consistent-tangent defect P1 named on H20 is not an H20
artefact: it fires at the IDENTICAL mesh location (element 534) under both a
Lagrange (TenNodeTetrahedron) and a Bernstein-std (BezierTet10) quadratic basis,
with the same qualitative signature (I1 just clears T, detAmin collapses 1e6-1e8x
against a -0.02..-0.08 ambient floor common to every plastic GP under non-associated
flow). What moves is WHEN: on this coarser tet mesh (h0=1.0, 3 elements/B less
resolved near the footing edge than H20's own grading) both non-relieved tet legs
wall at s/B ~ 0.0016-0.0017, roughly 7x earlier than H20uri's 0.0111 -- consistent
with a coarser/stiffer discretisation reaching the same mean-stress trigger sooner,
not with a different mechanism. `-bbar` (the tet's volumetric relief, playing the
role `uri`/`bbar` play on the hex family) removes the corner GP entirely over the
tested range and reaches its imposed target still hardening -- the same qualitative
protection bbar/uri buy on H20/H8, transferred to a completely different basis and
element family. This is independent, out-of-family confirmation of P1's H1 finding
and of P4's fix-lane framing (the defect is in the DruckerPrager corner tangent,
not in any one element's geometry).

## 4. What was NOT done

- `beziertet10` was run ("if time permits"); it FLOORed as fast as `tet10`, so no
  leg in this campaign reached s/B=0.008 except the bbar one -- the "matched s/B"
  table above has two N/A rows by necessity, not by omission.
- No attempt to push `beziertet10bbar` past s/B=0.02 to find ITS wall (if any) --
  out of scope for P3 (that is P4/P5 territory); reported as an allowance
  (TARGET at a chosen cap), never as a capacity (tail 6.1%, not a plateau).
- Bit-identical repeat (rule in force for a number entering a conclusion) was not
  run for any tet leg in this pass -- flagged for the owner review, same as P1's
  own repeat was done in a follow-up rather than the first pass.
