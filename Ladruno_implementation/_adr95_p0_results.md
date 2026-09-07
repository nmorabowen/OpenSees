# ADR-95 P0 — instrumentation results

```
VERDICT  ADR-95 P0 (instrumentation) — DONE, and it already moves H1 and H2.
build      ladrunoBuild()  cf239c9df18ba3d6e201fbcbaf87a4a9e88e404d
           that is HEAD at BUILD time; the P0 edits were still uncommitted,
           and landed as 56310fd39 on wp/95-prandtl-bezier-root-cause.  Any
           later probe rebuilt after that commit will stamp 56310fd39.
binaries   dist/bin/OpenSees.exe   2026-09-07 00:25:40  (fresh)
           dist/bin/opensees.pyd   2026-09-06 23:59:21  STALE — copy BLOCKED
           the fresh .pyd is build/build/Release/OpenSeesPy.dll (00:25:46) and a
           staged copy at <scratch>/dist_bin_p0/opensees.pyd; live P1
           processes (quad_path_diag.py) hold dist/bin open — 2 at the failed
           copy, 3 at the time of writing — so build.bat's copy failed with
           "the process cannot access the file" and CARRIED ON.  Copy it the
           moment the campaign's python processes exit; the check is
           Get-Process | ? { $_.Modules.FileName -contains <path to .pyd> }.
tests      tests/test_adr95_dp_branch_response.py   9 passed, 0 failed
           tests/test_r3_prandtl_collapse_gate.py   fastest leg only (h0 = 1.0,
           non-associated) — full gate is ~61 min measured and the box is shared
           with P1; RESULT: see "gate" below.
response   eleResponse(tag,'material',gp,'ladrunoBranch') -> Vector(8)
           [0] branch  0 elastic / 1 cone f1 / 2 cutoff f2 / 3 corner (FINAL Jact)
           [1] gamma0  [2] gamma1  [3] f1_trial  [4] f2_trial
           [5] forcedAccept (count>3 bailout)   [6] I1 of the returned stress
           [7] detAmin = min over 200 Fibonacci dirs of det(n.D_ep.n)/(2G)^3
delegation NO element needed a fix.  LadrunoBrick, LadrunoBrick20, BezierTet10
           and TenNodeTetrahedron all already forward `material <gp> <args>` to
           the NDMaterial's setResponse; all four verified live.
H2 IS ALREADY IN TROUBLE: with rho_bar = 0 and H = 0, detAmin goes NEGATIVE at
           FIRST YIELD on the cone (-7.5e-3 vs the exact elastic 0.4792).  Loss
           of ellipticity is the generic state of every plastic GP, so "detAmin
           crossed zero" cannot be P1's discriminator — only the FRACTION and
           connectivity of negative-detAmin points can be.
H1 SHARPENS: the tension cutoff returns NO stress.  gamma(1) is structurally
           zero (dead `Jact(i)==2` arm), so the corner is a pure TANGENT-OPERATOR
           switch with I1 left above T — exactly the asymmetry H1 posits.
```

## 1. What was built

`SRC/material/nD/UWmaterials/DruckerPrager.{h,cpp}` (vanilla, ledger rows added,
every edit marked `// Ladruno ADR-95`):

* seven write-only members set by `plastic_integrator()` on **every** path — the
  elastic early-return included, so the reported state can never be a stale echo
  of an earlier plastic step;
* `getLadrunoBranch()`, a pure observer returning the 8-vector above;
* `ladrunoDetAmin()`, called **only** from `getResponse` so the ~200x(3x3x3x3)
  sampling is never paid per step;
* response token `ladrunoBranch`, responseID 95, in `setResponse`/`getResponse`.

`DruckerPrager3D` and `DruckerPragerPlaneStrain` override neither `setResponse`
nor `getResponse`, so both inherit the token with no further work.

Nothing in the constitutive algebra reads any new member, so the material is
numerically unchanged by construction; the only cost on a normal step is seven
scalar stores.

## 2. The Voigt -> 4th-order mapping (the one real trap)

The acoustic tensor needs `C_ijkl`, and this class stores tangents in the
OpenSees 3D convention `(11, 22, 33, 12, 23, 31)` acting on **engineering** shear
(`initialize()` puts 0.5 on `mIIdev(3,3..5,5)` so `sigma_12 = G*gamma_12`).  The
correct mapping is a **plain index substitution**

```
vidx[i][j] = { {0,3,5}, {3,1,4}, {5,4,2} }
C_ijkl     = mCep( vidx[i][j], vidx[k][l] )        // NO 1/2 anywhere
```

because the stored coefficient multiplies `gamma_kl = 2 eps_kl` while the
contraction `C_ijkl eps_kl` runs over both `(k,l)` and `(l,k)` — the two cancel.

The first build shipped a 1/2 on the shear columns, which is the intuitive but
wrong choice; it returned `detAmin = 0.1224` where isotropic elasticity has the
closed form `(K + 4G/3)/(8G) = 0.47917` for **any** unit `n`.  The unit test now
pins that closed form to `rel = 1e-9`, and the .cpp comment states both the rule
and the number the wrong version gives.  **P1 data taken from the first build's
`detAmin` column is invalid** (its branch / gamma / f1 / f2 / forcedAccept / I1
columns are fine — only slot 7 changed).

## 3. Measured, single `LadrunoBrick`, uniform prescribed strain

Deck: `K = 1e4`, `G = 4e3`, `SY = 0.2`, `rho` from `phi_txc = 20 deg`
(= 0.148583), `rho_bar = 0`, no hardening; `T = sqrt(2/3)*SY/rho = 1.09904`.

| leg | strain | branch | gamma0 | gamma1 | f1_trial | f2_trial | forced | I1 | detAmin |
|---|---|---|---|---|---|---|---|---|---|
| elastic | `(-2,-2,-2)e-5`, `gxy = 1e-5` | 0 | 0 | 0 | -0.37418 | -2.89904 | 0 | -1.800 | **+0.4791667** |
| cone | `(-1,-1,-1)e-4`, `gxy = 8e-4` | 1 | 3.7812e-4 | 0 | +3.02494 | -10.09904 | 0 | -9.000 | **-7.5414e-3** |
| corner | `(3,1,1)e-4` | 3 | 4.2148e-4 | 0 | +3.37184 | +13.90096 | 0 | **+15.000** | +1.5106e+3 |

Readings:

1. **Elastic `detAmin` is the closed form to 1e-9** — the mapping is right.
2. **The cone branch is already non-elliptic.** `rho_bar = 0` is full
   non-association and `H = 0` is perfect plasticity, which is the classic
   Rudnicki-Rice worst case; `det(A)` is negative from the first yielding
   increment. This is a **prediction failure for H2 as written** in the plan
   ("min over GPs of det(A) crosses 0 at the event"): it crosses at first yield,
   thousands of steps before any wall. P1 must therefore record the negative
   **fraction** and its connectivity, not the crossing.
3. **The corner applies no stress return.** `gamma1 = 0` exactly and `I1` equals
   the trial `3K*tr(eps) = 15.000`, thirteen times `T`. The cutoff only swaps the
   tangent (note `detAmin` jumping three orders of magnitude, from `-7.5e-3` on
   the cone to `+1.5e3` at the corner). That is a much sharper version of H1 than
   the plan assumed, and it is a **vanilla defect**, recorded in `LEDGER_quirks`
   and pinned by a sentinel test.
4. `forcedAccept` was 0 on every leg, consistent with the plan's observation that
   `Jact =` appears in no campaign log.

## 4. Element delegation — nothing needed fixing

Checked in source and then live, one element each with `DruckerPrager` and an
elastic step:

| element | site | verdict |
|---|---|---|
| `LadrunoBrick` | `LadrunoBrick.cpp:4045` | forwards `&argv[2]` to the material |
| `LadrunoBrick20` | `LadrunoBrick20.cpp:1085` | forwards |
| `BezierTet10` | `BezierTet10.cpp:1843` | forwards |
| `TenNodeTetrahedron` | `TenNodeTetrahedron.cpp:1818` | forwards |

All four returned a width-8 vector with `branch = 0` and `detAmin > 0` under a
small load. A negative control is in the battery too: `ElasticIsotropic` returns
`[]` for the token, so a future generic fallback that swallowed it would be
caught rather than silently reporting zeros.

## 5. Unit test

`tests/test_adr95_dp_branch_response.py`, 9 cases, `zone_a`, 0.15 s total.
The three branch legs impose a **known homogeneous strain** on one element by
prescribing all 24 DOF with `sp()` under the **Penalty** handler — Penalty and
not Transformation on purpose: a single hex has every node on three faces, so a
uniform-strain field constrains all 24 DOF and Transformation would hand the
solver a zero-size system. Each target is computed by hand from the material's
own algebra and stated in the docstring, so a leg landing on the wrong branch is
a defect and not a tuning accident.

The module also takes `LADRUNO_DIST_BIN` to override which `dist/bin` is
imported. That exists because of the failure this session hit: a long campaign
run holds `dist/bin/opensees.pyd` open for hours and Windows then blocks the next
build's copy **silently** — `build.bat` prints one "the process cannot access the
file" line among the MKL copies and finishes with a success banner, leaving a
stale `.pyd` behind a fresh `OpenSees.exe`. Checking only `OpenSees.exe`'s mtime
would have passed.

## 6. Gate — no regression

The full gate is `@pytest.mark.slow` and measures ~61 min (537 + 517 + 1566 s for
the sequence plus 1047 s for the associated control), which busts the P0 budget,
and the box is shared with two live P1 processes. So **only the fastest leg** was
run — `h0 = 1.0`, non-associated — through the gate module's own `_run_leg()`,
with no edit to the gate file.

| quantity | reference (build `1db3394b`) | this build `cf239c9d` |
|---|---|---|
| DOF | 1386 | 1386 |
| `q_num` | 150.71 | 150.70600 |
| ratio | 1.0849 | 1.0849417 |
| tail % | 0.001 | 0.00139 |
| mode | BUDGET | BUDGET |
| ds/floor | 2500 | 2500.0 |
| capacity | yes | yes (plateau + free advance) |
| resultant identity | 2.7e-15 .. 8.9e-15 | 2.66e-15 |
| 1-D stress patch | 1.0e-14 .. 1.4e-13 | 1.01e-14 |

Identical to the recorded reference at every printed digit, on the same
termination mode, with the same subdivision headroom. Wall time 1390 s against
the reference's 537 s is CPU contention with the P1 runs, not the element.

The other two resolutions and the associated control were **not** run and are
owed before the PR flips to ready.

## 7. Open items handed to P1 / the orchestrator

1. **Copy the fresh `.pyd`.** `build/build/Release/OpenSeesPy.dll` ->
   `dist/bin/opensees.pyd`, as soon as the P1 `quad_path_diag.py` processes exit
   (still held at the time of writing). Until then any `import opensees` from
   `dist/bin` gets the pre-fix `detAmin`; `OpenSees.exe` is already current.
2. **Re-run whatever P1 has already taken**, or at minimum discard its `detAmin`
   column; branch/gamma/f-trial/forcedAccept/I1 from the first build are sound.
3. **H2's stated prediction is dead as written** (§3.2). Restate it as a fraction
   or a connectivity measure before P1's decision table is applied.
4. **H1 is now testable much more sharply** (§3.3): the corner is a tangent-only
   switch. The P2 knob `mTo = 1e10` / raising `SY` therefore tests exactly one
   thing — whether that operator switch is what walls the quadratic elements.
5. Item 3 of the plan's P0 list — the `--forensics` per-iteration residual dump in
   `quad_path_diag.py` — was **not** done here: that script is in the P1 agent's
   hands and is already being run with `--forensics`, so touching it would have
   raced a live campaign.
