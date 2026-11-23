---
title: model testNDMaterial — the nD material tester (WP-248)
project: Ladruno
tags:
  - guide
  - material
  - tester
---

# `model testNDMaterial` — the nD material tester

A classic-Tcl model builder that drives one `nDMaterial` point by point, with no
element and no domain. It sits beside the stock `testUniaxial` and
`testPlaneStress`, and differs from them in two ways:

- **Explicit commit.** The stock testers commit on every N-th strain call and
  cannot revert. Here nothing commits until `ndTest commit`, so any number of
  trials can be made from the same committed state.
- **Mixed control.** Any component can be given as a strain or as a stress; the
  free strains are found by Newton on the material tangent.

It was written for apeGmsh's material tester (apeGmsh ADR 0135 D9), which runs
`OpenSees.exe` on a generated script. openseespy has no model-builder testers, so
the command is Tcl only.

## Commands

```tcl
model testNDMaterial                  ;# takes no count
nDMaterial ElasticIsotropic 1 200000.0 0.3

ndTest material 1 ThreeDimensional    ;# copy the material in a view; returns the order n (6 here)
ndTest strain 1e-3 0 0 0 0 0          ;# trial strain, all n components
ndTest mixed esssss 1e-3 0 0 0 0 0    ;# uniaxial stress: e11 given, the other stresses zero
ndTest getStrain                      ;# n numbers
ndTest getStress                      ;# n numbers
ndTest getTangent                     ;# n*n numbers, row by row
ndTest commit                         ;# commitState
ndTest revert                         ;# back to the last commit
```

| Subcommand | Does | Returns |
|---|---|---|
| `material $tag $view` | `getCopy($view)` of `nDMaterial $tag`, replacing any material under test | the order n |
| `order` | — | n |
| `strain e1 … en` | `setTrialStrain` | — |
| `mixed $kinds v1 … vn ?-tol t? ?-maxIter k?` | strain control on the `e` components, stress control on the `s` ones | Newton iterations |
| `getStrain`, `getStress`, `getTangent` | read the trial state | numbers, `%.17g` |
| `commit` | `commitState`, and remembers the committed strain | — |
| `revert` | `revertToLastCommit`, then sets the committed strain as the trial | — |

The view is whatever the material's `getCopy(const char*)` accepts:
`ThreeDimensional`, `PlaneStrain`, `PlaneStress`, `AxiSymmetric`, `BeamFiber`,
`BeamFiber2d`, `PlateFiber`, and so on. A material without that view is refused
by name. Components follow OpenSees' order with engineering shear strains, e.g.
3D `11 22 33 12 23 31`.

## Mixed control

`$kinds` has one letter per component: `e` means the value is a strain, `s`
means it is a stress. Newton starts from the current trial strain, imposes the
`e` components, and iterates on the `s` block of the tangent until
`max |σ_s − target_s| ≤ tol · max(|σ|∞, |target_s|∞)` (default `tol` 1e-10,
`maxIter` 25, at most 1000). A singular block (e.g. stress control on a perfectly
plastic plateau) or no convergence is a Tcl error; the material stays at its last
trial, and `ndTest revert` returns it to the commit.

Typical states:

| State | `kinds` (3D) | values |
|---|---|---|
| uniaxial stress | `esssss` | `ε11 0 0 0 0 0` |
| triaxial, confinement p | `essess` | `ε11 −p −p 0 0 0` |
| simple shear, free normals | `sssess` | `0 0 0 γ12 0 0` |
| full stress control | `ssssss` | six stresses |

## Errors

Every refusal is a Tcl error whose message starts with `ndTest:` and names the
token at fault: an unknown tag or view, a strain count other than n, a token
that is not a finite number, a `kinds` string of the wrong length or with
another letter, an unknown option, a singular block, no convergence, a material
that refuses a trial strain, commit or revert. Queries before
`ndTest material` are refused.

## Tests

`tests/tcl/wp248_nd_material_tester.tcl`, driven by
`tests/test_wp248_nd_material_tester_tcl.py` (it runs when `dist\bin\OpenSees.exe`
exists or `LADRUNO_TCL_EXE` names one). Closed forms: 3D Hooke and its tangent;
uniaxial stress (σ11 = Eε11, ε22 = ε33 = −νε11, in one Newton step); triaxial
(σ1 = Eε1 − 2νp); the plane-stress view; the J2 perfectly plastic plateau at σy;
elastic unloading at E after yield; revert; repeatable trials; the refusals.
