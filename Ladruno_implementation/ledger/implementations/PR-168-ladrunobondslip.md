---
wp: PR-168
title: "LadrunoBondSlip"
pr: "#168"
status: "shipped — CI green (Zone-A Ubuntu build + classTag/manifest gates); Zone-A batt…"
section: "table"
legacy_seq: 78
---
| **LadrunoBondSlip** — 1D bond-slip `tau`-`s` `uniaxialMaterial` (CEB-FIP Model Code 2010 backbone) for the bar-axial slot of the embedded-rebar element ([[20_ladruno_embedded_reinforcement_adr]] D4); **v2 of the RC-3D plan, first slice**. v1 = MONOTONIC backbone + elastic (k0) unload/reload, sign-symmetric. Backbone: linear `s<s0` (D4.1 initial-slip regularization — kills the power-law `dτ/ds→∞` at s→0; `s0` default `0.1·s1`, `k0=τmax(s0/s1)^α/s0`) → ascending power law `τmax(s/s1)^α` → plateau `τmax` → softening to `τf` (NEGATIVE tangent ⇒ caller uses disp/arc/IMPLEX control, D4.2) → residual `τf`. **D4.3 fracture-energy reg.:** `-Gf` overrides `s3` so the softening triangle dissipates `Gf` per unit interface area (`s3=s2+2Gf/(τmax−τf)`); the element scales `τ` by `perimeter·L_trib`. `getInitialTangent=k0` (finite). Parse `uniaxialMaterial LadrunoBondSlip tag τmax s1 s2 s3 τf α <-Gf Gf> <-s0 s0>`; setResponse `bondStress`/`slip`/`tangent`. Registered: classTags.h, uniaxial CMake, Python + Tcl command maps, broker, header-stamp glob. | Material | **MAT_TAG 33002** (Ladruno *uniaxial* band, after LadrunoUniaxialJ2=33000, LadrunoRebarBuckling=33001) | `SRC/material/uniaxial/LadrunoBondSlip.{cpp,h}`, `tests/test_ladrunoBondSlip_material.py` | shipped — CI green (Zone-A Ubuntu build + classTag/manifest gates); Zone-A battery (backbone-vs-oracle, k0-finite, sign-symmetry, Gf-overrides-s3, elastic-unload) | [#168](https://github.com/nmorabowen/OpenSees/pull/168) |
