---
wp: LEGACY
title: "Symmetric-storage solvers silently DROP one u-p coupling block on LadrunoUP models — ProfileSPD (the no-system-command DEFAULT) returns rc = 0 with garbage por…"
legacy_seq: 170
---
## Symmetric-storage solvers silently DROP one u-p coupling block on LadrunoUP models — ProfileSPD (the no-`system`-command DEFAULT) returns rc = 0 with garbage pore pressures

**Symptom** (ADR-71 P1, 2026-07-11): a LadrunoUP transient run without a
`system` command (or with `system ProfileSPD` / any symmetric-storage SOE)
"succeeds" — `analyze()` returns 0 every step, displacements look plausible —
but the pore-pressure field is wrong by ~87 orders of magnitude. Measured on
the identical Terzaghi 1×10 Q4 column script: `system UmfPack` gives
p ∈ [0, 5] kPa (matches the Terzaghi series); `system ProfileSPD` gives
p ~ 1e88–1e89 with every `analyze()` still returning 0. Nothing fails loudly.

**Cause**: the honest-p contract (ADR-71 §3.2 ⟨FW-F1⟩) makes the effective
transient tangent UNSYMMETRIC — −Q lives in `getTangentStiff()` (u-rows) and
+Qᵀ in `getDamp()` (p-rows), so c₁K + c₂C never has matching off-diagonal
pairs. Symmetric-profile assembly stores only upper-triangle-in-profile
entries: one Q block is silently discarded at `addA()` time, and the solver
then factors and solves the mutilated (still well-conditioned-looking) system
cleanly. No framework hook lets an element reject an SOE.

**Rule**: every LadrunoUP model MUST name a general solver — `system UmfPack`
/ `SuperLU` / `FullGeneral` / `BandGeneral` (serial), MUMPS with `SYM=0` in
the MPI targets. The parser prints a one-line notice at element creation; the
Zone-B battery pins the divergence
(`tests/test_ladruno_up_element_analytic.py::test_wrong_solver_divergence_profilespd_vs_umfpack`).
