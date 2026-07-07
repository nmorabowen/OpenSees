---
title: "LadrunoComplexEigen — complex / state-space modal analysis (user guide)"
project: Ladruno
type: guide
related:
  - "[[46_ladruno_complex_modal_adr]]"
  - "[[45_ladruno_modal_family_roadmap_adr]]"
  - "[[LEDGER_quirks]]"
updated: 2026-07-07
---

# `complexEigen` — complex modal analysis for non-classically damped models

Real modes answer "at what frequency does the structure vibrate?"; they cannot
answer **"what is the damping ratio of the isolation mode?"** — the question
every base-isolation / supplemental-damper / SSI review asks. `complexEigen`
answers it: true per-mode damping ratios ζ_k, damped frequencies ω_d,k, and
**phased (complex) mode shapes** ψ_k for models whose damping matrix is not
proportional (localized dashpots, bearings, radiation damping).

## Quick start (openseespy)

```python
ops.eigen(12)                      # ordinary real eigen first (any solver)
res = ops.complexEigen()           # flat list, 7 numbers per physical mode:
# [omega0, omegaD, zeta, Re(lambda), Im(lambda), kind, resid] * nModes
#   kind: 0 underdamped (one entry per conjugate pair), 1 overdamped, 2 rigid
#   resid: ||(l^2 M + l C + K) z|| — per-mode quality metric
```

Tcl: `eigen 12; set res [complexEigen]` — same layout.

Options: `complexEigen [numModes] [-numModes p] [-tol eps] [-closedForm]`.

## The two projection routes

- **DEFAULT (assembled, "Route B")** — projects the model's **actual** M and C
  element-by-element from `getDamp()`/`getMass()` plus nodal mass/`alphaM`
  terms. This is *exactly the C your transient analysis feels*: scoped
  `region -rayleigh`, per-element factors, `betaKinit`/`betaKcomm`, material
  dampers (`Viscous` in `zeroLength`/`twoNodeLink`), and each element's
  `-doRayleigh` switch are all honored.
- **`-closedForm` ("Route A")** — the diagonal Rayleigh closed form
  `C̃ = diag(alphaM + betaK*w^2)` from the last **global** `rayleigh` call.
  Fast oracle; refuses `betaKinit`/`betaKcomm` and warns when scoped Rayleigh
  is present (it cannot see it).

## Complex mode shapes

`complexEigen` pushes ψ_k = Φ z_k onto the nodes. Record them per node/DOF:

```python
ops.recorder("Node", "-file", "re1.out", "-node", 1, 2, "-dof", 1, 2,
             "complexEigenRe1")     # real part of mode 1
ops.recorder("Node", "-file", "im1.out", "-node", 1, 2, "-dof", 1, 2,
             "complexEigenIm1")     # imaginary part of mode 1
ops.complexEigen(); ops.record()
```

Magnitude/phase are post-processing: `mag = hypot(re, im)`,
`phase = atan2(im, re)`. Phase convention: each mode is normalized so its
largest-magnitude entry is real-positive — relative phases are meaningful and
reproducible. For classical (proportional) damping the shapes collapse to the
real eigenvectors (Im ≈ 0): a nonzero phase spread is the fingerprint of
genuinely non-classical damping and shows *where* energy is dissipated.

## Contract & traps (read before trusting ζ)

1. **C = whatever `getDamp()` returns.** Damping that does not flow through
   `getDamp()` is invisible: `modalDamping` (integrator-side — the command
   WARNS), element `damping`-framework objects, HHT-α numerical damping, and
   any element whose `-doRayleigh` flag is left at its default **0** —
   `Truss` and the `zeroLength` family default OFF ([[LEDGER_quirks]]).
2. **Re-run `eigen` after any model edit.** M̃/C̃ are read live; the spectrum
   is a snapshot — and the reported ω₀ echoes the snapshot, masking staleness.
3. **Completeness is yours.** The projection spans only the retained p real
   modes; retain enough to cover the band of interest.
4. **Serial only** (parallel re-hosting = ADR 43 complex contours). Requires a
   constraint handler that distributes eigenvectors to slave DOFs — the
   default Transformation handler does; `constraints Plain` + MP does not.
5. `eigen` on tiny models: use `ops.eigen("-fullGenLapack", n)` (ARPACK needs
   numModes < nDOF).

## Validation

`tests/test_ladrunoComplexEigen.py` (21): independent companion-form numpy
oracle @1e-10; analytic SDOF branches; 2-DOF non-classical exact @1e-8;
`-doRayleigh` quirk pins; equalDOF slave damper; base-isolated stick gate —
complexEigen ζ_iso/ω_d vs exact state-space AND vs a Newmark free-decay
log-decrement fit (the G-A program gate, [[45_ladruno_modal_family_roadmap_adr]] §6).
