# `LadrunoModalResponse` — modal-superposition response toolkit (ADR 44)

Commands: `modalResponseHistory` (P1a, transient) · `responseSpectrumAnalysis
-combine` (P1b, CQC/SRSS) · `frequencyResponse` / `steadyStateDynamics` (P2, FRF/SSD).

> **Availability.** All three commands work in **openseespy** *and* in the classic
> Tcl exes (`OpenSees.exe` / `OpenSeesSP.exe` / `OpenSeesMP.exe`) — wired into
> `SRC/tcl/commands.cpp` dispatch. In Tcl the `{f, …}` result table is returned as
> a Tcl list *and* written to `-out <file>` (identical values); e.g.
> `set frf [frequencyResponse -freq 0.1 4.0 20 -baseAccel -dir 1 -damp 0.04 -node 3 -dof 1 -out frf.out]`.
> Smoke test: `Ladruno_scripts/verify_modal_classic_tcl.tcl`.

## `modalResponseHistory` — exact modal-superposition transient (ADR 44 P1a)

`LadrunoModalResponse` (classTag 33024) computes a **linear** time-history response
by exact modal superposition instead of a direct Newmark/HHT run. After you have a
mode basis (`eigen`) and the modal properties (`modalProperties`), each retained
mode becomes a decoupled SDOF oscillator that is advanced by a closed-form,
**exact-for-piecewise-linear-load** recurrence — no iteration, no factorization,
no per-step linear solve. The physical response is reconstructed as the modal sum
and **one domain step is committed per time station**, so your ordinary recorders
(node disp/vel/accel, element forces, reactions) capture the history exactly as in
a direct run.

This is dramatically cheaper than direct integration for genuinely linear models
(fast linear seismic, equipment, floor vibration) because the cost rides the
one-time eigensolve you already paid for.

> **Linear only.** Superposition is exact only while the model stays linear. Any
> material or geometric nonlinearity invalidates it — use a direct integrator then.

## Command

```
modalResponseHistory -dt $dt -nsteps $n -baseAccel $tsTag -dir $dir
    (-damp $xi | -rayleigh $a0 $a1 | -modalDamp $xi1 $xi2 ...)
    [-modes $m1 $m2 ...] [-t0 $t0]
```

openseespy:

```python
ops.eigen('-fullGenLapack', nModes)      # or the default solver for large N
ops.modalProperties()
ops.recorder('Node', '-file', 'roof.out', '-time', '-precision', 15,
             '-node', roof, '-dof', 1, 'disp')
ops.modalResponseHistory('-dt', dt, '-nsteps', n,
                         '-rayleigh', a0, a1,
                         '-baseAccel', tsTag, '-dir', 1)
```

| flag | meaning |
|---|---|
| `-dt $dt` | time step (must equal the recurrence step; sampled from `-baseAccel`) |
| `-nsteps $n` | number of steps; the run commits `n+1` stations (`t0 … t0+n·dt`) |
| `-baseAccel $tsTag` | `timeSeries` tag giving the ground acceleration `ü_g(t)` |
| `-dir $dir` | global excitation direction (1…ndf), matching `modalProperties` |
| `-damp $xi` | one modal damping ratio applied to every mode |
| `-rayleigh $a0 $a1` | Rayleigh factors → per-mode `ξ_a = a0/(2ω_a) + a1·ω_a/2` |
| `-modalDamp $xi1 …` | explicit per-mode ratios (absolute mode order) |
| `-modes $m1 …` | 1-based subset of modes to include (default: all extracted) |
| `-t0 $t0` | start time; base accel is sampled at `t0 + k·dt` (default 0) |

## What it computes

For a classically-damped linear system `M ü + C u̇ + K u = −M R ü_g(t)`, substituting
`u = Σ_a ψ_a q_a` (with `ψ_a` the eigenvector scaled exactly as `responseSpectrum`
uses it) decouples into `p` SDOFs

```
q_a'' + d_a q_a' + ω_a² q_a = −Γ_a ü_g(t),   d_a = c_a/m_a,  Γ_a = participation factor.
```

Each is advanced by a constant per-mode recurrence exact for a load that varies
linearly within the step (Nigam–Jennings / first-order hold), unified across the
under-, critically-, and over-damped branches (and the rigid `ω=0` case) by the
discriminant `Δ = ω² − (d/2)²`. Response is recovered `u = Σ_a ψ_a q_a` and likewise
`u̇`, `ü`.

## Gotchas

- **`-fullGenLapack` for small models.** The default ARPACK eigen solver needs
  `NEV < N`; use `ops.eigen('-fullGenLapack', N)` for tiny verification models.
- **Recorder precision.** Add `-precision 15` if you compare against an analytic
  reference — the default text output is ~6 significant figures.
- **Ground-motion length.** Make the `-baseAccel` record extend at least one sample
  past `t0 + nsteps·dt`; a `Path` series sampled at exactly its final abscissa can
  return 0 for the last station (floating-point round-off).
- **Do not validate against a direct Newmark run using *stiffness-proportional*
  (`betaK`) Rayleigh damping** — OpenSees's element-level `betaK` on some elements
  (e.g. Truss) does not reproduce the assembled `a1·K` and differs by a few percent.
  Mass-proportional (`alphaM`) Rayleigh and `modalDamping` agree with the exact
  modal solution.

## Response-spectrum combination — `responseSpectrumAnalysis -combine` (P1b)

OpenSees's `responseSpectrumAnalysis` computes per-mode modal displacements and
commits one step per mode, **leaving the SRSS/CQC/… combination to you**. P1b adds an
opt-in combination stage that writes a single combined nodal-displacement field:

```python
ops.eigen('-fullGenLapack', nModes)
ops.modalProperties()
# CQC (closely-spaced modes) — needs modal damping:
ops.responseSpectrumAnalysis(dir, '-Tn', *Tn, '-Sa', *Sa, '-combine', 'CQC', '-damp', 0.05)
# or SRSS / ABS / TenPercent (no damping needed for these):
ops.responseSpectrumAnalysis(dir, '-Tn', *Tn, '-Sa', *Sa, '-combine', 'SRSS')
# read the combined design displacements:
u = ops.nodeDisp(node, dof)
```

| rule | formula | use |
|---|---|---|
| `ABS` | `Σ|R_a|` | most conservative upper bound |
| `SRSS` | `√(Σ R_a²)` | well-separated modes |
| `TenPercent` | SRSS + `2Σ|R_iR_j|` over pairs within 10% period (Reg. Guide 1.92) | some closely-spaced |
| `CQC` | `√(ΣΣ ρ_ij R_iR_j)`, Der Kiureghian ρ_ij | closely-spaced modes (needs `-damp`/`-modalDamp`) |

- `-combine` absent ⇒ **byte-identical** to the stock per-mode behavior.
- `-combine` and `-mode` are mutually exclusive (combination needs all modes).
- **Combination is per-quantity and nonlinear**: this combines the nodal
  **displacement** field. Combined element forces/drifts are *not* obtained by
  deriving them from the combined displacements — combine those quantities' own
  per-mode peaks instead (record per-mode with `-mode k`, then combine).
- Known stock quirk: the `-scale` factor is ignored by the per-mode recovery
  (`solveMode` never applies it); the combined path matches that behavior.

## Frequency response & steady-state dynamics — `frequencyResponse` / `steadyStateDynamics` (P2)

For a **harmonic** base acceleration `ü_g(t) = amp·e^{iΩt}`, the steady response of
a linear structure is a small dense post-processor on the same mode basis — no time
stepping. Each mode has the complex transfer function

```
H_a(Ω) = 1 / (ω_a² − Ω² + i·Ω·d_a),   d_a = c_a/m_a,
```

and the physical (relative) displacement FRF is the modal sum
`û(Ω) = amp · Σ_a (ψ_a)(−Γ_a) H_a(Ω)` — the frequency-domain image of the P1a
recovery, so it inherits the same participation / eigenvector-scale normalization.

```python
ops.eigen('-fullGenLapack', nModes)
ops.modalProperties()
# full complex FRF of node `roof` dof 1 over 0.1..20 Hz, resonance-biased grid:
rows = ops.frequencyResponse('-freq', 0.1, 20.0, 400, '-biased',
                             '-baseAccel', '-dir', 1, '-rayleigh', a0, a1,
                             '-node', roof, '-dof', 1, '-out', 'frf.out')
# rows[k] = [f_Hz, Re, Im]

# steady-state harmonic amplitude (|response|) instead of the complex value:
amp = ops.steadyStateDynamics('-freq', 0.1, 20.0, 400, '-biased',
                              '-baseAccel', '-dir', 1, '-damp', 0.02,
                              '-node', roof, '-dof', 1, '-out', 'ssd.out')
# amp[k] = [f_Hz, |response|]
```

| flag | meaning |
|---|---|
| `-freq $fmin $fmax $nf` | sweep bounds **in Hz** and point count (`Ω = 2π f`) |
| `-lin` / `-log` / `-biased` | grid: uniform in `f` / geometric / linear + a ±5% cluster of points around each in-band modal frequency to resolve sharp peaks (default `-lin`) |
| `-baseAccel -dir $dir` | uniform harmonic base acceleration along global dir (no `timeSeries` — the sweep is per unit/`-amp` amplitude) |
| `-amp $a` | base-acceleration amplitude (default 1 → `frequencyResponse` is the pure transfer function) |
| `-damp`/`-rayleigh`/`-modalDamp` | same damping channels as P1a |
| `-node $tag -dof $dof` | the response DOF whose FRF is reported |
| `-resp disp\|vel\|accel` | response quantity (relative); `vel = iΩ·û`, `accel = −Ω²·û` (default `disp`) |
| `-modes $m1 …` | 1-based mode subset (default all) |
| `-out $file` | write the table to an ASCII file (whitespace columns); the table is also returned to the interpreter |

- `frequencyResponse` returns `[f, Re, Im]` per row; `steadyStateDynamics` returns
  `[f, |·|]`. Both also write `-out` if given.
- **Frequencies are Hz throughout** (in `-freq` and in the output `f` column). The
  FRF magnitude peaks at each damped modal frequency `f_a = ω_a/2π`; at `Ω→0` the
  response equals the static displacement under the load amplitude.
- **Sign convention is `e^{+iΩt}`** — the response lags 90° at resonance. Verified
  end-to-end against the direct complex solve `(K−Ω²M+iΩC)⁻¹(−MR)` in
  `modal_response_p2_spike/frf_oracle.py`.
- **Relative response.** For base excitation the reported disp/vel/accel are relative
  to the moving base (the P1a relative formulation); absolute acceleration would add
  the input back.
- Low-damping resonances are sharp — use `-biased` (or a fine `-log`) so the peak is
  not stepped over; an *undamped* mode sampled exactly at `Ω=ω_a` gives an infinite
  FRF by construction.

## Scope

P1a covers **uniform base-acceleration** transient excitation; P2 (`frequencyResponse`
/ `steadyStateDynamics`) covers **harmonic base-acceleration** frequency response —
both with classical (diagonal) modal damping. General load-pattern / nodal-force modal
loads, non-classical damping (see ADR 46 `complexEigen`), CQC/SRSS response-spectrum
combination (P1b), and random/PSD→RMS (P3) are separate phases of ADR 44.

Validation: `tests/test_ladrunoModalResponse.py` (SDOF exact vs numpy PWL across all
damping branches; multi-mode vs OpenSees Newmark and vs an independent full-matrix
Newmark; mode subset; damping-channel equivalence; guards). The closed-form
recurrence is cross-checked to ~1e-14 against a van-Loan `expm` first-order-hold
oracle in `modal_response_p1a_spike/pwl_recurrence_oracle.py`.
