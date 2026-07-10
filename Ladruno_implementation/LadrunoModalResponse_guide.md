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
| `-baseAccel -dir $dir` | uniform harmonic base acceleration along global dir (no `timeSeries` — the sweep is per unit/`-amp` amplitude); **relative** response |
| `-load $patternTag` | harmonic **nodal forces** `amp·P·e^{iΩt}` — `P` is the pattern's plain `NodalLoad` reference values (the pattern's own `timeSeries` is **ignored**; patterns carrying eleLoads/sp constraints, thermal nodal loads, or loads on nodes without eigenvectors are refused); **absolute** response. Mutually exclusive with `-baseAccel` |
| `-amp $a` | excitation amplitude (default 1 → `frequencyResponse` is the pure transfer function) |
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
- **Relative vs absolute.** For base excitation the reported disp/vel/accel are
  relative to the moving base (the P1a relative formulation); absolute acceleration
  would add the input back. For `-load` the response is absolute.
- **`-load` normalization.** The modal load is `f_a = ψ_aᵀP/m̃_a` with
  `ψ_a = evec·Vscale` and `m̃_a = generalizedMasses()(a)` — the exact basis
  `modalProperties` computes its products in, so the result is identical under
  default and `-unorm` normalization (pinned by test + the spike
  `modal_response_p3_spike/load_frf_oracle.py`; `-baseAccel` is the special case
  `P = −M·R`). Static limit: `u(Ω→0) = K⁻¹P` with all modes retained.
- Low-damping resonances are sharp — use `-biased` (or a fine `-log`) so the peak is
  not stepped over; an *undamped* mode sampled exactly at `Ω=ω_a` gives an infinite
  FRF by construction.

## Random response — `randomResponse` (P3)

For a **stationary random** base acceleration given as a **one-sided PSD `G(f)` in
Hz** ((accel)²/Hz — the engineering-spec convention: `σ_üg² = ∫₀^∞ G df`; vs the
textbook *two-sided rad/s* PSD, `G(f) = 4π·S(Ω)`), the response auto-PSD rides the
P2 FRF unchanged:

```
G_xx(f) = |H_x(f)|² · G(f),
m_k = ∫ f^k G_xx df  →  RMS = √m0,   ν₀ = √(m2/m0) [Hz],
E[peak] = (√(2·ln(ν₀T)) + 0.5772/√(2·ln(ν₀T))) · RMS   (with -duration T)
```

```python
ops.eigen('-fullGenLapack', nModes)
ops.modalProperties()
# input PSD as a timeSeries SAMPLED AT f IN Hz — e.g. band-limited flat:
ops.timeSeries('Path', 1, '-time', 0.0, 0.999*f1, f1, f2,
               '-values', 0.0, 0.0, G0, G0)
rms = ops.randomResponse('-freq', f1, f2, 800, '-biased',
                         '-baseAccel', '-dir', 1, '-inputPSD', 1,
                         '-damp', 0.03, '-node', roof, '-dof', 1)

# with statistics (and the Davenport expected peak over a 600 s exposure):
rms, nu0, m0, m2, peak = ops.randomResponse(..., '-stats', '-duration', 600.0)
```

| flag | meaning |
|---|---|
| `-inputPSD $tsTag` | **required** — one-sided PSD `G(f)` [Hz], a `timeSeries` sampled at `f` (Path with `f→G` breakpoints, Constant for white noise). With `-baseAccel`: PSD of the base acceleration. With `-load $patternTag`: PSD of the scalar `s(t)` multiplying the pattern's nodal-load shape (`P(t)=s(t)·P`, fully correlated; SDOF anchor `σ_u²=G_F/(8ξω³m²)`) |
| `-stats` | return `[rms, ν₀, m0, m2]` instead of the scalar RMS |
| `-duration $T` | append the Davenport expected peak as a 5th entry (implies the `-stats` list form even without `-stats`; the entry is **NaN** when `ν₀·T ≤ 1`, and the estimate is flagged unreliable for `ν₀·T < 2`) |
| `-out $file` | write `{f, G_in, G_xx}` rows |
| *(rest)* | `-freq/-lin/-log/-biased`, `-baseAccel -dir`, damping channels, `-node/-dof`, `-resp`, `-modes` — same as P2 (no `-amp`: the PSD carries the excitation scale) |

- **Convention pin.** White-noise SDOF anchor in this convention:
  `σ_x² = G0/(8ξω³)` (and `σ_v = ω·σ_x`). Pinned analytically AND by a Monte-Carlo
  synthetic realization in `modal_response_p3_spike/psd_rms_oracle.py` — a
  one-sided/two-sided mixup reads ~41% (√2), an Hz/rad mixup ~150% (√2π).
- **Use `-biased`.** The RMS is a band integral; a 200-pt `-lin` grid mis-integrates
  a ξ=0.5% resonance by ~14% where the biased grid reads 0.05%.
- **Refusals.** `nf ≥ 2` is required (a single point has zero measure); `G(f) < 0`
  anywhere on the grid is refused; a **zero-damped mode inside the band is refused**
  (its variance integral diverges — either damp it or exclude its resonance from
  the band; undamped modes strictly *outside* the band are legal). A **rigid-body
  mode (ω=0) with `fmin = 0`** is refused for *any* damping — its FRF diverges at
  its own `f=0` even with Rayleigh `a0 > 0` (use `fmin > 0` or `-modes` to drop it).
- The band is `[fmin, fmax]` — energy outside it is *not* counted. Make the band
  cover both the input's support and every resonant peak that carries response
  power (the white-noise tests use `fmax ≈ 30·fn` for sub-0.2% truncation).
- `ν₀` is the mean zero-upcrossing rate (narrow-band response: `ν₀ ≈ f_n`); the
  Davenport factor assumes a stationary Gaussian process over the exposure `T`.

## Scope

P1a covers **uniform base-acceleration** transient excitation; P2 (`frequencyResponse`
/ `steadyStateDynamics`) covers **harmonic** frequency response and P3
(`randomResponse`) **stationary random** response (single auto-PSD input →
RMS/ν₀/peak) — each for **base acceleration or a nodal-force `-load` pattern**
(the `-load` follow-up), all with classical (diagonal) modal damping.
Remaining follow-ups: `-load` excitation for the P1a *transient*, cross-PSD
(matrix `S_PP`) input; non-classical damping → ADR 46 `complexEigen`.

Validation: `tests/test_ladrunoModalResponse.py` (SDOF exact vs numpy PWL across all
damping branches; multi-mode vs OpenSees Newmark and vs an independent full-matrix
Newmark; mode subset; damping-channel equivalence; guards). The closed-form
recurrence is cross-checked to ~1e-14 against a van-Loan `expm` first-order-hold
oracle in `modal_response_p1a_spike/pwl_recurrence_oracle.py`.
