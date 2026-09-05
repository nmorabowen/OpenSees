---
title: "ADR 90 leg A0 — results: does a generic two-track Duvaut–Lions wrapper regularize localization width?"
project: Ladruno
type: measurement note (phase P0, WP-A)
status: "COMPLETE 2026-09-04 — the numbers below decide D2; the ADR text is written against them"
owner: nmora
related:
  - "[[_adr90_regularization_planning_brief]]"     # §2 F1/F2/F4/F6, §5 leg A0, §7 D2/D3, §8 OQ3/OQ5, §9 R1/R4
  - "[[_adr90_orchestration_plan]]"                # WP-A row, checkpoint CP1
  - "[[31_ladruno_concrete3d_adr]]"                # the shipped material-specific -eta and its PV gate
tags: [adr, regularization, viscoplastic, duvaut-lions, localization, oracle, measurement]
updated: 2026-09-04
---

# ADR 90 leg A0 — results

> [!abstract] One line
> The generic two-track wrapper is not an approximation to true Duvaut–Lions — in 1‑D under
> monotonic loading it **is** true Duvaut–Lions, exactly, for any hardening law (proved and
> measured to 1e‑13). It stops being so on **non‑proportional** and **unloading** paths, where
> the two differ by 4.4 e‑2 … 3.3 e‑1. On the A0 bar the two are indistinguishable (≤ 1.3 e‑5),
> the τ = 0 negative control fails objectivity exactly as it must (width ≡ h, work halves with
> h), and τ > 0 at fixed De converges — to a width that is a **steep** function of De, not a
> material length.

## 0. Provenance

| | |
|---|---|
| Oracle | `tests/_testbed/duvaut_lions_ref.py` (numpy only; imports no OpenSees) |
| Driver | `tests/_testbed/run_a0_sweep.py` — regenerates every table here (~2 min) |
| Gate | `tests/test_duvaut_lions_oracle.py` — 9 `zone_a` cases, **22.6 s measured** |
| Engine | numpy 2.4.6, python3.12; **no C++ build exists for this work package** |
| Git | `cc7c7f7a5` on `wp/90a-adr-oracle` | 
| Date | 2026-09-04 |

Every number below was produced by that driver at that hash. Nothing here is quoted from memory
or from the literature.

## 1. The A0 deck, and why these parameters

1‑D bar, `L = 100 mm`, `A = 1 mm²`, `N ∈ {20, 40, 80, 160}` two-node linear elements,
quasi-static (no inertia), end-displacement ramp `u(t) = u_max·t/T` with `u_max = 2.0 mm` and a
**uniform** `Δt = T/nsteps`. Newton with the consistent (blended) tangent; a consistent-tangent
predictor at each step start; backtracking safeguard. `De = τ/T`.

Material: `E = 20000 MPa`, `σ_Y = 20 MPa`, linear softening `H = −50 MPa`, residual `0.02 σ_Y`
carried on a nearly-flat branch `H_res = E/2e5`.

**Why `H = −E/400`.** A softening bar under end-displacement control snaps back — and Newton
cannot follow it — when the band is shorter than `L|H|/E`; from
`du = ℓ_b dσ/E_t + (L−ℓ_b) dσ/E` with `E_t = EH/(E+H) < 0`, `du/dσ` changes sign at
`ℓ_b = L|H|/E`. Here that threshold is **0.25 mm**, a factor 2.5 below the finest element
(`h = 0.625 mm` at `N = 160`), so **no leg of the sweep snaps back** — not just `N = 20`. This
is stronger than the brief asked for and it is what lets the τ = 0 negative control actually run
at every mesh, which is the leg that has to work for H1 to mean anything.

**Why `H_res > 0`.** A perfectly-plastic element in a 1‑D chain makes the assembled structural
tangent *exactly singular* (the chain extends at constant force). `H_res = E/2e5` removes that
without touching the pre-residual branch; the return map stays exact (two linear branches, the
active one chosen at the updated α). Default in the oracle is `H_res = 0` — a non-zero default
silently turned a "perfectly plastic" bar into a hardening one and broke PV3 once.

**Three imperfection conventions** (10 % strength reduction in all three):

| tag | what | mesh behaviour |
|---|---|---|
| `one_element` | de Borst's convention (a): one central element at 0.9 σ_Y | the defect SHRINKS with N |
| `fixed_flat` | the brief's convention (b) read literally: a fixed 10 % of L (2/4/8/16 elements) all at 0.9 σ_Y | fixed size, but every element in the zone is EXACTLY as weak as its neighbours |
| `fixed_graded` | convention (b′): the same fixed physical zone with a parabolic notch, 10 % at the centre tapering to 0 at the zone edge, centre offset by 1.3 % of the half-width to break the even-N tie | fixed size AND mesh-convergent as a field |

## 2. The three band-width definitions

Measured on the **post-peak** plastic-strain increment profile `ε_p(end) − ε_p(peak)`:

* **w1** — threshold count: `h × #{e : p_e > 1e-3·max p}`.
* **w2** — threshold-**free** second moment: `√(12·Var)`,
  `Var = Σ p_e[(x_e−x̄)² + h²/12] / Σ p_e`. The `h²/12` term is the within-element variance of the
  piecewise-constant profile; it makes a one-element band read **exactly h** and a k-element top
  hat read **exactly k·h**, so w2 is directly comparable across meshes. Without it a one-element
  band reads 0. This is the brief's **OQ5** answer, and `test_band_width_metric_is_calibrated`
  pins it.
* **w3** — FWHM: `h × #{e : p_e ≥ 0.5 max p}`.

w1 and w3 differ only in their threshold, and the sweep shows them disagreeing by up to a factor
of 40 on the same run (risk **R4** confirmed) — w2 is the one to quote.

## 3. PV1–PV6, both variants, both point models

| gate | claim | measured | verdict |
|---|---|---|---|
| PV1 | τ = 0 is **byte-identical** to the inviscid material (1‑D and J2; H < 0, = 0, > 0; loading and unloading) | `max|Δσ| = 0.0` exactly | PASS |
| PV2 | discrete → continuous, first order, mid-transient | order 0.91 (TT), 0.91 (TDL) | PASS |
| PV3 | steady overstress `= E·ε̇·τ`, independent of Δt ∈ {1, 0.25, 0.05} | max rel err `4.26e-15` | PASS |
| PV4 | the update **is** `(1−β)σ_tr + β σ̄` | `0.0` | PASS |
| PV5 | blended tangent vs central FD (1‑D and 6×6 J2) | max rel `3.31e-10` | PASS |
| PV6 | overstress monotone in τ, and > 0 (probe asserted plastic) | monotone over τ ∈ {0.1, 0.5, 2, 10} | PASS |

## 4. TT vs TDL at the material point — **a theorem, not a tolerance**

### 4.1 The result

**Theorem (1‑D, from rest, monotonic).** For the 1‑D associated model with *any* hardening
function K(α), the two-track wrapper and true Duvaut–Lions produce the **identical** stress path.

*Proof.* The TDL update gives `σ_{n+1} = σ_tr − βE(ᾱ−α_n)` and `α_{n+1} = α_n + β(ᾱ−α_n)`, hence
`σ_{n+1} + Eα_{n+1} = σ_tr + Eα_n`; so `ψ_n = σ_n + Eα_n` advances by exactly `EΔε` every step,
elastic steps included. From rest `ψ_0 = 0`, therefore **α_n = ε_n − σ_n/E**: TDL's relaxed
internal variable is exactly the plastic strain of its own stress. Substituting into the
projection equation `σ_tr − E(ᾱ−α_n) = K(ᾱ)` gives `Eε_{n+1} − Eᾱ = K(ᾱ)`, i.e.
`σ̄ = K(ᾱ)` with `ᾱ = ε_{n+1} − σ̄/E` — the definition of the **inviscid** stress at the total
strain `ε_{n+1}`. So TDL's projection target *is* `inner.getStress()` on the total strain path,
which is precisely what the two-track wrapper blends toward. Both return
`(1−β)σ_tr + β σ_inviscid(ε_{n+1})`. ∎

### 4.2 Measured

| case | max relative |σ_TT − σ_TDL| | |
|---|---|---|
| perfect plasticity | `3.20e-14` | identity |
| 1‑D **linear**, H/E ∈ {+0.10, +0.02, 0, −0.02, −0.10} × De ∈ {0, 0.003 … 0.3} | `9.23e-14` | identity |
| 1‑D **exponential** (nonlinear) law, softening and saturation-hardening, α_end/α_f ≈ 4.75 | `2.80e-13` | identity |
| J2, **proportional** path, H/E = ±0.05 | `4.72e-15` | identity |
| J2, **non-proportional** (axial 120 steps, then shear 80) | **`4.42e-02`** | **differs** |
| 1‑D linear, load / **unload** / reload | **`3.33e-01`** | **differs** |

Selected rows (the full tables are in the driver output):

```
   H/E   path      De   |s_inv|    |s_TT|   |s_TDL|        rel
  0.05   prop    0.10   155.638   159.516   159.516   4.37e-15
  0.05   nonprop 0.01    93.712    93.890    93.892   4.02e-03
  0.05   nonprop 0.10    93.712    95.787    95.782   4.42e-02
 -0.05   nonprop 0.01    92.634    92.723    92.725   8.99e-04
 -0.05   nonprop 0.10    92.634    94.024    94.051   2.62e-02

  cyclic (load/unload/reload, 1-D linear)
   H/E      De  sig_TT end  sig_TDL end  max rel diff
  0.10   0.010     33.4544      33.4545      4.63e-02
  0.10   0.100     45.1909      49.0359      3.14e-01
 -0.02   0.010     19.4152      19.5593      9.23e-02
 -0.02   0.100     35.7801      37.4076      3.33e-01
```

The difference on the non-proportional path grows monotonically with De (8.99e-4 → 2.62e-2 as
De goes 0.01 → 0.10), which is the signature of a genuine constitutive difference rather than
round-off.

> **A trap worth recording.** The first exponential-law run used `α_f = 0.5` while the path only
> reached `α ≈ 5e-3`; over `α/α_f ≈ 0.01` the exponential *is* linear, so it reproduced the
> linear-law identity to 1e-14 and would have been reported as "nonlinear laws agree too". The
> legs in the shipped oracle drive `α_end/α_f ≈ 4.75` (≈ 99 % of the law's range). Same family as
> the `-eta` gate's own tautological-fixture trap (`LEDGER_quirks.md`, `-eta` entry).

## 5. The A0 bar

### 5.1 H1 — negative control, τ = 0

```
convention        N        h  nsteps        w1        w2        w3    w2/h      W_50  W50 ratio   P_peak
one_element      20   5.0000    2000    5.0000    5.0000    5.0000    1.00   12.1500        n/a  18.0000
one_element      40   2.5000    2000    2.5000    2.5000    2.5000    1.00    6.0750      2.000  18.0000
one_element      80   1.2500    2000    1.2500    1.2500    1.2500    1.00    3.0375      2.000  18.0000
one_element     160   0.6250    2000    0.6250    0.6250    0.6250    1.00    1.5187      2.000  18.0000
fixed_graded     20   5.0000    4000    5.0000    5.0000    5.0000    1.00   12.7994        n/a  18.4730
fixed_graded     40   2.5000   16000    2.5000    2.5000    2.5000    1.00    6.1514      2.081  18.1109
fixed_graded     80   1.2500   32000    1.2500    1.2500    1.2500    1.00    3.0460      2.020  18.0250
fixed_graded    160   0.6250   64000    0.6250    0.6250    0.6250    1.00    1.5196      2.004  18.0040
fixed_flat       20   5.0000    2000   10.0000   10.0000   10.0000    2.00   24.3000        n/a  18.0000
fixed_flat       40   2.5000    2000    5.0000   13.2288    5.0000    5.29   12.1500      2.000  18.0000
fixed_flat       80   1.2500    2000    6.2500   10.1026    6.2500    8.08   15.1875      0.800  18.0000
fixed_flat      160   0.6250    2000    2.5000    8.7332    2.5000   13.97    6.0750      2.500  18.0000
```

`W_50` is the dissipated work at the instant the load first falls to **50 % of its peak** — a
common point on the softening curve, and therefore a mesh-*fair* energy measure. Comparing `W` at
a common end displacement is not fair: at fixed u a coarse band has softened far less than a fine
one.

* `one_element` and `fixed_graded`: **w1 = w2 = w3 = h exactly** at every mesh, and `W_50` halves
  with h to 3 significant figures. Textbook loss of objectivity.
* `fixed_flat` is **not usable as a negative control**: the tie among the weak elements means
  round-off picks the sub-band, and w2/h wanders 2.0 → 14.0 with no pattern.

### 5.2 H2 — at fixed De, does it converge? (`fixed_graded`, TT; TDL identical to 6 digits)

```
      De    N        w2  w2 ratio       w3   P_peak    P_end     W_end  curveL2 vs N/2  eps_p max
 1.0e-04   20    5.9257       n/a   5.0000  18.8619   1.2078   20.0105             n/a     0.3843
 1.0e-04   40    3.2703     0.552   2.5000  18.8177   2.0363   13.8671           0.217     0.7439
 1.0e-04   80    1.9700     0.602   1.2500  18.8314   3.7004   13.7058           0.112     1.4120
 1.0e-04  160    1.6910     0.858   1.2500  18.8184   3.7052   13.2451           0.040     1.5145
 3.0e-04   20    8.2485       n/a   5.0000  19.6026   4.7171   26.5897             n/a     0.3255
 3.0e-04   40    4.8026     0.582   2.5000  19.7240   5.2184   22.4539           0.132     0.5711
 3.0e-04   80    3.6286     0.756   2.5000  19.7519   5.2276   21.1199           0.071     0.6917
 3.0e-04  160    3.5266     0.972   2.5000  19.7538   5.2270   21.4887           0.021     0.7512
 5.0e-04   20   12.2547       n/a  10.0000  20.0426  10.0789   30.5649             n/a     0.2519
 5.0e-04   40   11.7026     0.955   5.0000  20.0558   4.8216   26.4621           0.155     0.3938
 5.0e-04   80   11.9144     1.018   2.5000  20.0590   8.4133   27.4292           0.056     0.5469
 5.0e-04  160   12.0320     1.010   3.1250  20.0597   7.0253   27.0102           0.024     0.5372
 1.0e-03   20   38.8560       n/a  10.0000  20.2418  14.2054   35.2609             n/a     0.1729
 1.0e-03   40   41.1880     1.060   5.0000  20.2543  13.5267   35.5675           0.012     0.2625
 1.0e-03   80   41.8736     1.017   5.0000  20.2575  13.4820   35.3967           0.011     0.2754
 1.0e-03  160   42.0407     1.004   6.2500  20.2583  13.4343   35.4125           0.002     0.2798
 3.0e-03   20   89.9698       n/a  10.0000  21.0303  20.1003   39.0775             n/a     0.0524
 3.0e-03   40   90.8211     1.009  10.0000  21.0428  20.1128   39.1000           0.001     0.0616
 3.0e-03   80   91.0308     1.002   7.5000  21.0460  20.1159   39.1056           0.000     0.0638
 3.0e-03  160   91.0831     1.001   8.1250  21.0467  20.1167   39.1070           0.000     0.0643
```

`w2 ratio` is w2(N)/w2(N/2) — **0.500 exactly** would be the τ = 0 behaviour. `curveL2` is the
normalized L2 distance between the load–displacement curve at N and at N/2.

* De ≥ 3e‑4 **converges**: ratios 0.972 / 1.010 / 1.004 / 1.001 at the finest pair, and the
  load–displacement curves converge (`curveL2` falls by ~3× per refinement).
* De = 1e‑4 does **not** converge in width (ratios still tracking 0.5–0.86) although its *curve*
  does (0.217 → 0.040). Width and curve convergence are **not** the same gate.
* The converged width is a **steep, monotone** function of De:
  `De = 3e-4 → 3.53 mm`, `5e-4 → 12.0 mm`, `1e-3 → 42.0 mm`, `3e-3 → 91.1 mm` (the bar is
  100 mm). A factor 3 in De moves the width by a factor 12. It is **not** the imperfection width
  (10 mm) — it crosses it.
* With `one_element` (convention a) the width does **not** converge at any De: the peak stays one
  element wide (`w3 = h` at every mesh) and w2 drifts toward 100 mm. The *curve* does converge for
  De ≥ 1e‑3 (`curveL2` 0.116 → 0.027 → 0.001). This is the conventions' whole point: an
  imperfection field that does not converge cannot give a converging width.
* `fixed_flat` is **bit-identical across N at every De > 0** — a tautological pass (see §5.6).

### 5.3 H3 — De collapse

(a) With the **same step count** the collapse is *exact by construction* and therefore worthless
as evidence: `β = Δt/(τ+Δt) = 1/(1 + nsteps·De)` and the strain increment `u_max/nsteps` carry no
other dependence on τ or T, so equal `(De, nsteps)` runs are bit-identical whatever `(τ, T)`
produced that De.

```
       tau         T         De                beta                w2            P_peak             W_end
   1.0e-03         1    1.0e-03    3.3333333333e-01   41.873586502200   20.257474397731   35.396701347229
   1.0e-02        10    1.0e-03    3.3333333333e-01   41.873586502200   20.257474397731   35.396701347229
   1.0e-04       0.1    1.0e-03    3.3333333333e-01   41.873586502200   20.257474397731   35.396701347229
```

(b) The substantive version — matched De with **different** step counts, so β differs:

```
       tau         T  nsteps        beta        w2    P_peak     W_end
   1.0e-03         1    1000  5.0000e-01   41.7324   20.2562   35.3909
   1.0e-02        10    2000  3.3333e-01   41.8736   20.2575   35.3967
   1.0e-04       0.1    4000  2.0000e-01   41.9292   20.2581   35.3996
```

Spread: 0.47 % in w2, 0.009 % in peak load, 0.025 % in dissipated work — and the entire residue
is the Δt-transient effect of §5.5, not a De violation.

### 5.4 H4 — TT vs TDL on the bar

Maximum relative differences over `De ∈ {0, 3e-4, 1e-3, 3e-3}` × `N ∈ {40, 80, 160}` × both
imperfection conventions:

| quantity | max relative difference |
|---|---|
| band width w2 | **1.31e-05** |
| peak load | 1.8e-16 |
| dissipated work | 1.25e-06 |
| load–displacement curve (L2) | 4.82e-06 |

Exactly 0 at τ = 0 and at De = 3e‑3. The ~1e‑5 residue comes from elements that yield and then
unload (where the §4 theorem's hypotheses fail) plus the Newton tolerance.

### 5.5 Δt sensitivity at fixed De (N = 80, De = 1e‑3, `fixed_graded`)

```
  nsteps          dt        beta        w2       w3    P_peak     P_end     W_end
     250   4.000e-03  8.0000e-01   41.1095   5.0000   20.2503   13.4356   35.3539
     500   2.000e-03  6.6667e-01   41.5724   5.0000   20.2538   13.4628   35.3791
    1000   1.000e-03  5.0000e-01   41.7324   5.0000   20.2562   13.4755   35.3909
    2000   5.000e-04  3.3333e-01   41.8736   5.0000   20.2575   13.4820   35.3967
    8000   1.250e-04  1.1111e-01   41.9568   5.0000   20.2585   13.4869   35.4011
```

The transient **is** Δt-dependent even though PV3's steady overstress is not: w2 moves 2.06 %
over a 32× refinement, converging at first order (successive differences 0.463, 0.160, 0.141,
0.056, 0.028). Peak load moves 0.004 %. Practical consequence for the C++ wrapper and for WP-D:
**a reported width must carry its step count**, not only τ and T.

### 5.6 Two gates that pass tautologically — recorded, not hidden

1. **De collapse at fixed nsteps** (§5.3a). Cannot fail. `test_a0_two_gates_that_pass_tautologically`
   asserts the bit-identity so the ADR cannot later cite it as evidence.
2. **A flat fixed-length imperfection.** With every element of the weak zone exactly as weak as
   its neighbours, the continuum solution is piecewise constant with the zone boundary on a mesh
   line, so **every mesh represents it exactly** and the answer is bit-identical in N at every
   De > 0. Convention (b) must be run **graded**, not flat.

### 5.7 The brief's requested De ∈ {0.01, 0.1, 1} — a non-tautology check

```
     De    N    P_peak  eps_p max  inelastic        w2       w3
   0.01   20   23.7884    0.02627     31.128   97.8650 100.0000
   0.01  160   23.8048    0.02890     31.105   98.0932 100.0000
   0.10   20   59.2479    0.01771     11.086   99.8105 100.0000
   0.10  160   59.2643    0.01795     11.083   99.8305 100.0000
   1.00   20  264.8369    0.00682      0.757    0.0000   0.0000
   1.00  160  264.8470    0.00685      0.757    0.0000   0.0000
```

All three are far past the point where this bar still localizes: at De = 0.01 and 0.1 the "band"
is the **whole 100 mm bar** (w3 = 100). At **De = 1 the load never peaks at all** — it rises
monotonically to 265 N, 14.7× the rate-independent peak of 18 N — so there is no post-peak branch
to measure and "the width converges" is a statement about an essentially elastic bar. The De
range where this deck is both regularized *and* still localizing is **3e‑4 … 1e‑3**, two orders
of magnitude below the brief's guess.

### 5.8 Steps-to-converge — the cost of *not* regularizing (`fixed_graded`, TT)

```
    N     tau=0   De=3e-4   De=1e-3   De=3e-3
   20      4000       250       250       250
   40     16000       250       250       250
   80     32000       250       250       250
  160     64000       250       250       250
```

The τ = 0 problem needs a step count that **grows ×4 then ×2 per refinement, 16× in all across
the sweep** — the numerical signature of ill-posedness — while every regularized run completes at
the coarsest step count tried (250), independent of mesh. This is arguably the cleanest single measurement in the note, and it
is free.

## 6. Hypothesis verdicts

| # | Hypothesis | Verdict | Carrying number |
|---|---|---|---|
| **H1** | τ = 0 ⇒ width ∝ h, dissipated work → 0 with h | **CONFIRMED** (conventions a and b′); **not testable** on convention b-flat | w1 = w2 = w3 = h *exactly* at all four meshes; `W_50` ratios 2.000 / 2.000 / 2.000 |
| **H2** | at fixed De, width and curve converge with N, for both TT and TDL | **CONFIRMED for De ≥ 3e‑4 with a mesh-convergent imperfection field; REFUTED for a one-element imperfection and for De ≤ 1e‑4** | w2(160)/w2(80) = 0.972 (De 3e‑4), 1.010 (5e‑4), 1.004 (1e‑3), 1.001 (3e‑3); vs 0.858 at De 1e‑4 and 1.356 for `one_element` at De 3e‑4 |
| **H3** | matched (τ, T) at equal De give the same result | **CONFIRMED, but the naive form is TAUTOLOGICAL**: identical at fixed nsteps by construction (bit-identical); 0.47 % apart when Δt differs, and that residue is the Δt effect | see §5.3 |
| **H4** | TT vs TDL differ measurably; is TT adequate for a softening law? | **TT is adequate — and the reason is a theorem, not a tolerance.** They are *identical* in 1‑D under monotonic loading for any hardening law, and indistinguishable on the whole A0 bar. They differ only on **non-proportional** and **unloading** paths | bar: ≤ 1.31e‑5; point, non-proportional J2: 4.42e‑2; point, cyclic: 3.33e‑1 |

## 7. Recommendation for decision **D2**

> **Build the generic two-track wrapper (brief §7 D2, first option). Do NOT pre-emptively build a
> material-specific `-tau` inside `LadrunoSANISAND` (WP-F). Scope the ADR's correctness claim to
> what §4 proves, and put the non-proportional/unloading gap on the record as the wrapper's
> declared approximation.**

The evidence line that carries it, in one sentence: *the two-track blend is not an approximation
to Duvaut–Lions but an exact reformulation of it whenever the internal variable is recoverable
from the total strain and the current stress (proved in §4.1, measured at 1e‑13 for linear,
nonlinear-exponential and proportional-J2 laws, and at ≤ 1.3e‑5 on the whole A0 structure) — so
the wrapper's error is confined to non-proportional and unloading paths, where it is 4.4e‑2 …
3.3e‑1 and must be disclosed, not to hardening or softening per se.*

Two riders the ADR must carry with that recommendation:

1. **The claim must be stated over paths, not over material classes.** The brief's F2 framed the
   two-track/true-DL split as "perfect plasticity vs hardening". That framing is **wrong** and
   §4 refutes it: hardening, softening and even strongly nonlinear hardening are all in the exact
   regime. The real boundary is *proportional-and-monotonic* vs *non-proportional-or-unloading*.
   SANISAND sits on the wrong side of that boundary — its α-tensor, fabric `z` and ψ-driven `M^b`
   make even a "monotonic" radial pushover non-proportional in the relevant sense — so the
   wrapper over SANISAND is an approximation whose size is **not yet measured** (this oracle
   measured a J2 surrogate: 0.09 % at De = 0.01, 4.4 % at De = 0.10). Measuring it on the real
   material is a WP-C/WP-D item, not a reason to block WP-C.
2. **WP-F stays parked, with a trigger.** Fire it only if the A2/A3 legs measure a
   wrapper-vs-inner discrepancy on SANISAND that exceeds the band TIMs declares. The A0 bar gives
   no reason to fire it now.

## 8. What this changes in the planning brief

| brief item | status after A0 |
|---|---|
| **OQ3** (is two-track adequate?) | **ANSWERED** — yes, with the path-scoped caveat of §7. |
| **OQ5** (threshold-free width metric) | **ANSWERED** — the calibrated second moment of §2, pinned by a unit test. |
| **R1** (BLOCKING: two-track ≠ DL for hardening, claim silently over-reaches) | **RETIRED as stated, REPLACED**: the risk is not hardening, it is non-proportional/unloading paths. Re-word in the ADR. |
| **R4** (width metric threshold-dependent, reverses across meshes) | **CONFIRMED as a real risk** — w1 and w3 disagree by up to 40× on the same run; mitigated by quoting w2. |
| **F2** ("for hardening / state-dependent models they do not [coincide]") | **Partially wrong** — hardening alone does not separate them. Correct the sentence. |
| **F4** (de Borst: quasi-static rate regularization is weak, no internal length) | **CORROBORATED and quantified**: there is no intrinsic length here; the converged width is set by De and moves a factor 12 for a factor 3 in De, and the *imperfection field's* mesh-convergence is a precondition for any width convergence at all. |
| brief §5 leg A0 gate wording ("τ = 0 → width ∝ h") | keep, and add: the negative control must use a **unique-weakest-point** imperfection; a flat weak zone cannot fail. |
| brief §5 "De × {½, 1, 2} matched pairs collapse" | keep, but **require different step counts** in the matched pairs, else the gate is bit-identical by construction. |

## 9. What this oracle does *not* settle

* It is 1‑D. The §4 theorem is a 1‑D statement; the multiaxial evidence is a J2 surrogate, not
  SANISAND. The non-associated Drucker–Prager plane-strain point model was scoped as optional and
  **was not built** — it is the natural first item if WP-C wants a second multiaxial datapoint
  before the A1 deck.
* It has no wave speed, so it cannot say anything about the dynamic internal length `ℓ = 2mc_e/E`
  that lane (a) of the brief's §3 would introduce.
* The A0 widths are specific to this deck's softening ductility (`σ_Y/|H| = 0.4` strain to
  residual). The *shape* of the De → width curve should transfer; the numbers should not be
  quoted as targets for A1/A2/A3.
* Nothing here has been run through the C++ engine. Every number is numpy.
