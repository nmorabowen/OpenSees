# WP-128 — SANISAND ring-state trace: can any integrator take the TIMs ring point? (TIMs F18(e) + the F18(a) baseline)

Investigation only: no C++ change, no build. Everything runs on WP-127's binary
(`ladrunoSANISANDReplay` + `substepStats`, build `234a7575`), in the WP-128 worktree.
Every number below is printed by a committed script in
`Ladruno_files/testbed/sanisand_ring_trace/` and saved under `out/` (§8 lists the
commands). Source lines are the WP-127 tree (`e8fb51cdb`) unless stated.

> **Correction (2026-09-27, after WP-134, PR #872).** The independent reference
> integrator ([[134_sanisand_reference_integrator]], exact Radau integration of the
> DM04 equations, not a port of the C++) overturns this report's ranking of
> mechanism **F** below. F (a negative Λ denominator taken as ELASTIC, plus the
> uncapped step factor) is not merely a compounder: it drives the substep error to
> **exactly 0**, so no tolerance can catch it; it accounts for **all 25**
> campaign-ModifiedEuler ring escapes from admissible starts; and it gives 20–65 %
> stress errors on benign 20–100 kPa states even at TolE 1e-8. This report's
> method could not see that: `md_port.py` shares the C++ error measure, so it
> judged F only by switching it off inside the same flawed estimator. WP-134 also
> found two defects this report does not list: **U9** (ModifiedEuler never
> re-evaluates K, G inside the increment: 0.6/6/24 % of the stress increment at
> δ = 1e-5/1e-4/1e-3, invisible to the error test) and **U10** (the loading test
> uses n:Δσ, not ∂f/∂σ:Δσ). G (trigger) and E (enabler) stand as reported. The
> WP-129 SAS-ME spec addresses all of them.

## 0. Answers in one page

1. **Finding B — how α gets ~6× outside the bounding surface.** Reproduced from a
   benign start in **one** accepted ModifiedEuler substep. The launch state is the
   one `Stress_Correction`'s low-p branch leaves (`σ = p_min·I`, `α = 0`), followed by
   a loading reversal that sets `α_in := α = 0`. On the next compression increment
   `(α − α_in):n = 0`, so `h` is `GetStateDependent`'s `1e10` sentinel. The first Heun
   stage then moves α by about `Δs/p_start`, evaluated at the floor pressure while
   the increment raises p roughly 20×. The second stage takes the `dγ < 0` branch,
   which is elastic in stress. Both stages have the same stress increment, so
   ModifiedEuler's **stress-only** error test (finding E) passes at `dT = 1`.
   **Smallest reproducer:** `σ = 0.0101·I, α = α_in = z = 0`, one plane-strain
   `dε_yy = +1e-4`. It returns `η = 10.79`, `α/α^b = 5.14` and `rc = 0` in one
   substep. The model's own answer, integrated with α in the error, is
   `α/α^b = 0.27` (§2).
   - The trigger is the pressure ratio of one increment, not its size. The escape
     sets in between `Δp/p = 3.3` (α/α^b 0.76) and `5.2` (2.16). At the floor that
     is a volumetric increment between 1e-5 and 3e-5.
   - The ring carries the signature. All four b8 rows with α outside the bounding
     surface have `α_in ≡ 0` exactly, and no row with `α_in ≠ 0` is outside.
   - The dT_min forced accept (finding C) is **not** needed: refusing it changes
     nothing. `Stress_Correction` supplies the launch state, but it is not the jump.
   - An error test that also measures α removes the escape on every path tried:
     the port with α in the error, and RK45 in C++. Measured in §2.4 and §6.
   - With the orchestrator's candidates F and G (§2.5), the ranking is:
     - **G** is the trigger: α_in is re-seated once per increment, so
       `(α − α_in):n` goes to 0 (h = 1e10) and then negative (h < 0) inside the
       substep. 37 of 38 crossing substeps have it, and a Macaulay bracket on it alone
       keeps α inside.
     - **E** is the enabler: the stress-only error.
     - **F**, the negative-denominator elastic drag, compounds the escape but is
       not necessary. The uncapped q, the dT_min forced accept and the
       `Stress_Correction` exit do not produce it.
2. **The f > 0 anomaly is a defect.** On an unloading increment through the tiny
   cone (`m = 0.005`), the same pair of stages is accepted with α moved. Then
   `Stress_Correction` cannot reduce `|f|` and takes its silent return: the
   "Couldn't decrease the yield function" branch, `:3054-3058`, which prints only
   under `debugFlag`. It hands back the uncorrected state with `f > 0` and `rc = 0`.
   This breaks the material checklist's "a non-converged return map must FAIL"
   (§3). Worst measured: `f = 11.2 kPa` committed at `p = 0.58` (q1_threshold).
3. **F18(e), plainly: no integrator takes the b8 worst point.**
   - The dumped state is **inadmissible**: α is 6.26× past the model's own `α^b_θ`,
     and `b:n = −8.19`.
   - Loading-type probes are "taken" trivially: the stress rides the cone around
     the bad α.
   - Every unloading-type probe either commits `f > 0` as success (ME at both
     tolerances) or **teleports** η from 12.9 to 1.33 through the dT_min
     forced-accept `Mc` clamp. That happens under ME at 1e-8, CPPM's ME fallback,
     the α-aware port, and RK45 even on a `1e-7` increment.
   - A 12.87 → 1.33 jump in one increment is not an integration.
   - An integrator should **refuse** such a state with a named code, not project it
     (§4).
4. **F18(a) baseline.** Today's norm is exactly `‖dσ₂−dσ₁‖ / max(2‖σ‖, 1 kPa)`:
   the "0.5 kPa switch" is a continuous 1 kPa floor, bit-identical on 640/640
   increments. So `σ_ref ∈ {0, 0.1, 1}` is inert wherever `‖σ‖ > 0.5`, which is
   every constant-p state and all but a handful of ring states. `σ_ref = 5` acts
   only where `‖σ‖ < 2.5 kPa`.
   - The substep count is **stability-limited**, not accuracy-limited. At
     constant p = 2 kPa and η/M^b = 1.00 it is 60 / 54 / 52 substeps at TolE
     1e-6 / 1e-5 / 1e-4. An accuracy-limited Heun would give 10×.
   - A 20 kPa floor saves 2 substeps there (52 → 50) and doubles the error.
   - Measured errors per 1e-5 increment:
     - constant-p chains: 1–2e-4 kPa at p = 2 and 5, 9e-4–2e-3 at p = 20
       (2e-5 relative);
     - ring states: median 1e-4, p95 3e-2, max 0.84 kPa, against an α-blind
       reference that under-states them (§5.4).
   - In the act's NormUnbalance, the smooth-state error sits below a
     1e-5·18.4 kN/m tolerance, about 2e-3 kPa at one Gauss point. The ring tail
     does not (§5.5).
5. **WP-129.** The error floor is the wrong fix: inert or a palliative, it buys
   nothing where the cost is. The minimal α-stability fix, each item with its
   proving test in §6:
   - (a) put α and z in ModifiedEuler's error, Sloan–Abbo–Sheng style;
   - (b) replace the dT_min forced accept with a refusal;
   - (c) make `Stress_Correction`'s give-up a refusal;
   - (d) check the committed state's admissibility on entry.

   (a) alone removes the escapes and all `f > 0` commits but leaves the Mc-clamp
   teleports, so it must ship with (b).

Side note (§7): the two flip-determinism failures are real, not an MKL or machine
effect. The first push step fails identically under `system FullGeneral`. At a
looser `NormDispIncr 1e-6` it converges to 9.7313 kN/m, against the test's pinned
9.659111 on `48c0e99bc`, so the deck's numerics moved after that build. The drift
was not bisected.

## 1. Harness

- **C++ replay** (`_boot.py`): four prototypes of the campaign set (`nu 0.312885`,
  `-Pmin 0.0101 -Presidual 0 -flipAlphaIn init`):
  - `ME(campaign)`: IntScheme 1, `-honorTolR 0` (TolE 1e-4), cap 20000;
  - `ME(1e-8,honor)`;
  - `RK45`: IntScheme **45**. The intake's "IntScheme 4" is `INT_MAXENE_FE`.
    `#define INT_RungeKutta45 45`.
  - `CPPM`: IntScheme 2, TolR 1e-10.

  Attachments are replayed `-convention compressionPositive` (finding A).
- **Chain driver** (`drive.py`): a single-point strain history is a chain of
  committed replays, each started from the previous return with `prevIncrNorm` =
  the previous increment. WP-127 pinned replay(k+1 from committed k) = analysis to
  1e-9. `run_const_p` adds a mixed-control chain: `dε_yy` is prescribed, and
  `dε_xx` is solved by secant so p stays at p₀ (miss ≤ 2e-9 kPa).
- **Python port of the IntScheme-1 path** (`md_port.py`). It is a line-by-line
  port of `integrate` (reversal test) → P2-5 guard → `explicit_integrator` →
  `IntersectionFactor(_Unloading)` → `ModifiedEuler` → `Stress_Correction`. It
  exists so that ModifiedEuler's substeps can be opened up and the error norm
  varied offline.
  - **Validated against the C++** on 80 ring rows × 20 probes (1600 increments,
    `validate_port.py`): 1581 identical substep census and state to **5e-11**.
  - The other 19 are round-off-chaotic in the C++ itself: perturbing σ by 1e-14
    changes the C++ substep count. One of them, 19717 vs a 20000 cap, is a
    cap-threshold near-miss.
  - Counterfactual switches, OFF by default, attribution only:
    - `alpha_err` / `fabric_err`: α / z in the substep error, RK45's form
      (absolute below norm 0.5, `/2‖·‖` above);
    - `drag="frozen"`: the `dγ<0` stage leaves α alone;
    - `forced_policy="refuse"`;
    - `correction=False`;
    - `err_floor=σ_ref`: the F18(a) norm.

## 2. Q1 — finding B: how a committed α gets outside the bounding surface

### 2.1 Search (`q1_search.py`)

K0 starts at p₀ = 2 / 5 / 20 kPa. Five plane-strain directions (passive, active,
shear, extension+shear, vertical unload), each run as 3 cycles of (n forward, n/2
back) at δ = 1e-6 / 1e-5 / 1e-4. α/α^b uses the model's own `α^b_θ` at the point's
Lode angle.

| path | δ = 1e-6 | δ = 1e-5 | δ = 1e-4 |
|---|---|---|---|
| passive / active / shear (all p₀) | ≤ 0.99 | ≤ 0.97 | ≤ 1.00 |
| `extShear` (all p₀) | 0.63 | 0.97 | **7.07** (first at k = 20) |
| `vertUnload` (all p₀) | 0.59 | 1.19 (k = 237, p 80 kPa) | **5.22** (first at k = 20) |

k = 20 is the **first reversal**. The first 20 increments push p to the floor,
where `Stress_Correction`'s low-p branch (`:2939-3010`) sets `σ = (p_min + p_r)·I`
and `α = 0` on every increment (q1_attrib's table: p_after = 0.0101 on k = 0..19,
the ~1930-substep lines). None of the three smaller-δ legs escapes.

### 2.2 Attribution (`q1_attrib.py`, port = C++ to 3e-3 on the chain)

The escaping increment (vertUnload, k = 20), opened:

- the reversal test sets `α_in := α_n = 0`, so `(α − α_in):n = 0` and **`h = 1e10`**
  (`GetStateDependent`, `:5398-5399`);
- **one substep**, `T = 0 → 1`, accepted, stage kinds `('plastic', 'elasticDrag')`:
  - stage 1 has `h = 1e10` and `Kp = 1.2e8`, so `dγ ≈ 0` and `dσ₁` is the elastic
    increment. But `dα₁ = (2/3) h dγ b`, and with `h → ∞`,
    `h·dγ → (2G n:dε_dev − K dε_v n:r)/((2/3) p b:n)`. That is the consistency
    motion `dα ≈ Δs/p` evaluated at **p = 0.0101**, the start of the substep,
    while the increment takes p to 0.20;
  - stage 2 has `(α + dα₁ − α_in):n = −13`, so `h < 0` and `dγ < 0`. This is the
    `dγ < 0` branch (`:1903-1910`): elastic stress plus the "drag"
    `Δα = Δ(s/p)`;
  - error: both stages' `dσ` are the elastic increment, so `‖dσ₂−dσ₁‖ ≈ 0`, and
    **the error never looks at α** (`:1966-1973`). α/α^b goes 0 → 5.64 in the
    substep.

Counterfactuals, same chain (port, one switch each):

| variant | max α/α^b | first escape | refused incr. |
|---|---|---|---|
| as built (= C++) | 5.18 | k = 20 | 0 |
| `dγ<0` stage leaves α alone (drag frozen) | 5.73 | k = 20 | 27 |
| forced accept at dT_min REFUSES | 5.18 | k = 20 | 0 |
| `Stress_Correction` off | 2.4e7 | k = 4 | 0 |
| **α also in the substep error** | **0.475** | none | 0 |
| TolE 1e-8 (α-blind) | 0.617 | none | 44 (cap) |
| C++, same total strain at δ = 1e-5 / 1e-6 | 1.19 / 0.57 | k = 237 / none | — |

The table verdicts, candidate by candidate:

- **(i) overshoot**: confirmed, in this precise form. The Euler stage of the α
  law, which scales as 1/p, is evaluated at the floor pressure across a substep
  that multiplies p by about 20. The stress-only error test then lets it through.
- **(ii) forced accept**: not the mechanism. Refusing it changes nothing.
- **(iii) `Stress_Correction` / `p_min`**: supplies the **launch state**
  (`α = 0`, and `α_in = 0` after the next reversal). Turning correction off is far
  worse, so it is not the jump.
- **(iv) the `dγ<0` "drag" branch**: not the jump either. Freezing it still
  escapes.
- **(v) finding E**, the α-blind error test: the enabling defect. Measuring α
  alone removes the escape.

### 2.3 Smallest reproducer and threshold (`q1_threshold.py`)

Start from `σ = p_s·I, α = α_in = z = 0` and apply one `dε = (0, δ, 0, 0, 0, 0)`.
C++ ME(campaign) versus the port with α in the error:

| p_s | δ | ME: substeps / η / α/α^b / f | Δp/p | α-aware: substeps / α/α^b |
|---|---|---|---|---|
| 0.0101 | 1e-5 | 1 / 1.19 / 0.80 / 2e-12 | 2.9 | 108 / 0.145 |
| 0.0101 | 3e-5 | 1 / 3.20 / **2.15** / 2e-9 | 6.7 | 152 / 0.204 |
| 0.0101 | 1e-4 | 1 / 10.79 / **5.14** / 2e-8 | 19.8 | 200 / 0.266 |
| 0.0101 | 3e-4 | 1 / 0.84 / **16.4** / **11.2** | 57.7 | 242 / 0.315 |
| 0.1 | 1e-4 | 1 / 3.38 / **1.62** / 8e-8 | 6.8 | 154 / 0.209 |
| 1 | 3e-4 | 1 / 3.14 / **2.16** / −3e-8 | 5.2 | 152 / 0.209 |
| 5 | 3e-4 | 1 / 1.54 / 0.76 / −9e-8 | 3.3 | 119 / 0.166 |

The escape tracks **Δp/p per increment**. α/α^b exceeds 1 once Δp/p lies between
3.3 (0.76) and 5.2 (2.16), whatever p_s is, and the increment is always one
accepted substep. RK45 (C++) independently agrees with the α-aware port: 0.145 / 0.200 /
0.250 / 0.286 on the four floor rows (`q1e_alpha_error.txt`), which confirms the
one-substep ME answer is wrong. Even below the escape threshold the α-blind answer
is 3–5× off in α. The **ring signature** (`q1_attrib`/q3 diagnostics):

- b8 rows with α outside `α^b_θ`: 1950/3 (6.26), 1950/2 (5.86), 1859/2 (1.03),
  1898/4 (1.01). All four have `α_in ≡ 0` exactly.
- The one other `α_in ≡ 0` row, 1950/4 (0.92), is just inside.
- No row with `α_in ≠ 0` is outside. b16: no row outside.

An `α_in` of exactly zero is only produced by a reversal *after* α was zeroed, and
only the low-p resets zero α (`Stress_Correction :3008-3009`,
`explicit_integrator`'s `p_n < p_r` reset).

### 2.4 Finding E — α-aware error control, measured (`q1e_alpha_error.py`)

| chain (p₀ 2, cycles) | ME: max α/α^b | RK45/1e-4 | RK45/1e-7 | ME + α error (port) |
|---|---|---|---|---|
| vertUnload δ 1e-4 | **5.22** (115k substeps) | 0.88 (21 clamp, 67 f>0) | 0.88 (20 clamp, 26 f>0) | **0.475** (798k substeps, 0 f>0) |
| extShear δ 1e-4 | **7.07** (142k) | 0.75 (41 clamp, 67 f>0) | 0.93 (63 clamp, 15 f>0) | **0.509** (863k, 0) |
| shear δ 1e-4 | 0.69 (19k) | 0.69 (90 f>0) | 0.76 (89 clamp) | 0.70 (19k) |
| vertUnload δ 1e-5 | 1.19 (76k) | **9541** (888 f>0) | 0.89 (366 clamp, 300 f>0) | **0.574** (686k) |

Ring rows × 8 probes (640 increments), new escapes (α inside → outside in one
increment) and `f > 1e-6` returned with `rc = 0`:

| | ME | RK45/1e-4 | RK45/1e-7 | ME + α error |
|---|---|---|---|---|
| inside → outside | 9 | 10 | 6 | **1** |
| f > 1e-6 commits | 4 | 512 | 66 | **0** |

- **Yes: α in the error keeps α inside where ME lets it escape.** It does so in
  the port, which differs from ME only in that one line. RK45 keeps it inside only
  by **clamping**: its `dT_min` is hard-coded 1e-3 (`:2249`) and it force-accepts
  with the same Mc clamp (`:2494-2508`).
- RK45 also runs **without `Stress_Correction`** (commented out, `:2521`), hence
  hundreds of `f > 0` commits. It blows up (α/α^b 9541) on the δ 1e-5 chain at
  TolR 1e-4. RK45 as it stands is not a candidate.
- The cost of α in the error is real on the floor-bouncing chains: 6–9× the
  substeps (§6, q5).

A milder, slower excursion also exists: vertUnload δ 1e-5 reaches 1.19 at
p ≈ 80 kPa after 223 increments. That chain is round-off chaotic (C++ vs port
diverge by 0.9 in α/α^b over 900 increments), so it is **not attributed**. The
α-aware port holds it to 0.574.

### 2.5 Candidates F and G (orchestrator's survey), and the ranking (`q1fg_candidates.py`)

- **F — negative-denominator misclassification.** `temp4 = Kp + 2G(B − C tr n³) −
  K D n:r < 0` makes a LOADING stage (numerator `n:Ce:dε > 0`) take the `dγ < 0`
  elastic+drag branch. On top of that, `q = max(0.8√(TolE/err), 0.5)` has no upper
  cap after an `err = 0` acceptance (`:2038`).
- **G — α_in staleness.** α_in is re-seated only once per global increment
  (`integrate()`, `:1049-1054`), so `(α − α_in):n` can go through 0 or negative
  inside the substeps. `h = b0/((α−α_in):n)` (`:5398-5401`) then blows up or
  flips sign.

The stage-by-stage signs, logged by the port:

| case | stage 1 | stage 2 |
|---|---|---|
| reproducer (floor, `dε_yy` 1e-4) | plastic: temp4 + (1.2e8), aain **= 0** (h = 1e10) | drag: temp4 + , numerator **−**, aain **−13** (h < 0, Kp < 0) |
| 1950/3 `shear+ 1e-5` (the f > 0 case) | plastic: temp4 +, aain = 0 | drag: temp4 +, numerator −, aain −0.16 |
| 1950/3 `isoExt 1e-6` (f = 0.118) | **drag by F**: temp4 **−1.9e10**, numerator + (loading), aain = 0 | **F** again: temp4 −2e10, numerator + |
| 1950/3 committed state | Kp = −368 (b:n = −8.19 because α is already outside), h = +192 | |

Every accepted substep that carried α/α^b across 1 on the two reversal chains,
classified by the signs of its two stages:

| chain | stage pattern | count |
|---|---|---|
| vertUnload | plastic(aain = 0) + drag(temp4 +, num −, aain −) | 3 (all) |
| extShear | plastic(aain = 0) + drag(temp4 +, num −, aain −) | 5 |
| extShear | plastic(aain = 0) + drag(**temp4 −, num +**, aain −), i.e. **F** | 9 |
| extShear | drag(temp4 −, num +, aain −) + plastic, i.e. **F** in stage 1 | 8 |
| extShear | plastic + plastic, with aain − in at least one stage (G, both stages "plastic" but h < 0) | 10 |
| extShear | plastic + plastic, aain + in both | 1 |
| extShear | drag + drag (num −: genuine unloading classification) | 2 |

Every crossing happened in an increment whose α_in was reset at its start, and
**37 of the 38 crossing substeps (3 vertUnload + 35 extShear) have `(α − α_in):n ≤ 0` in at least one stage**.
The exception is one small substep (dT 0.0017) in a chain already driven far out.

Counterfactuals, one switch each (port, same chains + the reproducer):

| variant | vertUnload | extShear | reproducer α/α^b | reproducer substeps |
|---|---|---|---|---|
| as built | 5.18 | 7.01 | 5.14 | 1 |
| **G-fix: h = b0/⟨(α−α_in):n⟩** (1e10 when ≤ 0) | **0.59** | **0.93** | **0.55** | 1 |
| F-q: cap q ≤ 2 | 5.18 | 7.01 | 5.14 | 1 |
| F-drag: `dγ<0` stage leaves α alone | 5.73 | 7.52 | 4.90 | 1 |
| G-fix + F-drag | 0.59 | 0.93 | 0.55 | 1 |
| **E-fix: α in the substep error** | **0.475** | **0.509** | **0.266** (≈ RK45 0.25) | 200 |

**Ranking of the B candidates** (which of them, removed alone, removes the escape):

| candidate | reproduces α outside? | smallest reproducer | verdict |
|---|---|---|---|
| **G — stale α_in, h < 0 / h = 1e10 mid-increment** | **yes, it is the trigger.** 37/38 crossing substeps have `(α−α_in):n ≤ 0`, and a Macaulay bracket on it alone keeps α inside (0.55–0.93) at no substep cost | §2.3 floor reproducer: `σ = 0.0101 I, α = α_in = z = 0, dε_yy = 1e-4` | primary **cause** |
| **E — stress-only error test** | **yes, it is the enabler.** It is why the G-corrupted stages are accepted at `dT = 1`. Adding α to the error alone keeps α inside and gives the right answer (0.27, = RK45) | same | primary **fix** (accuracy) |
| overshoot of the α law at the start-pressure (i) | the mechanism by which G's stage-1 moves α (`h·dγ` finite, evaluated at p_start) | same | the arithmetic of G+E, not separate |
| **F — negative denominator → elastic drag** | **contributes, not necessary.** It produces 17 of the 35 extShear crossings and the 1950/3 `isoExt` f > 0, but freezing the drag does NOT remove the escape (5.7 / 7.5). It also arises on its own once α is outside: b:n < 0 → Kp < 0 → temp4 < 0 | 1950/3 `isoExt 1e-6`: loading misread as elastic in both stages, err 0, f = 0.118 | a **consequence** that then compounds |
| F — uncapped q | no effect: the escapes already happen at `dT = 1` or `0.1`; capping q ≤ 2 changes nothing | — | not a mechanism here |
| dT_min forced accept (C) | no. Refusing it changes nothing on the escape chains (§2.2). It teleports η back to Mc on inadmissible states (§4) | — | separate defect |
| `Stress_Correction` silent exit | not the escape. It is the reason a misclassified increment returns **f > 0** as success (§3). Its low-p branch supplies the α = 0 launch state | 1950/3 `shear+ 1e-5` (f = 0.0139, rc 0) | separate defect (Q2) |
| start state with f > TolF not corrected before substepping | not observed as a trigger: all escapes start on or inside the cone (f ≤ 2e-8) | — | not reproduced |

What this means for WP-129: G and E are both needed.

- **G's fix is the model-level one.** Dafalias–Manzari's `α_in` is the
  back-stress at the last reversal, so `(α − α_in):n ≥ 0` must hold. Either
  re-seat α_in inside the substep loop when it goes negative, or at least
  Macaulay-bracket h.
- **Alone, G's fix is cheap but inaccurate.** It returns 0.55 on the reproducer
  where the α-aware answer is 0.27.
- **E's fix gives accuracy and the substep control that exposes the stiffness.**

## 3. Q2 — why an outside-the-yield-surface point returns success

At b8 1950/3, `shear+ 1e-5` (γ_xy):

- **The path is unload-then-plastic (4).** The elastic predictor is `G·γ ≈ 0.048 kPa`
  of shear, against a cone radius of `√(2/3)·m·p = 0.0014 kPa`. It crosses the cone
  and exits the far side. `IntersectionFactor_Unloading` finds the exit at 3.5 % of
  the increment, and ModifiedEuler runs the remaining 96.5 % from there.
- **The reversal sets `α_in := α`.** The trial direction opposes `α − α_in`, so
  stage 1 has `(α − α_in):n = 0` and `h = 1e10`. Stage 2 has
  `(α − α_in):n = −0.16`, so `h = −1.2e4` and `dγ < 0`: the elastic drag. Both
  `dσ` equal the elastic increment, the error is 4e-9, and the substep is
  **accepted at dT = 1**. α moves by the average of the two stages, which does not
  lie on the cone around the new stress: `f = +0.0139`.
- **`Stress_Correction` gives up silently.** It is called on the accepted substep.
  `|fr| > TolF`, so it enters the correction loop.
  - The first attempt, λ along `C:R` with α along `(2/3)h b`, does not reduce
    `|f|`. Here `h < 0` and `b:n` is inconsistent.
  - The fallback, λ along `∂f/∂σ = n − (n:r)/3·I`, does not either. At η ≈ 13,
    `n:r` is O(10), so the linear step mostly moves p, and the tiny cone makes the
    linearisation useless.
  - It then takes **`return;` without touching NextStress/NextAlpha**
    (`:3054-3058`, the "Couldn't decrease the yield function" message behind
    `debugFlag`). Port counter: `corrGiveUp = 1`.
- **Nothing downstream checks f.** `ModifiedEuler` sets `T = 1` and returns.
  `explicit_integrator` and `integrate` return nothing, and the wrapper's status
  sees no cap. **rc = 0 with f = 0.0139**, about 10× the cone radius.

The same path with `isoExt 1e-6` gives `f = 0.118` with err = 0. Starting from the
floor with `δ = 3e-4` gives `f = 11.2 kPa` (§2.3). On the ring, 4 of 640 ME
increments return `f > 1e-6` as success (q1e).

**This is a defect under the fork's own checklist.**
`.claude/skills/ladruno-new-material/SKILL.md`: "A non-converged return map must
FAIL (return < 0), never commit `f > 0` as success". It is a vanilla defect, since
the branch is upstream code. Fix direction: the give-up and the i == maxIter "Still
outside" branch must raise a refusal flag that ModifiedEuler turns into a failed
update. Quirks row added.

## 4. Q3 — F18(e): can ANY integrator take the b8 worst point?

The state is b8 el 1950 gp 3: `p' = 0.352, η = 12.87, f = 7.5e-10` (on its cone).
Its admissibility:

| quantity | value | admissible? |
|---|---|---|
| √(3/2)‖α‖ / α^b_θ (model's own, θ from n) | **6.26** (α^b_θ = 2.056, M^b_θ = 2.061) | no — b:n = −8.19 |
| (α − α_in):n | 9.87, α_in ≡ 0 | the §2.3 signature |
| tr(α), tr(z) | −1e-10, −6e-10 (‖z‖ = 10.5, zmax 12.5) | yes (round-off) |
| 1950/2 | α/α^b 5.86, tr(z) = −3.8e-4 on ‖z‖ 9.7 | no |
| 1859/2 | α/α^b 1.03, **tr(α) = 2.3e-3** | trace defect (WP-127 minor) |

Replay of 1950/3, probes ± isotropic and ± simple shear at 1e-7 / 1e-6 / 1e-5
(`q3_worst_point.py` → `q3_table.py`). Each cell reads substeps / returned η /
returned f. The flags mean:

- **C**: the dT_min forced accept **with the Mc clamp** fired (η teleported to
  about Mc);
- **+**: `f > 1e-6` returned with rc = 0;
- RK45 substeps are n/a because RK45 is not instrumented;
- a CPPM cell with substeps > 0 means CPPM failed and fell back to ME.

| probe | ME(campaign) | ME(1e-8,honor) | RK45(1e-10) | CPPM(1e-10) | ME+α err (port) |
|---|---|---|---|---|---|
| isoComp 1e-7 | 1 / 12.8 / 8e-10 | 1 / 12.8 / 8e-10 | 1.32 | 0 / 12.9 / −2e-17 | 1 / 12.8 / 8e-10 |
| isoExt 1e-7 | 1 / 13.0 / 1e-11 | 1 / 13.0 / 1e-11 | 13.0 / 0.024 + | 1 / 13.0 / 1e-11 | 7 / 13.0 / 8e-10 |
| shear− 1e-7 | 1 / 12.9 / 3e-10 | 58 / 12.9 / 8e-10 | 1.33 C | 0 / 12.9 / 7e-13 | 1 / 12.9 / 3e-10 |
| isoComp 1e-6 | 1 / 12.1 / 8e-10 | 1 / 12.1 / 8e-10 | 1.25 | 0 / 12.9 / 3e-16 | 1 / 12.1 / 8e-10 |
| isoExt 1e-6 | 1 / 13.8 / **0.12 +** | 1 / 13.8 / **0.12 +** | 13.7 / 0.23 + | 19 / 1.14 C | 144 / 1.33 C |
| shear+ 1e-6 | 1 / 12.9 / **0.0029 +** | 471 / 1.33 C | 1.32 | 4404 / 1.33 C | 215 / 12.8 / 2e-10 |
| shear− 1e-6 | 8 / 12.9 / 5e-12 | 225 / 12.9 / 8e-10 | 1.33 C | 0 / 12.9 / 1e-13 | 8 / 12.9 / 5e-12 |
| isoComp 1e-5 | 1 / 7.86 / 1e-9 | 1 / 7.86 / 1e-9 | 0.79 | 1 / 7.86 / 1e-9 | 1 / 7.86 / 1e-9 |
| isoExt 1e-5 | 356 / 1.33 C | 843 / 1.33 C | 27.8 / **2.0 +** | 6657 / 1.33 C | 404 / 1.33 C |
| shear+ 1e-5 | 1 / 12.7 / **0.014 +** | 1 / 12.7 / **0.014 +** | 1.24 | 5780 / 1.33 C | 388 / 1.33 C |
| shear− 1e-5 | 59 / 12.9 / 7e-9 | 289 / 12.9 / 8e-10 | 1.33 C | 2420 / 12.9 / 8e-10 | 59 / 12.9 / 7e-9 |

(`shear+ 1e-7` is elastic for all: f = −4.8e-4. All rc = 0: nothing refuses.)

**Plainly: no integrator takes this point.**

- **Loading probes** (isoComp, shear−) are "taken" by every scheme except RK45.
  The stress stays on the cone around the inadmissible α (η 12–13, or 7.9 when the
  compression raises p). The state stays inadmissible and no bounding-surface
  physics acts.
- **Unloading probes** (isoExt, shear+) end one of two ways. ME(campaign) and ME
  at 1e-8 commit **f > 0 as success** (§3). Everything else, including ME at 1e-8,
  CPPM (via its ME fallback), and the α-aware port, runs into dT_min and
  **teleports η 12.9 → 1.33** through the forced-accept Mc clamp. That is the
  clamp's arithmetic, not an integration of the model.
- **RK45** teleports on 7 of 11 non-zero probes, even at `1e-7`, because its
  dT_min is 1e-3.
- **CPPM converges** (0 substeps, `f ≈ 1e-16`) on the small loading probes only.
  It returns the same inadmissible α.

**What an integrator should do with an inadmissible committed state: refuse, don't
project.**

- The state cannot be reached by the continuous model. α is 6× past a bounding
  surface that only lets α exceed it by the contraction of M^b(ψ) in a step. It
  also carries a stale `α_in ≡ 0`. Every downstream quantity (h, Kp, b:n, D)
  answers garbage there.
- A projection (α back to √(2/3)·α^b·n, trace-free α/z) would silently rewrite
  committed history. It would also have to be defined consistently with σ, since
  f must stay 0, and it would hide the upstream defect that produced the state.
- The fix belongs **before** the commit (§6 (a)–(c)), so these states are never
  produced. The entry check is the backstop:
  - refuse a trial whose committed state has `α/α^b_θ > 1 + κ`. κ must admit the
    legitimate post-peak `b:n < 0` softening, about 1.03 here on 1859/2; κ = 0.25
    would separate the 1.01–1.03 rows from the 5.9–6.3 ones;
  - or refuse when `|tr α|, |tr z| > 1e-6·max(‖·‖, m)`.

  In both cases refuse with a code the WP-99 roster forwards, so the global step
  is cut, and log it once per point.
- The replay tool already projects traces and warns. That is right for a
  **diagnostic** tool and wrong for the integrator.

## 5. Q4 — the F18(a) baseline

Scripts: `q4_baseline.py`, `q4_stiffness.py`, `q4_refcheck.py`. In this section,
error means `‖σ_integrator − σ_reference‖` (contravariant Voigt norm, kPa) for
**one** increment from the same committed state.

The reference is C++ ME at TolR 1e-8 with `-honorTolR 1`. A case is excluded if
that reference force-accepted, capped, or disagrees with ME at 1e-9 by more than
10 % of the campaign error. The two obvious alternatives fail:

- RK45 is unusable as a reference: its dT_min is 1e-3 and it clamps (§4).
- An α-aware tight reference is too slow in Python on the 1e-5 ring set. §5.4
  checks the α-blind reference against it on the 1e-6 set instead.

### 5.1 Today's norm is already a 1 kPa floor

`:1966-1973`: `err = ‖dσ₂−dσ₁‖` if `‖σ‖ < 0.5`, else `/(2‖σ‖)`. That is
`‖dσ₂−dσ₁‖ / max(2‖σ‖, 1 kPa)`, continuous at `‖σ‖ = 0.5`. The port with
`σ_ref = 1` reproduces today **bit for bit on 640/640** ring increments.
`‖σ‖` is the full stress norm, hydrostatic part included, so `‖σ‖ ≥ √3·p`. The
consequences:

- `σ_ref ≤ 1` acts only where `‖σ‖ < 0.5 kPa`, i.e. `p ≲ 0.29`;
- `σ_ref = 5` acts only where `p ≲ 1.4`;
- `σ_ref` must exceed `2‖σ‖ ≈ 12 kPa` before it touches the act's
  p' ≈ 3.5 kPa ring (§1.3 of the intake).

The intake's "at p' ≈ 3.5 kPa that asks for an absolute stress error of about
1e-3" is right. But the floor proposed at `σ_ref ∈ {0.1, 1, 5}` does not move it.

### 5.2 Substeps and error per increment

Today = C++. The σ_ref columns are the port with the proposed norm.

| set | today substeps med/p95/max | today error kPa med/p95/max | σ_ref = 0 / 0.1 / 1 | σ_ref = 5 | σ_ref = 20 |
|---|---|---|---|---|---|
| ring 80×8 (523 with reference) | 5 / 71 / 1336 | 1.1e-4 / 2.9e-2 / 0.84 | = today (max 1340) | 5 / 70 / 860; err p95 3.6e-2 | 4 / 63 / 223; err p95 4.0e-2, max 0.74 |
| ring @1e-6 (274) | 1 / 12 / 138 | 6.6e-5 / 1.8e-3 / 4.9e-2 | = today | max 128 | max 30; err med 1.2e-4 |
| ring @1e-5 (249) | 15 / 96 / 1336 | 1.8e-4 / 0.19 / 0.84 | = today | max 860 | max 223 |
| constant-p p₀ = 2, to η/M^b = 1.004 (599 incr.) | 51 / 52 / 52 | 1.1e-4 / 1.5e-4 / 3.0e-4 | = today | = today | 50 / 51 / 51; err med **2.8e-4** |
| constant-p p₀ = 5, to 0.985 (600) | 34 / 34 / 36 | 2.1e-4 / 2.6e-4 / 0.16 | = today | = today | = today |
| constant-p p₀ = 20, to 0.945 (595) | 18 / 19 / 20 | 8.9e-4 / 1.9e-3 / 2.2e-3 | = today | = today | = today |
| `active` path p₀ 2 / 5 / 20 (p rises to 12.5 / 18 / 38) | 28/46/47, 23/31/36, 15/17/19 | 2.8e-4, 4.6e-4, 1.3e-3 (med) | = today | = today | = today |

Relative error on the smooth chains is 2e-5 per 1e-5 increment, i.e. TolE 1e-4
times 0.2.

### 5.3 The cost is stability-limited — a norm cannot buy it back (`q4_stiffness.txt`)

Substeps for the peak increment of each constant-p chain, by TolE (rows) and
σ_ref (columns):

| TolE | p₀ 2: σ_ref 1 / 20 / 100 / 1e4 | p₀ 20: σ_ref 1 / 20 / 100 / 1e4 |
|---|---|---|
| 1e-6 | 60 / 57 / 55 / 41 | 26 / 26 / 26 / 19 |
| 1e-5 | 54 / 53 / 51 / 8 | 20 / 20 / 20 / 14 |
| **1e-4 (today)** | **52** / 50 / 41 / 1 | **19** / 19 / 19 / 1 |
| 1e-3 | 42 / 31 / 8 / 1 | 15 / 15 / 14 / 1 |
| 1e-2 | 8 / 5 / 1 / 1 | 1 / 1 / 1 / 1 |

An accuracy-limited Heun would scale with TolE^−½, i.e. 10× from 1e-6 to 1e-4.
The measured change is 1.15× at p = 2 and 1.4× at p = 20: the substep size sits on
the explicit **stability** limit, which scales as 1/√p (51 / 34 / 18 substeps at
p = 2 / 5 / 20). Loosening the norm does nothing until it is so loose (σ_ref ≥ 100,
TolE ≥ 1e-3) that the estimator stops seeing the instability. Past that point the
answer is whatever one unstable step returns.

On the ring the tail behaves the same way: the max goes 1428 → 1336 from TolE
1e-6 to 1e-4. Only the median is accuracy-like: 72 / 27 / 14.

### 5.4 Caveat: the reference is α-blind (`q4_refcheck.txt`)

On the ring @1e-6 set, compare the C++ ME 1e-8 reference (ref A) with an α+z-aware
port at 1e-8 (ref B):

| comparison (313 cases, 7 excluded because a reference forced/refused) | median kPa | p95 kPa | max kPa |
|---|---|---|---|
| \|ref A − ref B\| | 5.6e-10 | 6.1e-3 | 7.0e-2 |
| \|today − ref A\| (what §5.2 reports) | 5.6e-5 | 5.1e-3 | 7.0e-2 |
| \|today − ref B\| | 1.1e-4 | 1.3e-2 | 7.8e-2 |

On most states the two references agree to round-off. On the tail (p95) the α-blind
reference sits 6e-3 kPa from the α-aware one, and today's error against the α-aware
reference is about 2× (median) to 2.6× (p95) what §5.2 reports.

The ring errors in §5.2 therefore **understate** the true integration error, by
about 2–3× on the 1e-6 set: they
measure ME against a tighter ME that makes the same α-blind acceptances. The
constant-p numbers are unaffected. There h is finite and both stages are plastic,
so the stress error does couple to α.

### 5.5 Against the global Newton's tolerance

The act runs `NormUnbalance` at `1e-5 × ‖reference load‖`. One Gauss point of a
B/8 quad (h = 0.1875 m, 2×2 rule) with a stress error δσ contributes about
`‖r‖ ≈ 0.4·h·δσ ≈ 0.08·δσ` kN/m. The act has to say which vector is "the
reference load"; the two ends:

| reference load | Newton tol (kN/m) | δσ at one GP that uses the whole tolerance |
|---|---|---|
| footing weight 18.4 kN/m | 1.8e-4 | **≈ 2e-3 kPa** |
| bearing load at the wall, ≈ 650 kPa × 1.5 m | 9.8e-3 | ≈ 0.12 kPa |

Where the proposed floor's admitted error is below the Newton tolerance:

- **Holds** on every constant-p / `active` state at p ≥ 2 kPa, today and at
  σ_ref = 20 alike. Those errors are 1e-4–2e-3 kPa per increment, below even the
  tight end.
- **Fails** on the ring tail: today p95 2.9e-2, max 0.84 kPa, and σ_ref = 20
  gives max 0.74. That tail exceeds the tight end at p95 and the loose end at max,
  and §5.4 says the true tail is larger.
- **Does not matter** in any case: the floor is inert where the error is small,
  and it cannot touch the tail. The tail is α-blind acceptance and
  stability-limited substepping, not the norm.

### 5.6 What WP-129 must measure once `-errFloor` exists

1. C++ byte-identity at the default: the default `σ_ref = 1`, **not 0**. §5.1: 0 is
   *not* today.
2. The port's predictions reproduced in C++: this section's ring and constant-p
   tables with `-errFloor 0 / 0.1 / 5 / 20`.
3. The same tables with α (and z) in the error. That, not the floor, is where the
   ring numbers move.

## 6. Q5 — recommendation for WP-129

**Verdict on F18(a): wrong fix.** It is not harmful: it is inert at σ_ref ≤ 5 on
every state measured, and at σ_ref = 20 it saves 2 of 52 substeps on the
constant-p peak while doubling the error. It targets neither the cost, which is
stability-limited (§5.3), nor the failure, which is α-blind acceptance plus
silent `f > 0` plus the Mc-clamp teleport (§2–§4). If the owner still wants the
flag for the act, ship it byte-identical with **default 1 kPa** and document §5.1.

**The minimal α-stability fix**, in the order that matters. Measured on the port
(`q5_error_variants.txt`):

| variant | reversal chains: max α/α^b (vertUnload / extShear) | substeps (vU / eS) | forced / refused (chains) | ring 640: substeps med/p95/max; in→out; forced; refused; f>1e-6 | constant-p, one increment per committed state: substeps med (p₀ 2 / 20); max \|Δσ\| vs today |
|---|---|---|---|---|---|
| today | 5.18 / 7.01 | 115k / 141k | 2 / 0 | 4/81/1336; 10; 2; 0; **5** | 51 / 18; — |
| + α in error | **0.475 / 0.509** | 798k / 863k | 0 / 18 (all cap hits) | 15/120/1369; **1**; **5**; 0; **0** | 51 / 18 (max 23); 2e-4 / **0.27** kPa |
| + α + z in error | identical to + α | | | identical | identical |
| + α + z, refuse at dT_min | identical to + α | | 0 / 18 | 15/120/1369; 1; 5 → **5 refused**; 0 | identical |

Reading it:

- **Adding z to the error changes nothing measured.** The fabric only moves when
  `D < 0` (dilation), and it is already bounded by `zmax`.
- **α alone removes 9 of 10 ring escapes and every `f > 0`.** It leaves 5 Mc-clamp
  teleports on the ring. Only refusing at dT_min turns those into refusals.
- **At p₀ = 20, α in the error moved one smooth-chain increment by 0.27 kPa.**
  Today's α-blind acceptance was not harmless there either.
- **The 18 refusals on extShear are `-maxSubsteps` 20000 cap hits.** On
  floor-bouncing paths the honest cost exceeds the campaign cap.

1. **(a) α and z in ModifiedEuler's substep error, Sloan–Abbo–Sheng style.** Use
   `err = max(‖dσ₂−dσ₁‖/max(2‖σ‖,1), ‖dα₂−dα₁‖/max(2‖α‖,1), ‖dz₂−dz₁‖/max(2‖z‖,1))`,
   RK45's form. This is the one change that removes the escapes and all `f > 0`
   commits on every path measured. **Alone it does not suffice.** It still meets
   dT_min on inadmissible or near-singular states and then force-accepts with the
   Mc clamp: 1950/3 isoExt / shear+ (§4); and see the refused column.
2. **(b) Refuse instead of force-accepting at dT_min** (finding C): return a named
   refusal, counted by the existing `forcedAtDTmin`. Without (a) this refuses
   nothing on the escape paths, where the escape happens at dT = 1, so it is not a
   substitute. With (a) it converts the remaining teleports into step cuts.
2b. **(a') Keep `(α − α_in):n ≥ 0` inside the substeps (candidate G).** Re-seat
   α_in at the substep where it goes negative, which is the model's definition of a
   reversal, or at minimum use `h = b0/⟨(α−α_in):n⟩`. Measured alone, the
   Macaulay form keeps α inside (0.55–0.93) at zero substep cost but is
   inaccurate (0.55 vs 0.27). Ship it **with** (a), not instead of it. Test: T1/T2
   with (a) disabled must still show α/α^b ≤ 1.
3. **(c) `Stress_Correction`'s give-up and "still outside" branches must fail**
   the update (Q2). Cheap, vanilla-additive (`// Ladruno`), and the checklist
   requires it.
4. **(d) An admissibility check on entry** (§4): refuse a committed state with
   α/α^b_θ > 1 + κ or a non-zero trace. It is a backstop and cannot fire once
   (a)–(c) hold.
5. **Not needed as separate fixes**, measured:
   - capping `h·Λ` per substep, and the `1e10` sentinel: with α in the error they
     only cost substeps;
   - freezing the `dγ<0` drag: it does not change the escape (§2.2).

   A hard α-bound (projection) is **not** recommended. It hides (a)'s failures
   instead of refusing them.

**Cost warning for a Sloan–Abbo–Sheng direction.** On the floor-bouncing chains,
(a) costs 6–9× the substeps (115k → 798k). That is the stiffness the α-blind test
was hiding. On smooth constant-p states the cost is small (table). An α-aware
*explicit* scheme is still stability-limited at low p (§5.3), so if the ring cost
matters after correctness is restored, the lever is implicit or stiffly-stable
stages (F18(c) CPPM), not the norm.

**Tests that would prove the fix (C++, WP-129):**

- **T1 (reproducer).** `ladrunoSANISANDReplay` from `σ = 0.0101·I, α = α_in = z = 0`,
  `dε_yy = 1e-4`. Today it returns α/α^b 5.14 in 1 substep. Fixed, it must return
  α/α^b < 1 (port: 0.266) or refuse. Parametrise over q1_threshold's grid: no
  α/α^b > 1 with rc = 0 anywhere.
- **T2 (chains).** q1_attrib's vertUnload / extShear δ 1e-4 cycles: max α/α^b ≤ 1
  + κ, zero `f > 1e-6` with rc = 0, and every refusal counted.
- **T3 (ring).** 80 rows × 8 probes: zero `f > 1e-6` with rc = 0, zero Mc-clamp
  teleports (a `forcedClampMc` > 0 must now be rc ≠ 0), and 1950/3's unloading
  probes refused.
- **T4 (Q2).** 1950/3 `shear+ 1e-5` must not return rc 0 with f = 0.0139.
- **T5 (no regression).** WP-127's byte-identity decks change only where the new
  error test binds. Report which. The constant-p chains' substeps and errors stay
  within the q5 port prediction.

## 7. Side note — the two pre-existing flip-determinism failures

The failing tests are
`tests/test_ladruno_sanisand_flip_determinism.py::test_first_ten_push_steps_bit_identical_across_mkl_threads`
and `::test_default_first_step_is_immune_to_the_hold_lottery`. Time-boxed look
(`flipdet_probe.py`, the test's own child deck, 1 thread):

| variant | first push step |
|---|---|
| as tested (Pardiso, NormDispIncr 1e-8, 100 it.) | **−3**, "failed to converge", load factor 9.67011 |
| `system FullGeneral` instead of Pardiso | **−3**, identical |
| Pardiso, NormDispIncr **1e-6** | converges: 9.731325, 13.148652, …, 35.943014 (10/10) |

- **Not environment.** The failure is solver-independent, so it is not MKL or
  Pardiso threading and not this machine.
- **The numerics moved.** The deck still reaches the F14 state (223/284
  round-off `α − α_in` after the holds), but the converged first step is 9.7313,
  where the test's docstring pinned 9.659111 on `48c0e99bc`. That drift is not a
  tolerance effect. Material commits in `48c0e99bc..c03a1bd4b` include
  `dee04dbe3` (WP-110 F15 tangent fix; this deck runs **TanType 2**, whose Newton
  path that fix changes) and `25136af9c` (WP-112).
- **Not bisected**, since that needs builds. The test's pinned constants and its
  1e-8 / 100-iteration convergence need re-deriving on the current tree. Not
  fixed here.

## 8. Reproduce

The runner is CPython 3.12 with `-S`; `_boot.py` asserts the worktree's pyd.

```
cd Ladruno_files/testbed/sanisand_ring_trace
PY=C:/Users/nmora/AppData/Local/Python/pythoncore-3.12-64/python.exe
$PY -S validate_port.py      # port vs C++  -> out/validate_port.txt            (~3 min)
$PY -S q1_search.py          # Q1 search    -> out/q1_search.txt                (~2 min)
$PY -S q1_attrib.py          # Q1 attribution + counterfactuals                 (~10 min)
$PY -S q1_threshold.py       # smallest reproducer + threshold                  (<1 min)
$PY -S q1_forced.py          # forced-accept counterfactual, delta 1e-5 chain   (~3 min)
$PY -S q1e_alpha_error.py    # finding E: ME vs RK45 vs alpha-aware port        (~8 min)
$PY -S q3_worst_point.py && $PY -S q3_table.py   # F18(e) table                 (<1 min)
$PY -S q4_baseline.py        # F18(a) baseline                                  (~10 min)
$PY -S q4_stiffness.py       # stability-limited substeps                       (<1 min)
$PY -S q4_refcheck.py        # alpha-blind reference caveat                     (~5 min)
$PY -S q5_error_variants.py  # WP-129 variants                                  (~15 min)
$PY -S q1fg_candidates.py    # orchestrator candidates F and G + ranking        (~5 min)
$PY -S flipdet_probe.py      # section 7                                        (~2 min)
```

## 9. Verified / not verified

- **Verified (measured):** everything in the tables above, on WP-127's binary and
  the validated port. The port's agreement with the C++ is 1581/1600 increments to
  5e-11, and the rest is shown to be round-off chaos of the C++ itself.
- **Not verified:**
  - that the act's explicit lane actually took an increment with `Δp/p ≳ 4` at the
    floor for 1950/2-3. The signature matches (`α_in ≡ 0`, η 12–13, p 0.35), but
    the lane's history is not available;
  - the 1.19 slow excursion (round-off chaotic, not attributed);
  - any C++ implementation of the §6 fixes. The cost and effect numbers for them
    are port predictions for WP-129 to reproduce;
  - the Newton-tolerance map's constant 0.08 kN/m per kPa, which is an
    order-of-magnitude element estimate, and which vector the act calls "the
    reference load";
  - the flip-determinism drift's cause (not bisected).
