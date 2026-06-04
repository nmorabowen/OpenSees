---
title: "LadrunoBondSlip — 1D Bond-Slip (τ–s) UniaxialMaterial"
project: Ladruno
type: reference / user guide
status: shipped (v1, MC2010 monotonic backbone + energy channels)
classTag: 33002
material: uniaxialMaterial LadrunoBondSlip
related:
  - "[[20_ladruno_embedded_reinforcement_adr]]"
  - "[[LadrunoRebarBuckling_guide]]"
  - "[[LadrunoUniaxialJ2_guide]]"
  - "[[Ladruno_materials_guide]]"
  - "[[11_brick_asdconcrete_integration]]"
tags:
  - material
  - uniaxial
  - bond-slip
  - bond
  - rebar
  - rc
  - embedded-reinforcement
---

# LadrunoBondSlip — 1D Bond-Slip (τ–s) UniaxialMaterial

> [!abstract] One-line summary
> `LadrunoBondSlip` (classTag **33002**) is a one-dimensional **bond-slip**
> `UniaxialMaterial`: it relates the local **bond stress** $\tau$ between a
> reinforcing bar and the surrounding concrete to the **slip** $s$ (relative
> bar-to-host displacement). It is the *axial constitutive slot* of the
> embedded-reinforcement element [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]],
> giving the fork a **DIANA-style** explicit bond-slip model for pull-out,
> development-length, anchorage, and lap-splice studies in reinforced concrete.
> v1 ships the **fib Model Code 2010 / CEB-FIP monotonic backbone** (cyclic
> degradation → v2).

This document is the **descriptive reference**: the constitutive theory exactly
as coded, the three must-fix numerical subtleties baked into v1, the OpenSees
wiring (material *and* the element that consumes it), and the intended use cases.
For the chronological design record, the validation gate plan, and the wider
embedded-reinforcement architecture (Mode P/T, the inverse-map, the penalty tie)
see [[20_ladruno_embedded_reinforcement_adr]] §D4.

> [!note] Code is the source of truth
> Where this guide and the ADR's *design-stage* signature differ, the **shipped
> source wins**. Two notable shipped deviations from the ADR sketch: (a) the bond
> fracture energy $G_F$ is an **optional `-Gf` flag** that *overrides* $s_3$, not a
> mandatory 8th positional argument; (b) the element converts stress→force with a
> single scalar `-bondScale` $(=\pi d\,L_\text{trib})$ rather than separate
> `-perimeter`/`-ltrib` flags. Both are flagged inline below.

---

## Contents

- [[#1 Intended use cases|1. Intended use cases]]
- [[#2 Why it exists — what drove the implementation|2. Why it exists]]
- [[#3 The bond-slip backbone (MC2010)|3. The backbone]]
- [[#4 The three must-fix subtleties|4. Must-fix subtleties]]
- [[#5 State update, unload/reload, and tangent|5. State update & tangent]]
- [[#6 The unit contract — stress to force|6. Unit contract]]
- [[#7 Command syntax — material|7. Command syntax (material)]]
- [[#8 Integration with the embedded-rebar element|8. The embedded element]]
- [[#9 Recording and energy channels|9. Recording & energy]]
- [[#10 Verification and validation|10. Verification]]
- [[#11 Limitations and boundaries|11. Limitations]]
- [[#12 References|12. References]]

---

## 1. Intended use cases

`LadrunoBondSlip` is built for **reinforced-concrete modeling where bar–concrete
slip matters** — i.e. you cannot assume the rebar is perfectly bonded into the
concrete bulk. The canonical targets:

| Use case | Why LadrunoBondSlip |
|---|---|
| **Pull-out tests** | The τ–s backbone *is* the pull-out constitutive law; a single embedded bar driven under displacement control reproduces the CEB-FIP pull-out response directly. |
| **Development-length / anchorage** | The bond stress integrated over the embedded length gives the developable bar force; the softening + residual branches capture loss of anchorage. |
| **Lap splices** | Two overlapping bars sharing a bond interface, each with its own τ–s law; well-confined splices are in scope (cover-controlled splitting is **not** — see §11). |
| **Tension stiffening** | The slip between cracks redistributes force between bar and concrete; the bond law sets how stiff that exchange is. |

The model it implements is the **DIANA-style** bond-slip reinforcement: an
**own-DOF rebar element** plus a **1D τ–s interface law**, the natural choice
when you need *explicit slip* rather than a tied (perfect-bond) embedding. The
contrast is the **perfect-bond** path (`-perfect` on the element, or the deferred
Mode T), which ties bar translations to the host and admits **no slip at all**.

**It is NOT:**
- **A steel material.** It carries *bond* stress (interface shear per unit bar
  surface), not *axial* bar stress. The bar's own σ–ε comes from a separate
  steel material (e.g. [[LadrunoUniaxialJ2_guide|LadrunoUniaxialJ2]] /
  `Steel02`) on the rebar element itself.
- **A rebar-buckling law.** That is the compression-side geometric overlay
  [[LadrunoRebarBuckling_guide|LadrunoRebarBuckling]] — a different axis of RC
  fiber behaviour, applied to the bar material, not the bond interface.
- **A splitting-bond model.** Cover-controlled splitting failure is out of scope
  for v1 (deferred to v1.1); see §11.

## 2. Why it exists — what drove the implementation

Embedded reinforcement in OpenSees historically meant **perfect bond**: tie the
bar's translations to the host element's displacement field and let the bar
carry whatever force the strain demands. That is correct for *well-anchored*
bars far from cracks, but it **over-smears** cracking and cannot represent
pull-out, anchorage loss, or the slip concentration at a crack — the exact
phenomena RC detailing studies care about.

The embedded-reinforcement ADR ([[20_ladruno_embedded_reinforcement_adr]])
introduces **Mode P** of `LadrunoEmbeddedRebar`: a penalty-coupled, own-DOF rebar
whose *axial* slot is a 1D bond-slip law. `LadrunoBondSlip` is that law. The
design constraints that shaped it:

- **A standard, parameterised backbone.** The fib Model Code 2010 / CEB-FIP τ–s
  curve is the field standard (codified, calibrated, four explicit regimes). The
  six shape parameters $\{\tau_{\max}, s_1, s_2, s_3, \tau_f, \alpha\}$ are
  **exposed**, not hard-coded, so users can dial in confined vs unconfined,
  good vs poor bond conditions per MC2010 Table 6.1-1.
- **Numerical robustness from the first iteration.** The raw power-law backbone
  has an *infinite* initial tangent (§4.1); left unregularized it kills the
  first Newton step. v1 bakes in the fix.
- **Mesh objectivity.** Bond softening localizes; a node-lumped bond law is *not*
  mesh-objective without fracture-energy regularization (§4.3). v1 exposes a
  bond fracture energy and ties the softening branch to it.
- **House hygiene.** No `exit()`; real `revert*`, `getCopy`, `sendSelf`/`recvSelf`
  round-trip, `Print`, and a `dissipatedEnergy`/`energy` response channel.

## 3. The bond-slip backbone (MC2010)

The constitutive law maps slip $s$ to bond stress $\tau$. It is **sign-symmetric**,
$\tau(-s)=-\tau(s)$, so the backbone is defined on the magnitude $a=|s|\ge 0$ and
the sign of $s$ is carried through. The four MC2010 regimes — *ascending power
law → plateau → linear softening → residual friction* — plus the v1 **initial
linear segment** (the D4.1 regularization, §4.1) give five pieces:

$$
\tau(a) =
\begin{cases}
k_0\,a, & 0 \le a < s_0 \quad\text{(initial linear segment, D4.1)}\\[4pt]
\tau_{\max}\left(\dfrac{a}{s_1}\right)^{\!\alpha}, & s_0 \le a < s_1 \quad\text{(ascending power law)}\\[10pt]
\tau_{\max}, & s_1 \le a < s_2 \quad\text{(plateau)}\\[4pt]
\tau_{\max} - (\tau_{\max}-\tau_f)\dfrac{a-s_2}{s_3-s_2}, & s_2 \le a < s_3 \quad\text{(linear softening)}\\[10pt]
\tau_f, & a \ge s_3 \quad\text{(residual friction)}
\end{cases}
$$

with the corresponding tangent

$$
\frac{\partial\tau}{\partial a} =
\begin{cases}
k_0, & 0\le a<s_0\\[2pt]
\dfrac{\alpha\,\tau_{\max}}{s_1}\left(\dfrac{a}{s_1}\right)^{\alpha-1}, & s_0\le a<s_1\\[8pt]
0, & s_1\le a<s_2\\[2pt]
-\dfrac{\tau_{\max}-\tau_f}{s_3-s_2}, & s_2\le a<s_3 \quad(<0,\ \text{softening})\\[8pt]
0, & a\ge s_3 .
\end{cases}
$$

**Parameters (all stresses; see §6 for the unit contract):**

| Symbol | Meaning |
|---|---|
| $\tau_{\max}$ | Peak bond stress (top of the ascending branch / plateau). |
| $s_1$ | Slip at which the plateau begins (peak slip). |
| $s_2$ | Slip at which softening begins (plateau end). |
| $s_3$ | Slip at which the residual plateau begins (softening end). |
| $\tau_f$ | Residual (frictional) bond stress, $a\ge s_3$. |
| $\alpha$ | Power-law exponent of the ascending branch (MC2010 typically $0.3$–$0.4$). |
| $G_F$ | *Optional* bond fracture energy; if given, **overrides $s_3$** (§4.3). |
| $s_0$ | *Optional* initial-linear-segment end; defaults to $0.1\,s_1$ (§4.1). |

The constructor sanitises the inputs (`setDerived`): $s_1$ is floored away from
zero, $s_0$ defaults to $0.1\,s_1$ and is clamped to $\le s_1$, and the regime
slips are kept ordered ($s_2\!\ge\!s_1$, $s_3\!\ge\!s_2$).

## 4. The three must-fix subtleties

These are the ADR-D4 "must-fix" items, all handled in v1.

### 4.1 Initial-slip regularization (singular first tangent)

The power-law branch $\tau=\tau_{\max}(a/s_1)^\alpha$ has, for $\alpha<1$,
$\partial\tau/\partial a \to \infty$ as $a\to 0^+$. The very first Newton
iteration (slip $\approx 0$) would then see an unbounded tangent and fail to
converge.

**Fix:** replace the power law on $[0, s_0)$ with a **linear segment** of finite
stiffness, secant to the power-law value at $s_0$:

$$
k_0 = \frac{\tau_{\max}\,(s_0/s_1)^\alpha}{s_0}, \qquad s_0 = 0.1\,s_1\ \text{(default)} .
$$

So `getInitialTangent()` returns the finite $k_0$, and a first tiny step gives
$\tau=k_0\,s$ — no blow-up. $s_0$ is tunable via `-s0` (and is used as the
unload/reload stiffness everywhere, §5).

### 4.2 Solution control on the softening branch (negative tangent)

Past $\tau_{\max}$, the softening branch $s_2\le a<s_3$ has a **negative
tangent**. Under **load control** a negative-stiffness branch diverges (the
structure snaps through).

> [!warning] Use displacement-like control past the peak
> Once any bond interface enters softening, drive the analysis with
> **`DisplacementControl`**, **`ArcLength`**, or an **IMPL-EX** scheme — *not*
> load control. This is the same regularization caveat that applies to every
> softening material in the fork (ASDConcrete3D, LadrunoRebarBuckling). v1 also
> nudges the *plateau and residual* tangents off a hard zero by a vanishing
> $\pm\varepsilon\,k_0$ to keep the assembled tangent non-singular at the seams.

### 4.3 Mesh objectivity (fracture-energy regularization)

A node-lumped bond-softening law is **not** mesh-objective: refine the rebar
discretisation and the dissipated bond energy changes, because softening
localizes into a single segment. The cure mirrors ASDConcrete3D's crack-band
regularization (see [[11_brick_asdconcrete_integration]]): fix the **energy
dissipated per unit interface area** rather than the slip at residual.

The bond fracture energy is the area of the softening triangle,

$$
G_F = \tfrac{1}{2}(\tau_{\max}-\tau_f)\,(s_3 - s_2)
\quad\Longrightarrow\quad
s_3 = s_2 + \frac{2\,G_F}{\tau_{\max}-\tau_f}.
$$

> [!note] Shipped semantics: `-Gf` overrides `s_3`
> In the shipped code, passing **`-Gf` $G_F$** (with $G_F>0$ and
> $\tau_{\max}>\tau_f$) **recomputes $s_3$** from the formula above, so the
> softening branch dissipates exactly $G_F$ per unit interface area regardless of
> the $s_3$ you typed. Omit `-Gf` (or pass $0$) to use the $s_3$ you supplied
> verbatim. The ADR's "7th explicit input / independent post-peak slope" sketch
> was simplified to this **`s_3`-override** so the constructor contract stays
> unambiguous — the code is the truth here.

For full mesh objectivity the *element* scales $G_F$ by the characteristic
length ($\text{lch}_\text{ref}/\text{lch}$); v1 is first-order $O(h^2)$ and wants
**≳ 6–8 segments per development length**. A distributed-Gauss bond integration
(rather than node-lumped) is the v2 upgrade.

## 5. State update, unload/reload, and tangent

The state update is a **damage-memory backbone with elastic unload/reload**:

- **Virgin loading** ($|s| \ge$ the committed peak $|s|$): the stress rides the
  envelope, $\tau = \operatorname{sgn}(s)\,\tau(|s|)$, with the branch tangent
  from §3. The peak-slip magnitude `CslipMaxAbs` is the **damage memory**.
- **Unload / reload** ($|s| <$ committed peak): the response is **elastic with
  stiffness $k_0$**, on a line through the committed-peak envelope point — i.e.
  the bond unloads stiffly toward the origin and reloads back up to where it left
  the envelope, then resumes softening. The tangent on this branch is $k_0$
  (finite, positive — no negative-stiffness surprise on unload).

This is a **monotonic-backbone** model: a purely monotonic history walks the full
backbone exactly, and reversals unload/reload elastically. **Cyclic
degradation** (bond deterioration, frictional pinching under load reversals) is
**deferred to v2**; v1's reversal behaviour is the simplest physically defensible
choice, not a calibrated cyclic law.

## 6. The unit contract — stress to force

`LadrunoBondSlip` returns a **bond stress** $\tau$ (force per unit *bar surface
area*). The embedded-rebar element converts that to a nodal **force** by
multiplying by the bar perimeter and the tributary bonded length:

$$
F \;=\; \tau \cdot (\pi d)\cdot L_\text{trib}
$$

where $d$ is the bar diameter and $L_\text{trib}$ the length of bar tributary to
the slip DOF. **Consequence:** $\tau_{\max}, \tau_f$ are **stresses** and the
*element* supplies the geometry $(\pi d\,L_\text{trib})$ — the material itself
knows nothing about bar size or spacing.

> [!note] Shipped wiring: one `-bondScale` scalar
> The shipped `LadrunoEmbeddedRebar` rolls the whole geometric factor into a
> single scalar **`-bondScale`** $bs = \pi d\,L_\text{trib}$ (default $1.0$), so
> $F = bs\cdot\tau$. The ADR's separate `-perimeter`/`-ltrib` flags were folded
> into this one number; supply the product. With `-bondScale 1` the element
> consumes $\tau$ directly as a force (handy for unit/validation tests where the
> τ–s law is exercised straight).

## 7. Command syntax — material

```tcl
uniaxialMaterial LadrunoBondSlip $tag $tauMax $s1 $s2 $s3 $tauf $alpha \
    <-Gf $Gf> <-s0 $s0>
```

Positional, in order: `$tag $tauMax $s1 $s2 $s3 $tauf $alpha` (seven required).
Then optional flags in any order:

- **`-Gf $Gf`** — bond fracture energy. If $>0$ (and $\tau_{\max}>\tau_f$),
  **overrides $s_3$** so the softening branch dissipates $G_F$ per unit area
  (§4.3). Omit to keep the typed $s_3$.
- **`-s0 $s0`** — end of the initial linear segment (§4.1). Defaults to
  $0.1\,s_1$; clamped to $\le s_1$.

```python
# OpenSeesPy — a representative MC2010 bond law (consistent units)
ops.uniaxialMaterial("LadrunoBondSlip", 1,
                     10.0,      # tau_max
                     1.0e-3,    # s1
                     3.0e-3,    # s2
                     10.0e-3,   # s3
                     2.0,       # tau_f
                     0.4,       # alpha
                     "-Gf", 0.0, "-s0", 0.0)   # both default if 0
```

```tcl
# Tcl — let Gf set the post-peak slope (overrides the s3 above):
uniaxialMaterial LadrunoBondSlip 1  10.0 1.0e-3 3.0e-3 10.0e-3 2.0 0.4 -Gf 0.024
```

## 8. Integration with the embedded-rebar element

`LadrunoBondSlip` is consumed by [[20_ladruno_embedded_reinforcement_adr|LadrunoEmbeddedRebar]]
(classTag **33005**), the penalty-coupled own-DOF rebar element. The bond-slip
law drives the bar **axial** slot; transverse coupling to the host is a separate
penalty tie (`-kt`). The shipped grammar (ADR §9 — it evolved past the original
`-mode P` sketch):

```tcl
element LadrunoEmbeddedRebar $tag $rebarNode {$nHost $h1..$hN | -host $eleTag} \
    {-shape $N1..$NN | -xi $x1..$x_ndm}  -dir $dx $dy [$dz] \
    ( -bond $matTag [-bondScale $bs] | -perfect $kAxial )  [-kt {$kt | auto}]
```

- **`-bond $matTag`** — the `LadrunoBondSlip` tag; this selects the **bond-slip**
  (DIANA-style, slip-admitting) path.
- **`-bondScale $bs`** — the $\pi d\,L_\text{trib}$ stress→force factor (§6).
- **`-perfect $kAxial`** — the *alternative* to `-bond`: a stiff penalty axial
  tie (**perfect bond, no slip**). This is the contrast case — choose one.
- **`-host $eleTag` / `-xi`** — let the host element supply its nodes and shape
  weights (the rebar stays host-agnostic); `-shape` works for any host.

The rebar element can be a truss (default, axial-only, explicit-clean) or a beam
(opt-in, adds transverse stiffness); **bond-slip acts axially only** in either
case. See the ADR §D5–D6 for the element-type trade-off and the Mode-T/Mode-P
history (Mode T is deferred indefinitely; all forward work is the Mode P
`-bond`/`-perfect` family).

## 9. Recording and energy channels

**Recordable responses** (material response on the embedded element, or the
direct uniaxial harness):

| Response | Returns |
|---|---|
| `stress` / `bondStress` | current bond stress $\tau$ |
| `strain` / `slip` | current slip $s$ |
| `tangent` | current $\partial\tau/\partial s$ |
| `energy` / `dissipatedEnergy` | cumulative bond work per unit area, $W=\int\tau\,\mathrm ds$ |

The **energy channel** (shipped via PR #177) accumulates the bond work
$W=\int\tau\,\mathrm ds$ by the trapezoidal rule, committed each step. This is
the **physical** bond energy the embedded element nets against the artificial
penalty energy in the energy-balance audit (ADR §10.2b) — a path-dependent
quantity, so it is committed in `commitState` and restored on `revertToLastCommit`.

The material also implements the full state cycle cleanly: `getCopy`,
`commitState` / `revertToLastCommit` / `revertToStart`, a `sendSelf`/`recvSelf`
round-trip (13 doubles: tag + 8 params + 4 committed state, including the work),
and `Print`.

## 10. Verification and validation

**Status:** shipped & CI-green on `ladruno` — `LadrunoBondSlip` (MAT 33002,
PR #168); energy channels (PR #177). classTags 33002 (material) and 33005
(element) are collision-free.

**Material battery** (`tests/test_ladrunoBondSlip_material.py`, Zone-A) verifies
the C++ backbone against a pure-Python MC2010 oracle through the direct uniaxial
harness (`testUniaxialMaterial` / `setStrain` / `getStress` / `getTangent`):

- **Backbone vs oracle** — a monotonic walk hits all five regimes (linear
  initial / power law / plateau / softening / residual) and matches the oracle to
  $10^{-6}$.
- **Plateau = $\tau_{\max}$**, **residual = $\tau_f$** — the named values at the
  obvious slips.
- **D4.1 initial tangent** — a first tiny step inside $[0, s_0)$ returns the
  **finite** $k_0$ (not a blow-up) and the right stress $k_0 s$.
- **Sign symmetry** — $\tau(-s)=-\tau(s)$ across representative slips.
- **D4.3 `-Gf` overrides $s_3$** — with $G_F>0$ the softening reaches residual
  *exactly* at $s_3=s_2+2G_F/(\tau_{\max}-\tau_f)$ (residual just above, still
  softening just below).
- **Elastic unload at $k_0$** — load into the plateau, unload one step: stress
  drops by $k_0\,\Delta s$ and the tangent is $k_0$.

The pull-out / development-length / lap-splice *element-level* validation legs
(DisplacementControl pull-out vs the CEB-FIP backbone, the rebar-refinement
mesh-objectivity leg, and the well-confined splice) live in the embedded-rebar
element's validation plan (ADR §6); the small-cover splice is **xfail'd** by
design (splitting is out of scope, §11).

## 11. Limitations and boundaries

- **Monotonic only (v1).** Reversals unload/reload elastically at $k_0$; there is
  **no cyclic bond degradation / pinching**. Cyclic deterioration → v2.
- **Well-confined bond only.** v1 targets pull-out where the surrounding concrete
  confines the bar. **Cover-controlled splitting** bond failure is **out of
  scope** (→ v1.1); a small-cover splice test is deliberately **xfail'd** so the
  model does not *silently over-predict* bond capacity for thin-cover details.
- **Softening ⇒ localization ⇒ control.** Past $\tau_{\max}$ the tangent is
  negative; use DisplacementControl / ArcLength / IMPL-EX, not load control
  (§4.2).
- **Mesh objectivity is first-order.** Node-lumped bond is $O(h^2)$ even with
  fracture-energy regularization; want ≳ 6–8 segments per development length.
  Distributed-Gauss bond → v2.
- **It is bond stress, not bar stress.** Pair it with a steel material on the
  rebar element; do not confuse the two slots (§1).
- **Small-strain bond assumption.** The element's tributary-length stretch
  ($L_\text{trib}\,\lambda$) under large host strain is deferred; v1 documents a
  host-axial-strain ≲ 2 % assumption (ADR §D6).

## 12. References

- **fib Model Code 2010** (fédération internationale du béton), §6.1 *Bond of
  embedded steel reinforcement* — the τ–s backbone and the parameter tables.
- **CEB-FIP Model Code 1990** — the original four-regime bond-slip law this
  generalises.
- **DIANA** bond-slip reinforcement (own-DOF rebar + 1D τ–s) — the modeling
  paradigm Mode P + this material implement.
- ASDConcrete3D crack-band / fracture-energy regularization — the mesh-objectivity
  pattern mirrored in §4.3 (see [[11_brick_asdconcrete_integration]]).
- ADR: [[20_ladruno_embedded_reinforcement_adr]] §D4 (this material), §D5–D6 and
  §9 (the consuming element), §10.2b (the energy channel).
