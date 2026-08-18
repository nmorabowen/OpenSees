---
title: ADR-85 T2 — 2D NTS friction, phase design notes
project: Ladruno
status: in progress (2026-08-18)
owner: nmora
tags: [implementation, contact, 2d, friction, design-note]
---

# ADR-85 T2 — 2D NTS friction: phase design notes

Companion to [[85_ladruno_contact_2d_adr]] (§How/4, §How/6, gate G-T2, deferrals
D3/D5/D6). This file records **the pre-T2 experiment and its decision**, the
measured baselines the gates diff against, and the failing-build stamps if the
pre-committed escalation fires.

## 0. Baseline at the phase tip (measured, not assumed)

| Item | Value |
|---|---|
| Branch | `claude/adr85-t2-2d-friction` off `origin/ladruno` |
| Tip (T1b merge, PR #751) | `bce163ccf8c24272c4f65f3dc21bb60f8799a9a6` |
| `ops.ladrunoBuild()` | `bce163ccf8c24272c4f65f3dc21bb60f8799a9a6` — **== HEAD** |
| Artifacts | `dist/bin/OpenSees.exe`, `dist/bin/opensees.pyd`, both 2026-08-18 11:32 |
| `contact_dump` baseline | `da73b6f89990c10fd5dc4cd3c1a3d5ce97d8ddb02bd50f0d9dd86ad58d5d7782` |
| ... reproducibility | 184 lines, 0 decks failed, **two runs byte-identical** |
| 2D suite `tests/test_adr85_contact2d_*.py` | **38 passed** |
| 3D battery `test_adr39_contact*` + `test_adr41_*` + `test_adr57_*` | **142 passed** (adr39 70 / adr41 50 / adr57 22) |

**Battery-count discrepancy (recorded, not silently adopted).** ADR-85's T0 log
and the T2 brief both cite the 3D battery as **121**. At this tip the
adr39/41/57 globs collect and pass **142**, and `git diff 87a882464..HEAD --
tests/` shows **no adr39/41/57 test was added or changed since T0** — so 142 was
also the count at T0 and the "121" refers to a narrower, unreconstructible
T0-era glob (no subset of 70/50/22 sums to 121). There is no
`tests/test_adr78_*.py` file at all, though the gate language names an "adr78
glob". **Resolution:** the gate that actually binds is *N unchanged across the
change at a pinned glob*, so T2 pins the glob above and its N = 142. The stale
121 / nonexistent adr78 glob should be corrected in ADR-85 (flagged to the lead).

## 1. THE PRE-T2 EXPERIMENT (§How/4) — reuse-vs-scalar for the friction return map

The panel established that the plane-constrained parity test has **zero
discriminating power** as a referee (on a plane-constrained configuration a
correct scalar path, a degenerate 3D path, and a *wrong-metric* 3D path all
pass), so the decision is made here, before any C++, and the parity test is
demoted to a regression gate.

### Method

Rather than transcribe the 3D kernel into numpy (which would only test the
transcription), the harness **includes the shipped header and calls it**:

- `contact_prototypes/t2_reuse_vs_scalar_check.cpp` — `#include
  "LadrunoFrictionKernel.h"`; evaluates **Path A** = the shipped
  `frictionReturnMap` / `frictionTangentBlock` on plane-constrained configs
  (`n = (cosθ, sinθ, 0)`, `t̂ = perp(n)`, slip vectors in-plane) against
  **Path B** = an independent scalar (`returnMap1D`) on the same inputs; dumps
  both at `%.17g`. 200 000 randomized configs (log-spaced `kt`, `kn`, `N`, `μ`,
  cohesion, `τmax`; 34 % of cases parked **astride the cone radius** perturbed by
  a few ulp, to probe the stick/slip threshold).
- `contact_prototypes/proto_t2_reuse_vs_scalar.py` — the referee. Partitions
  rows by branch agreement *first* (comparing a stick tangent against a slip
  tangent is meaningless), then adjudicates every threshold disagreement in
  **exact rational arithmetic** (`fractions.Fraction` — exact, stronger than the
  60-digit float used in T1a) to decide which path took the *correct* branch.
- `contact_prototypes/proto_t2_slip_cancellation.py` — the mechanism probe for
  the dominant finding below.

Reproduce:

```bash
g++ -O2 -std=c++17 -ISRC/domain/contact \
  Ladruno_implementation/contact_prototypes/t2_reuse_vs_scalar_check.cpp -o exp.exe
./exp.exe > dump.txt          # stderr carries the branch-disagreement count
python3.12 Ladruno_implementation/contact_prototypes/proto_t2_reuse_vs_scalar.py
python3.12 Ladruno_implementation/contact_prototypes/proto_t2_slip_cancellation.py
```

### Measured results

**(1) Stick/slip threshold correctness — the gate-relevant one.** 8 991 / 200 000
configs (4.5 %; 13.3 % of the threshold probes) take a **different stick/slip
branch** in the two paths. Adjudicated in exact rational arithmetic over all 167
disagreements in the dumped window:

| | wrong branch |
|---|---|
| Path A (degenerate 3D) | **149 (89.2 %)** |
| Path B (scalar) | **18 (10.8 %)** |

Path B's residual errors are the *single unavoidable* rounding of
`fl(kt·(s−sp))`; Path A adds avoidable `sqrt(tx²+ty²)`-vs-`fabs` and
normalization rounding on top. **G-T2(a) requires the threshold to be exact at
`tanφ = μ`** — Path A is measurably the wrong tool for that gate.

**(2) The dominant finding — catastrophic cancellation in a vector-stored slip.**
The 3D kernel takes the slip as a **global-frame 3-vector**. Storing it that way
means the tangential origin and the current slip are each projected onto the
global axes and rounded *there*, so their componentwise difference cancels:

| slip-increment relative error | `|ds|/|s|` 1e-16..1e-14 | 1e-14..1e-12 | 1e-12..1e-10 | 1e-10..1e-8 |
|---|---|---|---|---|
| Path A (global-frame vector store) | **1.5e-2** | 2.9e-4 | 3.2e-6 | 4.3e-8 |
| Path B (tangential scalar store) | **0** | 0 | 0 | 0 |

(medians; Path A max relative error **1.0**, i.e. 100 %.) Path B is *exactly*
zero-error by Sterbenz. In the regime a converging Newton iteration actually
lives in — small slip increments on a finite slip origin — **Path A loses the
slip increment's relative accuracy entirely**, and the slip increment is
precisely what the return map consumes. This is the ADR-41 C3
displacement-not-position crux resurfacing in a new frame: it is not enough to
measure slip from displacements, the slip must also be **kept as a scalar in the
tangential frame**, never as a global-frame vector differenced componentwise.

**(3) Slip tangent structure.** In a 1-D tangent space the slip-branch tangent
`(Pt − n̂⊗n̂)` is **identically zero**. Path A reaches exact zero in only 28.9 %
of slip cases (residual ≤ 4.5e-16·kt); Path B is structurally zero. Stick
tangent: Path A equals `kt` exactly in 54.2 % (≤ 2 ulp otherwise); Path B
exactly.

**(4) Leakage.** Out-of-plane leakage `|tF_z| + |gp_z|` is **exactly 0** — the
one thing reuse gets structurally right. But friction force leaks into the
**normal** direction, `|tF·n|` up to 3.85e-8 *relative* to `|tF|`: a degenerate
3D path silently perturbs the normal equilibrium the frictionless T1b lane
already gates.

### Decision

**Scalar `returnMap1D`, beside the 3D kernel in `LadrunoFrictionKernel.h`,
sharing the constants and the unified-cone clamps.**

Exactness and simplicity point the same way, so there is no trade-off to
arbitrate: the scalar path is exact where the gate demands exactness (1), is the
only one that preserves the slip increment (2), gives the structurally correct
slip tangent (3), and cannot leak into the normal direction (4) — while being
~20 lines against the 3-vector machinery. Reuse was the a-priori-cheaper option
and it loses on its own merits, measured.

**Binding implementation consequence:** the 2D path stores the tangential origin
and plastic slip as **scalars** in `FrictionState` (slot 0; slots 1–2 ≡ 0, per
§How/4), and forms the slip by projecting the **relative displacement** onto `t̂`
once. It must NOT store a 2-vector and difference it componentwise — that
reintroduces finding (2) exactly. `FrictionState` needs no new members: `gpT[0]`,
`gpTtrial[0]`, `gT0[0]`, `gT0committed[0]` carry the scalars and the existing
`commit`/`revertToLastCommit` double-buffer (contact-review P2) carries over.

## 2. Scope corrections found during recon

- **There is no `-mu` refusal row to invert.** The T2 brief and rev-2 phase plan
  describe "the retired `-mu` row" inverting into an acceptance test on the
  `pair_2d_now_live` precedent. The `-mu` 2D FATAL exists in C++
  (`LadrunoContactHandler.cpp:995-999`) but **no test anywhere in
  `tests/test_adr85_contact2d_*.py` asserts it** — the refusal-battery CASE rows
  are `mixed_dim`, `ndf3_frame_node_on_nts`, `nps2_on_3d_deck`,
  `pair_2d_now_live`, `mortar_2d_not_yet_supported`, `cross_dim_pair`,
  `plane_2d_before_surface`, `nps2_on_2d_deck_declares_ok`. So T2 **adds** a
  friction acceptance row (still following the `pair_2d_now_live` precedent and
  still asserting force *transfer*, not merely "it runs"); nothing is retired.
  Consequence: removing the FATAL cannot break a test, and nothing ever proved
  the FATAL fired — the new acceptance row is the first coverage either way.
- **`softKt`'s 2D branch is already coded and provably dead**
  (`LadrunoContactFE.cpp:773-783`, "Unreachable until T2 wires 2D friction").
  T2 **verifies it goes live**, per the brief — it is not rewritten.
- **D5's hazard is concrete.** The vertex adapter stores raw `Node*`
  `nts2dPrevFar`/`nts2dNextFar` (`LadrunoContactFE.cpp:255-256`, armed at
  `LadrunoContactHandler.cpp:1349-1401`) that are **outside the pair's own
  connectivity** (a vertex pair's nodes are slave + vtx). `pruneMissingNodes`
  (`LadrunoContactSurface.cpp:47`, driven from
  `LadrunoContactDomain.cpp:1324`) prunes surface segments; it does not
  revalidate those adapter-held pointers. The G-T2(f) vertex case is a
  use-after-free gate, not a bookkeeping gate.

## 3. Orchestrator review findings on the delegated implementation

**FIXED — `kn²` in the consistent friction cross term (`-consistanttan` only).**
`LadrunoContactFE::addFrictionTang2D` assembled
`Kss_s·(t̂⊗t̂) + dTN_s·kn·(t̂⊗n)`, but `tangentBlock1D` already returns
`dTN = −(∂cap/∂N)·kn·sgn` with `kn` **included** — exactly as the shipped 3D
`frictionTangentBlock` carries it once in `(−dCap_dN·kn·nh[i])·n[j]`. The extra
factor made the cross term scale as **kn²**; at a realistic `kn ~ 1e6..1e9` that
is a catastrophically wrong tangent, not a subtle one. The symmetric **default**
was unaffected (`dTN_s ≡ 0` unless `consistent == true`), so this would have
escaped every gate that does not exercise `-consistanttan`. Fixed at the call
site (`kn` appears exactly once) and both misleading doc comments corrected.

**CHECKED AND NOT A DEFECT — `cohesion`/`tauMax` on the 2D NTS lane.** The 2D
adapter calls `returnMap1D` without them (defaults 0, 0), and `-cohesion` /
`-tauMax` are parsed with no dimension gate
(`OpenSeesOutputCommands.cpp:747,754`) — which looks like silent wrongness. It is
not: the **shipped 3D NTS** call omits them identically
(`LadrunoContactFE.cpp:1187-1188`), because the unified cone's cohesion/τmax
clamps are an ADR-41 C3.1 **MORTAR** feature, not an NTS one. The 2D lane
therefore matches its 3D twin exactly, which is the parity contract. Whether the
NTS lane should *refuse* those flags is a pre-existing 3D question, not a T2
regression — not in scope, flagged only.

**VERIFIED — no shipped 3D path touched.** `LadrunoFrictionKernel.h` is a pure
append (zero deletions ⇒ the 3D kernel is byte-identical by construction). Every
deletion in `LadrunoContactFE.{h,cpp}` / `LadrunoContactHandler.cpp` sits inside a
2D branch, a 2D helper (`segment2DActive`), or the **2D-only** SEGMENT ctor
overload — the one taking `ndm2`, called only from the two 2D sites
(`LadrunoContactHandler.cpp:1382,1431`); the 3D SEGMENT ctor is a separate
overload used at `:1809`. This is *structural* grounds for the μ=0 byte-identity
expectation, independent of the gate that measures it.

### Byte-identity refuter on the shipped 3D friction path (pre-merge condition)

The phase brief required a refuter pass "if the implementation touches ANY shipped
3D friction path (`LadrunoFrictionKernel` or the 3D SEGMENT friction arms)". The 3D
*functions* were not modified (pure append), but the *file* was edited, so the
condition was discharged rather than argued away.

`contact_prototypes/t2_3d_byteidentity_refuter.cpp` compiles the **pre-PR (HEAD~1)**
and **post-PR** kernels into two separate namespaces and adversarially tries to make
`frictionReturnMap` / `frictionTangentBlock` / `frictionCap` disagree **bit-for-bit**
over 400 000 randomized configs, deliberately including every branch boundary
(`cap <= 0` free slip, the exact stick/slip threshold, the `tauMax` cap binding,
separated `N <= 0`, both `consistent` settings, `kn`/`kt` over 1e1..1e9).

**Result: 0 differing bits / branches — refutation FAILED, i.e. the 3D friction
kernel is bit-identical across this PR.** Corroborated independently by the
`contact_dump` harness (whose 7 canonical decks include a 3D NTS-friction deck)
hashing `da73b6f8...7782`, unchanged, and by the 3D battery at 142 passed.

Reproduce:

```bash
git show HEAD~1:SRC/domain/contact/LadrunoFrictionKernel.h > old.h   # + strip guards
g++ -O2 -std=c++17 -I. Ladruno_implementation/contact_prototypes/t2_3d_byteidentity_refuter.cpp -o refute.exe
```

## 4. Failing STAMPED builds (pre-committed escalation)

Escalation rule: the **second** failing stamped build of the parity/gate tests
hands the session to Opus-high. Stamps recorded here as they occur.

| # | Build stamp | What failed | Notes |
|---|---|---|---|
| — | — | none yet | — |
