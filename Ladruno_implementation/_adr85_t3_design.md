---
title: ADR-85 T3 — 2D mortar + ALM + tie + mortar friction + thickness, phase design notes
project: Ladruno
status: in progress (2026-08-18)
owner: nmora
tags: [implementation, contact, 2d, mortar, design-note]
---

# ADR-85 T3 — 2D mortar: phase design notes

Companion to [[85_ladruno_contact_2d_adr]] (§How/3, §How/7, gate G-T3, deferrals
D4/D5/D6). Records the measured baselines the gates diff against, the oracle
findings, the carried-debt decisions, and the failing-build stamps if the
pre-committed escalation fires.

## 0. Baseline at the phase tip (measured, not assumed)

| Item | Value |
|---|---|
| Branch | `claude/adr85-phase-t3-68a12a` off `origin/ladruno` (fresh worktree, zero extra commits) |
| Tip (T2 merge, PR #752) | `f33486853504c507b4caf194da2f0873e5b8b7be` |
| `ops.ladrunoBuild()` | `f33486853504c507b4caf194da2f0873e5b8b7be` — **== HEAD** |
| Artifacts | `dist/bin/OpenSees.exe` + `dist/bin/opensees.pyd`, 2026-08-18 13:31 (fresh full build, MUMPS included) |
| `contact_dump` baseline | SHA256 `B0F8F770267CB6AA42EC23C2342CCB8C8C936370F2FC445DD8E8ECA1655881E4` |
| ... reproducibility | 184 lines, 0 decks failed, **two runs byte-identical**; artifact preserved in session scratchpad as `contact_dump_baseline_f334868.txt` |
| 3D battery (pinned glob `test_adr39_contact*` + `test_adr41_*` + `test_adr57_*`, 28 files) | **142 passed** (matches T2's N) |
| 2D suite `tests/test_adr85_contact2d_*.py` (11 files) | **54 passed** (matches T2's N) |

**Cross-session dump-hash finding (recorded, not silently adopted).** T2's design
note records the dump at this same source tree as
`da73b6f8…7782`; this session's fresh build of the identical tree (harness
unchanged since T0, `git log -- tests/_testbed/contact_dump.py` stops at
`87a882464`) produces `B0F8F770…81E4`, twice, byte-identical. T2's worktree and
binaries no longer exist, so the two artifacts cannot be diffed; the causes
cannot be narrowed below {independent MUMPS/conan rebuild, toolchain drift
between builds}. **Consequence: `contact_dump` byte-identity is a
within-toolchain-session observable only.** The T3 gates diff against THIS
session's baseline, both `ladrunoBuild` stamps recorded — which is what the ADR
already prescribes. Quirks-ledger row owed.

## 1. Inherited rules that bind T3 (settled by measurement — do not re-litigate)

1. **Scalar tangential slip store** (T2 pre-experiment, `_adr85_t2_design.md` §1):
   a global-frame vector slip store loses the slip increment to catastrophic
   cancellation (median 1.5e-2, max 100% rel. error at |ds|/|s|~1e-15); a
   tangential scalar store is exactly zero-error. ADR-41 C3's 2D mortar friction
   therefore reuses `returnMap1D`/`tangentBlock1D` from
   `LadrunoFrictionKernel.h`, slip formed by projecting the relative
   DISPLACEMENT once onto t̂. No second scalar map (staleness trap).
2. **kn is already folded into `tangentBlock1D`'s dTN.** Reapplying kn at the
   scatter site gives a kn² tangent, invisible unless `-consistanttan` is
   exercised (T2 orchestrator-review catch, fixed at the call site there). The
   T3 mortar consistent tangent is gated explicitly on this.
3. **`-thickness h` applied ONCE where each density enters the residual** (ε_N,
   ε_T, μ_c, -tauMax/-cohesion clamps, tie stiffness); **λ_N is NOT separately
   scaled** (it inherits h from ε_N·ḡ in the Uzawa update; scaling again is an
   h² error — rev-2 panel catch).

## 2. Recon findings (measured against HEAD f33486853)

- **`LadrunoContact2DKernel.h` has NO mortar content** — pure NTS
  projection/vertex/B. T3 appends the 2D interval integrator to this header
  (ADR §Where names it there); pure append keeps the NTS section byte-stable.
- **3D mirror points** (`LadrunoMortarKernel.h`): `PairResult{status,ngp,area,
  gapL1,D,M,g}` (:99–108); status −1 refusal sites incl. `cos_t < 1e-12`
  (:420) — the exact 2D "perpendicular pair" analogue; sliver guard relative
  to `A_min` (:440–444); `A_phys = A_aux/cos_t` (:458–460); gN inside the
  integral (:484); `linearizePair` deliberately a stub — **the 2D mortar
  tangent stays the penalty Gram block only** (staleness-trap doctrine).
- **`returnMap1D`/`tangentBlock1D`** (`LadrunoFrictionKernel.h:240/:275`):
  scalar s/sp in the t̂ frame; `tFric` pre-negated; `spTrial` SET never `+=`;
  **`dTN` already carries kn** (comment :270–271; the T2 kn² catch at
  `LadrunoContactFE.cpp:589–596`).
- **The four `"T3"` refusal call sites**: mortar loop
  `LadrunoContactHandler.cpp:1870/:1895` (flip to live 2D branch), edge-edge
  `:2094/:2106` — edge-edge has NO 2D analogue (ADR scope fence), so these
  become a PERMANENT named scope-fence refusal, not a lane flip ("lands in
  T3" would become a lie).
- **`ndm(3)` hardcoded in the MORTAR ctor init list**
  (`LadrunoContactFE.cpp:154`) + ~14 stride-3 sites (resid :1273–1783, tang
  :1998/:2026/:2265). Template: the T1b/T2 `ndm2` SEGMENT ctor overload
  (`.cpp:221–267`) + separate 2D assembly functions (3D byte-identity by
  construction).
- **3D mortar friction is a VECTOR treatment** (`addMortarFriction`
  `.cpp:1431`: `gbarT[3]`, `st.gpT[3]`, `st.lambdaT[3]`, ALM offset
  `gTeff = (gbarT − gT0) + λ_T/kt`). In 2D this collapses to ONE scalar per
  slave node in slot 0 of the same `MortarNormalState` arrays (slots 1–2 ≡ 0)
  — the T2 cancellation finding makes the scalar store mandatory.
- **Uzawa site**: `LadrunoContactDomain.cpp:1176–1204`; tie branch :1185–1190
  (no clamp), normal `λ_N = min(0, λ_N + ε_N·ḡ)` :1201–1203. Untouched by T3.
- **`accumulateMortarTie` takes `const double rFacet[3]`** — 2D tie z-pads
  (the T1b `n[3]` idiom), no overload.
- **`WARN_2D_PERP_NO_NTS`**: zero stubs anywhere; append to `WarnTopic` after
  `WARN_VTX2D_SIDE_FLIP` (`LadrunoContactDomain.h:576`). No latch-listing
  query exists — the gate test asserts the stderr string via the child-process
  pattern (`_assert_refused`-style capture).
- **`-thickness`**: zero contact stubs. Parse in `OpenSeesOutputCommands.cpp`
  (⇒ new vanilla-ledger row), store on `MortarContact` (+1 param on
  `addMortarContact`, currently 33). h-multiply sites, ONCE each, at the 2D
  injection: `epsUse`, `epsTuse`, `muc`, `cohesion`, `tauMax`, `epsTie`;
  NEVER λ_N/λ_T/λ_Tie (they inherit h via `st.epsN`).
- **Ownership protocol**: the ADR's `:1179-1184` cite is STALE (that range is
  now the T1b vote magnitude gate); the real 3D machinery is
  `ladrunoEdgeEdgeOwns` (:555–577) + `continue` stand-down (:1951–1953). ADR
  correction owed.
- **Refusal-battery rows to invert**: `mortar_2d_not_yet_supported`
  (`test_adr85_contact2d_t0_refusals.py:137`, deck :349–363) + the same case
  referenced at `test_adr85_contact2d_t1b_nts.py:990`. There is NO tie row —
  T3 adds a tie acceptance case fresh. Precedent: the T1b
  `pair_2d_now_live` swap (asserts force TRANSFER, not "it runs").
- **D6 remainder is partially fictional**: RigidPlane/Mortar have NO per-node
  mass prescan at all (only NTS has one, `:864–901`); the only `getNodalMass`
  sites are `:460` and `:888`. The T2 note's "mis-attributing order" for
  those lanes referred to loops that do not exist.
- Housekeeping found: `stamp_headers.py` has a duplicate GLOBS entry for
  `LadrunoContact2DKernel.h` (:109, :117) — dedupe in T3.
- Banner line 52 tail to replace: `2D mortar/tie lands in T3 -- refused by
  name` → the T3 live clause; regen via `patch_banner.py` only.

## 2b. Binding design decisions (orchestrator)

1. **Kernel**: append `mortarSegmentPair2D(...)` (interval clip + D/M/g̃ on
   the slave-trace measure, 2-pt Gauss, `PairResult2D{status,ngp,len,gapL1,
   D[2][2],M[2][2],g[2]}`) to `LadrunoContact2DKernel.h`. Status codes mirror
   the 3D kernel: 0 accumulated, −1 empty/sliver/perpendicular (skip-with-
   routing), no −2 case (straight 2-node segments cannot be invalid the way
   quadratic facets can; zero-length is refused by τ_seg as in NTS).
2. **FE adapter**: new 2D MORTAR ctor overload (`ndm2` disambiguator, NSDMI
   members), separate `mortarActive2D`/`addMortarForce2D`/`addMortarTie2D`/
   `addMortarFriction2D`/2D tangent+visc branches; stride-2; 3D functions
   untouched byte-for-byte.
3. **Friction**: per-slave-node scalar on `MortarNormalState` slots [0]
   (`gpT[0]`, `gT0[0]`, `lambdaT[0]`), slip from projecting the weighted
   relative DISPLACEMENT once onto t̂; `returnMap1D`/`tangentBlock1D` reused;
   kn(=ε_N) appears exactly once in the consistent cross term.
4. **Routing**: handler computes `ntsCoDeclared` per mortar pair and passes
   it into the 2D MORTAR ctor. **Rule (corrected during T3 — the tag-equality
   form was unimplementable because NTS slave surfaces are node sets while
   mortar requires SLAVE_SEGMENTS, so the lanes can never share a slave
   tag):** an NTS Contact counts iff its MASTER surface tag equals the
   MortarContact's master tag AND its slave node set covers BOTH nodes of the
   mortar pair's slave segment (union over all matching NTS contacts,
   computed once per MortarContact). At evaluation, an ALIGN refusal
   (|cosθ| < τ_align) is silent iff ntsCoDeclared, else fires
   `WARN_2D_PERP_NO_NTS` once naming the pair. No double-count is structural:
   mortar contributes nothing on the refused pair; NTS enforces it. The
   master-only match was rejected as over-broad (NTS on unrelated slave nodes
   of the same master would silently swallow the latch — the ADR-78 P0
   shape).
5. **`-thickness h`** (default 1.0): parse-time h > 0 validation +
   requires `-mortar`; handle-time named FATAL if the pair is 3D and
   h ≠ 1.0. Applied once at the 2D injection site (list above).
6. **2D mortar + `-soft`**: SOFT=2 is fenced out of 2D by the ADR — named
   not-supported refusal (deferral row), never silent.
7. **`-ngp`**: 2D always uses 2-pt Gauss (exact for straight segments —
   oracle-proven); the flag is accepted and ignored in 2D with a doc note.
8. **Edge-edge on 2D**: permanent scope-fence refusal (new message text), both
   call sites.

## 3. Oracle findings (`proto_t3_mortar2d.py` — 70/70 PASS, verified by the
orchestrator: exit 0, seeded, repeat-run identical)

| CHECK | n | max rel err |
|---|---|---|
| 1 interval clip / refusals / τ_ov | 12 | 3.33e-16 |
| 2 2-pt Gauss exactness (vs exact rational quadrature, 3000 pairs) | 4 | 8.43e-15 |
| 3 tilted case + missing-Jacobian falsifier | 5 | 4.33e-16 |
| 4 mortar patch test (+ far-field twin at \|x\|~5e4) | 19 | 2.54e-16 |
| 5 tie / T-junction (ADR-41 C4) | 10 | 3.33e-16 |
| 6 weighted-gap sign convention | 7 | 1.74e-16 |
| 7 ALM / Uzawa (value = augTol residual by design) | 7 | 2.42e-06 |
| 8 `-thickness` scaling (bitwise) | 6 | 0.0 |

- The tilted twin without the alignment Jacobian fails by exactly cosθ
  (\|broken/correct − cosθ\| ≤ 2.2e-16 at cosθ ∈ {0.9, 0.7, 0.5}); at cosθ = 1
  correct and broken are bit-identical — the flat patch test is PROVABLY blind
  to the missing Jacobian. A 1-pt rule errs at 8.3e-2, so 2-pt "exact" is a
  result, not a tautology.
- **τ_ov default = 1e-9** (relative to min(L_s,L_m)/L_m — that normalizer is
  the largest attainable overlap, the 2D A_min analogue). Driver: far-field
  TRANSLATION noise re-represents [a,b] with ~9.0e-12 wobble; 1e-9 gives
  111× headroom. Companion **τ_align = 1e-6** (caps the Jacobian 1/\|cos\| at
  1e6). **Ordering τ_ov ≪ τ_align is a HARD, gated requirement**: for a
  contained slave b−a = L_s\|cos\|/L_m vs threshold τ_ov·L_s/L_m, so the
  sliver guard bites exactly when \|cos\| < τ_ov — reversed, it would swallow
  legitimate strongly-tilted overlaps instead of routing them.
- **A2 (amends §2b.1):** status codes stay {0, −1}, but `PairResult2D` adds a
  **`refusal` field (NONE/EMPTY/DEGEN/ALIGN) and `cos_t`**, set by the kernel
  that already computed them — otherwise the handler must re-derive cos_t to
  route `WARN_2D_PERP_NO_NTS` (the staleness trap on the exact quantity the
  gate turns on).
- **A5:** the tie patch is exact for fields with a TANGENTIAL gradient; a
  normal-gradient component differs by exactly g·n (physics). The tie gate
  deck therefore uses u = c₀ + c₁·x on a horizontal master.
- **A6:** the correct `-thickness 2` twin's committed multiplier is h·λ₁
  (λ inherits h); the bitwise-2× traction claim holds; double-scaling is the
  gated h² error.
- **A7:** ALM normalization adopted verbatim from shipped `proto_c2_alm.py`
  (ḡ = g̃/w, t = min(0, λ+ε_N·ḡ), compression negative).
- **A8 — C++ finding:** a noise-sized ACCEPTED interval is benign in the
  kernel (D/M/g̃ all scale with b−a) but the enforcement layer divides
  ḡ = g̃/w by the tributary — the active-set test needs its own tributary
  floor (the 3D `aGlobal <= 1e-300` guard is the precedent; use a relative
  floor).
- CHECK 8 reproduces the §How/7 panel catch: an UNSCALED cohesion/τ_max clamp
  flips the stick/slip decision under `-thickness 2`; h-scaled clamps leave
  it invariant.

## 4. Carried-forward debt decisions (D4/D5/D6) — decided, not rotting

- **D4 (`NTS2D_END_SLACK = 1e-3`) — RE-DEFERRED to T4, with reason.** The
  honest permanent fix (radial end-cap vertex pair at open terminals) edits
  the shipped, gate-protected 2D NTS pair path. Doing that inside the mortar
  phase couples two lanes' risk in one PR; T4 owns the benchmark battery
  (Hertz + full sweep) that would catch an end-cap regression, so the fix and
  its natural gate land together there. The window's cliff at 1+1e-3 remains
  disclosed in LEDGER_quirks.
- **D5 (`WARN_VTX2D_SIDE_FLIP` ungated) — FORMALLY RE-DEFERRED, demand-
  driven.** Four instrumented rigs failed to reproduce the cross-type flip
  (kinematic drive → singular; tethered → slave exits the wedge; untethered →
  coincidence floor / instability). The latch is wired defense-in-depth at a
  reachable site; what is missing is a REACHABLE deck, and constructing one
  has now failed 4× at bounded cost. Re-test trigger: the first collapse-deck
  application (the reachable case named by the T1b panel). Not a T3 blocker —
  T3 does not touch the vertex path.
- **D6 (RigidPlane/Mortar `-soft` prescan order) — CLOSED in T3, by
  evidence.** Recon measured that the "mis-attributing mass-check loops" for
  RigidPlane/Mortar DO NOT EXIST (`getNodalMass` appears only at :460 and the
  NTS scan :888); the T2 carry-forward described dead code. The serial
  mis-attribution for the mortar lane evaporates anyway once T3 makes the 2D
  mortar lane live (same mechanism as the NTS lane in T2). The MPI
  partitioning refusal loops (:809–837) are ADR-78 lane property (serial-only
  subsystem), unchanged. Quirks row records the correction.

## 4b. Review round (8-angle /code-review at high + orchestrator verification)

Refuter decision: **not required** — neither `LadrunoFrictionKernel.h` nor
`LadrunoMortarKernel.h` is touched (diff-stat proof); the 3D-path line edits
(preflight args, domain stream) are the exact class the `contact_dump`
observable certifies, and it is bit-identical (×2) with the battery at 142.
The orientation-vote extraction was verified **token-identical** by an
independent mechanical diff (finder B).

10 findings reported (5 correctness, 2 efficiency, 1 conventions, 2 docs) —
all fixed in-PR except two disclosed:

1. **FIXED — stale accumulator on the A8 floor-skip** (found independently by
   4 finders): the delta-keyed `accumulateMortarGap/Tie` now gets a ZERO
   report when a node falls below the tributary floor, mirroring the !active
   branch; a bare `continue` froze the previous eval's (g̃,a) into the Uzawa
   update and the penetration/tie queries. Reachable only via the T3 relative
   floor (3D's 1e-300 is unreachable).
2. **FIXED — `-epsN auto` × `-thickness h` h² double-count** (3 finders):
   auto-resolved penalties come from `getInitialStiff()` which already folds
   the element thickness (the §How/7 NTS absorption argument), so auto values
   are no longer h-multiplied; a defaulted ε_T inherits the provenance of the
   value it derives from (`epsTFromEpsN`).
3. **FIXED — tie ALIGN silent stand-down**: an NTS contact cannot own an
   equality bond; a TIE pair now always warns on ALIGN (tie-specific
   guidance) regardless of `ntsCoDeclared`.
4. **FIXED — vote slave-centroid double-count**: `ladruno2DOrientationVote`
   now dedups the slave list too (arithmetic NO-OP at the NTS call site —
   node sets have no duplicates — so NTS behavior is unchanged; the mortar
   SLAVE_SEGMENTS flat list duplicates interior vertices).
5. **FIXED — `LCD_FMT_VERSION` 1→2** (mortar record 37→38 slots): a v1
   stream now draws the named version refusal instead of a misalignment
   misdiagnosis.
6. **FIXED (efficiency)** — auto-penalty resolver memoized per master segment
   (its 2D inputs are invariant across slave segments, unlike 3D's
   per-facet orientDir); per-FE `mortar2dPerpWarned` latch spares the
   warnOnce set probe on permanently perpendicular pairs (latched on the
   ATTEMPT, win or lose the race).
7. **FIXED (docs/conventions)** — `// Ladruno` marker on the vanilla
   call-site edit; `mu` removed from the "h-scaled" ctor contract lists (it
   is dimensionless — an h-scaled μ is the §How/7 threshold error); SOFT=2
   fence comments corrected from "command surface" to the handle()-time FATAL
   (dimension unknown at parse); the preflight "no extra branch" comments now
   record the true reason (slave-loop-first ordering, not a per-node gate).
8. **DISCLOSED, not fixed** — the 3D `-thickness` backstop tests
   `hThickness != 1.0` exactly, so an explicit `-thickness 1.0` on a 3D deck
   passes silently; h=1 is force-neutral, and storing a `hasThickness` flag
   would grow the stream again for a diagnostic-only gain.
9. **DEFERRED (recorded debt)** — master-side routing match is declaration-
   tag identity while the slave side is node coverage (a per-contact-surface
   preprocessor pattern defeats it; failure mode is a LOUD spurious warn,
   not silence) — T4 hardening candidate. The friction tangent's unguarded
   `gT0[0]` read before first engagement mirrors the 3D twin's identical
   debt. Structural dedup of the 2D gather/Gram/Lref blocks is declined:
   the duplication mirrors the 3D twin's own structure (fork doctrine:
   twin-parity over DRY).

## 5. Failing STAMPED builds (pre-committed escalation)

Escalation rule: the **second** failing stamped build of the gate tests hands
the session to Opus-high. Stamps recorded here as they occur.

| # | Build stamp | What failed | Notes |
|---|---|---|---|
| — | **CLOSED GREEN** (artifacts 2026-08-18 15:09, review-fix build) | dump ×2 bit-identical `B0F8F770…81E4`; 3D battery **142 passed, N unchanged**; 2D suite **80 passed** (54 → 80; all G-T3 letters (a)–(f) + the post-review tie-ALIGN contract) | Escalation never fired: the single failing stamped build (#1 below) was root-caused test-side; the first C++-complete gate run was green. |
| 1 | artifacts 2026-08-18 14:19 (source = f33486853 + uncommitted T3) | 2D suite 12F/67P (all failures T3-new; 54 pre-T3 green; dump bit-identical, 3D 142 unchanged) | **ROOT-CAUSED: test-side, C++ exonerated by instrumented diagnosis.** (1) The mortar patch deck declared `-slave-segments 2` with a bare 4-node list `[200,201,202,203]` — the flat stride-2 pair convention makes that TWO disjoint segments with a hole where `[201,202]` should be; the parser legally accepts it (no shared node ⇒ the chain scan cannot object). The kernel-parity driver (MSVC, `t3_kernel_parity.cpp`) matched analytic D/M/g̃ exactly on all probe configs, and the "wrong" master forces [55.556, 88.889, 55.556] are the EXACT discrete solution of the holed surface (hand-solve: p = [0, 300] on the surviving left segment; M̄ᵀp gives 0.185185·300 = 55.556). (2) Tie/thickness decks used non-homogeneous `ops.sp`, unenforced under `constraints("LadrunoContact")` (warned in-run). Both defects systematic across the T3 files; fix delegated to the test author. **Quirks-ledger row owed: a holed 2-node segment listing is silently legal (indistinguishable from intentional disjoint patches) — the T4 user guide must state the chained-pairs convention loudly.** |
