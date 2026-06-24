---
title: Ladruno contact — session handoff (post C4 + D2)
project: Ladruno
status: handoff — mortar lane DONE (C1→C4) + viscous stabilization DONE (D2.1 NTS + D2.2 mortar)
owner: nmora
tags:
  - contact
  - handoff
  - roadmap
---

# Ladruno contact — session handoff

> Pick-up doc for the next session. The **authoritative** status/roadmap is the capstone
> [[48_ladruno_contact_capstone_adr]] (status-of-record table + unified roadmap); this is a
> point-in-time snapshot + recommended next step. Deferred set: [[47_ladruno_contact_deferrals_adr]].

## State (as of this handoff)

- **B3 SHIPPED** ([#389](https://github.com/nmorabowen/OpenSees/pull/389)) — the consistent ∂n/∂u
  geometric NORMAL tangent (`contact … -geomtan`) + the `ladrunoContactForce` query + a Hertz
  benchmark. Closes the geometric-tangent deferral the NTS lane carried since P2b (C2/C3/C4 all
  deferred it). `K_geom = kn·gN·∂²gN/∂u²` ⇒ quadratic Newton on curved/large-sliding; symmetric;
  gated off-default (indefinite → opt-in like `-consistanttan`). Battery **83→91/91**.
  NOTE: quantitative elliptic-Hertz `p(r)` is resolution-bound by NTS-penalty pressure recovery —
  4 new [[LEDGER_quirks]] entries (handler ignores imposed-disp; contact force not in nodeReaction;
  geom tangent indefinite; curved-indenter penalty bootstrap). Finer-mesh + displacement-control /
  D1 augmentation is the documented follow-up.
- **Contact battery: 91/91** — `tests/test_adr39_contact_p*.py` (+ `_p2b2c_hertz.py`) +
  `test_adr41_mortar_c2_{0,1,2}.py` + `_c3_{1,2,3}.py` + `_c4_{1,2}.py` + `test_adr41_viscous_d2.py`.

### Shipped this session
| PR | What |
|---|---|
| #381 / #382 / #383 | **C4 mortar mesh-tying** — permanent bond (zero-gap limit of contact): full 3-vec weighted relative-displacement `r_I → 0`, no-clamp ALM `λ_tie`. Resolved the shared-node MAJOR-1 quirk for the tie (linear `r` ⇒ the order-independent global accumulator `λ_N` uses). |
| #384 | Capstone status-of-record sync (C4). |
| #385 / #386 | **D2.1 viscous stabilization (NTS)** — `-visc μ_c`, RIGID_PLANE + SEGMENT. |
| #387 / #388 | **D2.2 mortar viscous** — extends `-visc` to the mortar contact ⇒ **D2 COMPLETE**. |

Workflow fix: recorded the `gh pr checks --watch` premature-exit trap in
[[../Ladruno_internal/WORKFLOW_GOTCHAS]] §3 (verify the actual Zone-A run by id, not the watch exit).

### Build / test recipe (this Windows box)
- Build: `Ladruno_scripts\build.bat OpenSeesPy` (exit 0).
- Run: `python3.12` with `PYTHONPATH=dist/bin LADRUNO_OPENSEES_QUIET=1 PYTHONIOENCODING=utf-8`.
- Oracles (numpy-only): `Ladruno_implementation/contact_prototypes/proto_*.py`.
- Zone-A CI is the merge gate; **verify the actual run id**, not `gh pr checks --watch` (it can exit
  before Zone-A registers): `RID=$(gh run list --branch <br> --limit 1 --json databaseId -q '.[0].databaseId'); gh run watch $RID; gh run view $RID --json conclusion`.
- Merge: `gh pr merge <n> --admin --squash`, base `ladruno`.

## Roadmap — what's left

### ✅ Done
- **Track A** (shared kernels): A1 friction extract #367, A2 projection extract #365.
- **Track C** (mortar lane): C1 kernel → C2 commit-cycle ALM → C3 friction → C4 mesh-tying. The whole
  implicit/accuracy lane.
- **Track D2** (viscous stabilization): NTS #385 + mortar #387.

### 📋 Remaining — COMMITTED (scoped, oracle-first ready)

**Track B — NTS explicit-stability hardening** (the explicit-first lane; independent of mortar):
- **B1** — P4 `SOFT=1` Courant-stable penalty (explicit `dt_cr` no longer throttled by `kₙ`).
  Gate: explicit stability + energy balance.
- **B2** — P5 `SOFT=2` segment-based penalty (corner / edge / T-intersection robustness).
- ~~**B3** — P2b-2c consistent `∂n/∂u` normal tangent + Hertz~~ **DONE [#389]** (see State above).
  Follow-up spun off: a quantitative FE Hertz harness (contact-force recorder ✅ via
  `ladrunoContactForce`; a robust curved-indentation driver — displacement-control or D1
  within-step augmentation — for the elliptic-`p(r)` / compliant-radius match).

**Track D1 — within-step augmentation refinement**: the held-load `analyze_augmented` proc is already
shipped + used by the C2/C4 tests, so D1 is largely delivered; what remains is a formal sign-off gate
(within-step `‖g̃‖→augTol` without recorder/load corruption) and folding the proc into a first-class
tested recipe.

### 🚫 Deferred — ADR-47 ([[47_ladruno_contact_deferrals_adr]]; each gated behind a re-open trigger that has NOT fired)
1. Dual / biorthogonal mortar basis (diagonal `D`, cheap nodal λ, LBB-optimal).
2. True-LM / saddle-point enforcement + inf-sup stabilization.
3. Self-contact (a surface vs itself).
4. Full slide-line Hermite smoothing (4) / averaged nodal-normal smoothing (4a) — faceted-normal chatter.
5. Anisotropic / elliptic friction (two principal μ).
6. Pressure/velocity-dependent `μ(N,v)` wired into contact (machinery exists in `frictionModel/`, unwired).
7. Custom global-Uzawa `LadrunoAugmentedNewton` `EquiSolnAlgo`.
8. NTN / NTS-via-mortar-weights unification.
9. Coupled-multiphysics surfaces (pore / thermal / acoustic) + pressure penetration.

### Disclosed limitations (not deferrals — won't change)
- Full ALM is **implicit-only** (explicit = single-pass penalty under the mass-only CDL).
- Both lanes are **penalty / Uzawa**, not exact-Lagrange stick.

## Recommended next step

**B1/B2** (explicit `SOFT` penalty tiers) — the remaining committed NTS hardening, for explicit
pounding/uplift robustness (`dt_cr` no longer throttled by `kₙ`). Or **D1** (within-step
augmentation sign-off — the `analyze_augmented` proc is shipped + used by C2/C4; what remains is a
formal `‖g̃‖→augTol` gate without recorder/load corruption, which ALSO unlocks the quantitative
Hertz follow-up via displacement-control-free curved indentation). After B3, the NTS lane's
deferred geometric tangent is CLOSED; everything else is ADR-47 (deferred behind triggers).
Process per fork standard: oracle-first → C++ → build → adversarial gate → PR on `ladruno`,
keep the capstone + ledgers current in the same PR.
