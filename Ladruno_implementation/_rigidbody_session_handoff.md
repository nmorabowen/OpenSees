---
title: "RigidBody (ADR-58) — session handoff"
project: Ladruno
type: handoff / next-session pointer
related:
  - "[[58_ladruno_rigid_body_adr]]"
updated: 2026-06-24
---

# RigidBody (ADR-58) — session handoff

Next-session pointer for `LadrunoRigidBody`. Full detail lives in the ADR and its
**Implementation log**: [[58_ladruno_rigid_body_adr]]. This is the short version.

## Status

| Piece | State | PR |
|---|---|---|
| **P1** — translation/ballistic (zero-stiffness Element + internal CoM node + rigid-link MP slaving + condensed mass) | ✅ merged on `ladruno` | [#426](https://github.com/nmorabowen/OpenSees/pull/426) |
| G9 manifest row for `ELE_TAG_LadrunoRigidBody` (33015) | ✅ merged | [#427](https://github.com/nmorabowen/OpenSees/pull/427) |
| **P2 Stage 1** — free-spin body-frame SO(3) rotation (`commitState` side channel) | ✅ merged | [#428](https://github.com/nmorabowen/OpenSees/pull/428) |
| **P2 Stage 2** — slave-following + moment gather + Housner | ⏭ **NEXT (not started)** | — |

The element today: a 6-DOF rigid body that **translates** correctly (ballistic gate)
and **rotates freely** (Dzhanibekov gate: `‖L‖` exact, energy 2nd-order, intermediate-
axis flips). It is **wip** — `status: planned` in `testbed/manifest.yaml`, **no banner
line** until the rocking MVP lands. Applied moments are currently **ignored** (free
spin only), and the **slaves do not yet track the body's rotation** (they follow only
translation — the frozen P1 transport block).

## Where to work (IMPORTANT — read [[ladruno-build-in-worktree-not-shared-checkout]])

- **Worktree:** `C:\Users\nmb\Documents\Github\OpenSees-rigidbody` on branch
  `feat/adr58-rigidbody-p2-stage2` (off the latest `ladruno`). MUMPS is junctioned in
  from the main checkout so the build skips the 15-min MUMPS rebuild. **Do all
  build+commit work here, NOT in the shared `…/OpenSees` checkout** (its branch
  switches under you — that already caused a commit-on-wrong-branch + a binary leak
  this session).
- **Build:** `& cmd.exe /c "cd /d <worktree> && Ladruno_scripts\build.bat OpenSeesPy"`
  via **PowerShell** (the Bash tool's `cmd /c … | tail` silently no-ops — only the
  banner prints, nothing recompiles). **Always verify `dist/bin/opensees.pyd` mtime**
  is newer than your edits before trusting a test run.
- **Run tests:** `python3.12` against the worktree build; the module imports as
  `import opensees` from `dist/bin` (add it via `os.add_dll_directory`), NOT
  `openseespy`. Both `rigidbody_tests/*.py` do this already.
- **Pre-PR:** run `python ci/check_classtags.py` + `ci/check_manifest.py` locally — the
  fork **merges even on a red gate** (no branch protection), so catch G9/manifest
  failures before the PR, not after.

## The element's current design (so you don't re-derive it)

- **Files:** `SRC/element/ladrunoRigidBody/{LadrunoRigidBody.{h,cpp},OPS_*.cpp,CMakeLists.txt}`.
  `ELE_TAG_LadrunoRigidBody = 33015`. Created via `element LadrunoRigidBody tag N
  s1..sN [-mass m] [-internalNode tag] [-omega wx wy wz]`.
- **Architecture (D1b):** a **zero-stiffness `Element`** (`K=M=P=C=0`) owning a private
  internal 6-DOF CoM `Node` (Joint3D pattern). Slaves tied by rigid-link
  `MP_Constraint`s. Mass condensed to the internal node; slaves zeroed to massless
  followers. (DomainComponent was rejected: >15-file cost + `Domain::commit` skips
  DomainComponents.)
- **Rotation (P2 S1):** integrated in `commitState`, **off the global solve** (D3 — a
  `DiagonalSOE` can't carry the dense inertia/gyroscopic coupling). State = **spatial
  angular momentum `L`** + orientation `Versor q`. Each step:
  `ω_body = Ibody⁻¹·(qᵀL)`; **2nd-order midpoint** `q ← q·exp(ω_mid·dt)`; `L ← L+m·dt`
  (`m=0` for now). Reuses `SRC/matrix/{GroupSO3.h,Versor.h}`; a local `rotateVec`
  replaces `Versor::rotate` (which needs left-scalar mult `VectorND<3>` lacks).
  `dt = getCurrentTime() − lastTime`, with a `started` flag (no giant first-`dt`).
- **Recorder (D9):** `eleResponse(tag, …)` → `orientation` (quaternion),
  `omega` (spatial), `omegaBody` (the flip metric), `angularMom` (spatial `L`),
  `comDisp`/`comVel`/`comAccel`, `mass`.
- **Tests:** `Ladruno_scripts/rigidbody_tests/test_p1_ballistic.py`,
  `test_p2_dzhanibekov.py` (both gate-passing). These are standalone scripts, **not**
  Zone-A pytests yet — the manifest `pytest` row is `PENDING` until a `tests/` port.

## P2 Stage 2 — the plan (start with a design pass)

The crux is **architectural**, so do a source-grounded design pass (Plan agent) BEFORE
coding, like P1/P2-S1 — see ADR §6 open questions and §7 P2 gate.

1. **Slave-following under finite rotation (the hard part).** The side-channel `q` must
   drive the slaves' positions, but `MP_Constraint` is **linear/homogeneous**
   (`u_c = Ccr·u_r`, no constant term) and the internal node's rotation DOFs are not
   the orientation carrier. Candidate mechanisms to evaluate:
   - a **time-varying `MP_RigidLink`** whose `getConstraint()` rebuilds the transport
     block from the current lever arm `R·d_i⁰` (the Transformation handler re-reads it
     each step when `isTimeVarying()==true` — verified at
     `TransformationDOF_Group.cpp:896-908`); vs.
   - the element **imposing slave trial displacements** directly each step; vs.
   - making the internal node's rotation DOFs the live orientation (gives up the clean
     side channel). Pick one with the design pass; the time-varying MP is the leading
     candidate and avoids Domain churn.
2. **Force/moment gather** — `gatherCoMTorque()` is stubbed to `0`. Implement
   `m_body = Σ (R·d_i⁰) × f_i` from slave loads/reactions (D3 dual of slaving) so an
   applied moment / toe-contact reaction actually drives the spin. Until this lands,
   any rotational load is silently ignored.
3. **Housner rocking-block gate** — the headline use case: a rigid block on a
   **compression-only toe** contact. The Plan pass recommended a `zeroLength` ENT/gap
   spring at the bottom-corner nodes (a *test-harness* contact, not part of the element;
   ADR D5 fences the contact engine out). Check the analytic Housner half-period +
   post-impact ω-ratio `(1−(3/2)sin²α)`.

**Deferred beyond Stage 2 (ADR):** parallel `sendSelf/recvSelf` (serial-only v1 — the
SO(3) state is self-contained but not serialized), implicit tangent (explicit-only v1),
and the AEM debris hook (D5, fenced out).

## Review posture

P1 and P2-S1 each had a focused adversarial review (findings folded in). P2-S2's
slave-following is the genuinely hard architectural call — **warrants the same
design-pass + review treatment**. The "run it, don't just review it" lesson held twice
(the P1 `setDomain` valid-flag bug and the P2 dissipative-scheme bug both surfaced only
on execution, exactly as the reviews predicted).
