---
wp: ADR-46
title: "wipe() does NOT recreate the Domain — new domain-level state MUST be reset in Domain::clearAll() (ADR 46 P1)"
legacy_seq: 158
---
## `wipe()` does NOT recreate the Domain — new domain-level state MUST be reset in `Domain::clearAll()` (ADR 46 P1)

`ops.wipe()` calls `Domain::clearAll()` on the SAME Domain object; nothing is
reconstructed. Any new domain-level member you add therefore leaks across model
generations unless you reset it in `clearAll()` yourself (the ADR-30 `theEQs` and
ADR-39 contact-engine cleanups are the same lesson). ADR 46 P1 hit this twice in
one PR: (1) the new Rayleigh-factor domain copy survived wipe → the NEXT model
reported phantom damping; (2) the UPSTREAM `theEigenvalues`/`theEigenvalueSetTime`
also survive wipe — latent forever because nothing read the spectrum across a wipe
until `complexEigen` did (a stale spectrum from the previous model silently answers
for the new one). Both now reset in `clearAll()` (`// Ladruno ADR46`).

Related trap: `Domain::getEigenvalues()` **exit(-1)s** (kernel-killer, MPCO class)
when the spectrum was never set — any graceful-failure caller must probe via the
additive `Domain::getNumEigenvalues()` (ADR 46) BEFORE calling it.

Related trap: `region ... -rayleigh` and per-element `-rayleigh` write factors
straight onto elements via `MeshRegion::setRayleighDampingFactors` / element
`setRayleighDampingFactors` and NEVER touch `Domain::setRayleighDampingFactors` —
the domain-level copy reflects only the last GLOBAL `rayleigh` call. `complexEigen`
scans elements and warns on mismatch (D1 policy: warn, never silently absorb).

ADR 46 P2 pre-warning (from the P1 Opus gate): once Route B projects a full
`ΦᵀCΦ` and emits `ψ=Φz`, REPEATED eigenvalues need the basis M-orthonormalized
WITHIN each repeated eigenspace (ARPACK does not guarantee it); the P1 closed form
is immune (uses eigenvalues only, never Φ). **RESOLVED BY DESIGN at P2:** the
assembler projects the FULL `M̃=ΦᵀMΦ` (never assumes `I`) and synthesizes
`K̃=sym(M̃·diag(ω²))` (exact: `KΦ=MΦΛ` column-wise, and `M̃` commutes with `Λ`
within a repeated eigenspace), so the QZ pencil is exact in ANY normalization and
for repeated eigenvalues — no re-orthonormalization step exists to get wrong.

Related trap (ADR 46 P2, OpenSees classic): `getDamp()`/`getMass()` on Elements
AND Nodes habitually return the SAME per-class/per-size static scratch matrix —
never hold references to both at once; deep-copy the first before calling the
second (`const Matrix Ce(el->getDamp());` then consume `el->getMass()` directly).

Route-B validation of the `-doRayleigh` contract (ADR 46 P2 test discovery):
**`Truss` also defaults `-doRayleigh` to 0** (the [[project_damping_channels]]
"OFF for the rest" family, not just zeroLength) — a default Truss feels NO global
`rayleigh` betaK in a transient, and the assembled complexEigen (correctly!)
reports zeta=0 for such a model while `-closedForm` claims bK·w/2. When a test or
model expects stiffness-Rayleigh on trusses, build them with `-doRayleigh 1`. The
divergence between the two routes on default elements is the P1-vs-P2 contract
difference, not a bug: Route B = the C the transient feels; closed form = the
global factors as if every element carried them.
is immune (uses eigenvalues only, never Φ).
