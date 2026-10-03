---
wp: LEGACY
title: "A token that answers to the wrong quantity can permanently SHADOW a later branch in the same strcmp ladder"
legacy_seq: 247
---
### A token that answers to the wrong quantity can permanently SHADOW a later branch in the same `strcmp` ladder
- **Bites:** `LadrunoEmbeddedNode::setResponse` listed `localForce` twice — once as a second spelling of the global-component tie traction (branch 1) and once, 16 branches later, for the genuine D9 local/interface-frame force (branch 17). First match wins, so branch 17's `localForce` was **dead code from the day it was written**, and anyone recording `localForce` silently got global components. `LadrunoIMKBeam`/`2d` had the same lie without the shadow: `force`, `globalForce` and `localForce` all returned `getResistingForce()`.
- **Tell:** a spelling appearing in two branches of one ladder is always a bug — either a copy/paste or a name that means two things. `grep -c '"localForce"'` per file is a cheap audit.
- **Workaround/status:** ✅ FIXED — `localForce` now resolves to the local-frame branch on `LadrunoEmbeddedNode`, is withdrawn from `LadrunoEmbeddedRebar` (no local frame exists there), and is a REAL element-frame force on both IMK beams (basic-`q` mapping, as DispBeamColumn). **Behaviour change for anyone who recorded `localForce` on those elements** — noted in [[LEDGER_implementations]]. *2026-07-28 (recorder-token consistency sweep).*
