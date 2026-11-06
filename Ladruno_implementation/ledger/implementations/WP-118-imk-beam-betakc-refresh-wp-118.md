---
wp: WP-118
title: "IMK beam betaKc refresh (WP-118)"
pr: "#853"
status: "fixed on the branch"
section: "table"
legacy_seq: 20
---
| **IMK beam `betaKc` refresh (WP-118)** — `LadrunoIMKBeam(2d)::commitState()` did not chain to `Element::commitState()`, so `Kc` (the committed stiffness behind `betaKc` Rayleigh damping) was never refreshed: captured correctly when `rayleigh` ran, then frozen at the initial stiffness — `betaKc` behaved like `betaK0`. Fix: chain to the base first, as `ElasticBeam2d/3d` do. Audit: no other fork element that uses Rayleigh skips the base call. **Behaviour change** for IMK + `betaKc` once the tangent changes. Gate `tests/test_imk_betakc.py`: betaKc-only differential vs `elasticBeamColumn`; the Corotational legs depart by 45% of peak pre-fix. | bug fix | — (existing IMK tags) | `SRC/element/ladrunoIMKBeam/LadrunoIMKBeam{,2d}.cpp`, `tests/test_imk_betakc.py` | **fixed on the branch** | #853 |
