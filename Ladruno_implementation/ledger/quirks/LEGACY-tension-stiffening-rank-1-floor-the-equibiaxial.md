---
wp: LEGACY
title: "Tension-stiffening rank-1 floor: the equibiaxial self-consistency coefficient is 0.5, not 1"
legacy_seq: 82
---
### Tension-stiffening rank-1 floor: the equibiaxial self-consistency coefficient is 0.5, not 1
- **Bites:** the RC-shell tension-stiffening floor (RC stack Phase 3a, `-tensStiff`) pins the principal tensile stress to `σ_ts(ε1)` by injecting `Δ·(p1⊗p1)` and measuring `n^Tσn`. For a true unit eigenvector the self-consistency coefficient `a²+b²+2(ab)² = (p1x²+p1y²)² = 1`, so one injection pins the normal exactly. In the DEGENERATE (equibiaxial) branch the natural fallback `(a,b,ab)=(0.5,0.5,0)` is **rank-2 (isotropic), not rank-1**, so its coefficient is `0.5²+0.5² = 0.5` — injecting `Δ·(0.5,0.5,0)` raises the mean normal by only `0.5·Δ`, so the floor reaches **halfway** to `σ_ts`, never `σ_ts`.
- **Fix (proven):** in the degenerate branch use distinct injection vs measurement vectors — inject `g=σ_ts−mean` to BOTH in-plane normals with full weight (`ts_inj=(1,1,0)`), measure the mean (`ts_meas=(0.5,0.5,0)`); then `ts_meas·ts_inj=1` and each in-plane normal reaches `σ_ts` exactly (verified to 1e-16 in the standalone g++ equibiaxial gate). General rule: the floor's self-consistency requires `ts_meas·ts_inj==1`, which holds for a rank-1 `p⊗p` but NOT for the rank-2 isotropic `½(I)`. A pure uniaxial gate (the original T1) never exercises the degenerate branch, so this only surfaces under an EQUIBIAXIAL test. Learned 2026-06-18 (caught by a 3-agent adversarial review of [[19_ladruno_rc_shell_adr|LadrunoRCConcrete]] Phase 3a).
