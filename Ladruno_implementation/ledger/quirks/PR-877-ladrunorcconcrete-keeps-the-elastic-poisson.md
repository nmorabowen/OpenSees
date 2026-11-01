---
wp: PR-877
title: "LadrunoRCConcrete keeps the ELASTIC Poisson ratio after cracking — the transverse crack strain unloads the strut and costs 8–10 % of a panel's shear strength (…"
pr: "#877"
date: 2026-10-10
---
### LadrunoRCConcrete keeps the ELASTIC Poisson ratio after cracking — the transverse crack strain unloads the strut and costs 8–10 % of a panel's shear strength (C2 adds `-crackedNu`)
- **Bites:** Vecchio–Collins PV panels with the full MCFT flag set (`-beta -lublinerReduced -tensStiff vc -tensStiffC 200`) and ν = 0.2 peak 7–9 % below the MCFT (PV20 4.25, PV19 3.71, PV27 5.93 MPa); the PV20 match with the test (4.26) is a cancellation of this against the MCFT's own +7 %. The pre-peak curve is soft from the first crack on (γ = 2.3e-3: τ = 2.41 vs 3.03 MPa MCFT).
- **Why:** the spine is an isotropic effective-stress model, `σ̃ = σ̃_n + C0(E,ν):Δε`; the strut's effective stress carries `ν ε1` of the transverse tensile strain (ε1 ≈ 2–4e-3), so opening cracks unload the compression strut. MCFT/DSFM take ν = 0 once cracked.
- **Fix/status (2026-09-26, Ladruno C2):** `-crackedNu 0` (irreversible latch at `ε1 ≥ ε_cr`) reproduces the global-ν = 0 run to 0.1 % while keeping ν = 0.2 before cracking. Off by default. See [[LadrunoRCConcrete_guide]] §5.1.
