---
wp: ADR-86b
title: "TIMs F4 — psi / yieldDistance responses on LadrunoSANISAND"
status: "shipped — wp/86c"
section: "table"
legacy_seq: 159
---
| **TIMs F4 — `psi` / `yieldDistance` responses on `LadrunoSANISAND`** (request `_tims_proposed_model_requests_2026-09-07.md` F4; [[86_ladruno_sanisand_adr]]) — two read-only scalar material responses in the fork's `setResponse` override: `psi` (alias `stateParameter`) = the model's own `GetPSI(e, p')` with `p' = p + p_residual` floored at `small`, i.e. the psi that fed `M^b`/`M^d` at the last commit; `yieldDistance` (alias `yieldFunction`) = `GetF(sigma_n, alpha_n)`, the signed distance to the cone (negative inside, ~`mTolF` on it). Both read the COMMITTED state. Touches no upstream file; the vanilla `ManzariDafalias` still answers neither (asserted). | material response (diagnostic) | **none new** — response ids 33094 / 33095 in the ADR-86b 3308x band | SRC: `SRC/material/nD/LadrunoSANISAND.cpp` (setResponse/getResponse only). tests: `tests/test_ladruno_sanisand_responses.py` (psi and f reconstructed from `stress`/`state[24]`/`alpha` to 1e-12 relative on the confine-first deck, both `p_r` legs; on-surface drift gate on the plastic leg; vanilla-absence). docs: `LadrunoSANISAND_implex_guide.md` §6.1; banner line. | **shipped — wp/86c** | PR pending |
