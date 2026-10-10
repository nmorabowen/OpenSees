---
wp: PR-715
title: "715 -- upstreamable-table row(s)"
pr: "#715"
files: ["`SRC/element/tetrahedron/TenNodeTetrahedron.cpp`", "`SRC/element/tetrahedron/FourNodeTetrahedron.cpp`"]
table: "upstreamable"
legacy_seq: [411, 412]
---
| `SRC/element/tetrahedron/TenNodeTetrahedron.cpp` | `// Ladruno`: **`setResponse` accepted ONLY the plural token spellings.** The two continuum branches were bare `strcmp(argv[0], "stresses")` / `"strains"`, so `-E stress` (singular) returned a null `Response` and the element was dropped from the recorder — silently, because the orchestrator had no diagnostic for a null response (fixed in the same PR). Neighbouring elements (`LadrunoBrick`, `FourNodeQuad`) take both spellings, so a mixed model recorded some element classes and not others from the SAME request. Both sites now route through `LadrunoResp::is(argv[0], "stress"/"strain")` (`SRC/element/LadrunoResponseTokens.h`, the canonical alias table introduced in cb292918b), which matches the singular AND the plural. **Alias only — nothing about the emitted response changes**: same `responseID` (3/4), same `ResponseType` labels, same `Vector(6*4)` width, so COLUMN_MAP layouts and MPCO/STKO readers are untouched. Gated by `tests/test_recorder_silent_drop.py::test_tet10_accepts_the_singular_token` (asserts 24 values for the singular spelling AND payload-identity with the plural); `tests/test_tet10_response_size.py` re-passes unchanged. | [#715](https://github.com/nmorabowen/OpenSees/pull/715) |
| `SRC/element/tetrahedron/FourNodeTetrahedron.cpp` | `// Ladruno`: same singular-token defect and same fix as the `TenNodeTetrahedron` row above — `strcmp(argv[0], "stresses")`/`"strains"` → `LadrunoResp::is(argv[0], "stress")`/`"strain"`, plus the `LadrunoResponseTokens.h` include. Alias only; the branch still emits `responseID` 3/4 over `Vector(6)` (this element has a single Gauss point) with unchanged labels. | [#715](https://github.com/nmorabowen/OpenSees/pull/715) |
