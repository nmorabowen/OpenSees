# ZFP re-encode benchmark (ADR 08 gate)

Offline spike answering the one number ADR
[`08_zfp_compression_research.md`](../../Ladruno_implementation/08_zfp_compression_research.md)
hinges on: **does ZFP fixed-accuracy beat the current `f32 + shuffle + gzip(4)`
baseline by ≥2× at acceptable error, on the frozen `[T×nIds×nComp]` layout?**
Zero recorder code, zero build risk — everything runs in Python against files
written by the *installed* Ladruno build.

## Pipeline

1. **Inputs** — three models per the ADR's benchmark spec, run with the venv
   python (installed `opensees.pyd` via boot `.pth`), all `-precision f64`:
   - `model_a_nonlinear2d.py` — nonlinear transient, many steps: 2D J2 cantilever
     plate, 30×8 quads, **broadband** random-phase 0.5–25 Hz tip load (seeded),
     3000 steps. Disp/vel/accel/reactions + stress/strain.
   - `model_b_fiber_column.py` — fiber-section element history: 3D RC cantilever
     column, 5 forceBeamColumn × 5 Lobatto IPs, Concrete02+Steel02, growing
     cyclic push into yield, 2000 steps (adaptive dt). localForce + section
     force/deformation + fiber stress/strain. **The pessimistic arm.**
   - `model_c_brick_field.py` — large nodal field: 16×16×8 stdBrick block
     (2601 nodes), Ricker base excitation, 600 steps. Disp/vel/accel.
     **The optimistic (spatially-correlated) arm.**
2. **`check_inputs.py`** — gates the inputs: fiber arm must actually have
   yielded (max |fiber strain| > ε_y), else the pessimistic case is vacuous.
3. **`reencode_bench.py`** — the measurement. Full 2×2 codec×layout matrix:
   gzip/ZFP (+SZ3) × frozen-3D/time-major, with **recorder-true chunking**
   (chunk budget in on-disk bytes, so f32 chunks span 2× the slabs — mirrors
   `Ladruno_Hdf5.h::chunkSlabsPerBlock`). Reports aggregate AND worst-family
   marginal ratio vs `f32_gzip4`; per-dataset breakdown in the JSON.

## Provenance

Harness hardened per an adversarial methodology review (Opus 4.8, 2026-07-06)
before any result was interpreted. The three blockers it caught and their
fixes: (1) f32 baseline was being re-encoded with f64 chunk geometry →
recorder-true chunking; (2) time-major prototype had no gzip control → full
codec×layout matrix; (3) aggregate ratio could mask per-family losses →
worst-family gate. Timing numbers are end-to-end Python wall-clock, NOT codec
throughput — they must not be quoted against the ADR's literature "20–60×
faster than gzip" claim.

## Decision rule (from the ADR) — and outcome

Adopt `-compression zfp` only if the **worst representative family** clears
≳2× marginal over `f32_gzip4` at acceptable range-relative error.

**Outcome (2026-07-06): REJECTED.** No config beat the baseline at matched
fidelity on either layout (best: 0.65× frozen / 0.97× time-major). The
"<1.5× ⇒ layout is the bottleneck" presumption was *refuted* — the frozen
interleaved layout favors shuffle+gzip (fiber families ~13:1), so
v2-for-compression is closed too. Full table, eight caveats, and reopening
criteria: ADR 08 Implementation log. Raw numbers: `results.json`.
