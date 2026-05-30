# Zone B — fork-local validation (NOT upstreamable)

Everything here needs heavy deps (apeGmsh, gmsh, LS-DYNA refs, perf baselines)
and **must never be included in a PR back to OpenSees/OpenSees**. Zone A
(`tests/`, `EXAMPLES/verification/`) is what travels upstream.

Contract: [`Ladruno_implementation/testbed/00_canonical_testbed.md`](../../Ladruno_implementation/testbed/00_canonical_testbed.md).

```
Ladruno_files/testbed/
├── perf/                 Axis 3 — state-determination ns/call micro-bench
│   ├── runner.py         median-of-7, MKL threads=1, warn +10% / fail +25%
│   └── baselines/        per-machine JSON baselines (re-baseline via explicit commit)
├── <feature>/            per-feature T2b (apeGmsh sweeps) + T3 (cross-code) cases
└── README.md
```

**Run:** `pytest -m zone_b tests/` (the cases live under `tests/` with the
`zone_b` marker; their fixtures/data live here). They are auto-skipped when
gmsh/apeGmsh are absent and run with `APEGMSH_SKIP_VIEWER=1` (conftest).

**T3 author checklist (so a case can't false-green):**
- re-declare loads/masses/SP on the bridge (`ops.fix`/`ops.mass`/`p.load`/`p.sp`)
  — `apeSees` does NOT auto-emit them; MP constraints it DOES.
- end every case with `assert_equilibrium(...)`.
- beam internal forces → native HDF5 path only (no beam vecxz in MPCO on this build).
- recorder `file=` → a `tmp_path`, never the OneDrive-synced worktree.
- assert on PG-scoped scalars (tip disp, max von Mises, period), never raw gmsh
  node ids; pin the gmsh version.
