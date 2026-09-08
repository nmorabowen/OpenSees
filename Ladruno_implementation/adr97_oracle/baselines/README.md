# ADR-97 gate-4 baselines — `Backward_Euler` inertness

`be_secant_baseline_3622d6214.json` is the committed stress/strain history of
**23 decks** (282 committed-stress rows) taken from `dist/bin/opensees.pyd`
built at `3622d6214ef4cdeb8cf65a102ee35f6cd9973337` — i.e. **before** any ADR-97
C++ edit, on the ASDP source that `wp/97a-plan-oracles` inherited unchanged.

ADR-97 **D1** keeps `Backward_Euler` byte-identical. Gate 4 is that promise made
checkable: `tests/test_adr97_p4_inertness.py` re-runs the same decks in FRESH
SUBPROCESSES against the current binary and compares bit for bit.

Fresh subprocesses are not optional: the per-tag `INT_OPT_*` / `DBL_OPT_*` option
maps are `std::map<int,...>` keyed by material tag and shared across every
instance of that tag *by design*, and the pytest heap leaks state across tests in
one process (ADR-94 `capfd` / `WinError 6`).

## Decks covered

* 9 `TenNodeTetrahedron` decks — MohrCoulomb (default / `strict_convergence 1` /
  `n_max_iterations 200`), VonMises (default / strict), MohrCoulombTensionCutoff,
  MC-from-the-MCTC-file, HoekBrown, DruckerPrager.
* 10 `stdBrick` VonMises decks — `Backward_Euler` x {Continuum, Secant, Elastic,
  Numerical_Algorithmic_FirstOrder} x {plastic leg, elastic leg}, plus the
  softening (`H = -120000`) and perfectly plastic (`H = 0`) legs.
* 4 `stdBrick` VonMises decks on the explicit integrators (`Forward_Euler`,
  `Forward_Euler_Subincrement`, `Modified_Euler_Error_Control`,
  `Runge_Kutta_45_Error_Control`) — these share the YF/PF/hardening headers that
  ADR-97 P1 adds members to, so they are part of the inertness claim.

## Regenerating

From `<worktree>/tests`:

```
PYTHONPATH=../dist/bin LADRUNO_OPENSEES_QUIET=1 \
    python3.12 ../Ladruno_implementation/adr97_oracle/baselines/dump_hist.py out.json
python3.12 ../Ladruno_implementation/adr97_oracle/baselines/cmp_hist.py \
    be_secant_baseline_3622d6214.json out.json
```

`cmp_hist.py` prints `IDENTICAL` / `DIFFERS` per deck with the worst absolute and
relative difference. Regenerate the baseline **only** when a deliberate,
documented change to `Backward_Euler` lands — never to make a red gate green.
