# ADR-97 gate-4 baselines — `Backward_Euler` inertness

`be_secant_baseline_3622d6214.json` is the committed stress/strain history of
**23 decks** (282 committed-stress rows) taken from `dist/bin/opensees.pyd`
built at `3622d6214ef4cdeb8cf65a102ee35f6cd9973337` — i.e. **before** any ADR-97
C++ edit, on the ASDP source that `wp/97a-plan-oracles` inherited unchanged.

**The bits belong to the dumping host as well as to the commit (WP-175).**
The file was re-dumped on 2026-10-06 on the current desktop (AMD Ryzen AI 7 PRO
350, MSVC 14.44.35207, Windows SDK 10.0.26100.0) from local builds of
`3622d6214` (18 decks) and `7e93e4381` (the wp/97e deliberate entry plus the four
explicit-integrator decks, which need wp/97f's `experimental 1`). The previous
copy came from the workstation the fork used until 2026-09-18, and `3622d6214`
rebuilt here missed it by up to 2.0e-09 relative in 14 decks, with unchanged
codes and step counts. If the `==` gate fails on a different Windows host,
first rebuild `3622d6214` there and compare it against this file. If that build
misses too, the cause is the host, not `Backward_Euler`. See the quirk
`WP-175-adr97-gate4-baseline-is-host-specific`.

## Adding a host (WP-177, #960)

Byte identity holds per **host**, not only per commit (WP-175). `hosts.json`
maps `hostkey.host_key()` (platform + CPU brand; run `python hostkey.py` to see
it) to that host's own baseline file. On a registered host the gate is `==`; on
any other host, Linux included, it compares against
`be_secant_baseline_3622d6214.json` at a 1e-8 relative floor, and
`test_gate4_this_host_has_its_own_baseline` reports the host as SKIPPED with
these instructions.

To register a new Windows machine, once:

1. Build `3622d6214` there (all 5 targets via `build.bat`, in its own
   worktree) and dump the 18 pre-ADR-97 decks; build `7e93e4381` and dump the
   wp/97e entry and the four `experimental 1` decks. This is the WP-175 recipe.
2. Merge the two into `be_secant_baseline_3622d6214.<short-host>.json`, where
   `<short-host>` is a slug of the CPU (e.g. `intel-i9-13900k`). Every record
   carries `host`, so check that it names this machine.
3. Check that a build of current `ladruno` matches the new file `==` on all 23
   decks; ADR-97 D1 on this host is exactly that.
4. Add the row to `hosts.json`: `file`, `dumped`, `toolchain`, `by`.

`test_gate4_hosts_registry_is_consistent` fails if a registered file is
missing or covers a different deck set. A toolchain upgrade does not change a
host's key. If it moves the bits, re-dump that host deliberately and say so in
the commit; never re-dump to make a red gate green.

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
documented change to `Backward_Euler` lands, or (WP-175) when the dev host
changed and a rebuild of the baseline's own commit proves the move belongs to
the host. Never regenerate just to make a red gate green.

The children are started as plain `python3.12 dump_hist.py`, without `-S`. A
`site-packages` boot `.pth` that imports `opensees` therefore wins over
`PYTHONPATH`. Check that the `ladrunoBuild` recorded in each child dump names
the build you meant, or run the decks from a `-S` driver with pinned paths.
