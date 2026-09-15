# `adr92_f7/` — the commit-time companion failure, reproduced and cured

WP-99 (F7). Two small decks, both of which are **mechanism demonstrations, not
capacity measurements**: the mesh is coarse on purpose and the `-implex` reading
hazard (ADR 92 §8) applies to every number they print.

## `strip_quad_implex.py` — the reproduction

A plane-strain `LadrunoQuad -formulation bbar` half-strip footing on SANISAND
(CP1 / Gorini parameters), self-weight on, in the TIMs Workbench configuration:
`-implex` **without** `-implexControl`, `-maxSubsteps 1000`, `-Pmin 0.0101`,
prescribed-settlement push (`sp` + `LoadControl(+ds)`, so `ops_Dt` is
proportional to settlement — the campaign deck's own shape).

144 elements (12 × 12, `h = 0.5 m`), 6 m × 6 m half-domain, `B/2 = 1 m`.

```
PYTHONPATH=<worktree>/dist/bin python3.12 strip_quad_implex.py out.csv [nsteps]
```

Add `--control` to run the same deck with `-implexControl 0.02` instead.

| binary | result |
|---|---|
| `9c2f964ea` (before WP-99) | **40 / 40 steps "converged" and committed**, with **567** commit-time companion refusals swallowed. `strip_quad_implex_OLD.csv`. |
| WP-99 | the first capped commit refuses, commits nothing, latches — and the **next** step's update is refused, so `analyze()` returns nonzero and the run stops. `strip_quad_implex_NEW.csv`. |

The OLD column is the defect in one line: the counter climbs from step 1 and the
load–settlement curve keeps rising anyway, because `Domain::commit()` drops every
element's `commitState()` return.

## `brick_implex_trajectory.py` — the bit-identity control

A free-DOF `LadrunoBrick` drained triaxial under `-implex` with an **adequate**
cap (`-maxSubsteps 20000`), i.e. a deck where the commit-time companion never
fails and WP-99's latch is therefore never armed. Writes one CSV row per
committed step (6 stresses, 6 strains, the six `implexDetail` slots, the four
`implexRefusals` counters) at `%.17g`.

Run it on both binaries and `fc` / `cmp` the files: byte-identical is the claim,
and the diff is the evidence. `brick_implex_trajectory_before.csv` /
`_after.csv` are the recorded pair.
