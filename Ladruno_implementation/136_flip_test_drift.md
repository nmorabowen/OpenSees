# WP-136 — why `test_ladruno_sanisand_flip_determinism.py`'s push-step tests went red

Status: DONE (draft PR #870). Test-only change; no C++ touched. Hand-over from
WP-128 (`128_sanisand_ring_trace.md` §7, draft #869).

## Verdict

**Fix consequence, not a regression.** WP-110's tangent fix `dee04dbe3` moved the
Newton path of a deck whose first push step has **no reachable equilibrium**, so
the "converged" value pinned by WP-112 was wherever the pre-fix tangent happened to
stop. The `-3` on the first push step is expected under the fixed tangent. It is
not a separate regression.

## Evidence (all on this worktree's builds, CPython 3.12 `-S`, 1 MKL thread unless stated)

1. **Reproduced** on the `ladruno` tip `877ee112`: the two tests fail with rc -3 on
   the first push step (0..3 holds). The third test in the file passes.
2. **Merge skew.** The pins were measured on `48c0e99bc` (2026-09-16), which
   contains neither WP-110 nor WP-112. #847 (WP-110) merged at 20:24 on 09-18. WP-112
   merged `ladruno` into its branch at 20:26 and landed (#849) at 20:39 without
   re-measuring.
3. **Candidates** in `48c0e99bc..ladruno` on the deck path. LadrunoQuad, Pardiso,
   KrylovNewton, DisplacementControl and Transformation are untouched. That leaves
   `7b81e7fde` (WP-104, diagnostic reset), `dee04dbe3` (WP-110) and `25136af9c`
   (WP-112, default token only; the pins were taken with `init` explicit).
4. **Attribution: tip with `dee04dbe3`'s C++ reverted** (one-file rebuild; this
   replaces building `dee04dbe3` and its parent). The first step returns rc 0 at
   **9.659111** (`0x1.351770cc471a5p+3`) in **96 of 100** iterations. Steps 1–4 match
   the docstring to 1e-6. Steps 5–10 drift about 0.3% from it even here, which is
   ULP-level binary differences amplified by the chatter (the docstring's own caveat).
5. **WP-110 cannot move the equilibrium.** Under IntScheme 1, `GetElastoPlasticTangent`
   only feeds `aCep1/aCep2 → aCep_Consistent` (TanType 2), never the stress or α update.
6. **The step has no equilibrium under any tangent.** First push step, `NormUnbalance`:

   | variant | residual after 300–400 its |
   |---|---|
   | tip, TanType 2 | 0.20–0.35 kN, rc -3 |
   | tip, TanType 1 | 0.27–1.07 kN, rc -3 |
   | tip without WP-110, TanType 2 | stuck, rc -3 |
   | tip, TanType 2, `-honorTolR 1` TolR 1e-8 | 0.12–0.13 kN, rc -3 |

   `NormDispIncr 1e-8` "converges" wherever a tangent's increments get small:
   9.626 (TanType 0), 9.647 (TanType 1), 7.239 (TanType 2 at TolR 1e-8), and
   9.731 at `NormDispIncr 1e-6`. The floor is ModifiedEuler's non-smooth σ(ε)
   (WP-134 mechanism F). The TolR sensitivity is its discretisation error (WP-134 U9).

## Change

- The push step runs `test FixedNumIter 20`, and no load value is asserted. F14 is a
  determinism property: identical arithmetic gives identical bits whether or not
  Newton converged.
  - Measured on `877ee112`: 10 push steps bit-identical at MKL threads 1/2/4/8.
  - The default's first step spreads 2e-14 over 0..3 holds; vanilla's reads
    4.084 / 6.994 / 8.264 / 9.530 (0.57).
- `_run_child` asserts the child's `ladrunoBuild()` for every leg. WP-124's battery
  launched the children with no `PYTHONPATH`, got `ModuleNotFoundError`, and read it
  as three failures.
- **Break-test:** mapping the default to `vanilla`:
  - fails the hold-count leg (spread 0.5715);
  - fails the warning leg;
  - does not fail the thread leg. That is the limit the docstring already states:
    this ~170-DOF deck's Pardiso arithmetic is thread-invariant.

## Follow-ups (not done here)

- Once WP-129's SAS-ME (#871) ships, re-run this deck under it. If the residual floor
  disappears, the push can go back to a residual-converged step.
- WP-110's tangent defects and ModifiedEuler's unassigned `aCep` are also on
  upstream `master` (`93f7e8e58`). Upstreaming is out of scope (owner decision).
