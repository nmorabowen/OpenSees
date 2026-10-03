---
wp: LEGACY
title: "Turning -implexControl ON is what walls a self-weight SANISAND push that TERMINATES fine without it — and control-off buys reach, not a confirmed curve"
legacy_seq: 432
---
### Turning `-implexControl` ON is what walls a self-weight SANISAND push that TERMINATES fine without it — and control-off buys reach, not a confirmed curve
- **Bites:** a self-weight SANISAND strip footing under `-implex -implexControl 0.05 0.01`
  refuses from the first push increments, the refusal counter runs into the
  hundreds (TIMs measured 204 120 on 9 720 Gauss points by step 3), every refusal
  halves the step, and the harness step falls below its floor at `s/B ~ 0.01`. It
  reads as a material wall — a dilatancy/fabric/`alpha_in` pathology, or the
  ADR-93 zero-confinement ring.
- **Measured (ADR-92 F10, `Ladruno_files/testbed/hypo_bearing/adr92_f10/`):** the
  SAME deck and mesh, identical in every other respect, with `-implexControl`
  simply **removed** — bare `-implex`, the doubling controller untouched —
  reaches the target `s/B = 0.05` in **104 steps, 0 subdivisions, 0 failed
  attempts, 58 s**, refusal ledger **0/0/0/0**. With the control on it walls at
  `s/B = 0.0085` with 724 refusals in 142 s. The ADR-95 campaign that reached
  `s/B = 0.15` on this material was itself control-OFF
  (`sanisand_path_diag.py:91`-`:104` passes only `-implex`).
- **What that does and does not establish.** It establishes TERMINATION, and it
  establishes it by measurement: the commit-time companion — the thing guide
  section 3 requires be able to FAIL rather than force-accept at `dT_min`, which
  is what `-maxSubsteps` buys — integrated every increment it was handed.
  `implexRefusals[3] = 0` on legs B, C, D, E, K, L, N and N1, and `<= 42`
  anywhere in the campaign (M 42, F1 29, I 12, H 6, F3 3, J 2). It does **not**
  touch the ACCURACY claim `-implexControl` also carries — ADR-92 section 8, and
  the engine's own control-off constructor echo at
  `LadrunoSANISAND.cpp:2092`-`:2095`: *"P0 measured IMPL-EX unusable from
  `d_eps = 5e-4` at `p0 = 5 kPa`, so at a low-confinement corner the control is a
  requirement, not an option."* **This deck sits inside that range, not outside
  it:** its minimum `p'` is **6.374 kPa, 1.27x the P0 corner**, and the
  control-off leg's strain increment crosses `5e-4` at step 24 (`s/B = 0.0012`)
  and runs at 2.6-4x the corner from `s/B = 0.00248` to the target — with **no
  implicit anchor past `s/B = 0.00227`**. So do not read this row as a general
  "the control is not needed at low confinement"; read it as "on this deck the
  control is what stops the run, its termination job is done by `-maxSubsteps`
  plus a watched bucket, and its accuracy job is simply not tested here."
- **Workaround/status (2026-09-14):** run bare `-implex` on a deck of this kind,
  **read `implexRefusals[3]` (the companion bucket) at the end of every leg**,
  and **do not quote the resulting curve as a capacity**. **Since WP-99 / PR #838
  (merged as `c75edc95c`, 2026-09-15) a capped companion commit ABORTS the run**
  (`Domain::commit()` fails, `analyze()` returns `-4`), so on any build from that
  merge on the bucket read is a second line of defence; the latch and its
  post-latch refusals live in `implexRefusals[4]`/`[5]`, not in the companion
  bucket `[3]`. **On a PRE-#838 build the read is the only thing that would catch
  it** — the F10 campaign itself ran on `9c2f964`, where such a commit is silent,
  which is the scope of every "silent" claim in this group. If the control is wanted
  anyway: `tol = 0.5` reaches the target with 18 refusals and agrees with the
  control-off arm to **0.395 % mean / 1.625 % max** over `0.002 <= s/B <= 0.05`
  — i.e. it makes the control nearly inert, which is the honest description of
  what it buys.
