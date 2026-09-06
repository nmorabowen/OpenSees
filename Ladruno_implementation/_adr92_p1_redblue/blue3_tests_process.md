# BLUE TEAM 3 — response to RED-3 (test battery, BVP gate instrument, process compliance)

Worktree `C:\Users\nmb\Documents\Github\OpenSees\.claude\worktrees\adr92-p1-implex`, branch
`wp/92c-implex-p1`, range `235934592..70f6bf0e8`. Read-only; no build, no pytest.

Scorecard: **F3 CONFIRMED · F4 PARTIAL (central consequence REFUTED) · F5 CONFIRMED ·
F6 CONFIRMED · F7 CONFIRMED · F8 CONFIRMED · F9 CONFIRMED · F10 CONFIRMED (with the
mitigation red itself states) · F11 PARTIAL · F12 CONFIRMED.**
Coverage matrix: rows 2, 6, 10, 17, 18, 20 upheld; **row 9 REFUTED as written**.

---

## F1/F2 — source-reading check only (owned by another defender)

RED-3's quotations match the source, with one mislocation. `mImplexDtCommit > 0.0` is
`LadrunoSANISAND.cpp:1329-1332`; `mImplexDt0 <= 0.0 && dt > 0.0` is `:1326-1327`;
the `-implexControl` escape `mImplexDt0 <= 0.0 || ...` is `:1615` — all verbatim. But the
D2 guard red calls "line 1615-ish" is at **`:1511`**, not 1615; red conflated two sites in
its "two consumers disagree about the sign contract" sentence, though the disagreement
itself is real (`:1511` tests `!= 0.0`, `:1329` tests `> 0.0`). `grep -n mImplexDtCommit`
shows the only write outside init/recv is `:1677 mImplexDtCommit = mImplexDt;` — no sign
flip anywhere, so red's stated refutation path for F1 does not exist. The `alpha = 1`
echo is confirmed: `adr92_bvp/implex_ctl/a2_h1.0_e0.6944_engine.log:3`.

---

## F3 — gate metric blind to subdivided steps — **CONFIRMED**

All three sub-claims hold, and the third holds by a *different* mechanism than red gives.

1. **`steps` is converged-only.** `sanisand_tau0_band.py:755` `rows.append(...)` sits after
   the `if not ok:` block at `:730-746`, whose every branch either `break`s or `continue`s
   (`nsub += 1; ds *= 0.5; ... continue`). `:915 steps=len(rows)`. A step that fails all
   three rungs never becomes a row.
2. **`past1` is blind to it.** `adr92_bvp_gate.py:82`
   `past1=(st - rung1) / st * 100.0` with `rung1 = st - rung2 - rung3` (`:69`); `ns`
   appears only in `rung2` (`:68`) and in `fail_it` (`:76-78`). Recomputed by importing
   `decompose` on the committed payloads: `implex_ctl` → `rung1=110, rung2=0, rung3=0,
   nsub=25, past1=0.0`; `implex_noctl` → `past1=0.0`; `implex_OFF` → `past1=48.4`. Arm A
   seized at `s/B = 0.00504` on `mode="FLOOR"` after 25 subdivisions and scores the best
   possible value on the headline metric.
3. **PREDICTION MET is reachable by a leg that dies early — but not arm A's way.** With
   `refused = 0` (the shipped reality, see F8) `failshare` *does* charge subdivisions, so
   arm A returns 85.0 % and fails `PASS_FAILED_SHARE`. The real hole is that
   `decompose()`/`table()`/`main()` read **only** `steps, nfail, nsub, nrelax, implex,
   build, tag` — never `mode`, `s_end_over_B`, `sfrac_target` or `verdict`. A leg that
   converges three clean rung-1 steps and then stops on `mode="WALL"` scores
   `past1 = 0.0, failshare = 0.0` and prints `VERDICT: PREDICTION MET` (`:136-138`).
   Worse: `load()` at `:59` filters `if "tag" in d and d.get("steps")`, so a leg with
   **zero** converged steps is silently dropped from the table rather than counted against
   the arm — survivorship, not just blindness. (That is exactly the pre-fix `implex_ON`
   leg: `run.log` "only 0 converged steps", no `.json` emitted at all.)
   Also confirmed: no `rung1 >= 0` / `rung2 >= 0` guard, and `RUNG_CAPS = (25, 40, 60)`
   (`:44`) is a hardcoded duplicate of the driver's `ladder` literal
   (`sanisand_tau0_band.py:697-700`) with no cross-check.

**Does the memo rescue it? No — it concedes it twice.**
`_adr92_p1_bvp_gate_results.md` §2 disclaims arm A's 85.0 % *on the ground that*
"no step in that leg left rung 1", i.e. it uses the blind metric to discard the one metric
that saw the subdivisions, leaving `past1 = 0.0` as the reported result. The memo does
disclose `subdivisions 25/80` and `FLOOR @ 0.00504` in its own table, so it is not hiding
the fact — but the instrument's headline number is unearned, and the commit body's
"NOT ONE STEP in either IMPL-EX leg left rung 1" is a statement about 110 rows chosen by
the same filter that discarded the 25 hard steps.

---

## F4 — provenance `build:'any'` — **PARTIAL; the central consequence is REFUTED**

**What holds.** `sanisand_tau0_band.py:785` and `:887` write `build=EXPECTED_BUILD` —
the *expected* value from `os.environ.get("LADRUNO_A2_EXPECT_BUILD", ...)` at `:234-235`,
not the measured `ops.ladrunoBuild()`. `:290-292` skips the assertion entirely when the
value is `"any"`. Measured from the committed JSONs: `implex_OFF.build =
'80e65a4de9bd6741261332b6f92de3a84bb4c642'`, `implex_ctl.build = 'any'`,
`implex_noctl.build = 'any'`. A payload schema that records a configuration echo in a
field named `build` is a genuine defect, and the fail-loud assertion was off for both
IMPL-EX arms.

**What is refuted — and it is red's own stated refutation condition.** Red asks for "a
recorded `ladrunoBuild()` return value (not `EXPECTED_BUILD`) for the two implex arms".
It is committed, in the same commit, one directory away:

```
adr92_bvp/implex_ctl/run.log:1    *** LADRUNO_A2_EXPECT_BUILD=any: build hash NOT pinned
                                      (running 8fe4f5630be5097454751f568bd423c281b8bbd0) ***
adr92_bvp/implex_ctl/run.log:3    engine build (ladrunoBuild) : 8fe4f5630be5097454751f568bd423c281b8bbd0
adr92_bvp/implex_noctl/run.log:3  engine build (ladrunoBuild) : 8fe4f5630be5097454751f568bd423c281b8bbd0
adr92_bvp/implex_OFF/run.log:2    engine build (ladrunoBuild) : 80e65a4de9bd6741261332b6f92de3a84bb4c642
adr92_bvp/implex_ON/run.log:2     engine build (ladrunoBuild) : 80e65a4de9bd6741261332b6f92de3a84bb4c642
```

So (i) the arms are **identified**, not unidentified; (ii) the commit's provenance
sentence ("the binary ... REPORTS `8fe4f5630`") matches the two implex `run.log`s exactly
— there are not "three mutually inconsistent provenance claims", there are two arms whose
measured hash agrees with the memo and a control whose measured hash agrees with its own
payload; (iii) the `any` path is not silent — it prints a loud unpinned-run warning as
line 1 of every run.log. The `.pyd` path is also recorded on the next line of each log,
and it is this worktree's `dist\bin\opensees.pyd` in all four runs.

**"Verified behaviourally" is sufficient here, narrowly.** The discrimination the memo
needs is pre- vs post-`3c788778f`, and the committed logs supply it directly: the pre-fix
run (`implex_ON`, build `80e65a4de`) contains **34,272** identical negative-`dt` refusals
and 0 converged steps; the post-fix arms contain **zero** `LadrunoSANISAND ... REFUSES`
lines and reach 110/142 steps. A pyd stale at `8fe4f5630` would carry the pre-fix guard
and could not have reached TARGET on a monotone-negative deck, so the stale-`.pyd` failure
mode of `31c4b5de7` is ruled out by the data, not merely asserted. `3c788778f` is the last
C++ commit in the range, so there is no "later build" to confuse it with.
One residual: the memo's cited count, "the negative-`dt` refusal fired **4896** times
before", matches nothing in the committed evidence (the committed pre-fix log has 34,272,
exactly 7×4896); that specific number is unverifiable from the repo.

---

## F5 — the battery cannot detect a no-op operator — **CONFIRMED**

The mutant-survival table is right on every row I could check, and its two strongest legs
are conceded by the file itself.

* **Tangent-identity affinity — CONFIRMED, and decisive.** `sigma~` is built at
  `LadrunoSANISAND.cpp:1563-1568` as `mSigma_n + Ce*(d_eps - f*mImplexDEpsP)`. The test
  (`tests/test_ladruno_sanisand_implex.py:596-607`) forms
  `d_sigma = sp - sm`, `d_eps = ep - em` over two probes that share the same committed
  state, the same `mImplexDEpsP` and (by the file's own norm argument, `:445-458`) the same
  `f`. The `f*mImplexDEpsP` term is identical in both probes and cancels exactly, so the
  secant is `Ce` for **any** frozen `f`, `f = 0` included, and for a total-strain form
  `Ce:(eps - f*eps_p)` likewise. The file states this itself at `:456-458`:
  "EXACTLY affine in `d_eps` for FIXED `f`".
* **`implexDetail[5]` is unpacked and dropped.** `:376`
  `total, dev, vol, fired_last, count, f_last = detail`; `grep -n f_last` returns that
  line and nothing else. The one Python-visible observable of `f` is never asserted on.
* **Zero-free-DOF blindness.** `sani._build` (`tests/test_ladruno_sanisand.py:349-371`)
  constrains all 24 DOFs (`fix` on the `x==0 / y==0 / k==0` faces, `sp` on `x==1 / y==1 /
  k==1`), so tests 2, 3 and 13 compare an ON path to an OFF path on a deck the file argues
  in two places cannot distinguish them. The flagship test's docstring
  (`:190-197`) says so outright: "a claim about the CODE PATH, provable without a binary".
* **`80e65a4de` shipped mutant (a) in production and added no test.**
  `git show --stat 80e65a4de` = `LEDGER_quirks.md`, `_adr92_p1_execution_plan.md`,
  `LadrunoSANISAND.cpp`. No `tests/` file. The defect it fixed was "every retried step in
  an adaptive run silently ran with the plastic extrapolation switched off".

**The floor-clamp question, answered.** Red hedges "probably red"; the evidence says the
clamp test is the one live constraint on the operator but that its firing is driven by the
extrapolation term, so `f≡0` most likely turns it red — and the file *measured* that.
`_build_floor_seeking_deck`'s docstring (`:274-296`) records that the first-draft isochoric
deck measured `implexDetail` count == 0 on the real binary, and diagnoses it as exactly the
`f`-term-absent case: "`mImplexDEpsP` is exactly zero until the FIRST commit after the
stage flip ... so `sigma~` ... is a much MILDER elastic-only prediction ... mild enough
here to stay clear of the floor". That is a direct measurement that an elastic predictor
(the `f≡0` mutant) does **not** cross `Pmin` on that deck. The current deck differs only by
`lat 0.5 → 1.5` (net dilation), i.e. it was reshaped until the clamp fired
(`17ccc6075`), which is red's other point. So: the clamp test is the battery's only live
operator constraint, its discriminating power was obtained by tuning the deck after a
`count == 0` measurement, and no test pins the count or `f` itself. `PARTIAL` only on
red's own hedge — the mechanism is confirmed, the direction is confirmed by the file's own
measurement, but neither of us has run it.

---

## F6 — db roundtrip skipped and structurally blind — **CONFIRMED**

`_roundtrip_implex` (`:824-859`) calls `sani._build` at `:832` and `:847`; `ref` at `:884`
and `ref_off` at `:900` call `sani._drive`, which is `_build` + legs
(`test_ladruno_sanisand.py:394-399`). All are the zero-free-DOF deck (see F5). By the
file's own gate-5 argument the ON and OFF committed answers are bit-identical there, so
`reldiff(ref_off, ref) == 0`, the non-vacuity guard at `:902-906`
(`if gap < _SENSITIVITY_FLOOR: pytest.skip(...)`, floor `1e-3`) fires, and the only
assertion that could see a dropped flag never executes. The two assertions that do run
(`:890`, `:894`) compare against `ref`, which a fully `-implex`-stripped restore would also
reproduce. Run records confirm the skip: `8fe4f5630` "13 passed / 2 failed / 1 skipped",
`17ccc6075` "11 passed / 4 failed / 1 skipped"; the only other `pytest.skip` sites are
`:842`/`:845` (`database()` unsupported), which fire before the assertions. `Vector(5) →
Vector(22)` (17 new slots, `LadrunoSANISAND.cpp:924-926`, `:987-989`) ships with no test
that can fail. No MPI `sendSelf`/`recvSelf` test exists.

---

## F7 — D2 semantics changed after "13/15", no test — **CONFIRMED** (one sub-claim softened)

Times and footprints, all verified:
`8fe4f5630` 2026-09-06 01:15:46 (`tests/…implex.py` only, 25+/8-) → `3c788778f`
01:53:58 (`SRC/material/nD/LadrunoSANISAND.cpp` only, 27+/9-) → `31c4b5de7` 01:55:24
(`LEDGER_quirks.md` only) → `70f6bf0e8` 02:25:51 (data + gate + memo). **No commit after
`8fe4f5630` touches `tests/`.** So the "13/15" figure was measured two behaviour commits
before the shipped semantics, and no post-fix pytest result is recorded anywhere in the
range or in `_adr92_p1_execution_plan.md`'s Log (whose last entry is the 11/4/1 run).

The guard moved from `mImplexDt < 0.0` to `mImplexDtCommit != 0.0 && mImplexDt *
mImplexDtCommit < 0.0` (`:1511`). The only negative-clock probe in the battery is
`ops.integrator('LoadControl', -1.0)` at `:746`, taken **after** `_establish_plastic_history`
has committed positive-`dt` steps (`:530-534`) — so it is a sign change and still refuses
under the new semantics. Red's sub-claim 2 ("the test name and docstring are now false")
is therefore **PARTIAL**: the test still passes and still exercises a shipped refusal, but
its name and its docstring's stated contract ("a NEGATIVE pseudo-time increment (`ops_Dt <
0`) is refused") describe a rule that no longer exists. Sub-claims 1, 3 and 4 stand:
no test drives two consecutive negative steps and asserts `analyze() == 0`
(`grep "LoadControl'," ` returns one negative site, `:746`), and the ledger row still reads
"a **negative** pseudo-time refused by symptom".

---

## F8 — cost-model edit inert — **CONFIRMED, in full, including the reconstruction**

`git show 70f6bf0e8 -- .../adr92_bvp_gate.py` is the diff red quotes, inside the results
commit. Thresholds were genuinely pre-registered: `git log --diff-filter=A` and
`git log -S REFUTE_PAST_RUNG1` both return **`5e5ec4db7`** (00:44:59), 1h41m before the
results commit — that half of the process holds.

`refused` is sourced only at `:93` `int(d.get('n_material_refused', 0) or 0)`.
`grep -rn n_material_refused` over the repo returns three hits: that line and two
"Owed"/"still not recorded" lines in `_adr92_p1_bvp_gate_results.md` (`:87`, `:116`).
No producer, and the key is absent from all three leg JSONs (verified by dumping
`sorted(d.keys())`). So `refused = 0` and `ladder_fails = max(ns - 0, 0) = ns`: the new
expression reduces to the old one identically. Arm A's `failshare` is still **85.03 %**,
above `PASS_FAILED_SHARE = 10.0`.

The reconstruction is exact. `implex_ctl/..._engine.log` (394 lines) contains, by
line-shape count: 75 × `Domain::update - domain failed in update`, and **25 each** of
`NewtonRaphson::solveCurrentStep`, `NewtonLineSearch::solveCurrentStep`,
`AcceleratedNewton::solveCurrentStep` — 3 × 25 = 75 = `nfail`, with `nsub = 25`, so every
failed rung belongs to a subdivided step and `refused = 25` is the correct substitution.
`decompose(d, 25)` → `failshare = 53.19 %`. Red's 53.2 % is right; the number was
recoverable from committed evidence and was disclaimed instead.

---

## F9 — `getCopy` test vacuous — **CONFIRMED**

`LadrunoSANISAND::getCopy(const char *type)` (`:809-843`) constructs
`LadrunoSANISAND3D` / `LadrunoSANISANDPlaneStrain` from scalar parameters only and then
`clone->setLadrunoImplexOptions(mImplexOpt, false)`. There is no history-copying statement
to mutate. Ordering confirmed in the test: `ops.element('stdBrick', 1, ...)` and
`ops.element('stdBrick', 2, ...)` at `:941-942`; the first `ops.analyze` is at `:980`.
Both clones are taken at zero committed strain, so `bystander_pstrain == 0` is satisfied
by any correct material. The assertion the ledger cites as coverage of the seam — that the
implex **options** are transferred — is untested; the sibling
`test_getcopy_void_carries_the_settings` uses the stronger construction this test does not.

---

## F10 — warrant package — **CONFIRMED (four gaps), with red's own mitigation upheld**

1. **Status stale.** The added row (one line, `git diff 235934592..70f6bf0e8 --
   LEDGER_implementations.md` shows `1 +` and no modification) reads "**draft — C++
   written, NOT YET BUILT.** No gate has run: oracle parity, tangent identity, floor
   clamp, byte-identity, symmetry, `dt` guards and the BVP gate all wait on the owner's
   `build.bat OpenSees OpenSeesPy`." Six later commits report gate results on a built
   binary, two of which changed the C++.
2. **Files column incomplete.** It lists `SRC/material/nD/LadrunoSANISAND.{h,cpp}`,
   `LadrunoSANISAND3D.cpp`, `LadrunoSANISANDPlaneStrain.cpp` and nothing else. The range
   also adds `tests/test_ladruno_sanisand_implex.py` (998 lines),
   `Ladruno_files/testbed/hypo_bearing/adr92_bvp_gate.py`,
   `_adr92_p1_{execution_plan,bvp_gate_results}.md`, `LEDGER_quirks.md` (+136) and 37 data
   files. `CLAUDE.md` requires the row to name the files.
3. **Behaviour claim wrong** post-`3c788778f` — see F7.
4. **Manifest untouched, CI cannot notice.** `git log -3 -- Ladruno_implementation/testbed/
   manifest.yaml` → `c89005e2a`, on `ladruno` before this branch. `ci/check_manifest.py`
   iterates `ladruno_tags()`, which regexes `#define` lines in `SRC/classTags.h` whose
   trailing comment contains "ladruno" (`:23-30`), then requires a manifest row per
   *tag_symbol*. ADR-92 P1 adds no classTag (the row says "none new"), so the check passes
   without looking at the feature. `manifest.yaml:913`'s `ND_TAG_LadrunoSANISAND` row still
   names only `tests/test_ladruno_sanisand.py`, and its coverage evidence is
   mutation-proven throughout ("PROVEN by mutation: ... fails THIS TEST AND ONLY THIS
   TEST"). `Ladruno_scripts/mutation_gate.py` exists (15,067 bytes, mtime pre-branch) and
   `grep -ril implex Ladruno_scripts/ ci/` returns only two unrelated recorder-harness
   files — no ADR-92 mutation run.

**Mitigation upheld, exactly as red states it.** `gh pr view 798`:
`isDraft: true`, `baseRefName: ladruno`, `createdAt 2026-09-06T04:28:02Z`; first commit
`02b1a7c8b` at `2026-09-05T23:27:34-05:00` = `04:27:34Z` → **28 s**, day-one draft,
compliant. Banner untouched (`git diff --stat` over the range shows no
`Ladruno_scripts/banner_features.txt`, `SRC/tcl/tclMain.cpp` or
`SRC/interpreter/PythonModule.cpp`), matching "No banner line until CP2". Vanilla footprint
zero: the range touches only `SRC/material/nD/LadrunoSANISAND.{h,cpp}`,
`LadrunoSANISAND3D.cpp`, `LadrunoSANISANDPlaneStrain.cpp` — all fork-owned — so no
`LEDGER_vanilla_files.md` row is owed. ADR-87 requires the complete warrant package at the
ready-flip, not on a draft.

---

## F11 — red tests unmarked, 16 vs 15 — **PARTIAL**

**Confirmed.** AST walk of `tests/test_ladruno_sanisand_implex.py` returns 12
`test_*` functions, of which one carries
`@pytest.mark.parametrize('scheme', [3, 5, 7, 8, 9])` → **16 collected**, not 15.
`grep "xfail"` returns nothing; `grep "pytest.skip"` returns only `:842`, `:845`, `:903`.
`8fe4f5630`'s message reports "13 passed / 2 failed / 1 skipped", so the published "13/15"
counts against a denominator that omits the skipped test — which, per F6, is the one
hiding a structural blindness. The file carries `pytestmark = [pytest.mark.zone_a]`
(`:88`), so those tests are in the Zone-A selection, and `.github/workflows/ladruno.yml:54-59`
skips the Zone-A build on drafts (`github.event.pull_request.draft == false`), so nothing
runs them.

**Softened on two points.** (a) Red cites `.github/workflows/build_cmake.yml:43-48` as the
pytest gate; that workflow triggers only on `master` (`:4`, `:7`) and this PR targets
`ladruno`, so it never runs on this branch at all. The correct citation is
`ladruno.yml:104-109` (`pytest -v -m "zone_a"`). (b) Red names the two failures as
`test_tangent_identity_...` and `test_implexcontrol_refuses_past_tolerance...`. The record
does not say: `8fe4f5630` states only "Both remaining failures are under diagnosis and
attributed to neither side yet", and that same commit *deleted* the three `ops.reset()`
calls that `_adr92_p1_execution_plan.md:218-226` diagnoses as the tangent test's cause.
Which two are red is unrecorded — which is the same evidentiary hole F7 identifies.

---

## F12 — tuned constants, unchecked precondition, 12.4 MB log — **CONFIRMED**

1. `_CAP_ADEQUATE = 20000` at `:144`, under a comment block (`:120-143`) that states it was
   read off "this build (2026-09-06)". It *is* applied as a controlled variable on every
   ON/OFF pair; no test pins the substep count, so a regression moving the required count
   above 20000 silently restores the unequal-cap comparison `5e5ec4db7` fixed.
2. `_PROBE_P0 = 100.0` at `:466-467`, "away from the p_min floor" by comment only. The
   tangent test reads `'stress'`, `'strain'` and `'tangent'` and never `implexDetail`, so
   it never checks that `implexDetail[3]`/`[4]` (clamp fired / count) stayed zero across
   the two probes — despite the new `LEDGER_quirks` entry stating "a test that reports
   'identity to machine precision' on a clamped path is measuring nothing", and despite
   `test_floor_clamp_fires_...` proving the observable is reachable.
3. `adr92_bvp/implex_ON/a2_h1.0_e0.6944_engine.log` = **12,407,322 bytes / 34,396 lines**;
   digit-normalised line-shape count gives **34,272** copies of one string, the D2
   negative-`dt` warning. It is the superseded pre-fix run (`run.log` "only 0 converged
   steps", `FAILED.txt` at `s/B = 0.00000`, build `80e65a4de`). Sibling post-fix logs are
   394 and 8 lines. No `.gitignore` rule covers `*_engine.log` (only
   `Ladruno_scripts/_*.log` and `build_p1b*.log`). The asymmetry is confirmed at source:
   the D2 warning (`:1512-1518`) has no throttle counter, while the `-implexControl`
   refusal (`:1615-1616`) returns `LADRUNO_MATERIAL_REFUSED` with **no `opserr` at all` —
   the `implex_ctl` log contains zero `LadrunoSANISAND ... REFUSES` lines and only
   `LadrunoBrick`'s own message, throttled at 10 with one suppression notice — which is
   precisely why `n_material_refused` (F8) is unrecoverable.

---

# COVERAGE MATRIX — spot checks

| row | red's verdict | blue |
|---|---|---|
| 2 — `-implex` on does not corrupt the committed answer | leak detector only | **UPHELD.** `test_implex_on_matches_off_on_a_zero_free_dof_deck` runs `sani._drive_confined`, zero free DOFs; its docstring (`:190-197`) calls the claim "provable without a binary". Survives every operator mutant. |
| 6 — frozen tangent `d(σ~)/dε == Ce` | RED (unmarked) + survives `f≡0` and total form; precondition unchecked | **UPHELD on substance, PARTIAL on "RED".** The affinity argument is exact (see F5); the precondition is unchecked (F12.2); no `xfail`. Whether *this* test is one of the two still-failing is unrecorded (F11b). |
| 9 — D2 sign-change refusal, **shipped** semantics | NOT COVERED | **REFUTED as written.** `test_implex_negative_pseudo_dt_is_refused...` drives `LoadControl(-1.0)` (`:746`) *after* positive committed steps (`:530-534`), so `mImplexDtCommit > 0`, `mImplexDt < 0`, product `< 0` — it exercises the shipped `:1511` guard directly and asserts `rc != 0`. The defect is that its name and docstring describe the retired rule, not that the shipped rule is uncovered. |
| 10 — D2 must **not** refuse a monotone negative clock | NOT COVERED | **UPHELD.** Only one negative `LoadControl` site in the file, and it is preceded by positive steps. Nothing drives two consecutive negative increments and asserts `analyze() == 0`. This is the regression `3c788778f` fixed. |
| 17 — PlaneStrain lane | NOT COVERED | **UPHELD.** `grep "'-ndm'"` returns four hits, all `3` (`:304, :483, :685, :931`); zero `-ndm 2`; `PlaneStrain` appears only in prose. `LadrunoSANISANDPlaneStrain.cpp` is +20 lines in the range. |
| 18 — first-step arming / `f` frozen once per step | NOT COVERED | **UPHELD.** `f_last` is bound at `:376` and referenced nowhere else. `80e65a4de` (the arm-consumption defect) added no `tests/` file. |
| 20 — companion commit-time failure (`-maxSubsteps` exhausted inside `ladrunoImplexCommit`) | NOT COVERED | **UPHELD.** The path is `LadrunoSANISAND.cpp:1641-1660` ("There is no good move here"), reachable only with `-implexControl` OFF. `_CAP_TIGHT = 200` is used at exactly one site, `:633`, inside `test_implex_refuses_unsupported_schemes` — a parser refusal, not a commit-time companion failure. No test pairs a tight cap with control off on a free-DOF deck. |

Rows not spot-checked are not endorsed; rows 22 and 24 were incidentally confirmed
(`grep implexAlpha` returns nothing in the battery; `'-implexDt'` is passed once, `:490`,
with `'strain'`).
