# RED TEAM 3 — test battery, gate instrument, process compliance
## ADR-92 P1 `-implex` on `LadrunoSANISAND` — PR #798, `wp/92c-implex-p1`, range `235934592..70f6bf0e8`

Read-only audit. 12 findings, most severe first, then the coverage matrix.

---

## F1 — BLOCKER — the headline BVP leg ran with the extrapolation factor pinned at α, because the campaign's own deck has a negative clock

**Where:** `SRC/material/nD/LadrunoSANISAND.cpp:1325-1332` (introduced `56c382ea1`, blessed by `3c788778f`), consumed by `Ladruno_implementation/_adr92_p1_bvp_gate_results.md` / commit `70f6bf0e8`.

**CLAIM.** Commit `3c788778f` relaxed D2 so that a *consistently negative* pseudo-clock is legal, on the stated ground that "two consecutive NEGATIVE increments give a POSITIVE, correct ratio". That ratio is never computed. Fifteen lines above the guard it changed:

```c
    mImplexDt = dt;
    if (mImplexDt0 <= 0.0 && dt > 0.0)
        mImplexDt0 = dt;

    if (mImplexDtCommit > 0.0)
        mImplexFactor = (dt / mImplexDtCommit) * mImplexOpt.alpha;
    else
        mImplexFactor = mImplexOpt.alpha;
```

**EVIDENCE.** The campaign driver drives settlement downward — `sanisand_tau0_band.py:721`, `ops.integrator("LoadControl", -ds)` — so `dt < 0` on every step of every BVP leg (the commit message says exactly this: "the pseudo-time increment is negative on EVERY step, by design"). Therefore `mImplexDtCommit > 0.0` is **false on every step for the life of the leg**, and `mImplexFactor` is the constant `mImplexOpt.alpha` (= 1, per the echo line in `adr92_bvp/implex_ctl/..._engine.log:3`, "alpha = 1"). The `dt_{n+1}/dt_n` adaptation — the thing that makes this IMPL-EX and not a fixed elastic predictor — never runs on any leg in `adr92_bvp/`.

**CONSEQUENCE.** The headline result of `70f6bf0e8` ("TARGET @ 0.25000", "the ladder is GONE") was produced by an operator that is *not the one the ADR describes*. Every conclusion in `_adr92_p1_bvp_gate_results.md` about IMPL-EX's cost behaviour is a statement about `f ≡ 1`, not about `f = (dt_{n+1}/dt_n)·α`. `3c788778f`'s own justification for permitting the deck is falsified by the code it left unchanged.

**WHAT WOULD REFUTE IT.** A leg record or engine echo showing `implexDetail[5]` (`mImplexFactor`) ≠ 1.0 on a `LoadControl(-ds)` deck; or a source path I missed where `mImplexDtCommit` is stored with a sign flip. Note `mImplexDtCommit > 0.0` is a strict positive test, not `!= 0.0` — unlike the D2 guard's own `mImplexDtCommit != 0.0` on line 1615-ish, so the two consumers of the same variable disagree about its sign contract.

---

## F2 — BLOCKER — `-implexControl` can never stop refusing on a negative clock; the memo blames the tolerance and books a sweep

**Where:** `SRC/material/nD/LadrunoSANISAND.cpp:1326-1327` and `:1615`; `_adr92_p1_bvp_gate_results.md` / `70f6bf0e8` commit body.

**CLAIM.** `mImplexDt0` (header: "the first non-zero dt seen -- the `-implexControl` floor") is armed **only from a positive dt**:

```c
    if (mImplexDt0 <= 0.0 && dt > 0.0)
        mImplexDt0 = dt;
```

On a negative-clock deck it stays `0.0` forever. The refusal escape then reads:

```c
        if (mImplexError > mImplexOpt.errorTol) {
            // Below the reduction limit there is nothing left to cut, so refusing
            // again would only turn a bounded inaccuracy into a dead analysis.
            if (mImplexDt0 <= 0.0 || mImplexDt >= mImplexOpt.reductionLimit * mImplexDt0)
                return LADRUNO_MATERIAL_REFUSED;
        }
```

`mImplexDt0 <= 0.0` short-circuits **true**, so the "nothing left to cut" brake is permanently disabled and the material refuses without limit.

**EVIDENCE.** Arm A (`adr92_bvp/implex_ctl/a2_h1.0_e0.6944.json`): `steps=110, nfail=75, nsub=25, nrelax=0, mode="FLOOR"` — exactly 3 failed rungs per subdivided step, 25 subdivisions, terminating on the `DS_MIN` floor at `s/B = 0.00504`. That is the signature of a refusal that never stops. The `70f6bf0e8` commit body instead concludes: *"With control ON at tol = 0.05 the leg SEIZES ... The two arms bracket it: the mechanism works, the tolerance as set does not"*, and books as owed next *"find the `-implexControl` operating point (0.05 is unusable; sweep)"*.

**CONSEQUENCE.** A code defect is recorded as a tuning finding, and the next work item is a parameter sweep that cannot succeed at any tolerance on this deck class. `-implexControl` is described in `LEDGER_implementations.md` as **mandatory at the low-p corner** ("unusable from `5e-4` at `p0 = 5`"), so the only accuracy police the feature has is inoperative on exactly the deck the ADR was opened for.

**WHAT WOULD REFUTE IT.** A leg with `dt < 0` throughout that shows `-implexControl` ceasing to refuse as `ds` shrinks; or a second arming site for `mImplexDt0` — `grep -n mImplexDt0` returns only `:1108` (zero), `:1326-1327` (positive-only), `:925`/`:988` (wire), `:1615` (this consumer), `:2014` (print).

---

## F3 — BLOCKER — the gate's headline metric is structurally blind to the failure mode that killed arm A

**Where:** `Ladruno_files/testbed/hypo_bearing/adr92_bvp_gate.py:66-79` vs `sanisand_tau0_band.py:708-748`.

**CLAIM.** `past rung 1` is computed over **converged steps only** and therefore cannot see a step that failed all three rungs and subdivided.

**EVIDENCE.** Driver: `rows` (hence `steps = len(rows)`) is appended **only after** `if not ok:` has been skipped; a step that fails all three rungs does `nsub += 1; ds *= 0.5; continue` and never becomes a row. Gate:

```python
    rung3 = nr
    rung2 = nf - 3 * ns - 2 * nr
    rung1 = st - rung2 - rung3
    ...  past1=(st - rung1) / st * 100.0
```

`nsub` appears in `rung2`'s correction but **never in `past1`'s numerator or denominator**. Verified on the committed payloads: arm A (`implex_ctl`) decomposes to `rung1=110, rung2=0, rung3=0 → past1 = 0.0 %` — a leg that seized at a tenth of the control's depth after 25 subdivisions reports the *best possible score on the gate's headline metric*. Arm B likewise 0.0 %.

**CONSEQUENCE.** The commit's load-bearing sentence — *"NOT ONE STEP in either IMPL-EX leg left rung 1"* — is true and meaningless: the metric awards 0 % to a leg every one of whose hard steps was thrown away rather than climbed. The gate's own `VERDICT: PREDICTION MET` branch (`ip <= 10 and ifs <= 10`) is reachable by a leg that dies immediately. There is no guard on `rung1 >= 0`/`rung2 >= 0` either, and the hardcoded `RUNG_CAPS = (25, 40, 60)` is not checked against the driver's `ladder` literal, so a CP1 baseline produced by any other ladder is silently mis-decomposed.

**WHAT WOULD REFUTE IT.** A definition of `steps` in the driver that includes failed attempts (it does not — `rows.append` is after the `if not ok` block), or a second metric in the gate that charges subdivisions to the ladder (only `failshare` does, and see F8).

---

## F4 — BLOCKER — provenance: both IMPL-EX arms were run with the fail-loud build check **switched off**, in the same session that documented the silent-stale-`.pyd` trap

**Where:** `sanisand_tau0_band.py:234, 290-292, 785, 887`; payloads under `Ladruno_files/testbed/hypo_bearing/adr92_bvp/`.

**CLAIM.** The `build` field in a leg record is the **expected** build read from the environment, not the measured `ladrunoBuild()`, and `"any"` disables the assertion entirely:

```python
EXPECTED_BUILD = os.environ.get(...)
...
        if EXPECTED_BUILD.lower() != "any":
            assert build == EXPECTED_BUILD, (...)
...
                build=EXPECTED_BUILD, driver=os.path.abspath(__file__),
```

**EVIDENCE (measured from the committed JSON):**

| dir | `build` | `implex` | mode |
|---|---|---|---|
| `implex_OFF` (control) | `80e35...` → `80e65a4de9bd6741261332b6f92de3a84bb4c642` | `False` | WALL |
| `implex_ctl` (arm A) | **`any`** | `True` | FLOOR |
| `implex_noctl` (arm B) | **`any`** | `True` | TARGET |

The three-arm comparison table in `70f6bf0e8` therefore compares one identified binary against two unidentified ones. The commit's own provenance sentence — *"the binary CONTAINS the D2 sign-guard fix but REPORTS `8fe4f5630`"* — matches **neither** payload: the control says `80e65a4de`, the implex arms say `any`. Three mutually inconsistent provenance claims for one table. And `31c4b5de7`, committed 30 minutes earlier in the same session, documents that this build system *silently ships a stale `.pyd` and exits 0* — the memory-recorded fork rule is "open every probe/harness with `ladrunoBuild()` (fail-loud wrong-build check + stale-pyd detector)".

**CONSEQUENCE.** The only evidence that the arms ran the intended code is the commit message's "verified BEHAVIOURALLY, not by hash". That verification ("the negative-dt refusal fired 4896 times before and 0 after") distinguishes pre- from post-`3c788778f`; it does **not** distinguish post-`3c788778f` from any later or intermediate build, and it cannot rule out the exact stale-`.pyd` failure the same session recorded.

**WHAT WOULD REFUTE IT.** A recorded `ladrunoBuild()` return value (not `EXPECTED_BUILD`) for the two implex arms, or the `dist/bin/opensees.pyd` mtime the quirk entry itself prescribes checking.

---

## F5 — BLOCKER — the battery cannot detect `-implex` being a no-op, nor the two operator errors P0 named; exactly one live test constrains the operator, and its deck was reshaped until it fired

**Where:** `tests/test_ladruno_sanisand_implex.py`, whole file.

**CLAIM.** Walk the 16 collected tests against three mutants — (a) the extrapolation inert (`sigma~ = sigma_n + Ce:d_eps`, i.e. `f ≡ 0`), (b) the **total**-strain form P0 measured 78 % wrong at `p0 = 5 kPa`, (c) `mImplexFactor` frozen at the wrong value:

| test | survives (a) `f≡0`? | survives (b) total form? |
|---|---|---|
| `..._unset_reproduces_recorded_sensitivity` | yes (no `-implex`) | yes |
| `..._on_matches_off_on_a_zero_free_dof_deck` | **yes — by its own docstring** | **yes** |
| `..._stage0_inertness_...` | yes | yes |
| `..._floor_clamp_fires_and_implexDetail_reports_it` | probably red | yes (clamp fires *more*) |
| `..._tangent_identity_...` | **yes** — `σ~ = Ce:(ε − f·ε_p)` is affine in `dε` with the *same* `Ce` for **any** frozen `f`, including 0, and for the total form too | **yes** |
| 5× `..._refuses_unsupported_schemes`, `..._requires_maxsubsteps`, `..._maxsubsteps_zero...` | yes (parser only) | yes |
| `..._negative_pseudo_dt_is_refused...` | yes (see F7) | yes |
| `..._implexcontrol_refuses_past_tolerance...` | currently **red** | — |
| `..._db_roundtrip...` | **SKIPPED** (F6) | — |
| `..._getcopy_does_not_share_history...` | **vacuous** (F9) | yes |

**EVIDENCE.** The flagship test's own docstring proves it cannot fail on a broken operator: *"`ladrunoImplexTrial()` (the extrapolation) never writes `mEpsilon` … So the OFF path's single `integrate()` call and the ON path's commit-time `integrate()` call … must produce the same answer -- **this is a claim about the CODE PATH, provable without a binary**"*. A tautology proved without a binary is a leak detector, not a feature gate. The tangent test's own mechanism note concedes the same algebra: *"`sigma~ = sigma_n + Ce*(d_eps - f*d_eps_p)` is EXACTLY affine in `d_eps` for FIXED `f`"* — affine for **any** f, so the secant returns `Ce` whatever `f` is. And `implexDetail[5]` = `mImplexFactor` — the one Python-visible observable of `f` — is unpacked in the floor-clamp test as `f_last` and then **never asserted on**.

Independent confirmation from the team's own history: `80e65a4de` fixed a defect in which "*every retried step in an adaptive run silently ran with the plastic extrapolation switched off*" — i.e. mutant (a), in production — and the commit says the battery "could never have seen" it. **No regression test was added for it** (the file has no test that asserts `f ≠ 0` after a failed-and-reverted step).

**CONSEQUENCE.** The only test that would go red on a fully inert `-implex` is `test_floor_clamp_fires_and_implexDetail_reports_it`, and per `17ccc6075` that test's deck was **changed until it fired** (`lat 0.5 → 1.5`, isochoric → net-dilating) after measuring `count == 0` on the real binary. The battery's live green surface for the actual IMPL-EX operator is one test whose deck is a tuned knob. Gate 2 (oracle parity) — the only thing that could show `sigma~` is *numerically right* — is admitted "still not run" in `70f6bf0e8`.

**WHAT WOULD REFUTE IT.** A run of `Ladruno_scripts/mutation_gate.py` on this family showing that zeroing `mImplexFactor`, or switching `ladrunoImplexTrial` to the total form, turns some test other than the floor clamp red. See F10 — that gate has not been run for ADR-92.

---

## F6 — MAJOR — the wire grew `Vector(5) → Vector(22)` and has **zero green coverage**; the one test that aims at it is blind by construction *and* reports SKIPPED

**Where:** `tests/test_ladruno_sanisand_implex.py:862-912`; `LEDGER_implementations.md` ADR-92 P1 row ("Wire grows `Vector(5) -> Vector(22)`").

**CLAIM.** `test_implex_db_roundtrip_carries_flags_and_history` asserts `reldiff(ref, out) <= 1e-12` where `ref = sani._drive(..., opts_on)` and `out` is the post-restore run. `sani._drive` builds `sani._build`, which is **zero-free-DOF** (every DOF is either `ops.fix` or `ops.sp`, `test_ladruno_sanisand.py:349-371`). By this same file's gate-5 argument, `-implex` ON and OFF commit **bit-identical** stress on such a deck. Therefore a restore that silently dropped every `-implex` flag and the whole `mImplexDEpsP` history would land on **exactly** `ref` and the assertion would pass.

**EVIDENCE.** The file's own non-vacuity guard proves it:

```python
    gap = sani._reldiff(ref_off, ref)
    if gap < _SENSITIVITY_FLOOR:
        pytest.skip('this deck does not distinguish -implex on from off ...')
```

and the run records confirm it fires: `8fe4f5630` reports **"13 passed / 2 failed / 1 skipped"**, `17ccc6075` reports "11 passed / 4 failed / 1 skipped". One test skips on every run. It is this one (the only other `skip` sites are `database()` unsupported, which would skip *before* the two assertions).

**CONSEQUENCE.** The ADR-92 P1 ledger row claims the six-override discipline and the new history variable survive `sendSelf`/`recvSelf`, sourced to a test that (i) reports SKIPPED and (ii) could not detect the failure even if it reported PASSED. 17 new wire slots (including `mImplexDt0`, see F2, and `mImplexFactor`, see F1) ship untested. No MPI `sendSelf`/`recvSelf` test exists at all.

**WHAT WOULD REFUTE IT.** A deck in this file where `-implex` ON and OFF give different **committed** answers — the file argues in two separate places that no zero-free-DOF deck can, and every ON/OFF pair in it is zero-free-DOF.

---

## F7 — MAJOR — `3c788778f` changed shipped D2 semantics **after** the "13/15" measurement, with no test, no re-run, and no doc update; one test's name and docstring now describe behaviour that no longer exists

**Where:** commit `3c788778f` (2026-09-06 01:53:58) vs `8fe4f5630` (01:15:46); `tests/test_ladruno_sanisand_implex.py:720-762`; `LEDGER_implementations.md` ADR-92 P1 row.

**CLAIM.** The guard moved from `mImplexDt < 0.0` to `mImplexDtCommit != 0.0 && mImplexDt * mImplexDtCommit < 0.0`. No commit after `3c788778f` touches `tests/`, so:

1. **The "13/15" figure predates the change.** It was measured on `8fe4f5630`'s binary, two behaviour commits ago. No post-fix pytest result is recorded anywhere in the range.
2. **The test name and docstring are now false.** `test_implex_negative_pseudo_dt_is_refused_and_leaves_committed_state_unchanged`: *"a NEGATIVE pseudo-time increment (`ops_Dt < 0`) is refused"* — a negative increment is now explicitly **legal** and the fix's whole point is that it must not be refused.
3. **The new behaviour has no test.** `3c788778f` states the regression's cause exactly — *"Every unit test drove LoadControl POSITIVE, which is why 13/15 green did not catch this"* — and then adds no test that drives it negative. `grep` confirms: no `LoadControl(-` in the battery except the single sign-reversal probe.
4. **The docs still describe the old contract.** `LEDGER_implementations.md`'s P1 row: *"a **negative** pseudo-time refused by symptom (a load factor turned round past a limit point)"*.

**CONSEQUENCE.** The most consequential class of deck for this feature (monotone negative clock — the campaign's own) is unrepresented in the battery both before and after the fix, which is how F1 and F2 got through.

**WHAT WOULD REFUTE IT.** A pytest transcript dated after 01:53:58, or a test that constructs `dt < 0` on two consecutive steps and asserts `analyze() == 0`.

---

## F8 — MAJOR — the gate's cost model was edited in the same commit as the results, after the number came back unfavourable, and the edit is **inert**

**Where:** `git show 70f6bf0e8 -- Ladruno_files/testbed/hypo_bearing/adr92_bvp_gate.py`.

**CLAIM.** Thresholds *were* pre-registered (`5e5ec4db7`, 00:44:59, adds `REFUTE_PAST_RUNG1 = 40.0`, `PASS_PAST_RUNG1 = 10.0`, `PASS_FAILED_SHARE = 10.0`) — that part of the process holds. But the **cost model** was changed inside the results commit itself:

```diff
-    fail_it = ns * sum(RUNG_CAPS) + rung2 * RUNG_CAPS[0] + ...
+    ladder_fails = max(ns - refused, 0)
+    fail_it = (ladder_fails * sum(RUNG_CAPS) + refused * RUNG_CAPS[0] + ...)
```

with the comment *"measured 2026-09-06, arm A came back 85.0 % ... which is impossible as a ladder cost and is an artefact of this line."*

**EVIDENCE.** `refused` comes from `d.get('n_material_refused', 0)`. **No producer writes that key** — `grep -rn n_material_refused` across the repo returns only the gate itself and two "Owed"/"still not recorded" lines in `_adr92_p1_bvp_gate_results.md`; not one of the three committed leg JSONs contains it. So the correction evaluates `refused = 0` and reduces to the original expression, byte for byte.

**CONSEQUENCE.** A metric was reformulated post-hoc to answer an inconvenient number, shipped in the same commit as that number, and does nothing. Arm A's reported `failshare` is still 85.0 %, which is above `PASS_FAILED_SHARE`, so arm A does **not** pass the gate's own criteria — yet the commit headline is "the ladder is GONE". The datum can in fact be reconstructed from the *committed* evidence: `implex_ctl`'s engine log has exactly 75 `Domain::update - domain failed in update` lines = `nfail`, and `nsub = 25`, `3 × 25 = 75`, so every failed rung belongs to a subdivided step; substituting `refused = 25` gives `failshare = 53.2 %` — still failing, but a number, not a disclaimer.

**WHAT WOULD REFUTE IT.** Any producer writing `n_material_refused`, or a leg JSON containing it.

---

## F9 — MAJOR — the `getCopy` test is vacuous: the clones are made before any history exists, and the implementation copies no history for *any* material

**Where:** `tests/test_ladruno_sanisand_implex.py:914-998`; `SRC/material/nD/LadrunoSANISAND.cpp:809-843`.

**CLAIM.** The test claims to verify that `getCopy(const char*)` "deliberately does NOT transfer `mImplexDEpsP`". It cannot.

**EVIDENCE.** Two facts. (1) `getCopy(const char*)` constructs a brand-new `LadrunoSANISAND3D` from **scalar parameters only** and then calls `clone->setLadrunoImplexOptions(mImplexOpt, false)` — there is no history-copying code to break; a mutant that *tried* to transfer `mImplexDEpsP` would have to be written from scratch. (2) Both `ops.element('stdBrick', 1, ...)` and `ops.element('stdBrick', 2, ...)` run **before** the first `ops.analyze`, so at copy time the prototype's `mImplexDEpsP` is identically zero. The assertion `all(abs(v) == 0.0 for v in bystander_pstrain)` is satisfied by a material that sits at zero strain — which any correct material does, `-implex` or not.

**CONSEQUENCE.** Zero coverage of the getCopy seam's actual risk (the `LEDGER_implementations` claim that the implex options *are* transferred, `setLadrunoImplexOptions(mImplexOpt, false)`), and a green test that reads as coverage in the ledger. The stronger construction — drive the *prototype* through a plastic history, then build a second element from the same tag — is the one the sibling `test_getcopy_void_carries_the_settings` uses and this test does not.

**WHAT WOULD REFUTE IT.** A code path where `ops.element(...)` clones from the prototype *after* the prototype has committed plastic strain, on this deck.

---

## F10 — MAJOR — ADR-87 warrant package: `LEDGER_implementations.md` row is stale and inaccurate; `manifest.yaml` not updated; mutation gate not run; the manifest CI check passes vacuously

**Where:** `Ladruno_implementation/LEDGER_implementations.md` (the one added row); `Ladruno_implementation/testbed/manifest.yaml`; `ci/check_manifest.py`.

**CLAIM.** Four distinct compliance gaps.

1. **Status is stale.** The row reads: *"**draft — C++ written, NOT YET BUILT.** No gate has run: oracle parity, tangent identity, floor clamp, byte-identity, symmetry, `dt` guards and the BVP gate all wait on the owner's `build.bat`."* Six later commits report gate results on a built binary, including the BVP gate whose findings changed the C++ twice.
2. **Files column is incomplete.** It lists only the three `SRC/` files. Missing: `tests/test_ladruno_sanisand_implex.py`, `Ladruno_files/testbed/hypo_bearing/adr92_bvp_gate.py`, `Ladruno_implementation/_adr92_p1_{execution_plan,bvp_gate_results}.md`, `LEDGER_quirks.md`, and the ~40 committed data files. `CLAUDE.md` requires the row to name the files.
3. **Behaviour claim is wrong** post-`3c788778f` (see F7).
4. **Manifest not updated, and the CI gate cannot notice.** `ci/check_manifest.py` keys on *new Ladruno `#define`s in `SRC/classTags.h`*. ADR-92 P1 adds **no new classTag** (the row says so), so `check_manifest` passes without ever looking at the feature. `manifest.yaml`'s `ND_TAG_LadrunoSANISAND` row (line 913) still names only `tests/test_ladruno_sanisand.py` as its pytest, and its cited coverage evidence is all *mutation-proven* ("PROVEN by mutation: … fails THIS TEST AND ONLY THIS TEST"). **No mutation proof exists for any ADR-92 mechanism** — which is precisely the instrument that would have caught F5, F6 and F9.

**CONSEQUENCE.** The required-check surface (`check_manifest`, Zone-A) is green while the feature is unregistered and unmutation-tested. Compliance is asserted by a ledger row that describes a state five commits old.

**Mitigation the blue team is entitled to:** the PR is still a **draft**, and ADR-87 requires the complete warrant package only at the ready-flip. On the other draft-day-one requirement the branch is **compliant**: PR #798 `createdAt = 2026-09-06T04:28:02Z` versus first commit `02b1a7c8b` at `2026-09-05 23:27:34 -0500` = `04:27:34Z` — **28 seconds**. Banner: correctly untouched (`banner_features.txt`, `tclMain.cpp`, `PythonModule.cpp` all absent from the range), matching the row's "No banner line until CP2". `LEDGER_vanilla_files.md`: correctly untouched — all four `SRC/` files in the range (`LadrunoSANISAND.{h,cpp}`, `LadrunoSANISAND3D.cpp`, `LadrunoSANISANDPlaneStrain.cpp`) are fork-owned, no upstream footprint, so no row is owed.

**WHAT WOULD REFUTE IT.** A mutation-gate transcript for ADR-92 P1, or a `manifest.yaml` diff in this range.

---

## F11 — MAJOR — two known failures are committed **red**, not `xfail`/`skip`, and the draft-PR policy hides them; the reported denominator is wrong

**Where:** `tests/test_ladruno_sanisand_implex.py` (no `xfail` anywhere — `grep "xfail"` returns nothing); `.github/workflows/build_cmake.yml:43-48` (`pytest -v` over `tests/`); `CLAUDE.md` ("a draft … the full Zone-A build intentionally skips drafts").

**CLAIM.** `test_tangent_identity_frozen_ce_on_strain_dt_source` and `test_implexcontrol_refuses_past_tolerance_and_leaves_committed_state_unchanged` are known-failing and carry no marker, so the branch's Zone-A pytest run is **red by construction** — and because the PR is a draft, the full Zone-A build is skipped and nobody sees it. `CLAUDE.md`: *"Do not work around a red check — fix it."* Leaving them unmarked on a branch whose CI does not run them is the same effect as working around it.

**EVIDENCE.** AST count of the file: 11 single tests + one `@pytest.mark.parametrize('scheme', [3,5,7,8,9])` = **16 collected tests**, not 15. The reported "13/15" is 13 of the 15 non-skipped; the file's own commit messages give the true shape ("13 passed / 2 failed / 1 skipped"). The task brief and the ledger both propagate "15".

**CONSEQUENCE.** The headline pass rate is stated against a denominator that excludes the skipped test (F6) — the one that hides a structural blindness — so the published figure both overstates coverage and understates the test count.

**WHAT WOULD REFUTE IT.** An `xfail(strict=True)` marker on either failing test, or a Zone-A configuration that deselects this file.

---

## F12 — MINOR — tuned constants presented as controlled variables; a stated precondition never checked; 12.4 MB of superseded log committed

**Where:** `tests/test_ladruno_sanisand_implex.py:141-146, 341-420, 470`; `LEDGER_quirks.md` (new entry, item 2); `adr92_bvp/implex_ON/`.

Three smaller items, each concrete.

1. **`_CAP_ADEQUATE = 20000` is a knob that was read off this build.** The block comment is admirably honest — *"Measured on this build (2026-09-06) … this deck needs between 1000 and 5000 ModifiedEuler substeps"* — but that makes it a value derived from the binary under test, not from the physics. Making it a *controlled* variable (same on both arms of every ON/OFF pair) was the right fix and is correctly applied; the residual objection is that no test pins the substep count itself, so a regression that changed the required count from 1000-5000 to 25000 would silently convert every ON/OFF pair into an unequal-cap comparison again — the exact bug `5e5ec4db7` fixed.
2. **The tangent-identity gate never checks its own stated precondition.** The team's own new quirk entry says: *"A tangent-identity gate must therefore be run on a path where the clamp is idle, and a test that reports 'identity to machine precision' on a clamped path is measuring nothing."* `_PROBE_P0 = 100.0` is chosen "away from the p_min floor" by comment only — the test never reads `implexDetail[3]`/`[4]` during the two probes to assert the clamp stayed idle, although the floor-clamp test proves that observable is available.
3. **Repo hygiene.** `adr92_bvp/implex_ON/a2_h1.0_e0.6944_engine.log` is **12,407,322 bytes / 34,396 lines**, of which **34,272 lines are one repeated warning string** (verified by normalising digits and counting unique line shapes). It is the log of the **superseded pre-fix** run — the one `3c788778f` fixed — and its directory also holds a 2-line curve CSV and a `.FAILED.txt`. Committing campaign data *is* this repo's convention (`adr90_a2/` etc. are tracked), but a 12 MB dead-run log is not: the sibling post-fix arms' logs are 394 and 8 lines. A 40-line excerpt plus the count would carry the same evidentiary weight. Related and worth fixing at the source: the `-implexControl` refusal at `LadrunoSANISAND.cpp:1615` returns `LADRUNO_MATERIAL_REFUSED` with **no `opserr` at all** (the ctl arm's log contains zero material-refusal lines; only `LadrunoBrick`'s own throttled-at-10 message), which is why `n_material_refused` is unrecoverable — while the D2 guard, which needs no such census, shouts 34,272 times unthrottled.

**WHAT WOULD REFUTE IT.** A `.gitignore` rule or LFS pointer for `*_engine.log`; an `opserr` on the `-implexControl` refusal path.

---

# COVERAGE MATRIX — feature aspect × test

| # | Feature aspect | Covered by | Verdict |
|---|---|---|---|
| 1 | `-implex` off is byte-identical | `test_implex_unset_reproduces_recorded_sensitivity` | **real** — recorded ADR-86b number, genuine regression check |
| 2 | `-implex` on does not corrupt the **committed** answer | `test_implex_on_matches_off_on_a_zero_free_dof_deck` | **leak detector only** — tautological by its own docstring; survives every operator mutant (F5) |
| 3 | Stage-0 (`mElastFlag == 0`) inertness + `LoadControl(0.0)` hold | `test_stage0_inertness_gravity_and_hold_is_bit_identical` | **real**, but same zero-free-DOF blindness — passes with the operator inert |
| 4 | `p_min` clamp on `sigma~` fires | `test_floor_clamp_fires_and_implexDetail_reports_it` | **real but tuned** — deck reshaped `lat 0.5 → 1.5` after measuring `count == 0` (`17ccc6075`) |
| 5 | `implexDetail` split identity `‖dσ‖² = ‖ddev‖² + 3(dp)²` | same test | **weak** — a norm-decomposition identity true for *any* vector, satisfied by `d = 0` |
| 6 | Frozen tangent `d(σ~)/dε == Ce` on `-implexDt strain` | `test_tangent_identity_frozen_ce_on_strain_dt_source` | **RED (unmarked)** + would pass on `f≡0` and on the total form (F5); precondition unchecked (F12.2) |
| 7 | Parser refusal: IntScheme 3/5/7/8/9 | `test_implex_refuses_unsupported_schemes` ×5 | **real** |
| 8 | Parser refusal: `-maxSubsteps` required / `0` refused | `..._scheme1_requires_maxsubsteps`, `..._maxsubsteps_zero_is_also_refused` | **real** |
| 9 | D2 sign-change refusal (**shipped** semantics) | — | **NOT COVERED** — the one test present asserts the *old* `dt < 0` semantics (F7) |
| 10 | D2 must **not** refuse a monotone negative clock | — | **NOT COVERED** — the regression `3c788778f` fixed; still untested (F7) |
| 11 | `-implexControl` refuses past tolerance | `test_implexcontrol_refuses_past_tolerance...` | **RED (unmarked)**; even green it only asserts `rc != 0` with no positive control separating refusal from ordinary non-convergence |
| 12 | Refusal leaves committed state unchanged | both refusal tests, `before == after` | **vacuous** — a failed `StaticAnalysis` step already ran `revertToLastCommit()`; true for any failure cause |
| 13 | `sendSelf`/`recvSelf` — 17 new wire slots | `test_implex_db_roundtrip...` | **SKIPPED + structurally blind** (F6) |
| 14 | Parallel (MPI) `sendSelf`/`recvSelf` | — | **NOT COVERED** |
| 15 | `getCopy(const char*)` transfers the options / not the history | `test_implex_getcopy...` | **VACUOUS** (F9) |
| 16 | `getCopy(void)` under `-implex` | — | **NOT COVERED** (the sibling file covers it for the ADR-86 settings only) |
| 17 | **PlaneStrain** lane (`LadrunoSANISANDPlaneStrain.cpp`, +20 lines in range) | — | **NOT COVERED** — every deck in the file is `ndm 3` (4 occurrences of `'-ndm', 3`, zero of `'-ndm', 2`) |
| 18 | First-step arming / `f` frozen once per step | — | **NOT COVERED** — `implexDetail[5]` (`mImplexFactor`) is read as `f_last` and never asserted; the `80e65a4de` arm-consumption bug has no regression test |
| 19 | `revertToLastCommit` ghost-`update()` arm consumption | — | **NOT COVERED** (the defect that motivated `80e65a4de`) |
| 20 | Companion (commit-time) return **failure** — `-maxSubsteps` exhausted inside `ladrunoImplexCommit` | — | **NOT COVERED** — the case the source calls "there is no good move left" |
| 21 | Subdivision / retry interaction with the frozen `f` | — | **NOT COVERED** in pytest; exercised only by the BVP gate, whose metric is blind to it (F3) |
| 22 | `-implexAlpha` | — | **NOT COVERED** — flag parsed at `LadrunoSANISAND.cpp:318`, never passed by any test |
| 23 | `-implexDt user` | — | **NOT COVERED** |
| 24 | `-implexDt pseudo` explicitly (vs by default) | — | **NOT COVERED** (only `'strain'` is ever passed explicitly) |
| 25 | `dt == 0` hold ⇒ `f = 0` (the `80e65a4de` rule) | `test_stage0_...` reaches `LoadControl(0.0)` but at **stage 0**, where `-implex` is inert by design | **NOT COVERED at stage 1** |
| 26 | Oracle parity — `sigma~` is numerically *correct* | — | **NOT COVERED** — gate 2 admitted "still not run" (`70f6bf0e8`); and per `LEDGER_quirks` the P0 oracle has no clamp, so it is not a reference on a clamped path anyway |
| 27 | Tangent symmetry claim (`TanType` inert under `-implex`) | — | **NOT COVERED** — asserted in the echo line and the ledger; no test reads the 6×6 and checks `A == Aᵀ` |
| 28 | Non-`LadrunoBrick` elements silently swallow the refusal | — | **NOT COVERED** — documented in `LEDGER_quirks`/ADR-86b as making a capped run "silently INVALID"; no test asserts the warning or a refusal |

**Summary:** of 28 aspects, 6 are genuinely covered, 5 are covered but non-discriminating (tautological, vacuous, tuned, or skipped), 2 are red-unmarked, and **15 are not covered at all** — including both PlaneStrain, both parallel, the whole `f`-freezing mechanism the feature's own quirk entry calls the central trap, and the numerical correctness of the extrapolated stress.
