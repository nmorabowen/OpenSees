---
title: "WP-237 orchestration plan — work packages, roles, gates and protocol for certified limit analysis"
project: Ladruno
status: "Proposed (2026-11-18) — awaiting the owner's approval (§9). No program issue exists yet"
companion: 237_ladruno_limit_analysis_adr
updated: 2026-11-18
---

# WP-237 orchestration plan

ADR-87 requires that agent assignment and effort per work package be decided at kickoff and recorded.
`Ladruno_internal/ISSUES.md` requires that live state live on GitHub. This plan does the first and names the
homes for the second. Everything here is proposed until the owner approves it. The design is in
[[237_ladruno_limit_analysis_adr]]; where the two differ, the ADR wins.

## 0. Mechanics

- **Model and effort** are set per spawned agent when the tool exposes both. Check at each kickoff. If `effort`
  is gone, write "effort: harness default" in the handoff.
- **State is not kept here.** This file is static. Status, SHAs, job ids and "resume here" notes go to the GitHub
  homes in §4.
- **WP numbers are claimed by the program issue** (§8). Every number below is provisional until then.

## 1. Roles

| Role | Model · effort | Rule |
|---|---|---|
| Orchestrator: sequencing, verification, decisions, launching builds and long runs | current session · high | edits no `SRC/`, `tests/` or apeGmsh code |
| Contract author (237a) | Opus · high | the orchestrator freezes the result |
| Python reference author (237b) | Opus · high | does not read the checker |
| Independent checker author (237c) | Fable · high | does not read the reference or the engine |
| FEniCSx reference author (237e) | Opus · high | outside the fork; does not read fork limit code |
| C++ engine author (237f onward) | Opus · high | reads the interface and cases, not `tests/limit_reference.py` |
| Solver and scale measurements (237d) | Sonnet · medium | the orchestrator launches long runs |
| apeGmsh plumbing under an accepted ADR | Sonnet · medium, as T arranges | through T only |
| Critic at each ready-flip | the other model family from the author · high | briefed with the diff and the contract, not with conclusions |
| Read-only sweeps | Explore · very thorough | no edits |

The reference, checker, FEniCSx and C++ authors are four different agents. Whether ADR-87 D4 also needs an
out-of-family verdict is open (ADR O8).

`/code-review` levels: `low` for documents, `medium` for Python references, `high` for contracts and C++.

## 2. Work packages

### 2.1 Fork

| WP | Phase | Delivers | Owned paths | Agent | Starts when | Done when |
|---|---|---|---|---|---|---|
| 237 | setup | ADR and this plan | `Ladruno_implementation/237_*`, `_wp237_*` | Opus · high | owner go | owner accepts |
| 237a | P0 | frozen mathematics, data contract, cases, including the true-3D gate if A1 is accepted | `tests/limit_interface.md`, `tests/limit_cases.txt`, `tests/data/limit/`, contract tests | Opus · high | owner go | orchestrator freezes it, after O16 is answered |
| 237b | P0 | Python reference assembler | `tests/limit_reference.py`, `tests/test_limit_reference.py` | Opus · high | 237a frozen | L0 and L3 reference gates |
| 237c | P0 | independent checker and negative fixtures | `tests/limit_certificate.py`, `tests/test_limit_certificate.py`, fixtures | Fable · high | 237a frozen | L1 gates |
| 237d | P0 | solver viability, locking and scale study | `Ladruno_files/testbed/limit_analysis/` | Sonnet · medium | 237b and 237c usable | L2, 3D measurements, K1 and the locking gate decided |
| 237f | P1 | C++ kernel, assembler, command, export, writer, adapter | new limit directory, command row, CMake, ledger fragments, manifest row | Opus · high | I1; O2, O3, O4, O8, O9 settled | P1 gates and warrant package |
| 237g | P3 | V-H-M and 3D query workflow | query code and tests | Opus · high | 237f and bridge usable | P3 gates |
| 237h | P4 | adaptive 3D study | study driver and gates | re-decided at kickoff | I3 | P4 gates |
| 237i | P5 | semidefinite cones, facet output, and jumps if not already moved earlier | as carded | re-decided at kickoff | P4 stable, or owner reprioritises | P5 gates |
| 237j | P6 | performance and advanced meshing | as carded | re-decided at kickoff | measured need | measured gain over P4 |

### 2.2 Outside the fork

| WP | Delivers | Location | Starts when | Done when |
|---|---|---|---|---|
| 237e | FEniCSx assembly of the same cases, solved with open solvers | outside the fork and apeGmsh (O11) | 237a frozen; O11 answered | same problem class and exact cases agree within gates |

### 2.3 apeGmsh, through T

| Link | Delivers | Starts when | Done when |
|---|---|---|---|
| L1 | apeGmsh limit-analysis ADR, Proposed | fork contract draft exists; T accepts routing | maintainer ratifies |
| L2 | `ops.limit`, sidecar writing, stock refusal | L1 accepted; 237f contract stable | parity gates |
| L3 | `Results.limit` and viewer support | 237f output contract stable | reader and fixture gates |
| L4 | `FORK_PIN` move and `live-fork` proof | a fork release contains the engine | fork-sensitive tests green |
| C0 | ADR 0136 S0 with amendments A1 to A5 and estimator decoupling | owner and T route apeGmsh#3052 | maintainer accepts |
| C1 | S1: `results.error.indicator`, ZZ estimator | C0 | patch and Kirsch gates, indicator tests |
| C2 | S2: `g.mesh.adapt.size_field`, list-data writer, Min refresh | C0 | 2D and 3D size-field gates |
| C3 | S3: fresh-session driver | C2 and generic estimator ready | 2D gates, declaration gates, generic estimator test |

ADR 0136's A1 to A5 are not this ADR's A1 to A4. Always write "0136-A1" or "237-A1" in handoffs.

### 2.4 Restart lane

| Link | Slice | Delivers | Starts when | Done when |
|---|---|---|---|---|
| D0 | ADR 236 S0 | owner decisions on its §6 | now | ADR accepted |
| D1 | 236 S1 | `getHistory`, `setHistory`, commands, first three materials | D0 | its INV gates |
| D2 | 236 S2 | SANISAND family and PM4Sand | D1 | INV-1 on that family |
| D3 | 236 S3 | balancing-step measurements | D1 | numbers in ADR 236 |
| D4 | 236 S4 | bulk history in `.ladruno` | D1 | format gates |
| D5 | 236 S5 | footing restart with contact reset | D1 to D3 | same-mesh restart gate |
| D6 | 0136 S4 | `FORK_PIN` containing D1 and D2 | fork release | pin proof |
| D7 | 0136 S5 | transfer bundle and fork commands | C3, D6 | J2 restart gate |
| D8 | 0136 S6 | transient transfer | D7 | wave-energy gate |

D1 is coordinated with the MPM program, because ADR 232 phase 9 waits for the same primitive. No limit-analysis phase waits
for this lane.

## 3. Lane dependencies

```text
Docs:     237 ADR + this plan
Lane A:   237a -> 237b ----+
               -> 237c ----+-> 237d -> I1 -> 237f -> 237g
Lane E:   237a -> 237e (outside the fork; informs I1 if started)
Lane B:   237f contract -> I2 -> L1 -> L2, L3 -> L4 -> usable bridge for 237g
Lane C:   C0 -> C1
             -> C2 -> C3
Adaptive: 237g + C2 + C3 + 0136-A1 + 0136-A2 + 0136-A4 -> I3 -> 237h
Lane D:   D0 -> D1 -> D2, D3, D4; D1 to D3 -> D5; D1 + D2 -> D6; C3 + D6 -> D7 -> D8
Later:    237i; 237j
```

Safe in parallel once 237a is frozen: 237b, 237c and 237e. C0 and D0 are decision work and can proceed independently. C1 and
C2 can run in parallel if T gives disjoint files. Interacting fork PRs merge one at a time, with affected tests rerun on the
merged result.

## 4. GitHub homes for state

Nothing is posted until the owner approves the program issue.

| Kind | Home |
|---|---|
| Program status, WP claims, delegations, `human-gate` questions | one fork issue labelled `program` (to be created) |
| A WP's state and handoffs | its draft PR, using the header `### Handoff — WP-<nnn> — <date> — <model> @ <effort>` |
| Fork-wide open work | board #937, one row for the program |
| apeGmsh chain and rows | T coordinator #1202 and board #1203, arranged by T |
| A finding outside the WP | its own issue |
| A failure class seen twice | a `lesson` issue, then a lint rule or quirk within 7 days |
| Gotchas that must not be rediscovered | `Ladruno_implementation/ledger/quirks/` fragments |
| Long-run records | handoff comment: command, host, build, input hash, output path, start time |
| Private comparisons | outside public repositories; never in the fork or apeGmsh |

Private memory holds pointers only.

## 5. apeGmsh goes through T

- Every apeGmsh PR, merge, board row and label for this program is arranged with the T coordinator (#1202, #1203).
- If no T session is live, post on #1202 **after the owner approves the post**, and wait. Do not land work around it.
- Fork workers may draft apeGmsh text. They do not edit an apeGmsh checkout.
- ADR acceptance is the apeGmsh maintainer's act. No agent marks one accepted.
- Fork workers do not read apeGmsh source for this feature's fork code. ADR text may be read.

## 6. Gates

### 6.1 Integration

| Gate | Condition |
|---|---|
| I1, P0 to C++ | 237b, 237c and 237d passed; K1 did not fire; proposed addition under 237-A1: the true-3D locking question has an accepted answer |
| I2, fork to apeGmsh | problem and output contracts frozen, with fixtures both sides can read |
| I3, adaptive | 237g, C2, C3, 0136-A1, 0136-A2, 0136-A4 and a 3D size-field gate done |
| I4, restart | D1 and D2 in a released fork, and apeGmsh pins it |

0136-A3 joins I3 only if measurements show regeneration under the budget needs coarsening.

### 6.2 Every fork commit handed back

Run these from the worktree, using `python3.12`:

```bash
python3.12 ci/check_classtags.py
python3.12 ci/check_manifest.py
python3.12 ci/check_ladruno_commands.py
python3.12 ci/check_ledger_fragments.py
python3.12 ci/check_quirk_patterns.py
python3.12 ci/check_duplicate_defs.py
python3.12 -m ruff check --no-cache --select F811 ci tests Ladruno_scripts
python3.12 Ladruno_scripts/stamp_headers.py --check
python3.12 -m pytest -q ci/test_quirk_*.py ci/test_duplicate_defs.py
```

### 6.3 Ready-flip for a family PR

Oracles that are not self-referential; the mutation gate red under its mutants; the guide; both shipped engines or a recorded
absence; the manifest row, ledger fragments and banner line; a dispatched Zone-A run, because a ready-flip starts none; one critic
from the other family; and the ADR-87 D4 verdict if O8 requires it.

### 6.4 Results

A gate passes only on certified bounds. Every quoted result names build or commit, driver, output path, date, criterion, mesh, solver,
`domain_claim` and interface law. A certified Lower above a certified Upper stops the program. Kill criteria K1 to K5 are in ADR §18.

## 7. Protocol

1. Boot: read the program issue, the WP's last handoff and the board. Check `in-flight` and `WP-<n>` PR titles for duplicates.
2. Refresh refs from `origin/ladruno` and apeGmsh `origin/main`.
3. Dispatch one brief per worker, each with its own worktree and owned files. Independent workers never get each other's implementation.
4. Verify: read the named files and rerun the cheapest gate. A worker's report alone is not acceptance.
5. The orchestrator launches builds and long runs. Workers do not background builds.
6. A third fix round on one finding stops and goes to the orchestrator.
7. Land: the owner merges fork PRs unless the program issue records a delegation (O12). Merge commits only; never `--auto` or `--admin`.
8. Record: handoff on the PR, one line on the program issue, board row refreshed.

### Standing worker rules

1. Work only in your worktree and owned files. Report its path and head SHA.
2. The design is frozen. Objections are proposed amendments, not silent deviations.
3. Open code or ADR text before claiming that something exists, works or is missing. Cite `path::symbol`.
4. Do not open GPL, LGPL, CC-BY-SA or unlicensed implementation source for this feature, including DOLFINx, UFL, `dolfinx_optim`,
   `fenics_optim`, ECOS and `bearing-bounds`. Use papers and public mathematics.
5. No licence-requiring component in the fork. No code moves from the FEniCSx lane into the fork.
6. Use ledger fragments, never `LEDGER_*.md`.
7. Do not push, comment or open a PR unless the brief says the orchestrator has approved it.
8. Keep private project data, private paths and names out of public files. Use published analytic references.
9. Never tune an element, tolerance or mesh rule to reach a reference.

## 8. WP numbers

- **The program issue claims the numbers.** Its title carries the range, as earlier program issues do, for example
  "Program: certified limit analysis (WP-237…)".
- Until it exists, 237 and every suffix in this plan are provisional. If the number is taken first, the documents are renumbered before
  any push.
- Suffixes and new numbers are added through that issue, never by guessing the next free one.
- `Ladruno_internal/ISSUES.md` describes program issues but has no reservation rule. This rule is program practice, as in ADR 232.
- **This WP does not create the issue.** That needs the owner's approval.

## 9. Go / no-go items for the owner

None is settled by this plan.

1. Approve §1 and §2, or reassign.
2. Approve opening the program issue claiming WP-237 onward, its row on #937, and the first post to T on #1202.
3. Answer O16: accept, change or reject amendments 237-A1 to 237-A4. 237a should not freeze before this.
4. Confirm working defaults O1 to O5.
5. Decide O8 and O12.
6. Answer O11 before 237e starts.
7. Decide whether this ADR and plan get a critic's review before the first push.
