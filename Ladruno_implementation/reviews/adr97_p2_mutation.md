---
title: ADR-97 P2 — mutation gate (eigenprojection rotation term)
project: Ladruno
status: complete
owner: nmora
tags:
  - implementation
  - material
  - review
---

# ADR-97 P2 mutation gate

**PR** [#824](https://github.com/nmorabowen/OpenSees/pull/824) ·
**report** [[reviews/adr97_p2_report]] · **script**
`Ladruno_implementation/adr97_scripts/mutate_p2.py`

## The mutation

In the 6D back-transform of the principal-space Koiter tangent
(`ASDPlasticMaterial3D::cp_principal_return`), the eigenprojection **rotation
term** — the shear-slot coefficients

```
T[3+s, 3+s] = (y_i - y_j) / (x_i - x_j)
```

that carry the fact that the return moves the principal **values** while the
principal **directions** follow the trial stress — is replaced by `1.0`, the
value it takes for an elastic step (`y == x`).

This is a plausible, subtle bug rather than a wrecking ball. It touches **only
the tangent**: the residual, the region classification and the returned stress
are all untouched, so every committed state stays exact and `Backward_Euler` is
not involved at all.

Only the **non-degenerate** branch is mutated; the l'Hôpital limit used when two
trial eigenvalues coincide is left alone.

```
python3.12 Ladruno_implementation/adr97_scripts/mutate_p2.py apply
Ladruno_scripts\build.bat OpenSees OpenSeesPy
cd tests && PYTHONPATH=../dist/bin python3.12 -m pytest test_adr97_p2_principal.py -q
python3.12 Ladruno_implementation/adr97_scripts/mutate_p2.py revert
```

The script round-trips **byte-identically** against git (`git status SRC` clean
after `apply` + `revert`) before any build is launched.

## Result — 4 killed, 34 survivors

| test | clean binary | mutant |
|---|---|---|
| `test_gate2_algorithmic_fd_on_the_face_region[face, sheared]` | rel_err **1.08e-08** | rel_err **1.29e+00** — **KILLED** |
| `test_gate2_algorithmic_fd_on_the_face_region[face, axis aligned]` | rel_err **8.40e-09** | the rig no longer converges (`analyze -> -3`) — **KILLED** |
| `test_gate2_algorithmic_fd_on_the_degenerate_edge_region` | rel_err **2.88e-11** | rel_err **8.74e-02** — **KILLED** |
| `test_gate2_negative_control_backward_euler_cannot_do_the_face_rig` | CP converges (`rc = 0`), BE does not | CP no longer converges either — **KILLED** |

The 129 % error on the sheared face rig is the cleanest single number: the
rotation term is not a correction, it is most of the shear block of `C_alg`.

**The degenerate-edge test dies too, and that was not the prediction.** The
mutation was written to spare it (its `s1 == s2` gap takes the l'Hôpital branch,
which is untouched), but the oedometric state has `s1 == s2 != s3`: only Voigt
slot 3, the `(1,2)` pair, is degenerate — slots 4 and 5 pair a degenerate
eigenvalue with `s3` and take the MUTATED branch. So the degenerate branch is
covered by exactly one of the three shear slots on that rig, and no test in this
file isolates it alone. Recorded rather than papered over; a test that isolates
the l'Hôpital branch would need a fully hydrostatic trial, where the tangent is
identically zero (the apex) and there is nothing to measure. **The l'Hôpital
branch is therefore verified only in combination**, which is an honest limit of
this gate.

## The 34 survivors

Every one is a case the mutation provably cannot reach, and the pattern is the
point:

* **all 6 gate-1a region pins** — 1.5e-16 or better, unchanged to the digit;
* **all 3 gate-1b path audits** and the rotation guard;
* **all 3 gate-1c MCTC oracles**, including the two bit-identical
  CP-vs-`Backward_Euler` comparisons;
* **all 4 gate-4 step-size rows**, the `MC_ds` finding in both directions, and
  CP ≡ BE at the apex (2.05e-16);
* **all 12 gate-6 refusal / degeneracy / strict-inertness tests**.

The residual and the region classification are untouched, so the committed
stress must stay exact — and it does, which is the evidence that the mutation
really is confined to the tangent.

One survivor **degraded without failing**: the MCTC iteration contrast went from
**22** global Newton iterations (3, 3, 4, 4, 4, 4) to **26** (3, 3, 5, 5, 5, 5)
against `Backward_Euler`'s 74. Its assertion is `sum(it_cp) <= sum(it_be)`, which
26 still satisfies. That is the right assertion for a robustness gate but it
means the MCTC contrast is **not** a mutation detector; the three finite-
difference tests are.

## After the revert

`mutate_p2.py revert`, rebuild, re-run: `tests/test_adr97_p2_principal.py`
**38 passed**, and the full ASDPlasticMaterial3D battery (19 files:
`test_adr84_*`, `test_asdplastic_*`, `test_adr94*` except the doc-rewriting
`test_adr94_matrix.py`, `test_adr97_*`) **172 passed, 2 skipped**.
