# sanisand_reference — an independent SANISAND (DM04) reference integrator

WP-134. A Python oracle for the SANISAND model of Dafalias & Manzari (2004), *J. Eng.
Mech.* 130(6):622–634. It is written from the paper's equations, not from
`ManzariDafalias.cpp` or WP-128's `md_port.py`. It exists to answer one question for
WP-128/129: **what do the rate equations do from state X under dε**, with no
integrator error.

Full write-up: [`Ladruno_implementation/134_sanisand_reference_integrator.md`](../../Ladruno_implementation/134_sanisand_reference_integrator.md).

## What it does

- The DM04 rate equations are an ODE in pseudo-time t ∈ [0, 1] over one prescribed
  increment. They are integrated with SciPy's implicit **Radau IIA** at rtol 1e-10,
  with a scaled atol and a fixed-step finite-difference Jacobian.
- Every discontinuity of the right-hand side is an **event** that ends a smooth
  segment:
  - yield (f → 0⁺);
  - unloading (N → 0);
  - the loading-denominator sign;
  - the α_in reversal ((α−α_in):n → 0);
  - the two Macaulay kinks;
  - p → p_floor.

  The mode is decided again at each event.
- On the plastic branch the consistency condition holds identically. f at exit and the
  maximum |f| drift are reported.
- A point where N > 0 but H ≤ 0 has **no admissible rate solution**. The integrator
  stops there with `H_nonpositive`. It never maps that case to "elastic".
- It reports two measures of α against the bounding surface:
  - `rho_b` = ‖α‖ / √(2/3)·α^b(θ_n): WP-128's measure;
  - `rho_alpha`, measured in α's own direction: the geometric admissibility test.
- Drivers:
  - a single increment;
  - chains of increments;
  - mixed stress/strain control (drained or undrained triaxial, simple shear).

## Options (paper by default; each UW addition a separate switch)

| field | paper (default) | UW value | what |
|---|---|---|---|
| `d_factor` | False | True | U1: low-p dilatancy sigmoid (p < 0.05 P_atm) |
| `p_residual` | 0 | p_r | U2: p → p + p_r on the plastic side |
| `p_min` | 0 | 0.0101 | U3: G, K use max(p, p_min) |
| `g_void_ratio` | `current` | `initial` | U4: G with e_init instead of e |
| `void_ratio_law` | `current` | `initial` | U5: de = −(1+e_init) dε_v |
| `alpha_in_rule` | `paper` | `uw` | U6: reseat α_in on (α−α_in):n < 0, or once per increment on the elastic trial |
| `h_cap` | None | 1e10 | U7: h = 1e10 when \|(α−α_in):n\| < 1e-10 |
| `elastic_moduli` | `continuous` | `frozen` / `frozen_increment` | U8 / U9: K, G of the committed state on the elastic part, or on the whole increment (ModifiedEuler) |

Presets:
- `Options()`: the paper.
- `Options.uw()`: all UW additions, the RK45 comparator.
- `Options.uw_me()`: `uw` plus U9, the ModifiedEuler comparator.
- `ring.ring_variants()["uw_model"]`: the UW constitutive additions U1–U5 with the
  paper's α_in rule and continuous moduli. This is **the oracle a corrected C++
  integrator should reproduce**.

## Conventions

- Compression positive.
- Tensors are 6-vectors of tensor components in the order xx yy zz xy yz zx
  (σ, α, α_in, z).
- Strain input is Voigt with **engineering** shear. This is the
  `ladrunoSANISANDReplay -convention compressionPositive` convention.
- Units are kPa.

## Use

```python
import sys; sys.path.insert(0, "Ladruno_scripts")
from sanisand_reference import CAMPAIGN, Options, State, integrate
st = State.from_voigt(sigma, alpha, z, e, alpha_in)      # compression positive
res = integrate(st, [0, 1e-4, 0, 0, 0, 0], CAMPAIGN, Options())
res.status, res.f_end, res.max_rho_alpha, res.state, res.segments, res.reseats
```

CLI (prints JSON), run from `Ladruno_scripts/`:

```
python -m sanisand_reference reproducer                 # WP-128's smallest reproducer
python -m sanisand_reference ring --mesh b8 --row 0 --probe shear --delta 1e-5 [--options uw_me]
python -m sanisand_reference increment --sigma ... --alpha ... --e 0.7 --deps 0 1e-4 0 0 0 0
python -m sanisand_reference triaxial --set toyoura --p0 100 --e0 0.833 --undrained
```

Validation runners and their outputs are in `Ladruno_files/testbed/sanisand_reference/`:
- `run_paper.py`
- `run_crosscheck.py` + `analyse_crosscheck.py`
- `run_reproducer.py`
- `run_ring.py`

## Requirements

- numpy and scipy (tested with scipy 1.17, CPython 3.11).
- The C++ cross-check (`cxx.py`) runs the fork's `opensees.pyd` in a CPython 3.12
  subprocess (`_cxx_runner.py`, `-S`, manual paths, asserting `opensees.__file__`),
  because that interpreter has no scipy. The paths can be overridden with
  `SANISAND_REF_OPENSEES_BIN`, `SANISAND_REF_PY312` and `SANISAND_REF_SITE312`.
- The ring CSVs are read from the checkout. If they are absent, they are read with
  `git show` from `origin/wp/127-sanisand-replay-counters`.

Tests:
- `tests/test_sanisand_reference.py`: pure Python, about 15 s.
- `tests/test_sanisand_reference_crosscheck.py`: skipped without the WP-127 pyd.
