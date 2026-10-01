# norsand_oracle: WP-144 LadrunoNORSAND reference oracles and gate G1 (Zone B)

The two independent Python oracles for the LadrunoNORSAND material (NorSand in the Andrade & Borja 2006 form),
both written from [`144a_norsand_equation_sheet.md`](../../../Ladruno_implementation/144a_norsand_equation_sheet.md)
by different authors and models. Plan: [`144_ladruno_norsand_plan.md`](../../../Ladruno_implementation/144_ladruno_norsand_plan.md) §5–6.

| Folder | What it is | Role |
|---|---|---|
| `o1_rate/` | The continuum rate equations, integrated by SciPy Radau (rtol 1e-10), with no return map | **The truth** |
| `o2_algo/` | Backward-Euler spectral return map (AB06 Box 2) + the closed-form consistent tangent | **The C++ contract**. Its constants (`kernel.py`, README) bind the P1 kernel |
| `tests/` | The gate G1 suite: K1 closed forms, O2→O1 convergence, tangents (small and finite strain), K2 benchmark (AB06 §6.1) + sensitivity, identities, cap, the Gudehus–Argyris census | Expected values come only from the sheet, published numbers or stated convergence arguments |

`tests/conftest.py` maps the common parameter names onto each oracle (`make_params`, `K2_BASE`).

## Run

The suite needs numpy, scipy and sympy, so it is Zone B. Run it on Esmeralda (`~/wp144/venv`), from this folder:

```bash
python -m pytest tests -q -p no:cacheprovider
```

It takes about 6 minutes; G1 closed at 391 / 391 on 2026-10-01.
`tests/out/` holds the generated K2 sensitivity table and the elastic-energy convexity table.
