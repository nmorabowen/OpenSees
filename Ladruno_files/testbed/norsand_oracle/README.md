# norsand_oracle: WP-144 LadrunoNORSAND reference oracles and gate G1 (Zone B)

The two independent Python oracles for the LadrunoNORSAND material (NorSand in the Andrade & Borja 2006 form),
both written from [`144a_norsand_equation_sheet.md`](../../../Ladruno_implementation/144a_norsand_equation_sheet.md)
by different authors and models. Plan: [`144_ladruno_norsand_plan.md`](../../../Ladruno_implementation/144_ladruno_norsand_plan.md) §5–6.

| Folder | What it is | Role |
|---|---|---|
| `o1_rate/` | The continuum rate equations, integrated by SciPy Radau (rtol 1e-10), with no return map | **The truth** |
| `o2_algo/` | Backward-Euler spectral return map (AB06 Box 2) + the closed-form consistent tangent | **The C++ contract**. Its constants (`kernel.py`, README) bind the P1 kernel |
| `tests/` | The gate G1 suite: K1 closed forms, O2→O1 convergence, tangents (small and finite strain), K2 benchmark (AB06 §6.1) + sensitivity, identities, cap, the Gudehus–Argyris census | Expected values come only from the sheet, published numbers or stated convergence arguments |
| `kernel_parity/` | P1a: the C++ kernel `SRC/material/nD/LadrunoNorSandKernel.h` against O2, step by step, through a ctypes shim built with g++ (`ns_shim.cpp`, `ns_kernel.py`), plus the g++ self-check driver `tests/ladrunonorsand_kernel_check.cpp` (the two P1a gates). `mutate_kernel.sh` re-runs it against 13 kernel mutants, all killed (nine kernel mutants plus the four chain mutants of the chained substep tangent: `last_substep_tangent`, `chain_drop_Spi_carry`, `chain_drop_v_column`, `chain_assemble_1e-6_off`). `tools/` holds `syntax_check.sh` (shell TUs, `-Wall -Wextra -Wshadow`), `syntax_check_reg.sh` (registration TUs) and `voigt_check.cpp` (the Voigt/tensor mapping of the shell against the real kernel by finite differences) | Gate 1e-10 relative (the argument and the corner/vertex-band tangent exception are in the test's docstring) |

`tests/conftest.py` maps the common parameter names onto each oracle (`make_params`, `K2_BASE`).

## Run

The suite needs numpy, scipy and sympy, so it is Zone B. Run it on Esmeralda (`~/wp144/venv`), from this folder:

```bash
python -m pytest tests kernel_parity -q -p no:cacheprovider
```

It takes about 6 minutes; G1 closed at 391 / 391 on 2026-10-01. `kernel_parity` adds 37 tests (~10 s;
it needs `g++` on PATH and skips without it).
`tests/out/` holds the generated K2 sensitivity table and the elastic-energy convexity table.
