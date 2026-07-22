#!/usr/bin/env python3
"""Emit reproducible evidence for the four-transformation CMS oracle."""

from __future__ import annotations

import os
import sys

import numpy as np
from scipy import linalg as sla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ladruno_cms_reference import (  # noqa: E402
    assemble_global,
    build_spring_mass_elements,
    run_paper_hierarchical_cms,
)


N = 120
MODEL = build_spring_mass_elements(N, 1.0, 1.0, mass_type="consistent")
K_REF, M_REF = assemble_global(MODEL)
EIGVALS_REF, _ = sla.eigh(K_REF, M_REF)

print("=" * 88)
print("Yu et al. four-transformation CMS oracle | N=120 | consistent element mass")
print("=" * 88)

print("\n--- Full local bases: exactness of the four transformations ---")
full = run_paper_hierarchical_cms(MODEL, 4, 3, None, None)
for name, value in full["dimensions"].items():
    print(f"{name:>42}: {value:4d}")
print(f"{'max absolute eigenvalue error':>42}: {np.max(np.abs(full['eigenvalues'] - EIGVALS_REF)):.4e}")
print(f"{'max backward residual':>42}: {np.max(full['residuals']):.4e}")
print(f"{'max duplicate-interface jump':>42}: {full['max_duplicate_jump']:.4e}")
print(f"{'K congruence error':>42}: {np.max(np.abs(full['T'].T @ K_REF @ full['T'] - full['K_star'])):.4e}")
print(f"{'M congruence error':>42}: {np.max(np.abs(full['T'].T @ M_REF @ full['T'] - full['M_star'])):.4e}")

print("\n--- Truncation/enrichment evidence (first five requested modes) ---")
print(f"{'k_L2':>6} {'k_L1':>6} {'r_final':>8} {'reduction':>11} {'max rel eig err':>18} {'max residual':>14}")
for k_level2, k_level1 in [(5, 6), (8, 10), (10, 15), (12, 20), (15, 25)]:
    result = run_paper_hierarchical_cms(
        MODEL, 4, 3, k_level2, k_level1, num_modes=10
    )
    rel = np.abs(result["eigenvalues"][:5] - EIGVALS_REF[:5]) / EIGVALS_REF[:5]
    r_final = result["dimensions"]["final"]
    reduction = 100.0 * (1.0 - r_final / N)
    print(
        f"{k_level2:6d} {k_level1:6d} {r_final:8d} {reduction:10.1f}% "
        f"{np.max(rel):18.4e} {np.max(result['residuals'][:5]):14.4e}"
    )

print("\nThe report demonstrates algebraic correctness, not MPI scalability or a production tolerance.")
