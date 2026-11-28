"""WP-265: the mutation gate's IDENT registry cannot silently fall back to "no floor".

Pure Python, no binary. Pins what WP-265 put in Ladruno_scripts/mutation_gate.py:

  * every family with live LADRUNO_MUTATE_* call sites has an IDENT floor > 0
    (the pre-WP-265 0.0 floors were declared, never measured);
  * the IDENT suite of CONTINUUM and SHELLMOD is the ordinary suite PLUS its
    tangent-vs-FD gate file, and that gate file stays OUT of the ZERO/SCALE
    suites (where K and the FD are mutated consistently, it could only survive);
  * the CPML suite still holds the K-vs-FD tests that make its own paths an IDENT
    suite (if one is renamed or deleted, CPML's IDENT floor loses its warrant).
"""
import importlib.util
import os

import pytest

pytestmark = [pytest.mark.zone_a]

_HERE = os.path.dirname(os.path.abspath(__file__))
_GATE = os.path.normpath(os.path.join(_HERE, "..", "Ladruno_scripts", "mutation_gate.py"))

GATE_FILES = {"CONTINUUM": "test_continuum_tangent_gate.py",
              "SHELLMOD": "test_shellmod_tangent_gate.py"}
CPML_FD_TESTS = {"test_ladruno_cpml_waveguide_mms.py": "def test_group1_opensees_finite_difference_tangent",
                 "test_ladruno_cpml_up.py": "def test_up_gate0_oracle_14_cases_and_fd_tangent"}


@pytest.fixture(scope="module")
def mg():
    spec = importlib.util.spec_from_file_location("mutation_gate", _GATE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_every_gated_family_has_a_nonzero_ident_floor(mg):
    for fam in sorted(mg.GATED_FAMILIES):
        floor = mg.FAMILIES[fam]["min_score"].get("IDENT", 0.0)
        assert floor > 0.0, f"{fam}: IDENT floor is {floor}; WP-265 requires > 0 for a gated family"


@pytest.mark.parametrize("fam", sorted(GATE_FILES))
def test_ident_suite_is_the_family_suite_plus_its_tangent_gate(mg, fam):
    base = set(mg.resolve_paths(fam))
    ident = set(mg.resolve_paths(fam, "IDENT"))
    gate = GATE_FILES[fam]
    assert os.path.isfile(os.path.join(_HERE, gate))
    assert gate in ident, f"{fam}: {gate} is not in the IDENT suite"
    assert base <= ident, f"{fam}: the IDENT suite dropped ordinary tests: {sorted(base - ident)}"
    for mode in ("ZERO", "SCALE"):
        assert gate not in set(mg.resolve_paths(fam, mode)), (
            f"{fam}: {gate} leaked into the {mode} suite, where it can only survive")


def test_cpml_suite_keeps_its_fd_tangent_tests(mg):
    paths = set(mg.resolve_paths("CPML", "IDENT"))
    for fname, marker in CPML_FD_TESTS.items():
        assert fname in paths
        with open(os.path.join(_HERE, fname), encoding="utf-8") as fh:
            assert marker in fh.read(), f"{fname}: '{marker}' is gone (CPML IDENT warrant)"
