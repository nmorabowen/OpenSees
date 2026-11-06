"""Self-test for quirk-lint rule `revert` (alias L10, WP-153) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_revert.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


BASE_H = ("class IncrementalIntegrator { public: virtual int revertToLastStep(void); };\n"
          "class TransientIntegrator : public IncrementalIntegrator {};\n")
DR_H = kit.STAMP + ("class Relax : public TransientIntegrator\n{\n  public:\n"
                "    int newStep(double dt);\n    int update(const Vector &U);\n"
                "  private:\n    Vector *Ut, *Vhalf;\n};\n")
DR_CPP = kit.STAMP + ("int Relax::newStep(double dt) { Ut->addVector(1.0, *Vhalf, dt); return 0; }\n"
                  "int Relax::update(const Vector &U) { return 0; }\n")


def _run(tmp_path, files):
    files = dict(files)
    files.setdefault("SRC/analysis/integrator/IncrementalIntegrator.h", BASE_H)
    root = kit.tree(tmp_path, files)
    used = set()
    return cq.check_revert(root, kit.rel(root), used) + cq.check_stale_waivers(root, kit.rel(root), used)


def test_flags_an_integrator_that_inherits_the_no_op(tmp_path):
    # the LadrunoDynamicRelaxation incident (WP-153 #899)
    out = _run(tmp_path, {"SRC/analysis/integrator/Relax.h": DR_H,
                         "SRC/analysis/integrator/Relax.cpp": DR_CPP})
    assert len(out) == 1 and out[0].startswith("L10 ") and "Relax" in out[0], out


def test_passes_the_override(tmp_path):
    fixed = DR_H.replace("    int update(const Vector &U);\n",
                         "    int update(const Vector &U);\n    int revertToLastStep(void);\n")
    assert _run(tmp_path, {"SRC/analysis/integrator/Relax.h": fixed,
                          "SRC/analysis/integrator/Relax.cpp": DR_CPP}) == []


def test_passes_an_override_inherited_from_a_fork_parent(tmp_path):
    parent = kit.STAMP + ("class Leap : public TransientIntegrator\n{\n  public:\n"
                      "    int revertToLastStep(void);\n  private:\n    Vector *Ut;\n};\n")
    child = DR_H.replace("class Relax : public TransientIntegrator", "class Relax : public Leap")
    assert _run(tmp_path, {"SRC/analysis/integrator/Leap.h": parent,
                          "SRC/analysis/integrator/Relax.h": child,
                          "SRC/analysis/integrator/Relax.cpp": DR_CPP}) == []


def test_ignores_stateless_integrators_non_integrators_and_vanilla(tmp_path):
    stateless = DR_H.replace("    Vector *Ut, *Vhalf;\n", "    double dt;\n")
    assert _run(tmp_path / "a", {"SRC/analysis/integrator/Relax.h": stateless,
                                "SRC/analysis/integrator/Relax.cpp": DR_CPP}) == []
    other = DR_H.replace("public TransientIntegrator", "public Element")
    assert _run(tmp_path / "b", {"SRC/element/Relax.h": other,
                                "SRC/element/Relax.cpp": DR_CPP}) == []
    assert _run(tmp_path / "c", {"SRC/analysis/integrator/Relax.h": DR_H.replace(kit.STAMP, ""),
                                "SRC/analysis/integrator/Relax.cpp": DR_CPP.replace(kit.STAMP, "")}) == []


def test_waiver_and_stale_waiver(tmp_path):
    waived = DR_H.replace("class Relax", "// ladruno-lint: revert-ok retry-aware: commit() snapshots the march\nclass Relax")
    assert _run(tmp_path / "a", {"SRC/analysis/integrator/Relax.h": waived,
                                "SRC/analysis/integrator/Relax.cpp": DR_CPP}) == []
    short = DR_H.replace("class Relax", "// ladruno-lint: revert-ok ok\nclass Relax")
    out = _run(tmp_path / "b", {"SRC/analysis/integrator/Relax.h": short,
                               "SRC/analysis/integrator/Relax.cpp": DR_CPP})
    assert len(out) == 1 and "too short" in out[0], out
    stale = waived.replace("    int update(const Vector &U);\n",
                           "    int update(const Vector &U);\n    int revertToLastStep(void);\n")
    out = _run(tmp_path / "c", {"SRC/analysis/integrator/Relax.h": stale,
                               "SRC/analysis/integrator/Relax.cpp": DR_CPP})
    assert len(out) == 1 and "stale revert-ok" in out[0], out
