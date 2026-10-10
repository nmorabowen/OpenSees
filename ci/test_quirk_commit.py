"""Self-test for quirk-lint rule `commit` (alias L4, WP-118) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_commit.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


# the LadrunoIMKBeam shape (pre WP-118): commits materials + transform, never the base
IMK_COMMIT = kit.STAMP + """
int Beam::commitState(void)
{
  int ok = 0;
  for (int i = 0; i < 2; i++)
    ok += theMat[i]->commitState();
  ok += theCoordTransf->commitState();
  return ok;
}
"""


ELEMENT_H = "class Beam : public Element\n{\n};\n"


def _run(tmp_path, files):
    files = dict(files)
    files.setdefault("SRC/element/Beam.h", ELEMENT_H)
    root = kit.tree(tmp_path, files)
    used = set()
    return cq.check_commitstate(root, kit.rel(root), used) + cq.check_stale_waivers(root, kit.rel(root), used)


def test_ignores_non_elements(tmp_path):
    # materials / sections / transforms have a commitState() too, but own no Kc
    out = _run(tmp_path, {"SRC/element/Beam.cpp": IMK_COMMIT,
                         "SRC/element/Beam.h": "class Beam : public UniaxialMaterial\n{\n};\n"})
    assert out == []


def test_follows_indirect_inheritance(tmp_path):
    out = _run(tmp_path, {"SRC/element/Beam.cpp": IMK_COMMIT,
                         "SRC/element/Beam.h": "class Beam : public BeamBase\n{\n};\n",
                         "SRC/element/BeamBase.h": "class BeamBase : public Element\n{\n};\n"})
    assert len(out) == 1, out


def test_flags_a_commitstate_that_skips_the_base(tmp_path):
    out = _run(tmp_path, {"SRC/element/Beam.cpp": IMK_COMMIT})
    assert len(out) == 1 and "Beam::commitState does not call Element::commitState" in out[0], out


def test_passes_when_chaining_to_element(tmp_path):
    fixed = IMK_COMMIT.replace("  int ok = 0;", "  int ok = this->Element::commitState();")
    assert _run(tmp_path, {"SRC/element/Beam.cpp": fixed}) == []


def test_passes_when_chaining_to_a_parent(tmp_path):
    parent = IMK_COMMIT.replace("  int ok = 0;", "  int ok = BeamBase::commitState();")
    assert _run(tmp_path, {"SRC/element/Beam.cpp": parent}) == []


def test_exempts_a_class_that_overrides_rayleigh(tmp_path):
    # the LadrunoRigidBody shape: Rayleigh ignored, Kc never allocated
    in_cpp = IMK_COMMIT + "int Beam::setRayleighDampingFactors(double, double, double, double) { return 0; }\n"
    assert _run(tmp_path / "a", {"SRC/element/Beam.cpp": in_cpp}) == []
    in_h = {"SRC/element/Beam.cpp": IMK_COMMIT,
            "SRC/element/Beam.h": "class Beam : public Element\n{\n  int setRayleighDampingFactors(double, double, double, double);\n};\n"}
    assert _run(tmp_path / "b", in_h) == []


def test_ignores_member_commits_and_comments(tmp_path):
    # `theMat[i]->commitState()` is not a base call, and a commented base call does not count
    commented = IMK_COMMIT.replace("  int ok = 0;", "  int ok = 0;   // this->Element::commitState();")
    assert len(_run(tmp_path, {"SRC/element/Beam.cpp": commented})) == 1


def test_waiver_and_stale_waiver(tmp_path):
    waived = IMK_COMMIT.replace("int Beam::commitState(void)",
                                "// ladruno-lint: commit-ok Kc handled by an owned sub-element\n"
                                "int Beam::commitState(void)")
    assert _run(tmp_path / "a", {"SRC/element/Beam.cpp": waived}) == []
    stale = waived.replace("  int ok = 0;", "  int ok = this->Element::commitState();")
    out = _run(tmp_path / "b", {"SRC/element/Beam.cpp": stale})
    assert len(out) == 1 and "stale commit-ok" in out[0], out
