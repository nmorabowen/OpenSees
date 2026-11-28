"""Self-test for quirk-lint rule `reaction-scatter` in ci/check_quirk_patterns.py.

The threaded reaction loop (Domain::ladrunoThreadedReactions) scatters form-safe elements
with Element's base rule, so a class overriding addResistingForceToNodalReaction must opt
out with ladrunoThreadSafeReaction() == false.
Run: pytest -q ci/test_quirk_reaction_scatter.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402

OVERRIDE_H = """class Absorb : public Element
{
  public:
    int addResistingForceToNodalReaction(int flag);
};
"""


def _run(tmp_path, files):
    root = kit.tree(tmp_path, files)
    return cq.check_reaction_scatter(root, kit.rel(root))


def test_flags_an_override_without_the_opt_out(tmp_path):
    out = _run(tmp_path, {"SRC/element/Absorb.h": OVERRIDE_H})
    assert len(out) == 1 and "Absorb.h:4" in out[0] and "ladrunoThreadSafeReaction" in out[0], out


def test_passes_with_the_opt_out(tmp_path):
    fixed = OVERRIDE_H.replace(
        "    int addResistingForceToNodalReaction(int flag);",
        "    int addResistingForceToNodalReaction(int flag);\n"
        "    bool ladrunoThreadSafeReaction(void) const { return false; }")
    assert _run(tmp_path, {"SRC/element/Absorb.h": fixed}) == []


def test_ignores_the_base_and_comments(tmp_path):
    out = _run(tmp_path, {
        "SRC/element/Element.h": OVERRIDE_H,
        "SRC/domain/subdomain/Subdomain.h": OVERRIDE_H,
        "SRC/element/Plain.h": "class Plain : public Element\n{\n  // addResistingForceToNodalReaction(int)\n};\n",
    })
    assert out == [], out


def test_the_fork_tree_is_clean():
    root = Path(__file__).resolve().parent.parent
    assert cq.check_reaction_scatter(root, lambda p: str(p)) == []
