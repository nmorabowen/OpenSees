"""Self-test for quirk-lint rule `pointers` (alias L3, WP-115) in ci/check_quirk_patterns.py.

One file per rule (WP-162), so two PRs adding rules never edit the same test file.
Shared helpers live in ci/_quirk_testkit.py.
Run: pytest -q ci/test_quirk_pointers.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _quirk_testkit as kit  # noqa: E402
import check_quirk_patterns as cq  # noqa: E402


def test_flags_an_orphaned_pointer(tmp_path):
    root = kit.tree(tmp_path, {
        "Ladruno_implementation/LEDGER_quirks.md": "### `wipe()` does NOT recreate the Domain\n",
        ".claude/skills/g/SKILL.md": 'Quirks: "`wipe()` does NOT recreate the Domain", "renamed heading".\n',
    })
    out = cq.check_pointers(root, kit.rel(root))
    assert len(out) == 1 and "renamed heading" in out[0]
